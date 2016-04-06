from __future__ import division
import sys,os, time
sys.path.append('../../../')
from copy import deepcopy
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from shared_data_holder import shared_data_holder
from tools.userError import userError
from tools.globalVariables import *
from tools.pyomoSolverCreator import pyomoSolverCreator
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.fva import fva
from tools.fba.set_specific_bounds import set_specific_bounds
from imp import load_source
from coopr.pyomo import *
from coopr.opt import *
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class MUST_doubles(object):
    """
    This class defines the optimization problem to identify the MUST double sets 
    of reactions, i.e., MUST_LU (or MUST_UL), MUST_LL and MUST_UU

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: March 23, 2016
    """
    def __init__(self,model, product_exchrxn_id, MUST_double_type, flux_bounds_ref, flux_bounds_overprod, growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, product_targetYield_percent = 80, min_biomass_percent = 10, blocked_rxns = [], MUST_X = [], MUST_L = [], MUST_U = [], build_new_optModel = True, objective_thr = 5, validate_results = False, results_filename = '', create_new_results_file = True, optimization_solver = default_optim_solver, warnings = True, stdout_msgs = True): 
        """
        All inputs are a subset of those for OptForce

        INPUTS:
        ------
               shared_data: An instance of class shared_data_holder
          MUST_double_type: Tyype of MUST doulbe to find. Allowed choices include: MUST_LU 
                            (or MUST_UL), MUST_LL and MUST_UU
        build_new_optModel: If True ia new optimization model is created
             objective_thr: If the outer objective function is less than this threshold the
                            code stops
          validate_results: If True each obtained solution is validated before finding the next
                            one
          results_filename: A string containing the name of the file where the results should
                            be written in
          create_new_results_file: If True the results are written in a new file and if False they 
                            are appended to an existing file specified by results_filename
        
        MUST_LU:  Max_OP(v1 - v2) < Min_WT(v1 - v2)
                                           v1 - v2
        Wild-type                          |-----|
        over-producing             |-----|

        MUST_UL: Min_OP(v1 - v2) > Max_WT(v1 - v2)
                            v1 - v2
        Wild-type           |-----|
        over-producing               |-----|

        MUST_UU: Min_OP(v1 + v2) > Max_WT(v1 + v2)
                            v1 + v2
        Wild-type           |-----|
        over-producing               |-----|

        MUST_LL: Max_OP(v1 + v2) < Min_WT(v1 + v2)
                                           v1 + v2
        Wild-type                          |-----|
        over-producing             |-----|


        Therefore, the objective of function of the inner problem is maximized for MUST_LU and MUST_LL and it is 
        minimized for MUST_UL and MUST_UL

        Ali R. Zomorrodi - Segre lab @ BU
        Last updated: 01/04/2016
        """
        # model
        self.model = model
        self.biomass_reaction = model.biomass_reaction
        self.biomass_rxn_id = model.biomass_reaction.id

        # Flux bounds in the reference and overproducing strains
        self.flux_bounds_ref = flux_bounds_ref
        self.flux_bounds_overprod = flux_bounds_overprod
 
        # Growth condition
        self.growthMedium_flux_bounds = growthMedium_flux_bounds

        # Product and biomass yields
        self.product_exchrxn_id = product_exchrxn_id
        self.product_targetYield_percent = product_targetYield_percent
        self.min_biomass_percent = min_biomass_percent

        # List of blocked rxns
        self.blocked_rxns = blocked_rxns

        # MUST singles
        self.MUST_X = MUST_X
        self.MUST_L = MUST_L
        self.MUST_U = MUST_U

        # Type of MUST doubles to identify
        self.MUST_double_type = MUST_double_type

        # build_new_optModel
        self.build_new_optModel = build_new_optModel

        # outer objective threshold
        self.objective_thr = objective_thr

        # validate results
        self.validate_results = validate_results

        # results file name
        self.results_filename = results_filename

        # create_new_results_file
        self.create_new_results_file = create_new_results_file

        # optimization_solver
        self.optimization_solver = optimization_solver

        # warnings
        self.warnings = warnings

        # stdout_msgs
        self.stdout_msgs = stdout_msgs

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
         attr_name: Attribute name
        attr_value: Attribute value
        """
        # model 
        if attr_name == 'model' and not isinstance(attr_value,model):
            raise TypeError('model must be instance of class model')

        # product_exchrxn_id 
        if attr_name == 'product_exchrxn_id' and not isinstance(attr_value,str):
            raise TypeError('product_exchrxn_id must be a string')

        # flux_bounds_ref, flux_bounds_overprod
        if attr_name in ['flux_bounds_ref','flux_bounds_overprod'] and not isinstance(attr_value,dict):
            raise TypeError('{} must be a dictionary'.format())

        # min_biomass_percent and product_targetYield_percent 
        if attr_name in ['min_biomass_percent','product_targetYield_percent'] and (not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError('{} must be an integer or a float'.format(attr_name))
        if attr_name in ['min_biomass_percent','product_targetYield_percent'] and (attr_value < 0 or attr_value > 100):
            raise ValueError('{} must be between 0 and 100'.format(attr_name))

        # growthMedium_flux_bounds 
        if attr_name == 'growthMedium_flux_bounds' and not isinstance(attr_value,dict):
            raise TypeError('growthMedium_flux_bounds must be a dictionary')
        if attr_name == 'growthMedium_flux_bounds' and len([k for k in attr_value.keys() if k.lower() not in ['flux_bounds_filename','flux_bounds_dict']]) > 0:
            raise ValueError('Invalid key for growthMedium_flux_bounds. Allowed keys are flux_bounds_filename and flux_bounds_dict')

        # MUST_X, MUST_L and MUST_U and blocked rxns 
        if attr_name in ['MUST_X','MUST_L','MUST_U','blocked_rxns'] and not isinstance(attr_value,list):
            raise TypeError('{} must be a list of strings'.format(attr_name))
        elif attr_name in ['MUST_X','MUST_L','MUST_U','blocked_rxns'] and len([k for k in attr_value if not isinstance(k,str)]) > 0:
            raise TypeError('{} must be a list of string. Non-string objects were observed in the list'.format(attr_name))

        if attr_name == 'MUST_double_type' and not isinstance(attr_value,str):
            raise TypeError('MUST_double_type must be a string')
        elif attr_name == 'MUST_double_type' and attr_value.lower() not in ['must_lu','must_ul','must_uu','must_ll']: 
            raise TypeError('Invalid MUST_double_type value! Allowed choices are MUST_LU, MUST_UL, MUST_LL, MUST_UU')

        if attr_name == 'build_new_optModel' and not isinstance(attr_value,bool):
            raise TypeError('build_new_optModel must be either True or False')

        if attr_name == 'objective_thr' and not isinstance(attr_value,float) and not isinstance(attr_value,int):
            raise TypeError('objective_thr must be either a float or an integer')

        if attr_name == 'validate_results' and not isinstance(attr_value,bool):
            raise TypeError('validate_results must be eiher True or False')

        if attr_name == 'results_filename' and not isinstance(attr_value,str):
            raise TypeError('results_filename must be a string')

        if attr_name == 'create_new_results_file' and not isinstance(attr_value,bool):
            raise TypeError('create_new_results_file must be eiher True or False')

        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi','gurobi_ampl','cplexamp']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError('warnings must be eiher True or False')

        if attr_name == 'stoud_msgs' and not isinstance(attr_value,bool):
            raise TypeError('stdout_msgs must be eiher True or False')

        self.__dict__[attr_name] = attr_value

    def find_maxBiomass_flux(self):
        """
        Finds the maximum thoeretical biomass flux
        """
        for rxn in self.model.reactions:
            rxn.objective_coefficient = 0
        self.model.biomass_reaction.objective_coefficient = 1

        set_specific_bounds(model = self.model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)

        self.model.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False, warnings = True)
        if self.model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            self.max_biomass_flux = self.model.fba_model.solution['objective_value']
            if self.stdout_msgs:
                print 'OptForce: Theoretical maximum biomass formation flux = {}'.format(self.max_biomass_flux)
        else:
            raise userError('The FBA problem for finding the maximum theoretical biomass flux was not solved to optimality!')

    def find_prod_max_theor_yield(self):
        """
        Finds the maximum thoeretical yield for the product of interest (this is the max flux of the exchange reactions for the product)
        """
        for rxn in self.model.reactions:
            rxn.objective_coefficient = 0
        self.model.reactions_by_id[self.product_exchrxn_id].objective_coefficient = 1

        set_specific_bounds(model = self.model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)

        self.model.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False, warnings = True)
        if self.model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            self.product_max_theor_yield = self.model.fba_model.solution['objective_value']
            if self.stdout_msgs:
                print 'OptForce: Theoretical maximum product yield = {}'.format(self.product_max_theor_yield)
        else:
            raise userError('The FBA problem for finding the theoretical maximum product yield was not solved to optimality!')

    #------------------------------------------------------------------------
    #--- Define parameters needed for the optimizaiton problem ---
    #------------------------------------------------------------------------
    def set_model_flux_bounds(self):
        """
        Set the flux bounds for the growth medium and min biomass flux in the model
        This functions needs to be caleed before creating an optModel
        """
        # Find the maximum biomass flux and  maximum theoretical yield of the product
        self.find_maxBiomass_flux()
        self.find_prod_max_theor_yield()

        # Set the flux bounds for the model
        set_specific_bounds(model = self.model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)
        # Impost constraints on product yield and min biomass formation flux
        set_specific_bounds(model = self.model, flux_bounds = {self.biomass_rxn_id:[(self.min_biomass_percent/100)*self.max_biomass_flux,None], self.product_exchrxn_id:[(self.product_targetYield_percent/100)*self.product_max_theor_yield,None]}, reset_flux_bounds = False)


    def define_optModel_params(self):
        """
        Assigns model parameters
        """
        self._muLB_max = 1e4
        self._muUB_max = 1e4

        # Reactions that must be excluded (ignore) from MUST doube consideration
        self._ignore = []
        # Rxns in MUST single sets 
        # blocked rxns under the given uptake and aeration conditions
        # Rxns whose flux bounds have a lower and upper bound of zero
        # Rxns whose flux bounds in the reference and overproducing strains are the same
        # Transport reactions
        # Exchange reactions
        # Reactions with no gene association
        self._ignore = self.blocked_rxns
        for rxn in self.model.reactions:
            if rxn.id in self.MUST_X + self.MUST_L + self.MUST_U or \
               (abs(rxn.flux_bounds[0]) < 1e-6 and abs(rxn.flux_bounds[1]) < 1e-6) or \
               (self.flux_bounds_ref[rxn.id] == self.flux_bounds_overprod[rxn.id]) or \
               'transport' in rxn.name.lower() or \
               rxn.reversibility.lower() == 'exchange' or \
               (rxn.genes == [] and rxn.gene_reaction_rule == ''):
                self._ignore.append(rxn.id)

        # biomass
        self._ignore.append(self.model.biomass_reaction.id)

        # ATPM
        if 'ATPM' in self.model.reactions_by_id.keys(): 
            self._ignore.append('ATPM')

    def assignFluxBounds(self,optModel,j):
        """
        Define the flux bounds
        """
        # v_biomass >= v_biomass_min
        if j == self.biomass_rxn_id:
            LB_UB = [(self.min_biomass_percent/100)*self.max_biomass_flux,self.model.reactions_by_id[j].flux_bounds[1]]

        # v_prodict >= v_product_target
        elif j == self.product_exchrxn_id: 
            LB_UB = [(self.product_targetYield_percent/100)*self.product_max_theor_yield, self.model.reactions_by_id[j].flux_bounds[1]]

        else: 
            LB_UB = self.model.reactions_by_id[j].flux_bounds

        return LB_UB 

    #------------------------------------------------------------------------
    #--- Define rules for the objective function and various constraints ----
    #------------------------------------------------------------------------
    #--- Constraints of the outer-level problem ---
    # Objective function
    def outer_objectiveFunc_rule(self,optModel):
        """
        Objective function of the outer-level problem 
        """
        # MUST_LU:  Max_OP(v1 - v2) < Min_WT(v1 - v2)
        # z = sum(j,[LBw(j)*yL(j) - UBw(j)*yU(j)] - [v_yL(j) - v_yU(j)])
        if self.MUST_double_type.lower() == 'must_lu':
            obj = sum([(self.flux_bounds_ref[j][0]*optModel.yL[j] - self.flux_bounds_ref[j][1]*optModel.yU[j]) - (optModel.v_yL[j] - optModel.v_yU[j]) for j in optModel.J])

        # MUST_UL: Min_OP(v1 - v2) > Max_WT(v1 - v2)
        # z = sum(j,[v_yU(j) - v_yL(j)] - [UBw(j)*yU(j) - LBw(j)*yL(j)])
        elif self.MUST_double_type.lower() == 'must_ul':
            obj = sum([(optModel.v_yU[j] - optModel.v_yL[j]) - (self.flux_bounds_ref[j][1]*optModel.yU[j] - self.flux_bounds_ref[j][0]*optModel.yL[j]) for j in optModel.J])

        # MUST_UU: Min_OP(v1 + v2) > Max_WT(v1 + v2)
        # z = sum(j,[v_y1U(j) + v_y2U(j)] - [UBw(j)*y1U(j) + UBw(j)*y2U(j)])
        elif self.MUST_double_type.lower() == 'must_uu':
            obj = sum([(optModel.v_y1U[j] + optModel.v_y2U[j]) - (self.flux_bounds_ref[j][1]*optModel.y1U[j] + self.flux_bounds_ref[j][1]*optModel.y2U[j]) for j in optModel.J])

        # MUST_LL: Max_OP(v1 + v2) < Min_WT(v1 + v2)
        # z = sum(j,[LBw(j)*y1L(j) + LBw(j)*y2L(j)] - [v_y1L(j) + v_y2L(j)]);
        elif self.MUST_double_type.lower() == 'must_ll':
            obj = sum([(self.flux_bounds_ref[j][0]*optModel.y1L[j] + self.flux_bounds_ref[j][0]*optModel.y2L[j]) - (optModel.v_y1L[j] + optModel.v_y2L[j]) for j in optModel.J])

        return obj

    def onlyOne_yL_isOne_const_rule(self, optModel):
        """
        sum(j,yL(j)) = 1
        """
        return sum([optModel.yL[j] for j in optModel.J]) == 1

    def onlyOne_yU_isOne_const_rule(self, optModel):
        """
        sum(j,yU(j)) = 1
        """
        return sum([optModel.yU[j] for j in optModel.J]) == 1

    def each_rxn_yLyU_const_rule(self, optModel, j):
        """
        yL[j] + yU[j] = 1 for each j 
        """
        return (optModel.yL[j] + optModel.yU[j]) <= 1

    def onlyOne_y1U_isOne_const_rule(self, optModel):
        """
        sum(j,y1U(j)) = 1
        """
        return sum([optModel.y1U[j] for j in optModel.J]) == 1

    def onlyOne_y2U_isOne_const_rule(self, optModel):
        """
        sum(j,y2U(j)) = 1
        """
        return sum([optModel.y2U[j] for j in optModel.J]) == 1

    def each_rxn_y1Uy2U_const_rule(self, optModel, j):
        """
        y1U[j] + y2U[j] = 1 for each j 
        """
        return (optModel.y1U[j] + optModel.y2U[j]) <= 1

    def onlyOne_y1L_isOne_const_rule(self, optModel):
        """
        sum(j,y1L(j)) = 1
        """
        return sum([optModel.y1L[j] for j in optModel.J]) == 1

    def onlyOne_y2L_isOne_const_rule(self, optModel):
        """
        sum(j,y2L(j)) = 1
        """
        return sum([optModel.y2L[j] for j in optModel.J]) == 1

    def each_rxn_y1Ly2L_const_rule(self, optModel, j):
        """
        y1L[j] + y2L[j] = 1 for each j 
        """
        return (optModel.y1L[j] + optModel.y2L[j]) <= 1

    #--- Constraints of the primal problem ---
    def primal_objectiveFunc_rule(self,optModel):
        """
        Objective function of the primal problem 
        """
        if self.MUST_double_type.lower() == 'must_lu':
            obj = sum([optModel.v[j]*optModel.yL[j] - optModel.v[j]*optModel.yU[j] for j in optModel.J])

        elif self.MUST_double_type.lower() == 'must_ul':
            obj = sum([optModel.v[j]*optModel.yU[j] - optModel.v[j]*optModel.yL[j] for j in optModel.J])

        elif self.MUST_double_type.lower() == 'must_uu':
            obj = sum([optModel.v[j]*optModel.y1U[j] + optModel.v[j]*optModel.y2U[j] for j in optModel.J])

        elif self.MUST_double_type.lower() == 'must_ll':
            obj = sum([optModel.v[j]*optModel.y1L[j] + optModel.v[j]*optModel.y2L[j] for j in optModel.J])

        return obj

    def massBalance_const_rule(self,optModel,i):
        """
        Mass balance 
        """
        return sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0

    #--- Constraints of the dual problem ---
    def dual_objectiveFunc_rule(self,optModel):
        """
        Objective function of the dual problem
        """
        if self.MUST_double_type.lower() in ['must_lu','must_ll']:
            obj = sum([optModel.muUB[j]*self.model.reactions_by_id[j].flux_bounds[1] for j in optModel.J]) - \
                  sum([optModel.muLB*self.model.reactions_by_id[j].flux_bounds[0] for j in optModel.J if j not in [self.biomass_rxn_id,self.product_exchrxn_id]]) - \
                  optModel.muLB[self.biomass_rxn_id]*(self.min_biomass_percent/100)*self.max_biomass_flux - \
                  optModel.muLB[self.product_exchrxn_id]*(self.product_targetYield_percent/100)*self.product_max_theor_yield 

        elif self.MUST_double_type.lower() in ['must_ul','must_uu']:
            obj = sum([optModel.muLB*self.model.reactions_by_id[j].flux_bounds[0] for j in optModel.J if j not in [self.biomass_rxn_id,self.product_exchrxn_id]]) + \
                  optModel.muLB[self.biomass_rxn_id]*(self.min_biomass_percent/100)*self.max_biomass_flux + \
                  optModel.muLB[self.product_exchrxn_id]*(self.product_targetYield_percent/100)*self.product_max_theor_yield - \
                  sum([optModel.muUB[j]*self.model.reactions_by_id[j].flux_bounds[1] for j in optModel.J]) 

        return obj

    def dual_const_rule(self,optModel,j):
        """
        Constraints of the dual problem
        """
        rxn = self.model.reactions_by_id[j]  # rxn object with id j

        # MUST_LU
        if self.MUST_double_type.lower() == 'must_lu':
            const = sum(rxn.stoichiometry[i]*optModel.Lambda[i.id] for i in rxn.compounds) + optModel.muUB[j] - optModel.muLB[j] == optModel.yL[j] - optModel.yU[j] 

        # MUST_UL
        elif self.MUST_double_type.lower() == 'must_ul':
            const = sum(rxn.stoichiometry[i]*optModel.Lambda[i.id] for i in rxn.compounds) + optModel.muLB[j] - optModel.muUB[j] == optModel.yU[j] - optModel.yL[j] 

        # MUST_LL
        elif self.MUST_double_type.lower() == 'must_ll':
            const = sum(rxn.stoichiometry[i]*optModel.Lambda[i.id] for i in rxn.compounds) + optModel.muUB[j] - optModel.muLB[j] == optModel.y1L[j] + optModel.y2L[j] 

        # MUST_UU
        elif self.MUST_double_type.lower() == 'must_uu':
            const = sum(rxn.stoichiometry[i]*optModel.Lambda[i.id] for i in rxn.compounds) + optModel.muLB[j] - optModel.muUB[j] == optModel.y1U[j] + optModel.y2U[j] 

        return const

    # --- Strong duality constraint ---
    def strong_duality_const_rule(self,optModel):
        """
        Strong duality constraint: primal objective = dual objective 
        """
        if self.MUST_double_type.lower() == 'must_lu':
            const = sum([optModel.v_yL[j] - optModel.v_yU[j] for j in optModel.J]) == \
                    sum([optModel.muUB[j]*self.model.reactions_by_id[j].flux_bounds[1] for j in optModel.J]) - \
                    sum([optModel.muLB[j]*self.model.reactions_by_id[j].flux_bounds[0] for j in optModel.J if j not in [self.biomass_rxn_id,self.product_exchrxn_id]]) - \
                    optModel.muLB[self.biomass_rxn_id]*(self.min_biomass_percent/100)*self.max_biomass_flux - \
                    optModel.muLB[self.product_exchrxn_id]*(self.product_targetYield_percent/100)*self.product_max_theor_yield 

        elif self.MUST_double_type.lower() == 'must_ul':
            const = sum([optModel.v_yU[j] - optModel.v_yL[j] for j in optModel.J]) == \
                    sum([optModel.muLB[j]*self.model.reactions_by_id[j].flux_bounds[0] for j in optModel.J if j not in [self.biomass_rxn_id,self.product_exchrxn_id]]) + \
                    optModel.muLB[self.biomass_rxn_id]*(self.min_biomass_percent/100)*self.max_biomass_flux + \
                    optModel.muLB[self.product_exchrxn_id]*(self.product_targetYield_percent/100)*self.product_max_theor_yield - \
                    sum([optModel.muUB[j]*self.model.reactions_by_id[j].flux_bounds[1] for j in optModel.J]) 


        elif self.MUST_double_type.lower() == 'must_uu':
            const = sum([optModel.v_y1U[j] + optModel.v_y2U[j]*optModel.y2U[j] for j in optModel.J]) == \
                    sum([optModel.muLB[j]*self.model.reactions_by_id[j].flux_bounds[0] for j in optModel.J if j not in [self.biomass_rxn_id,self.product_exchrxn_id]]) + \
                    optModel.muLB[self.biomass_rxn_id]*(self.min_biomass_percent/100)*self.max_biomass_flux + \
                    optModel.muLB[self.product_exchrxn_id]*(self.product_targetYield_percent/100)*self.product_max_theor_yield - \
                    sum([optModel.muUB[j]*self.model.reactions_by_id[j].flux_bounds[1] for j in optModel.J]) 


        elif self.MUST_double_type.lower() == 'must_ll':
            const = sum([optModel.v_y1L[j] + optModel.v_y2L[j] for j in optModel.J]) == \
                    sum([optModel.muUB[j]*self.model.reactions_by_id[j].flux_bounds[1] for j in optModel.J]) - \
                    sum([optModel.muLB[j]*self.model.reactions_by_id[j].flux_bounds[0] for j in optModel.J if j not in [self.biomass_rxn_id,self.product_exchrxn_id]]) - \
                    optModel.muLB[self.biomass_rxn_id]*(self.min_biomass_percent/100)*self.max_biomass_flux - \
                    optModel.muLB[self.product_exchrxn_id]*(self.product_targetYield_percent/100)*self.product_max_theor_yield 

        return const

    #--- Constraints linearizing the product of reactions fluxes and binary variables ---
    #-- v[j]*yL[j] --
    def linearize_v_yL_const1_rule(self,optModel,j):
        """
        LB[j]*yL[j] <= v_yL[j] <= UB[j]*yL[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return rxn.flux_bounds[0]*optModel.yL[j] <= optModel.v_yL[j]

    def linearize_v_yL_const2_rule(self,optModel,j):
        """
        LB[j]*yL[j] <= v_yL[j] <= UB[j]*yL[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_yL[j] <= rxn.flux_bounds[1]*optModel.yL[j]

    def linearize_v_yL_const3_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - yL[j])*UB[j] <= v_yL[j] <= v[j] - (1 - yL[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v[j] - (1 - optModel.yL[j])*rxn.flux_bounds[1] <= optModel.v_yL[j]

    def linearize_v_yL_const4_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - yL[j])*UB[j] <= v_yL[j] <= v[j] - (1 - yL[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_yL[j] <= optModel.v[j] - (1 - optModel.yL[j])*rxn.flux_bounds[0]

    #-- v[j]*yU[j] --
    def linearize_v_yU_const1_rule(self,optModel,j):
        """
        LB[j]*yU[j] <= v_yU[j] <= UB[j]*yU[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return rxn.flux_bounds[0]*optModel.yU[j] <= optModel.v_yU[j]

    def linearize_v_yU_const2_rule(self,optModel,j):
        """
        LB[j]*yU[j] <= v_yU[j] <= UB[j]*yU[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_yU[j] <= rxn.flux_bounds[1]*optModel.yU[j]

    def linearize_v_yU_const3_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - yU[j])*UB[j] <= v_yU[j] <= v[j] - (1 - yU[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v[j] - (1 - optModel.yU[j])*rxn.flux_bounds[1] <= optModel.v_yU[j]

    def linearize_v_yU_const4_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - yU[j])*UB[j] <= v_yU[j] <= v[j] - (1 - yU[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_yU[j] <= optModel.v[j] - (1 - optModel.yU[j])*rxn.flux_bounds[0]

    #-- v[j]*y1U[j] --
    def linearize_v_y1U_const1_rule(self,optModel,j):
        """
        LB[j]*y1U[j] <= v_y1U[j] <= UB[j]*y1U[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return rxn.flux_bounds[0]*optModel.y1U[j] <= optModel.v_y1U[j]

    def linearize_v_y1U_const2_rule(self,optModel,j):
        """
        LB[j]*y1U[j] <= v_y1U[j] <= UB[j]*y1U[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y1U[j] <= rxn.flux_bounds[1]*optModel.y1U[j]

    def linearize_v_y1U_const3_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y1U[j])*UB[j] <= v_y1U[j] <= v[j] - (1 - y1U[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v[j] - (1 - optModel.y1U[j])*rxn.flux_bounds[1] <= optModel.v_y1U[j]

    def linearize_v_y1U_const4_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y1U[j])*UB[j] <= v_y1U[j] <= v[j] - (1 - y1U[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y1U[j] <= optModel.v[j] - (1 - optModel.y1U[j])*rxn.flux_bounds[0]

    #-- v[j]*y2U[j] --
    def linearize_v_y2U_const1_rule(self,optModel,j):
        """
        LB[j]*y2U[j] <= v_y2U[j] <= UB[j]*y2U[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return rxn.flux_bounds[0]*optModel.y2U[j] <= optModel.v_y2U[j]

    def linearize_v_y2U_const2_rule(self,optModel,j):
        """
        LB[j]*y2U[j] <= v_y2U[j] <= UB[j]*y2U[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y2U[j] <= rxn.flux_bounds[1]*optModel.y2U[j]

    def linearize_v_y2U_const3_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y2U[j])*UB[j] <= v_y2U[j] <= v[j] - (1 - y2U[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v[j] - (1 - optModel.y2U[j])*rxn.flux_bounds[1] <= optModel.v_y2U[j]

    def linearize_v_y2U_const4_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y2U[j])*UB[j] <= v_y2U[j] <= v[j] - (1 - y2U[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y2U[j] <= optModel.v[j] - (1 - optModel.y2U[j])*rxn.flux_bounds[0]

    #-- v[j]*y1L[j] --
    def linearize_v_y1L_const1_rule(self,optModel,j):
        """
        LB[j]*y1L[j] <= v_y1L[j] <= UB[j]*y1L[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return rxn.flux_bounds[0]*optModel.y1L[j] <= optModel.v_y1L[j]

    def linearize_v_y1L_const2_rule(self,optModel,j):
        """
        LB[j]*y1L[j] <= v_y1L[j] <= UB[j]*y1L[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y1L[j] <= rxn.flux_bounds[1]*optModel.y1L[j]

    def linearize_v_y1L_const3_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y1L[j])*UB[j] <= v_y1L[j] <= v[j] - (1 - y1L[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v[j] - (1 - optModel.y1L[j])*rxn.flux_bounds[1] <= optModel.v_y1L[j]

    def linearize_v_y1L_const4_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y1L[j])*UB[j] <= v_y1L[j] <= v[j] - (1 - y1L[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y1L[j] <= optModel.v[j] - (1 - optModel.y1L[j])*rxn.flux_bounds[0]


    #-- v[j]*y2L[j] --
    def linearize_v_y2L_const1_rule(self,optModel,j):
        """
        LB[j]*y2L[j] <= v_y2L[j] <= UB[j]*y2L[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return rxn.flux_bounds[0]*optModel.y2L[j] <= optModel.v_y2L[j]

    def linearize_v_y2L_const2_rule(self,optModel,j):
        """
        LB[j]*y2L[j] <= v_y2L[j] <= UB[j]*y2L[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y2L[j] <= rxn.flux_bounds[1]*optModel.y2L[j]

    def linearize_v_y2L_const3_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y2L[j])*UB[j] <= v_y2L[j] <= v[j] - (1 - y2L[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v[j] - (1 - optModel.y2L[j])*rxn.flux_bounds[1] <= optModel.v_y2L[j]

    def linearize_v_y2L_const4_rule(self,optModel,j):
        """
        x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        v[j] - (1 - y2L[j])*UB[j] <= v_y2L[j] <= v[j] - (1 - y2L[j])*LB[j]

        The functions returns (lb, expr, ub) for a constraint in the form lb <= expor <= ub
        """
        rxn = self.model.reactions_by_id[j]  # rxn object
        return optModel.v_y2L[j] <= optModel.v[j] - (1 - optModel.y2L[j])*rxn.flux_bounds[0]

    #--------------------------------------------------------
    #---------------- Create optimization models ------------        
    #--------------------------------------------------------
    def build_primal_optModel(self):
        """
        Creates a pyomo optimization model for the primal problem 
        """
        # Define parameters and scalars needed to define the optimizaiton problem
        self.define_optModel_params()

        # Create a pyomo model optimization model
        optModel = ConcreteModel()

        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Variables --- 
        # Reaction fluxes
        optModel.v = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        # Binary variables in for reactions that must be down-regulated (yL) or up-regulated (yU) 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.yL = Var(optModel.J, domain=Boolean)
            optModel.yU = Var(optModel.J, domain=Boolean)
        elif self.MUST_double_type.lower() == 'must_ll':
            optModel.y1L = Var(optModel.J, domain=Boolean)
            optModel.y2L = Var(optModel.J, domain=Boolean)
        elif self.MUST_double_type.lower() == 'must_uu':
            optModel.y1U = Var(optModel.J, domain=Boolean)
            optModel.y2U = Var(optModel.J, domain=Boolean)

        #--- Objective function ---
        if self.MUST_double_type.lower() in ['must_lu','must_ll']:
            optModel.objectiveFunc = Objective(rule=self.primal_objectiveFunc_rule, sense = maximize)
        if self.MUST_double_type.lower() in ['must_ul','must_uu']:
            optModel.objectiveFunc = Objective(rule=self.primal_objectiveFunc_rule, sense = minimize)

        #--- Constraints ----
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule=self.massBalance_const_rule)

        self.optModel = optModel

    def build_dual_optModel(self):
        """
        Creates a pyomo optimization model for the dual problem 
        """
        # Define parameters and scalars needed to define the optimizaiton problem
        self.define_optModel_params()

        # Create a pyomo model optimization model
        optModel = ConcreteModel()

        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Variables --- 
        # Dual variables associated with steady-state mass balance constraints
        # It's better to not provide upper and lower bound on lambda otherwise the dual objective 
        # will be slightly different from the primal's
        optModel.Lambda = Var(optModel.I, domain=Reals)

        # Dual variables associated with v_j >= LB_j and v_j <= UB_j
        optModel.muLB = Var(optModel.J, domain=Reals, bounds = (0,self._muLB_max))
        optModel.muUB = Var(optModel.J, domain=Reals, bounds = (0,self._muUB_max))

        # Product of binary variables and reaction fluxes 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.v_yL = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
            optModel.v_yU = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
        elif self.MUST_double_type.lower() == 'must_ll':
            optModel.v_y1L = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
            optModel.v_y2L = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
        elif self.MUST_double_type.lower() == 'must_uu':
            optModel.v_y1U = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
            optModel.v_y2U = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        # Binary variables in for reactions that must be down-regulated (yL) or up-regulated (yU) 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.yL = Var(optModel.J, domain=Boolean)
            optModel.yU = Var(optModel.J, domain=Boolean)
        elif self.MUST_double_type.lower() == 'must_ll':
            optModel.y1L = Var(optModel.J, domain=Boolean)
            optModel.y2L = Var(optModel.J, domain=Boolean)
        elif self.MUST_double_type.lower() == 'must_uu':
            optModel.y1U = Var(optModel.J, domain=Boolean)
            optModel.y2U = Var(optModel.J, domain=Boolean)

        #--- Objective function ---
        if self.MUST_double_type.lower() in ['must_lu','must_ll']:
            optModel.objectiveFunc = Objective(rule=self.dual_objectiveFunc_rule, sense = minimize)
        if self.MUST_double_type.lower() in ['must_ul','must_uu']:
            optModel.objectiveFunc = Objective(rule=self.dual_objectiveFunc_rule, sense = maximize)

        #-- Constraints of the dual problem --
        # dual constraints 
        optModel.dual_const = Constraint(optModel.J, rule=self.dual_const_rule)

        # Constraints linearizing the product of v_j and binary variables 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.linearize_v_yL_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const1_rule)
            optModel.linearize_v_yL_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const2_rule)
            optModel.linearize_v_yL_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const3_rule)
            optModel.linearize_v_yL_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const4_rule)

            optModel.linearize_v_yU_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yU_const1_rule)
            optModel.linearize_v_yI_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yU_const2_rule)
            optModel.linearize_v_yU_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yU_const3_rule)
            optModel.linearize_v_yU_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yU_const4_rule)

        elif self.MUST_double_type.lower() in 'must_uu':
            optModel.linearize_v_y1U_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const1_rule)
            optModel.linearize_v_y1U_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const2_rule)
            optModel.linearize_v_y1U_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const3_rule)
            optModel.linearize_v_y1U_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const4_rule)

            optModel.linearize_v_y2U_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const1_rule)
            optModel.linearize_v_y2U_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const2_rule)
            optModel.linearize_v_y2U_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const3_rule)
            optModel.linearize_v_y2U_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const4_rule)

        elif self.MUST_double_type.lower() in 'must_ll':
            optModel.linearize_v_y1L_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const1_rule)
            optModel.linearize_v_y1L_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const2_rule)
            optModel.linearize_v_y1L_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const3_rule)
            optModel.linearize_v_y1L_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const4_rule)

            optModel.linearize_v_y2L_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const1_rule)
            optModel.linearize_v_y2L_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const2_rule)
            optModel.linearize_v_y2L_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const3_rule)
            optModel.linearize_v_y2L_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const4_rule)

        self.optModel = optModel

    def build_bilevel_optModel(self):
        """
        Creates a pyomo optimization model for the bilevel problem 
        """
        # Define parameters and scalars needed to define the optimizaiton problem
        self.define_optModel_params()

        # Create a pyomo model optimization model
        optModel = ConcreteModel()

        #------------------ Sets -----------------------
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #----------------- Variables -------------------
        #-- Variables of the primal problem --
        # Reaction fluxes
        optModel.v = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        # Binary variables in for reactions that must be down-regulated (yL) or up-regulated (yU) 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.yL = Var(optModel.J, domain=Boolean)
            optModel.yU = Var(optModel.J, domain=Boolean)
        elif self.MUST_double_type.lower() == 'must_ll':
            optModel.y1L = Var(optModel.J, domain=Boolean)
            optModel.y2L = Var(optModel.J, domain=Boolean)
        elif self.MUST_double_type.lower() == 'must_uu':
            optModel.y1U = Var(optModel.J, domain=Boolean)
            optModel.y2U = Var(optModel.J, domain=Boolean)

        #-- Variables of the dual problem --
        # Dual variables associated with steady-state mass balance constraints
        # It's better to not provide upper and lower bound on lambda otherwise the dual objective 
        # will be slightly different from the primal's
        optModel.Lambda = Var(optModel.I, domain=Reals)

        # Dual variables associated with v_j >= LB_j and v_j <= UB_j
        optModel.muLB = Var(optModel.J, domain=Reals, bounds = (0,self._muLB_max))
        optModel.muUB = Var(optModel.J, domain=Reals, bounds = (0,self._muUB_max))

        # Product of binary variables and reaction fluxes 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.v_yL = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
            optModel.v_yU = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
        elif self.MUST_double_type.lower() == 'must_ll':
            optModel.v_y1L = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
            optModel.v_y2L = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
        elif self.MUST_double_type.lower() == 'must_uu':
            optModel.v_y1U = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
            optModel.v_y2U = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        #------ Objective function of the outer problem ------
        optModel.objectiveFunc = Objective(rule = self.outer_objectiveFunc_rule, sense = maximize)

        #---------------- Constraints ----------------------
        #-- Constraints of the outer problem --
        # Only one yL or yU can be one at a time and each rxn can participate in only L or U set 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.onlyOne_yL_isOne_const = Constraint(rule=self.onlyOne_yL_isOne_const_rule)
            optModel.onlyOne_yU_isOne_const = Constraint(rule=self.onlyOne_yU_isOne_const_rule)
            optModel.each_rxn_yLyU_const = Constraint(optModel.J, rule=self.each_rxn_yLyU_const_rule)
        elif self.MUST_double_type.lower() == 'must_uu':
            optModel.onlyOne_y1U_isOne_const = Constraint(rule=self.onlyOne_y1U_isOne_const_rule)
            optModel.onlyOne_y2U_isOne_const = Constraint(rule=self.onlyOne_y2U_isOne_const_rule)
            optModel.each_rxn_y1Uy2U_const = Constraint(optModel.J, rule=self.each_rxn_y1Uy2U_const_rule)
        elif self.MUST_double_type.lower() == 'must_ll':
            optModel.onlyOne_y1L_isOne_const = Constraint(rule=self.onlyOne_y1L_isOne_const_rule)
            optModel.onlyOne_y2L_isOne_const = Constraint(rule=self.onlyOne_y2L_isOne_const_rule)
            optModel.each_rxn_y1Ly2L_const = Constraint(optModel.J, rule=self.each_rxn_y1Ly2L_const_rule)

        # Integer cuts
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.integer_cuts = ConstraintList(noruleinit=True)
        elif self.MUST_double_type.lower() in ['must_uu','must_ll']:
            optModel.integer_cuts1 = ConstraintList(noruleinit=True)
            optModel.integer_cuts2 = ConstraintList(noruleinit=True)
         
        #-- Constraints of the primal problem --
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule=self.massBalance_const_rule)

        # Constraints linearizing the product of v_j and binary variables 
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            optModel.linearize_v_vL_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const1_rule)
            optModel.linearize_v_vL_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const2_rule)
            optModel.linearize_v_vL_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const3_rule)
            optModel.linearize_v_vL_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const4_rule)

            optModel.linearize_v_yU_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yU_const1_rule)
            optModel.linearize_v_yU_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yU_const2_rule)
            optModel.linearize_v_yU_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const3_rule)
            optModel.linearize_v_yU_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_yL_const4_rule)

        elif self.MUST_double_type.lower() in 'must_uu':
            optModel.linearize_v_y1U_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const1_rule)
            optModel.linearize_v_y1U_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const2_rule)
            optModel.linearize_v_y1U_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const3_rule)
            optModel.linearize_v_y1U_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1U_const4_rule)

            optModel.linearize_v_y2U_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const1_rule)
            optModel.linearize_v_y2U_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const2_rule)
            optModel.linearize_v_y2U_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const3_rule)
            optModel.linearize_v_y2U_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2U_const4_rule)

        elif self.MUST_double_type.lower() in 'must_ll':
            optModel.linearize_v_y1L_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const1_rule)
            optModel.linearize_v_y1L_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const2_rule)
            optModel.linearize_v_y1L_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const3_rule)
            optModel.linearize_v_y1L_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y1L_const4_rule)

            optModel.linearize_v_y2L_const1 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const1_rule)
            optModel.linearize_v_y2L_const2 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const2_rule)
            optModel.linearize_v_y2L_const3 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const3_rule)
            optModel.linearize_v_y2L_const4 = Constraint([j for j in optModel.J if j not in self._ignore], rule=self.linearize_v_y2L_const4_rule)

        #-- Constraints of the dual problem --
        # dual constraints 
        optModel.dual_const = Constraint(optModel.J, rule=self.dual_const_rule)

        #-- Strong duality -- 
        optModel.strong_duality_const = Constraint(rule=self.strong_duality_const_rule)

        self.optModel = optModel

    def fix_known_variables(self):
        """
        Fixies all known binary variables to zero/one
        """
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            for rxn in self._ignore:
                self.optModel.yL[rxn] = 0
                self.optModel.yL[rxn].fixed = True
                self.optModel.v_yL[rxn] = 0
                self.optModel.v_yL[rxn].fixed = True

                self.optModel.yU[rxn] = 0
                self.optModel.yU[rxn].fixed = True
                self.optModel.v_yU[rxn] = 0
                self.optModel.v_yU[rxn].fixed = True

        elif self.MUST_double_type.lower() in 'must_uu':
            for rxn in self._ignore:
                self.optModel.y1U[rxn] = 0
                self.optModel.y1U[rxn].fixed = True
                self.optModel.v_y1U[rxn] = 0
                self.optModel.v_y1U[rxn].fixed = True

                self.optModel.y2U[rxn] = 0
                self.optModel.y2U[rxn].fixed = True
                self.optModel.v_y2U[rxn] = 0
                self.optModel.v_y2U[rxn].fixed = True

        elif self.MUST_double_type.lower() in 'must_ll':
            for rxn in self._ignore:
                self.optModel.y1L[rxn] = 0
                self.optModel.y1L[rxn].fixed = True
                self.optModel.v_y1L[rxn] = 0
                self.optModel.v_y1L[rxn].fixed = True

                self.optModel.y2L[rxn] = 0
                self.optModel.y2L[rxn].fixed = True
                self.optModel.v_y2L[rxn] = 0
                self.optModel.v_y2L[rxn].fixed = True

    def run_single(self):
        """
        Performs a single run of the code
        """
        # Processing and wall time for pyomo
        start_preproc_pyomo_pt = time.clock()
        start_preproc_pyomo_wt = time.time()

        # Instantiate the optModel
        self.optModel.preprocess()

        elapsed_preproc_pyomo_pt = str(timedelta(seconds = time.clock() - start_preproc_pyomo_pt))
        elapsed_preproc_pyomo_wt = str(timedelta(seconds = time.time() - start_preproc_pyomo_wt))

        #---- Solve the model ----
        #- Solve the optModel (tee=True shows the solver output) -
        try:
            # Processing and wall time for the solver
            start_solver_pt = time.clock()
            start_solver_wt = time.time()
            optSoln = self._optSolver.solve(self.optModel,tee=False)
            solver_flag = 'normal'

        # In the case of an error switch the solver
        except  Exception, e:
            if self.warnings:
                print '**WARNING (MUST_doubles)! {} failed with the following error: \n{} \nAn alternative solver is tried'.format(self.optimization_solver,e)

            if self.optimization_solver.lower() == 'gurobi':
                new_solver = 'cplex'
            elif self.optimization_solver.lower() == 'cplex':
                new_solver = 'gurobi'

            # Try solving with the alternative solver
            optSolver = pyomoSolverCreator(new_solver)
            try:
                start_solver_pt = time.clock()
                start_solver_wt = time.time()
                optSoln = optSolver.solve(self.optModel,tee=False)
                solver_flag = 'normal'
            except   Exception, e:
                solver_flag = 'solverError'
                if self.warnings:
                    print '**WARNING (break_cycle.py)! {} failed with the following error: \n{} \nAn alternative solver is tried'.format(new_solver,e)

        elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
        elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))

        #----- Print the results in the output ------
        if solver_flag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':

            exit_flag = 'globallyOptimal'

            # Load the results
            self.optModel.load(optSoln)

            # Optimal value of the objective function
            opt_objValue = self.optModel.objectiveFunc()

            # List of reactions for which yL or yU is one 
            # NOTE: Instead of checking if binary variables are exactly equal to one, we check whether the difference 
            # between a binary variable and one is less than a threshold. This is because sometimes the solver stops with 
            # a binary variable that is not eactly one or zero according to the specified integrality tolerence used
            if self.MUST_double_type.lower() in ['must_lu','must_ul']:
                one_yLopt_rxns = [j for j in self.optModel.J if abs(self.optModel.yL[j].value - 1) <= mip_integrality_tol] 
                one_yUopt_rxns = [j for j in self.optModel.J if abs(self.optModel.yU[j].value - 1) <= mip_integrality_tol] 

                # self._yLopt_curr and self._yUopt_curr hold the optimal value of the binary variables at the
                # current interation
                self._yLopt_curr = {}
                self._yUopt_curr = {}
                for j in self.optModel.J:
                    self._yLopt_curr[j] = 0
                    self._yUopt_curr[j] = 0
                for j in one_yLopt_rxns:
                    self._yLopt_curr[j] = 1
                for j in one_yUopt_rxns:
                    self._yUopt_curr[j] = 1

            elif self.MUST_double_type.lower() in 'must_uu':
                one_y1Uopt_rxns = [j for j in self.optModel.J if abs(self.optModel.y1U[j].value - 1) <= mip_integrality_tol] 
                one_y2Uopt_rxns = [j for j in self.optModel.J if abs(self.optModel.y2U[j].value - 1) <= mip_integrality_tol] 

                # self._y1Uopt_curr and self._y2Uopt_curr hold the optimal value of the binary variables at the
                # current interation
                self._y1Uopt_curr = {}
                self._y2Uopt_curr = {}
                for j in self.optModel.J:
                    self._y1Uopt_curr[j] = 0
                    self._y2Uopt_curr[j] = 0
                for j in one_y1Uopt_rxns:
                    self._y1Uopt_curr[j] = 1
                for j in one_y2Uopt_rxns:
                    self._y2Uopt_curr[j] = 1

            elif self.MUST_double_type.lower() in 'must_ll':
                one_y1Lopt_rxns = [j for j in self.optModel.J if abs(self.optModel.y1L[j].value - 1) <= mip_integrality_tol] 
                one_y2Lopt_rxns = [j for j in self.optModel.J if abs(self.optModel.y2L[j].value - 1) <= mip_integrality_tol] 

                # self._y1Lopt_curr and self._y2Lopt_curr hold the optimal value of the binary variables at the
                # current interation
                self._y1Lopt_curr = {}
                self._y2Lopt_curr = {}
                for j in self.optModel.J:
                    self._y1Lopt_curr[j] = 0
                    self._y2Lopt_curr[j] = 0
                for j in one_y1Lopt_rxns:
                    self._y1Lopt_curr[j] = 1
                for j in one_y2Lopt_rxns:
                    self._y2Lopt_curr[j] = 1

        # If there was a solver error or if an optimal solution was not returned 
        else:
            if solver_flag == 'solverError':
                exit_flag = solver_flag
            else:
                exit_flag = str(optSoln.solver.termination_condition)

            opt_objValue = None
            if self.MUST_double_type.lower() in ['must_lu','must_ul']:
                one_yLopt_rxns = []
                one_yUopt_rxns = []
            elif self.MUST_double_type.lower() in 'must_uu':
                one_y1Uopt_rxns = []
                one_y2Uopt_rxns = []
            elif self.MUST_double_type.lower() in 'must_ll':
                one_y1Lopt_rxns = []
                one_y2Lopt_rxns = []

        #---- Solution ----
        if self.MUST_double_type.lower() in ['must_lu','must_ul']:
            self._curr_soln = {'exit_flag':exit_flag,'objective_value':opt_objValue,'L':one_yLopt_rxns,'U':one_yUopt_rxns}

        elif self.MUST_double_type.lower() in 'must_uu':
            self._curr_soln = {'exit_flag':exit_flag,'objective_value':opt_objValue,'U1':one_y1Uopt_rxns,'U2':one_y2Uopt_rxns}

        elif self.MUST_double_type.lower() in 'must_ll':
            self._curr_soln = {'exit_flag':exit_flag,'objective_value':opt_objValue,'L1':one_y1Lopt_rxns,'L2':one_y2Lopt_rxns}

        # Print the results on the screen 
        if self.stdout_msgs:
            print '\nObjective value = {}, Optimality status = {}, Solution status = {}, Solver run status = {}'.format(opt_objValue, optSoln.solver.termination_condition, optSoln.Solution.status, solver_flag)
            print 'Took (hh:mm:ss) {}/{} of processing/walltime to create a pyomo model, {}/{} to  preprcoess the model and {}/{} to solve the model\n'.format(self._elapsed_create_pyomo_pt, self._elapsed_create_pyomo_wt, elapsed_preproc_pyomo_pt,elapsed_preproc_pyomo_wt, elapsed_solver_pt,elapsed_solver_wt)

    def validate_soln(self):
        """
        This function validates the obtained results
        """
        # MUST_LU:  Max_OP(v1 - v2) < Min_WT(v1 - v2)
        if self.MUST_double_type.lower() == 'must_lu':
            if len(self._curr_soln['L']) != 1:
                raise userError('Length of one_yLopt_rxns is not one! len(one_yLopt_rxns) = {}'.format(len(self._curr_soln['L'])))
            else:
                L_rxn = self._curr_soln['L'][0]                   

                if (L_rxn in self.MUST_X) or (L_rxn in self.MUST_L) or (L_rxn in self.MUST_U):
                    raise userError('L_rxn {} appears in MUST_L or MUST_U or MUST_X'.format(L_rxn)) 

            if len(self._curr_soln['U']) != 1: 
                raise userError('Length of one_yUopt_rxns is not one! len(one_yUopt_rxns) = {}'.format(len(self._curr_soln['U'])))
            else:
                U_rxn = self._curr_soln['U'][0]   

                if (U_rxn in self.MUST_X) or (U_rxn in self.MUST_L) or (U_rxn in self.MUST_U):
                    raise userError('U_rxn {} appears in MUST_L or MUST_U or MUST_X'.format(U_rxn)) 

            if (self.flux_bounds_ref[L_rxn][0] - self.flux_bounds_ref[U_rxn][1]) - (self.optModel.v_yL[L_rxn].value - self.optModel.v_yU[U_rxn].value) < self.objective_thr:
                raise userError('Min_WT(v1 - v2) - Max_OP(v1 - v2) = ({}) - ({}) = {} is not less than objective_thr = {} for L_rxn = {}, U_rxn = {}'.format(self.flux_bounds_ref[L_rxn][0] - self.flux_bounds_ref[U_rxn][1], self.optModel.v_yL[L_rxn].value - self.optModel.v_yU[U_rxn].value, (self.flux_bounds_ref[L_rxn][0] - self.flux_bounds_ref[U_rxn][1]) - (self.optModel.v_yL[L_rxn].value - self.optModel.v_yU[U_rxn].value), self.objective_thr, L_rxn, U_rxn)) 

        # MUST_UL: Min_OP(v1 - v2) > Max_WT(v1 - v2)
        elif self.MUST_double_type.lower() == 'must_ul':
            if len(self._curr_soln['L']) != 1:
                raise userError('Length of one_yLopt_rxns is not one! len(one_yLopt_rxns) = {}'.format(len(self._curr_soln['L'])))
            else: 
                L_rxn = self._curr_soln['L'][0]                   

            if len(self._curr_soln['U']) != 1: 
                raise userError('Length of one_yUopt_rxns is not one! len(one_yUopt_rxns) = {}'.format(len(self._curr_soln['U'])))
            else: 
                U_rxn = self._curr_soln['U'][0]   

            if (self.optModel.v_yU[U_rxn].value - self.optModel.v_yL[L_rxn].value) - (self.flux_bounds_ref[U_rxn][1] - self.flux_bounds_ref[L_rxn][0]) < self.objective_thr:
                raise userError('Min_OP(v1 - v2) - Max_WT(v1 - v2) = ({}) - ({}) = {} is not less than objective_thr = {} for U_rxn = {} , L_rxn = {}'.format(self.optModel.v_yU[U_rxn].value - self.optModel.v_yL[L_rxn].value, self.flux_bounds_ref[U_rxn][1] - self.flux_bounds_ref[L_rxn][0], (self.optModel.v_yU[U_rxn].value - self.optModel.v_yL[L_rxn].value) - (self.flux_bounds_ref[U_rxn][1] - self.flux_bounds_ref[L_rxn][0]), self.objective_thr, U_rxn, L_rxn))

        # MUST_UU: Min_OP(v1 + v2) > Max_WT(v1 + v2)
        elif self.MUST_double_type.lower() == 'must_uu':
            if len(self._curr_soln['U1']) != 1:
                raise userError('Length of one_y1Uopt_rxns is not one! len(one_y1Uopt_rxns) = {}'.format(len(self._curr_soln['U1'])))
            else:
                U1_rxn = self._curr_soln['U1'][0]                   

                if (U1_rxn in self.MUST_X) or (U1_rxn in self.MUST_L) or (U1_rxn in self.MUST_U):
                    raise userError('U1_rxn {} appears in MUST_L or MUST_U or MUST_X'.format(U1_rxn)) 

            if len(self._curr_soln['U2']) != 1:
                raise userError('Length of one_y2Uopt_rxns is not one! len(one_y2Uopt_rxns) = {}'.format(len(self._curr_soln['U1'])))
            else:
                U2_rxn = self._curr_soln['U2'][0]   

                if (U2_rxn in self.MUST_X) or (U2_rxn in self.MUST_L) or (U2_rxn in self.MUST_U):
                    raise userError('U2_rxn {} appears in MUST_L or MUST_U or MUST_X'.format(U2_rxn)) 

            if (self.optModel.v_y1U[U1_rxn].value + self.optModel.v_y2U[U2_rxn].value) - (self.flux_bounds_ref[U1_rxn][1] + self.flux_bounds_ref[U2_rxn][1]) < self.objective_thr:
                raise userError('Min_OP(v1 + v2) - Max_WT(v1 + v2) = ({}) - ({}) = {} is not less than objective_thr = {} for U1_rxn = {} , U2_rxn = {}'.format(self.optModel.v_y1U[U1_rxn].value + self.optModel.v_y2U[U2_rxn].value, self.flux_bounds_ref[U1_rxn][1] + self.flux_bounds_ref[U2_rxn][1], (self.optModel.v_y1U[U1_rxn].value + self.optModel.v_y2U[U2_rxn].value) -  (self.flux_bounds_ref[U1_rxn][1] + self.flux_bounds_ref[U2_rxn][1]), self.objective_thr, U1_rxn, U2_rxn))    

        # MUST_LL: Max_OP(v1 + v2) < Min_WT(v1 + v2)
        elif self.MUST_double_type.lower() == 'must_ll':
            if len(self._curr_soln['L1']) != 1:
                raise userError('Length of one_y1Lopt_rxns is not one! len(one_y1Lopt_rxns) = {}'.format(len(self._curr_soln['L1'])))
            else:  
                L1_rxn = self._curr_soln['L1'][0]                   

                if (L1_rxn in self.MUST_X) or (L1_rxn in self.MUST_L) or (L1_rxn in self.MUST_U):
                    raise userError('L1_rxn {} appears in MUST_L or MUST_U or MUST_X'.format(L1_rxn)) 

            if len(self._curr_soln['L2']) != 1:
                raise userError('Length of one_y2Lopt_rxns is not one! len(one_y2Lopt_rxns) = {}'.format(len(self._curr_soln['L2'])))
            else: 
                L2_rxn = self._curr_soln['L2'][0]   

                if (L2_rxn in self.MUST_X) or (L2_rxn in self.MUST_L) or (L2_rxn in self.MUST_U):
                    raise userError('L2_rxn {} appears in MUST_L or MUST_U or MUST_X'.format(L2_rxn)) 

            if (self.flux_bounds_ref[L1_rxn][0] + self.flux_bounds_ref[L1_rxn][0]) - (self.optModel.v_y1L[L1_rxn].value + self.optModel.v_y2U[L2_rxn].value) < self.objective_thr:
                raise userError('Min_WT(v1 + v2) - Max_OP(v1 + v2) = ({}) - ({}) = {} is not less than objective_thr = {} for L1_rxn = {} , L2_rxn = {}'.format(self.flux_bounds_ref[L1_rxn][0] + self.flux_bounds_ref[L2_rxn][0], self.optModel.v_y1L[L1_rxn].value + self.optModel.v_y2L[L2_rxn].value), (self.flux_bounds_ref[L1_rxn][0] + self.flux_bounds_ref[L2_rxn][0]) - (self.optModel.v_y1L[L1_rxn].value + self.optModel.v_y2L[L2_rxn]).value, self.objective_thr, L1_rxn, L2_rxn)     

    def run(self, optModel_name = 'bilevel', max_solution_num = 100):
        """ 
        This method runs FBA. 

        INPUTS:
        -------
               optModel_name: Name of the optimizaiton model to solve
            max_solution_num: Total number of solutions to return. The default is a very large numbr meaning that it will return
                              all possible solutions
        OUTPUT:
        -------
        solution: A list of dictionaries with the following keys:
                        exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                                   or what is stored in optSoln.solver.termination_condition
                  objective_value: Optimal objective funtion value
                 one_yLopt_rxns: List of reactions for which the optimal value of yL is one
                 one_yUopt_rxns: List of reactions for which the optimal value of yU is one
        """
        # Total processing and wall time required to create the pyomo model, solve it and store the results 
        start_total_pt = time.clock()
        start_total_wt = time.time()

        # Set the model flux bounds
        self.set_model_flux_bounds()

        self.solution = []

        # Number of solutions found so far
        found_solutions_num = 0

        #---- Creating the pyomo optModel ----
        if self.build_new_optModel:
            start_pyomo_pt = time.clock()
            start_pyomo_wt = time.time()
            # Create the pyomo model optModel only if self.build_new_optModel == 1        
            if optModel_name.lower() == 'bilevel':
                self.build_bilevel_optModel()
            elif optModel_name.lower() == 'primal':
                self.build_primal_optModel()
            elif optModel_name.lower() == 'dual':
                self.build_dual_optModel()
            else:
                raise ValueError("Invalid optModel_name value. Allowed choices are 'bilevel', 'primal' and 'dual'")
            # Time to create a pyomo model
            self._elapsed_create_pyomo_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
            self._elapsed_create_pyomo_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))
        else:
            self._elapsed_create_pyomo_pt = 0 
            self._elapsed_create_pyomo_wt = 0 

        # Fix known variables
        self.fix_known_variables()

        # Create a solver and set the options
        self._optSolver = pyomoSolverCreator(self.optimization_solver)

        # results file
        if self.results_filename != '':
            if self.create_new_results_file:
                file_opener = 'w'
            else:
                file_opener = 'a'
            with open(self.results_filename,file_opener) as f:
                f.write('\n' + self.MUST_double_type + '_details = [\n')

        # A list of dictionaries holding the optimal solution in different iterations
        self.solutions = []

        done = False
        while not done:
            # Solve the optimizaiton model 
            self.run_single()

            # Check whether the current run was successful
            if self._curr_soln['exit_flag'] == 'globallyOptimal':

                if self._curr_soln['objective_value'] < self.objective_thr:
                    done = True
                    self.termination_condition = 'Objective function ({}) less than objective_thr ({})'.format(self._curr_soln['objective_value'],self.objective_thr)
                else:
                    # Validate the solution
                    if self.validate_results:
                        self.validate_soln()

                    if self.MUST_double_type.lower() in ['must_lu','must_ul']:
                        self._curr_soln['U'] = self._curr_soln['U'][0]
                        self._curr_soln['L'] = self._curr_soln['L'][0]
                    elif self.MUST_double_type.lower() == 'must_uu':
                        self._curr_soln['U1'] = self._curr_soln['U1'][0]
                        self._curr_soln['U2'] = self._curr_soln['U2'][0]
                    elif self.MUST_double_type.lower() == 'must_ll':
                        self._curr_soln['L1'] = self._curr_soln['L1'][0]
                        self._curr_soln['L2'] = self._curr_soln['L2'][0]

                    self.solutions.append(self._curr_soln)

                    # Add integer cut
                    if self.MUST_double_type.lower() in ['must_lu','must_ul']:
                        self.optModel.integer_cuts.add(sum([self.optModel.yL[j] for j in self.optModel.J if self._yLopt_curr[j] == 1]) + sum([self.optModel.yU[j] for j in self.optModel.J if self._yUopt_curr[j] == 1]) <= 1)
                    elif self.MUST_double_type.lower() == 'must_uu':
                        # Here, (y1U(j) , y2U(j2) and (y1U(j2) , y2U(j1) are equivalent solutions
                        self.optModel.integer_cuts1.add(sum([self.optModel.y1U[j] for j in self.optModel.J if self._y1Uopt_curr[j] == 1]) + sum([self.optModel.y2U[j] for j in self.optModel.J if self._y2Uopt_curr[j] == 1]) <= 1)
                        self.optModel.integer_cuts2.add(sum([self.optModel.y2U[j] for j in self.optModel.J if self._y1Uopt_curr[j] == 1]) + sum([self.optModel.y1U[j] for j in self.optModel.J if self._y2Uopt_curr[j] == 1]) <= 1) 
                    elif self.MUST_double_type.lower() == 'must_ll':
                        # Here, (y1L(j) , y2L(j2) and (y1L(j2) , y2L(j1) are equivalent solutions
                        self.self.optModel.integer_cuts1.add(sum([self.optModel.y1L[j] for j in self.optModel.J if self._y1Lopt_curr[j] == 1]) + sum([self.optModel.y2U[j] for j in self.optModel.J if self._y2Uopt_curr[j] == 1]) <= 1)
                        self.self.optModel.integer_cuts2.add(const = sum([self.optModel.y2L[j] for j in self.optModel.J if self._y1Lopt_curr[j] == 1]) + sum([self.optModel.y1U[j] for j in self.optModel.J if self._y2Uopt_curr[j] == 1]) <= 1)

                    found_solutions_num += 1

                    if self.stdout_msgs:
                       print '{} solutions found'.format(found_solutions_num)

                    # Write results into the file
                    if self.results_filename != '':
                        with open(self.results_filename,'a') as f:
                            # MUST_LU:  Max_OP(v1 - v2) < Min_WT(v1 - v2)
                            if self.MUST_double_type.lower() in 'must_lu':
                                L_rxn = self._curr_soln['L'] 
                                U_rxn = self._curr_soln['U']   
                                f.write("{{'L':'{}', 'U':'{}', 'Min_WT(v1 - v2)':{}, 'Max_OP(v1 - v2)':{}}},\n".format(L_rxn,U_rxn, self.flux_bounds_ref[L_rxn][0] - self.flux_bounds_ref[U_rxn][1], self.optModel.v_yL[L_rxn].value - self.optModel.v_yU[U_rxn].value))     
            
                            # MUST_UL: Min_OP(v1 - v2) > Max_WT(v1 - v2)
                            elif self.MUST_double_type.lower() in 'must_ul':
                                L_rxn = self._curr_soln['L']                   
                                U_rxn = self._curr_soln['U']   
                                f.write("{{'U':'{}', 'L':'{}', 'Min_OP(v1 - v2)':{}, 'Max_WT(v1 - v2)':{}}},\n".format(U_rxn,L_rxn, self.optModel.v_yU[U_rxn].value - self.optModel.v_yL[L_rxn].value, self.flux_bounds_ref[U_rxn][1] - self.flux_bounds_ref[L_rxn][0]))      
 
                            # MUST_UU: Min_OP(v1 + v2) > Max_WT(v1 + v2)
                            elif self.MUST_double_type.lower() in 'must_uu':
                                U1_rxn = self._curr_soln['U1']
                                U2_rxn = self._curr_soln['U2'] 
                                f.write("{{'U1':'{}', 'U2':'{}', 'Min_OP(v1 + v2)':{}, 'Max_WT(v1 + v2)':{}}},\n".format(U1_rxn,U2_rxn, self.optModel.v_y1U[U1_rxn].value + self.optModel.v_y2U[U2_rxn].value, self.flux_bounds_ref[U1_rxn][1] + self.flux_bounds_ref[U2_rxn][1]))    
            
                            # MUST_LL: Max_OP(v1 + v2) < Min_WT(v1 + v2)
                            elif self.MUST_double_type.lower() in 'must_ll':
                                L1_rxn = self._curr_soln['L1'] 
                                L2_rxn = self._curr_soln['L2']  
                                f.write("{{'L1':'{}', 'L2':'{}', 'Min_WT(v1 - v2)':{}, 'Max_OP(v1 - v2)':{}}},\n".format(L1_rxn,L2_rxn, self.flux_bounds_ref[L_rxn][0] + self.flux_bounds_ref[U_rxn][0], self.optModel.v_y1L[L1_rxn].value + self.optModel.v_y2U[L2_rxn].value))     

                    if found_solutions_num >= max_solution_num:
                        done = True
                        self.termination_condition = 'max_solution_num of {} was reached'.format(max_solution_num)

            else:
                done = True 
                self.termination_condition = 'Optimization problem was not solved to optimality: {}:'.format(self._curr_soln['exit_flag'])

        # Time required to perform FBA
        elapsed_total_pt = str(timedelta(seconds = time.clock() - start_total_pt))
        elapsed_total_wt = str(timedelta(seconds = time.time() - start_total_wt))

        if self.stdout_msgs:
           print '\nFinding {}s took {}/{} (hh:mm:ss) processing/walltime\n'.format(self.MUST_double_type, elapsed_total_pt,elapsed_total_wt)

        if self.results_filename != '':
            with open(self.results_filename,'a') as f:
                f.write(']\n')

        return (self.solutions, self.termination_condition)

#-------------------------
if __name__ == '__main__':
    pass 
 
