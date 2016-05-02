from __future__ import division
import sys,os, time
sys.path.append('../../../')
from copy import deepcopy
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from MUST_doubles import MUST_doubles
from tools.userError import userError
from tools.globalVariables import *
from tools.pyomoSolverCreator import pyomoSolverCreator
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.create_model import create_model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from imp import load_source
from coopr.pyomo import *
from coopr.opt import *
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class FORCE(object):
    """
    This class defines the optimization problem to identify the FORCE set of 
    of reactions, 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: March 28, 2016
    """
    def __init__(self,model, product_exchrxn_id, flux_bounds_ref, flux_bounds_overprod, growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, min_biomass_percent = 10, stopWith_product_yield_percent = 80, MUST_X = [], MUST_L = [], MUST_U = [], MUST_LU_L = [], MUST_LU_U = [], MUST_LL_L1 = [], MUST_LL_L2 = [], MUST_UU_U1 = [], MUST_UU_U2 = [], total_interven_num = 10, notMUST_total_interven_num = 0, notXrxns_interven_num = 0, notLrxns_interven_num = 0, notUrxns_interven_num = 0, fixed_X_rxns = [], ignored_X_rxns = [], fixed_L_rxns = [], ignored_L_rxns = [], fixed_U_rxns = {}, ignored_U_rxns = [], inSilico_essential_rxns = [], inVivo_essential_rxns= [], blocked_rxns = [], build_new_optModel = True, dual_formulation_type = 'simplified', validate_results = True, results_filename = '', optimization_solver = default_optim_solver, warnings = True, stdout_msgs = True, stdout_msgs_details = False): 
        """
        All inputs are a subset of those for OptForce

        INPUTS:
        ------
                total_interven_num: The maximum total number of interventions
        notMUST_total_interven_num: The maximum allowed number of interventions out of the MUST sets
             notXrxns_interven_num: Allowed number of knockouts outside of X_rxns
             notUrxns_interven_num: Allowed number of down_regulations outside of U_rxns
             notLrxns_interven_num: Allowed number of down_regulations outside of L_rxns
                build_new_optModel: If True ia new optimization model is created
          stopWith_product_yield_percent: Product target yield percent. Iterations stop if the objective function reaches this target yield
             dual_formulation_type: A string showing which type of dual formulation to use. 
                                      'standard': We deal with the product of binary variables and dual varibales, whcih ar elinearized
                                    'simplified': We analyze each term in the dual objective function considering the duality theory and
                                                  try to simplify the product of binary and dual variables. No linearization of the 
                                                  binary and dual variables is needed here.
        fixed_X_rxns, fixed_L_rxns, fixed_U_rxns: A list of reactions ids that must be knocked out, down-regulated or 
                                   up-regulated by fixing their corresponidng yX, yL or yU to one
        ignored_X_rxns, ignored_L_rxns, ignored_U_rxns: A list of reactions ids that must NOT be knocked out, down-regulate
                                   or up-regulated by fixing their corresponindg yX, yL or yU to zero
          inSilico_essential_rxns: List of in silico essential reaciton ids
            inVivo_essential_rxns: List of in vivo essential reaction ids
                     blocked_rxns: List of alwqys blocked reactions in the model
         growthMedium_flux_bounds: Information about growth media
                                   flux_bounds_filename: Name of and path of the file containing the flux bounds
                                                         file for exchange reactions corresponding to compounds in 
                                                         in the growth media
                                       flux_bounds_dict: A dictionary containing the flux bounds for other reactions (such as carbon
                                                         source uptake, oxygen uptake, etc), that have to be set for reactions not 
                                                         in flux_data_filenam or media_filename
                 validate_results: If True each obtained solution is validated before finding the next
                                   one
                 results_filename: A string containing the name of the file where the results should
                                   be written in

        Ali R. Zomorrodi - Segre lab @ BU
        Last updated: 04/01/2016
        """
        # model 
        self.model = model
        self.biomass_reaction = model.biomass_reaction
        self.biomass_rxn_id = model.biomass_reaction.id
        
        # Growth medium
        self.growthMedium_flux_bounds = growthMedium_flux_bounds

        # flux bounds
        self.flux_bounds_ref = flux_bounds_ref
        self.flux_bounds_overprod = flux_bounds_overprod

        # Product and biomass yields
        self.product_exchrxn_id = product_exchrxn_id
        self.stopWith_product_yield_percent = stopWith_product_yield_percent
        self.min_biomass_percent = min_biomass_percent

        # MUST sets
        self.MUST_X = MUST_X
        self.MUST_L = MUST_L
        self.MUST_U = MUST_U
        self.MUST_LU_L = MUST_LU_L
        self.MUST_LU_U = MUST_LU_U
        self.MUST_LL_L1 = MUST_LL_L1
        self.MUST_LL_L2 = MUST_LL_L2
        self.MUST_UU_U1 = MUST_UU_U1
        self.MUST_UU_U2 = MUST_UU_U2
        
        # Total number of interventions
        self.total_interven_num = total_interven_num

        # Allowed number of intervnetions out of the MUST sets
        self.notMUST_total_interven_num = notMUST_total_interven_num

        # Allowed number of interventions outside of X_rxns, L_rxns and U_rxns
        self.notXrxns_interven_num = notXrxns_interven_num
        self.notLrxns_interven_num = notLrxns_interven_num
        self.notUrxns_interven_num = notUrxns_interven_num

        # build_new_optModel
        self.build_new_optModel = build_new_optModel

        # List of in silico essential, in vivo essential and always blocked reactions
        self.inSilico_essential_rxns = inSilico_essential_rxns
        self.inVivo_essential_rxns = inVivo_essential_rxns
        self.blocked_rxns = blocked_rxns        

        # Type of the dual formulation
        self.dual_formulation_type = dual_formulation_type

        # Product target yield percent
        self.stopWith_product_yield_percent = stopWith_product_yield_percent

        # validate results
        self.validate_results = validate_results

        # results file name
        self.results_filename = results_filename

        # optimization solver
        self.optimization_solver = optimization_solver

        # warnings
        self.warnings = warnings

        # stdout_msgs
        self.stdout_msgs = stdout_msgs
        self.stdout_msgs_details = stdout_msgs_details

        # - Perform some preprocessing --
        # Reactions that must be up- or down-regulated according to MUST sets
        self.U_rxns = list(set(self.MUST_U + self.MUST_LU_U + self.MUST_UU_U1 + self.MUST_UU_U2))
        self.L_rxns = list(set(self.MUST_L + self.MUST_LU_L + self.MUST_LL_L1 + self.MUST_LL_L2))
        self.X_rxns = self.MUST_X 

        # Reactions that must be knocked out, down-regulated or up-regulated by fixing their 
        # corresponding binary variables yX, yL and yU to one
        self.fixed_X_rxns = fixed_X_rxns
        self.fixed_L_rxns = fixed_L_rxns
        self.fixed_U_rxns = fixed_U_rxns
        self.ignored_X_rxns = ignored_X_rxns
        self.ignored_L_rxns = ignored_L_rxns
        self.ignored_U_rxns = ignored_U_rxns

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

        # growthMedium_flux_bounds 
        if attr_name == 'growthMedium_flux_bounds' and not isinstance(attr_value,dict):
            raise TypeError('growthMedium_flux_bounds must be a dictionary')
        if attr_name == 'growthMedium_flux_bounds' and len([k for k in attr_value.keys() if k.lower() not in ['flux_bounds_filename','flux_bounds_dict']]) > 0:
            raise ValueError('Invalid key for growthMedium_flux_bounds. Allowed keys are flux_bounds_filename and flux_bounds_dict')

        # product_exchrxn_id 
        if attr_name == 'product_exchrxn_id' and not isinstance(attr_value,str):
            raise TypeError('product_exchrxn_id must be a string')

        # min_biomass_percent and stopWith_product_yield_percent 
        if attr_name in ['min_biomass_percent','stopWith_product_yield_percent'] and (not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError('{} must be an integer or a float'.format(attr_name))
        if attr_name in ['min_biomass_percent','stopWith_product_yield_percent'] and (attr_value < 0 or attr_value > 100):
            raise ValueError('{} must be between 0 and 100'.format(attr_name))

        # MUST_singles, MUST_doubles, essential and blocked rxns 
        if attr_name in ['MUST_X','MUST_L','MUST_U','MUST_LU_L', 'MUST_LU_U', 'MUST_LL_L1', 'MUST_LL_L2', 'MUST_UU_U1', 'MUST_UU_U2', 'blocked_rxns', 'inSilico_essential_rxns','inVivo_essential_rxns'] and not isinstance(attr_value,list):
            raise TypeError('{} must be a list of strings'.format(attr_name))
        elif attr_name in ['MUST_X','MUST_L','MUST_U','MUST_LU_L', 'MUST_LU_U', 'MUST_LL_L1', 'MUST_LL_L2', 'MUST_UU_U1', 'MUST_UU_U2', 'blocked_rxns', 'inSilico_essential_rxns','inVivo_essential_rxns'] and len([k for k in attr_value if not isinstance(k,str)]) > 0:
            raise TypeError('{} must be a list of string. Non-string objects were observed in the list'.format(attr_name))

        if attr_name == 'total_interven_num' and not isinstance(attr_value,int):
            raise TypeError('total_interven_num must be an integer')

        if attr_name == 'notMUST_total_interven_num' and not isinstance(attr_value,int):
            raise TypeError('notMUST_total_interven_num must be an integer')
        elif attr_name == 'notMUST_total_interven_num' and attr_value < 0:
            raise ValueError('notMUST_total_interven_num must be a non-negative integer') 

        if attr_name == 'notXrxns_interven_num' and not isinstance(attr_value,int):
            raise TypeError('notXrxns_interven_num must be an integer')
        elif attr_name == 'notXrxns_interven_num' and attr_value < 0:
            raise ValueError('notXrxns_interven_num must be a non-negative integer') 

        if attr_name == 'notLrxns_interven_num' and not isinstance(attr_value,int):
            raise TypeError('notLrxns_interven_num must be an integer')
        elif attr_name == 'notLrxns_interven_num' and attr_value < 0:
            raise ValueError('notLrxns_interven_num must be a non-negative integer') 

        if attr_name == 'notUrxns_interven_num' and not isinstance(attr_value,int):
            raise TypeError('notUrxns_interven_num must be an integer')
        elif attr_name == 'notUrxns_interven_num' and attr_value < 0:
            raise ValueError('notUrxns_interven_num must be a non-negative integer') 

        if attr_name in ['fixed_X_rxns','fixed_L_rxns','fixed_U_rxns','ignore_X_rxns','ignored_L_rxns','ignored_U_rxns'] and not isinstance(attr_value,list):
            raise TypeError('{} must be a list of strings'.format(attr_name))
        elif attr_name in ['fixed_X_rxns','fixed_L_rxns','fixed_U_rxns'] and len([k for k in attr_value if not isinstance(k,str)]) > 0: 
            raise TypeError('Elements of {} must be a string'.format(attr_name))

        if attr_name in ['inSilico_essential_rxns','inVivo_essential_rxns','blocked_rxns'] and not isinstance(attr_value,list):
            raise TypeError('{} must be a list of strings'.fomrat(attr_name))
        elif attr_name in ['inSilico_essential_rxns','inVivo_essential_rxns','blocked_rxns'] and len([k for k in attr_value if not isinstance(k,str)]) > 0: 
            raise TypeError('Element of {} must be a string'.format(attr_name))

        if attr_name == 'build_new_optModel' and not isinstance(attr_value,bool):
            raise TypeError('build_new_optModel must be either True or False')

        if attr_name == 'dual_formulation_type' and not isinstance(attr_value,str):
            raise TypeError('dual_formulation_type must be a string')
        elif attr_name == 'dual_formulation_type' and attr_value.lower() not in ['standard','simplified']: 
            raise TypeError('Invalid dual_formulation_type value! Allowed choices are standard and simplified')

        if attr_name == 'stopWith_product_yield_percent' and not isinstance(attr_value,float) and not isinstance(attr_value,int):
            raise TypeError('stopWith_product_yield_percent must be an integer or float')
        elif attr_name == 'stopWith_product_yield_percent' and (attr_value < 0 or attr_value > 100):
            raise TypeError('stopWith_product_yield_percent must be an integer or float between 0 and 100')

        if attr_name == 'validate_results' and not isinstance(attr_value,bool):
            raise TypeError('validate_results must be eiher True or False')

        if attr_name == 'results_filename' and not isinstance(attr_value,str):
            raise TypeError('results_filename must be a string')

        if attr_name in ['warnings','stoud_msgs','stdout_msgs_details'] and not isinstance(attr_value,bool):
            raise TypeError('{} must be eiher True or False'.format(attr_name))

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
        if not hasattr(self,'max_biomass_flux'):
            self.find_maxBiomass_flux()
        if not hasattr(self,'product_max_theor_yield'):
            self.find_prod_max_theor_yield()

        # Set the flux bounds for the model
        set_specific_bounds(model = self.model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)
        # Impost constraints on min biomass formation flux
        set_specific_bounds(model = self.model, flux_bounds = {self.biomass_rxn_id:[(self.min_biomass_percent/100)*self.max_biomass_flux,None]}, reset_flux_bounds = False)

    def define_optModel_params(self):
        """
        Assigns model parameters
        """
        # big M
        self._bigM_rxnflux = 1000
  
        # Max value of the dual variables
        self._bigM_dualvar = 1e4

        # Reactions that must be excluded (ignore) from any type of manupulaiton (X, L or U) 
        ignore_all = []
        for rxn in self.model.reactions:
            if rxn in self.blocked_rxns or \
               rxn.flux_bounds == [0,0] or \
               'transport' in rxn.name.lower() or \
               rxn.subsystem != '' and 'transport' in rxn.subsystem.lower() or \
               rxn.reversibility.lower() == 'exchange' or \
               (rxn.genes == [] and rxn.gene_reaction_rule == ''):
                ignore_all.append(rxn.id)

        # biomass
        ignore_all.append(self.model.biomass_reaction.id)

        # ATPM
        if 'ATPM' in self.model.reactions_by_id.keys():
            ignore_all.append('ATPM')

        #--- Reactions that must be ignored for knock outs --- 
        self.ignored_X_rxns += ignore_all

        # in silico and in vivo essential rxns
        self.ignored_X_rxns += list(set(self.inSilico_essential_rxns + self.inVivo_essential_rxns))

        # Reactions not in X_rxns 
        if self.notXrxns_interven_num == 0:
            self.ignored_X_rxns += [r.id for r in self.model.reactions if r.id not in self.X_rxns]

        #---- L_rxns ----
        self.ignored_L_rxns += ignore_all + self.U_rxns

        # Reactions not in L_rxns 
        if self.notLrxns_interven_num == 0:
            self.ignored_L_rxns +=  [r.id for r in self.model.reactions if r.id not in self.L_rxns]

        #---- U_rxns ----
        self.ignored_U_rxns += ignore_all + self.L_rxns

        # Reactions not in U_rxns 
        if self.notUrxns_interven_num == 0:
            self.ignored_U_rxns += [r.id for r in self.model.reactions if r.id not in self.U_rxns]

    #------------------------------------------------------------------------
    #--- Define rules for the objective function and various constraints ----
    #------------------------------------------------------------------------
    def strong_duality_const_rule(self,optModel):
        """
        Objective function of the dual problem
        """
        if self.dual_formulation_type.lower() == 'standard':
            const = \
                  optModel.v[self.product_exchrxn_id] ==  \
                  sum([self.model.reactions_by_id[j].flux_bounds[0]*(optModel.muLB[j] - optModel.muLB_yX[j]) for j in optModel.J]) - \
                  sum([self.model.reactions_by_id[j].flux_bounds[1]*(optModel.muUB[j] - optModel.muUB_yX[j]) for j in optModel.J]) + \
                  sum([(self.flux_bounds_overprod[j][0] + self._bigM_rxnflux)*optModel.thethaLB_yU[j] - self._bigM_rxnflux*optModel.thethaLB[j] for j in optModel.J if j not in self.ignored_U_rxns]) - \
                  sum([(self.flux_bounds_overprod[j][1] - self._bigM_rxnflux)*optModel.thethaUB_yL[j] + self._bigM_rxnflux*optModel.thethaUB[j] for j in optModel.J if j not in self.ignored_L_rxns])
                        
        elif self.dual_formulation_type.lower() == 'simplified':
            const = \
                  optModel.v[self.product_exchrxn_id] == \
                  sum([self.model.reactions_by_id[j].flux_bounds[0]*optModel.muLB[j] for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[0] != -self._bigM_rxnflux]) - \
                  sum([self.model.reactions_by_id[j].flux_bounds[1]*optModel.muUB[j] for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[1] != self._bigM_rxnflux]) + \
                  sum([self.flux_bounds_overprod[j][0]*optModel.thethaLB[j] for j in optModel.J if j not in self.ignored_U_rxns]) - \
                  sum([self.flux_bounds_overprod[j][1]*optModel.thethaUB[j] for j in optModel.J if j not in self.ignored_L_rxns])

        return const


    #--- Constraints of the dual problem ---
    def dual_objectiveFunc_rule(self,optModel):
        """
        Objective function of the dual problem
        """
        if self.dual_formulation_type.lower() == 'standard':
            # sum(j, muLB(j)*[LB(j)*(1 - yX(j)] + 
            #        muUB(j)*[UB(j)*(1 - yX(j)] +]
            #        thethaLB(j)*[(LBon(j) - (-M))*yU(j) + (-M)] +
            #        thethaUB(j)*[(UBon(j) - bigM)*yU(j) + bigM]
            obj = sum([self.model.reactions_by_id[j].flux_bounds[0]*(optModel.muLB[j] - optModel.muLB_yX[j]) for j in optModel.J if j not in self.ignored_X_rxns]) - \
                  sum([self.model.reactions_by_id[j].flux_bounds[1]*(optModel.muUB[j] - optModel.muUB_yX[j]) for j in optModel.J if j not in self.ignored_X_rxns]) + \
                  sum([(self.flux_bounds_overprod[j][0] + self._bigM_rxnflux)*optModel.thethaLB_yU[j] - self._bigM_rxnflux*optModel.thethaLB[j] for j in optModel.J if j not in self.ignored_U_rxns]) - \
                  sum([(self.flux_bounds_overprod[j][1] - self._bigM_rxnflux)*optModel.thethaUB_yL[j] + self._bigM_rxnflux*optModel.thethaUB[j] for j in optModel.J if j not in self.ignored_L_rxns])
                        
        elif self.dual_formulation_type.lower() == 'simplified':
            obj = sum([self.model.reactions_by_id[j].flux_bounds[0]*optModel.muLB[j] for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[0] != -self._bigM_rxnflux]) - \
                  sum([self.model.reactions_by_id[j].flux_bounds[1]*optModel.muUB[j] for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[1] != self._bigM_rxnflux]) + \
                  sum([self.flux_bounds_overprod[j][0]*optModel.thethaLB[j] for j in optModel.J if j not in self.ignored_U_rxns]) - \
                  sum([self.flux_bounds_overprod[j][1]*optModel.thethaUB[j] for j in optModel.J if j not in self.ignored_L_rxns])

        return obj

    def dual_const_rule(self,optModel,j):
        """
        Constraints of the dual problem
        """
        rxn = self.model.reactions_by_id[j]  # rxn object with id j
        if j == self.product_exchrxn_id:
            const = sum(rxn.stoichiometry[i]*optModel.Lambda[i.id] for i in rxn.compounds) + optModel.muLB[j] - optModel.muUB[j] + optModel.thethaLB[j] - optModel.thethaUB[j] == 1 
        else:
            const = sum(rxn.stoichiometry[i]*optModel.Lambda[i.id] for i in rxn.compounds) + optModel.muLB[j] - optModel.muUB[j] + optModel.thethaLB[j] - optModel.thethaUB[j] == 0
        return const

    #--------------------------------------------------------
    #---------------- Create optimization models ------------        
    #--------------------------------------------------------
    def build_primal_optModel(self):
        """
        Creates a pyomo optimization model for the primal problem 
        """
        # Create a pyomo model optimization model
        optModel = ConcreteModel()

        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Variables --- 
        # Reaction fluxes (NOTE: Make sure flux bounds are correctly assigned)
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds)

        # Binary variables in for reactions considered for up-regulation, down-regulation or knockout 
        optModel.yL = Var(optModel.J, domain=Boolean)
        optModel.yU = Var(optModel.J, domain=Boolean)
        optModel.yX = Var(optModel.J, domain=Boolean)

        #--- Objective function ---
        optModel.objectiveFunc = Objective(rule = lambda optModel: optModel.v[self.product_exchrxn_id], sense = minimize)

        #--- Constraints ----
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)

        # Constrain for U_rxns
        # v(j) =g= LBon(j)*yU(j) + (1-yU(j))*LB(j) or v(j) =g= (LBon(j) - LB(j))*yU(j) + LB(j).
        # Instead of LB(j) use -bigM to make the constraint inactive when yU(j) = 0:
        # v(j) =g= (LBon(j) - (-M))*yU(j) + (-M) 
        optModel.U_rxns_const = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.v[j] >= (self.flux_bounds_overprod[j][0] + self._bigM_rxnflux)*optModel.yU[j] - self._bigM_rxnflux)

        # Constraint for U_rxns
        # v(j) =l= UBon(j)*yL(j) + (1-yL(j))*UB(j) or v(j) =l= (UBon(j) - UB(j))*yL(j) + UB(j)
        # Instead of UB(j) use bigM to make the constraint inactive when yL(j) = 0:
        # v(j) =l= (UBon(j) - bigM)*yU(j) + bigM
        optModel.L_rxns_const = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optModel, j: optModel.v[j] <= (self.flux_bounds_overprod[j][1] - self._bigM_rxnflux)*optModel.yL[j] + self._bigM_rxnflux)

        # Cconstraints for X_rxns. This constraint is written over all reactions (not just thos ein X_rxns) to incorporate
        # their LB and UB into the calculations
        # LB(j)*(1 - yX(j)) <= v(j)  and v(j) <= UB(j)*(1 - yX(j)) 
        optModel.X_rxns_const1 = Constraint(optModel.J, rule = lambda optModel, j: optModel.v[j] >= self.model.reactions_by_id[j].flux_bounds[0]*(1 - optModel.yX[j]))
        optModel.X_rxns_const2 = Constraint(optModel.J, rule = lambda optModel, j: optModel.v[j] <= self.model.reactions_by_id[j].flux_bounds[1]*(1 - optModel.yX[j]))

        self.optModel = optModel

    def build_dual_optModel(self):
        """
        Creates a pyomo optimization model for the dual problem 
        """
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

        # Dual variables associated with v_j >= LB_j*(1-yX(j) and v_j <= UB_j*(1-yX(j))
        optModel.muLB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        optModel.muUB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        if self.dual_formulation_type.lower() == 'standard':
            optModel.muLB_yX = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
            optModel.muUB_yX = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))

        # Dual variables associated with U_rxns_const: v(j) >= (LBon(j) - (-M))*yU(j) + (-M) (thethaLB) and  
        # L_rxns_const: v(j) <= (UBon(j) - bigM)*yU(j) + bigM (thethaUB)
        optModel.thethaLB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        optModel.thethaUB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        if self.dual_formulation_type.lower() == 'standard':
            optModel.thethaLB_yU = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
            optModel.thethaUB_yL = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))

        # Binary variables
        optModel.yL = Var(optModel.J, domain=Boolean)
        optModel.yU = Var(optModel.J, domain=Boolean)
        optModel.yX = Var(optModel.J, domain=Boolean)

        #--- Objective function ---
        optModel.objectiveFunc = Objective(rule=self.dual_objectiveFunc_rule, sense = maximize)

        #-- Constraints of the dual problem --
        # dual constraints 
        optModel.dual_const = Constraint(optModel.J, rule=self.dual_const_rule)

        #- Constraints linearizing the product of v_j and binary variables -
        # x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        # muLB[j] - (1 - yX[j])*muLB_max[j] <= muLB_yX[j] <= muLB[j] - (1 - yX[j])*muLB_min[j] and muLB_min = 0
        if self.dual_formulation_type.lower() == 'standard':
            optModel.linearize_muLB_yX_const1 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muLB_yX[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.linearize_muLB_yX_const2 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muLB[j] - (1 - optModel.yX[j])*self._bigM_dualvar <= optModel.muLB_yX[j])
            optModel.linearize_muLB_yX_const3 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muLB_yX[j] <= optModel.muLB[j] - (1 - optModel.yX[j])*0)

            # muUB_min[j]*yX[j] <= muUB_yX[j] <= muUB_max*yX[j]  and muUB_min = 0
            optModel.linearize_muUB_yX_const1 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muUB_yX[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.linearize_muUB_yX_const2 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muUB[j] - (1 - optModel.yX[j])*self._bigM_dualvar <= optModel.muUB_yX[j])
            optModel.linearize_muUB_yX_const3 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muUB_yX[j] <= optModel.muUB[j] - (1 - optModel.yX[j])*0)

            # thethaLB_min[j]*yU[j] <= thethaLB_yU[j] <= thethaLB_max*yU[j]  and thethaLB_min = 0
            optModel.linearize_thethaLB_yU_const1 = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.thethaLB_yU[j] <= self._bigM_dualvar*optModel.yU[j])
            optModel.linearize_thethaLB_yU_const2 = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.thethaLB[j] - (1 - optModel.yU[j])*self._bigM_dualvar <= optModel.thethaLB_yU[j])
            optModel.linearize_thethaLB_yU_const3 = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.thethaLB_yU[j] <= optModel.thethaLB[j] - (1 - optModel.yU[j])*0)

            # thethaUB_min[j]*yL[j] <= thethaUB_yL[j] <= thethaUB_max*yL[j]  and thethaUB_min = 0
            optModel.linearize_thethaUB_yL_const1 = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optModel, j: optModel.thethaUB_yL[j] <= self._bigM_dualvar*optModel.yL[j])
            optModel.linearize_thethaUB_yL_const2 = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optModel, j: optModel.thethaUB[j] - (1 - optModel.yL[j])*self._bigM_dualvar <= optModel.thethaUB_yL[j])
            optModel.linearize_thethaUB_yL_const3 = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optMode, j: optModel.thethaUB_yL[j] <= optModel.thethaUB[j] - (1 - optModel.yL[j])*0)

        if self.dual_formulation_type.lower() == 'simplified':
            optModel.muLB_yX_const = Constraint([j for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[0] == -self._bigM_rxnflux], rule = lambda optModel,j:optModel.muLB[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.muUB_yX_const = Constraint([j for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[1] == self._bigM_rxnflux], rule = lambda optModel,j: optModel.muUB[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.thethaLB_yU_const = Constraint(optModel.J, rule = lambda optModel, j: optModel.thethaLB[j] <= self._bigM_dualvar*optModel.yU[j])
            optModel.thethaUB_yL_const = Constraint(optModel.J, rule = lambda optModel,j: optModel.thethaUB[j] <= self._bigM_dualvar*optModel.yL[j])

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
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds)

        # Binary variables in for reactions considered for up-regulation, down-regulation or knockout 
        optModel.yL = Var(optModel.J, domain=Boolean)
        optModel.yU = Var(optModel.J, domain=Boolean)
        optModel.yX = Var(optModel.J, domain=Boolean)

        #-- Variables of the dual problem --
        # Dual variables associated with steady-state mass balance constraints
        # It's better to not provide upper and lower bound on lambda otherwise the dual objective 
        # will be slightly different from the primal's
        optModel.Lambda = Var(optModel.I, domain=Reals)

        # Dual variables associated with v_j >= LB_j*(1-yX(j) and v_j <= UB_j*(1-yX(j))
        optModel.muLB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        optModel.muUB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        if self.dual_formulation_type.lower() == 'standard':
            optModel.muLB_yX = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
            optModel.muUB_yX = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))

        # Dual variables associated with U_rxns_const: v(j) >= (LBon(j) - (-M))*yU(j) + (-M) (thethaLB) and  
        # L_rxns_const: v(j) <= (UBon(j) - bigM)*yU(j) + bigM (thethaUB)
        optModel.thethaLB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        optModel.thethaUB = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
        if self.dual_formulation_type.lower() == 'standard':
            optModel.thethaLB_yU = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))
            optModel.thethaUB_yL = Var(optModel.J, domain=Reals, bounds = (0,self._bigM_dualvar))

        #------ Objective function of the outer problem ------
        optModel.objectiveFunc = Objective(rule = lambda optModel: optModel.v[self.product_exchrxn_id], sense = maximize)

        #---------------- Constraints ----------------------
        #-- Constraints of the outer problem --
        # Each rxn can participate in only one intervention: yL(j) + yU(j) + yX(j) =l= 1
        optModel.one_interven_forEachRxn_const = Constraint(optModel.J, rule = lambda optModel,j: optModel.yL[j] + optModel.yU[j] + optModel.yX[j] <= 1)

        # Restrict the total number of interventions: sum(j,yL(j) + yU(j) + yX(j)) =l= interven_num_curr;
        optModel.total_interven_const = Constraint(rule = lambda optModel: sum([optModel.yL[j] + optModel.yU[j] + optModel.yX[j] for j in optModel.J]) <= self._interven_num_curr)

        # Maximum allowed number of interventions r out of the must sets
        optModel.notMUST_interven_const = Constraint(rule = lambda optModel: sum([optModel.yX[j] for j in optModel.J if j not in self.X_rxns]) + sum([optModel.yL[j] for j in optModel.J if j not in self.L_rxns]) + sum([optModel.yU[j] for j in optModel.J if j not in self.U_rxns]) <= self.notMUST_total_interven_num)

        # Maximum allowed number of knockouts, up-regulations and down-regulations outside of X_rxns, U_rxns and L_rxns
        optModel.notXrxns_interven_const = Constraint(rule = lambda optModel: sum([optModel.yX[j] for j in optModel.J if j not in self.X_rxns]) <= self.notXrxns_interven_num) 
        optModel.notUrxns_interven_const = Constraint(rule = lambda optModel: sum([optModel.yU[j] for j in optModel.J if j not in self.U_rxns]) <= self.notUrxns_interven_num) 
        optModel.notLrxns_interven_const = Constraint(rule = lambda optModel: sum([optModel.yL[j] for j in optModel.J if j not in self.L_rxns]) <= self.notLrxns_interven_num) 

        # Constraint on the max product formation flux achieved in earlier iterations of the code
        optModel.max_product_prev = Constraint(rule = lambda optModel: optModel.v[self.product_exchrxn_id] >= 1.005*self._best_product_yield_soFar)

        # Integer cuts
        optModel.integer_cuts = ConstraintList(noruleinit=True)

        #-- Constraints of the primal problem --
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)

        # Constrain for U_rxns
        # v(j) =g= LBon(j)*yU(j) + (1-yU(j))*LB(j) or v(j) =g= (LBon(j) - LB(j))*yU(j) + LB(j).
        # Instead of LB(j) use -bigM to make the constraint inactive when yU(j) = 0:
        # v(j) =g= (LBon(j) - (-M))*yU(j) + (-M) 
        optModel.U_rxns_const = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.v[j] >= (self.flux_bounds_overprod[j][0] + self._bigM_rxnflux)*optModel.yU[j] - self._bigM_rxnflux)

        # Constrain for L_rxns
        # v(j) =l= UBon(j)*yL(j) + (1-yL(j))*UB(j) or v(j) =l= (UBon(j) - UB(j))*yL(j) + UB(j)
        # Instead of UB(j) use bigM to make the constraint inactive when yL(j) = 0:
        # v(j) =l= (UBon(j) - bigM)*yU(j) + bigM
        optModel.L_rxns_const = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optModel, j: optModel.v[j] <= (self.flux_bounds_overprod[j][1] - self._bigM_rxnflux)*optModel.yL[j] + self._bigM_rxnflux)

        # Cconstraints for X_rxns. This constraint is written over all reactions (not just thos ein X_rxns) to incorporate
        # their LB and UB into the calculations
        # LB(j)*(1 - yX(j)) <= v(j)  and v(j) <= UB(j)*(1 - yX(j)) 
        optModel.X_rxns_const1 = Constraint(optModel.J, rule = lambda optModel, j: optModel.v[j] >= self.model.reactions_by_id[j].flux_bounds[0]*(1 - optModel.yX[j]) )
        optModel.X_rxns_const2 = Constraint(optModel.J, rule = lambda optModel, j: optModel.v[j] <= self.model.reactions_by_id[j].flux_bounds[1]*(1 - optModel.yX[j]))

        #-- Constraints of the dual problem --
        # dual constraints 
        optModel.dual_const = Constraint(optModel.J, rule=self.dual_const_rule)

        # Constraints linearizing the product of v_j and binary variables
        # x - (1-y)*UpperBound_of_x <= s <= x- (1-y)*LowerBound_of_x 
        # muLB[j] - (1 - yX[j])*muLB_max[j] <= muLB_yX[j] <= muLB[j] - (1 - yX[j])*muLB_min[j] and muLB_min = 0
        if self.dual_formulation_type.lower() == 'standard':
            optModel.linearize_muLB_yX_const1 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muLB_yX[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.linearize_muLB_yX_const2 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muLB[j] - (1 - optModel.yX[j])*self._bigM_dualvar <= optModel.muLB_yX[j])
            optModel.linearize_muLB_yX_const3 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muLB_yX[j] <= optModel.muLB[j] - (1 - optModel.yX[j])*0)

            # muUB_min[j]*yX[j] <= muUB_yX[j] <= muUB_max*yX[j]  and muUB_min = 0
            optModel.linearize_muUB_yX_const1 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muUB_yX[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.linearize_muUB_yX_const2 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muUB[j] - (1 - optModel.yX[j])*self._bigM_dualvar <= optModel.muUB_yX[j])
            optModel.linearize_muUB_yX_const3 = Constraint([j for j in optModel.J if j not in self.ignored_X_rxns], rule = lambda optModel, j: optModel.muUB_yX[j] <= optModel.muUB[j] - (1 - optModel.yX[j])*0)

            # thethaLB_min[j]*yU[j] <= thethaLB_yU[j] <= thethaLB_max*yU[j]  and thethaLB_min = 0
            optModel.linearize_thethaLB_yU_const1 = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.thethaLB_yU[j] <= self._bigM_dualvar*optModel.yU[j])
            optModel.linearize_thethaLB_yU_const2 = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.thethaLB[j] - (1 - optModel.yU[j])*self._bigM_dualvar <= optModel.thethaLB_yU[j])
            optModel.linearize_thethaLB_yU_const3 = Constraint([j for j in optModel.J if j not in self.ignored_U_rxns], rule = lambda optModel, j: optModel.thethaLB_yU[j] <= optModel.thethaLB[j] - (1 - optModel.yU[j])*0)

            # thethaUB_min[j]*yL[j] <= thethaUB_yL[j] <= thethaUB_max*yL[j]  and thethaUB_min = 0
            optModel.linearize_thethaUB_yL_const1 = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optModel, j: optModel.thethaUB_yL[j] <= self._bigM_dualvar*optModel.yL[j])
            optModel.linearize_thethaUB_yL_const2 = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optModel, j: optModel.thethaUB[j] - (1 - optModel.yL[j])*self._bigM_dualvar <= optModel.thethaUB_yL[j])
            optModel.linearize_thethaUB_yL_const3 = Constraint([j for j in optModel.J if j not in self.ignored_L_rxns], rule = lambda optMode, j: optModel.thethaUB_yL[j] <= optModel.thethaUB[j] - (1 - optModel.yL[j])*0)

        if self.dual_formulation_type.lower() == 'simplified':
            optModel.muLB_yX_const = Constraint([j for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[0] == -self._bigM_rxnflux], rule = lambda optModel,j:optModel.muLB[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.muUB_yX_const = Constraint([j for j in optModel.J if self.model.reactions_by_id[j].flux_bounds[1] == self._bigM_rxnflux], rule = lambda optModel,j: optModel.muUB[j] <= self._bigM_dualvar*optModel.yX[j])
            optModel.thethaLB_yU_const = Constraint(optModel.J, rule = lambda optModel, j: optModel.thethaLB[j] <= self._bigM_dualvar*optModel.yU[j])
            optModel.thethaUB_yL_const = Constraint(optModel.J, rule = lambda optModel,j: optModel.thethaUB[j] <= self._bigM_dualvar*optModel.yL[j])

        #-- Strong duality -- 
        optModel.strong_duality_const = Constraint(rule=self.strong_duality_const_rule)

        self.optModel = optModel

    def fix_known_variables(self):
        """
        Fixies all known binary variables to zero/one
        """
        # Do not manipulate reactions in self.ignorex_X_rxns
        for rxn in self.ignored_X_rxns:
            self.optModel.yX[rxn] = 0
            self.optModel.yX[rxn].fixed = True
            if self.dual_formulation_type.lower() == 'standard':
                self.optModel.muLB_yX[rxn] = 0
                self.optModel.muLB_yX[rxn].fixed = True
                self.optModel.muUB_yX[rxn] = 0
                self.optModel.muUB_yX[rxn].fixed = True

        for rxn in self.ignored_L_rxns:
            self.optModel.yL[rxn] = 0
            self.optModel.yL[rxn].fixed = True
            if self.dual_formulation_type.lower() == 'standard':
                self.optModel.thethaUB_yL[rxn] = 0 
                self.optModel.thethaUB_yL[rxn].fixed = True

        for rxn in self.ignored_U_rxns:
            self.optModel.yU[rxn] = 0
            self.optModel.yU[rxn].fixed = True
            if self.dual_formulation_type.lower() == 'standard':
                optModel.thethaLB_yU[rxn] = 0 
                optModel.thethaLB_yU[rxn].fixed = True 


        # Fix binary variables for a priori imposed knockouts, down-regulations and up-regulations
        for rxn in self.fixed_X_rxns:
            self.optModel.yX[rxn] = 1 
            self.optModel.yX[rxn].fixed = True

        for rxn in self.fixed_L_rxns:
            self.optModel.yL[rxn] = 1 
            self.optModel.yL[rxn].fixed = True

        for rxn in self.fixed_U_rxns:
            self.optModel.yU[rxn] = 1 
            self.optModel.yU[rxn].fixed = True

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
            solver_flag = 'solverError'
            if self.warnings:
                print '**WARNING (FORCE)! {} failed with the following error: \n{} \n'.format(self.optimization_solver,e)

        elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
        elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))

        if solver_flag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':

            exit_flag = 'globallyOptimal'

            # Load the results
            self.optModel.load(optSoln)

            # Optimal value of the objective function
            opt_objValue = self.optModel.objectiveFunc()

            # Reactions that must be knocked out, up-regulated or down-regulated
            yX_one_rxns = [j for j in self.optModel.J if abs(self.optModel.yX[j].value - 1) <= mip_integrality_tol]
            yL_one_rxns = [j for j in self.optModel.J if abs(self.optModel.yL[j].value - 1) <= mip_integrality_tol]
            yU_one_rxns = [j for j in self.optModel.J if abs(self.optModel.yU[j].value - 1) <= mip_integrality_tol]

        # If there was a solver error or if an optimal solution was not returned 
        else:
            opt_objValue = None
            if solver_flag == 'solverError':
                exit_flag = solver_flag
            else:
                exit_flag = str(optSoln.solver.termination_condition)
 
            # Reactions that must be knocked out, up-regulated or down-regulated
            yX_one_rxns = []
            yL_one_rxns = []
            yU_one_rxns = []

        # Store the solution
        self._curr_soln = {'exit_flag':exit_flag,'objective_value':opt_objValue,'interven_num': len(yX_one_rxns) + len(yL_one_rxns) + len(yU_one_rxns),'X_rxns':yX_one_rxns,'L_rxns':yL_one_rxns,'U_rxns':yU_one_rxns}

        # Print the results on the screen 
        if self.stdout_msgs_details:
            print '\nObjective value = {}, Optimality status = {}, Solution status = {}, Solver run status = {}'.format(opt_objValue, optSoln.solver.termination_condition, optSoln.Solution.status, solver_flag)
            print 'Took (hh:mm:ss) {}/{} of processing/walltime to create a pyomo model, {}/{} to  preprcoess the model and {}/{} to solve the model\n'.format(self._elapsed_create_optModel_pt, self._elapsed_create_optModel_wt, elapsed_preproc_pyomo_pt,elapsed_preproc_pyomo_wt, elapsed_solver_pt,elapsed_solver_wt)

    def find_min_product_yield(self, validate = True):
        """
        If validate is True it Validates an obtained solution. 
        if validate is False, it finds the min product yield for a set of fixed binary
        variables provided by self.fixed_X_rxns, self.fixed_L_rxns and self.fixed_U_rxns
        """
        for rxn in self.model.reactions:
            rxn.objective_coefficient = 0
        self.model.reactions_by_id[self.product_exchrxn_id].objective_coefficient = 1                

        if validate:
            for j in self._curr_soln['X_rxns']:
                self.model.reactions_by_id[j].flux_bounds = [0,0]
            for j in self._curr_soln['L_rxns']:
                self.model.reactions_by_id[j].flux_bounds[1] = self.flux_bounds_overprod[j][1]
            for j in self._curr_soln['U_rxns']:
                self.model.reactions_by_id[j].flux_bounds[0] = self.flux_bounds_overprod[j][0]
        else:
            for j in self.fixed_X_rxns:
                self.model.reactions_by_id[j].flux_bounds = [0,0]
            for j in self.fixed_L_rxns:
                self.model.reactions_by_id[j].flux_bounds[1] = self.flux_bounds_overprod[j][1]
            for j in self.fixed_U_rxns:
                self.model.reactions_by_id[j].flux_bounds[0] = self.flux_bounds_overprod[j][0]

        if hasattr(self.model,'fba_model'):
            self.model.fba(build_new_optModel = False, maximize = False, stdout_msgs = False)
        else:
            self.model.fba(build_new_optModel = True, maximize = False, stdout_msgs = False)

        if self.model.fba_model.solution['exit_flag'] == 'globallyOptimal': 
            if validate and (self._curr_soln['objective_value'] - self.model.fba_model.solution['objective_value'] > 1e-6):
                raise userError('Validation failed for X_rxns = {} , L_rxns = {} , U_rxns = {} because the production flux of the product ({}) is less than the objective function value of the optimizaiton problem for identifhing FORCE sets ({})'.format(self._curr_soln['X_rxns'], self._curr_soln['L_rxns'], self._curr_soln['U_rxns'], self.model.fba_model.solution['objective_value'], self._curr_soln['objective_value']))
            elif not validate:
                self._best_product_yield_soFar_init = self.model.fba_model.solution['objective_value']
            else:
                if self.stdout_msgs_details:
                    print 'The following soluiton was successfully validates: X_rxns = {}\nL_rxns = {}\nU_rxns = {}'.format(self._curr_soln['X_rxns'],self._curr_soln['L_rxns'],self._curr_soln['U_rxns'])

        elif self.model.fba_model.solution['exit_flag'] != 'globallyOptimal':
            raise userError('Validaiton failed for X_rxns = {} , L_rxns = {} , U_rxns = {} because the fba problem to find the minimum production flux of the product was not solved to optimality and ended with an exit_flag of {}'.format(self._curr_soln['X_rxns'], self._curr_soln['L_rxns'], self._curr_soln['U_rxns'], self.model.fba_model.solution['exit_flag']))

        else:  # If validated, then reset the flux bounds for the model
            self.set_model_flux_bounds()

    def run(self):
        """ 
        This method runs FBA. 

        OUTPUT:
        -------
        solution: A list of dictionaries with the following keys:
        """
        # Total processing and wall time required to create the pyomo model, solve it and store the results 
        start_total_pt = time.clock()
        start_total_wt = time.time()

        # Set flux bounds for the model
        self.set_model_flux_bounds()

        # Number of solutions found so far
        found_solutions_num = 0

        # Current max number of interventions
        self._interven_num_curr = len(self.fixed_X_rxns) + len(self.fixed_L_rxns) + len(self.fixed_U_rxns) + 1

        # Find the min product yield, if a number of fixed interventions has already provided using fixed_X_rxns, 
        # fixed_L-rxns or fixed_U_rxns,
        if len(self.fixed_X_rxns) > 0 or len(self.fixed_L_rxns) > 0 or len(self.fixed_U_rxns) > 0:
            self.find_min_product_yield(validate = False)
        else:
            self._best_product_yield_soFar_init = 0.0
        if self.stdout_msgs:
            print '\nMin product yield = {:0.4} ({:0.3}% of theoretical maximum = {:0.4}) , biomass flux = {} ({:0.3}%) of max biomass = {:0.4})'.format(self._best_product_yield_soFar_init, 100*self._best_product_yield_soFar_init/self.product_max_theor_yield, self.product_max_theor_yield, self.model.fba_model.optModel.v[self.biomass_rxn_id].value, 100*self.model.fba_model.optModel.v[self.biomass_rxn_id].value/self.max_biomass_flux, self.max_biomass_flux)

        # Best product yield achieved so far
        self._best_product_yield_soFar = self._best_product_yield_soFar_init

        # The following parameter is essentially the same as self._best_product_yield_soFar but it is
        # reset to zero each time self._interven_num_curr increases. This parameter store the best
        # objective function with the current number of interventions. 
        best_product_yield_soFar = self._best_product_yield_soFar_init

        if self.results_filename != '':
            with open(self.results_filename,'w') as f:
                f.write('')

        # Creating the pyomo optModel 
        if self.build_new_optModel:
            start_pyomo_pt = time.clock()
            start_pyomo_wt = time.time()

            self.build_bilevel_optModel()

            self._elapsed_create_optModel_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
            self._elapsed_create_optModel_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))
        else:
            self._elapsed_create_optModel_pt = 0 
            self._elapsed_create_optModel_wt = 0 

        # Fix known variables
        self.fix_known_variables()

        # Create a solver and set the options
        self._optSolver = pyomoSolverCreator(self.optimization_solver)

        # A list of dictionaries holding the optimal solution in different iterations
        self.solutions = []
        done = False

        if self.stdout_msgs:
            print '\n-------- # of interventions = {} ----------\n'.format(self._interven_num_curr)

        # For a given number of interventions we keep finding alternative solutions all giving the same objective 
        # function value. If the problem becomes infeasible or the objective funciton is less than the best current
        # then we increase the number of iterventions 
        while not done: 

            # Solve the optimizaiton model 
            self.run_single()

            # Check whether the current run was successful
            if self._curr_soln['exit_flag'] == 'globallyOptimal':

                # If the product yield for the current number of iterations is zero, move on  
                # by increasing the number of interventions
                if abs(self._curr_soln['objective_value'] - self._best_product_yield_soFar_init) < 1e-6:
                    if self.stdout_msgs:
                        print '\n** Iterations with {} interventions stopped with the following solution, which results in an objective function of almost zero'.format(self._interven_num_curr)
                        self.print_results_summary()

                    self._interven_num_curr += 1
                    if self.stdout_msgs and self._interven_num_curr <= self.total_interven_num:
                        print '\n-------- # of interventions = {} ----------\n'.format(self._interven_num_curr)
                    best_product_yield_soFar = self._best_product_yield_soFar_init 
                    self._best_product_yield_soFar = self._best_product_yield_soFar_init 
         
                # If the product yield with the current number of iterations is less than best_product_yield_soFar
                # then move on by increasing the number of iterations  
                elif best_product_yield_soFar - self._curr_soln['objective_value'] > 1e-3:
                    # print a summary of results in the output
                    if self.stdout_msgs:
                        print '\n** Iterations with {} interventions stopped with the following solution, which is worse than the best solution achieved so far with this number of interventions ({})'.format(self._interven_num_curr, best_product_yield_soFar)
                        self.print_results_summary()

                    # Stop if the target yield has already been achieved
                    if best_product_yield_soFar >= (self.stopWith_product_yield_percent/100)*self.product_max_theor_yield:
                        done = True
                        if stdout_msgs:
                            print 'Iterations stopped because the best product yield achieved so far ({}) is satisfies the target product yield ({})'.format(best_product_yield_soFar, (self.stopWith_product_yield_percent/100)*self.product_max_theor_yield)
                    else:
                        self._interven_num_curr += 1
                        if self._interven_num_curr <= self.total_interven_num:
                            if self.stdout_msgs and self._interven_num_curr <= self.total_interven_num:
                                print '\n-------- # of interventions = {} ----------\n'.format(self._interven_num_curr)
                            self._best_product_yield_soFar = best_product_yield_soFar
                            #best_product_yield_soFar = 0  # reset to zero before moving to higher order interventions 
   
                # Otherwise store the results
                else:
                    found_solutions_num += 1

                    # If this is the first solution with the current number of iterations specified 
                    # by self._interven_num_curr, then set best_product_yield_soFar to the objective funciton
                    # value, if it is greater than zero
                    if self._curr_soln['objective_value'] > best_product_yield_soFar:
                        best_product_yield_soFar = self._curr_soln['objective_value']

                    # Validate the solution
                    if self.validate_results:
                        self.find_min_product_yield()

                    self.solutions.append(self._curr_soln)

                    # Add an integer cut excluding the current solution. 
                    # we check if their distance from one is less than mip_integrality_tol
                    # Total number of yX, yL and yU variables is 3*len(self.optModel.J)
                    self.optModel.integer_cuts.add(
                         sum([self.optModel.yX[j] for j in self._curr_soln['X_rxns']]) + \
                         sum([self.optModel.yL[j] for j in self._curr_soln['L_rxns']]) + \
                         sum([self.optModel.yU[j] for j in self._curr_soln['U_rxns']]) + \
                         sum([1 - self.optModel.yX[j] for j in self.optModel.J if j not in self._curr_soln['X_rxns']]) + \
                         sum([1 - self.optModel.yL[j] for j in self.optModel.J if j not in self._curr_soln['L_rxns']]) + \
                         sum([1 - self.optModel.yU[j] for j in self.optModel.J if j not in self._curr_soln['U_rxns']]) <= \
                         3*len(self.optModel.J) - 1)

                    # print a summary of results in the output
                    if self.stdout_msgs:
                        print '\n{})'.format(found_solutions_num)
                        self.print_results_summary()
                        print '\nAdded integer cut (not including the ones that were zero): {}'.format(sum([self.optModel.yX[j] for j in self._curr_soln['X_rxns']]) + sum([self.optModel.yL[j] for j in self._curr_soln['L_rxns']]) + sum([self.optModel.yU[j] for j in self._curr_soln['U_rxns']]) <= len(self.optModel.J) - 1)
                        print 'max_product_prev const: {}'.format(self.optModel.v[self.product_exchrxn_id] >= 1.005*self._best_product_yield_soFar)

                    # Write results into the file
                    if self.results_filename != '':
                       with open(self.results_filename,'a') as f:
                           f.write('X_rxns = [\n') 
                           for r_id in self._curr_soln['X_rxns']:
                               rxn = self.model.reactions_by_id[r_id] # reaction object
                               f.write("{{'id':'{}','name','{}', 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(r_id,rxn.name,rxn.gene_reaction_rule, self.flux_bounds_ref[r_id], self.flux_bounds_overprod[r_id]))
                           f.write(']\n') 

                           f.write('L_rxns = [\n') 
                           for r_id in self._curr_soln['L_rxns']:
                               rxn = self.model.reactions_by_id[r_id] # reaction object
                               f.write("{{'id':'{}','name','{}', 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(r_id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[r_id], self.flux_bounds_overprod[r_id]))
                           f.write(']\n') 

                           f.write('U_rxns = [\n') 
                           for r_id in self._curr_soln['U_rxns']:
                               rxn = self.model.reactions_by_id[r_id] # reaction object
                               f.write("{{'id':'{}','name','{}', 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(r_id,rxn.name,rxn.gene_reaction_rule, self.flux_bounds_ref[r_id], self.flux_bounds_overprod[r_id]))
                           f.write(']\n') 

            # If the optimization problem was not solved to optimality 
            else:
                if self.stdout_msgs:
                    print '**Iterations with {} interventions stopped because the optimization problem was not solved to optimality'.format(self._interven_num_curr)

                # Stop if the target yield has already been achieved
                if best_product_yield_soFar >= (self.stopWith_product_yield_percent/100)*self.product_max_theor_yield:
                    done = True
                    if stdout_msgs:
                        print 'Iterations stopped because the best product yield achieved so far ({}) is satisfies the target product yield ({})'.format(best_product_yield_soFar, (self.stopWith_product_yield_percent/100)*self.product_max_theor_yield)
                else: 
                    self._interven_num_curr += 1
                    if self._interven_num_curr <= self.total_interven_num:
                        if self.stdout_msgs:
                            print '\n-------- # of interventions = {} ----------\n'.format(self._interven_num_curr)
                        best_product_yield_soFar = 0 

            # Stop if the maximum number of iterations has exceeded
            if self._interven_num_curr > self.total_interven_num:
                if self.stdout_msgs and not done:
                    print '\n**Iterations stopped because the maximum number of interventions = {} has been reached\n'.format(self.total_interven_num)
                done = True

        elapsed_total_pt = str(timedelta(seconds = time.clock() - start_total_pt))
        elapsed_total_wt = str(timedelta(seconds = time.time() - start_total_wt))
        if self.stdout_msgs:
            print '\nFound a total of {} FORCE solutions'.format(found_solutions_num)
            print 'Finding FORCE sets took (hh:mm:ss) a total of {}/{} of processing/walltime\n'.format(elapsed_total_pt, elapsed_total_wt)


    def print_results_summary(self):
        """
        prints a summary of the current results in the output
        """
        print '\nMaxmin product yield = {:0.4} ({:0.3}% of theoretical maximum = {:0.4})  , biomass flux = {} ({:0.3}%) of max biomass = {:0.4})'.format(self.optModel.v[self.product_exchrxn_id].value, 100*self.optModel.v[self.product_exchrxn_id].value/self.product_max_theor_yield, self.product_max_theor_yield, self.optModel.v[self.biomass_rxn_id].value, 100*self.optModel.v[self.biomass_rxn_id].value/self.max_biomass_flux, self.max_biomass_flux)

        if len(self._curr_soln['X_rxns']) > 0:
             print 'Knockouts = \n{}'.format(self._curr_soln['X_rxns'])
        else:
             print 'Knockouts = \nNone'

        print '\nDown-regulations = '
        if len(self._curr_soln['L_rxns']) > 0:
            for rxn in self._curr_soln['L_rxns']:
                print '{}:\tflux_bounds_ref = {} ,  flux_bounds_overprod = {}'.format(rxn,self.flux_bounds_ref[rxn],self.flux_bounds_overprod[rxn])
        else:
            print 'None'

        print '\nUp-regulations ='
        if len(self._curr_soln['U_rxns']) > 0:
            for rxn in self._curr_soln['U_rxns']:
                print '{}:\tflux_bounds_ref = {} ,  flux_bounds_overprod = {}'.format(rxn,self.flux_bounds_ref[rxn],self.flux_bounds_overprod[rxn])
        else:
            print 'None' 

    def run_primal_dual(self):
        """ 
        Runs primal and dual problems to test if they are equal 
        """
        # Total processing and wall time required to create the pyomo model, solve it and store the results 

        # Set flux bounds for the model
        self.set_model_flux_bounds()
        self.define_optModel_params()
   
        # Empty ignored reacitons
        self.ignored_X_rxns = []
        self.ignored_L_rxns = []
        self.ignored_U_rxns = []

        # Create a solver and set the options
        self._optSolver = pyomoSolverCreator(self.optimization_solver)

        print '--------- Primal problem --------'
        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()
        self.build_primal_optModel()
        self._elapsed_create_optModel_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
        self._elapsed_create_optModel_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))

        for j in self.optModel.J:
            self.optModel.yX[j] = 0
            self.optModel.yL[j] = 0
            self.optModel.yU[j] = 0

        self.optModel.yU['PylB'] = 1

        for j in self.optModel.J:
            self.optModel.yX[j].fixed = True
            self.optModel.yL[j].fixed = True
            self.optModel.yU[j].fixed = True

        # Solve the optimizaiton model 
        self.run_single()

        print '--------- Dual problem --------'
        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()
        self.build_dual_optModel()
        self._elapsed_create_optModel_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
        self._elapsed_create_optModel_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))

        for j in self.optModel.J:
            self.optModel.yX[j] = 0
            self.optModel.yL[j] = 0
            self.optModel.yU[j] = 0

        self.optModel.yU['PylB'] = 1

        for j in self.optModel.J:
            self.optModel.yX[j].fixed = True
            self.optModel.yL[j].fixed = True
            self.optModel.yU[j].fixed = True

        # Solve the optimizaiton model 
        self.run_single()
#-------------------------
if __name__ == '__main__':
    pass 
 
