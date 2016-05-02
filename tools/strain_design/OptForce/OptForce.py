from __future__ import division
import sys,os, time
sys.path.append('../../../')
from copy import deepcopy
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from shared_data_holder import shared_data_holder
from MUST_doubles import MUST_doubles
from FORCE import FORCE
from tools.userError import userError
from tools.globalVariables import *
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.create_model import create_model
from tools.fba.fba import fba
from tools.fba.fva import fva
from tools.fba.set_specific_bounds import set_specific_bounds
from imp import load_source
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class OptForce(object):
    """
    Checks whether a given flux data leads to a feasible FBA for a given metabolic model
    """
    def __init__(self, model, product_exchrxn_id, product_targetYield_percent = 80, min_biomass_percent = 10, growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, read_exp_flux_bounds_ref_fromThisFile = '',  read_blocked_rxns_fromThisFile = '', read_inSilico_essential_rxns_fromThisFile = '', read_inVivo_essential_rxns_fromThisFile = '', save_flux_bounds_ref_toThisFile = '', read_flux_bounds_ref_fromThisFile = '', save_flux_bounds_overprod_toThisFile = '', read_flux_bounds_overprod_fromThisFile = '', save_MUST_singles_toThisFile = '', read_MUST_singles_fromThisFile = '', read_MUST_doubles_fromThisFile = '', MUST_singles_diff_thr = 1, MUST_doubles_params = {}, FORCE_params = {}, optimization_solver = default_optim_solver, simulation_condition = '', stdout_msgs = True, warnings = True):
        """
        INPUTS:
        ------
                          model: An instance of class model
             product_exchrxn_id: The id of the exchange reaction for the product of interest (a string)
                exp_flux_bounds: A dictionary with keys and values as follows:
                                   keys: Id of reactions with available flux data
                                 values: A list of two elements in the form [LB,UB]
                                         containing the LB and UB on reaction fluxes. 
                                         LB and UB can be None as well (if there is no
                                         experimental data for reaction fluxes)
        product_targetYield_percent: Desired target yield for the product of interest in the form of the percent of maximum theoretical yield
            min_biomass_percent: Min required biomass formation flux in the form of the percent of the max theoretical biomass formation flux
        growthMedium_flux_bounds: Information about growth media
                                 read_flux_bounds_ref_fromThisFile: Name of and path of the file containing the flux bounds
                                                      file for exchange reactions corresponding to compounds in 
                                                       in the growth media
                                 flux_bounds_dict: A dictionary containing the flux bounds for other reactions (such as carbon
                                                  source uptake, oxygen uptake, etc), that have to be set for reactions not 
                                                   in flux_data_filenam or media_filename
                   read_blocked_rxns_fromThisFile:
        read_inSilico_essential_rxns_fromThisFile:
          read_inVivo_essential_rxns_fromThisFile:
                                  These are the files in which the list of blocked reactions under the examined condition, list of 
                                  essential rxns in silico and essential rxns in vivo are stored 
               save_flux_bounds_ref_toThisFile: 
          save_flux_bounds_overprod_toThisFile: 
                  save_MUST_singles_toThisFile: 
                                  These are the files to save the the results in. If an empty string is provided
                                  the results will not be saved in files
             read_flux_bounds_ref_fromThisFile: 
        read_flux_bounds_overprod_fromThisFile: 
                read_MUST_singles_fromThisFile: 
                read_MUST_doubles_fromThisFile: 
                                 These are the files to read the flux bounds, MUST singles and MUST doulbes from.
                                 Check out function find_flux_ref_overprod_flux_bounds and class MUST_doubles > run
                                 for details of how the results shoud be stored in files.
            : A dictionary containing the specific input parameters for the MUST_doubles class
                   FORCE_params: A dictionary containing the specific input parameters for the FORCE class
               MUST_singles_diff_thr: The minimum value above which between the difference between LB/UB in
                                 the reference and overproducing strains are considered to be 
                                 signficant for declating a MUST single. For example, if 
                                 LB_overprod(j) - UB_ref(j) >= MUST_singles_diff_thr, then j is a MUST_U
            optimization_solver: A string containing showing which optimization solver to use
                                 Allowed choices are 'gurobi' and 'cplex'
           simulation_condition: A string representing the simulation conditon 
                       warnings: Can be True or False shwoing whether the warnings should be 
                                 writtten to the 
                                 screen or not. The default is True  
                    stdout_msgs: Prints a summary of the FBA and a number of other messages 
                                 in the output if 'on.
       
        Ali R. Zomorrodi - Segre's Lab @ BU
        Last updated: 04-04-2016
        """
        # model 
        self.model = model 
        self.biomass_reaction = model.biomass_reaction
        self.biomass_rxn_id = model.biomass_reaction.id

        # Growth medium
        self.growthMedium_flux_bounds = growthMedium_flux_bounds

        # Product and biomass yields
        self.product_exchrxn_id = product_exchrxn_id 
        self.product_targetYield_percent = product_targetYield_percent
        self.min_biomass_percent = min_biomass_percent

        # Parameters for MUST_singles, MUST_doubles and FORCE classes
        self.MUST_singles_diff_thr = MUST_singles_diff_thr
        self.MUST_doubles_params = MUST_doubles_params
        self.FORCE_params = FORCE_params

        # Other general parameters
        self.optimization_solver =  optimization_solver
        self.simulation_condition = simulation_condition
        self.warnings = warnings
        self.stdout_msgs =  stdout_msgs

        # Input and output files
        self.read_exp_flux_bounds_ref_fromThisFile = read_exp_flux_bounds_ref_fromThisFile
        self.read_blocked_rxns_fromThisFile = read_blocked_rxns_fromThisFile
        self.read_inSilico_essential_rxns_fromThisFile = read_inSilico_essential_rxns_fromThisFile
        self.read_inVivo_essential_rxns_fromThisFile = read_inVivo_essential_rxns_fromThisFile
        self.read_exp_flux_bounds_ref_fromThisFile = read_exp_flux_bounds_ref_fromThisFile
        self.save_flux_bounds_ref_toThisFile = save_flux_bounds_ref_toThisFile
        self.read_flux_bounds_ref_fromThisFile = read_flux_bounds_ref_fromThisFile
        self.save_flux_bounds_overprod_toThisFile = save_flux_bounds_overprod_toThisFile
        self.read_flux_bounds_overprod_fromThisFile = read_flux_bounds_overprod_fromThisFile
        self.save_MUST_singles_toThisFile = save_MUST_singles_toThisFile
        self.read_MUST_singles_fromThisFile = read_MUST_singles_fromThisFile
        self.read_MUST_doubles_fromThisFile = read_MUST_doubles_fromThisFile

        #-- Load some data from files --
        # Experimental flux bounds 
        if self.read_exp_flux_bounds_ref_fromThisFile != '':
            load_source('dataFile',self.read_exp_flux_bounds_ref_fromThisFile)
            import dataFile
            self.exp_flux_bounds_ref = dataFile.exp_flux_bounds 
        else: 
            self.exp_flux_bounds_ref = []
 
        # In vivo essential reactions
        if self.read_inVivo_essential_rxns_fromThisFile != '':
            load_source('dataFile',self.read_inVivo_essential_rxns_fromThisFile)
            import dataFile
            self.inVivo_essential_rxns = dataFile.inVivo_essential_rxns 
        else: 
            self.inVivo_essential_rxns = []
 
        # In silico essential reactions
        if self.read_inSilico_essential_rxns_fromThisFile != '':
            load_source('dataFile',self.read_inSilico_essential_rxns_fromThisFile)
            import dataFile
            self.inSilico_essential_rxns = dataFile.essential_rxns 
        else: 
            self.inSilico_essential_rxns = []
 
        # Blocked rxns under the examined uptake and aeration conditions 
        if self.read_blocked_rxns_fromThisFile != '':
            load_source('dataFile',self.read_blocked_rxns_fromThisFile)
            import dataFile
            self.blocked_rxns = dataFile.blocked_rxns 
        else: 
            self.blocked_rxns = []
 
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

        # product_targetYield_percent
        if attr_name == 'product_targetYield_percent' and (not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError('product_targetYield_percent must be an integer or a float')
        if attr_name == 'product_targetYield_percent' and (attr_value < 0 or attr_value > 100):
            raise ValueError('product_targetYield_percent must be between 0 and 100')

        # min_biomass_percent 
        if attr_name == 'min_biomass_percent' and (not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError('min_biomass_percent must be an integer or a float')
        if attr_name == 'min_biomass_percent' and (attr_value < 0 or attr_value > 100):
            raise ValueError('min_biomass_percent must be between 0 and 100')

        # growthMedium_flux_bounds 
        if attr_name == 'growthMedium_flux_bounds' and not isinstance(attr_value,dict):
            raise TypeError('growthMedium_flux_bounds must be a dictionary')
        if attr_name == 'growthMedium_flux_bounds' and len([k for k in attr_value.keys() if k.lower() not in ['flux_bounds_filename','flux_bounds_dict']]) > 0:
            raise ValueError('Invalid key for growthMedium_flux_bounds. Allowed keys are flux_bounds_filename and flux_bounds_dict')

        # file names
        if attr_name in ['read_exp_flux_bounds_ref_fromThisFile','read_blocked_rxns_fromThisFile','read_inSilico_essential_rxns_fromThisFile','read_inVivo_essential_rxns_fromThisFile', 'save_flux_bounds_ref_toThisFile','read_flux_bounds_ref_fromThisFile','save_flux_bounds_overprod_toThisFile','read_flux_bounds_overprod_fromThisFile','save_MUST_singles_toThisFile','read_MUST_singles_fromThisFile','read_MUST_doubles_fromThisFile'] and not isinstance(attr_value,str):
            raise TypeError('{} must be a string'.format(attr_name))

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi','gurobi_ampl','cplexamp']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # MUST_singles_diff_thr
        if attr_name == 'MUST_singles_diff_thr' and not isinstance(attr_value,float) and not isinstance(attr_value,int):
            raise TypeError('MUST_singles_diff_thr must be a non-negative integer or float')
        elif attr_name == 'MUST_singles_diff_thr' and attr_value < 0: 
            raise ValueError('MUST_singles_diff_thr must be a non-negative integer or float')

        # MUST_doubles params
        if attr_name == 'MUST_doubles_params' and not isinstance(attr_value,dict):
            raise TypeError('MUST_doubles_params must be a dictionary')
        elif attr_name == 'MUST_doubles_params' and len([k for k in attr_value.keys() if k.lower() not in ['growthMedium_flux_bounds', 'product_targetYield_percent', 'min_biomass_percent', 'blocked_rxns', 'MUST_X', 'MUST_L', 'MUST_U', 'build_new_optModel', 'objective_thr', 'validate_results', 'results_filename', 'create_new_results_file', 'optimization_solver','warnings', 'stdout_msgs']]) > 0:
            raise userError('Invalid key(s) for MUST_doubles: {}'.format([k for k in attr_value.keys() if k.lower() not in ['growthMedium_flux_bounds', 'product_targetYield_percent', 'min_biomass_percent', 'blocked_rxns', 'MUST_X', 'MUST_L', 'MUST_U', 'build_new_optModel', 'objective_thr', 'validate_results', 'results_filename', 'create_new_results_file', 'optimization_solver','warnings', 'stdout_msgs']]))
        elif attr_name == 'MUST_doubles_params':
            # Assign default values to MUST_doubles_params
            if 'product_targetYield_percent' not in attr_value.keys():
                attr_value['product_targetYield_percent'] = 80
            if 'min_biomass_percent' not in attr_value.keys():
                attr_value['min_biomass_percent'] = 10
            if 'blocked_rxns' not in attr_value.keys():
                attr_value['blocked_rxns'] = []
            if 'MUST_L' not in attr_value.keys():
                attr_value['MUST_L'] = []
            if 'MUST_U' not in attr_value.keys():
                attr_value['MUST_U'] = []
            if 'build_new_optModel' not in attr_value.keys():
                attr_value['build_new_optModel'] = True
            if 'objective_thr' not in attr_value.keys():
                attr_value['objective_thr'] = 1
            if 'validate_results' not in attr_value.keys():
                attr_value['validate_results'] = False
            if 'create_new_results_file' not in attr_value.keys():
                attr_value['create_new_results_file'] = True
            if 'optimization_solver' not in attr_value.keys():
                attr_value['optimization_solver'] = default_optim_solver
            if 'warnings' not in attr_value.keys():
                attr_value['warnings'] = True
            if 'stdout_msgs' not in attr_value.keys():
                attr_value['stdout_msgs'] = True

        if attr_name == 'FORCE_params' and not isinstance(attr_value,dict):
            raise TypeError('FORCE_params must be a dictionary')
        elif attr_name == 'FORCE_params' and len([k for k in attr_value.keys() if k not in ['growthMedium_flux_bounds', 'min_biomass_percent', 'stopWith_product_yield_percent', 'MUST_X', 'MUST_L', 'MUST_U', 'MUST_LU_L', 'MUST_LU_U', 'MUST_LL_L1', 'MUST_LL_L2', 'MUST_UU_U1', 'MUST_UU_U2', 'total_interven_num', 'notMUST_total_interven_num', 'notXrxns_interven_num', 'notLrxns_interven_num', 'notUrxns_interven_num', 'fixed_X_rxns', 'fixed_L_rxns', 'fixed_U_rxns', 'ignored_X_rxns', 'ignored_L_rxns', 'ignored_U_rxns', 'inSilico_essential_rxns', 'inVivo_essential_rxns', 'blocked_rxns', 'build_new_optModel', 'dual_formulation_type', 'validate_results', 'results_filename', 'warnings', 'stdout_msgs']]) > 0:
            raise ValueError('Invalid key(s) for FORCE_params: {}'.format([k for k in attr_value.keys() if k not in ['growthMedium_flux_bounds', 'min_biomass_percent', 'stopWith_product_yield_percent', 'MUST_X', 'MUST_L', 'MUST_U', 'MUST_LU_L', 'MUST_LU_U', 'MUST_LL_L1', 'MUST_LL_L2', 'MUST_UU_U1', 'MUST_UU_U2', 'total_interven_num', 'notMUST_total_interven_num', 'notXrxns_interven_num', 'notLrxns_interven_num', 'notUrxns_interven_num', 'fixed_X_rxns', 'fixed_L_rxns', 'fixed_U_rxns', 'ignored_X_rxns', 'ignored_L_rxns', 'ignored_U_rxns', 'inSilico_essential_rxns', 'inVivo_essential_rxns', 'blocked_rxns', 'build_new_optModel', 'dual_formulation_type', 'validate_results', 'results_filename', 'warnings', 'stdout_msgs']])) 
        elif attr_name == 'FORCE_params':
            # Assign default values to FORCE_params
            if 'growthMedium_flux_bounds' not in attr_value.keys():
                attr_value['growthMedium_flux_bounds'] = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}
            if 'stopWith_product_yield_percent' not in attr_value.keys():
                attr_value['stopWith_product_yield_percent'] = 80
            if 'min_biomass_percent' not in attr_value.keys():
                attr_value['min_biomass_percent'] = 10
            if 'blocked_rxns' not in attr_value.keys():
                attr_value['blocked_rxns'] = []
            if 'MUST_L' not in attr_value.keys():
                attr_value['MUST_L'] = []
            if 'MUST_U' not in attr_value.keys():
                attr_value['MUST_U'] = []
            if 'MUST_LU_L' not in attr_value.keys():
                attr_value['MUST_LU_L'] = []
            if 'MUST_LU_U' not in attr_value.keys():
                attr_value['MUST_LU_U'] = []
            if 'MUST_LL_L1' not in attr_value.keys():
                attr_value['MUST_LL_L1'] = []
            if 'MUST_LL_L2' not in attr_value.keys():
                attr_value['MUST_LL_L2'] = []
            if 'MUST_UU_U1' not in attr_value.keys():
                attr_value['MUST_UU_U1'] = []
            if 'MUST_UU_U2' not in attr_value.keys():
                attr_value['MUST_UU_U2'] = []
            if 'total_interven_num' not in attr_value.keys():
                attr_value['total_interven_num'] = 10
            if 'notMUST_total_interven_num' not in attr_value.keys():
                attr_value['notMUST_total_interven_num'] = 0
            if 'notXrxns_interven_num' not in attr_value.keys():
                attr_value['notXrxns_interven_num'] = 0
            if 'notLrxns_interven_num' not in attr_value.keys():
                attr_value['notLrxns_interven_num'] = 0
            if 'notUrxns_interven_num' not in attr_value.keys():
                attr_value['notUrxns_interven_num'] = 0
            if 'fixed_X_rxns' not in attr_value.keys():
                attr_value['fixed_X_rxns'] = []
            if 'fixed_L_rxns' not in attr_value.keys():
                attr_value['fixed_L_rxns'] = []
            if 'fixed_U_rxns' not in attr_value.keys():
                attr_value['fixed_U_rxns'] = []
            if 'ignored_X_rxns' not in attr_value.keys():
                attr_value['ignored_X_rxns'] = []
            if 'ignored_L_rxns' not in attr_value.keys():
                attr_value['ignored_L_rxns'] = []
            if 'ignored_U_rxns' not in attr_value.keys():
                attr_value['ignored_U_rxns'] = []
            if 'read_inSilico_essential_rxns_fromThisFile' not in attr_value.keys():
                attr_value['read_inSilico_essential_rxns_fromThisFile'] = ''
            if 'read_inVivo_essential_rxns_fromThisFile' not in attr_value.keys():
                attr_value['read_inVivo_essential_rxns_fromThisFile'] = ''
            if 'read_blocked_rxns_fromThisFile' not in attr_value.keys():
                attr_value['read_blocked_rxns_fromThisFile'] = ''
            if 'build_new_optModel' not in attr_value.keys():
                attr_value['build_new_optModel'] = True
            if 'dual_formulation_type' not in attr_value.keys():
                attr_value['dual_formulation_type'] = 'simplified'
            if 'taget_product_yield_percent' not in attr_value.keys():
                attr_value['taget_product_yield_percent'] = 80
            if 'validate_results' not in attr_value.keys():
                attr_value['validate_results'] = False
            if 'optimization_solver' not in attr_value.keys():
                attr_value['optimization_solver'] = default_optim_solver
            if 'warnings' not in attr_value.keys():
                attr_value['warnings'] = True
            if 'stdout_msgs' not in attr_value.keys():
                attr_value['stdout_msgs'] = True

        # Simulation conditions name
        if attr_name == 'simulation_condition' and (attr_value is not None and not isinstance(attr_value,str)):
            raise userError('Invalid simulation_condition for fba model. simulation_condition must be a string')

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("'stdout_msgs' must be either True or False")
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("'warnings' must be either True or False")

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

    def find_ref_flux_bounds(self):
        """
        Finds the flux bounds for the reference strain
        """
        # Growth medium
        set_specific_bounds(model = self.model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)
        # Experimental flux data
        set_specific_bounds(model = self.model, flux_bounds = self.exp_flux_bounds_ref, reset_flux_bounds = False)
        # The following is to make sure that 100 moles of glucose is taken up (was used in iAF12600_BoundsWT.gms for IPP)
        set_specific_bounds(model = self.model, flux_bounds = {'GLCptspp':[100,100]}, reset_flux_bounds = False)

        #--- Find max biomass after imposing experimental flux data ---
        for rxn in self.model.reactions:
            rxn.objective_coefficient = 0
        self.model.biomass_reaction.objective_coefficient = 1
        self.model.fba(stdout_msgs = False)
        max_biomass_exp = self.model.fba_model.solution['objective_value'] 
        if self.stdout_msgs:
            print 'max biomass in the presence of experimental flux data = {}\n'.format(max_biomass_exp)
        self.model.biomass_reaction.flux_bounds = [max_biomass_exp, max_biomass_exp]

        # Perform fba again to make sure that the problem is not infeasible when imposing max_biomass_exp
        self.model.fba(stdout_msgs = False)
        if self.model.fba_model.solution['exit_flag'] != 'globallyOptimal': 
            if self.warnings:
                print '\nWARNING (OptForce.py)! Infeasible fba problem after imposing the max biomass flux in the presence of experimental flux data. The LB on biomass reaction flux was set to 99% of max_biomass_exp'
            self.model.biomass_reaction.flux_bounds = [0.99*max_biomass_exp, max_biomass_exp]

            self.model.fba(stdout_msgs = False)
            if self.model.fba_model.solution['exit_flag'] != 'globallyOptimal': 
                raise userError('Infeasible fba problem after even after setting the LB on biomass flux to 99% of the max biomass flux in the presence of experimental flux data')  

        #--- Perform fva  ---
        if self.stdout_msgs:
            print 'OptForce: Finding flux bounds in the reference strain ...',
            start_refBnd_pt = time.clock()
            start_refBnd_wt = time.time()

        self.flux_bounds_ref = fva(model = self.model, save_to_model = False, results_filename = self.save_flux_bounds_ref_toThisFile, warmstart = True, optimization_solver = 'gurobi_ampl', stdout_msgs = False) 

        if self.stdout_msgs:
            elapsed_refBnd_pt = str(timedelta(seconds = time.clock() - start_refBnd_pt))
            elapsed_refBnd_wt = str(timedelta(seconds = time.time() - start_refBnd_wt))
            print ' took {} of processing time and {} of wall time'.format(elapsed_refBnd_pt, elapsed_refBnd_wt)

    def find_overprod_flux_bounds(self):
        """
        Finds the flux bounds for the overproducing strains 
        """
        if self.stdout_msgs:
            print 'OptForce: Finding flux bounds in the overproducing strain ...',
            start_opsBnd_pt = time.clock()
            start_opsBnd_wt = time.time()

        # Growth medium
        set_specific_bounds(model = self.model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)
        # Impost constraints on product yield and min biomass formation flux
        set_specific_bounds(model = self.model, flux_bounds = {self.biomass_rxn_id:[(self.min_biomass_percent/100)*self.max_biomass_flux,None], self.product_exchrxn_id:[(self.product_targetYield_percent/100)*self.product_max_theor_yield,None]}, reset_flux_bounds = False)
    
        self.flux_bounds_overprod = fva(model = self.model, save_to_model = False, results_filename = self.save_flux_bounds_overprod_toThisFile, warmstart = True, optimization_solver = 'gurobi_ampl', stdout_msgs = False) 

        if self.stdout_msgs:
            elapsed_opsBnd_pt = str(timedelta(seconds = time.clock() - start_opsBnd_pt))
            elapsed_opsBnd_wt = str(timedelta(seconds = time.time() - start_opsBnd_wt))
            print ' took {} of processing time and {} of wall time'.format(elapsed_opsBnd_pt, elapsed_opsBnd_wt)

    def find_MUST_singles(self):
        """
        Finds MUST single sets (MUST_L, MUST_U and MUST_X) 
        """
        if self.stdout_msgs:
            print 'OptForce: Finding MUST singles ...',
            start_mustSingl_pt = time.clock()
            start_mustSingl_wt = time.time()

        # A threshold below which the flux value is considered to be practically zero
        zero_flux_thr = 1e-4

        self.MUST_X, self.MUST_L, self.MUST_U = [], [], []
        for rxn in [r for r in self.model.reactions if r.reversibility.lower() != 'exchange']:
            # MUST X: abs(LB_overprod) < zero_flux_thr and abs(UB_overprod) < zero_flux_thr and (abs(LB_ref) > zero_flux_thr or abs(UB_ref) > zero_flux_thr) 
            if (self.flux_bounds_overprod[rxn.id][0] != None and self.flux_bounds_overprod[rxn.id][1] != None and abs(self.flux_bounds_overprod[rxn.id][0]) < zero_flux_thr and abs(self.flux_bounds_overprod[rxn.id][1] < zero_flux_thr)) and (abs(self.flux_bounds_ref[rxn.id][0]) > zero_flux_thr or abs(self.flux_bounds_ref[rxn.id][1]) > zero_flux_thr):
                self.MUST_X.append(rxn.id)

            # MUST_U: UB_ref < LB_oveprod
            elif self.flux_bounds_ref[rxn.id][1] != None and self.flux_bounds_overprod[rxn.id][0] != None and (self.flux_bounds_overprod[rxn.id][0] - self.flux_bounds_ref[rxn.id][1]) >= self.MUST_singles_diff_thr:
                self.MUST_U.append(rxn.id)

            # MUST_L: UB_overprod < LB_ref
            elif self.flux_bounds_overprod[rxn.id][1] != None and self.flux_bounds_ref[rxn.id][0] != None and (self.flux_bounds_ref[rxn.id][0] - self.flux_bounds_overprod[rxn.id][1]) >= self.MUST_singles_diff_thr:
                self.MUST_L.append(rxn.id)

        if self.stdout_msgs:
            elapsed_mustSingl_pt = str(timedelta(seconds = time.clock() - start_mustSingl_pt))
            elapsed_mustSingl_wt = str(timedelta(seconds = time.time() - start_mustSingl_wt))
            print ' found {} MUST_Ls, {} MUST_Us and {} MUST_Xs, took {} of processing time and {} of wall time'.format(len(self.MUST_L), len(self.MUST_U), len(self.MUST_X), elapsed_mustSingl_pt, elapsed_mustSingl_wt)

        #-- Save the MUST single results into a file --
        if self.save_MUST_singles_toThisFile != '':
            with open(self.save_MUST_singles_toThisFile,'w') as f:
                f.write('MUST_X_details = {\n')
                for r_id in self.MUST_X:
                    rxn = self.model.reactions_by_id[r_id]
                    if "'" not in rxn.name:
                        f.write("'{}':{{'name':'{}', 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(rxn.id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[rxn.id],self.flux_bounds_overprod[rxn.id]))               
                    else:
                        f.write("'{}':{{'name':\"{}\", 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(rxn.id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[rxn.id],self.flux_bounds_overprod[rxn.id]))               
                f.write('}\n')

                f.write('\nMUST_L_details = {\n')
                for r_id in self.MUST_L:
                    rxn = self.model.reactions_by_id[r_id]
                    if "'" not in rxn.name:
                        f.write("'{}':{{'name':'{}', 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(rxn.id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[rxn.id],self.flux_bounds_overprod[rxn.id]))               
                    else:
                        f.write("'{}':{{'name':\"{}\", 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(rxn.id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[rxn.id],self.flux_bounds_overprod[rxn.id]))               
                f.write('}\n')

                f.write('\nMUST_U_details = {\n')
                for r_id in self.MUST_U:
                    rxn = self.model.reactions_by_id[r_id]
                    if "'" not in rxn.name:
                        f.write("'{}':{{'name':'{}', 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(rxn.id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[rxn.id],self.flux_bounds_overprod[rxn.id]))               
                    else:
                        f.write("'{}':{{'name':\"{}\", 'gpr':'{}', 'ref_bounds':{}, 'overprod_bounds':{}}},\n".format(rxn.id,rxn.name, rxn.gene_reaction_rule, self.flux_bounds_ref[rxn.id],self.flux_bounds_overprod[rxn.id]))               
                f.write('}\n')

    def find_MUST_doubles(self):
        """
        Finds MUST double sets (MUST_LU, MUST_UU and MUST_LL) 
        """
        # MUST_LU
        if self.stdout_msgs:
            print 'OptForce: Finding MUST LUs ...',
            start_mustlu_pt = time.clock()
            start_mustlu_wt = time.time()

        must_LU_inst = MUST_doubles(model = self.model, product_exchrxn_id = self.product_exchrxn_id, flux_bounds_ref = self.flux_bounds_ref, flux_bounds_overprod = self.flux_bounds_overprod, growthMedium_flux_bounds = self.growthMedium_flux_bounds, product_targetYield_percent = self.product_targetYield_percent, min_biomass_percent = self.min_biomass_percent, blocked_rxns = self.blocked_rxns, MUST_X = self.MUST_X, MUST_L = self.MUST_L, MUST_U = self.MUST_U, MUST_double_type = 'MUST_LU', build_new_optModel = self.MUST_doubles_params['build_new_optModel'], objective_thr = self.MUST_doubles_params['objective_thr'], validate_results = self.MUST_doubles_params['validate_results'], results_filename = self.MUST_doubles_params['results_filename'], create_new_results_file = True, stdout_msgs = self.MUST_doubles_params['stdout_msgs'])
        (solutions, termination_condition) = must_LU_inst.run()
        self.MUST_LU_L = []
        self.MUST_LU_U = []
        for soln in solutions:
            self.MUST_LU_L.append(soln['L'])
            self.MUST_LU_U.append(soln['U'])

        if self.stdout_msgs:
            elapsed_mustlu_pt = str(timedelta(seconds = time.clock() - start_mustlu_pt))
            elapsed_mustlu_wt = str(timedelta(seconds = time.time() - start_mustlu_wt))
            print '\t{} MUST_LUs found, ended with "{}" and took {} of processing time and {} of wall time'.format(len(solutions), termination_condition, elapsed_mustlu_pt, elapsed_mustlu_wt)

        # MUST_UU
        if self.stdout_msgs:
            print 'OptForce: Finding MUST UUs ...',
            start_mustuu_pt = time.clock()
            start_mustuu_wt = time.time()

        must_UU_inst = MUST_doubles(model = self.model, product_exchrxn_id = self.product_exchrxn_id, flux_bounds_ref = self.flux_bounds_ref, flux_bounds_overprod = self.flux_bounds_overprod, growthMedium_flux_bounds = self.growthMedium_flux_bounds, product_targetYield_percent = self.product_targetYield_percent, min_biomass_percent = self.min_biomass_percent, blocked_rxns = self.blocked_rxns, MUST_X = self.MUST_X, MUST_L = self.MUST_L, MUST_U = self.MUST_U, MUST_double_type = 'MUST_UU', build_new_optModel = self.MUST_doubles_params['build_new_optModel'], objective_thr = self.MUST_doubles_params['objective_thr'], validate_results = self.MUST_doubles_params['validate_results'], results_filename = self.MUST_doubles_params['results_filename'], create_new_results_file = False, stdout_msgs = self.MUST_doubles_params['stdout_msgs'])
        (solutions, termination_condition) = must_UU_inst.run()
        self.MUST_UU_U1 = []
        self.MUST_UU_U2 = []
        for soln in solutions:
            self.MUST_UU_U1.append(soln['U1'])
            self.MUST_UU_U2.append(soln['U2'])

        if self.stdout_msgs:
            elapsed_mustuu_pt = str(timedelta(seconds = time.clock() - start_mustuu_pt))
            elapsed_mustuu_wt = str(timedelta(seconds = time.time() - start_mustuu_wt))
            print '\t{} MUST_UUs found, ended with "{}" and took {} of processing time and {} of wall time'.format(len(solutions), termination_condition, elapsed_mustuu_pt, elapsed_mustuu_wt)

        # MUST_LL
        if self.stdout_msgs:
            print 'OptForce: Finding MUST LLs ...',
            start_mustll_pt = time.clock()
            start_mustll_wt = time.time()

        must_LL_inst = MUST_doubles(model = self.model, product_exchrxn_id = self.product_exchrxn_id, flux_bounds_ref = self.flux_bounds_ref, flux_bounds_overprod = self.flux_bounds_overprod, growthMedium_flux_bounds = self.growthMedium_flux_bounds, product_targetYield_percent = self.product_targetYield_percent, min_biomass_percent = self.min_biomass_percent, blocked_rxns = self.blocked_rxns, MUST_X = self.MUST_X, MUST_L = self.MUST_L, MUST_U = self.MUST_U, MUST_double_type = 'MUST_LL', build_new_optModel = self.MUST_doubles_params['build_new_optModel'], objective_thr = self.MUST_doubles_params['objective_thr'], validate_results = self.MUST_doubles_params['validate_results'], results_filename = self.MUST_doubles_params['results_filename'], create_new_results_file = False, stdout_msgs = self.MUST_doubles_params['stdout_msgs'])
        (solutions, termination_condition) = must_LL_inst.run()
        self.MUST_LL_L1 = []
        self.MUST_LL_L2 = []
        for soln in solutions:
            self.MUST_LL_L1.append(soln['L1'])
            self.MUST_LL_L2.append(soln['L2'])

        if self.stdout_msgs:
            elapsed_mustll_pt = str(timedelta(seconds = time.clock() - start_mustll_pt))
            elapsed_mustll_wt = str(timedelta(seconds = time.time() - start_mustll_wt))
            print '\t{} MUST_LLs found, ended with "{}" and took {} of processing time and {} of wall time'.format(len(solutions), termination_condition, elapsed_mustll_pt, elapsed_mustll_wt)

    def run(self):
        """
        Runs the entire OptForce procedure 
        """
        if self.stdout_msgs:
            print '\nStart running OptForce ...'

        # Find the maximum biomass flux
        self.find_maxBiomass_flux() 

        # Find the maximum theoretical yield of the product
        self.find_prod_max_theor_yield()

        # Find flux bounds in the reference strain 
        if self.read_flux_bounds_ref_fromThisFile == '':
            self.find_ref_flux_bounds() 
        else:
            load_source('dataFile',self.read_flux_bounds_ref_fromThisFile)
            import dataFile
            self.flux_bounds_ref = dataFile.fva_flux_bounds

        # Find flux bounds in the overproducing strain 
        if self.read_flux_bounds_overprod_fromThisFile == '':
            self.find_overprod_flux_bounds() 
        else:
            load_source('dataFile',self.read_flux_bounds_overprod_fromThisFile)
            import dataFile
            self.flux_bounds_overprod = dataFile.fva_flux_bounds

        # Find MUST singles (MUST_L, MUST_U and MUST_X)
        if self.read_MUST_singles_fromThisFile == '':
            self.find_MUST_singles()
        else:
            load_source('dataFile',self.read_MUST_singles_fromThisFile)
            import dataFile
            self.MUST_X = dataFile.MUST_X_details.keys()
            self.MUST_L = dataFile.MUST_L_details.keys()
            self.MUST_U = dataFile.MUST_U_details.keys()
            if self.stdout_msgs:
                print 'OptForce: {} MUST_Ls, {} MUST_Us and {} MUST_Xs were read from the input file'.format(len(self.MUST_L), len(self.MUST_U), len(self.MUST_X))

        # Find MUST doubles (MUST_LU, MUST_UU and MUST_LU)
        if self.read_MUST_doubles_fromThisFile == '':
            self.find_MUST_doubles()
        else:
            load_source('dataFile',self.read_MUST_doubles_fromThisFile)
            import dataFile
            # MUST_LU
            MUST_LU_details = dataFile.MUST_LU_details
            self.MUST_LU_L = []
            self.MUST_LU_U = []
            for MUST_LU in MUST_LU_details:
                self.MUST_LU_L.append(MUST_LU['L'])
                self.MUST_LU_U.append(MUST_LU['U'])

            # MUST_UU
            MUST_UU_details = dataFile.MUST_UU_details
            self.MUST_UU_U1 = []
            self.MUST_UU_U2 = []
            for MUST_UU in MUST_UU_details:
                self.MUST_UU_U1.append(MUST_UU['U1'])
                self.MUST_UU_U2.append(MUST_UU['U2'])

            # MUST_LL
            MUST_LL_details = dataFile.MUST_LL_details
            self.MUST_LL_L1 = []
            self.MUST_LL_L2 = []
            for MUST_LL in MUST_LL_details:
                self.MUST_LL_L1.append(MUST_LL['L1'])
                self.MUST_LL_L2.append(MUST_LL['L2'])

            if self.stdout_msgs:
                print 'OptForce: {} MUST_LUs, {} MUST_UUs and {} MUST_LLs were read from the input file'.format(len(MUST_LU_details), len(MUST_UU_details), len(MUST_LL_details))


        # Find FORCE sets
        if self.stdout_msgs:
            print '\nOptForce: Finding FORCE sets ...'
        force_inst = FORCE(model = self.model, product_exchrxn_id = self.product_exchrxn_id, flux_bounds_ref = self.flux_bounds_ref, 
                           flux_bounds_overprod = self.flux_bounds_overprod, growthMedium_flux_bounds = self.growthMedium_flux_bounds, min_biomass_percent = self.min_biomass_percent, stopWith_product_yield_percent = self.FORCE_params['stopWith_product_yield_percent'], 
                           MUST_X = self.MUST_X, MUST_L = self.MUST_L, MUST_U = self.MUST_U, MUST_LU_L = self.MUST_LU_L, MUST_LU_U = self.MUST_LU_U, MUST_LL_L1 = self.MUST_LL_L1, MUST_LL_L2 = self.MUST_LL_L2, MUST_UU_U1 = self.MUST_UU_U1, MUST_UU_U2 = self.MUST_UU_U2, 
                           total_interven_num = self.FORCE_params['total_interven_num'], notMUST_total_interven_num = self.FORCE_params['notMUST_total_interven_num'], notXrxns_interven_num = self.FORCE_params['notXrxns_interven_num'], notLrxns_interven_num = self.FORCE_params['notLrxns_interven_num'], notUrxns_interven_num = self.FORCE_params['notUrxns_interven_num'], 
                           fixed_X_rxns = self.FORCE_params['fixed_X_rxns'], fixed_L_rxns = self.FORCE_params['fixed_L_rxns'], fixed_U_rxns = self.FORCE_params['fixed_U_rxns'], 
                           ignored_X_rxns = self.FORCE_params['ignored_X_rxns'], ignored_L_rxns = self.FORCE_params['ignored_L_rxns'], ignored_U_rxns = self.FORCE_params['ignored_U_rxns'], 
                           inSilico_essential_rxns = self.inSilico_essential_rxns, inVivo_essential_rxns = self.inVivo_essential_rxns, blocked_rxns = self.blocked_rxns, 
                           build_new_optModel = self.FORCE_params['build_new_optModel'], dual_formulation_type = self.FORCE_params['dual_formulation_type'], validate_results = self.FORCE_params['validate_results'], 
                           results_filename = self.FORCE_params['results_filename'], optimization_solver = self.FORCE_params['optimization_solver'], warnings = self.FORCE_params['warnings'], stdout_msgs = self.FORCE_params['stdout_msgs'])
        #force_inst.run_primal_dual()
        force_inst.run()
 
