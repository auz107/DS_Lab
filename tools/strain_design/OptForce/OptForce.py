from __future__ import division
import sys,os, time
sys.path.append('../../../')
from copy import deepcopy
from tools.userError import userError
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
from coopr.pyomo import *
from coopr.opt import *
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class OptForce(object):
    """
    Checks whether a given flux data leads to a feasible FBA for a given metabolic model
    """
    def __init__(self, model, product_exchrxn_id, exp_flux_bounds, product_targetYield_percent = 80, min_biomass_percent = 10, growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, flux_bounds_filename = '', MUST_singles_filename = '', MUST_doubles_filename = '', results_filename_base = '', optimization_solver = 'gurobi', simulation_conditions = '', stdout_msgs = True, warnings = True):
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
                                 flux_bounds_filename: Name of and path of the file containing the flux bounds
                                                      file for exchange reactions corresponding to compounds in 
                                                       in the growth media
                                 flux_bounds_dict: A dictionary containing the flux bounds for other reactions (such as carbon
                                                  source uptake, oxygen uptake, etc), that have to be set for reactions not 
                                                   in flux_data_filenam or media_filename
            flux_bounds_filename: A string containing the filename (including path) storinring the 
                                 flux bounds for the reference and overproducing strains (see function
                                 find_ref_overprod_flux_bounds for how the flux bounds are stored in
                                 files).
          MUST_singles_filename: A string containing the filename (including path) storing the 
                                 information for reactons appearing in MUST_X, MUST_L and MUST_U sets.
                                 (see function find_MUST_singles for how MUST single are stored in files)
           MUST_double_filename: A string containing the filename (including path) storing the 
                                 information for reactons appearing in MUST_LL, MUST_UU and MUST_LU 
                                 sets (see function find_MUST_soubles for how MUST single are stored 
                                 in files)
            optimization_solver: A string containing showing which optimization solver to use
                                 Allowed choices are 'gurobi' and 'cplex'
          simulation_conditions: A string representing the simulation conditon 
          results_filename_base: The base name of the results file including it path. Various results are written into files 
                                 results_filename_base + '_' + [specific results filename]
                       warnings: Can be True or False shwoing whether the warnings should be 
                                 writtten to the 
                                 screen or not. The default is True  
                    stdout_msgs: Prints a summary of the FBA and a number of other messages 
                                 in the output if 'on.
       
        Ali R. Zomorrodi - Segre's Lab @ BU
        Last updated: 03-16-2016
        """
        self.shared_data = shared_data_holder(model = model, product_exchrxn_id = product_exchrxn_id, exp_flux_bounds = exp_flux_bounds, product_targetYield_percent = product_targetYield_percent, min_biomass_percent = min_biomass_percent, growthMedium_flux_bounds = growthMedium_flux_bounds, flux_bounds_filename =  flux_bounds_filenam, MUST_singles_filename =  MUST_singles_filenam, MUST_doubles_filename =  MUST_doubles_filenam, results_filename_base =  results_filename_bas, optimization_solver =  optimization_solve, simulation_conditions = ,imulation_condition, warnings = warnings, stdout_msgs =  stdout_msg):

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
         attr_name: Attribute name
        attr_value: Attribute value
        """
        if attr_name == 'shared_data' and not isinstance(attr_value,shared_data_holder):
            raise TypeError('shared_data must be instance of class shared_data_holder')

        self.__dict__[attr_name] = attr_value

    def find_maxBiomass_flux(self):
        """
        Finds the maximum thoeretical biomass flux
        """
        for rxn in self.shared_data.model.reactions:
            rxn.objective_coefficient = 0
        self.shared_data.model.biomass_reaction.objective_coefficient = 1

        set_specific_bounds(model = self.shared_data.model, file_name = self.shared_data.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.shared_data.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)

        self.shared_data.model.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False, warnings = True)
        if self.shared_data.model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            self.shared_data.max_biomass_flux = self.shared_data.model.fba_model.solution['objective_value']
            if self.shared_data.stdout_msgs:
                print 'OptForce: Theoretical maximum biomass formation flux = {}'.format(self.shared_data.max_biomass_flux)
        else:    
            raise userError('The FBA problem for finding the maximum theoretical biomass flux was not solved to optimality!')

    def find_prod_max_theor_yield(self):
        """
        Finds the maximum thoeretical yield for the product of interest (this is the max flux of the exchange reactions for the product)
        """
        for rxn in self.shared_data.model.reactions:
            rxn.objective_coefficient = 0
        self.shared_data.model.reactions_by_id[self.shared_data.product_exchrxn_id].objective_coefficient = 1

        set_specific_bounds(model = self.shared_data.model, file_name = self.shared_data.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.shared_data.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)

        self.shared_data.model.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False, warnings = True)
        if self.shared_data.model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            self.shared_data.product_max_theor_yield = self.shared_data.model.fba_model.solution['objective_value']
            if self.shared_data.stdout_msgs:
                print 'OptForce: Theoretical maximum product yield = {}'.format(self.shared_data.product_max_theor_yield)
        else:    
            raise userError('The FBA problem for finding the theoretical maximum product yield was not solved to optimality!')

    def find_ref_overprod_flux_bounds(self):
        """
        Finds the flux bounds for the reference and overproducing strains 
        """
        #-- Finds flux bounds for the reference strain -- 
        if self.shared_data.stdout_msgs:
            print 'OptForce: Finding flux bounds in the reference strain ...'

        # Growth medium
        set_specific_bounds(model = self.shared_data.model, file_name = self.shared_data.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.shared_data.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)
        # Experimental flux data
        set_specific_bounds(model = self.shared_data.model, flux_bounds = self.shared_data.exp_flux_bounds, reset_flux_bounds = False)
    
        self.shared_data.flux_bounds_ref = fva(model = self.shared_data.model, store_fva_flux_bounds = False, stdout_msgs = False) 

        #-- Finds flux bounds for the overproducing strain -- 
        if self.shared_data.stdout_msgs:
            print 'OptForce: Finding flux bounds in the overproducing strain ...'

        # Growth medium
        set_specific_bounds(model = self.shared_data.model, file_name = self.shared_data.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.shared_data.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)
        # Impost constraints on product yield and min biomass formation flux
        set_specific_bounds(model = self.shared_data.model, flux_bounds = {self.shared_data.model.biomass_reaction.id:[(self.shared_data.min_biomass_percent/100)*self.shared_data.max_biomass_flux,None], self.shared_data.product_exchrxn_id:[(self.shared_data.product_targetYield_percent/100)*self.shared_data.product_max_theor_yield,None]}, reset_flux_bounds = False)
    
        self.shared_data.flux_bounds_overprod = fva(model = self.shared_data.model, store_fva_flux_bounds = False, stdout_msgs = False) 

        #-- Save the flux bounds results into a file --
        if self.shared_data.results_filename_base != '':
            with open(self.shared_data.results_filename_base + '_flux_bounds.py','w') as f:
                f.write('flux_bounds_ref = {\n')
                for rxn in self.shared_data.model.reactions:
                    f.write("'{}':{},\n".format(rxn.id,self.shared_data.flux_bounds_ref[rxn.id]))
                f.write('}\n')

                f.write('\nflux_bounds_overprod = {\n')
                for rxn in self.shared_data.model.reactions:
                    f.write("'{}':{},\n".format(rxn.id,self.shared_data.flux_bounds_overprod[rxn.id]))
                f.write('}\n')

        if self.shared_data.stdout_msgs:
            print '\nFlux bounds for the reference and overproducing strains were written into {}\n'.format(self.shared_data.results_filename_base + '_flux_bounds.py')

    def find_MUST_singles(self):
        """
        Finds MUST single sets (MUST_L, MUST_U and MUST_X) 
        """
        if self.shared_data.stdout_msgs:
            print 'OptForce: Finding MUST singles ...\n'

        # A threshold below which the flux value is considered to be practically zero
        zero_flux_thr = 1e-4

        self.shared_data.MUST_X, self.shared_data.MUST_L, self.shared_data.MUST_U = [], [], []
        for rxn in [r for r in self.shared_data.model.reactions if r.reversibility.lower() != 'exchange']:
            # MUST X: abs(LB_overprod) < zero_flux_thr and abs(UB_overprod) < zero_flux_thr and (abs(LB_ref) > zero_flux_thr or abs(UB_ref) > zero_flux_thr) 
            if (self.shared_data.flux_bounds_overprod[rxn.id][0] != None and self.shared_data.flux_bounds_overprod[rxn.id][1] != None and abs(self.shared_data.flux_bounds_overprod[rxn.id][0]) < zero_flux_thr and abs(self.shared_data.flux_bounds_overprod[rxn.id][1] < zero_flux_thr)) and (abs(self.shared_data.flux_bounds_ref[rxn.id][0]) > zero_flux_thr or abs(self.shared_data.flux_bounds_ref[rxn.id][1]) > zero_flux_thr):
                self.shared_data.MUST_X.append(rxn)

            # MUST_U: UB_ref < LB_oveprod
            elif self.shared_data.flux_bounds_ref[rxn.id][1] != None and self.shared_data.flux_bounds_overprod[rxn.id][0] != None and self.shared_data.flux_bounds_ref[rxn.id][1] < self.shared_data.flux_bounds_overprod[rxn.id][0]:
                self.shared_data.MUST_U.append(rxn)

            # MUST_L: UB_overprod < LB_ref
            elif self.shared_data.flux_bounds_overprod[rxn.id][1] != None and self.shared_data.flux_bounds_ref[rxn.id][0] != None and self.shared_data.flux_bounds_overprod[rxn.id][1] < self.shared_data.flux_bounds_ref[rxn.id][0]:
                self.shared_data.MUST_L.append(rxn)

        #-- Save the MUST single results into a file --
        if self.shared_data.results_filename_base != '':
            with open(self.shared_data.results_filename_base + '_MUST_singles.py','w') as f:
                f.write('MUST_X_details = {\n')
                for rx in self.shared_data.MUST_X:
                    f.write("'{}':{{'name':{}, 'gpr':{}, 'ref_bounds':{}, 'overprod_bounds':{}}},\n':".format(rx.id,rx.name, rx.gene_reaction_rule, self.shared_data.flux_bounds_ref[rx.id],self.shared_data.flux_bounds_overprod[rx.id]))               
                f.write('}\n')

                f.write('\nMUST_L_details = {\n')
                for rl in self.shared_data.MUST_L:
                    f.write("'{}':{{'name':{}, 'gpr':{}, 'ref_bounds':{}, 'overprod_bounds':{}}},\n':".format(rl.id,rl.name, rl.gene_reaction_rule, self.shared_data.flux_bounds_ref[rl.id],self.shared_data.flux_bounds_overprod[rl.id]))               
                f.write('}\n')

                f.write('\nMUST_U_details = {\n')
                for ru in self.shared_data.MUST_U:
                    f.write("'{}':{{'name':{}, 'gpr':{}, 'ref_bounds':{}, 'overprod_bounds':{}}},\n':".format(ru.id,ru.name, ru.gene_reaction_rule, self.shared_data.flux_bounds_ref[ru.id],self.shared_data.flux_bounds_overprod[ru.id]))               
                f.write('}\n')

        if self.shared_data.stdout_msgs:
            print '\nMUST single reactions were written into {}\n'.format(self.shared_data.results_filename_base + '_MUST_singles.py')

    def run(self):
        """
        Runs the entire OptForce procedure 
        """
        if self.shared_data.stdout_msgs:
            print '\nStart running OptForce ...'

        # Find the maximum biomass flux
        self.find_maxBiomass_flux() 

        # Find the maximum theoretical yield of the product
        self.find_prod_max_theor_yield()

        # Find flux bounds in the wild-type and overproducing strains
        if self.flux_bounds_filename == '':
            self.find_ref_overprod_flux_bounds() 
        else:
            load_source('dataFile',self.shared_data.flux_bounds_filename)
            import dataFile
            self.shared_data.flux_bounds_ref = dataFile.flux_bounds_ref
            self.shared_data.flux_bounds_overprod = dataFile.flux_bounds_overprod

        # Find MUST singles (MUSTL_L, MUST_U and MUST_X)
        if self.shared_data.MUST_singles_filename == '':
            self.find_MUST_singles()
        else:
            load_source('dataFile',self.shared_data.MUST_singles_filename)
            import dataFile
            self.shared_data.MUST_X = dataFile.MUST_X_details.keys()
            self.shared_data.MUST_L = dataFile.MUST_L_details.keys()
            self.shared_data.MUST_U = dataFile.MUST_U_details.keys()

class MUST_doubles(object):
    """
    This class contains defines the optimization problem to identify the MUST double sets 
    of reactions, i.e., MUST_LU (or MUST_UL), MUST_LL and MUST_UU

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: March 19, 2016
    """

    def __init__(self,shared_data, MUST_double_type): 
        """
        All inputs are a subset of those for OptForce

        INPUTS:
        ------
             shared_data: An instance of class shared_data_holder
        MUST_double_type: Tyype of MUST doulbe to find. Allowed choices include: MUST_LU 
        (or MUST_UL), MUST_LL and MUST_UU
        """
        # Shared_data 
        self.shared_data = shared_data

        # Type of MUST doubles to identify
        self.MUST_double_type = MUST_double_type

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
         attr_name: Attribute name
        attr_value: Attribute value
        """
        if attr_name == 'shared_data' and not isinstance(attr_value,shared_data_holder):
            raise TypeError('shared_data must be instance of class shared_data_holder')

        if attr_name == 'MUST_double_type' and not isinstance(attr_value,str):
            raise TypeError('MUST_double_type must be a string')
        elif attr_name == 'MUST_double_type' and attr_value.lower() not in ['must_lu','must_ul','must_uu','must_ll']: 
            raise TypeError('Invalid MUST_double_type value! Allowed choices are MUST_LU, MUST_UL, MUST_LL, MUST_UU')


        self.__dict__[attr_name] = attr_value


class shared_data(object):
    """
    A class storing the data that should be shared among the main OptForce class and the
    classes finding MUST doubles and FORCE sets of reactions. This class just facilitates
    data sharing among the three aformentioned classes.
    """
    def __init__(self, model, product_exchrxn_id, exp_flux_bounds, product_targetYield_percent = 80, min_biomass_percent = 10, growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, flux_bounds_filename = '', MUST_singles_filename = '', MUST_doubles_filename = '', results_filename_base = '', optimization_solver = 'gurobi', simulation_conditions = '', stdout_msgs = True, warnings = True):

        # Metabolic model
        self.model = model

        # Id of the exchange reaction for the product of interest
        self.product_exchrxn_id = product_exchrxn_id

        # Target yeild
        self.product_targetYield_percent = product_targetYield_percent

        # Minimum biomass formation
        self.min_biomass_percent = min_biomass_percent

        # Experimental flux data
        self.exp_flux_bounds = exp_flux_bounds

        # Flux bounds for exchange reactions corresponding to compounds present in the growth medium
        self.growthMedium_flux_bounds = growthMedium_flux_bounds

        # Solver name
        self.optimization_solver = optimization_solver

        # Simulation conditions
        self.simulation_conditions = simulation_conditions

        # results_filename_base
        self.results_filename_base = results_filename_base

        # File names holding the flux bounds, MUST single and MUST doulbe sets
        self.flux_bounds_filename = flux_bounds_filename
        self.MUST_singles_filename = MUST_singles_filename
        self.MUST_doubles_filename = MUST_doubles_filename

        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.warnings = warnings

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

        # exp_flux_bounds
        if attr_name == 'exp_flux_bounds' and not isinstance(attr_value,dict):
            raise TypeError('exp_flux_bounds must be a dictionary')

        # growthMedium_flux_bounds 
        if attr_name == 'growthMedium_flux_bounds' and not isinstance(attr_value,dict): 
            raise TypeError('growthMedium_flux_bounds must be a dictionary')
        if attr_name == 'growthMedium_flux_bounds' and len([k for k in attr_value.keys() if k.lower() not in ['flux_bounds_filename','flux_bounds_dict']]) > 0: 
            raise ValueError('Invalid key for growthMedium_flux_bounds. Allowed keys are flux_bounds_filename and flux_bounds_dict')

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Simulation conditions name
        if attr_name == 'simulation_conditions' and (attr_value is not None and not isinstance(attr_value,str)): 
            raise userError('Invalid simulation_conditions for fba model. simulation_conditions must be a string')

        # results_filename_base 
        if attr_name == 'results_filename_base' and not isinstance(attr_value,str): 
            raise TypeError('results_filename_base must be a string')
         
        # flux_bounds_filename, MUST_singles_filename and MUST_doubles_filename
        if attr_name == 'flux_bounds_filename' and not isinstance(attr_value,str): 
            raise TypeError('flux_bounds_filename must be a string')
        if attr_name == 'MUST_singles_filename' and not isinstance(attr_value,str): 
            raise TypeError('MUST_singles_filename must be a string')
        if attr_name == 'MUST_dboules_filename' and not isinstance(attr_value,str): 
            raise TypeError('MUST_doubles_filename must be a string')

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("'stdout_msgs' must be either True or False")
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("'warnings' must be either True or False")

#-------------------------
if __name__ == '__main__':
    pass 
 
