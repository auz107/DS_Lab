from __future__ import division
import sys,os, time
sys.path.append('../../')
from copy import deepcopy
from tools.userError import userError
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.fva import fva
from tools.fba.create_model import create_model
from tools.fba.set_specific_bounds import set_specific_bounds
from imp import load_source
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class check_fluxData(object):
    """
    Checks whether a given flux data leads to a feasible FBA for a given metabolic model
    """

    def __init__(self, model, exp_flux_bounds, optimization_solver = 'gurobi', simulation_conditions = '', stdout_msgs = True, warnings = True):
        """
        INPUTS:
        ------
                        model: An instance of class model
                exp_flux_bounds: A dictionary with keys and values as follows:
                                 keys: Id of reactions with available flux data
                               values: A list of two elements in the form [LB,UB]
                                       containing the LB and UB on reaction fluxes. 
                                       LB and UB can be None as well (if there is no
                                       experimental data for reaction fluxes)
          optimization_solver: A string containing showing which optimization solver to use
                               Allowed choices are 'gurobi' and 'cplex'
        simulation_conditions: A string representing the simulation conditon 
                     warnings: Can be True or False shwoing whether the warnings should be 
                               writtten to the 
                               screen or not. The default is True  
                  stdout_msgs: Prints a summary of the FBA and a number of other messages 
                               in the output if 'on.
       
        Ali R. Zomorrodi - Segre's Lab @ BU
        Last updated: 03-14-2016
        """
        # Metabolic model
        self.model = model

        # Experimental flux data
        self.exp_flux_bounds = exp_flux_bounds

        # Solver name
        self.optimization_solver = optimization_solver

        # Simulation conditions
        self.simulation_conditions = simulation_conditions

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

        # exp_flux_bounds
        if attr_name == 'exp_flux_bounds' and not isinstance(attr_value,dict):
            raise TypeError('exp_flux_bounds must be a dictionary')

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Simulation conditions name
        if attr_name == 'simulation_conditions' and (attr_value is not None and not isinstance(attr_value,str)): 
            raise userError('Invalid simulation_conditions for fba model. simulation_conditions must be a striing')

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("'stdout_msgs' must be either True or False")
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("'warnings' must be either True or False")

        self.__dict__[attr_name] = attr_value

    def check_feasbility(self):
        """
        Check the feasibility of an FBA problem for a given flux data  
        """
        # Assign the objective function coefficients
        for rxn in self.model.reactions:
            rxn.objective_coefficient = 0
        self.model.biomass_reaction.objective_coefficient = 1
    
        # Growth medium
        set_specific_bounds(model = self.model, flux_bounds = self.exp_flux_bounds, reset_flux_bounds = False)
    
        # Perform FBA 
        self.model.fba(assign_wildType_max_biomass = False, stdout_msgs = self.stdout_msgs)

    def fit_max_exp_data(self):
        """
        Solves an MILP to find the maximum data that allow a feasible FBA problem. In other words,
        it finds the experimental LB and UB on fluxes that make the FBA probelm infeasible
        """
        pass

#-------------------------------------------------------------------------------
def Ecoli_minimalAerobic(): 
    use_iAF1260 = False
    use_iJO1366 = True

    # Import the flux bounds in flux_data_filename
    load_source('dataFile','/usr2/postdoc/alizom/work/tools/strain_design/flux_data/Escherichia_coli/iJO1366_minimal_aerobic_glc_PMID_23036703.py')
    import dataFile
    exp_flux_bounds = dataFile.exp_flux_bounds

    if use_iAF1260:
        print '---- iAF1260 model ---'
        # Model path
        model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iAF1260/'
        print '\nFBA before adding flux data: '
        iAF1260 = create_model(organism_info = {'id':'Ecoli', 'name':'Escherichia_coli'}, model_info = {'id':'iAF1260_minimal_aerobic_glc_PMID_23036703', 'sbml_filename':model_path + 'iAF1260_updated.xml', 'biomassrxn_id':'Ec_biomass_iAF1260_core_59p81M'}, growthMedium_fluxBounds = {'media_filename':model_path + 'iAF1260_minimal_glucose_aerobic.py', 'other_flux_bounds': {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}}, stdout_msgs = True, warnings = True)

        print '\nFBA after adding flux data: '
        checkFlux_iAF1260 = check_fluxData(model = iAF1260, exp_flux_bounds = exp_flux_bounds, simulation_conditions = 'minimal aetobic glucose', stdout_msgs = True, warnings = True)
        checkFlux_iAF1260.check_feasbility()

    if use_iJO1366:
        print '---- iJO1366 model ---'
        model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'
        iJO1366 = create_model(organism_info = {'id':'Ecoli', 'name':'Escherichia_coli'}, model_info = {'id':'iJO1366_minimal_aerobic_glc_PMID_23036703', 'sbml_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_fluxBounds = {'fluxBounds_filename':model_path + 'iJO1366_minimal_glucose_aerobic.py', 'fluxBounds_dict': {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}}, stdout_msgs = True, warnings = True)

        print '\nFBA after adding flux data: '
        checkFlux_iJO1366 = check_fluxData(model = iJO1366, exp_flux_bounds = exp_flux_bounds, simulation_conditions = 'minimal aetobic glucose', stdout_msgs = True, warnings = True)
        checkFlux_iJO1366.check_feasbility()

        print '\nFBA after incorporating the biomass flux obtained with using exp flux data on top of the flux data themselves: '
        exp_flux_bounds[iJO1366.biomass_reaction.id] = [5.04,5.04]
        checkFlux_iJO1366 = check_fluxData(model = iJO1366, exp_flux_bounds = exp_flux_bounds, simulation_conditions = 'minimal aetobic glucose', stdout_msgs = True, warnings = True)
        checkFlux_iJO1366.check_feasbility()

        #fva_flux_bounds = fva(model = iJO1366, store_fva_flux_bounds = False, stdout_msgs = True)
        #print '\n\n',fva_flux_bounds
        #with open('res.txt','w') as f:
        #    for r in fva_flux_bounds.keys():
        #        f.write('{}:{}\n'.format(r,fva_flux_bounds[r]))

#-------------------------
if __name__ == '__main__':
    Ecoli_minimalAerobic()
 
