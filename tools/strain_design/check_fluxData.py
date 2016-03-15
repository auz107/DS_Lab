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
def create_model(organism_info = {'id':'', 'name':''}, model_info = {'id':'', 'sbml_filename':None, 'biomassrxn_id':None}, media_info = {'media_filename':None, 'other_flux_bounds': {}}, stdout_msgs = True, warnings = True):
    """
    Creates a metabolic model and assigns the flux bounds for a desired growth media and performs FBA.

    INPUTS:
    ------
     organism_info: A dictionary containing the information about he organism (id and anme)
        model_info: A dictionary containing the metabolic model information with the following keys:
                                'id': model id
                     'sbml_filename': Name of and path of the SBML file containing the model
                     'biomassrxn_id': The if of biomass reaction
        media_info: Information about growth media
                       media_filename: Name of and path of the file containing the flux bounds
                                       file for exchange reactions corresponding to compounds in 
                                       in the growth media
                    other_flux_bounds: A dictionary containing the flux bounds for other reactions (such as carbon
                                       source uptake, oxygen uptake, etc), that have to be set for reactions not 
                                       in flux_data_filenam or media_filename
    """
    # Define the organism
    model_organism = organism(id = organism_info['id'], name = organism_info['name'])

    model = read_sbml_model(file_name = model_info['sbml_filename'], model_id = model_info['id'], model_organism = model_organism, model_type = 'metabolic',import_params = False)

    model.biomass_reaction = model.reactions_by_id[model_info['biomassrxn_id']]

    # Assign the objective function coefficients
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.biomass_reaction.objective_coefficient = 1

    # Growth medium
    set_specific_bounds(model = model, file_name = media_info['media_filename'], flux_bounds = media_info['other_flux_bounds'])

    # Perform FBA for the wild-type
    model.fba(assign_wildType_max_biomass = True, stdout_msgs = stdout_msgs)

    return model

#-------------------------------------------------------------------------------
def Ecoli_minimalAerobic(): 
    use_iAF1260 = True
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
        iAF1260 = create_model(organism_info = {'id':'Ecoli', 'name':'Escherichia_coli'}, model_info = {'id':'iAF1260_minimal_aerobic_glc_PMID_23036703', 'sbml_filename':model_path + 'iAF1260_updated.xml', 'biomassrxn_id':'Ec_biomass_iAF1260_core_59p81M'}, media_info = {'media_filename':model_path + 'iAF1260_minimal_glucose_aerobic.py', 'other_flux_bounds': {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}}, stdout_msgs = True, warnings = True)

        print '\nFBA after adding flux data: '
        checkFlux_iAF1260 = check_fluxData(model = iAF1260, exp_flux_bounds = exp_flux_bounds, simulation_conditions = 'minimal aetobic glucose', stdout_msgs = True, warnings = True)
        checkFlux_iAF1260.check_feasbility()

    if use_iJO1366:
        print '---- iJO1366 model ---'
        model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'
        iJO1366 = create_model(organism_info = {'id':'Ecoli', 'name':'Escherichia_coli'}, model_info = {'id':'iJO1366_minimal_aerobic_glc_PMID_23036703', 'sbml_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, media_info = {'media_filename':model_path + 'iJO1366_minimal_glucose_aerobic.py', 'other_flux_bounds': {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}}, stdout_msgs = True, warnings = True)

        print '\nFBA after adding flux data: '
        checkFlux_iJO1366 = check_fluxData(model = iJO1366, exp_flux_bounds = exp_flux_bounds, simulation_conditions = 'minimal aetobic glucose', stdout_msgs = True, warnings = True)
        checkFlux_iJO1366.check_feasbility()

#-------------------------
if __name__ == '__main__':
    Ecoli_minimalAerobic()
 
