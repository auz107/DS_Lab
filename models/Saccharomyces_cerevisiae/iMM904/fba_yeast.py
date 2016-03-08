import sys
sys.path.append('/data/alizom')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
import cobra

# Needed to read the Pf model
sys.path.append('/data/alizom/models/Pfluorescens/iSB1139')

if __name__ == "__main__":

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Model path
    model_path = '/data/alizom/models/Saccharomyces_cerevisiae/iMM904/'

    # Define the organism
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '')

    # Load the original model
    iMM904 = read_sbml_model(file_name = model_path + 'iMM904.xml', model_name = 'iMM904',model_organism = model_organism, model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    iMM904.validate()

    # Set the growth medium
    set_specific_bounds(model = iMM904,file_name = model_path + 'iMM904_minimal_glucose.py',flux_bounds = {'EX_glc_e_':[-10,1000],'EX_o2_e_':[-2,1000]},simulation_condition = 'minimal_glucose_aerobic')

    for rxn in iMM904.reactions:
        if rxn.id == 'biomass_SC5_notrace':
            rxn.objective_coefficient = 1
        else:
            rxn.objective_coefficient = 0

    # FBA  
    iMM904.fba(optimization_solver = optimization_solver) 

    #--- Simulate histidone defect ---
    print '\n--- Simulation histidine synthesis defect ---'
    set_specific_bounds(model = iMM904,file_name = model_path + 'iMM904_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-2,1000],'IGPDH':[0,0]})
    iMM904.fba(optimization_solver = optimization_solver, create_model = True)

    set_specific_bounds(model = iMM904,file_name = model_path + 'iMM904_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_his_L_e_':[-0.5,1000],'EX_o2_e_':[-2,1000],'IGPDH':[0,0]})
    iMM904.fba(optimization_solver = optimization_solver, create_model = True)

