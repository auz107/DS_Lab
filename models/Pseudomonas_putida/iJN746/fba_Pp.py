import sys
sys.path.append('/data/alizom/')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
import cobra


if __name__ == "__main__":

    # Model path
    model_path = '/data/alizom/models/Pseudomonas_putida/iJN746/'

    # Define the organism
    model_organism = organism(id = 'Pputida', name = 'Pseudomonas putida',domain = 'Bacteria', genus = 'Pseudomonas', species = 'putida', strain = 'KT2440')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Read the sbml file
    iJN746 = read_sbml_model(file_name = model_path + 'iJN746.xml', model_name = 'iJN746',model_organism = model_organism,model_type = 'metabolic',import_bounds = False)

    # Check the model after the modifications 
    iJN746.validate()

    iJN746.print_reactions()

    """
    # Set the growth medium
    set_specific_bounds(model = iJN746,file_name = model_path + 'iJN746_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')

    # Assign and objective function coefficients
    for rxn in iJN746.reactions:
        if rxn.id == 'EX_EC9324':
            rxn.objective_coefficient = 1
        else:
            rxn.objective_coefficient = 0

    # FBA  
    iJN746.fba(optimization_solver = optimization_solver) 
    """

