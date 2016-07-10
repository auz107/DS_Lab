import sys
sys.path.append('../../../')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
from tools.utilities.export_model_selected import export_model_selected
import cobra, re, copy

if __name__ == "__main__":

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Pseudomonas_putida/iJN746/'

    # Define the organism
    model_organism = organism(id = 'Pputida', name = 'Pseudomonas putida',domain = 'Bacteria', genus = 'Pseudomonas', species = 'putida', strain = 'KT2440')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Read the sbml file
    iJN746 = read_sbml_model(file_name = model_path + 'iJN746.xml', model_name = 'iJN746',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    iJN746.validate()

    #--- Save biomass reaction as a pydict model to a file ---
    iJN746.biomass_reaction = iJN746.reactions_by_id['BiomassKT_TEMP2']
    iJN746.compounds_by_id['h2o_c'].ModelSEED_id = ['cpd00001'] # To avoid assigning cpd15275 

    biomass_rxn = copy.deepcopy(iJN746.biomass_reaction)
    biomass_rxn.id = 'Pseudomonas_putida_iJN746_biomass_rxn'
    
    # Change compartment of compounds that is c to co and also p to c0
    cpts_by_id = dict([(cpt.id,cpt) for cpt in biomass_rxn.compartments])
    cpts_by_id['c'].id = 'c0' # Change c to c0
    for cpd in [c for c in biomass_rxn.compounds if c.compartment.id == 'c0']:
        # Replace _p in their id with _c0
        cpd.id = re.sub('_c$','_c0',cpd.id)
    for cpd in [c for c in biomass_rxn.compounds if c.compartment.id == 'p']:
        cpd.compartment = cpts_by_id['c']
        # Replace _p in their id with _c0
        cpd.id = re.sub('_p$','_c0',cpd.id)
    biomass_rxn.assign_props()   
    
    export_model_selected(rxns_list = [biomass_rxn], output_format = 'pydict', output_filename = 'iJN746_biomass_rxn_model.py', exported_model_id = 'iJN746_biomass_rxn', organism_id = 'Pseudomonas_putida', obtain_ModelSEED_ids = True, replace_with_ModelSEED_id = True)

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

