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
    model_path = '/data/alizom/models/Serratia_marcescens/kbase/'

    # Define the organism
    model_organism = organism(id = 'Enterobacter_aerogenes', name = 'Enterobacter aerogenes putida',domain = 'Bacteria', genus = 'Enterobacter', species = 'aerogenes', strain = 'KCTC 2190')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Read the sbml file
    Sm_kbase = read_sbml_model(file_name = model_path + 'Serratia_marcescens.12.xml', model_id = 'Sm_kbase',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    Sm_kbase.validate()

    # Assign and objective function coefficients
    for rxn in Sm_kbase.reactions:
        rxn.objective_coefficient = 0
    Sm_kbase.get_reactions({'biomass0':'id'}).objective_coefficient = 1

    #--- List of carbon sources it should grow on according to experimental data ---
    carbon_sources = {
    #'B-methyl-D-glucoside':'EX_cpd11603_e0', 
    #'D-glucose-1-phosphate':'EX_cpd00089_e0', 
    'D-mannitol':'EX_cpd00314_e0',  
    'glycyl-L-glutamic acid':'EX_cpd11592_e0',
    'galacturonic acid':'EX_cpd00280_e0',  
    'glucose':'EX_cpd00027_e0',
    #'i-erythritol':'EX_cpd00392_e0',
    #'L_asparagine':'EX_cpd00132_e0',  
    'L_serine':'EX_cpd00054_e0', 
    'N-acetyl-D-glucosamine':'EX_cpd00122_e0',
    'putrescine':'EX_cpd00118_e0',  
    #'tween 80':'EX_cpd13392_e0'
    } 

    print '\n--- Growth with no carbon source ---'
    # Set the growth medium
    set_specific_bounds(model = Sm_kbase,flux_bounds = {'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'kbase_minimal_medium.py')

    # FBA  
    Sm_kbase.fba(optimization_solver = optimization_solver) 

    print '\n--- Growth on glucose ---'
    # Set the growth medium
    set_specific_bounds(model = Sm_kbase,flux_bounds = {carbon_sources['glucose']:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'kbase_minimal_medium.py')

    # FBA  
    Sm_kbase.fba(optimization_solver = optimization_solver) 

    #--- Growth on other carbon sources ---
    for carbon_src in sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose']):
        print '\n--- Growth on {}  ---'.format(carbon_src)
        set_specific_bounds(model = Sm_kbase,flux_bounds = {carbon_sources[carbon_src]:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'kbase_minimal_medium.py')
        Sm_kbase.fba(optimization_solver = optimization_solver) 




