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
    model_path = '/data/alizom/models/Enterobacter_aerogenes/kbase/'

    # Define the organism
    model_organism = organism(id = 'Enterobacter_aerogenes', name = 'Enterobacter aerogenes putida',domain = 'Bacteria', genus = 'Enterobacter', species = 'aerogenes', strain = 'KCTC 2190')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Read the sbml file
    Ea_kbase = read_sbml_model(file_name = model_path + 'Enterobacter_aerogenes.20.xml', model_id = 'Ea_kbase',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    Ea_kbase.validate()

    # Assign and objective function coefficients
    for rxn in Ea_kbase.reactions:
        rxn.objective_coefficient = 0
    Ea_kbase.get_reactions({'biomass0':'id'}).objective_coefficient = 1

    #--- List of carbon sources it should grow on according to experimental data ---
    carbon_sources = {
    '4_hydroxy_benzoic_acid':'EX_cpd00136_e0', 
    'alpha-D-lactose':'EX_cpd00208_e0', 
    #'B-methyl-D-glucoside':'EX_cpd11603_e0', 
    #'D-glucose-1-phosphate':'EX_cpd00089_e0', 
    'D-mannitol':'EX_cpd00314_e0',  
    'D,L-alpha-glycerol phosphate':'EX_cpd00080_e0', 
    'D-cellobiose':'EX_cpd00158_e0', 
    #'D-galactonic acid gamma lactone':'EX_cpd02143_e0',
    #'D-glucosaminic acid':'EX_cpd02351_e0', 
    #'D_malic_acid':'EX_cpd00386_e0',    
    'D-xylose':'EX_cpd00154_e0',  
    'galacturonic acid':'EX_cpd00280_e0',  
    'glycyl-L-glutamic acid':'EX_cpd11592_e0',
    'glucose':'EX_cpd00027_e0',   
    'L_arginine':'EX_cpd00051_e0',
    #'ctionL_asparagine':'EX_cpd00132_e0',  
    'L_serine':'EX_cpd00054_e0', 
    'N-acetyl-D-glucosamine':'EX_cpd00122_e0',
    'putrescine':'EX_cpd00118_e0',  
    #'phenylethylamine':'EX_cpd03161_e0', 
    } 

    print '\n--- Growth with no carbon source ---'
    # Set the growth medium
    set_specific_bounds(model = Ea_kbase,flux_bounds = {'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'kbase_minimal_medium.py')

    # FBA  
    Ea_kbase.fba(optimization_solver = optimization_solver) 

    print '\n--- Growth on glucose ---'
    # Set the growth medium
    set_specific_bounds(model = Ea_kbase,flux_bounds = {carbon_sources['glucose']:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'kbase_minimal_medium.py')

    # FBA  
    Ea_kbase.fba(optimization_solver = optimization_solver) 

    #--- Growth on other carbon sources ---
    for carbon_src in sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose']):
        print '\n--- Growth on {}  ---'.format(carbon_src)
        set_specific_bounds(model = Ea_kbase,flux_bounds = {carbon_sources[carbon_src]:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'kbase_minimal_medium.py')
        Ea_kbase.fba(optimization_solver = optimization_solver) 




