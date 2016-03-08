import sys
sys.path.append('/data/alizom')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
from copy import deepcopy
import cobra
# Needed to read the Pf model
sys.path.append('/data/alizom/models/Pfluorescens/iSB1139')
from read_excel_model import create_model

if __name__ == "__main__":
 
    # Path to the model
    model_path = '/data/alizom/models/Pseudomonas_fluorescens/iSB1139/'

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    model_organism = organism(id = 'Pfluorescens', name = 'Pseudomonas fluorescens',domain = 'Bacteria', genus = 'Pseudomonas', species = 'fluorescens', strain = 'SBW25')

    print '\n---- Model from the original SBML file -------\n'
    #iSB1139_origSBML = read_sbml_model(file_name = model_path + 'iSB1139_orig.xml', model_name = 'iSB1139',model_organism = model_organism,model_type = 'metabolic',import_bounds = False)
    iSB1139_origSBML = read_sbml_model(file_name = model_path + 'pFluorescens.xml', model_name = 'iSB1139',model_organism = model_organism,model_type = 'metabolic',import_bounds = False)

    # Remove all compounds not invovled in a reaction
    iSB1139_origSBML.remove_compounds([m for m in iSB1139_origSBML.compounds if len(m.reactions) == 0])

    # Remove all reactions with no reactants and products 
    iSB1139_origSBML.remove_reactions([r for r in iSB1139_origSBML.reactions if len(r.compounds) == 0])
    iSB1139_origSBML.validate()

    # Export the updated model to a new SBML file using COBRA toolbox
    cobra.io.write_sbml_model(iSB1139_origSBML.export_to_cobra(),model_path + 'iSB1139_origSBML_cobra.xml',use_fbc_package=False)

    # Set the objective function coefficients
    iSB1139_origSBML.biomass_reaction = iSB1139_origSBML.get_reactions({'BIOMASSRXN':'id'})
    for rxn in iSB1139_origSBML.reactions:
        rxn.objective_coefficient = 0
    iSB1139_origSBML.biomass_reaction.objective_coefficient = 1

    for rxn in [r for r in iSB1139_origSBML.reactions if r.type.lower() == 'exchange']:
        rxn.flux_bounds[0] = 0

    iSB1139_origSBML.fba(optimization_solver = optimization_solver) 

    print '\n---- Model from the Excel file -------\n'
    # Read the model using the customized model reader 
    iSB1139_excel = create_model(python_model_file_name = model_path + 'iSB1139.py',organism_name = 'Pfluorescens', model_type = 'metabolic', model_name = 'iSB1139')  

    iSB1139_excel.validate()

    iSB1139_excel.biomass_reaction = iSB1139_excel.get_reactions({'BIOMASSRXN':'id'})

    # Set the objective function coefficients
    for rxn in iSB1139_excel.reactions:
        rxn.objective_coefficient = 0
    iSB1139_excel.biomass_reaction.objective_coefficient = 1

    # Export the updated model to a new SBML file
    cobra.io.write_sbml_model(iSB1139_excel.export_to_cobra(),model_path + 'iSB1139_excel_cobra.xml',use_fbc_package=False)

    # Set the growth medium
    set_specific_bounds(model = iSB1139_excel,file_name = model_path + 'iSB1139_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')

    iSB1139_excel.fba(optimization_solver = optimization_solver) 
         
