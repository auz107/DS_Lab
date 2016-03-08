import sys
sys.path.append('../../')
from tools.userError import *
from tools.core.organism import organism
from tools.core.compartment import compartment
from tools.core.gene import gene
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
import cobra

def write_gams_model(model,file_name = None)
    """
    Writes a model into text to be used as inputs for constructing a GAMS optimization model

    INPUTS:
    -------
        model: Input model (an instance of object model)
    file_name: If a directory and/or file name is given it is added to the begining of all output 
               files, otherwise general file names (reactions.txt, metabolites.txt, reaction_types.txt,
               stoichiometric_matrix.txt and genes.txt) are written in the current directory
               Example: If file_name = 'results/Ecoli' then metabolite names are written into
               results/Ecoli_metabolites.txt
 
    Ali R. Zomorrodi, Daneil Segre Lab @ BU 
    Last updated: April 06, 2015
    """
    if file_name != None:
        # Name of the file for metabolite names
        metab_file_name = file_name + '_metabolites.txt' 

        # Name of the file for reaction names
        rxn_file_name = file_name + '_reactions.txt' 

        # Name of the file for reaction types
        rxntype_file_name = file_name + '_reaction_types.txt' 

        # Name of the file for stoichiometric matrix 
        stoic_file_name = file_name + '_stoichiometric_matrix.txt' 

        # Name of the file for genes 
        gene_file_name = file_name + '_genes.txt' 

    else:
        # Name of the file for metabolite names
        metab_file_name = 'metabolites.txt' 

        # Name of the file for reaction names
        rxn_file_name = 'reactions.txt' 

        # Name of the file for reaction types
        rxntype_file_name = 'reaction_types.txt' 

        # Name of the file for stoichiometric matrix 
        stoic_file_name = 'stoichiometric_matrix.txt' 

        # Name of the file for genes 
        gene_file_name = 'genes.txt' 


    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Model path
    model_path = '/data/alizom/models/Saccharomyces_cerevisiae/iAZ900/'

    # Define the organism
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '')

    #------- Original iAZ900 model -----
    print '\n----- Original iAZ900 model -----\n'
    # Load the original model
    iAZ900 = read_sbml_model(file_name = model_path + 'iAZ900.xml', model_name = 'iAZ900',model_organism = model_organism, model_type = 'metabolic',import_bounds = False)

    # Check the model after the modifications 
    iAZ900.validate()

    # Set the growth medium
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')

    for rxn in iAZ900.reactions:
        if rxn.id == 'biomass_wildType':
        #if rxn.id == 'biomass_core':
            rxn.objective_coefficient = 1
        else:
            rxn.objective_coefficient = 0

    # FBA  
    iAZ900.fba(optimization_solver = optimization_solver) 

    #------- iAZ900 model modified with cycles removed (by Anupam) -----
    print '\n----- iAZ900_noCyles model modified with cycles removed (by Anupam) -----\n'
    # Load the original model
    iAZ900_noCyles = read_sbml_model(file_name = model_path + 'iAZ900_Anupam_March_25_2015.xml', model_name = 'iAZ900_noCyles',model_organism = model_organism, model_type = 'metabolic',import_bounds = False)

    # Check the model after the modifications 
    iAZ900_noCyles.validate()

    with open('metabolites.txt','w') as f:
        f.write('\/\n')
        for cmp in iAZ900_noCyles.compounds:
            f.write("'" + cmp.id + "'\n")
        f.write('\/\n')

    with open('reactions.txt','w') as f:
        f.write('\/\n')
        for rxn in iAZ900_noCyles.reactions:
            f.write("'" + rxn.id + "'\n")
        f.write('\/\n')

    with open('reaction_types.txt','w') as f:
        f.write('\/\n')
        for rxn in iAZ900_noCyles.reactions:
            if rxn.type.lower() == 'irreversible':
                f.write("'" + rxn.id + "'\t0\n")
            elif rxn.type.lower() == 'reversible':
                f.write("'" + rxn.id + "'\t1\n")
            elif rxn.type.lower() == 'exchange':
                f.write("'" + rxn.id + "'\t3\n")
            else:
                raise userError('unknonw reaction name')
        f.write('\/\n')

    with open('stoichiometric_matrix.txt','w') as f:
        for cmp in iAZ900_noCyles.compounds:
            for rxn in cmp.reactant_reactions:
                f.write("S('" + cmp.id + "','" + rxn.id + "') = " + rxn.stoichiometry[cmp] +';\n')
            for rxn in cmp.product_reactions:
                f.write("S('" + cmp.id + "','" + rxn.id + "') = " + rxn.stoichiometry[cmp] +';\n')
            f.write('\n')

    # Set the growth medium
    set_specific_bounds(model = iAZ900_noCyles,file_name = model_path + 'iAZ900_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')

    for rxn in iAZ900_noCyles.reactions:
        if rxn.id == 'biomass_wildType':
        #if rxn.id == 'biomass_core':
            rxn.objective_coefficient = 1
        else:
            rxn.objective_coefficient = 0

    # FBA  
    iAZ900_noCyles.fba(optimization_solver = optimization_solver) 

