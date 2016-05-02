from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from tools.core.organism import organism
from tools.core.compartment import compartment
from tools.core.gene import gene
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
from imp import load_source
import re

def read_pydict_model(file_name, import_params = False, validate = True, warnings = True, stdout_msgs = True): 
    """
    A class to read a sbml model into a model object. The sbml model is 
    first converted into a COBRApy model, which is used as an input to create
    the corresponding model object. 

    INPUTS: 
    ------
           file_name: The name of the python file containing the model (String)
                      This string can contain the full path to the file  
            model_id: Name of the model (string). If empty, model_id is read from the input file 
          model_type: Type of the model (string, e.g., 'metabolic'). If empty, it is read from the input file
       import_params: Indicates whether the model parameters (flux bounds, objective coefficients)  
                      should be imported from the input file (True) or not (False). The default is False
            warnings: Can be 'on' or 'off' shwoing whether the warnings should be writtten to the 
                      screen or not
         stdout_msgs: Can be 'on' or 'off' shwoing whether any messages should be written to the
                      screen

    OUTPUT:
    -------
               model: A model object

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 04-18-2016
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError("stdout_msgs must be True or False")
    if not isinstance(warnings,bool):
        raise TypeError("warnings must be True or False")


    # Load the python dictionary containing the model information
    load_source('dataFile',file_name)
    import dataFile
    model_pydict = dataFile.model
 

    #--- organism ---
    model_organism = organism(id = model_pydict['organism']['id'], name = model_pydict['organism']['name'], domain = model_pydict['organism']['domain'], genus = model_pydict['organism']['genus'], species = model_pydict['organism']['species'], strain = model_pydict['organism']['strain'], ModelSEED_type = model_pydict['organism']['ModelSEED_type'], gDW_per_cell = model_pydict['organism']['gDW_per_cell'], gWW_per_cell = model_pydict['organism']['gWW_per_cell'], cells_per_ml = model_pydict['organism']['cells_per_ml'], gDW_per_ml = model_pydict['organism']['gDW_per_ml'], mu = model_pydict['organism']['mu'], random_mortality_rate = model_pydict['organism']['random_mortality_rate'], notes = model_pydict['organism']['notes'])
    for attr in [k for k in model_pydict['organism'].keys() if k not in ['id','name', 'domain', 'genus', 'species', 'strain', 'ModelSEED_type', 'gDW_per_cell', 'gDW_per_cell', 'cells_per_ml', 'gDW_per_ml', 'mu', 'random_mortality_rate', 'notes']]:
        model_organism.__setattr__(attr, model_pydict['organism'][attr])
        

    #--- compartments ---
    compartments = []
    compartments_by_id = {}
    for cpt_id in model_pydict['compartments'].keys():
        cpt = compartment(id = model_pydict['compartments'][cpt_id]['id'], name = model_pydict['compartments'][cpt_id]['name'], name_aliases = model_pydict['compartments'][cpt_id]['name_aliases'], notes = model_pydict['compartments'][cpt_id]['notes'])
        for attrs in [k for k in model_pydict['compartments'][cpt_id].keys() if k not in ['id', 'name', 'synonyms', 'model', 'notes']]:
            cpt.__setattr__(attr, model_pydict['compartments'][cpt_id][attr])
        compartments_by_id[model_pydict['compartments'][cpt_id]['id']] = cpt
    compartments.append(cpt)

    #--- compounds ---
    compounds = []
    compounds_by_id = {}
    for cpd_id in model_pydict['compounds'].keys():
        cpd = compound(id = model_pydict['compounds'][cpd_id]['id'], compartment = compartments_by_id[model_pydict['compounds'][cpd_id]['compartment']], name = model_pydict['compounds'][cpd_id]['name'], name_aliases = model_pydict['compounds'][cpd_id]['name_aliases'], KEGG_id = model_pydict['compounds'][cpd_id]['KEGG_id'], ModelSEED_id = model_pydict['compounds'][cpd_id]['ModelSEED_id'], BiGG_id = model_pydict['compounds'][cpd_id]['BiGG_id'], formula = model_pydict['compounds'][cpd_id]['formula'], molecular_weight = model_pydict['compounds'][cpd_id]['molecular_weight'], charge = model_pydict['compounds'][cpd_id]['charge'], reactions = [], reactant_reactions = [], product_reactions = [], concentration = model_pydict['compounds'][cpd_id]['concentration'], deltaG = model_pydict['compounds'][cpd_id]['deltaG'], deltaG_uncertainty = model_pydict['compounds'][cpd_id]['deltaG_uncertainty'], deltaG_range = model_pydict['compounds'][cpd_id]['deltaG_range'], notes = model_pydict['compounds'][cpd_id]['notes'], warnings = model_pydict['compounds'][cpd_id]['warnings'])

        for attr in [k for k in model_pydict['compounds'][cpd_id].keys() if k not in ['id', 'compartment', 'name', 'name_aliases', 'KEGG_id', 'ModelSEED_id', 'BiGG_id', 'formula', 'molecular_weight', 'charge', 'reactions', 'reactant_reactions', 'product_reactions', 'concentration', 'deltaG', 'deltaG_uncertainty', 'deltaG_range', 'notes', 'warnings']]:
            cpd.__setattr__(attr, model_pydict['compounds'][cpd_id][attr])
        compounds.append(cpd)
        compounds_by_id[model_pydict['compounds'][cpd_id]['id']] = cpd

    #--- genes ---
    genes = []
    genes_by_id = {}
    for gen_id in model_pydict['genes'].keys():
        gen = gene(id = model_pydict['genes'][gen_id]['id'], name = model_pydict['genes'][gen_id]['name'], name_aliases = model_pydict['genes'][gen_id]['name_aliases'], compartment = [compartments_by_id[c] for c in model_pydict['genes'][gen_id]['compartment']], reactions = [], locus_pos = model_pydict['genes'][gen_id]['locus_pos'], expression_level = model_pydict['genes'][gen_id]['expression_level'], notes = model_pydict['genes'][gen_id]['notes'])

        for attr in [k for k in model_pydict['genes'][gen_id].keys() if k not in ['id','name','name_aliases','compartment','reactions','locus_pos', 'expression_level', 'notes']]:
            gen.__setattr__(attr, model_pydict['genes'][gen_id][attr])
        genes.append(gen)
        genes_by_id[model_pydict['genes'][gen_id]['id']] = gen

    #--- reactions ---
    reactions = []
    reactions_by_id = {}
    for rxn_id in model_pydict['reactions'].keys():
        r_stoic = dict([(compounds_by_id[c],model_pydict['reactions'][rxn_id]['stoichiometry'][c]) for c in model_pydict['reactions'][rxn_id]['stoichiometry'].keys()])
        
        if import_params:
            r_flux_bounds = model_pydict['reactions'][rxn_id]['flux_bounds']
        else:
            r_flux_bounds = []

        rxn = reaction(id = model_pydict['reactions'][rxn_id]['id'], stoichiometry = r_stoic, reversibility = model_pydict['reactions'][rxn_id]['reversibility'], name = model_pydict['reactions'][rxn_id]['name'], name_aliases = model_pydict['reactions'][rxn_id]['name_aliases'], KEGG_id = model_pydict['reactions'][rxn_id]['KEGG_id'], ModelSEED_id = model_pydict['reactions'][rxn_id]['ModelSEED_id'], BiGG_id = model_pydict['reactions'][rxn_id]['BiGG_id'], EC_numbers = model_pydict['reactions'][rxn_id]['EC_numbers'], subsystem = model_pydict['reactions'][rxn_id]['subsystem'], pathways = model_pydict['reactions'][rxn_id]['pathways'], compartment = [compartments_by_id[c] for c in model_pydict['reactions'][rxn_id]['compartment']], genes = [genes_by_id[g] for g in model_pydict['reactions'][rxn_id]['genes']], gene_reaction_rule = model_pydict['reactions'][rxn_id]['gene_reaction_rule'], objective_coefficient = model_pydict['reactions'][rxn_id]['objective_coefficient'], flux = model_pydict['reactions'][rxn_id]['flux'], store_flux = model_pydict['reactions'][rxn_id]['store_flux'], flux_bounds = r_flux_bounds, deltaG = model_pydict['reactions'][rxn_id]['deltaG'], deltaG_uncertainty = model_pydict['reactions'][rxn_id]['deltaG_uncertainty'], deltaG_range = model_pydict['reactions'][rxn_id]['deltaG_range'], kinetics = model_pydict['reactions'][rxn_id]['kinetics'], kinetic_compounds = model_pydict['reactions'][rxn_id]['kinetic_compounds'], confidence_level = model_pydict['reactions'][rxn_id]['confidence_level'], notes = model_pydict['reactions'][rxn_id]['notes'], warnings = model_pydict['reactions'][rxn_id]['warnings'])

        for attr in [k for k in model_pydict['reactions'][rxn_id].keys() if k not in ['id', 'stoichiometry', 'reversibility', 'name', 'name_aliases', 'KEGG_id', 'ModelSEED_id', 'BiGG_id', 'EC_numbers', 'subsystem', 'pathways', 'compartment', 'genes', 'gene_reaction_rule', 'objective_coefficient', 'flux', 'store_flux', 'flux_bounds', 'deltaG', 'deltaG_uncertainty', 'deltaG_range', 'kinetics', 'kinetic_compounds', 'confidence_level', 'notes', 'warnings', 'compounds', 'reactants', 'products']]:
            rxn.__setattr__(attr, model_pydict['reactions'][rxn_id][attr])
        reactions.append(rxn)
        reactions_by_id[model_pydict['reactions'][rxn_id]['id']] = rxn

    # biomass rxn
    if model_pydict['model_instance_attrs']['biomass_reaction'] != None:
        biomass_rxn = reactions_by_id[model_pydict['model_instance_attrs']['biomass_reaction']], 
    else:
        biomass_rxn = None

    # ATPM rxn
    if model_pydict['model_instance_attrs']['atpm_reaction'] != None:
        atpm_rxn = reactions_by_id[model_pydict['model_instance_attrs']['atpm_reaction']], 
    else:
        atpm_rxn = None

    #--- model ---
    model_class = model(id = model_pydict['model_instance_attrs']['id'], 
                        type = model_pydict['model_instance_attrs']['type'], 
                         organism = model_organism, 
                         reactions = reactions, 
                         compounds = compounds, 
                         genes = genes, 
                         compartments = compartments, 
                         name = model_pydict['model_instance_attrs']['name'], 
                         biomass_reaction = biomass_rxn,  
                         atpm_reaction = atpm_rxn,  
                         notes = model_pydict['model_instance_attrs']['notes'], 
                         validate = validate, 
                         stdout_msgs = model_pydict['model_instance_attrs']['stdout_msgs'], 
                         warnings = model_pydict['model_instance_attrs']['warnings'])

    for attr in [k for k in model_pydict['model_instance_attrs'].keys() if k not in ['id',' organis', 'compounds', 'reactions', 'genes', 'compartments', 'name', 'biomass_reaction', 'atpm_reaction', 'notes', 'stdout_msgs', 'warnings']]:
       model_class.__setattr__(attr, model_pydict['model_instance_attrs'][attr]) 
    
    return model_class


