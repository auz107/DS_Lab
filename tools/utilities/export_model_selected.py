import sys
sys.path.append('../../')
from tools.userError import userError
from tools.core.organism import organism
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.gene import gene
from tools.core.model import model
from get_ModelSEED_ids import get_ModelSEED_ids
from copy import deepcopy

def export_model_selected(cpds_list = [], rxns_list = [], genes_list = [], output_format = 'sbml', output_filename = '', exported_model_id = '', organism_id = '', obtain_ModelSEED_ids = True, replace_with_ModelSEED_id = True):
    """
    This function exports a selected list of compounds and/or reactions to a format 
    acceptable by export() function (see the inputs of function model.export() for details)

    INPUTS:
    -------
    cpds_list:
    A list of selected compounds

    rxns_list: 
    A list of selected reactions to export

    exported_model_id:
    An id for the imported model

    obtain_ModelSEED_ids:
    If True, ModelSEED ids are obtained before export

    replace_with_ModelSEED_id:
    If True, the id of any compound or reaction with a unique ModelSEED id is replaced 
    with its corresponding ModelSEED id and a new attribute called original_model_id is added
    to the compound or reaction object.

    NOTE: If not input is provided for both cpds_list and rxns_list then the entire
          model will be exported

    Ali R. Zomorrodi - Segre lab @ BU
    Last updated: 06/03/2016
    """
    # Compartments to export
    cpts_list = list(set([cpd.compartment for cpd in cpds_list] + [cpt for rxn in rxns_list for cpt in rxn.compartments] + [cpt for gen in genes_list for cpt in gen.compartment]))
    cpts_to_export = [deepcopy(cpt) for cpt in cpts_list]
    cpts_by_id = dict([(cpt.id,cpt) for cpt in cpts_list])

    # Compounds
    cpds_to_export = []
    for cpd in list(set(cpds_list + [c for rxn in rxns_list for c in rxn.compounds])):
         cpd_copy = deepcopy(cpd)
         cpd_copy.compartment = cpts_by_id[cpd.compartment.id]
         cpd_copy.reset_props()
         cpds_to_export.append(cpd_copy)
    cpds_by_id = dict([(cpd.id,cpd) for cpd in cpds_to_export]) 

    # Genes
    genes_to_export = []
    for gen in list(set(genes_list + [g for rxn in rxns_list for g in rxn.genes])):
         gen_copy = deepcopy(gen)
         gen_copy.compartment = [cpts_by_id[cpt.id] for cpt in gen.compartment]
         genes_to_export.append(gen_copy)
    genes_by_id = dict([(gen.id,gen) for gen in genes_to_export]) 

    # reactions
    rxns_to_export = []
    for rxn in rxns_list:
        rxn_copy = deepcopy(rxn)
        rxn_copy.reset_props()
        rxn_stoic = {}
        for cpd in rxn.stoichiometry.keys():
            rxn_stoic[cpds_by_id[cpd.id]] = rxn.stoichiometry[cpd]
        rxn_copy.set_stoichiometry(stoichiometry = rxn_stoic)
        rxns_to_export.append(rxn_copy)

    # Exported model id 
    if exported_model_id == '':
        exported_model_id = 'Selected parts of a model'

    # Organism name
    if organism_id == '':
        org = None
    else:
        org = organism(id = organism_id) 
    
    model_to_export = model(id = exported_model_id, organism = org, compounds = tuple(cpds_to_export), reactions = tuple(rxns_to_export), genes = tuple(genes_to_export), compartments = cpts_to_export)  
    # Obtain ModelSEED ids
    if obtain_ModelSEED_ids:
        get_ModelSEED_ids(model = model_to_export, stdout_msgs = False)

    # Replace the model ids with the corresponding ModelSEED ids
    if replace_with_ModelSEED_id:
        for cpd in [c for c in model_to_export.compounds if len(c.ModelSEED_id) == 1]:
            cpd.original_model_id = cpd.id
            cpd.id = cpd.ModelSEED_id[0] + '_' + cpd.compartment.id

        for rxn in [r for r in model_to_export.reactions if len(r.ModelSEED_id) == 1]:
            rxn.original_model_id = rxn.id
            rxn.id = rxn.ModelSEED_id[0]

    # Validate the model
    model_to_export.validate()

    # Export the model
    model_to_export.export(output_format = output_format, output_filename = output_filename) 


