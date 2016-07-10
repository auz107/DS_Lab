from __future__ import division
import sys, time
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
sys.path.append('../../')
from tools.userError import *
from tools.core.organism import organism
from tools.core.compartment import compartment
from tools.core.gene import gene
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
from tools.ancillary.get_ModelSEED_ids import get_ModelSEED_ids, remove_compartment
from models.ModelSEED.ModelSEED_cpds_master import cpds_master as ModelSEED_cpds
from models.ModelSEED.ModelSEED_rxns_master import rxns_master as ModelSEED_rxns
import itertools

def integrate_models(models_list, add_biomass_rxns = False, replace_cpt_ids = True, obtain_ModelSEED_ids = True, warnings = True, stdout_msgs = True): 
    """
    INPUTS: 
    ------
    models_list: 
    A list of dictionaries containing the information about the models to be integrates. 
    The keys of this dictionary are as follows:
        'model': The model object

        'standard_to_model_compartID_map': 
         A dictionary where keys are one of the letters below (standard compartment ids and 
         values are the corresponding compartment ids in the model id in the model. 
         c: Cytosol (cytoplasm),   e: Extracellular,  g: Golgi,     m: Mitochondria   
         n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,   x: Peroxisome
         For example, if a model has two compartments c and e, one can provide  {'c':'c id in 
         the model', 'e':'e id in the model'}. One can also provide {'c':'', 'e': ''} in which 
         the code searches for these two compartments in the model

         'biomass_rxn_id': 
         Id of the biomass reaction in the model (optional). Enter an empty string('') if you 
         don't wish to provide

    add_biomass_rxns: 
    If True, biomass reactions are all added to the integrated model

    replace_cpt_ids:
    If True, compartment ids in the compound ids are replaced with the corresponding standard
    compartment ids. Set the value to False if all models have the same compartment ids (for
    example if you are integrating multiples models created on KBase or ModelSEED)
 
    OUTPUT:
    -------
    integrated_model: The integrated model 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06-01-2016
    """
    if not isinstance(models_list,dict):
        raise TypeError('models_list must be a list of dictionaries. An object of type {} was provided instead'.format(type(models_list)))
    elif len([m.id for m in models if len([k for k in m.keys() if k.lower() not in ['model','standard_to_model_compartID_map', 'biomass_rxn_id']]) > 0]) > 0:
        raise userError('Non-authorized keys were observed for the following models: {}'.format([m.id for m in models if len([k for k in m.keys() if k.lower() not in ['model','standard_to_model_compartID_map', 'biomass_rxn_id']]) > 0]))

    if not isinstance(add_biomass_rxns,bool):
        raise TypeError('add_biomass_rxns must be either True or False')
    if not isinstance(obtain_ModelSEED_ids,bool):
        raise TypeError('obtain_ModelSEED_ids must be either True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True or False')
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True or False')

    # First obtain ModelSEED ids for all compounds and reactions in all models
    if obtain_ModelSEED_ids:
        if stdout_msgs:
            print 'Getting ModelSEED ids for all models in the list ...'
        for m in [mo['model'] for mo in models_list]:
            get_ModelSEED_ids(model = m, warnings = warnings, stdout_msgs = False)

    # For each model create a dictionary mapping the compartment ids in the model to standard
    # compartment ids
    for mdl in models_list:
        mdl['model_to_standard_compartID_map'] = dict([(v,k) for (k,v) in mdl['standard_to_model_compartID_map'].items()])

    # ----- organism ------
    models_with_organism = [mdl for mdl in models if mdl['model'].organism 1= None]
    if len(models_with_organism) > 0:
        integMdl_org = models_with_organism[0]['model'].organism
    
        if len(models_with_organism) > 1:
            for mdl in models_with_organism[1:]:
                integMdl_org.name_aliases = list(set(integMdl_org.name_aliases + [mdl['model'].organism.name] + mdl['model'].organism.name-aliases))

                if integMdl_org == '' and mdl['model'].organism.domain != '':
                    integMdl_org.domain = mdl['model'].organism.domain
                if integMdl_org == '' and mdl['model'].organism.genus != '':
                    integMdl_org.genus = mdl['model'].organism.genus
                if integMdl_org == '' and mdl['model'].organism.species != '':
                    integMdl_org.species = mdl['model'].organism.species
                if integMdl_org == '' and mdl['model'].organism.strain != '':
                    integMdl_org.strain = mdl['model'].organism.strain
                if integMdl_org == '' and mdl['model'].organism.ModelSEED_type != '':
                    integMdl_org.ModelSEED_type = mdl['model'].organism.ModelSEED_type

    # ----- compartments ------
    # HEre, we first find out all compartments present in all models and then create new
    # compartment objects used for the integrated model
    # A dicitonary with names being standard comprtment ids and values being their names
    standard_cpt_id_name_map = {'c':'Cytosol', 'e': 'Extracellular', 'g':'Golgi', 'm':'Mitochondria', 'n':'Nucleus', 'p':'Periplasm', 'r':'Endoplasmic reticulum', 'x':'Peroxisome'} 

    comparts = []
    comparts_by_id = {}

    # If replace cpt ids with standard ones
    if replace_cpt_ids:
        cpts_list = list(set([c for mdl in models_list for c in mdl['standard_to_model_compartID_map'].keys()]))
        standard_to_integModel_cpt_map = dict([(c,c) for c in cpts_list])

    # Typically, if we don't want to replace compartment ids with standard compartment ids, all models
    # have the same comaprtment ids (e.g., integrating multiple ModelSEED models)
    else:
        cpts_list = list(set([c for mdl in models_list for c in mdl['standard_to_model_compartID_map'].values()]))
        standard_to_integModel_cpt_map =  dict(set([c for mdl in models_list for c in mdl['standard_to_model_compartID_map'].items()]))


    for cpt in cpts_list: 
        compart = compartment(id = cpt, name = standard_cpt_id_name_map[cpt]) 
        comparts.append(compart)
        comparts_by_id[cpt] = compart
    comparts = tuple(comparts)

    # ----- genes ------
    genes = []
    genes_by_id = {}

    # A dictionary showing the map between genes in different models that are not added to 
    # genes (for integrated model) because there is already a match for that gene included in
    # genes. The keys are a tuple with the first element being the model id and the second element
    # being a genes object in that model and values are a gene object in genes (for 
    # integrated model). Example: {('model1',gen1):gen0, ...}
    not_included_genes_map = {}

    # First add all compounds in the first model
    # First add all compounds in the first model
    for gen in models_list[0]['model'].genes:
        gen.reset_props()
        genes.append(gen)
        genes_by_id[gen.id] = gen

    # Next, compare the rest of the models one-by-one with integrated_model and add any
    # gene for which a unique map was not found
    for mdl in models_list[1:]:
        if stdout_msgs:
            print '\nAdding genes from model {} to integrated model ...\n'.format(mdl['model'].id)

        genInegMdl_to_genMdl_map, genMdl_to_genIntegMdl_map = compare_genes(genes_list1 = genes, genes_list2 = mdl['model'].genes, genes_list1_id = 'integrate_model', genes_list2_id = mdl['model'].id, warnings = warnings, stdout_msgs = False)[0:2]

        # Update genes attributes for those genes in mdl with a unique match in 
        # integrated_model
        for gen in [g for g in genMdl_to_genIntegMdl_map if len(genMdl_to_genIntegMdl_map[gen]) == 1]:

            not_included_gens_map[(mdl['model'].id,gen)] = genMdl_to_genIntegMdl_map[gen][0]            

            genMdl_to_genIntegMdl_map[gen][0].name_aliases = list(set(genMdl_to_genIntegMdl_map[gen][0].name_aliases + [gen.name] + gen.name_aliases))

        # Consinder only genes in mdl for which there is no unique map in integrated_model
        for gen in [g for g in genMdl_to_genIntegMdl_map if len(genMdl_to_genIntegMdl_map[gen]) != 1]:
            gen.reset_props()
            gen.compartment = [comparts_by_id[mdl['model_to_standard_compartID_map'][cpt.id]] for cpt in gen.compartment]
            genes.append(gen)
            genes_by_id[gen.id] = gen

    # ----- Compounds ------
    compounds = []
    compounds_by_id = {}

    # A dictionary showing the map between compounds in different models that are not added to 
    # compounds (for integrated model) because there is already a match for that compound included in
    # compounds. The keys are a tuple with the first element being the model id and the second element
    # being a compounds object in that model and values are a compound object in compounds (for 
    # integrated model). Example: {('model1',cpd1):cpd0, ...}
    not_included_cpds_map = {}

    # First add all compounds in the first model
    for cpd in models_list[0]['model'].compounds:
        cpd.reset_props()
        cpd.compartment = comparts_by_id[models_list[0]['model_to_standard_compartID_map'][cpd.compartment.id]]
        if replace_cpt_ids:
            # Replace compartment id in the compound id with _standard_comaprtment_id
            cpd.id = re.sub('\[' + cpd.compartment.id + '\]|_' + cpd.compartment.id, '_' + models_list[0]['model_to_standard_compartID_map'][cpd.compartment.id], cpd.id)
        compounds.append(cpd)
        compounds_by_id[cpd.id] = cpd

    # Next, compare the rest of the models one-by-one with integrated_model and add any
    # compound for which a unique map was not found
    for mdl in models_list[1:]:
        if stdout_msgs:
            print '\nAdding compounds from model {} to integrated model ...\n'.format(mdl['model'].id)

        cpdInegMdl_to_cpdMdl_map, cpdMdl_to_cpdIntegMdl_map = compare_compounds(cpds_list1 = compounds, cpds_list2 = mdl['model'].compounds, standard_to_cpds_list1_cpt_ids_map = , standard_to_cpds_list2_cpt_ids_map = mdl['standard_to_model_compartID_map'], cpds_list1_id = 'integrate_model', cpds_list2_id = mdl['model'].id, stdout_msgs = False)[0:2]


        # Update compounds attributes for those compounds in mdl with a unique match in 
        # integrated_model
        for cpd in [c for c in cpdMdl_to_cpdIntegMdl_map if len(cpdMdl_to_cpdIntegMdl_map[cpd]) == 1]:

            not_included_cpds_map[(mdl['model'].id,cpd)] = cpdMdl_to_cpdIntegMdl_map[cpd][0]            

            cpdMdl_to_cpdIntegMdl_map[cpd][0].name_aliases = list(set(cpdMdl_to_cpdIntegMdl_map[cpd][0].name_aliases + [cpd.name] + cpd.name_aliases))

            cpdMdl_to_cpdIntegMdl_map[cpd][0].KEGG_id = list(set(cpdMdl_to_cpdIntegMdl_map[cpd][0].KEGG_id + cpd.KEGG_id)) 
            cpdMdl_to_cpdIntegMdl_map[cpd][0].ModelSEED_id = list(set(cpdMdl_to_cpdIntegMdl_map[cpd][0].ModelSEED_id + cpd.ModelSEED_id)) 
            cpdMdl_to_cpdIntegMdl_map[cpd][0].BiGG_id = list(set(cpdMdl_to_cpdIntegMdl_map[cpd][0].BiGG_id + cpd.BiGG_id)) 

            if cpd.formula != '' and cpdMdl_to_cpdIntegMdl_map[cpd][0].formulat == '':
                cpdMdl_to_cpdIntegMdl_map[cpd][0].formula = cpd.formula 

            if cpd.molecular_weight != None and cpdMdl_to_cpdIntegMdl_map[cpd][0].molecular_weightt == None:
                cpdMdl_to_cpdIntegMdl_map[cpd][0].molecular_weight = cpd.molecular_weight 

            if cpd.charge != None and cpdMdl_to_cpdIntegMdl_map[cpd][0].charget == None:
                cpdMdl_to_cpdIntegMdl_map[cpd][0].charge = cpd.charge 

            if cpd.deltaG != None and cpdMdl_to_cpdIntegMdl_map[cpd][0].deltaGt == None:
                cpdMdl_to_cpdIntegMdl_map[cpd][0].deltaG = cpd.deltaG 

            if cpd.deltaG_uncertainty != None and cpdMdl_to_cpdIntegMdl_map[cpd][0].deltaG_uncertaintyt == None:
                cpdMdl_to_cpdIntegMdl_map[cpd][0].deltaG_uncertainty = cpd.deltaG_uncertainty 

            if cpd.deltaG_range != [] and cpdMdl_to_cpdIntegMdl_map[cpd][0].deltaG_ranget == []:
                cpdMdl_to_cpdIntegMdl_map[cpd][0].deltaG_range = cpd.deltaG_range 

        # Consinder only compounds in mdl for which there is no unique map in integrated_model
        for cpd in [c for c in cpdMdl_to_cpdIntegMdl_map if len(cpdMdl_to_cpdIntegMdl_map[cpd]) != 1]:
            cpd.reset_props()
            cpd.compartment = comparts_by_id[mdl['model_to_standard_compartID_map'][cpd.compartment.id]]
            if replace_cpt_ids:
                # Replace compartment id in the compound id with _ssntadard_comaprtment_id
                cpd.id = re.sub('\[' + cpd.compartment.id + '\]|_' + cpd.compartment.id, '_' + mdl['model_to_standard_compartID_map'][cpd.compartment.id], cpt.id)
            compounds.append(cpd)
            compounds_by_id[cpd.id] = cpd

    # ----- Reactions ------
    reactions = []
    reactions_by_id = {}

    # First add all reactions in the first model
    for rxn in models_list[0]['model'].reactions:
        rxn.reset_props()
        rxn.compartment = [comparts_by_id[models_list[0]['model_to_standard_compartID_map'][rxn.cpt.id]] for cpt in rxn.compartment]

        reactions.append(rxn)
        reactions_by_id[rxn.id] = rxn

    # Next, compare the rest of the models one-by-one with integrated_model and add any
    # reaction for which a unique map was not found
    for mdl in models_list[1:]:
        if stdout_msgs:
            print '\nAdding reactions from model {} to integrated model ...\n'.format(mdl['model'].id)

        rxnInegMdl_to_rxnMdl_map, rxnMdl_to_rxnIntegMdl_map = compare_reactions(rxns_list1 = reactions, rxns_list2 = mdl['model'].reactions, standard_to_rxns_list1_cpt_ids_map = standard_to_integModel_cpt_map, standard_to_rxns_list2_cpt_ids_map = mdl['standard_to_model_compartID_map'], rxns_list1_id = 'integrate_model', rxns_list2_id = mdl['model'].id, cpd1_to_cpd2_map = cpdInegMdl_to_cpdMdl_map[mdl['model'].id], stdout_msgs = False)[0:2]

        # Update reactions attributes for those reactions in mdl with a unique match in 
        # integrated_model
        for rxn in [r for r in rxnMdl_to_rxnIntegMdl_map if len(rxnMdl_to_rxnIntegMdl_map[rxn]) == 1]:

            rxnMdl_to_rxnIntegMdl_map[rxn][0].name_aliases = list(set(rxnMdl_to_rxnIntegMdl_map[rxn][0].name_aliases + [rxn.name] + rxn.name_aliases)) 

            rxnMdl_to_rxnIntegMdl_map[rxn][0].KEGG_id = list(set(rxnMdl_to_rxnIntegMdl_map[rxn][0].KEGG_id + rxn.KEGG_id)) 
            rxnMdl_to_rxnIntegMdl_map[rxn][0].ModelSEED_id = list(set(rxnMdl_to_rxnIntegMdl_map[rxn][0].ModelSEED_id + rxn.ModelSEED_id)) 
            rxnMdl_to_rxnIntegMdl_map[rxn][0].BiGG_id = list(set(rxnMdl_to_rxnIntegMdl_map[rxn][0].BiGG_id + rxn.BiGG_id)) 

            rxnMdl_to_rxnIntegMdl_map[rxn][0].EC_numbers = list(set(rxnMdl_to_rxnIntegMdl_map[rxn][0].EC_numbers + rxn.EC_numbers))

            if rxn.reversibility != '' and rxnMdl_to_rxnIntegMdl_map[rxn][0].reversibility != '' and rxn.reversibility.lower() != rxnMdl_to_rxnIntegMdl_map[rxn][0].reversibility.lower():
                print "**WARNING! Reversiblity of reaction {} in model {} is '{}' but this does not match the reversbility of its corresponding reaction {} in the integrated model, which is '{}'".format(rxn.id, mdl['model'].id, rxn.reversibility, mdl['model'].id, rxnMdl_to_rxnIntegMdl_map[rxn][0].id, rxnMdl_to_rxnIntegMdl_map[rxn][0].reversibility)
            elif rxn.reversibility != '' and rxnMdl_to_rxnIntegMdl_map[rxn][0].reversibility == '':
                rxnMdl_to_rxnIntegMdl_map[rxn][0].reversibility = rxn.reversibility

            if rxn.gene_reaction_rule != '' and rxnMdl_to_rxnIntegMdl_map[rxn][0].gene_reaction_rule != '' and rxn.gene_reaction_rule.lower() != rxnMdl_to_rxnIntegMdl_map[rxn][0].gene_reaction_rule.lower():
                print "**WARNING! gene_reaction_rule of reaction {} in model {} is '{}' but this does not match gene_reaction_rule of its corresponding reaction {} in the integrated model, which is '{}'".format(rxn.id, mdl['model'].id, rxn.gene_reaction_rule, mdl['model'].id, rxnMdl_to_rxnIntegMdl_map[rxn][0].id, rxnMdl_to_rxnIntegMdl_map[rxn][0].gene_reaction_rule)
            elif rxn.gene_reaction_rule != '' and rxnMdl_to_rxnIntegMdl_map[rxn][0].gene_reaction_rule == '':
                rxnMdl_to_rxnIntegMdl_map[rxn][0].gene_reaction_rule = rxn.gene_reaction_rule

            # Genes
            rxnMdl_to_rxnIntegMdl_map[rxn][0].genes = list(set([g for g in rxn.genes if (mdl['model'].id,g) not in not_included_gens_map.keys()] + [not_included_gens_map[(mdl['model'].id,g)] for g in rxn.genes if (mdl['model'].id,g) in not_included_gens_map.keys()]))

            if rxn.deltaG != None and rxnMdl_to_rxnIntegMdl_map[rxn][0].deltaGt == None:
                rxnMdl_to_rxnIntegMdl_map[rxn][0].deltaG = rxn.deltaG 

            if rxn.deltaG_uncertainty != None and rxnMdl_to_rxnIntegMdl_map[rxn][0].deltaG_uncertaintyt == None:
                rxnMdl_to_rxnIntegMdl_map[rxn][0].deltaG_uncertainty = rxn.deltaG_uncertainty 

            if rxn.deltaG_range != [] and rxnMdl_to_rxnIntegMdl_map[rxn][0].deltaG_ranget == []:
                rxnMdl_to_rxnIntegMdl_map[rxn][0].deltaG_range = rxn.deltaG_range 

        # Consinder only reactions in mdl for which there is no unique map in integrated_model
        for rxn in [c for c in rxnMdl_to_rxnIntegMdl_map if len(rxnMdl_to_rxnIntegMdl_map[rxn]) != 1]:
            rxn.reset_props()
            rxn.compartment = [comparts_by_id[mdl['model_to_standard_compartID_map'][cpt.id]] for cpt in rxn.compartment]

            # Replace any compound in the reaction stoichiometry with a unique match in integrated model
            # with its match
            rxn_stoic = dict(rxn.stoichiometry)
            for cpd in [c for c in rxn.stoichiometry.keys() if (mdl['model'].id,c) in not_included_cpds_map.keys()]:
                rxn_stoic[not_included_cpds_map[(mdl['model'].id,cpd)]] = rxn.stoichiometry[cpd]           
                del rxn_stoic[cpd]
            rxn.set_stoichiometry(self, stoichiometry = rxn_stoic, replace = True)

            reactions.append(rxn)
            reactions_by_id[rxn.id] = rxn


