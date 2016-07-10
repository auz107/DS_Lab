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
from tools.ancillary.get_ModelSEED_ids import remove_compartment
from models.ModelSEED.ModelSEED_cpds_master import cpds_master as ModelSEED_cpds
from models.ModelSEED.ModelSEED_rxns_master import rxns_master as ModelSEED_rxns
import itertools

def add_ATPM_rxn(model, ATPM_cpds_ids = {'h2o_c':'', 'atp_c':'', 'adp_c':'', 'h_c':'', 'pi_c':''}):
    """
    Add ATPM (ATP maintenance) reaction to an input model 

    INPUTS:
    -------
            model: The input metabolic model
    ATPM_cpds_ids: The id of compounds participating in the ATPM reaction in the model
    """
    #- Add ATP maintenance reaction to the model ATPM   [c] : atp + h2o --> adp + h + pi - 
    # atp: cpd00002   h2o: cpd00001   adp: cpd00008   h: cpd00067    pi: cpd00012
    atp_c = model.compounds_by_id[ATPM_cpds_ids['atp_c']] 
    h2o_c = model.compounds_by_id[ATPM_cpds_ids['h2o_c']]
    adp_c = model.compounds_by_id[ATPM_cpds_ids['adp_c']]
    h_c = model.compounds_by_id[ATPM_cpds_ids['h_c']]
    pi_c = model.compounds_by_id[ATPM_cpds_ids['pi_c']]

    ATPM = reaction(id = 'ATPM', name = 'Non-growth associated maintenance ATP', stoichiometry = {atp_c:-1, h2o_c:-1, adp_c:1, h_c:1, pi_c:1}) 

    # Export the updated model to a new SBML file
    model.add_reactions([ATPM])

    return model 

def add_biomass_rxn(model, biomass_filename, standard_to_biomassModel_cpt_ids_map, standard_to_model_cpt_ids_map, biomass_file_format = 'sbml', biomass_rxn_id = '', obtain_ModelSEED_ids = True, warnings = True, stdout_msgs = True):
    """
    Adds a biomass reaction to a metabolic model

    INPUTS:
    -------
                  model: Metabolic model to which the biomass reaction should be added
       biomass_filename: Name (and path) of the file containing the biomass reaction
    standard_to_biomassModel_cpt_ids_map:
           standard_to_model_cpt_ids_map:
                          A dictionary where keys are standard compartment ids as follows:
                          c: Cytosol (cytoplasm),   e: Extracellular,   g: Golgi,     m: Mitochondria
                          n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,    x: Peroxisome
                          and values are corresponding compartment ids in biomass_model and model

    biomass_file_format: The format of the file holding the biomass inforation (sbml, pydict, pickle, json)
         biomass_rxn_id: Id of the biomass reaction in the input file. If no input is provided, the first reaction
                         in the model is taken as the biomass reaction
    obtain_ModelSEED_ids: If True, ModelSEED ids are obtained for model and biomass_model
    """
    if not isinstance(filename,str):
        raise TypeError('Invalid format for filename. filename must be a string')
    if not isinstance(biomass_rxn_id,str):
        raise TypeError('Invalid format for biomass_rxn_id. biomass_rxn_id must be a string')
    if file_format.lower() not in ['sbml', 'pydic', 'pickle', 'json']: 
        raise userError('Invalid input for file_format. Allowed choices are: sbml, pydict, pickle and json')

    # Import biomass reaction as a model containing the biomass reaction and all participating compounds
    biomass_model = create_model(model_info = {'id':'biomass_model', 'file_format':biomass_file_format, 'model_filename':biomass_filename}, validate = True, stdout_msgs = True, warnings = True)

    if biomass_rxn_id == '':
        biomass_rxn = model.reactions[0]
    else:
        biomass_rxn = biomass_model.reactions_by_id[biomass_rxn_id]

    # Remove any reactions other than the biomass and compounds except for those participating in the
    # biomass reaction from biomass_model
    biomass_model.del_reactions([r for r in biomass_model.reactions if r != biomass_rxn])
    biomass_model.del_compounds([c for c in biomass_model.compounds if c not in biomass_rxn.compounds])

    biomass_rxn.reset_props()
    for cpd in biomass_model.compounds:
        cpd.reset_props()

    # Get ModelSEED ids for compounds and reactions in biomass_model
    get_ModelSEED_ids(biomass_model)

    # Compare biomass model with original model to find common compounds 
    cpdBiomass_to_cpdModel_map, cpdModel_to_cpdBiomass_map = compare_compounds(cpds_list1i = biomass_model1.compounds, cpds_list2 = model.compounds, standard_to_model1_cpt_ids_map = standard_to_biomassModel_cpt_ids_map, standard_to_model2_cpt_ids_map = standard_to_model_cpt_ids_map, cpds_list1_id = biomass_model.id, cpds_list2_id = model.id, warnings = warnings, stdout_msgs = stdout_msgs)[0:2] 

    # Do the following just for testing compare_models
    cpdBiomass_to_cpdModel_map, cpdModel_to_cpdBiomass_map, rxnBiomass_to_rxnModel_map, rxnModel_to_rxnBiomass_map = compare_models(model1 = biomass_model, model2 = model, standard_to_model1_cpt_ids_map = standard_to_biomassModel_cpt_ids_map, standard_to_model2_cpt_ids_map = standard_to_model_cpt_ids_map, obtain_ModelSEED_ids = obtain_ModelSEED_ids, warnings = False, stdout_msgs = False)

    # Now add the biomass reaction and all participating compounds to the model
    




