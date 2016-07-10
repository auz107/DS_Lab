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
from tools.utilities.get_ModelSEED_ids import remove_compartment, get_ModelSEED_ids
from models.ModelSEED.ModelSEED_cpds_master import cpds_master as ModelSEED_cpds
from models.ModelSEED.ModelSEED_rxns_master import rxns_master as ModelSEED_rxns
import itertools

def create_superModel_from_ModelSEED(original_model, standard_to_model_compartID_map, add_c_cpds_to_otherComparts = False, obtain_ModelSEED_ids = True, validate = True,  warnings = True, stdout_msgs = True): 
    """
    INPUTS: 
    ------
    original_model: 
    An instance of object model. The input model should have ModelSEED ids assigned to any 
    compound or reaction if possible 

    standard_to_model_compartID_map: 
    A dictionary where keys are one of the letters below (standard compartment
    ids and values are the corresponding compartment ids in the model
    id in the original model. 
    c: Cytosol (cytoplasm),   e: Extracellular,   g: Golgi,     m: Mitochondria
    n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,    x: Peroxisome
    For example, if a model has two compartments c and e, one can provide 
    {'c':'c id in the model', 'e':'e id in the model'}. One can also provide
    {'c':'', 'e': ''} in which the code searches for these two compartments 
    in the model

    original_model: 
    The original metabolic model

    obtain_ModelSEED_ids: 
    If True ModelSEED ids of reactions and compounds are obtained

    validate: 
    If True, the constructed model undergoes the validation process (this can be very
    time-consuming take 2-3 h)

    add_c_cpds_to_otherComparts: 
    If True, adds a copy of any compound that is in cytosol but not in  model
    compartments other than [p] and [e] along with their transport reactions. Note that a 
    copy of any compound in cytosol that is not present in [p] or [e] compartments or vice
    versa will be always added to the super_model no matter if this parameter is True or False
                    
 
    OUTPUT:
    -------
    super_model: 
    The integrated model containing both the original model and the compounds and reactions 
    in the ModelSEED

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06-20-2016
    """
    start_pt = time.clock()
    start_wt = time.time()

    if not isinstance(original_model,model):
        raise TypeError('Original model must be an instance of class model. An object of type {} was provided instead'.format(type(original_model)))
    if not isinstance(add_c_cpds_to_otherComparts,bool):
        raise TypeError('add_c_cpds_to_otherComparts must be either True or False')

    if not isinstance(standard_to_model_compartID_map,dict):
        raise TypeError('A dicitonary expected for compart_id but an object of type {} was provided instead'.format(type(compart_id)))
    elif not set(standard_to_model_compartID_map.keys()).issubset(['c','e','g','m','n','p','r','x']):
        raise ValueError('Invalid key for standard_to_model_compartID_map. Eligible keys are: c, e, g, m, n, p, r, x')
    elif len([v for v in standard_to_model_compartID_map.values() if not isinstance(v,str)]) > 0:
        raise TypeError('A string expected for the values of dictionary compart_id but objects of type {} were provided instead'.format(list(set([type(v) for v in standard_to_model_compartID_map.values() if not isinstance(v,str)]))))

    if not isinstance(obtain_ModelSEED_ids, bool):
        raise TypeError('obtain_ModelSEED_ids must be either True or False')

    # Reactions for specific types of organisms
    if original_model.organism.ModelSEED_type == '':
        raise userError('ModelSEED_type has not been assigned for the organism object of the original model')
    elif original_model.organism.ModelSEED_type.lower() == 'bacteria_grampositive':
        from models.ModelSEED.ModelSEED_rxns_GramPositive import rxns_GramPositive as organismType_specific_ModelSEED_rxns
    elif original_model.organism.ModelSEED_type.lower() == 'bacteria_gramnegative':
        from models.ModelSEED.ModelSEED_rxns_GramNegative import rxns_GramNegative as organismType_specific_ModelSEED_rxns
    elif original_model.organism.ModelSEED_type.lower() == 'human':
        from models.ModelSEED.ModelSEED_rxns_Human import rxns_Human as organismType_specific_ModelSEED_rxns

    #---- Compartments in the model ----
    # Map from compartment ids in the model to standard ones
    model_to_standard_compartID_map = dict([(standard_to_model_compartID_map[k],k) for k in standard_to_model_compartID_map.keys()])

    #------ Reset model-related attributes of reactions and compartments -----
    # First reset the properties of all compounds and reactions related to the model they belong to 
    # Also assign a new attribute external to each comound and reaction showing if it is an external 
    # reaciton (True) or not (False)
    for cpd in original_model.compounds:
        cpd.reset_props()
        cpd.external = False
    for rxn in original_model.reactions:
        rxn.reset_props()
        rxn.external = False        
    for cpt in original_model.compartments:
        cpt.model = None
    original_model.organism.model = None

    # Get ModelSEED ids
    if obtain_ModelSEED_ids:
        get_ModelSEED_ids(model = original_model, stdout_msgs = False)

    # ModelSEED if of compounds and reactions in the model with a unique ModelSEED id
    orig_model_cpds_ModelSEED_ids = []
    orig_model_rxns_ModelSEED_ids = [mid for r in original_model.reactions if len(r.ModelSEED_id) > 1 for mid in r.ModelSEED_id]

    # ModelSEED ids corresponding to more than one compound in the model
    all_orig_model_cpds_ModelSEED_ids = list(set([mid for c in original_model.compounds if len(c.ModelSEED_id) > 0 for mid in c.ModelSEED_id]))
    ModelSEED_id_cpd_map = dict([(mid,[]) for mid in all_orig_model_cpds_ModelSEED_ids])
    for mid in all_orig_model_cpds_ModelSEED_ids:
        mid_cpds = [c for c in original_model.compounds if mid in c.ModelSEED_id]
        # Compartment id appears in the names of compounds and reactions for models constructed
        # on ModelSEED or KBase, so we remove them from the names to find out whether the 
        # compound appears in more than one compartment based on their name 
        mid_cpds_names = [remove_compartment(input_string = c.name, compartments_info = [c.compartment.id]) for c in mid_cpds]
        mid_cpds_ids_co_cpt = [remove_compartment(input_string = c.id, compartments_info = [c.compartment.id]) for c in mid_cpds]
        ModelSEED_id_cpd_map[mid] = [c.id for c in mid_cpds if mid_cpds_ids_co_cpt.count(remove_compartment(input_string = c.id, compartments_info = [c.compartment.id])) == 1  or mid_cpds_names.count(remove_compartment(input_string = c.name, compartments_info = [c.compartment.id])) == 1]
    ModelSEED_ids_non_unique_cpds = [mid for mid in ModelSEED_id_cpd_map.keys() if len(ModelSEED_id_cpd_map[mid]) > 1]

    # List of compounds and reactions in the super model
    organism = original_model.organism
    cpds_superModel = [] 
    rxns_superModel = []

    # Access reactions in the super_model by id. For compounds or reactions int he model, this id is the
    # ModelSEED id if that reaction or cmpound has at least one ModelSEED id. Otherwise, it is its id in the
    # original model. For compounds or reactions with more than one ModelSEED id we consider the first in the 
    # list assuming that these are equivalent due to having the same name, forumula, stoichioemtry, etc
    cpds_superModel_by_id = {}
    rxns_superModel_by_id = {}

    for cpd in original_model.compounds:
        cpds_superModel.append(cpd)
        if len(cpd.ModelSEED_id) > 0:
           ModelSEED_id_cpt = cpd.ModelSEED_id[0] + '_' + model_to_standard_compartID_map[cpd.compartment.id] + '0'
           if ModelSEED_id_cpt not in cpds_superModel_by_id.keys(): 
                cpds_superModel_by_id[ModelSEED_id_cpt] = [cpd]
           elif cpd.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and ModelSEED_id_cpt in cpds_superModel_by_id.keys():
                cpds_superModel_by_id[ModelSEED_id_cpt].append(cpd)
                orig_model_cpds_ModelSEED_ids.append(ModelSEED_id_cpt)
        else:
            cpds_superModel_by_id[cpd.id] = [cpd]
            orig_model_cpds_ModelSEED_ids.append(cpd.id)

    for rxn in original_model.reactions:
        rxn.base_cost = 0
        if len(rxn.ModelSEED_id) > 1 and not rxn.is_exchange and rxn.ModelSEED_id[0] in organismType_specific_ModelSEED_rxns.keys():
            rxn.forward_cost = organismType_specific_ModelSEED_rxns[rxn.ModelSEED_id[0]]['forward_cost']
            rxn.backward_cost = organismType_specific_ModelSEED_rxns[rxn.ModelSEED_id[0]]['backward_cost']
            rxn.probability_dG_lessThanZero = ModelSEED_rxns[rxn.ModelSEED_id[0]]['probability_dG_lessThanZero']
            rxns_superModel_by_id[rxn.ModelSEED_id[0]] = rxn
        else:
            # Reaction is reversible in the original model
            if rxn.reversibility == 'reversible' and not rxn.is_exchange:
                rxn.forward_cost = 0 
                rxn.backward_cost = 0
            # Reaction is irreversible (forward) in the original model
            if rxn.reversibility == 'irreversible' and not rxn.is_exchange:
                rxn.forward_cost = 0 
                rxn.backward_cost = 5 
            if rxn.is_exchange: # Exchange reactions
                rxn.forward_cost = 0 
                rxn.backward_cost = 0 
            rxn.probability_dG_lessThanZero = None
            rxns_superModel_by_id[rxn.id] = rxn
        rxns_superModel.append(rxn)

    #-------- Modifications from the original model
    # All transport reactions in the ModelSEED database
    ModelSEED_transport_rxns = [r for r in ModelSEED_rxns if ModelSEED_rxns[r]['is_transport']]

    #--- Add exchange reactions to the model for compounds in the [e] compartment with no exchange reaction ---
    for cpd_e in [c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and len([r for r in original_model.reactions if (r.is_exchange or 'exchange' in r.name or 'EX_' in r.name) and c in r.reactants]) == 0]:
        exch_rxn = create_exchrxn(cpd_e) 
        rxns_superModel.append(exch_rxn)
        rxns_superModel_by_id[exch_rxn.id] = exch_rxn

    # --- Add missing transport reactions between copies of the same compound in different comaprtments ---
    # (A) If the model has no [p] compartment, 
    #     Add transport reactions for any compounds present in cytosol and other compartments 
    #     if no such transport reaction already exists in the model
    if 'p' not in standard_to_model_compartID_map.keys():
        # All compartments in the model other than [c]
        for cpt in [ct for ct in standard_to_model_compartID_map.keys() if ct != 'c']: 
            # Compounds in [c] and in cpt 
            for cpd_c in [c_c for c_c in original_model.compounds if c_c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map[cpt]] and c.name == c_c.name]) == 1]:
                # Corresponding compound in other compartment
                cpd_cpt = [c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map[cpt]] and c.name == cpd_c.name][0]
   
                # Add a transport reaction transport reaction between [c] and cpt compartment. if
                # 1. If no transport reaction exists in the model
                # 2. If there is no such transport reaction in ModelSEED when cpt == 'e' 
                # in the ModelSEED
                if len([r for r in original_model.reactions if cpd_c in r.stoichiometry.keys() and cpd_cpt in r.stoichiometry.keys()]) == 0 and (cpt != 'e' or (cpt == 'e' and not (len(cpd_c.ModelSEED_id) > 1 and len([r for r in ModelSEED_transport_rxns if (cpd_c.ModelSEED_id[0],'c0') in ModelSEED_rxns[r]['stoichiometry'] and (cpd_c.ModelSEED_id[0],'e0') in ModelSEED_rxns[r]['stoichiometry']]) > 0))):
                    transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_cpt)
                    if transport_rxn.id not in rxns_superModel_by_id.keys():
                        rxns_superModel.append(transport_rxn)
                        rxns_superModel_by_id[transport_rxn.id] = transport_rxn
    
    # (B) If the model has an [p] compartment,
    #     (1) Add transport reactions for any compound present in cytosol and all other compartments 
    #         other than [e] if no such transport reaction already exists in the model
    #     (2) Add transport reactions for any compound present in [e] and [p] 
    #         if no such transport reaction already exists in the model
    else:
        # (1)
        # All compartments in the model other than [c] and [e]
        for cpt in [ct for ct in standard_to_model_compartID_map.keys() if ct != 'c' and ct != 'e']: 
            # Compounds in [c] and in cpt 
            for cpd_c in [c_c for c_c in original_model.compounds if c_c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map[cpt]] and c.name == c_c.name]) == 1]:
                # Corresponding compound in other compartment
                cpd_cpt = [c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map[cpt]] and c.name == cpd_c.name][0]
                if len([r for r in original_model.reactions if cpd_c in r.stoichiometry.keys() and cpd_cpt in r.stoichiometry.keys()]) == 0: 
                    transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_cpt)
                    if transport_rxn.id not in rxns_superModel_by_id.keys():
                        rxns_superModel.append(transport_rxn)
                        rxns_superModel_by_id[transport_rxn.id] = transport_rxn
    
        # (2)
        for cpd_p in [c_p for c_p in original_model.compounds if c_p.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['p']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and c.name == c_p.name]) == 1]:
            # Corresponding compound in other compartment
            cpd_e = [c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and c.name == c_p.name][0] 
            if len([r for r in original_model.reactions if cpd_p in r.stoichiometry.keys() and cpd_e in r.stoichiometry.keys()]) == 0: 
                transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_p, cpd_cpt2 = cpd_e)
                if transport_rxn.id not in rxns_superModel_by_id.keys():
                    rxns_superModel.append(transport_rxn)
                    rxns_superModel_by_id[transport_rxn.id] = transport_rxn

    # --- Add a secretion pathway for all compounds in the model ---
    # (A) If the model has no [p] compartment, 
    #     For any compounds present in cytosol but not in [e] compartment
    #     (1) Add a new compound in [e] 
    #     (2) Add a transport reaction between [c] and [e] 
    #     (3) Add an exchange reaction 
    #
    #     For any compounds present in [e] but not in [c] compartment
    #     (4) Add a new compound in [c] 
    #     (5) Add a transport reaction between [c] and [e] 
    if 'p' not in standard_to_model_compartID_map.keys():
        # For any compounds present in cytosol but not in [e] compartment
        for cpd_c in [c_c for c_c in original_model.compounds if c_c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and c.name == c_c.name]) == 0]:
            # (1)
            # Create a new compound in compartment cpt
            cpd_e = copy_cpd_new_cpt(cpd = cpd_c, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['e']])

            cpd_e_added = False 
            if len(cpd_e.ModelSEED_id) >= 1:
               if cpd_e.ModelSEED_id[0] + '_e0' not in cpds_superModel_by_id.keys(): 
                    cpds_superModel_by_id[cpd_e.ModelSEED_id[0] + '_e0'] = [cpd_e]
                    cpds_superModel.append(cpd_e)
                    cpd_e_added = True
               elif cpd_e.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_e.ModelSEED_id[0] + '_e0' in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_e.ModelSEED_id[0] + '_e0'].append(cpd_e)
                    cpds_superModel.append(cpd_e)
                    cpd_e_added = True
            elif cpd_e.id not in cpds_superModel_by_id.keys():
                cpds_superModel_by_id[cpd_e.id] = [cpd_e]
                cpds_superModel.append(cpd_e)
                cpd_e_added = True
    
            if cpd_e_added:
                # (2) 
                # Add a transport reaction between [c] and [e] only if there is no such transport reaciton 
                # in the ModelSEED
                if not (len(cpd_c.ModelSEED_id) > 1 and len([r for r in ModelSEED_transport_rxns if (cpd_c.ModelSEED_id[0],'c0') in ModelSEED_rxns[r]['stoichiometry'] and (cpd_c.ModelSEED_id[0],'e0') in ModelSEED_rxns[r]['stoichiometry']]) > 0):
                    transport_rxn = create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_e)
                    if transport_rxn.id not in rxns_superModel_by_id.keys():
                        rxns_superModel.append(transport_rxn)
                        rxns_superModel_by_id[transport_rxn.id] = transport_rxn
        
                # (3)
                exch_rxn = create_exchrxn(cpd_e) 
                rxns_superModel.append(exch_rxn)
                rxns_superModel_by_id[exch_rxn.id] = exch_rxn

        # For any compounds present in [e] but not in [c] compartment
        for cpd_e in [c_e for c_e in original_model.compounds if c_e.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and c.name == c_e.name]) == 0]:
            # (4)
            # Create a new compound in compartment cpt
            cpd_c = copy_cpd_new_cpt(cpd = cpd_e, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['c']])
            
            cpd_c_added = False
            if len(cpd_c.ModelSEED_id) >= 1:
                if cpd_c.ModelSEED_id[0] + '_c0' not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_c.ModelSEED_id[0] + '_c0'] = [cpd_c]
                    cpd_c_added = True
                    cpds_superModel.append(cpd_c)
                elif cpd_c.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_c.ModelSEED_id[0] + '_c0' in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_c.ModelSEED_id[0] + '_c0'].append(cpd_c)
                    cpd_c_added = True
                    cpds_superModel.append(cpd_c)
            elif cpd_c.id not in cpds_superModel_by_id.keys():
                cpds_superModel_by_id[cpd_c.id] = [cpd_c]
                cpds_superModel.append(cpd_c)
                cpd_c_added = True
    
            # (5) 
            # Add a transport reaction between [c] and [e] only if there is no such transport reaciton 
            # in the ModelSEED
            if cpd_c_added and not (len(cpd_c.ModelSEED_id) > 1 and len([r for r in ModelSEED_transport_rxns if (cpd_c.ModelSEED_id[0],'c0') in ModelSEED_rxns[r]['stoichiometry'] and (cpd_c.ModelSEED_id[0],'e0') in ModelSEED_rxns[r]['stoichiometry']]) > 0):
                transport_rxn = create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_e)
                if transport_rxn.id not in rxns_superModel_by_id.keys():
                    rxns_superModel.append(transport_rxn)
                    rxns_superModel_by_id[transport_rxn.id] = transport_rxn

    # (B) If the model has a compartment [p], 
    #     For any compounds present in cytosol but not in [p] compartment
    #     (1) Add a new compound in [p] 
    #     (2) Add a transport reaction between [c] and [p] 
    #
    #     For any compounds present in [p] but not in cytosol 
    #     (3) Add a new compound in [c] 
    #     (4) Add a transport reaction between [c] and [p] 
    #
    #     For any compounds present in [p] but not in [e] compartment
    #     (5) Add a new compound in [e] 
    #     (6) Add a transport reaction between [p] and [e] 
    #     (7) Add an exchange reaction
    #
    #     For any compounds present in [e] but not in [p] compartment
    #     (7) Add a new compound in [p] 
    #     (8) Add a transport reaction between [p] and [e] 
    else:
        # For any compounds present in cytosol but not in [p] compartment
        for cpd_c in [c_c for c_c in original_model.compounds if c_c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['p']] and c.name == c_c.name]) == 0]:

            # (1)
            cpd_p = copy_cpd_new_cpt(cpd = cpd_c, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['p']])

            cpd_p_added = False
            if len(cpd_p.ModelSEED_id) >= 1:
                if cpd_p.ModelSEED_id[0] + '_p0' not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_p.ModelSEED_id[0] + '_p0'] = [cpd_p]
                    cpds_superModel.append(cpd_p)
                    cpd_p_added = True
                elif cpd_p.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_p.ModelSEED_id[0] + '_p0' in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_p.ModelSEED_id[0] + '_p0'].append(cpd_p)
                    cpds_superModel.append(cpd_p)
                    cpd_p_added = True
            elif cpd_p.id not in cpds_superModel_by_id.keys():
                cpds_superModel_by_id[cpd_p.id] = [cpd_p]
                cpds_superModel.append(cpd_p)
                cpd_p_added = True
    
            # (2) 
            if cpd_p_added:
                # Add a transport reaction between [c] and [p] 
                transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_p)
                if transport_rxn.id not in rxns_superModel_by_id.keys():
                    rxns_superModel.append(transport_rxn)
                    rxns_superModel_by_id[transport_rxn.id] = transport_rxn
    
        # For any compounds present in [p] but not in [c] compartment
        for cpd_p in [c_p for c_p in original_model.compounds if c_p.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['p']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and c.name == c_p.name]) == 0]:
            # (3)
            cpd_c = copy_cpd_new_cpt(cpd = cpd_p, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['c']])

            cpd_c_added = False
            if len(cpd_c.ModelSEED_id) >= 1:
                if cpd_c.ModelSEED_id[0] + '_c0' not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_c.ModelSEED_id[0] + '_c0'] = [cpd_c]
                    cpds_superModel.append(cpd_c)
                    cpd_c_added = True
                elif cpd_c.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_c.ModelSEED_id[0] + '_c0' in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_c.ModelSEED_id[0] + '_c0'].append(cpd_c)
                    cpds_superModel.append(cpd_c)
                    cpd_c_added = True
            elif cpd_c.id in cpds_superModel_by_id.keys():
                cpds_superModel_by_id[cpd_c.id] = [cpd_c]
                cpds_superModel.append(cpd_c)
                cpd_c_added = True

            # (4) 
            if cpd_c_added:
                # Add a transport reaction between [c] and [p] 
                transport_rxn = create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_p)
                if transport_rxn.id not in rxns_superModel_by_id.keys():
                    rxns_superModel.append(transport_rxn)
                    rxns_superModel_by_id[transport_rxn.id] = transport_rxn

        # For any compounds present in [p] but not in [e] compartment
        for cpd_p in [c_p for c_p in original_model.compounds if c_p.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['p']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and c.name == c_p.name]) == 0]:
            # (5)
            cpd_e = copy_cpd_new_cpt(cpd = cpd_p, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['e']])

            cpd_e_added = False
            if len(cpd_e.ModelSEED_id) >= 1:
                if cpd_e.ModelSEED_id[0] + '_e0' not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_e.ModelSEED_id[0] + '_e0'] = [cpd_e]
                    cpds_superModel.append(cpd_e)
                    cpd_e_added = True
                elif cpd_e.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_e.ModelSEED_id[0] + '_e0' in cpds_superModel_by_id.keys(): 
                    cpds_superModel_by_id[cpd_e.ModelSEED_id[0] + '_e0'].append(cpd_e)
                    cpds_superModel.append(cpd_e)
                    cpd_e_added = True
            else:
                cpds_superModel_by_id[cpd_e.id] = [cpd_e]
                cpds_superModel.append(cpd_e)
                cpd_e_added = True
    
            if cpd_e_added:
                # (6) 
                # Add a transport reaction between [c] and [p] 
                transport_rxn = create_transport_rxn(cpd_cpt1 = cpd_p, cpd_cpt2 = cpd_e)
                if transport_rxn.id not in rxns_superModel_by_id.keys():
                    rxns_superModel.append(transport_rxn)
                    rxns_superModel_by_id[transport_rxn.id] = transport_rxn
    
                # (7)
                exch_rxn = create_exchrxn(cpd_e) 
                rxns_superModel.append(exch_rxn)
                rxns_superModel_by_id[exch_rxn.id] = exch_rxn

        # For any compounds present in [e] but not in [p] compartment
        for cpd_e in [c_e for c_e in original_model.compounds if c_e.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['e']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['p']] and c.name == c_e.name]) == 0]:
            # (8)
            # Create a new compound in compartment cpt
            cpd_p = copy_cpd_new_cpt(cpd = cpd_e, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['p']])

            cpd_p_added = False
   
            if len(cpd_p.ModelSEED_id) >= 1:
                if cpd_p.ModelSEED_id[0] + '_p0' not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_p.ModelSEED_id[0] + '_p0'] = [cpd_p]
                    cpds_superModel.append(cpd_p)
                    cpd_p_added = True
                elif cpd_p.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_p.ModelSEED_id[0] + '_p0' in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_p.ModelSEED_id[0] + '_p0'].append(cpd_p)
                    cpds_superModel.append(cpd_p)
                    cpd_p_added = True
            elif cpd_p.id not in cpds_superModel_by_id.keys():
                cpds_superModel_by_id[cpd_p.id] = [cpd_p]
                cpds_superModel.append(cpd_p)
                cpd_p_added = True
  
            # (9) 
            if cpd_p_added:
                transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_p, cpd_cpt2 = cpd_e)
                if transport_rxn.id not in rxns_superModel_by_id.keys():
                    rxns_superModel.append(transport_rxn)
                    rxns_superModel_by_id[transport_rxn.id] = transport_rxn

    #--- Add copies of compounds in cytosol missing in other comaprtments and their transport reactions ---
    if add_c_cpds_to_otherComparts:
        # For any compound in the cytosol for which there is not any corresponding compound in all other compartments
        # other than [e] (and [p] if the model has compartment [p]) 
        # 1. Add a new corresponding compound in that compartment 
        # 2. Add a transport reaction between cytosol and that compartment
        #
        # For any compound in any model compartment than [e] (and [p], if the model has a [p] compartment) for which there 
        # is not corresponding compound in [c] 
        # 3. Add a new corresponding compound [c] 
        # 4. Add a transport reaction between cytosol and that compartment
        if 'p' not in standard_to_model_compartID_map.keys():
            excluded_cpts = ['e'] 
        else:
            excluded_cpts = ['e','p']

        # All compartments in the model other than [c] and excluded compounds
        for cpt in [ct for ct in standard_to_model_compartID_map.keys() if ct not in ['c'] + excluded_cpts]: 
            # Compounds in [c] but not in cpt 
            for cpd_c in [c_c for c_c in original_model.compounds if c_c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map[cpt]] and c.name == c_c.name]) == 0]:
                # (1)
                cpd_cpt = copy_cpd_new_cpt(cpd = cpd_c, cpt = original_model.compartments_by_id[standard_to_model_compartID_map[cpt]])

                cpd_cpt_added = False
                if len(cpd_cpt.ModelSEED_id) >= 1:
                    if cpd_cpt.ModelSEED_id[0] + '_' + cpt + '0' not in cpds_superModel_by_id.keys():
                        cpds_superModel_by_id[cpd_cpt.ModelSEED_id[0] + '_' + cpt + '0'] = [cpd_cpt]
                        cpds_superModel.append(cpd_cpt)
                        cpd_cpt_added = True
                    elif cpd_cpt.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_cpt.ModelSEED_id[0] + '_' + cpt + '0' in cpds_superModel_by_id.keys(): 
                        cpds_superModel_by_id[cpd_cpt.ModelSEED_id[0] + '_' + cpt + '0'].append(cpd_cpt)
                        cpds_superModel.append(cpd_cpt)
                        cpd_cpt_added = True
                elif cpd_cpt.id not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_cpt.id] = [cpd_cpt]
                    cpds_superModel.append(cpd_cpt)
                    cpd_cpt_added = True
    
                # (2) 
                if cpd_cpt_added:
                    transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_cpt)
                    if transport_rxn.id not in rxns_superModel_by_id.keys():
                        rxns_superModel.append(transport_rxn)
                        rxns_superModel_by_id[transport_rxn.id] = transport_rxn
    
            # Compounds in cpt but not in [c] 
            for cpd_cpt in [c_cpt for c_cpt in original_model.compounds if c_cpt.compartment == original_model.compartments_by_id[standard_to_model_compartID_map[cpt]] and len([c for c in original_model.compounds if c.compartment == original_model.compartments_by_id[standard_to_model_compartID_map['c']] and c.name == c_cpt.name]) == 0]:
                # (3)
                cpd_c = copy_cpd_new_cpt(cpd = cpd_cpt, cpt = original_model.compartments_by_id[standard_to_model_compartID_map['c']])

                cpd_c_added = False                
                if len(cpd_c.ModelSEED_id) >= 1:
                    if cpd_c.ModelSEED_id[0] + '_c0' not in cpds_superModel_by_id.keys():
                        cpds_superModel_by_id[cpd_c.ModelSEED_id[0] + '_c0'] = [cpd_c]
                        cpds_superModel.append(cpd_c)
                        cpd_c_added = True
                    elif cpd_c.ModelSEED_id[0] in ModelSEED_ids_non_unique_cpds and cpd_c.ModelSEED_id[0] + '_c0' in cpds_superModel_by_id.keys():
                        cpds_superModel_by_id[cpd_c.ModelSEED_id[0] + '_c0'].append(cpd_c)
                        cpds_superModel.append(cpd_c)
                        cpd_c_added = True
                elif cpd_c.id not in cpds_superModel_by_id.keys():
                    cpds_superModel_by_id[cpd_c.id] = [cpd_c]
                    cpds_superModel.append(cpd_c)
                    cpd_c_added = True
  
                # (4) 
                if cpd_c_added:
                    transport_rxn =  create_transport_rxn(cpd_cpt1 = cpd_c, cpd_cpt2 = cpd_cpt)
                    if transport_rxn.id not in rxns_superModel_by_id.keys():
                        rxns_superModel.append(transport_rxn)
                        rxns_superModel_by_id[transport_rxn.id] = transport_rxn
    

    #--- ModelSEED compounds ---
    # Exchange reactions in the ModelSEED
    ModelSEED_exchrxns = [r for r in ModelSEED_rxns if 'exchange' in ModelSEED_rxns[r]['name'].lower() or len([n for n in ModelSEED_rxns[r]['name_aliases'] if 'exchange' in n]) > 0 or (ModelSEED_rxns[r]['BiGG_id'] != None and  'EX_' in ModelSEED_rxns[r]['BiGG_id'])]

    # List of compounds in the ModelSEED participating in the ModelSEED reactions
    # ** Consider only those in c0 and e0 compartments **
    cpds_in_ModelSEED_rxns = list(set([(cpd,cpt) for r in ModelSEED_rxns.keys() if ModelSEED_rxns[r]['stoichiometry'] != None and r not in ModelSEED_exchrxns + orig_model_rxns_ModelSEED_ids and 'a0' not in [cpd_cpt[1] for cpd_cpt in ModelSEED_rxns[r]['stoichiometry'].keys()] for (cpd,cpt) in ModelSEED_rxns[r]['stoichiometry'].keys() if cpt in ['c0','e0']]))

    # Compounds in the ModelSEED except those in the original model. Also excludes compounds that do not participate
    # in any ModelSEED reactions
    for (cid,cpt) in [(c,ct) for (c,ct) in cpds_in_ModelSEED_rxns if c + '_' + ct not in orig_model_cpds_ModelSEED_ids and c + '_' + ct not in cpds_superModel_by_id.keys()]:
        if cpt == 'c0':
           compart = original_model.compartments_by_id[standard_to_model_compartID_map['c']]
        elif cpt == 'e0':
            compart = original_model.compartments_by_id[standard_to_model_compartID_map['e']]
        else:
            raise userError('Unknonw compartment: {} for compound {}'.format(cpt,cid))

        ModelSEED_cpd = compound(id = cid + '_' + cpt, compartment = compart, name = ModelSEED_cpds[cid]['name'], name_aliases = ModelSEED_cpds[cid]['name_aliases'], KEGG_id = ModelSEED_cpds[cid]['KEGG_id'], ModelSEED_id = cid, BiGG_id = ModelSEED_cpds[cid]['BiGG_id'], formula = ModelSEED_cpds[cid]['formula'], deltaG = ModelSEED_cpds[cid]['deltaG'], deltaG_uncertainty = ModelSEED_cpds[cid]['deltaG_uncertainty'])

        ModelSEED_cpd.external = True
        ModelSEED_cpd.external_type = 'ModelSEED compound'

        cpds_superModel.append(ModelSEED_cpd)
        cpds_superModel_by_id[ModelSEED_cpd.id] = [ModelSEED_cpd]


    #--- ModelSEED reactions ---
    # Maximum cost of reactions whose template_rxn_type is NOT conditional
    max_cost_nonCond = max([organismType_specific_ModelSEED_rxns[rid]['base_cost'] + max(organismType_specific_ModelSEED_rxns[rid]['forward_cost'], organismType_specific_ModelSEED_rxns[rid]['backward_cost']) for rid in organismType_specific_ModelSEED_rxns.keys() if organismType_specific_ModelSEED_rxns[rid]['template_rxn_type'] != 'conditional']) 

    # Maximum cost of reactions whose template_rxn_type is conditional
    max_cost_cond = max([organismType_specific_ModelSEED_rxns[rid]['base_cost'] + max(organismType_specific_ModelSEED_rxns[rid]['forward_cost'], organismType_specific_ModelSEED_rxns[rid]['backward_cost']) for rid in organismType_specific_ModelSEED_rxns.keys() if organismType_specific_ModelSEED_rxns[rid]['template_rxn_type'] == 'conditional']) 

    counter = 0

    if stdout_msgs:
        print 'create_superModel: Started reading reactions ...',

    # Do not consider reactions in the original model and those having a compartment 'a0', which is related to 
    # template Human (the only compartments appearing in gram positive/negative are c0 and e0)
    for rid in [r for r in ModelSEED_rxns.keys() if ModelSEED_rxns[r]['stoichiometry'] != None and r not in ModelSEED_exchrxns + orig_model_rxns_ModelSEED_ids and 'a0' not in [cpd_cpt[1] for cpd_cpt in ModelSEED_rxns[r]['stoichiometry'].keys()] and r not in ModelSEED_exchrxns]:
        counter += 1
        if stdout_msgs: 
            if counter/5000 == int(counter/5000):
                print '{} '.format(counter),
            sys.stdout.flush()

        # Reversibility
        if rid in organismType_specific_ModelSEED_rxns and organismType_specific_ModelSEED_rxns[rid]['reversibility_ModelSEED_curated_template'].lower() != 'unknown':
            r_rev = organismType_specific_ModelSEED_rxns[rid]['reversibility_ModelSEED_curated_template'] 
        else:
            r_rev = ModelSEED_rxns[rid]['reversibility_ModelSEED_curated_master']

        # If an organism-type specific reaction 
        if rid in organismType_specific_ModelSEED_rxns.keys():
            external_type = 'ModelSEED template rxn'
 
            # Base cost of adding this reaction
            if organismType_specific_ModelSEED_rxns[rid]['template_rxn_type'] != 'conditional':
                base_cost = organismType_specific_ModelSEED_rxns[rid]['base_cost'] 
            else:
                # Add max + 1 to the base_cost of conditional reactions. This is to make sure that the cost of adding 
                # non-conditoinal reactions is lower than those of conditional ones
                base_cost = organismType_specific_ModelSEED_rxns[rid]['base_cost'] + (max_cost_nonCond + 1) 

            # forward and backward costs of adding this reaction
            forward_cost = organismType_specific_ModelSEED_rxns[rid]['forward_cost'] 
            backward_cost = organismType_specific_ModelSEED_rxns[rid]['backward_cost'] 

            # template reaction type
            template_rxn_type = organismType_specific_ModelSEED_rxns[rid]['template_rxn_type'] 

        else: 
            external_type = 'ModelSEED non-template rxn'
            template_rxn_type = None
 
            # The cost of adding reactions not in the template should be higher than the maximum cost of condiitonal
            # reactions in the template
            base_cost = max_cost_cond + (max_cost_nonCond + 1) + 1

            # forward and backward costs are not defined for reactions not in the template, i.e., for these reactions
            # we cosider only the cost of adding them to the model in the direction specified by reversiblity 
            forward_cost = None 
            backward_cost = None 
        
        probability_dG_lessThanZero = ModelSEED_rxns[rid]['probability_dG_lessThanZero'] 

        # Reaction compartments
        r_compart = []
        for (cpd,cpt) in ModelSEED_rxns[rid]['stoichiometry'].keys():
            if cpt == 'c0':
                r_compart.append(original_model.compartments_by_id[standard_to_model_compartID_map['c']])
            elif cpt == 'e0':
                r_compart.append(original_model.compartments_by_id[standard_to_model_compartID_map['e']])
            else:
                raise userError('Unknown compartment: {} for compound {} in reaction {}'.format(cpt,cpd, rid))


        #-- Reaction stoichiometry and compartments --
        # The following is a list of lists where each inner list is a list of tuples composed of the compound 
        # objects corresponding to a compound ModelSEED id and their stoichiometric coefficient. Example,
        # reaction stoichiometry: stoic = {'a':-1,'b':2,'c':3} 
        # cpds_superModel_by_id = {'a':['a1','a2'],'b':['b1','b2'], 'c':['c']}
        # r_stoic_cpds_list = [[('a1', -1), ('a2', -1)], [('c', 3)], [('b1', 2), ('b2', 2)]]
        # Except that 'a1', 'a2', 'b1', 'b2' and 'c' are replaced by compound objects on our case
        r_stoic_cpds_list = [[(cpd_obj,s_coeff) for cpd_obj in cpds_superModel_by_id[cpd + '_' + cpt]] for ((cpd,cpt),s_coeff) in ModelSEED_rxns[rid]['stoichiometry'].items()]

        # The folloiwng is a list of reaction stoichioemtries. For the example above
        # r_all_stoic = [[('a1', -1), ('c', 3), ('b1', 2)],
        #                [('a1', -1), ('c', 3), ('b2', 2)]
        #                [('a2', -1), ('c', 3), ('b1', 2)]
        #                [('a2', -1), ('c', 3), ('b2', 2)]]
        r_all_stoic = [list(a) for a in itertools.product(*r_stoic_cpds_list)] 
      
        if len(r_all_stoic) > 1:
            # Ids of compounds participating in ModelSEED ids in the reaction stoichiometry corresponding to more
            # than one compound object. For the example above this is ['a1','a2','b1','b2']
            duplated_cpd_ids = [cpd_obj.id for cpd_obj in [c for (cpd,cpt) in ModelSEED_rxns[rid]['stoichiometry'].keys() if len(cpds_superModel_by_id[cpd + '_' + cpt]) > 1 for c in cpds_superModel_by_id[cpd + '_' + cpt]]] 

        for r_stoic in r_all_stoic: 
            r_stoic = dict(r_stoic)

            if len(r_all_stoic) == 1:
                r_id = rid
            else:
                r_id = rid + '_' + '_'.join([c_id for c_id in [c.id for c in r_stoic.keys()] if c_id in duplated_cpd_ids]) 

            # Reaction 
            ModelSEED_rxn = reaction(id = r_id, stoichiometry = r_stoic, reversibility = r_rev, name = ModelSEED_rxns[rid]['name'], name_aliases = [], KEGG_id = ModelSEED_rxns[rid]['KEGG_id'], ModelSEED_id = rid, BiGG_id = ModelSEED_rxns[rid]['BiGG_id'], EC_numbers = [], subsystem = '', pathways = [], genes = [], gene_reaction_rule = '', deltaG = ModelSEED_rxns[rid]['deltaG_ModelSEED'], deltaG_uncertainty = ModelSEED_rxns[rid]['deltaG_uncertainty_ModelSEED'], deltaG_range = [ModelSEED_rxns[rid]['deltaG_min_ModelSEED'],ModelSEED_rxns[rid]['deltaG_max_ModelSEED']])

            ModelSEED_rxn.external = True
            ModelSEED_rxn.external_type = external_type
            ModelSEED_rxn.base_cost = base_cost 
            ModelSEED_rxn.forward_cost = forward_cost 
            ModelSEED_rxn.backward_cost = backward_cost 
            ModelSEED_rxn.template_rxn_type = template_rxn_type
            ModelSEED_rxn_probability_dG_lessThanZero = probability_dG_lessThanZero 

            rxns_superModel.append(ModelSEED_rxn)
            rxns_superModel_by_id[ModelSEED_rxn.id] = ModelSEED_rxn

    if stdout_msgs:
        print '{}\n'.format(counter)
        sys.stdout.flush()

    #---------- Create the model -----------
    super_model =  model(id = 'super_model', type = 'metabolic', organism = original_model.organism, reactions = rxns_superModel, compounds = cpds_superModel, compartments = [cpt for cpt in original_model.compartments], validate = validate, warnings = warnings, stdout_msgs = stdout_msgs)

    elapsed_pt = str(timedelta(seconds = time.clock() - start_pt))
    elapsed_wt = str(timedelta(seconds = time.time() - start_wt))
    if stdout_msgs:
        print 'superModel: Super model was created. It took {}/{} of professing/wall time '.format(elapsed_pt,elapsed_wt)

    return super_model


def create_transport_rxn(cpd_cpt1, cpd_cpt2):
    """
    Creates a transport reaction between two compounds 

    INPUTS:
    -------
    cpd_cpt1: The compound object in compartment 1
    cpd_cpt2: The compound object in compartment 2
    """
    # Id with compartment name removed
    no_cpt_id = remove_compartment(input_string = cpd_cpt1.id, compartments_info = [cpd_cpt1.compartment.id])

    transport_rxn = reaction(id = no_cpt_id + '_' + '_'.join(sorted([cpd_cpt1.compartment.id, cpd_cpt2.compartment.id])) + '_transport', stoichiometry = {cpd_cpt1:-1, cpd_cpt2:1}, reversibility = 'reversible', name = cpd_cpt1.name + ' transport')
    transport_rxn.base_cost = 0 
    transport_rxn.forward_cost = 0 
    transport_rxn.backward_cost = 0 
    transport_rxn.probability_dG_lessThanZero = None
        
    transport_rxn.external = True 
    transport_rxn.external_type = 'transport'
    
    return transport_rxn 

def create_exchrxn(cpd_e):
    """
    Creates an exchange reaction for a compound in the e compartment

    INPUTS:
    ------
    cpd_e: Compound object in the [e] compartment
    """
    # Id with compartment name removed
    no_cpt_id = remove_compartment(input_string = cpd_e.id, compartments_info = [cpd_e.compartment.id])

    exch_rxn = reaction(id = 'EX_' + no_cpt_id + '_e0', stoichiometry = {cpd_e:-1}, reversibility = 'reversible', is_exchange = True, name = cpd_e.name + ' exchange', ModelSEED_id = ['EX_' + mid + '_e0' for mid in cpd_e.ModelSEED_id], KEGG_id = ['EX_' + mid + '_e0' for mid in cpd_e.KEGG_id], BiGG_id = ['EX_' + mid + '_e0' for mid in cpd_e.BiGG_id])
    exch_rxn.base_cost = 0
    exch_rxn.forward_cost = 0 
    exch_rxn.backward_cost = 0 
    exch_rxn.probability_dG_lessThanZero = None
    
    exch_rxn.external = True 
    exch_rxn.external_type = 'exchange'

    return exch_rxn   

def copy_cpd_new_cpt(cpd, cpt):
    """
    Copies a compound object to a new compartment

    INPUTS:
    ------
    cpd: Compound object 
    cpt: The compartment object to which cpd should be copied
    """
    # Id with compartment name removed
    no_cpt_id = remove_compartment(input_string = cpd.id, compartments_info = [cpd.compartment.id])
    
    # Create a new compound in compartment cpt
    cpd_cpt = compound(id = no_cpt_id + '_' + cpt.id, compartment = cpt, name = cpd.name, name_aliases = cpd.name_aliases, KEGG_id = cpd.KEGG_id, ModelSEED_id = cpd.ModelSEED_id, BiGG_id = cpd.BiGG_id, formula = cpd.formula, deltaG = cpd.deltaG, deltaG_uncertainty = cpd.deltaG_uncertainty)
    cpd_cpt.external = True
    cpd_cpt.external_type = 'original model compound in new comaprtment'

    return cpd_cpt


 
