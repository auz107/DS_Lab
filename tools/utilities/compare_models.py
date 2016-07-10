from __future__ import division
import sys, time
sys.path.append('../../')
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
from get_ModelSEED_ids import get_ModelSEED_ids, remove_compartment
from remove_non_alphanumeric import remove_non_alphanumeric

#--------------------------------------
# Ali R. Zomorrodi - Segre lab @ BU
# Last updated: 04-28-2016
#--------------------------------------

#--------------------------------------
#--------- Ancillary functions --------
#--------------------------------------
def create_ModelSEED_id_map(obj_list):
    """ 
    Creates a dictionary mapping the ModelSEED ids for componds or reactions to objects of
    type compound or reaction in a list. This is applied to objects with a ModelSEED id 

    INPUTS:
    -------
    obj_list: A list of objects of type compounds or reactions

    OUTPUTS:
    -------
    ModelSEED_id_map: A dictionary where keys are ModelSEED ids and values are 
                      objects of type compounds or reactions   
    """ 
    ModelSEED_id_map = {}
    for obj in [o for o in obj_list if o.ModelSEED_id != []]:
        for sid in obj.ModelSEED_id:
            if sid in ModelSEED_id_map.keys():
                ModelSEED_id_map[sid] += [obj]
            else:
                ModelSEED_id_map[sid] = [obj]
    return  ModelSEED_id_map

def create_KEGG_id_map(obj_list):
    """ 
    Creates a dictionary mapping the KEGG ids for componds or reactions to objects of
    type compound or reaction in a list. This is applied to objects with a KEGG id 

    INPUTS:
    -------
    obj_list: A list of objects of type compound or reaction

    OUTPUTS:
    -------
    KEGG_id_map: A dictionary where keys are KEGG ids and values are 
                 objects of type compound or reaction   
    """ 
    KEGG_id_map = {}
    for obj in [o for o in obj_list if o.KEGG_id != []]:
        for kid in obj.KEGG_id:
            if kid in KEGG_id_map.keys():
                KEGG_id_map[kid] += [obj]
            else:
                KEGG_id_map[kid] = [obj]
    return KEGG_id_map

def create_BiGG_id_map(obj_list):
    """ 
    Creates a dictionary mapping the BiGG ids for componds or reactions to objects of
    type compound or reaction in a list. This is applied to objects with a BiGG id 

    INPUTS:
    -------
    obj_list: A list of objects of type compounds or reactions

    OUTPUTS:
    -------
    BiGG_id_map: A dictionary where keys are BiGG ids and values are 
                      objects of type compounds or reactions   
    """ 
    BiGG_id_map = {}
    for obj in [o for o in obj_list if o.BiGG_id != []]:
        for bid in obj.BiGG_id:
            if bid in BiGG_id_map.keys():
                BiGG_id_map[bid] += [obj]
            else:
                BiGG_id_map[bid] = [obj]
    return  BiGG_id_map

def create_cleanNames_map(obj_list,compart_list):
    """ 
    Creates a dictionary mapping the non-alphanumeric names of componds or reactions to 
    objects of type compound or reaction in a list

    INPUTS:
    -------
        obj_list: A list of objects of type compound or reaction
    compart_list: A list of strings containing the ids of compartments for obj_list 

    OUTPUTS:
    -------
    clean_names_map: A dictionary where keys are non-alphanumeric names and values are 
    """ 
    clean_names_map = {}
    for obj in obj_list:
        for name in [obj.name] + obj.name_aliases:
            clean_name = remove_non_alphanumeric(remove_compartment(name,compartments_info = compart_list)).lower()
            if clean_name in clean_names_map.keys():
                clean_names_map[clean_name] += [obj]
            else:
                clean_names_map[clean_name] = [obj]
    # Names of objects in model2 where all non-alpahnumeric characters are removed 
    return clean_names_map

def create_cpd_formulas_map(cpd_list):
    """ 
    Create a dictionary mapping the non-alphanumeric formula of the compounds to  
    compound objects in a given list of compounds. This is applied only to objects
    with a given formula

    INPUTS:
    -------
    cpd_list: A list of objects of type compound

    OUTPUTS:
    -------
    clean_formulas_map: A dictionary where keys are non-alphanumeric names and values are 
    """ 
    formulas_map = {}
    for cpd in [c for c in cpd_list if c.formula != None]:
        if cpd.formula in formulas_map.keys():
            formulas_map[cpd.formula)] += [cpd]
        else:
            formulas_map[cpd.formula] = [cpd]
    return formulas_map

#--------------------------------------
#--------- Main functions --------
#--------------------------------------
def compare_compounds(cpds_list1,cpds_list2, standard_to_cpds_list1_cpt_ids_map, standard_to_cpds_list2_cpt_ids_map, cpds_list1_id = None, cpds_list2_id = None, warnings = True, stdout_msgs = True):
    """
    Find common compounds in two given list of compounds

    INPUTS:
    -------
                cpds_list1: An instance of object compound or a list of such objects
                cpds_list2: An instance of object compound or a list of such objects
             cpds_list1_id: An id for cpds_list1. This is requried if not all compounds in 
                            cpds_list1 belong to the same model or if no model is assigned 
                            to all reactions of list 1
             cpds_list2_id: An id for cpds_list2. This is requried if not all compounds in 
                            cpds_list1 belong to the same model or if no model is assigned 
                            to all reactions of list 2
    standard_to_cpds_list1_cpt_ids_map:
    standard_to_cpds_list2_cpt_ids_map: 
                            A dictionary where keys are standard compartment ids as follows:
                            c: Cytosol (cytoplasm),   e: Extracellular,   g: Golgi,     m: Mitochondria
                            n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,    x: Peroxisome
                            and values are corresponding compartment ids in cpds_list1 and cpds_list2

    OUTPUTS:
    --------
    cpd1_to_cpd2_map: A dictionary where the keys are the compounds objects from cpds_list1 and values 
                      are a list of compound objects corresponding to this compound from cpds_list2 
    cpd2_to_cpd1_map: A dictionary where the keys are the compounds objects from cpds_list2 and values 
                      are a list of compound objects corresponding to this compound from cpds_list1 
    cpd1_to_cpd2_noCopart_map:
    cpd2_to_cpd1_noCopart_map:
                      Same as cpd1_to_cpd2_map and cpd2_to_cpd1_map except that the mataches here are found
                      without considering the compartments
    cpd1_to_cpd2_match_found_by:
    cpd2_to_cpd1_match_found_by:
                      A dictionary with keys being ids of compounds in cpds_list1 or cpds_list2 and values being a 
                      string showing how the matches in cpd1_to_cpd2_map and cpd2_to_cpd1_map were found 
    cpd1_to_cpd2_noCompart_match_found_by:
    cpd2_to_cpd1_noCompart_match_found_by:
                      Same cpd1_to_cpd2_match_found_by and cpd2_to_cpd1_match_found_by for cpd1_to_cpd2_noCopart_map 
                      and cpd2_to_cpd1_noCopart_map
    """
    if not isinstance(strdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be True or False')

    if cpds_list1_id == None:
        model1_id = list(set([c.model.id for c in cpds_list1 if 'model' in dir(c)]))
        if len(model1_id) == 1:
            cpds_list1_id = model1_id[0]
        else:
            raise userError('cpds_list1_id is needed because not all compounds in cpds_list1 have the same model id (' + str(model1_id) + ')')

    if cpds_list2_id == None:
        model2_id = list(set([c.model.id for c in cpds_list2 if 'model' in dir(c)]))
        if len(model2_id) == 1:
            cpds_list2_id = model2_id[0]
        else:
            raise userError('cpds_list2_id is needed because not all compounds in cpds_list2 have the same model id (' + str(model2_id) + ')')

    # Dictionaries where keys are compartment ids in cpds_list1 and cpds_list2 and values are corresponding standard 
    # comaprtment ids
    cpds_list1_cpt_ids_to_standard_map = {}
    cpds_list2_cpt_ids_to_standard_map = {}
    for cpt_id in standard_to_cpds_list1_cpt_ids_map.values():
        cpds_list1_cpt_ids_to_standard_map[cpt_id] = [sid for sid in standard_to_cpds_list1_cpt_ids_map.keys() if standard_to_cpds_list1_cpt_ids_map[sid] == cpt_id][0]
    for cpt_id in standard_to_cpds_list2_cpt_ids_map.values():
        cpds_list2_cpt_ids_to_standard_map[cpt_id] = [sid for sid in standard_to_cpds_list2_cpt_ids_map.keys() if standard_to_cpds_list2_cpt_ids_map[sid] == cpt_id][0]


    # A dictionary mapping the ids of componds in cpds_list2 to compound objects
    cpds_list2_id_map = dict([(c.id,c) for c in cpds_list2]) 
    # ModelSEED ids for compounds in model2
    cpds_list2_ids = cpds_list2_id_map.keys()

    # A dictionary mapping the ModelSEED ids of componds to compounds in cpds_list2
    cpds_list2_ModelSEED_id_map = create_ModelSEED_id_map(cpds_list2)
    # ModelSEED ids for compounds in model2
    cpds_list2_ModelSEED_ids = cpds_list2_ModelSEED_id_map.keys()

    # A dictionary mapping the KEGG ids of componds to compounds in model2
    cpds_list2_KEGG_id_map = create_KEGG_id_map(cpds_list2)
    # ModelSEED ids for compounds in model2
    cpds_list2_KEGG_ids = cpds_list2_KEGG_id_map.keys()

    # A dictionary mapping the BiGG ids of componds to compounds in model2
    cpds_list2_BiGG_id_map = create_BiGG_id_map(cpds_list2)
    # ModelSEED ids for compounds in model2
    cpds_list2_BiGG_ids = cpds_list2_BiGG_id_map.keys()

    # A dictionary, which maps the non-alphanumeric names of the compounds to 
    # compounds in model2 
    cpds_list2_clean_names_map = create_cleanNames_map(cpds_list2,compart_list = cpds_list2_cpt_ids_to_standard_map.keys())
    # Names of compounds in model2 where all non-alpahnumeric characters are removed 
    cpds_list2_clean_names = cpds_list2_clean_names_map.keys()

    # Create a dictionary which maps the non-alphanumeric formul of the compounds and 
    # their ids in model2 
    cpds_list2_formulas_map = create_cpd_formulas_map(cpds_list2)
    # Formulas of compounds in model2 where all non-alpahnumeric characters are removed 
    cpds_list2_formulas = cpds_list2_formulas_map.keys()

    #-------- Find common compounds ----
    if stdout_msgs:
        print 'Checking for common compounds ...'

    cpd1_to_cpd2_map = dict((cpd1,[]) for cpd1 in cpds_list1) 
    cpd1_to_cpd2_noCopart_map = dict((cpd1,[]) for cpd1 in cpds_list1) 
    cpd1_to_cpd2_match_found_by = dict((cpd1.id,'No match') for cpd1 in cpds_list1) 
    cpd1_to_cpd2_noCompart_match_found_by = dict((cpd1.id,'No match') for cpd1 in cpds_list1) 

    # counters showing how many of common reactions were found by comparing their
    # ModelSEED id, KEGG id, name or model id
    id_counter = 0
    ModelSEED_counter = 0
    KEGG_counter = 0
    BiGG_counter = 0
    name_counter = 0
    formula_counter = 0

    for cpd1 in cpds_list1:

        # A list containing compounds form cpds_list2 matching cpd1 by id, BiGG id, KEGG id, name and formula
        match_byID = []
        match_byModelSEED = []
        match_byBiGG = []
        match_byKEGG = []
        match_byName = []
        match_byFormula = []

        # Search by id, If a match found by id no need to search by other criteria
        if cpd1.id in cpds_list2_ids:
            cpd1_to_cpd2_map[cpd1] += cpds_list2_id_map[cpd1.id]
            id_counter += 1

        else:
            # Compre ModelSEED ids. Consider ModelSEED id only if it is unique for cpd1
            ModelSEED_intersect = list(set(cpd1.ModelSEED_id).intersection(set(cpds_list2_ModelSEED_ids)): 
            if len(ModelSEED_intersect) > 0:
                ModelSEED_counter += 1
                for mid in ModelSEED_intersect:
                    match_byModelSEED += cpds_list2_ModelSEED_id_map[mid] 
                match_byModelSEED = list(set(match_byModelSEED))
  
            # Compre KEGG ids. Consider KEGG id only if it is unique for cpd1
            KEGG_intersect = list(set(cpd1.KEGG_id).intersection(set(cpds_list2_KEGG_ids)): 
            if len(KEGG_intersect) > 0:
                KEGG_counter += 1
                for kid in KEGG_intersect:
                    match_byKEGG += cpds_list2_KEGG_id_map[kid] 
                match_byKEGG = list(set(match_byKEGG))

            # Compre BiGG ids. Consider BiGG id only if it is unique for cpd1
            BiGG_intersect = list(set(cpd1.BiGG_id).intersection(set(cpds_list2_BiGG_ids)): 
            if len(BiGG_intersect) > 0:
                BiGG_counter += 1
                for bid in BiGG_intersect:
                    match_byBiGG += cpds_list2_BiGG_id_map[bid] 
                match_byBiGG = list(set(match_byBiGG))    

            # Compare names
            if cpd1.name != '':
                clean_names = [remove_non_alphanumeric(remove_compartment(cpd1.name,compartments_info = cpds_list1_compart_ids)).lower()] + [remove_non_alphanumeric(remove_compartment(n,compartments_info = cpds_list1_compart_ids)).lower() for n in cpd1.name_aliases]
                clean_names_intersect = list(set(clean_names).intersection(set(cpds_list2_clean_names)): 
                if len(clean_names_intersect) > 0:
                    clean_name_counter += 1
                    match_byName += clean_names_intersect 
                match_byName = list(set(match_byName))    

            # Compare formula
            if cpd1.formula != '' and cpd.formula in cpds_list2_formulas: 
                formula_counter += 1
                match_byFormula = cpds_list2_formulas_map[cpd1.formula] 
             
        #--- Assign an accurate ModelSEED id according to searches above ---
        # If match_byID is not '', you've found a unique match id
        if match_byID != []:
            cpd1_to_cpd2_map[cpd1] = match_byID
            cpd1_to_cpd2_noCompart_map[cpd1] = match_byID
            cpd1_to_cpd2_match_found_by[cpd.id] = 'model id'
            cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'model id'

        else:
            # If only one match was found by ModelSEED_id, you're done
            if len(match_byModelSEED_id) == 1:
                cpd1_to_cpd2_noCompart_map[cpd1] = match_byModelSEED_id
                cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'ModelSEED_id id'
    
            # If only one match was found by KEGG_id, you're done
            elif len(match_byKEGG) == 1:
                cpd1_to_cpd2_noCompart_map[cpd1] = match_byKEGG
                cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'KEGG id'
    
            # If only one match was found by BiGG_id, you're done
            elif len(match_byBiGG) == 1:
                cpd1_to_cpd2_noCompart_map[cpd1] = match_byBiGG
                cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'BiGG id'
    
            # If only one match was found by formula, you're done
            elif len(match_byFormula) == 1:
                cpd1_to_cpd2_noCompart_map[cpd1]  = match_byFormula
                cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'formula'
            else:
                # Find the intersection of all searches
                if len([lst for lst in [match_byModelSEED, match_byKEGG, match_byBiGG, match_byName, match_byFormula] if len(lst) > 0]) > 0:
                    match_intersect = list(set.intersection(*[set(lst) for lst in [match_byModelSEED, match_byKEGG, match_byBiGG, match_byName, match_byFormula] if len(lst) > 0]))
                    if len(match_intersect) == 0:
                        cpd1_to_cpd2_noComaprt_map[cpd1] = match_intersect
                        cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'intersection'
    
                    # Otherwise find their union
                    else:
                        cpd1_to_cpd2_noCompart_map[cpd1] = list(set(match_byKEGG + match_byBiGG + match_byName))
                        cpd1_to_cpd2_noCompart_match_found_by[cpd.id] = 'union'

    # If any match has found by ModelSEED_id, kEGG_id, BiGG_id, Name or Formula, check whether the compartments match too
    if len(cpd1_to_cpd2_map[cpd1]) == 0 and len(cpd1_to_cpd2_noCompart_map[cpd1]) > 0:   
        cpd1_to_cpd2_map[cpd1] = [cpd2 for cpd2 in cpd1_to_cpd2_noCompart_map[cpd] if cpds_list1_cpt_ids_to_standard_map[cpd1.compartment.id] = cpds_list2_cpt_ids_to_standard_map[cpd2.compartment.id]]        

    cpd2_to_cpd1_map = dict((cpd2,[]) for cpd2 in cpds_list2) 
    cpd2_to_cpd1_noCompart_map = dict((cpd2,[]) for cpd2 in cpds_list2) 
    cpd2_to_cpd1_match_found_by = dict((cpd2.id,'No match') for cpd2 in cpds_list1) 
    cpd2_to_cpd1_noCompart_match_found_by = dict((cpd2.id,'No match') for cpd2 in cpds_list1) 
    for cpd1 in [c for c in cpds_list1 if len(cpd1_to_cpd2_map[c]) > 0]:
        for cpd2 in cpd1_to_cpd2_map[cpd1]
            cpd2_to_cpd1_map[cpd2].append(cpd1)
        cpd2_to_cpd1_match_found_by[cpd2.id] = cpd1_to_cpd2_match_found_by[cpd1]
    for cpd1 in [c for c in cpds_list1 if len(cpd1_to_cpd2_noCompart_map[c]) > 0]:
        for cpd2 in cpd1_to_cpd2_noCompart_map[cpd1]
            cpd2_to_cpd1_noCompart_map[cpd2].append(cpd1)
        cpd2_to_cpd1_noCompart_match_found_by[cpd2.id] = cpd1_to_cpd2_noCompart_match_found_by[cpd1]

    # Compounds with no match
    no_match_cpds_list1 = [c for c in cpds_list1 if len(cpd1_to_cpd2_map) == 0] 
    no_match_cpds_list2 = [c for c in cpds_list2 if len(cpd2_to_cpd1_map) == 0] 

    # Compounds with more than one match
    more_than_one_match_cpds_list1 = [(c1.id,[c2.id for c2 in cpd1_to_cpd2_map[c1]]) for c1 in cpds_list1 if len(cpd1_to_cpd2_map) > 1] 
    more_than_one_match_cpds_list2 = [(c2.id,[c1.id for c1 in cpd2_to_cpd1_map[c2]]) for c2 in cpds_list2 if len(cpd2_to_cpd1_map) > 1] 

    if warnings:
        print '\n**WARNING! No match found for {} compounds in cpds_list1 and for {} compounds in cpds_list2'.format(len(no_match_cpds_list1),len(no_match_cpds_list2))
        print '\n**WARNING! More than one match found for the following {} compounds in cpds_list1: {}'.format(len(more_than_one_match_cpds_list1.keys()1),more_than_one_match_cpds_list1)
        print '\n**WARNING! More than one match found for the following {} compounds in cpds_list2: {}'.format(len(more_than_one_match_cpds_list2.keys()1),more_than_one_match_cpds_list2)
 
    if stdout_msgs:
        print '\nSummary of the compare_compounds:'
        print '\tTotal # of compounds in cpds_list1 with a match = {} (out of {})'.format(len([c for c in cpds_list1 if  cpd1_to_cpd2_map[c] != []]), len(cpds_list1))
        print '\tTotal # of compounds in cpds_list2 with a match = {} (out of {})'.format(len([c for c in cpds_list2 if  cpd2_to_cpd1_map[c] != []]), len(cpds_list2))
        print '\t\t# of compounds with a matched model id = ',id_counter
        print '\t\t# of compounds with a match KEGG id = ',KEGG_counter
        print '\t\t# of compounds with a matched BiGG id = ',BiGG_counter
        print '\t\t# of compounds with a matched name = ',name_counter
        print '\t\t# of compounds with a matached formula = ',formula_counter

        print '\n\t\t# of compounds in cpds_list1 with matches found by model id = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'model id'])
        print '\t\t# of compounds with matches found by ModelSEED id = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'ModelSEED id'])
        print '\t\t# of compounds with matches found by KEGG id = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'KEGG id'])
        print '\t\t# of compounds with matches found by BiGG id = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'BiGG id'])
        print '\t\t# of compounds with matches found by formula = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'formula id'])
        print '\t\t# of compounds with matches found by intersecton of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'intersection'])
        print '\t\t# of compounds with matches found by union of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([c for c in cpd1_to_cpd2_match_found_by.keys() if cpd1_to_cpd2_match_found_by[c] == 'union'])

    return cpd1_to_cpd2_map, cpd2_to_cpd1_map, cpd1_to_cpd2_noCopart_map, cpd2_to_cpd1_noCopart_map, cpd1_to_cpd2_match_found_by, cpd2_to_cpd1_match_found_by, cpd1_to_cpd2_noCompart_match_found_by, cpd2_to_cpd1_noCompart_match_found_by


def compare_reactions(rxns_list1, rxns_list2, standard_to_cpds_list1_cpt_ids_map, standard_to_cpds_list2_cpt_ids_map, rxns_list1_id = None, rxns_list2_id = None, cpd1_to_cpd2_map = None, warnings = True, stdout_msgs = True):
    """
    Find common reactions in two given list of reactions

    INPUTS:
    -------
        rxns_list1: A list of reaction objects 
        rxns_list2: A list of reaction objects
     rxns_list1_id: An id for rxns_list1. This is requried if not all reactions in 
                   rxns_list1 belong to the same model or if not model is assigned 
                   to all reactions of list 1
     rxns_list2_id: An id for rxns_list2. This is requried if not all reactions in 
                   rxns_list1 belong to the same model or if not model is assigned 
                   to all reactions of list 2
    standard_to_rxns_list1_cpt_ids_map:
    standard_to_rxns_list2_cpt_ids_map: 
                            A dictionary where keys are standard compartment ids as follows:
                            c: Cytosol (cytoplasm),   e: Extracellular,   g: Golgi,     m: Mitochondria
                            n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,    x: Peroxisome
                            and values are corresponding compartment ids in rxns_list1 and rxns_list2
    cpd1_to_cpd2_map: The output from compare compounds, which is a dictionary mapping compounds participating 
                      in rxns_list1 to those in rxns_list2. If this input is provided it can be used to 
                      compare reactions by matching their stoichiometries

    OUTPUTS:
    --------
    rxn1_to_rxn2_map: A dictionary where the keys are the reactions objects from rxns_list1 and values 
                      are a list of reaction objects corresponding to this reaction from rxns_list2 
    rxn2_to_rxn1_map: A dictionary where the keys are the reactions objects from rxns_list2 and values 
                      are a list of reaction objects corresponding to this reaction from rxns_list1 
    rxn1_to_rxn2_noCopart_map:
    rxn2_to_rxn1_noCopart_map:
                      Same as rxn1_to_rxn2_map and rxn2_to_rxn1_map except that the mataches here are found
                      without considering the compartments
    rxn1_to_rxn2_match_found_by:
    rxn2_to_rxn1_match_found_by:
                      A dictionary with keys being ids of reactions in rxns_list1 or rxns_list2 and values being a 
                      string showing how the matches in rxn1_to_rxn2_map and rxn2_to_rxn1_map were found 
    rxn1_to_rxn2_noCompart_match_found_by:
    rxn2_to_rxn1_noCompart_match_found_by:
                      Same rxn1_to_rxn2_match_found_by and rxn2_to_rxn1_match_found_by for rxn1_to_rxn2_noCopart_map 
                      and rxn2_to_rxn1_noCopart_map
    """
    if not isinstance(strdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be True or False')

    if rxns_list1_id == None:
        model1_id = list(set([r.model.id for r in rxns_list1 if 'model' in dir(r)]))
        if len(model1_id) == 1:
            rxns_list1_id = model1_id[0]
        else:
            raise userError('rxns_list1_id is needed because not all reactions in rxns_list1 have the same model id (' + str(model1_id) + ')')

    if rxns_list2_id == None:
        model2_id = list(set([r.model.id for r in rxns_list2 if 'model' in dir(r)]))
        if len(model2_id) == 1:
            rxns_list2_id = model2_id[0]
        else:
            raise userError('rxns_list2_id is needed because not all reactions in rxns_list2 have the same model id (' + str(model2_id) + ')')

    # Dictionaries where keys are compartment ids in cpds_list1 and cpds_list2 and values are corresponding standard 
    # comaprtment ids
    rxns_list1_cpt_ids_to_standard_map = {}
    rxns_list2_cpt_ids_to_standard_map = {}
    for cpt_id in standard_to_rxns_list1_cpt_ids_map.values():
        rxns_list1_cpt_ids_to_standard_map[cpt_id] = [sid for sid in standard_to_rxns_list1_cpt_ids_map.keys() if standard_to_rxns_list1_cpt_ids_map[sid] == cpt_id][0]
    for cpt_id in standard_to_rxns_list2_cpt_ids_map.values():
        rxns_list2_cpt_ids_to_standard_map[cpt_id] = [sid for sid in standard_to_rxns_list2_cpt_ids_map.keys() if standard_to_rxns_list2_cpt_ids_map[sid] == cpt_id][0]

    # A dictionary mapping the ids of componds in rxns_list2 to compound objects
    rxns_list2_id_map = dict([(c.id,c) for c in rxns_list2])
    # ModelSEED ids for compounds in model2
    rxns_list2_ids = rxns_list2_id_map.keys()

    # A dictionary mapping the ModelSEED ids of reactions in rxns_list2 to reaction objects in that list 
    rxns_list2_ModelSEED_id_map = create_ModelSEED_id_map(rxns_list2)
    # ModelSEED ids for reactions in rxns_list2
    rxns_list2_ModelSEED_ids = rxns_list2_ModelSEED_id_map.keys()

    # A dictionary mapping the KEGG ids of reactions in rxns_list2 to reaction objects in that list 
    rxns_list2_KEGG_id_map = create_KEGG_id_map(rxns_list2)
    # KEGG ids for reactions in rxns_list2
    rxns_list2_KEGG_ids = rxns_list2_KEGG_id_map.keys()

    # A dictionary mapping the BiGG ids of reactions in rxns_list2 to reaction objects in that list 
    rxns_list2_BiGG_id_map = create_BiGG_id_map(rxns_list2)
    # BiGG ids of reactions in rxns_list2
    rxns_list2_BiGG_ids = rxns_list2_BiGG_id_map.keys()

    # A dictionary mapping the clean names of reactions in rxns_list2 to reaction objects in that list 
    rxns_list2_clean_names_map = create_cleanNames_map(rxns_list2,compart_list = rxns_list2_cpt_ids_to_standard_map.keys())
    # Clean names of reactions in rxns_list2m
    rxns_list2_clean_names = rxns_list2_clean_names_map.keys()

    # A dictionary mapping the stoichiometry of reactions in rxns_list2 (as a tuple of tuples) to the corresponding reaction object
    rxns_list2_stoic_map = dict([(tuple(sorted(r.stoichiometry.items(), key = lambda x:x[0].id)),r) for r in rxns_list2]) 
    # Stoichiometries of reactions in rxns_list2
    rxns_list2_stoics = rxns_list2_stoic_map.keys()

    rxn1_to_rxn2_map = dict((rxn1,[]) for rxn1 in rxns_list1)
    rxn1_to_rxn2_noCopart_map = dict((rxn1,[]) for rxn1 in rxns_list1)
    rxn1_to_rxn2_match_found_by = dict((rxn1.id,'No match') for rxn1 in rxns_list1)
    rxn1_to_rxn2_noCompart_match_found_by = dict((rxn1.id,'No match') for rxn1 in rxns_list1)

    id_counter = 0
    ModelSEED_counter = 0
    KEGG_counter = 0
    BiGG_counter = 0
    name_counter = 0
    stoic_counter = 0

    for rxn1 in rxns_list1:

        # A list containing reactions form rxns_list2 matching rxn1 by id, BiGG id, KEGG id, name and stoic
        match_byID = []
        match_byModelSEED = []
        match_byBiGG = []
        match_byKEGG = []
        match_byName = []
        match_byStoic = []

        # Search by id, If a match found by id no need to search by other criteria
        if rxn1.id in rxns_list2_ids:
            rxn1_to_rxn2_map[rxn1] += rxns_list2_id_map[rxn1.id]
            id_counter += 1

        else:
            # Compre ModelSEED ids. Consider ModelSEED id only if it is unique for rxn1
            ModelSEED_intersect = list(set(rxn1.ModelSEED_id).intersection(set(rxns_list2_ModelSEED_ids)): 
            if len(ModelSEED_intersect) > 0:
                ModelSEED_counter += 1
                for mid in ModelSEED_intersect:
                    match_byModelSEED += rxns_list2_ModelSEED_id_map[mid] 
                match_byModelSEED = list(set(match_byModelSEED))
  
            # Compre KEGG ids. Consider KEGG id only if it is unique for rxn1
            KEGG_intersect = list(set(rxn1.KEGG_id).intersection(set(rxns_list2_KEGG_ids)): 
            if len(KEGG_intersect) > 0:
                KEGG_counter += 1
                for kid in KEGG_intersect:
                    match_byKEGG += rxns_list2_KEGG_id_map[kid] 
                match_byKEGG = list(set(match_byKEGG))

            # Compre BiGG ids. Consider BiGG id only if it is unique for rxn1
            BiGG_intersect = list(set(rxn1.BiGG_id).intersection(set(rxns_list2_BiGG_ids)): 
            if len(BiGG_intersect) > 0:
                BiGG_counter += 1
                for bid in BiGG_intersect:
                    match_byBiGG += rxns_list2_BiGG_id_map[bid] 
                match_byBiGG = list(set(match_byBiGG))    

            # Compare names
            if rxn1.name != '':
                clean_names = [remove_non_alphanumeric(remove_compartment(rxn1.name,compartments_info = rxns_list1_compart_ids)).lower()] + [remove_non_alphanumeric(remove_compartment(n,compartments_info = rxns_list1_compart_ids)).lower() for n in rxn1.name_aliases]
                clean_names_intersect = list(set(clean_names).intersection(set(rxns_list2_clean_names)): 
                if len(clean_names_intersect) > 0:
                    clean_name_counter += 1
                    match_byName += clean_names_intersect 
                match_byName = list(set(match_byName))    

            # Compare stoic only if all compoudns participating in rxns have at least one match in reactions
            # participating in rxns_list2
            if cpd1_to_cpd2_map != None and len([c for c in rxn1.stoichiometry.keys() if len(cpd1_to_cpd2_map[c]) > 0]) == len(rxn1.stoichiometry.keys()):
                # Stoichiometry of rxn1 as tuple with reactions replaced with their matches in rxns_list2
                rxn1_stoic = tuple(sorted([(cpd1_to_cpd2_map[c][0],rxn1.stoichiometry[c]) for c in rxn1.stoichiometry.keys()], key = lambda x:x[0].id))
                if rxn1_stoic in rxns_list2_stoics: 
                    stoic_counter += 1
                    match_byStoic = rxns_list2_stoics_map[rxn1_stoic] 
             
        #--- Assign an accurate ModelSEED id according to searches above ---
        # If match_byID is not '', you've found a unique match id
        if match_byID != []:
            rxn1_to_rxn2_map[rxn1] = match_byID
            rxn1_to_rxn2_noCompart_map[rxn1] = match_byID
            rxn1_to_rxn2_match_found_by[rxn.id] = 'model id'
            rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'model id'

        # If only one match was found by stoic, you're done
        elif len(match_byStoic) == 1:
            rxn1_to_rxn2_map[rxn1] = match_byStoic
            rxn1_to_rxn2_noCompart_map[rxn1] = match_byStoic
            rxn1_to_rxn2_match_found_by[rxn.id] = 'stoic'
            rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'stoic'

        else:
            # If only one match was found by ModelSEED_id, you're done
            if len(match_byModelSEED_id) == 1:
                rxn1_to_rxn2_noCompart_map[rxn1] = match_byModelSEED_id
                rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'ModelSEED_id id'
    
            # If only one match was found by KEGG_id, you're done
            elif len(match_byKEGG) == 1:
                rxn1_to_rxn2_noCompart_map[rxn1] = match_byKEGG
                rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'KEGG id'
    
            # If only one match was found by BiGG_id, you're done
            elif len(match_byBiGG) == 1:
                rxn1_to_rxn2_noCompart_map[rxn1] = match_byBiGG
                rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'BiGG id'
    
            else:
                # Find the intersection of all searches
                if len([lst for lst in [match_byModelSEED, match_byKEGG, match_byBiGG, match_byName] if len(lst) > 0]) > 0:
                    match_intersect = list(set.intersection(*[set(lst) for lst in [match_byKEGG, match_byBiGG, match_byName] if len(lst) > 0]))
                    if len(match_intersect) == 0:
                        rxn1_to_rxn2_noComaprt_map[rxn1] = match_intersect
                        rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'intersection'
    
                    # Otherwise find their union
                    else:
                        rxn1_to_rxn2_noCompart_map[rxn1] = list(set(match_byKEGG + match_byBiGG + match_byName))
                        rxn1_to_rxn2_noCompart_match_found_by[rxn.id] = 'union'

    # If no match was found by model id or stoic by there are matches found by ModelSEED_id, kEGG_id, BiGG_id or 
    # Name check whether the compartments match too
    if len(rxn1_to_rxn2_map[rxn1]) == 0 and len(rxn1_to_rxn2_noCompart_map[rxn1]) > 0:   
        rxn1_to_rxn2_map[rxn1] = [rxn2 for rxn2 in rxn1_to_rxn2_noCompart_map[rxn] if set([rxns_list1_cpt_ids_to_standard_map[cpt.id] for cpt in rxn1.compartment]) = set([rxns_list2_cpt_ids_to_standard_map[cpt.id] for cpt in rxn2.compartment])]        

    rxn2_to_rxn1_map = dict((rxn2,[]) for rxn2 in rxns_list2) 
    rxn2_to_rxn1_noCompart_map = dict((rxn2,[]) for rxn2 in rxns_list2) 
    rxn2_to_rxn1_match_found_by = dict((rxn2.id,'No match') for rxn2 in rxns_list1) 
    rxn2_to_rxn1_noCompart_match_found_by = dict((rxn2.id,'No match') for rxn2 in rxns_list1) 
    for rxn1 in [r for r in rxns_list1 if len(rxn1_to_rxn2_map[c]) > 0]:
        for rxn2 in rxn1_to_rxn2_map[rxn1]
            rxn2_to_rxn1_map[rxn2].append(rxn1)
        rxn2_to_rxn1_match_found_by[rxn2.id] = rxn1_to_rxn2_match_found_by[rxn1]
    for rxn1 in [r for r in rxns_list1 if len(rxn1_to_rxn2_noCompart_map[c]) > 0]:
        for rxn2 in rxn1_to_rxn2_noCompart_map[rxn1]
            rxn2_to_rxn1_noCompart_map[rxn2].append(rxn1)
        rxn2_to_rxn1_noCompart_match_found_by[rxn2.id] = rxn1_to_rxn2_noCompart_match_found_by[rxn1]

    # Compounds with no match
    no_match_rxns_list1 = [r for r in rxns_list1 if len(rxn1_to_rxn2_map) == 0] 
    no_match_rxns_list2 = [r for r in rxns_list2 if len(rxn2_to_rxn1_map) == 0] 

    # Compounds with more than one match
    more_than_one_match_rxns_list1 = [(c1.id,[c2.id for c2 in rxn1_to_rxn2_map[c1]]) for c1 in rxns_list1 if len(rxn1_to_rxn2_map) > 1] 
    more_than_one_match_rxns_list2 = [(c2.id,[c1.id for c1 in rxn2_to_rxn1_map[c2]]) for c2 in rxns_list2 if len(rxn2_to_rxn1_map) > 1] 

    if warnings:
        print '\n**WARNING! No match found for {} reactions in rxns_list1 and for {} reactions in rxns_list2'.format(len(no_match_rxns_list1),len(no_match_rxns_list2))
        print '\n**WARNING! More than one match found for the following {} reactions in rxns_list1: {}'.format(len(more_than_one_match_rxns_list1.keys()1),more_than_one_match_rxns_list1)
        print '\n**WARNING! More than one match found for the following {} reactions in rxns_list2: {}'.format(len(more_than_one_match_rxns_list2.keys()1),more_than_one_match_rxns_list2)
 
    if stdout_msgs:
        print '\nSummary of the compare_reactions:'
        print '\tTotal # of reactions in rxns_list1 with a match = {} (out of {})'.format(len([r for r in rxns_list1 if  rxn1_to_rxn2_map[c] != []]), len(rxns_list1))
        print '\tTotal # of reactions in rxns_list2 with a match = {} (out of {})'.format(len([r for r in rxns_list2 if  rxn2_to_rxn1_map[c] != []]), len(rxns_list2))
        print '\t\t# of reactions with a matched model id = ',id_counter
        print '\t\t# of reactions with a match KEGG id = ',KEGG_counter
        print '\t\t# of reactions with a matched BiGG id = ',BiGG_counter
        print '\t\t# of reactions with a matched name = ',name_counter
        print '\t\t# of reactions with a matached stoic = ',stoic_counter

        print '\n\t\t# of reactions in rxns_list1 with matches found by model id = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'model id'])
        print '\t\t# of reactions with matches found by ModelSEED id = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'ModelSEED id'])
        print '\t\t# of reactions with matches found by KEGG id = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'KEGG id'])
        print '\t\t# of reactions with matches found by BiGG id = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'BiGG id'])
        print '\t\t# of reactions with matches found by stoic = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'stoic id'])
        print '\t\t# of reactions with matches found by intersecton of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'intersection'])
        print '\t\t# of reactions with matches found by union of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([r for r in rxn1_to_rxn2_match_found_by.keys() if rxn1_to_rxn2_match_found_by[c] == 'union'])

    return rxn1_to_rxn2_map, rxn2_to_rxn1_map, rxn1_to_rxn2_noCopart_map, rxn2_to_rxn1_noCopart_map, rxn1_to_rxn2_match_found_by, rxn2_to_rxn1_match_found_by, rxn1_to_rxn2_noCompart_match_found_by, rxn2_to_rxn1_noCompart_match_found_by


def compare_models(model1, model2, standard_to_model1_cpt_ids_map, standard_to_model2_cpt_ids_map, obtain_ModelSEED_ids = True, warnings = True, stdout_msgs = True):
    """
    This compares two metabolic models. It comapres both reactions and compounds. 
    by comparing their id, KEGG_id and ModelSEED_id, names or formula (for compunds only)

    INPUTS:
    ------
                  model1: Model 1
                  model2: Model 2
    obtain_ModelSEED_ids: If True ModelSEED ids are obtained first for model1 and model2
    standard_to_model1_cpt_ids_map:
    standard_to_model2_cpt_ids_map: 
                           A dictionary where keys are standard compartment ids as follows:
                           c: Cytosol (cytoplasm),   e: Extracellular,   g: Golgi,     m: Mitochondria
                           n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,    x: Peroxisome
                           and values are corresponding compartment ids in model1 and model2


    OUTPUTS:
    --------
    """
    if not isinstance(strdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
       
    # Try to get the ModelSEED ids 
    if stdout_msgs:
        print '\n---------- Compre two metabolic models -----------\n'

    # Get ModelSEED ids for reactions and metabolites
    if obtain_ModelSEED_ids:
        get_ModelSEED_ids(model1, warnings = warnings, stdout_msgs = stdout_msgs)
        get_ModelSEED_ids(model2, warnings = warnings, stdout_msgs = stdout_msgs)

    # Compare compounds
    cpd1_to_cpd2_map, cpd2_to_cpd1_map = compare_compounds(cpds_list1i = model1.compounds, cpds_list2 = model2.compounds, standard_to_model1_cpt_ids_map, standard_to_model2_cpt_ids_map, cpds_list1_id = model1.id, cpds_list2_id = model2.id, warnings = warnings, stdout_msgs = stdout_msgs)[0:2]

    # Compare reactions
    rxn1_to_rxn2_map, rxn2_to_rxn1_map = compare_reactions(rxns_list1 = model1.reactions, rxns_list2 = model2.reactions, standard_to_cpds_list1_cpt_ids_map = standard_to_model1_cpt_ids_map, standard_to_cpds_list2_cpt_ids_map = standard_to_model2_cpt_ids_map, rxns_list1_id = model1.id, rxns_list2_id = model2.id, cpd1_to_cpd2_map = cpd1_to_cpd2_map, warnings = warnings, stdout_msgs = stdout_msgs)[0:2]

    if stdout_msgs:
        print '\n---------- Done with comparing models -----------\n'

    return cpd1_to_cpd2_map, cpd2_to_cpd1_map, rxn1_to_rxn2_map, rxn2_to_rxn1_map 

