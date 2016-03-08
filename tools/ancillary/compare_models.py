from __future__ import division
import sys, time
sys.path.append('../../')
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
from get_ModelSeed_ids import get_ModelSeed_ids, remove_compartment
from remove_non_alphanumeric import remove_non_alphanumeric

#--------------------------------------
# Ali R. Zomorrodi - Segre lab @ BU
# Last updated: 05-13-2015
#--------------------------------------

#--------------------------------------
#--------- Ancillary functions --------
#--------------------------------------
def create_ModelSeedID_map(obj_list):
    """ 
    Creates a dictionary mapping the ModelSeed ids for componds or reactions to objects of
    type compound or reaction in a list. This is applied to objects with a ModelSeed id 

    INPUTS:
    -------
    obj_list: A list of objects of type compounds or reactions

    OUTPUTS:
    -------
    ModelSeed_id_map: A dictionary where keys are ModelSeed ids and values are 
                      objects of type compounds or reactions   
    """ 
    ModelSeed_id_map = {}
    for obj in [o for o in obj_list if o.ModelSeed_id != None]:
        if isinstance(obj.ModelSeed_id,str):
            if obj.ModelSeed_id in ModelSeed_id_map.keys():
                ModelSeed_id_map[obj.ModelSeed_id] += [obj]
            else:
                ModelSeed_id_map[obj.ModelSeed_id] = [obj]
        elif isinstance(obj.ModelSeed_id,list):
            for sid in obj.ModelSeed_id:
                if sid in ModelSeed_id_map.keys():
                    ModelSeed_id_map[sid] += [obj]
                else:
                    ModelSeed_id_map[sid] = [obj]
    return  ModelSeed_id_map

def create_KeggID_map(obj_list):
    """ 
    Creates a dictionary mapping the Kegg ids for componds or reactions to objects of
    type compound or reaction in a list. This is applied to objects with a Kegg id 

    INPUTS:
    -------
    obj_list: A list of objects of type compound or reaction

    OUTPUTS:
    -------
    Kegg_id_map: A dictionary where keys are Kegg ids and values are 
                 objects of type compound or reaction   
    """ 
    Kegg_id_map = {}
    for obj in [o for o in obj_list if o.Kegg_id != None]:
        if isinstance(obj.Kegg_id,str):
            if obj.Kegg_id in Kegg_id_map.keys():
                Kegg_id_map[obj.Kegg_id] += [obj]
            else:
                Kegg_id_map[obj.Kegg_id] = [obj]
        elif isinstance(obj.Kegg_id,list):
            for sid in obj.Kegg_id:
                if sid in Kegg_id_map.keys():
                    Kegg_id_map[sid] += [obj]
                else:
                    Kegg_id_map[sid] = [obj]
    return Kegg_id_map

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
        for name in obj.name:
            clean_name = remove_non_alphanumeric(remove_compartment(name,compartments_info = compart_list)).lower()
            if clean_name in clean_names_map.keys():
                clean_names_map[clean_name] += [obj]
            else:
                clean_names_map[clean_name] = [obj]
    # Names of objects in model2 where all non-alpahnumeric characters are removed 
    return clean_names_map

def create_cmp_cleanFormulas_map(cmp_list):
    """ 
    Create a dictionary mapping the non-alphanumeric formula of the compounds to  
    compound objects in a given list of compounds. This is applied only to objects
    with a given formula

    INPUTS:
    -------
    cmp_list: A list of objects of type compound

    OUTPUTS:
    -------
    clean_formulas_map: A dictionary where keys are non-alphanumeric names and values are 
    """ 
    clean_formulas_map = {}
    for cmp in [c for c in cmp_list if c.formula != None]:
        clean_formula = remove_non_alphanumeric(cmp.formula).lower()
        if clean_formula in clean_formulas_map.keys():
            clean_formulas_map[clean_formula] += [cmp]
        else:
            clean_formulas_map[clean_formula] = [cmp]
    return clean_formulas_map

#--------------------------------------
#--------- Main functions --------
#--------------------------------------
def compare_compounds(cmp_list1,cmp_list2,cmp_list1_id = None, cmp_list2_id = None, compart_list1 = None, compart_list2 = None, warnings = True, stdout_msgs = True):
    """
    Find common compounds in two given list of compounds

    INPUTS:
    -------
        cmp_list1: An instance of object compound or a list of such objects
        cmp_list2: An instance of object compound or a list of such objects
     cmp_list1_id: An id for cmp_list1. This is requried if not all compounds in 
                   cmp_list1 belong to the same model or if no model is assigned 
                   to all reacitons of list 1
     cmp_list2_id: An id for cmp_list2. This is requried if not all compounds in 
                   cmp_list1 belong to the same model or if no model is assigned 
                   to all reacitons of list 2
    compart_list1: A list of string containing the ids of compartments for cmp_list1 
    compart_list2: A list of string containing the ids of compartments for cmp_list2 
         warnings: Can be True or False showing whether warnings should be written  
                   in the output
      stdout_msgs: Can be True or False showing whether details of comparison should be
                   written in the output

    OUTPUTS:
    --------
    common_compounds: A list of tuples where the first element is an object of type compound
                      referring to a compound in cmp_list 1 and the second element is an object 
                      of list of objects of type compound containing the corresponding compounds
                      in cmp_list2 
               match: This function also assigns to each compound and reaction object a new field 
                      termed match. If match is None it means that reaction or compound has no 
                      match in other models. Otherwise it is a list of objects of type compound
                      (for compounds) that match to this compound in other models
    """
    if not isinstance(strdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be True or False')

    if cmp_list1_id == None:
        model1_id = list(set([c.model.id for c in cmp_list1 if 'model' in dir(c)]))
        if len(model1_id) == 1:
            cmp_list1_id = model1_id[0]
        else:
            raise userError('cmp_list1_id is needed because not all compounds in cmp_list1 have the same model id (' + str(model1_id) + ')')

    if cmp_list2_id == None:
        model2_id = list(set([c.model.id for c in cmp_list2 if 'model' in dir(c)]))
        if len(model2_id) == 1:
            cmp_list2_id = model2_id[0]
        else:
            raise userError('cmp_list2_id is needed because not all compounds in cmp_list2 have the same model id (' + str(model2_id) + ')')

    # List of compartments
    if compart_list1 == None:
       compart_list1 = list(set([c.compartment.id for c in cmp_list1])) 
    if compart_list2 == None:
       compart_list2 = list(set([c.compartment.id for c in cmp_list2])) 

    # A dictionary mapping the ModelSeed ids of componds to compounds in cmp_list2
    cmp_list2_ModelSeed_id_map = create_ModelSeedID_map(cmp_list2)
    # ModelSeed ids for compounds in model2
    cmp_list2_ModelSeed_ids = cmp_list2_ModelSeed_id_map.keys()

    # A dictionary mapping the ModelSeed ids of componds to compounds in model2
    cmp_list2_Kegg_id_map = create_KeggID_map(cmp_list2)
    # ModelSeed ids for compounds in model2
    cmp_list2_Kegg_ids = cmp_list2_Kegg_id_map.keys()

    # A dictionary which maps the non-alphanumeric names of the compounds to 
    # compounds in model2 
    cmp_list2_clean_names_map = create_cleanNames_map(cmp_list2,compart_list = compart_list2)
    # Names of compounds in model2 where all non-alpahnumeric characters are removed 
    cmp_list2_clean_names = cmp_list2_clean_names_map.keys()

    # Create a dictionary which maps the non-alphanumeric formul of the compounds and 
    # their ids in model2 
    cmp_list2_clean_formulas_map = create_cmp_cleanFormulas_map(cmp_list2)
    # Formulas of compounds in model2 where all non-alpahnumeric characters are removed 
    cmp_list2_clean_formulas = cmp_list2_clean_formulas_map.keys()

    # Add a new field termed hasMatch to all reactions, which indicates whether
    # it has a match in another model
    for cmp in cmp_list1:
        cmp.match = None

    for cmp in cmp_list2:
        cmp.match = None

    #-------- Find common compounds ----
    print 'Checking for common compounds ...'
    common_compounds = [] 

    # counters showing how many of common reactions were found by comparing their
    # ModelSeed id, Kegg id, name or model id
    ModelSeed_counter = 0
    Kegg_counter = 0
    name_counter = 0
    id_counter = 0
    formula_counter = 0

    for cmp1 in cmp_list1:
        # Remove non-alphanumeric characters and compartment ids from the id and replace 
        # all upper case letters with lowercase for name, id and formula
        clean_id = remove_non_alphanumeric(remove_compartment(cmp1.id,compartments_info = compart_list1)).lower()
        if cmp1.name != None:
            if isinstance(cmp1.name,str):
                clean_name = remove_non_alphanumeric(remove_compartment(cmp1.name,compartments_info = compart_list1)).lower()
            else: # if it is a list
                clean_name = [remove_non_alphanumeric(remove_compartment(n,compartments_info = compart_list1)).lower() for n in cmp1.name]
           
        if cmp1.formula != None:
            clean_formula = remove_non_alphanumeric(cmp1.formula).lower()

        # Compounds in model 2 that match to this compound
        match_m2 = []

        # Compre ModelSeed ids. Consider ModelSeed id only if it is unique for cmp1
        if isinstance(cmp1.ModelSeed_id,str) and cmp1.ModelSeed_id in cmp_list2_ModelSeed_ids: 
            ModelSeed_counter += 1
            match_m2 += cmp_list2_ModelSeed_id_map[cmp1.ModelSeed_id] 

        # Compare Kegg ids
        elif isinstance(cmp1.Kegg_id,str) and cmp1.Kegg_id in cmp_list2_Kegg_ids: 
            Kegg_counter += 1
            match_m2 += cmp_list2_Kegg_id_map[cmp1.Kegg_id] 

        # Compare ids in the model
        elif len([cmp2 for cmp2 in cmp_list2 if cmp2.id.lower() == cmp1.id.lower()]) > 0:
            id_counter += 1
            match_m2 += [cmp2 for cmp2 in cmp_list2 if cmp2.id.lower() == cmp1.id.lower()]

        # Compare names
        elif cmp1.name != None:
            if isinstance(cmp1.name,str) and clean_name in cmp_list2_clean_names:  
                name_counter += 1
                match_m2 += cmp_list2_clean_names_map[clean_name]

            # cmp1 name is list and cmp2 name is string
            elif isinstance(cmp1.name,list) and len([n for n in cmp1.name if n in cmp_list2_clean_names]) > 0: 
                name_counter += 1
                match_m2 = []
                for n in [n for n in cmp1.name if n in cmp_list2_clean_names]:
                    match_m2 += cmp_list2_clean_names_map[n] 

        # Compare formula
        elif cmp1.formula != None and clean_formula in cmp_list2_clean_formulas: 
            formula_counter += 1
            match_m2 = cmp_list2_clean_formulas_map[clean_formula]

        if len(match_m2) >= 1:
            # Check the compartments
            compartment_match = [c for c in match_m2 if c.compartment.id.lower() == cmp1.compartment.id.lower() or c.compartment.name.lower() == cmp1.compartment.name.lower()]
            if len(compartment_match) == 1:
                match_m2 = compartment_match[0]
                print 'hello 1'
            elif len(compartment_match) > 1:
                match_m2 = compartment_match
                print 'hello 2'
            else:
                match_m2 = []

            if len(match_m2) >= 1:
                cmp1.match = match_m2
                match_m2.match = cmp1
                common_compounds.append((cmp1,match_m2))           

    # Compounds in cmp_list1 with more than one match in cmp_list2
    more_than_one_match_cmps = [cmp_pair for cmp_pair in common_compounds if isinstance(cmp_pair[1],list)]
    if len(more_than_one_match_cmps) > 0 and warnings:
         print 'WARNING! ',len(more_than_one_match_cmps),' compounds in ',cmp_list1_id,' have more than one match in ',cmp_list2_id
         for cmp_pair in  more_than_one_match_cmps:
             print '\t', (((cmp_list_id,cmp_pair[0].id),(cmp_list2_id,cmp_pair[1].id)))
    
    if stdout_msgs:
        print '\nSummary of the compare_compounds:'
        print '\tTotal # of compounds in %s = %i '%(cmp_list1_id,len(cmp_list1)) 
        print '\tTotal # of compounds in %s = %i '%(cmp_list2_id,len(cmp_list2)) 
        print '\tTotal # of common compounds = ',len(common_compounds)
        print '\t\t# of compounds matched by ModelSeed id = ',ModelSeed_counter
        print '\t\t# of compounds matched by Kegg id = ',Kegg_counter
        print '\t\t# of compounds matched by name = ',name_counter
        print '\t\t# of compounds matached by model id = ',id_counter
        print '\t\t# of compounds matached by formula = ',formula_counter
        print '\tTotal # of compounds unique to %s = %i '%(cmp_list1_id,len([c for c in cmp_list1 if c.match != None])) 
        print '\tTotal # of compounds unique to %s = %i\n'%(cmp_list2_id,len([c for c in cmp_list2 if c.match != None])) 

    return common_compounds

def compare_reactions(rxn_list1,rxn_list2,rxn_list1_id = None, rxn_list2_id = None, compart_list1 = None, compart_list2 = None, warnings = True, stdout_msgs = True):
    """
    Find common reactions in two given list of reactions

    INPUTS:
    -------
        rxn_list1: An instance of object reaction or a list of such objects
        rxn_list2: An instance of object reaction or a list of such objects
     rxn_list1_id: An id for rxn_list1. This is requried if not all reactions in 
                   rxn_list1 belong to the same model or if not model is assigned 
                   to all reacitons of list 1
     rxn_list2_id: An id for rxn_list2. This is requried if not all reactions in 
                   rxn_list1 belong to the same model or if not model is assigned 
                   to all reacitons of list 2
    compart_list1: A list of string containing the ids of compartments for rxn_list 
    compart_list2: A list of string containing the ids of compartments for rxn_list 
         warnings: Can be True or False showing whether warnings should be written  
                   in the output
      stdout_msgs: Can be True or False showing whether details of comparison should be
                   written in the output


    OUTPUTS:
    --------
    common_reactions: A list of tuples where the first element is an object of type reaction
                      referring to a compound in rxn_list 1 and the second element is an object 
                      of list of objects of type reaction containing the corresponding reactions
                      in rxn_list2 
    """
    if not isinstance(strdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be True or False')

    if rxn_list1_id == None:
        model1_id = list(set([r.model.id for r in rxn_list1 if 'model' in dir(r)]))
        if len(model1_id) == 1:
            rxn_list1_id = model1_id[0]
        else:
            raise userError('rxn_list1_id is needed because not all reactions in rxn_list1 have the same model id (' + str(model1_id) + ')')

    if rxn_list2_id == None:
        model2_id = list(set([r.model.id for r in rxn_list2 if 'model' in dir(r)]))
        if len(model2_id) == 1:
            rxn_list2_id = model2_id[0]
        else:
            raise userError('rxn_list2_id is needed because not all reactions in rxn_list2 have the same model id (' + str(model2_id) + ')')

    # List oc compartments
    if compart_list1 == None:
       compart_list1 = list(set([compart.id for r in rxn_list1 for compart in r.compartments]))
    if compart_list2 == None:
       compart_list2 = list(set([compart.id for r in rxn_list2 for compart in r.compartments]))


    # A dictionary mapping the ModelSeed ids of reactions ito reactions in model2
    rxn_list2_ModelSeed_id_map = create_ModelSeedID_map(rxn_list2)
    # ModelSeed ids for reactions in model2
    rxn_list2_ModelSeed_ids = rxn_list2_ModelSeed_id_map.keys()

    # A dictionary mapping the ModelSeed ids of reactions ito reactions in model2
    rxn_list2_Kegg_id_map = create_KeggID_map(rxn_list2)
    # ModelSeed ids for reactions in model2
    rxn_list2_Kegg_ids = rxn_list2_Kegg_id_map.keys()

    # A dictionary which maps the non-alphanumeric names of the reactions to 
    # ids in model2 
    rxn_list2_clean_names_map = create_cleanNames_map(rxn_list2,compart_list = compart_list2)
    # Names of reactions in model2 where all non-alpahnumeric characters are removed 
    rxn_list2_clean_names = rxn_list2_clean_names_map.keys()

    common_reactions = [] 
    common_reaction = []

    # Add a new field termed hasMatch to all reactions, which indicates whether
    # it has a match in another model
    for rxn in rxn_list1:
        rxn.match = None

    for rxn in rxn_list2:
        rxn.match = None

    # counters showing how many of common reactions were found by comparing their
    # ModelSeed id, Kegg id, name or model id
    ModelSeed_counter = 0
    Kegg_counter = 0
    name_counter = 0
    id_counter = 0

    #--- Start the search ----
    for rxn1 in rxn_list1:
        # Remove non-alphanumeric characters and compartment ids from the id and replace 
        # all upper case letters with lowercase for name, id and formula
        clean_id = remove_non_alphanumeric(remove_compartment(rxn1.id,compartments_info = compart_list1)).lower()
        if rxn1.name != None:
            if isinstance(rxn1.name,str):
                clean_name = remove_non_alphanumeric(remove_compartment(rxn1.name,compartments_info = compart_list1)).lower()
            else: # if it is a list
                clean_name = [remove_non_alphanumeric(remove_compartment(n,compartments_info = compart_list1)).lower() for n in rxn1.name]

        # Reactions in model 2 that match to this reaction
        match_m2 = []

        # Comprae ModelSeed ids. Consider the ModelSeed id if it is unique for rxn1
        if isinstance(rxn1.ModelSeed_id,str) and rxn1.ModelSeed_id in rxn_list2_ModelSeed_ids: 
            ModelSeed_counter += 1
            match_m2 += rxn_list2_ModelSeed_id_map[rxn1.ModelSeed_id] 

        # Compare kegg ids
        elif rxn1.Kegg_id != None and rxn1.Kegg_id in rxn_list2_Kegg_ids: 
            Kegg_counter += 1
            match_m2 += rxn_list2_Kegg_id_map[rxn1.Kegg_id] 

        # Compare ids
        elif len([rxn2 for rxn2 in rxn_list2 if rxn2.match == None and rxn2.id.lower() == rxn1.id.lower()]) > 0:
            id_counter += 1
            match_m2 += [rxn2 for rxn2 in rxn_list2 if rxn2.match == None and rxn2.id.lower() == rxn1.id.lower()]

        # Compare names
        elif rxn1.name != None and clean_name in rxn_list2_clean_names: 
            name_counter += 1
            match_m2 += rxn_list2_clean_names_map[clean_name] 

        if len(match_m2) >= 1:
            # Check the compartments
            compartment_match = [r for r in match_m2 if set([c.id.lower() for c in r.compartments]) == set([c.id.lower() for c in rxn1.compartments]) or set([c.name.lower() for c in r.compartments]) == set([c.name.lower() for c in rxn1.compartments])]
            if len(compartment_match) == 1:
                match_m2 = compartment_match[0]
            elif len(compartment_match) > 1:
                match_m2 = compartment_match
            else:
                match_m2 = []

            if len(match_m2) >= 1:
                rxn1.match = match_m2
                match_m2.match = rxn1
                common_reactions.append((rxn1,match_m2))           

    # Reactions in rxn_list1 with more than one match in rxn_list2
    more_than_one_match_rxns = [rxn_pair for rxn_pair in common_reactions if isinstance(rxn_pair[1],list)]
    if len(more_than_one_match_rxns) > 0 and warnings:
         print 'WARNING! ',len(more_than_one_match_rxns),' reactions in ',rxn_list1_id,' have more than one match in ',rxn_list2_id
         for rxn_pair in  more_than_one_match_rxns:
             print '\t', (((rxn_list1_id,rxn_pair[0].id),(rxn_list2_id,rxn_pair[1].id)))

    if stdout_msgs:
        print '\nSummary of the compare_reactiongs:'
        print '\tTotal # of reactiongs in %s = %i '%(rxn_list1_id,len(rxn_list1)) 
        print '\tTotal # of reactiongs in %s = %i '%(rxn_list2_id,len(rxn_list2)) 
        print '\tTotal # of common reactions = ',len(common_reactions)
        print '\t\t# of reactiongs matched by ModelSeed id = ',ModelSeed_counter
        print '\t\t# of reactiongs matched by Kegg id = ',Kegg_counter
        print '\t\t# of reactiongs matched by name = ',name_counter
        print '\t\t# of reactiongs matached by model id = ',id_counter
        print '\tTotal # of reactiongs unique to %s = %i '%(rxn_list1_id,len([c for c in rxn_list1 if c.match != None])) 
        print '\tTotal # of reactiongs unique to %s = %i\n'%(rxn_list2_id,len([c for c in rxn_list2 if c.match != None])) 
    
    return common_reactions


def compare_models(model1,model2,stdout_msgs = True):
    """
    This compares two metabolic models. It comapres both reactions and compounds. 
    by comparing their id, Kegg_id and ModelSeed_id, names or formula (for compunds only)

    INPUTS:
    ------
         model1: Model 1
         model2: Model 2
    stdout_msgs: Can be True or False showing whether details of comparison should be
                 written in the output


    OUTPUTS:
    --------
    common_compounds: A list of tuples where the first element is an object of type reaction
                      referring to a compound in cmp_list 1 and the second element is an object 
                      of list of objects of type compound containing the corresponding compounds
                      in cmp_list2 
    common_reactions: A list of tuples where the first element is an object of type reaction
                      referring to a compound in rxn_list 1 and the second element is an object 
                      of list of objects of type reaction containing the corresponding reactions
                      in rxn_list2 
               match: This function also assigns to each compound and reaction object a new field 
                      termed match. If match is None it means that reaction or compound has no 
                      match in other models. Otherwise it is a list of objects of type compound
                      (for compounds) or reactions (for reactions) that match to this compound
                      or reaction in other models 
    """
    if not isinstance(strdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
       
    # Try to get the ModelSeed ids 
    if stdout_msgs:
        print '\n---------- Compre two metabolic models -----------\n'

    # Get ModelSeed ids for reactions and metabolites
    get_ModelSeed_ids(model1,stdout_msgs = stdout_msgs)
    get_ModelSeed_ids(model2,stdout_msgs = stdout_msgs)

    if stdout_msgs:
        print '\n'
        print len([r for r in model1.reactions if r.ModelSeed_id != None]),' (out of ',len(model1.reactions),') reactions have ModelSeed id in ',model1.id
        print len([c for c in model1.compounds if c.ModelSeed_id != None]),' (out of ',len(model1.compounds),') compounds have ModelSeed id in ',model1.id
        print len([r for r in model2.reactions if r.ModelSeed_id != None]),' (out of ',len(model2.reactions),') reactions have ModelSeed id in ',model2.id
        print len([c for c in model2.compounds if c.ModelSeed_id != None]),' (out of ',len(model2.compounds),') compounds have ModelSeed id in ',model2.id
        print '\n'

    # List of compartment ids
    model1_comparts = [compart.id for compart in model1.compartments]
    model2_comparts = [compart.id for compart in model2.compartments]

    # Compare compounds
    common_compounds = compare_compounds(cmp_list1 = model1.compounds,cmp_list2 = model2.compounds,compart_list1 = model1_comparts, compart_list2 = model2_comparts, stdout_msgs = stdout_msgs)

    # Compare reactions
    common_reactions = compare_reactions(rxn_list1 = model1.reactions,rxn_list2 = model2.reactions,compart_list1 = model1_comparts, compart_list2 = model2_comparts, stdout_msgs = stdout_msgs)

    if stdout_msgs:
        print '\n---------- Done with comparing models -----------\n'

    return [common_compounds,common_reactions]

