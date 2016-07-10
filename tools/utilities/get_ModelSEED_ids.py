from __future__ import division
import sys, time
sys.path.append('../../')
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
from remove_non_alphanumeric import remove_non_alphanumeric
from models.ModelSEED.ModelSEED_cpds_master import cpds_master as ModelSEED_cpds 
from models.ModelSEED.ModelSEED_cpds_name_aliases import cpds_by_clean_name_underline as ModelSEED_cpds_by_clean_name 
from models.ModelSEED.ModelSEED_cpds_KEGG_aliases import cpds_by_KEGG_id as ModelSEED_cpds_by_KEGG_id 
from models.ModelSEED.ModelSEED_cpds_formula_aliases import  cpds_by_formula as ModelSEED_cpds_by_formula 
from models.ModelSEED.ModelSEED_cpds_BiGG_aliases import  cpds_by_clean_BiGG_id  as ModelSEED_cpds_by_clean_BiGG_id 
from models.ModelSEED.ModelSEED_rxns_master import rxns_master as ModelSEED_rxns
from models.ModelSEED.ModelSEED_rxns_name_aliases import rxns_by_clean_name_underline as ModelSEED_rxns_by_clean_name
from models.ModelSEED.ModelSEED_rxns_KEGG_aliases import rxns_by_KEGG_id as ModelSEED_rxns_by_KEGG_id
from models.ModelSEED.ModelSEED_rxns_BiGG_aliases import rxns_by_clean_BiGG_id as ModelSEED_rxns_by_clean_BiGG_id
from models.ModelSEED.ModelSEED_rxns_EC_number_aliases import rxns_by_EC_number as ModelSEED_rxns_by_EC_number
import re

"""
This file contains the folloiwng functions:

FUNCTIONS:
---------
   remove_compartment: Removes compartment ids from metabolite or reaction names or ids
   match_rxn_eqn(rxn): Checks if a reaction equaton matches with those in the ModelSEED. 
get_cpds_ModelSEED_id: Tries to get the ModelSEED ids for one or more compounds
get_rxns_ModelSEED_id: Tries to get the ModelSEED ids for one or more reactions

Ali R. Zomorrodi - Segre Lab @ Boston University
Last updated: 04-30-2016 
"""

#---------- Global variables -------------
#--- compounds ---
# All ModelSEED ids in kegg
ModelSEED_cpd_ids = ModelSEED_cpds.keys()

# Names of compounds in the ModelSEED where all non-alpahnumeric characters are removed 
ModelSEED_cpd_clean_names = ModelSEED_cpds_by_clean_name.keys() 

ModelSEED_cpd_names_non_unique = [n for n in ModelSEED_cpd_clean_names if len(ModelSEED_cpds_by_clean_name[n]) > 1]
if len(ModelSEED_cpd_names_non_unique) > 0:
    print '**WARNING! {} compound names in ModelSEED correspond to more than one ModelSEED id'.format(len(ModelSEED_cpd_names_non_unique))

# KEGG ids for compounds in the ModelSEED 
ModelSEED_cpd_KEGG_ids = ModelSEED_cpds_by_KEGG_id.keys() 

# Formulas of compounds in the ModelSEED where all non-alpahnumeric characters are removed 
ModelSEED_cpd_formulas = ModelSEED_cpds_by_formula.keys() 

ModelSEED_cpd_formulas_non_unique = [f for f in ModelSEED_cpd_formulas if len(ModelSEED_cpds_by_formula[f]) > 1]
if len(ModelSEED_cpd_formulas_non_unique) > 0:
    print '**WARNING! {} compound formulas in ModelSEED correspond to more than one ModelSEED id'.format(len(ModelSEED_cpd_formulas_non_unique))

# BiGG ids of ModelSEED compounds 
ModelSEED_cpds_clean_BiGG_ids = ModelSEED_cpds_by_clean_BiGG_id.keys() 

#--- reactions ---
# All ModelSEED ids in kegg
ModelSEED_rxn_ids = ModelSEED_rxns.keys()

# Names of reactions in the ModelSEED where all non-alpahnumeric characters are removed 
ModelSEED_rxn_clean_names = ModelSEED_rxns_by_clean_name.keys() 

ModelSEED_rxn_names_non_unique = [n for n in ModelSEED_rxn_clean_names if len(ModelSEED_rxns_by_clean_name[n]) > 1]
if len(ModelSEED_rxn_names_non_unique) > 0:
    print '**WARNING! {} reaction names in ModelSEED correspond to more than one ModelSEED id'.format(len(ModelSEED_rxn_names_non_unique))

# KEGG ids for reactions in the ModelSEED 
ModelSEED_rxn_KEGG_ids = ModelSEED_rxns_by_KEGG_id.keys() 

# BiGG ids of ModelSEED compounds 
ModelSEED_rxns_clean_BiGG_ids = ModelSEED_rxns_by_clean_BiGG_id.keys()

# EC numbers of ModelSEED compounds 
ModelSEED_rxns_EC_numbers = ModelSEED_rxns_by_EC_number.keys()

# ModelSEED reactions by the list of compounds
ModelSEED_rxns_by_cpdsList = dict([(k,[]) for k in [tuple(sorted([c for c in ModelSEED_rxns[rxn]['stoichiometry'].keys()], key = lambda x: (x[0],x[1]))) for rxn in ModelSEED_rxns.keys() if ModelSEED_rxns[rxn]['stoichiometry'] != None]])
for rxn in [r for r in ModelSEED_rxns.keys() if ModelSEED_rxns[r]['stoichiometry'] != None]:
    cpds_list = tuple(sorted(ModelSEED_rxns[rxn]['stoichiometry'].keys(), key = lambda x: (x[0],x[1])))
    ModelSEED_rxns_by_cpdsList[cpds_list].append(rxn)
print '**WARNING! {} compounds list in the ModelSEED database  correspond to more than ModelSEED rxn'.format(len([c for c in ModelSEED_rxns_by_cpdsList.keys() if len(ModelSEED_rxns_by_cpdsList[c]) != 1]))

# ModelSEED reactions by the list of compounds except for h (cpd00067) and pi (cpd00009). 
# This is to search for reactions in the model that may not be balanced (may have missing 
# h and pi)
# Keys of the dictionary
ModelSEED_rxns_by_cpdsList_nohp = dict([(k,[]) for k in [tuple(sorted([c for c in ModelSEED_rxns[rxn]['stoichiometry'].keys() if c[0] not in ['cpd00067','cpd00009']], key = lambda x: (x[0],x[1]))) for rxn in ModelSEED_rxns.keys() if ModelSEED_rxns[rxn]['stoichiometry'] != None]])
for rxn in [r for r in ModelSEED_rxns.keys() if ModelSEED_rxns[r]['stoichiometry'] != None]:
    cpds_list = tuple(sorted([c for c in ModelSEED_rxns[rxn]['stoichiometry'].keys() if c[0] not in ['cpd00067','cpd00009']], key = lambda x: (x[0],x[1])))
    ModelSEED_rxns_by_cpdsList_nohp[cpds_list].append(rxn)
print '**WARNING! {} compounds list excluding h and pi in the ModelSEED database correspond to more than rxn'.format(len([c for c in ModelSEED_rxns_by_cpdsList_nohp.keys() if len(ModelSEED_rxns_by_cpdsList_nohp[c]) > 1]))

# Check how many stoichiometires correspond to more than one reaction
ModelSEED_rxns_by_stoic = dict([(tuple(ModelSEED_rxns[rxn]['stoichiometry'].items()),[]) for rxn in ModelSEED_rxns.keys() if ModelSEED_rxns[rxn]['stoichiometry'] != None])
for rxn in [r for r in ModelSEED_rxns.keys() if ModelSEED_rxns[r]['stoichiometry'] != None]:
    stoic = tuple(ModelSEED_rxns[rxn]['stoichiometry'].items())
    ModelSEED_rxns_by_stoic[stoic].append(rxn)
ModelSEED_rxn_stoic_non_unique = [s for s in  ModelSEED_rxns_by_stoic.keys() if len(ModelSEED_rxns_by_stoic[s]) > 1]
if len(ModelSEED_rxn_stoic_non_unique) > 0:
    print 'WARNING! {} stoichiometries in ModelSEED database correspond to more than ModelSEED rxn\n'.format(len(ModelSEED_rxn_stoic_non_unique))

def remove_compartment(input_string,compartments_info = None):
    """
    This functions removes the compartment names and/or ids from the input strings,
    which can be compound or reaction names or ids 

    INPUTS:
    ------
         input_string: A string or a list of strings
    compartments_info: A list of strings containing the names or ids of compartments that
                       appear in the list of strings..
                       If not input is proivded a number of pre-set compartment ids
                       (including _c, _e, _m, _n, _p, _r, _v, _x, [c], [e], [m], [n], [r], [v], 
                       [x]]) are used.

    OUTPUTS:
    --------
    A list of strings from which compartments have been removed 
    """
    converted_string = []

    # If the input is a string convert it to a list with one element
    if isinstance(input_string,str):
        input_string = [input_string]

    elif not isinstance(input_string,list):
        raise userError('Invluad input for remove_compartment! String or list of strings expected.' + type(input_string) + ' entered')

    for s in input_string:
        # First some pre-set patterns 
        pattern = '_[c,e,g,m,n,p,r,v,x]0$|_[c,e,g,m,n,p,r,v,x]$|\[[c,e,p,m,x,n,v,g]\]$|\[[c,e,g,m,n,p,r,v,x]0\]$|\[[c,e,p,m,x,n,v,g]\]$|\[[c,e,g,m,n,p,r,v,x]\]$'
        if compartments_info != None:
            for compt in compartments_info:
                pattern += '|_' + compt + '$|\[' + compt + '\]$'          

        s = re.sub(pattern,'',s)

        converted_string.append(s)

    if len(converted_string) == 1:
        converted_string = converted_string[0]

    return converted_string   

def match_rxn_eqn(rxn):
    """
    Checks if a reaction equaton matches with those in the ModelSEED. This function
    works only for reactions where each compound has a unique ModelSEED id
  
    INPUTS:
    -------
                              rxn: A reaction object

    INPUTS (As global variables)
    ---------------------------:
                     ModelSEED_rxns: ModelSEED rxns by id
         ModelSEED_rxns_by_cpdsList: A dictionary with keys being the list of compounds 
                                     participating in that reaction and values being the 
                                     ModelSEED id of reaction
    ModelSEED_rxns_by_cpdsList_nohp: Same as ModelSEED_rxns_by_cpdsList except that h and pi
                                     have been removed from the list of compounds

    OUTPUTS:
    --------
         rxn_ModelSEED_id: A list containing the identified ModelSEED ids.
    ModelSEED_id_found_by: A string showing whether the identified ModelSEED id was 
                           found by matching the equation or by matching the equation
                           with removed h and pi

    The reactions in a metabolic model can have any of the following compartments:
    1. id: e, [e], e0  ,  Name: Extracellular, Extra_cellular, Extra_organism
    2. id: p, [p], p0  ,  Name: Periplasm
    3. id: c, [c], c0  ,  Name: Cytosol, Cytoplasm
    4. id: m, [m], m0  ,  Name: Mitochondria
    5. id: r, [r], r0  ,  Name: Endoplasmic Reticulum
    6. id: g, [g], g0  ,  Name: Golgi, Golgi Apparatus
    7. id: n, [n], n0  ,  Name: Nucleus
    8. id: x, [x], x0  ,  Name: Peroxisome
    9. id: v, [v], v0  ,  Name: Vacuole

    Since in the compounds participating in reactions of the ModelSEED appear with only tow compartments of 
    either c0 or eo, we follow the rules below to convert the stoichiomery of the reaction in the model to 
    a form consistent with those in ModelSEED (e.g., 'stoichiometry': {('cpd00794', 'c0'): -1.0, ('cpd00001', 'c0'): 1.0) 
    a) If all reactions occurs in the same compartment (no matter what that compartment is), the the compartment 
       is converted to c0
    b) If some compounds in the reaciton have one compartment and the others have another compartment (typically there 
       are only two compartments in total):
       - If one compartment is cytosol (c0), convert the other to e0
       - If one compartment is p0 and the other e0, convert p0 to c0
    """
    rxn_ModelSEED_id = []
    ModelSEED_id_found_by = 'No match'

    # If each participating compound in this reaction has a unique ModelSEED id
    if len(rxn.ModelSEED_id) != 1 and len([c for c in rxn.compounds if len(c.ModelSEED_id) == 1]) == len(rxn.compounds):

        # A map converting the compartment of each compound in the reaciton to c0 or e0
        compart_map = {}
    
        # Compartments of compounds participating in this reaction
        rxn_comparts = list(set([(ct.id.lower(),ct.name.lower()) for ct in rxn.compartments]))
        if len(rxn_comparts) == 1:
            compart_map[rxn_comparts[0][0]] = 'c0' 
     
        elif len(rxn_comparts) == 2:
            # If one compartment is cytosol, set the other to e0
            c_compart = list(set([c[0] for c in rxn_comparts if c[0] in ['c', '[c]', 'c0'] or c[1].lower() in ['cytosol','cytoplasm']]))
            if len(c_compart) == 1: 
                # Id of the cytosl compartment in the model
                c_compart = c_compart[0][0] 
                compart_map[c_compart] = 'c0'
    
                # The other compartment
                compart_map[[c[0] for c in rxn_comparts if c[0] != c_compart][0]] = 'e0'
    
            # If one compartment is e0 and the other to p0, set p0 to c0
            else:
                e_compart = list(set([c[0] for c in rxn_comparts if c[0] in ['e', '[e]','e0'] or 'extra' in c[1].lower()]))
                p_compart = list(set([c[0] for c in rxn_comparts if c[0] in ['p','[p]', 'p0'] or c[1].lower() ==  'periplasm']))
    
                if len(e_compart) == 1 and len(p_compart) == 1:
                    compart_map[e_compart[0][0]] = 'e0'
                    compart_map[p_compart[0][0]] = 'c0'
    
        #print 'rxn id = {} , rxn_comparts = {} , compart_map = {}'.format(rxn.id, rxn_comparts, compart_map)
    
        if compart_map != {}:
    
            # stoichiometry of the rxn in the model with a ModelSEED consistent format
            rxn_stoic = {}
            for cpd in rxn.compounds:
                rxn_stoic[(cpd.ModelSEED_id[0], compart_map[cpd.compartment.id])] = rxn.stoichiometry[cpd] 
         
            # stoichiometry of the rxn in the model with a ModelSEED consistent format while 
            # removing h and pi from the list of compounds
            rxn_stoic_nohp = {}
            for cpd in [c for c in rxn.compounds if c.ModelSEED_id[0] not in ['cpd00067','cpd00009']]:
                rxn_stoic_nohp[(cpd.ModelSEED_id[0], compart_map[cpd.compartment.id])] = rxn.stoichiometry[cpd] 
         
            # Keys of rxn_stoic sorted first with cpd id and next with compartment id
            sorted_cpds = tuple(sorted(rxn_stoic.keys(), key = lambda x: (x[0],x[1])))
            sorted_cpds_nohp = tuple(sorted(rxn_stoic_nohp.keys(), key = lambda x: (x[0],x[1])))
        
            # Now search this in ModelSEED rxns
            if sorted_cpds in ModelSEED_rxns_by_cpdsList.keys():
                # Now check if the stoichiometric coefficients match. Note that here we check
                # check both forward and backward stoichiometries as the same reaction may be 
                # written in reverse directions in the model and in the ModelSEED
                for ModelSEED_rxn in ModelSEED_rxns_by_cpdsList[sorted_cpds]:
                    if len([cpd for cpd in sorted_cpds if rxn_stoic[cpd] == ModelSEED_rxns[ModelSEED_rxn]['stoichiometry'][cpd]]) == len(sorted_cpds) or len([cpd for cpd in sorted_cpds if -rxn_stoic[cpd] == ModelSEED_rxns[ModelSEED_rxn]['stoichiometry'][cpd]]) == len(sorted_cpds):
                        rxn_ModelSEED_id.append(ModelSEED_rxn)
                        ModelSEED_id_found_by = 'equation'

            # Search in ModelSEED with h and pi removed form reactions (As the reaction in 
            # the model may not be balanced and have missing h and pi)
            elif sorted_cpds_nohp in ModelSEED_rxns_by_cpdsList_nohp.keys():
                for ModelSEED_rxn in ModelSEED_rxns_by_cpdsList_nohp[sorted_cpds_nohp]:
                    if len([cpd for cpd in sorted_cpds_nohp if rxn_stoic_nohp[cpd] == ModelSEED_rxns[ModelSEED_rxn]['stoichiometry'][cpd]]) == len(sorted_cpds) or len([cpd for cpd in sorted_cpds_nohp if -rxn_stoic[cpd] == ModelSEED_rxns[ModelSEED_rxn]['stoichiometry'][cpd]]) == len(sorted_cpds_nohp):
                        rxn_ModelSEED_id.append(ModelSEED_rxn)
                        ModelSEED_id_found_by = 'equation nohp'
        
    else: 
        print '**WARNING! No matching equation foound for reaction {} becaue not all participating compounds have a unique ModelSEED id'.format(rxn.id)

    return rxn_ModelSEED_id, ModelSEED_id_found_by

def get_cpds_ModelSEED_id(cpds_list, compart_list = None, warnings = True, stdout_msgs = True):
    """
    This functions finds the ModelSEED id for a compound or a list of compounds cpd

    INPUTS:
    -------
    cpds_list: 
    A list (tuple) of instances of type compound 

    compart_list: 
    A list of string containing the compartment ids for compounds  in cpds_list

    warnings: 
    Can be True or False showing whether warnings should be written in the output

    stdout_msgs: 
    Can be True or False showing whether details of comparison should be written in the output

    OUTPUTS:
    --------
    If a ModelSEED id is found it is added to the ModelSEED_id field of the object compound
    """
    if not isinstance(stdout_msgs,bool): 
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool): 
        raise TypeError('warnings must be True or False')

    # List of all compartments in the input compounds
    if compart_list == None:
        compart_list = list(set([c.compartment.id for c in cpds_list if c.compartment != None]))

    # counters showing how many of common reactions were found by comparing their
    # KEGG id, name or model id
    id_counter = 0
    KEGG_counter = 0
    name_counter = 0
    formula_counter = 0
    BiGG_counter = 0

    # The following parameter shows how the ModelSEED id was found for each compounds
    ModelSEED_id_found_by = {}

    for cpd in [c for c in cpds_list if len(c.ModelSEED_id) > 0]:
        ModelSEED_id_found_by[cpd.id] = 'already assigned' 
    for cpd in [c for c in cpds_list if len(c.ModelSEED_id) == 0]:
        ModelSEED_id_found_by[cpd.id] = 'No match'

    for cpd in [c for c in cpds_list if c.ModelSEED_id == []]:
        # Remove non-alphanumeric characters and compartment ids from the id and replace 
        # all upper case letters with lowercase for name, id and formula
        clean_id = remove_non_alphanumeric(remove_compartment(cpd.id,compartments_info = compart_list)).lower()
        if cpd.name != '':
            clean_name = remove_non_alphanumeric(remove_compartment(cpd.name,compartments_info = compart_list)).lower() 

        ModelSEED_id_byID = []
        ModelSEED_id_byBiGG = []
        ModelSEED_id_byKEGG = []
        ModelSEED_id_byName = []
        ModelSEED_id_byFormula = []
 
        # Search by id
        if clean_id in ModelSEED_cpd_ids:
            ModelSEED_id_byID = [clean_id] 
            id_counter += 1

        else:
            # Search by BiGG id
            if clean_id in ModelSEED_cpds_clean_BiGG_ids:
                ModelSEED_id_byBiGG = ModelSEED_cpds_by_clean_BiGG_id[clean_id]
                BiGG_counter += 1
  
            # Search by KEGG_id, if it is available 
            if cpd.KEGG_id != []:
                for ckid in cpd.KEGG_id: 
                    if ckid in ModelSEED_cpd_KEGG_ids: 
                        ModelSEED_id_byKEGG += ModelSEED_cpds_by_KEGG_id[cpd.KEGG_id]
                ModelSEED_id_byKEGG = list(set(ModelSEED_id_byKEGG))
                if len(ModelSEED_id_byKEGG) > 0:
                    KEGG_counter += 1

            # Search by name where all non-alphanumeric characters are removed 
            if cpd.name != '':
                clean_names = [remove_non_alphanumeric(remove_compartment(cpd.name,compartments_info = compart_list)).lower()] + [remove_non_alphanumeric(remove_compartment(n,compartments_info = compart_list)).lower() for n in cpd.name_aliases]
                clean_names_intersect = list(set(clean_names).intersection(set(ModelSEED_cpd_clean_names)))
                if len(clean_names_intersect) > 0:
                    name_counter += 1
                    for cname in clean_names_intersect:
                        ModelSEED_id_byName = ModelSEED_cpds_by_clean_name[cname]

            # Search by formula 
            if cpd.formula != '' and cpd.formula in ModelSEED_cpd_formulas: 
                ModelSEED_id_byFormula = ModelSEED_cpds_by_formula[cpd.formula] 
                formula_counter += 1

        #--- Assign an accurate ModelSEED id according to searches above ---
        # If ModelSEED_id_byID is not '', you've found a unique ModelSEED id
        if ModelSEED_id_byID != []:
            cpd.ModelSEED_id = ModelSEED_id_byID
            ModelSEED_id_found_by[cpd.id] = 'model id' 

        # If only one KEGG id was found, you're done
        elif len(ModelSEED_id_byKEGG) == 1:
            cpd.ModelSEED_id = ModelSEED_id_byKEGG
            ModelSEED_id_found_by[cpd.id] = 'KEGG id' 
  
        # If only one BiGG id was found, you're done
        elif len(ModelSEED_id_byBiGG) == 1:
            cpd.ModelSEED_id = ModelSEED_id_byBiGG
            ModelSEED_id_found_by[cpd.id] = 'BiGG id' 
    
        # If a ModelSEED id by formula has found, you're done
        elif len(ModelSEED_id_byFormula) == 1:
            cpd.ModelSEED_id = ModelSEED_id_byFormula
            ModelSEED_id_found_by[cpd.id] = 'formula' 
        else:
            # Find the intersection of all searches
            if len([set(lst) for lst in [ModelSEED_id_byKEGG,ModelSEED_id_byBiGG,ModelSEED_id_byName] if len(lst) > 0]) > 0:
                ModelSEED_id_intersect = list(set.intersection(*[set(lst) for lst in [ModelSEED_id_byKEGG,ModelSEED_id_byBiGG,ModelSEED_id_byName] if len(lst) > 0])) 
                if len(ModelSEED_id_intersect) == 0:
                    cpd.ModelSEED_id = ModelSEED_id_intersect
                    ModelSEED_id_found_by[cpd.id] = 'intersection' 

                # Otherwise find their union
                else:
                    cpd.ModelSEED_id = list(set(ModelSEED_id_byKEGG + ModelSEED_id_byBiGG + ModelSEED_id_byName))
                    ModelSEED_id_found_by[cpd.id] = 'union' 
 
        # Replace name, KEGG_id and formula from the ModelSEED, if only one ModelSEED_id is found 
        if len(cpd.ModelSEED_id) == 1:
            # If the compound does not already have a name, use the first name
            # as the name and the rest as name_aliases
            if cpd.name == '':
                cpd.name = ModelSEED_cpds[cpd.ModelSEED_id]['name']
                cpd.name_aliases = ModelSEED_cpds[cpd.ModelSEED_id[0]]['name_aliases']
            # Otherwise assign the ModelSEED names as name_aliases
            else:
                cpd.name_aliases = [ModelSEED_cpds[cpd.ModelSEED_id[0]]['name']] + ModelSEED_cpds[cpd.ModelSEED_id[0]]['name_aliases']

            # Assign KEGG id, BiGG_id and formula
            if cpd.KEGG_id == [] and ModelSEED_cpds[cpd.ModelSEED_id[0]]['KEGG_id'] != None:
                cpd.KEGG_id = ModelSEED_cpds[cpd.ModelSEED_id[0]]['KEGG_id']
            if cpd.BiGG_id == [] and ModelSEED_cpds[cpd.ModelSEED_id[0]]['BiGG_id'] != None:
                cpd.BiGG_id = ModelSEED_cpds[cpd.ModelSEED_id[0]]['BiGG_id']
            if cpd.formula == '' and ModelSEED_cpds[cpd.ModelSEED_id[0]]['formula'] != None:
                cpd.formula = ModelSEED_cpds[cpd.ModelSEED_id[0]]['formula']

    for c in cpds_list:
        c.ModelSEED_id_found_by = ModelSEED_id_found_by[c.id]

    # Compounds with more than one ModelSEED id
    more_than_one_ModelSEED_id_cpds = [c for c in cpds_list if len(c.ModelSEED_id) > 1]

    # Compounds with more than one ModelSEED id
    no_ModelSEED_cpds = [c.id for c in cpds_list if c.ModelSEED_id == []]

    # ModelSEED ids corresponding to more than one compound
    all_ModelSEED_ids = list(set([mid for c in cpds_list if len(c.ModelSEED_id) > 0 for mid in c.ModelSEED_id]))
    ModelSEED_id_cpd_map = dict([(mid,[]) for mid in all_ModelSEED_ids])
    for mid in all_ModelSEED_ids:
        mid_cpds = [c for c in cpds_list if mid in c.ModelSEED_id]

        # Remove compartment id from name (compartmentid appears in cpd names in ModelSEED
        # kbase sbml model)
        mid_cpds_names = [remove_compartment(c.name, compartments_info = compart_list) for c in mid_cpds]

        # Consider only the ones whose name appears only nce in mid_rxns_names. This is to 
        # avoid reporting the same compound that appears in more than one compartment.
        ModelSEED_id_cpd_map[mid] = [c.id for c in mid_cpds if mid_cpds_names.count(remove_compartment(c.name, compartments_info = compart_list)) == 1]
    mid_non_unique_cpds = [(mid,ModelSEED_id_cpd_map[mid]) for mid in ModelSEED_id_cpd_map.keys() if len(ModelSEED_id_cpd_map[mid]) > 1]

    if warnings:
        if len(more_than_one_ModelSEED_id_cpds) > 0 and warnings:
            print '\nWARNING! More than one ModelSEED id was found for {} compounds including: {}'.format(len(more_than_one_ModelSEED_id_cpds), [(c.id,c.ModelSEED_id) for c in more_than_one_ModelSEED_id_cpds])

        if len(no_ModelSEED_cpds) > 0:
            print '\n**WARNINGS! No ModelSEED id was found for the following {} compounds: {}\n'.format(len(no_ModelSEED_cpds),no_ModelSEED_cpds)

        if len(mid_non_unique_cpds) > 0:
            print '**WARNING! The following {} ModelSEED ids correspond to more than one compound: {}'.format(len(mid_non_unique_cpds),mid_non_unique_cpds)

    if stdout_msgs:
        print '\nSummary of the search for the ModelSEED id of compounds:'
        print '\tTotal # of compounds with a ModelSEED id = {} (out of {})'.format(len([c for c in cpds_list if c.ModelSEED_id != []]), len(cpds_list))
        print '\t\t# of compounds with a matached model id = ',id_counter
        print '\t\t# of compounds with a match KEGG id = ',KEGG_counter
        print '\t\t# of compounds with a matched BiGG id = ',BiGG_counter
        print '\t\t# of compounds with a matched name = ',name_counter
        print '\t\t# of compounds with a matached formula = ',formula_counter

        print '\n\t\t# of compounds with ModelSEED ids found by model id = ',len([c for c in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[c] == 'model id'])
        print '\t\t# of compounds with ModelSEED ids found by KEGG id = ',len([c for c in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[c] == 'KEGG id'])
        print '\t\t# of compounds with ModelSEED ids found by BiGG id = ',len([c for c in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[c] == 'BiGG id'])
        print '\t\t# of compounds with ModelSEED ids found by formula = ',len([c for c in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[c] == 'formula id'])
        print '\t\t# of compounds with ModelSEED ids found by intersecton of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([c for c in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[c] == 'intersection'])
        print '\t\t# of compounds with ModelSEED ids found by union of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([c for c in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[c] == 'union'])

        print '\tTotal # of compounds with no ModelSEED id = {} (out of {})\n'.format(len([c for c in cpds_list if c.ModelSEED_id == []]),len(cpds_list))



def get_rxns_ModelSEED_id(rxns_list, compart_list = None, warnings = True, stdout_msgs = True): 
    """
    This functions finds the ModelSEED id for a reaction or a list of reactions

    INPUTS:
    -------
    rxns_list: 
    A list (tuple) of type reaction or a list of reactions

    compart_list: 
    A list of string containing the compartment ids for reactions in rxns_list

    warnings: 
    Can be True or False showing whether warnings should be written  in the output

    stdout_msgs: 
    Can be True or False showing whether details of comparison should be  written in the output

    OUTPUTS:
    --------
    If a ModelSEED id is found it is added to the ModelSEED_id field of the object rxn
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be True or False')

    # List of all compartments in the input compounds
    if compart_list == None:
        compart_list = list(set([compart.id for r in rxns_list for compart in r.compartments if r.compartments != None]))

    # counters showing how many of common reactions were found by comparing their
    # ModelSEED id, KEGG id, BiGG id, name, model id or by their equation
    id_counter = 0
    KEGG_counter = 0
    name_counter = 0
    BiGG_counter = 0
    EC_counter = 0 
    eqn_counter = 0
    exch_counter = 0

    # The following parameter shows how the ModelSEED id was found for each compounds
    ModelSEED_id_found_by = {}
    for rxn in [r for r in rxns_list if len(r.ModelSEED_id) > 0]: 
        ModelSEED_id_found_by[rxn.id] = 'already assigned'
    for rxn in [r for r in rxns_list if len(r.ModelSEED_id) == 0]: 
        ModelSEED_id_found_by[rxn.id] = 'No match'

    # Examine non-exchange reactions first
    for rxn in [r for r in rxns_list if r.ModelSEED_id == [] and 'exchange' not in r.reversibility.lower()]:
        # Remove non-alphanumeric characters and compartment ids from the id and replace 
        # all upper case letters with lowercase for name, id and formula
        clean_id = remove_non_alphanumeric(remove_compartment(rxn.id,compartments_info = compart_list)).lower()
        clean_name = remove_non_alphanumeric(remove_compartment(rxn.name,compartments_info = compart_list)).lower()

        ModelSEED_id_byID = []
        ModelSEED_id_byBiGG = []
        ModelSEED_id_byKEGG = []
        ModelSEED_id_byName = []
        ModelSEED_id_byECNumber = []
        ModelSEED_id_byEqn = []

        # Search by id
        if clean_id in ModelSEED_rxn_ids:
            ModelSEED_id_byID = [clean_id]
            id_counter += 1

        else:
            # Search by BiGG id
            if clean_id in ModelSEED_rxns_clean_BiGG_ids:
                ModelSEED_id_byBiGG = ModelSEED_rxns_by_clean_BiGG_id[clean_id]
                BiGG_counter += 1
 
            # Search by KEGG_id, if it is available 
            if rxn.KEGG_id != []:
                for r_kid in rxn.KEGG_id: 
                    if r_kid in ModelSEED_rxn_KEGG_ids: 
                        ModelSEED_id_byKEGG += ModelSEED_rxns_by_KEGG_id[r_kid]
                ModelSEED_id_byKEGG = list(set(ModelSEED_id_byKEGG))
                if len(ModelSEED_id_byKEGG) > 0:
                    KEGG_counter += 1

            # Search by EC_number, if it is available 
            if rxn.EC_numbers != []:
                for r_ec in rxn.EC_numbers: 
                    if r_ec in ModelSEED_rxns_EC_numbers: 
                        ModelSEED_id_byECNumber += ModelSEED_rxns_by_EC_number[r_ec]
                ModelSEED_id_byECNumber = list(set(ModelSEED_id_byECNumber))
                if len(ModelSEED_id_byECNumber) > 0:
                    ECNumber_counter += 1

            # Search by name where all non-alphanumeric characters are removed 
            if rxn.name != '' and clean_name in ModelSEED_rxn_clean_names:
	        ModelSEED_id_byName = ModelSEED_rxns_by_clean_name[clean_name]
                name_counter += 1

            # If no or more than one ModelSEED id was found, then search by equation. This search is performed only if all
            # compounds participating in this reaction have a unique ModelSEED_id
            if len(rxn.ModelSEED_id) != 1 and len([c for c in rxn.compounds if len(c.ModelSEED_id) == 1]) == len(rxn.compounds):
                ModelSEED_id_byEqn, eqn_found_by = match_rxn_eqn(rxn) 
                if ModelSEED_id_byEqn != []:
                    eqn_counter += 1

        #--- Assign an accurate ModelSEED id according to searches above ---
        # If ModelSEED_id_byID is not '', you've found a unique ModelSEED id
        if ModelSEED_id_byID != []:
            rxn.ModelSEED_id = ModelSEED_id_byID
            ModelSEED_id_found_by[rxn.id] = 'model id' 

        else:
            # If only one ModelSEED_id by KEGG id was found, you're done
            if len(ModelSEED_id_byKEGG) == 1:
                rxn.ModelSEED_id = ModelSEED_id_byKEGG
                ModelSEED_id_found_by[rxn.id] = 'KEGG id' 
    
            # If only one ModelSEED by BiGG id was found, you're done
            elif len(ModelSEED_id_byBiGG) == 1:
                rxn.ModelSEED_id = ModelSEED_id_byBiGG
                ModelSEED_id_found_by[rxn.id] = 'BiGG id' 

            # If only one ModelSEED_id by reaction equation was found, you're done
            elif len(ModelSEED_id_byEqn) == 1:
                rxn.ModelSEED_id = ModelSEED_id_byEqn
                ModelSEED_id_found_by[rxn.id] = eqn_found_by 

            # If one ModelSEED id by EC number has found, you're done
            elif len(ModelSEED_id_byECNumber) == 1:
                rxn.ModelSEED_id = ModelSEED_id_byECNumber
                ModelSEED_id_found_by[rxn.id] = 'EC number' 

            else:
                # Find the intersection of all searches
                if len([lst for lst in [ModelSEED_id_byKEGG, ModelSEED_id_byBiGG, ModelSEED_id_byEqn, ModelSEED_id_byECNumber, ModelSEED_id_byName] if len(lst) > 0]) > 0:
                    ModelSEED_id_intersect = list(set.intersection(*[set(lst) for lst in [ModelSEED_id_byKEGG, ModelSEED_id_byBiGG, ModelSEED_id_byEqn, ModelSEED_id_byECNumber, ModelSEED_id_byName] if len(lst) > 0])) 
                    if len(ModelSEED_id_intersect) > 0:
                        rxn.ModelSEED_id = ModelSEED_id_intersect
                        ModelSEED_id_found_by[rxn.id] = 'intersection' 

                    # Otherwise find their union
                    else: 
                        rxn.ModelSEED_id = list(set(ModelSEED_id_byKEGG + ModelSEED_id_byBiGG + ModelSEED_id_byName))
                        ModelSEED_id_found_by[rxn.id] = 'union' 

                    
        # Replace name, KEGG_id and formula from the ModelSEED, if only one ModelSEED_id is found 
        if len(rxn.ModelSEED_id) == 1:
            # If the compound does not already have a name, use the first name
            # as the name and the rest as name_aliases
            if rxn.name == '':
                rxn.name = ModelSEED_rxns[rxn.ModelSEED_id]['name']
                rxn.name_aliases = ModelSEED_rxns[rxn.ModelSEED_id[0]]['name_aliases']
            # Otherwise assign the ModelSEED names as name_aliases
            else:
                rxn.name_aliases = [ModelSEED_rxns[rxn.ModelSEED_id[0]]['name']] + ModelSEED_rxns[rxn.ModelSEED_id[0]]['name_aliases']

            # Assign KEGG id, BiGG_id and EC number 
            if rxn.KEGG_id == [] and ModelSEED_rxns[rxn.ModelSEED_id[0]]['KEGG_id'] != None:
                rxn.KEGG_id = ModelSEED_rxns[rxn.ModelSEED_id[0]]['KEGG_id']
            if rxn.BiGG_id == [] and ModelSEED_rxns[rxn.ModelSEED_id[0]]['BiGG_id'] != None:
                rxn.BiGG_id = ModelSEED_rxns[rxn.ModelSEED_id[0]]['BiGG_id']
            if rxn.EC_numbers == [] and ModelSEED_rxns[rxn.ModelSEED_id[0]]['EC_numbers'] != None:
                rxn.EC_numbers = ModelSEED_rxns[rxn.ModelSEED_id[0]]['EC_numbers']

    # Now examine exchange reactions by checking whether the metabolite participating in 
    # this exchange reaction has a kegg id (cpd.ModelSEEDid). If it did, then assign 
    # EX_cpdModelSEEDid_e0 as the ModelSEED id for the exchange reaction
    for rxn in [r for r in rxns_list if r.ModelSEED_id == [] and 'exchange' in r.reversibility.lower()]:
        # Check whether the compound participating in this exchange reaction has a ModelSEED id
        if len(rxn.compounds[0].ModelSEED_id) == 1:
            exch_counter += 1
            rxn.ModelSEED_id = ['EX_' + rxn.compounds[0].ModelSEED_id[0] + '_e0']
            exch_counter += 1
            ModelSEED_id_found_by[rxn.id] = 'exchange'

            if rxn.name == '':
                rxn.name = ModelSEED_cpds[rxn.compounds[0].ModelSEED_id[0]]['name'] + ' exchange'
            if rxn.name_aliases == []:
                rxn.name_aliases = [na + ' exchange' for na in ModelSEED_cpds[rxn.compounds[0].ModelSEED_id[0]]['name_aliases']]
            if rxn.KEGG_id == [] and ModelSEED_cpds[rxn.compounds[0].ModelSEED_id[0]]['KEGG_id'] != None:
                for kid in ModelSEED_cpds[rxn.compounds[0].ModelSEED_id[0]]['KEGG_id']:
                    rxn.KEGG_id.append('EX_' + kid + '_e0')
            if rxn.BiGG_id == [] and ModelSEED_cpds[rxn.compounds[0].ModelSEED_id[0]]['BiGG_id'] != None:
                for bid in ModelSEED_cpds[rxn.compounds[0].ModelSEED_id[0]]['BiGG_id']:
                    rxn.BiGG_id.append('EX_' + bid + '_e0')

        elif len(rxn.compounds[0].ModelSEED_id) > 1:
            rxn.ModelSEED_id = []
            for sid in rxn.compounds[0].ModelSEED_id:
                rxn.ModelSEED_id.append('EX_' + sid + '_e0')
            exch_counter += 1
            ModelSEED_id_found_by[rxn.id] = 'exchange'

    for r in rxns_list:
        rxn.ModelSEED_id_found_by = ModelSEED_id_found_by[r.id]

    # Reactions with more than one ModelSEED id
    more_than_one_ModelSEED_id_rxns = [c for c in rxns_list if len(c.ModelSEED_id) > 1]

    # Reactions with no ModelSEED id
    no_ModelSEED_rxns = [r.id for r in rxns_list if r.ModelSEED_id == []]

    # ModelSEED ids corresponding to more than one compound
    all_ModelSEED_ids = list(set([mid for r in rxns_list if len(r.ModelSEED_id) > 0 for mid in r.ModelSEED_id]))
    ModelSEED_id_rxn_map = dict([(mid,[]) for mid in all_ModelSEED_ids])
    for mid in all_ModelSEED_ids:
        mid_rxns = [r for r in rxns_list if mid in r.ModelSEED_id]

        # Remove compartment ids from reaction names (compartment ids appear in reaciton 
        # names in the sbml files of the ModelSEED and kbase model)
        mid_rxns_names = [remove_compartment(r.name, compartments_info = compart_list) for r in mid_rxns]

        # Consider only the ones whose name appears only nce in mid_rxns_names. This is to 
        # avoid reporting the same compound that appears in more than one compartment.
        ModelSEED_id_rxn_map[mid] = [r.id for r in mid_rxns if mid_rxns_names.count(remove_compartment(r.name, compartments_info = compart_list)) == 1]

    mid_non_unique_rxns = [(mid,ModelSEED_id_rxn_map[mid]) for mid in ModelSEED_id_rxn_map.keys() if len(ModelSEED_id_rxn_map[mid]) > 1]

    if warnings:
        if len(more_than_one_ModelSEED_id_rxns) > 0:
            print '\nWARNING! More than one ModelSEED id was found for the following {} reactions: {}'.format(len(more_than_one_ModelSEED_id_rxns), [(r.id,r.ModelSEED_id) for r in more_than_one_ModelSEED_id_rxns])

        if len(no_ModelSEED_rxns) > 0 and warnings:
            print '\n**WARNINGS! No ModelSEED id was found for the following {} reactions: {}\n'.format(len(no_ModelSEED_rxns),no_ModelSEED_rxns)

        if len(mid_non_unique_rxns) > 0:
            print '**WARNING! The following {} ModelSEED ids correspond to more than one compound: {}'.format(len(mid_non_unique_rxns),mid_non_unique_rxns)

        if stdout_msgs:
            print '\nSummary of the search for the ModelSEED id of reactions:'
            print '\tTotal # of reactions with a ModelSEED id = {} (out of {})'.format(len([c for c in rxns_list if c.ModelSEED_id != []]), len(rxns_list))
            print '\t\t# of reactions with a matached model id = ',id_counter
            print '\t\t# of reactions with a match KEGG id = ',KEGG_counter
            print '\t\t# of reactions with a matched BiGG id = ',BiGG_counter
            print '\t\t# of reactions with a matched name = ',name_counter
            print '\t\t# of reactions with a matached equation = ',eqn_counter
            print '\t\t# of exchange reactions with a matached compound = ',exch_counter
    
            print '\n\t\t# of reactions with ModelSEED ids found by model id = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'model id'])
            print '\t\t# of reactions with ModelSEED ids found by KEGG id = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'KEGG id'])
            print '\t\t# of reactions with ModelSEED ids found by BiGG id = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'BiGG id'])
            print '\t\t# of reactions with ModelSEED ids found by equation = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'equation'])
            print '\t\t# of reactions with ModelSEED ids found by equation with removed h and pi = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'equation nohp'])
            print '\t\t# of exchange reactions with ModelSEED ids found by their participating compound  = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'exchange'])
            print '\t\t# of reactions with ModelSEED ids found by intersecton of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'intersection'])
            print '\t\t# of reactions with ModelSEED ids found by union of ModelSEED ids matching their KEGG id, BiGG id and names = ',len([r for r in ModelSEED_id_found_by.keys() if ModelSEED_id_found_by[r] == 'union'])
            print '\t\t# of reactions with ModelSEED ids found by intersecton of ModelSEED ids matching their KEGG id, BiGG id and names and those found by matching rxn eqn = ',len([r for r in ModelSEED_id_found_by.keys() if 'intersection with' in ModelSEED_id_found_by[r]])
            print '\t\t# of reactions with ModelSEED ids found by union of ModelSEED ids matching their KEGG id, BiGG id and names and those found by matching rxn eqn = ',len([r for r in ModelSEED_id_found_by.keys() if 'union with' in ModelSEED_id_found_by[r]])
    
            print '\tTotal # of reactions with no ModelSEED id = {} (out of {})\n'.format(len([r for r in rxns_list if r.ModelSEED_id == []]),len(rxns_list))
    

def get_ModelSEED_ids(model,stdout_msgs = True):
    """
    This finds the ModelSEED ids for all reactions and metabolites ia a model 
    INPUTS:
    ------
    model: An object of type model

    OUTPUTS:
    --------
    The ModelSEED ids are assigned to the field ModelSEED_id of a compound or 
    reaction in the model
    """
    print '    Getting ModelSEED ids for compounds ...'
    get_cpds_ModelSEED_id(cpds_list = model.compounds, compart_list = [c.id for c in model.compartments],stdout_msgs = stdout_msgs)
    print '    Getting ModelSEED ids for reactions ...'
    get_rxns_ModelSEED_id(rxns_list = model.reactions, compart_list = [c.id for c in model.compartments],stdout_msgs = stdout_msgs)
 
