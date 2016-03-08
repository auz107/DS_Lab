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
  remove_compartment: Removes compartment ids from metabolite or reaciton names or ids
get_cpd_ModelSEED_id: Tries to get the ModelSEED ids for one or more compounds
get_rxn_ModelSEED_id: Tries to get the ModelSEED ids for one or more reactions

Ali R. Zomorrodi - Segre Lab @ Boston University
Last updated: 02-29-2016 
"""

def remove_compartment(input_string,compartments_info = None):
    """
    This functions removes the compartment names and/or ids from the input strings,
    which can be compound or reaction names or ids 

    INPUTS:
    ------
         input_string: A string or a list of strings
    compartments_info: A list of strings containing the names and ids of compartments.
                       If not input is proivded a number of pre-set compartment ids
                       (including _c, _e, _m, _x, _v, _n, [c], [e], [m], [x], [n]) are used.

    OUTPUTS:
    a string with all non-alaphabetical and non-numerical characters replaced with an underline
    """
    converted_string = []

    # If the input is a string convert it to a list with one element
    if isinstance(input_string,str):
        input_string = [input_string]

    elif not isinstance(input_string,list):
        raise userError('Invluad input for remove_compartment! String or list of strings expected.' + isinstance(input_string) + ' entered')

    for s in input_string:
        # First some pre-set patterns 
        pattern = '_[c,e,p,m,x,n,v,g]0$|_[c,e,p,m,x,n,v,g]$|\[[c,e,p,m,x,n,v,g]\]$|\[[c,e,p,m,x,n,v,g]0\]$'
        if compartments_info != None:
            for compt in compartments_info:
                pattern += '|_' + compt + '$|\[' + compt + '\]$'          

        s = re.sub(pattern,'',s)

        converted_string.append(s)

    if len(converted_string) == 1:
        converted_string = converted_string[0]

    return converted_string   

def get_cpd_ModelSEED_id(cpd_list, compart_list = None, warnings = True, stdout_msgs = True):
    """
    This functions finds the ModelSEED id for a compound or a list of compounds cpd

    INPUTS:
    -------
         cpd_list: An instance or a list of instances of type compound 
     compart_list: A list of string containing the compartment ids for compounds 
                   in cpd_list
         warnings: Can be True or False showing whether warnings should be written  
                   in the output
      stdout_msgs: Can be True or False showing whether details of comparison should be
                   written in the output

    OUTPUTS:
    --------
    If a ModelSEED id is found it is added to the ModelSEED_id field of the object compound
    """
    if not isinstance(stdout_msgs,bool): 
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool): 
        raise TypeError('warnings must be True or False')

    if not isinstance(cpd_list,list):
        cpd_list = [cpd_list]

    # List of all compartments in the input compounds
    if compart_list == None:
        compart_list = list(set([c.compartment.id for c in cpd_list if c.compartment != None]))

    # All ModelSEED ids in kegg
    ModelSEED_cpd_ids = ModelSEED_cpds.keys()

    # Names of compounds in the ModelSEED where all non-alpahnumeric characters are removed 
    ModelSEED_cpd_clean_names = ModelSEED_cpds_by_clean_name.keys() 

    # KEGG ids for compounds in the ModelSEED 
    ModelSEED_cpd_KEGG_ids = ModelSEED_cpds_by_KEGG_id.keys() 

    # Formulas of compounds in the ModelSEED where all non-alpahnumeric characters are removed 
    ModelSEED_cpd_formulas = ModelSEED_cpds_by_formula.keys() 

    # BiGG ids of ModelSEED compounds 
    ModelSEED_cpds_clean_BiGG_ids = ModelSEED_cpds_by_clean_BiGG_id.keys() 

    # counters showing how many of common reactions were found by comparing their
    # KEGG id, name or model id
    KEGG_counter = 0
    name_counter = 0
    formula_counter = 0
    BiGG_counter = 0
    id_counter = 0

    for cpd in [c for c in cpd_list if c.ModelSEED_id == []]:
        # Remove non-alphanumeric characters and compartment ids from the id and replace 
        # all upper case letters with lowercase for name, id and formula
        clean_id = remove_non_alphanumeric(remove_compartment(cpd.id,compartments_info = compart_list)).lower()
        if cpd.name != '':
            clean_name = remove_non_alphanumeric(remove_compartment(cpd.name,compartments_info = compart_list)).lower() 

        ModelSEED_id_byID = ''
        ModelSEED_id_byBiGG = []
        ModelSEED_id_byKEGG = []
        ModelSEED_id_byName = []
        ModelSEED_id_byFormula = []
 
        # Search by id
        if clean_id in ModelSEED_cpd_ids:
            cpd.ModelSEED_id = clean_id 
            id_counter += 1

        else:
            # Search by BiGG id
            if clean_id in ModelSEED_cpds_clean_BiGG_ids:
                ModelSEED_id_byBiGG = ModelSEED_cpds_by_clean_BiGG_id[clean_id] 
                BiGG_counter += 1
  
            # Search by KEGG_id, if it is available 
            if cpd.KEGG_id != [] and cpd.KEGG_id in ModelSEED_cpd_KEGG_ids: 
                ModelSEED_id_byKEGG = ModelSEED_cpds_by_KEGG_id[cpd.KEGG_id] 
                KEGG_counter += 1

            # Search by name where all non-alphanumeric characters are removed 
            if cpd.name != None and clean_name in ModelSEED_cpd_clean_names:
	        ModelSEED_id_byName = ModelSEED_cpds_by_clean_name[clean_name]
                name_counter += 1

            # Search by formula 
            if cpd.formula != None and cpd.formula in ModelSEED_cpd_formulas: 
                ModelSEED_id_byFormula = ModelSEED_cpds_by_formula[clean_formula] 
                formula_counter += 1

        #--- Assign an accurate ModelSEED id according to searches above ---
        # If ModelSEED_id_byID is not '', you've found a unique ModelSEED id
        if ModelSEED_id_byID != '':
            cpd.ModelSEED_id = [ModelSEED_id_byID]
        
        else:
            # If only one KEGG id was found, you're done
            if len(ModelSEED_id_byKEGG) == 1:
                cpd.ModelSEED_id = ModelSEED_id_byKEGG
            else:
                # If only one BiGG id was found, you're done
                if len(ModelSEED_id_byBiGG) == 1:
                    cpd.ModelSEED_id = ModelSEED_id_byBiGG
                else:
                    # Find the intersection of all searches
                    ModelSEED_id_intersect = list(set.intersection(*[set(lst) for lst in [ModelSEED_id_byKEGG,ModelSEED_id_byBiGG,ModelSEED_id_byName,ModelSEED_id_byFormula] if len(lst) > 0])) 

                    # If the intersection is not empty assign it as the ModelSEED id
                    if len( ModelSEED_id_intersec) > 0:
                        cpd.ModelSEED_id =  ModelSEED_id_intersec

                    # Otherwise find their union
                    else:
                        cpd.ModelSEED_id = set(ModelSEED_id_byKEGG + ModelSEED_id_byBiGG + ModelSEED_id_byName + ModelSEED_id_byFormula)

 
        # If a ModelSEED id is found replace the name, KEGG_id and formula from the ModelSEED 
        if cpd.ModelSEED_id != []:
            # If more than one ModelSEED is was found
            if isinstance(cpd.ModelSEED_id,list) and len(cpd.ModelSEED_id) > 1:
                # If more than one ModelSEED id is found, then don't assign the name because it
                # can be confusing. Just specify synonyms as a set of sets
                cpd.synonyms = []
                cpd.formula = []
                cpd.KEGG_id = []
                for seed_id in cpd.ModelSEED_id:
                    cpd.synonyms += ModelSEED_cpds[seed_id]['name']
                    if ModelSEED_cpds[seed_id]['formula'] != None:
                        cpd.formula += ModelSEED_cpds[seed_id]['formula']
                    if ModelSEED_cpds[seed_id]['KEGG_id'] != None:
                        cpd.KEGG_id += ModelSEED_cpds[seed_id]['KEGG_id']

            # If only ModelSEED id was found
            elif isinstance(cpd.ModelSEED_id,list) and len(cpd.ModelSEED_id) == 1:
                cpd.ModelSEED_id = cpd.ModelSEED_id[0]      
                # If the compound does not already have a name, use the first name
                # as the name and the rest as synonyms
                if cpd.name == None or (cpd.name != None and len(cpd.name) == 0):
                    cpd.name = ModelSEED_cpds[cpd.ModelSEED_id]['name'][0]
                    cpd.synonyms = ModelSEED_cpds[cpd.ModelSEED_id]['name'][1:]
                # Otherwise assign the ModelSEED names as synonyms
                else:
                    cpd.synonyms = ModelSEED_cpds[cpd.ModelSEED_id]['name']
                if ModelSEED_cpds[cpd.ModelSEED_id]['formula'] != None:
                    cpd.formula = ModelSEED_cpds[cpd.ModelSEED_id]['formula']
                if ModelSEED_cpds[cpd.ModelSEED_id]['KEGG_id'] != None:
                    cpd.KEGG_id = ModelSEED_cpds[cpd.ModelSEED_id]['KEGG_id']

            # If a ModelSEED has already been assigned (as a string) in the original model 
            else:
                if cpd.name == None or (cpd.name != None and len(cpd.name) == 0):
                    cpd.name = ModelSEED_cpds[cpd.ModelSEED_id]['name'][0]
                    cpd.synonyms = ModelSEED_cpds[cpd.ModelSEED_id]['name'][1:]
                else:
                    cpd.synonyms = ModelSEED_cpds[cpd.ModelSEED_id]['name']
                if ModelSEED_cpds[cpd.ModelSEED_id]['formula'] != None:
                    cpd.formula = ModelSEED_cpds[cpd.ModelSEED_id]['formula']
                if ModelSEED_cpds[cpd.ModelSEED_id]['KEGG_id'] != None:
                    cpd.KEGG_id = ModelSEED_cpds[cpd.ModelSEED_id]['KEGG_id']

    # Compounds with more than one ModelSEED id
    more_than_one_ModelSEED_id_cpds = [c for c in cpd_list if c.ModelSEED_id != [] and isinstance(c.ModelSEED_id,list)]
    if stdout_msgs:
        if len(more_than_one_ModelSEED_id_cpds) > 0 and warnings:
            print '\nWARNING! More than one ModelSEED id was found for %i compounds including:'%(len(more_than_one_ModelSEED_id_cpds))
        for cpd in more_than_one_ModelSEED_id_cpds:
            print '\t',cpd.id,'\t',cpd.ModelSEED_id

        print '\nSummary of the search for the ModelSEED id of compounds:'
        print '\tTotal # of compounds with a ModelSEED id = %i (out of %i)'%(KEGG_counter + name_counter + id_counter + formula_counter,len(cpd_list))
        print '\t\t# of compounds matched by KEGG id = ',KEGG_counter
        print '\t\t# of compounds matched by name = ',name_counter
        print '\t\t# of compounds matached by model id = ',id_counter
        print '\t\t# of compounds matached by formula = ',formula_counter
        print '\tTotal # of compounds with no ModelSEED id = %i (out of %i)\n'%(len(cpd_list) - (KEGG_counter + name_counter + id_counter + formula_counter),len(cpd_list))

def get_rxn_ModelSEED_id(rxn_list, compart_list = None, warnings = True, stdout_msgs = True): 
    """
    This functions finds the ModelSEED id for a reaction or a list of reactions

    INPUTS:
    -------
         rxn_list: An instance of type reaction or a list of reactions
     compart_list: A list of string containing the compartment ids for reactions 
                   in rxn_list
         warnings: Can be True or False showing whether warnings should be written  
                   in the output
      stdout_msgs: Can be True or False showing whether details of comparison should be
                   written in the output

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
        compart_list = list(set([compart.id for r in rxn_list for compart in r.compartments if r.compartments != None]))

    # All ModelSEED ids in kegg
    ModelSEED_rxn_ids = ModelSEED_rxns.keys()

    # Names of reactions in the ModelSEED where all non-alpahnumeric characters are removed 
    ModelSEED_rxn_clean_names = ModelSEED_rxns_by_clean_name.keys() 

    # KEGG ids for reactions in the ModelSEED 
    ModelSEED_rxn_KEGG_ids = ModelSEED_rxns_by_KEGG_id.keys() 

    # counters showing how many of common reactions were found by comparing their
    # ModelSEED id, KEGG id, name or model id
    KEGG_counter = 0
    name_counter = 0
    id_counter = 0
    exch_counter = 0

    if not isinstance(rxn_list,list):
        rxn_list = [rxn_list]

    # Examine non-exchange reacitons first
    for rxn in [r for r in rxn_list if r.ModelSEED_id == [] and 'exchange' not in r.type.lower()]:
        # Remove non-alphanumeric characters and compartment ids from the id and replace 
        # all upper case letters with lowercase for name, id and formula
        clean_id = remove_non_alphanumeric(remove_compartment(rxn.id,compartments_info = compart_list)).lower()
        if rxn.name != '':
            clean_name = remove_non_alphanumeric(remove_compartment(rxn.name,compartments_info = compart_list)).lower() 
        # Search by id 
        if clean_id in ModelSEED_rxn_ids:
            rxn.ModelSEED_id = clean_id 
            id_counter += 1

        # Search by KEGG_id, if it is available 
        elif rxn.KEGG_id != [] and rxn.KEGG_id in ModelSEED_rxn_KEGG_ids: 
            rxn.ModelSEED_id = ModelSEED_rxns_by_KEGG_id[rxn.KEGG_id] 
            KEGG_counter += 1

        # Search by name where all non-alphanumeric characters are removed 
        elif rxn.name != None and clean_name in ModelSEED_rxns_by_clean_name:
	    rxn.ModelSEED_id = ModelSEED_rxns_by_clean_name[clean_name]
            name_counter += 1

        # If a ModelSEED id is found replace the name, KEGG_id and EC_number from the ModelSEED 
        if rxn.ModelSEED_id != []:
            # If more than one ModelSEED id was found
            if isinstance(rxn.ModelSEED_id,list) and len(rxn.ModelSEED_id) > 1:
                # If there are more than one ModelSEED id don't assign the names from ModelSEED
                # because it can be confusing. Just specify the synonyms as a list of lists
                rxn.synonyms = []
                rxn.EC_number = []
                rxn.KEGG_id = []
                for seed_id in rxn.ModelSEED_id:
                    rxn.synonyms += ModelSEED_rxns[seed_id]['name']
                    if ModelSEED_rxns[seed_id]['EC_number'] != None:
                        rxn.EC_number += ModelSEED_rxns[seed_id]['EC_number']
                    if ModelSEED_rxns[seed_id]['KEGG_id'] != None:
                        rxn.KEGG_id += ModelSEED_rxns[seed_id]['KEGG_id']
    
            # If only one ModelSEED id was found
            elif isinstance(rxn.ModelSEED_id,list) and len(rxn.ModelSEED_id) == 1:
                rxn.ModelSEED_id = rxn.ModelSEED_id[0]      
                # If a name has not already been assigned, assign the first element from
                # ModelSEED names to name and the rest to synonyms
                if rxn.name == None or (rxn.name != None and len(rxn.name) == 0):
                    rxn.name = ModelSEED_rxns[rxn.ModelSEED_id]['name'][0]
                    rxn.synonyms = ModelSEED_rxns[rxn.ModelSEED_id]['name'][1:]
                # Otherwise assign the ModelSEED names to synonyms
                else:
                    rxn.synonyms = ModelSEED_rxns[rxn.ModelSEED_id]['name']
                if ModelSEED_rxns[rxn.ModelSEED_id]['EC_number'] != None:
                    rxn.EC_number = ModelSEED_rxns[rxn.ModelSEED_id]['EC_number']
                if ModelSEED_rxns[rxn.ModelSEED_id]['KEGG_id'] != None:
                    rxn.KEGG_id = ModelSEED_rxns[rxn.ModelSEED_id]['KEGG_id']
    
            # If a ModelSEED id has already been assigned (as a string) in the original model
            else: 
                if rxn.name == None or (rxn.name != None and len(rxn.name) == 0):
                    rxn.name = ModelSEED_rxns[rxn.ModelSEED_id]['name'][0]
                    rxn.synonyms = ModelSEED_rxns[rxn.ModelSEED_id]['name'][1:]
                else:
                    rxn.synonyms = ModelSEED_rxns[rxn.ModelSEED_id]['name']
                if ModelSEED_rxns[rxn.ModelSEED_id]['EC_number'] != None:
                    rxn.EC_number = ModelSEED_rxns[rxn.ModelSEED_id]['EC_number']
                if ModelSEED_rxns[rxn.ModelSEED_id]['KEGG_id'] != None:
                    rxn.KEGG_id = ModelSEED_rxns[rxn.ModelSEED_id]['KEGG_id']

    # Now examine exchange reactions by checking whether the metabolite participating in 
    # this exchange reaction has a kegg id (cpdModelSEEDid). If it did, then assign 
    # EX_cpdModelSEEDid_e0 as the ModelSEED id for the exchange reaction
    for rxn in [r for r in rxn_list if r.ModelSEED_id == [] and 'exchange' in r.type.lower()]:
        # Check whether the compound participating in this exchange reaction has a ModelSEED id
        if rxn.compounds[0].ModelSEED_id != []:
            exch_counter += 1
            if isinstance(rxn.compounds[0].ModelSEED_id,str):
                rxn.ModelSEED_id = 'EX_' + rxn.compounds[0].ModelSEED_id + '_e0'
            elif isinstance(rxn.compounds[0].ModelSEED_id,list):
                rxn.ModelSEED_id = []
                for sid in rxn.compounds[0].ModelSEED_id:
                    rxn.ModelSEED_id.append('EX_' + sid + '_e0')

    # Reactions with more than one ModelSEED id
    more_than_one_ModelSEED_id_rxns = [c for c in rxn_list if c.ModelSEED_id != [] and isinstance(c.ModelSEED_id,list)]
    if stdout_msgs:
        if len(more_than_one_ModelSEED_id_rxns) > 0 and warnings:
            print '\nWARNING! More than one ModelSEED id was found for %i reaction(s) including:'%(len(more_than_one_ModelSEED_id_rxns))
        for rxn in more_than_one_ModelSEED_id_rxns:
            print '\t',rxn.id,'\t',rxn.ModelSEED_id

        print '\nSummary of the search for the ModelSEED id of reactions:'
        print '\tTotal # of reactions with a ModelSEED id = %i (out of %i)'%(KEGG_counter + name_counter + id_counter + exch_counter,len(rxn_list))
        print '\t\t# of reactions matched by KEGG id = ',KEGG_counter
        print '\t\t# of reactions matched by name = ',name_counter
        print '\t\t# of reactions matached by model id = ',id_counter
        print '\t\t# of exchange reactions matached by their compound = ',exch_counter
        print '\tTotal # of reactions with no ModelSEED id = %i (out of %i)\n'%(len(rxn_list) - (KEGG_counter + name_counter + id_counter + exch_counter),len(rxn_list))

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
    get_cpd_ModelSEED_id(model.compounds, compart_list = [c.id for c in model.compartments],stdout_msgs = stdout_msgs)
    print '    Getting ModelSEED ids for reactions ...'
    get_rxn_ModelSEED_id(model.reactions, compart_list = [c.id for c in model.compartments],stdout_msgs = stdout_msgs)
 
