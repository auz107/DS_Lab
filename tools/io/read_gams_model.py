from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.compartment import compartment
from tools.core.model import model
from models import *
from imp import load_source

def read_gams_model(file_name,model_name,model_organism, model_type, warnings = True, stdout_msgs = True): 
    """
    A class to read a gams model in python-formatted data
    and to convert it into a format we need. 
        

    INPUTS: 
    ------
          file_name: The name of the sbml file containing the model (String)
                     This string can contain the full path to the file  
         model_name: Name of the model (string) 
     model_organism: Can be either a string containing the name of the organism or it can be 
                     an instance of the object organism
         model_type: Type of the model (string, e.g., 'metabolic')
           warnings: Can be 'on' or 'off' shwoing whether the warnings should be writtten to the 
                     screen or not
        stdout_msgs: Can be 'on' or 'off' shwoing whether any messages should be written to the
                     screen

    OUTPUT:
               model: A model object

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 08-10-2015
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError("stdout_msgs must be True or False")
    if not isinstance(warnings,bool):
        raise TypeError("warnings must be True or False")

    # Import the data in module stored in file_name
    if type(file_name) == str:
        load_source('gams_model',file_name)
    import gams_model

    # organism
    if type(model_organism) is str:             # If name is provided
        model_organism = organism(id = model_organism)
    elif type(model_organism) is not organism:  # if an object of type organism is not provided
        raise userError('Model_organism must be either a string containing the name (id) of the organism or an instance of the class "organism"')

    # compartments
    compartments = []
    compartments_by_id = []
    comparts = {'c':'cytosol','p':'preplasm','e':'extracellular'}
    for compartment_id in ['c','e','p']: 
        comp = compartment(id = compartment_id, name = comparts[compartment_id])
        compartments.append(comp)
        compartments_by_id.append((compartment_id,comp))
    compartments_by_id = dict(compartments_by_id)

    # Compounds
    get_cmp_by_id = []
    compounds = []
    for cmp_id in sorted(gams_model.cmp_names):
        if '[c]' in cmp_id:
            compart = compartments_by_id['c']
        elif '[p]' in cmp_id:
            compart = compartments_by_id['p']
        elif '[e]' in cmp_id:
            compart = compartments_by_id['e']
        else:
            raise userError('Unknown compartment for metaabolite ' + cmp_id)

        m = compound(id = cmp_id,compartment = compart)
        compounds += [m]
        get_cmp_by_id += [(cmp_id,m)]

    get_cmp_by_id = dict(get_cmp_by_id)

    # reactions
    stoic_keys = gams_model.stoic_matrix.keys() # keys of the stoichiomteric matrix
    get_rxn_by_id = []
    reactions = []
    for rxn_id in sorted(gams_model.rxn_names):
        stoichiometry = []

        # Reaction reversibility
        if gams_model.rxn_types[rxn_id] == 0:
            reaction_rev = 'irreversible'
        elif gams_model.rxn_types[rxn_id] == 1 and 2 not in gams_model.rxn_types.values():
            reaction_rev = 'reversible'
        elif gams_model.rxn_types[rxn_id] == 1 and 2 in gams_model.rxn_types.values():
            reaction_rev = 'reversible_forward'
        elif gams_model.rxn_types[rxn_id] == 2: 
            reaction_rev = 'reversible_backward'
        elif gams_model.rxn_types[rxn_id] == 3 and 4 not in gams_model.rxn_types.values():
            reaction_rev = 'exchange'
        elif gams_model.rxn_types[rxn_id] == 3 and 4 in gams_model.rxn_types.values():
            reaction_rev = 'exchange_forward'
        elif gams_model.rxn_types[rxn_id] == 4:
            reaction_rev = 'exchange_backward'


        # Reaction stoichiometry
        for m_id in [m[0] for m in stoic_keys if m[1] == rxn_id]:
            stoichiometry.append((get_cmp_by_id[m_id],gams_model.stoic_matrix[(m_id,rxn_id)]))   

        if len(stoichiometry) == 0:
            raise userError("**Error! The field 'stoichiometry' was not assigned for reaction " + str(rxn_id))

        reactions.append(reaction(id = rxn_id, stoichiometry = dict(stoichiometry), reversibility = reaction_rev))
 
    # biomass reaction
    biomass_reaction = [r for r in reactions for b in gams_model.biomass_rxn_names if r.id == b]
    if len(biomass_reaction) == 0 and warnings:
        print 'WARNING! No biomass reactions found. Check the model manually\n'
        biomass_reaction = None
    elif len(biomass_reaction) > 1 and warnings:
        print 'WARNING! More than one biomass reactions found: ',[r.id for r in biomass_reaction],'\n'

    # model
    return model(id = model_name, type = model_type, organism = model_organism, reactions = reactions, compounds = compounds, compartments = compartments, biomass_reaction = biomass_reaction, warnings = warnings, stdout_msgs = stdout_msgs)

#------------ Sample run -------------
if __name__ == "__main__":
    pass
