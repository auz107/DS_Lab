from __future__ import division
import sys
sys.path.append('../../')
import reaction
from compartment import compartment
from tools.userError import userError
from copy import deepcopy

class gene(object):
    """
    A class holding the information about a gene 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 11-24-2014
    """
    def __init__(self, id, name = '', name_aliases = [], compartment = [], reactions = [], model_id = '', locus_pos = None, expression_level = None, notes = None): 
    
        # Gene id 
        self.id = id

        # Gene name (case insensitive string)
        self.name = name

        # Name name_aliases
        self.name_aliases = name_aliases

        # Gene compartment (case insensitive string)
        self.compartment = compartment

        # Reactions coded for by this gene 
        self.reactions = reactions

        # model id
        self.model_id = model_id

        # A tuple of type (start,end) indicating the start and end locus position of the gene 
        self.locus_pos = locus_pos

        # Expression level of the gene 
        self.expression_level = expression_level
  
        # Notes and comments
        self.notes = notes

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # id 
        if attr_name == 'id' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'id' for gene " + str(attr_value) + "! 'id' must be a string. A " + str(id) + "! 'name_aliases' must be a list of strings. A " + str(type(attr_value)) + " type object was entered instead")

        # Name
        if attr_name == 'name' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'name' for gene " + self.id + "! 'name' must be a string. A " + str(id) + "! 'name_aliases' must be a list of strings. A " + str(type(attr_value)) + " type object was entered instead")

        # Name aliases
        if attr_name == 'name_aliases' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'name_aliases' for gene " + str(id) + "! 'name_aliases' must be a list of strings. A " + str(id) + "! 'name_aliases' must be a list of strings. A " + str(type(attr_value)) + " type object was entered instead")
        if attr_name == 'name_aliases' and len([n for n in attr_value if not isinstance(r,str)]) > 0:
            raise TypeError("Invalid 'name_aliases' for gene " + str(id) + "! 'name_aliases' must be a list of strings. Objects that are not string found in the list: " + str([n for n in attr_value if not isinstance(r,str)]))

        # compartment 
        if attr_name == 'compartments' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'compartments' format for reaction " + self.id + "! compartments for a model must be a list of objects of type compartment. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'compartments' and len([n for n in attr_value if not isinstance(n,compartment)]) > 0:
            raise TypeError("Invalid 'compartments' format for reaction " + self.id + "! compartments for a model must be a list of objects of type compartment. Objects that are not of type compartment found in the list: " + str([r for r in attr_value if not isinstance(n,compartment)]))

        if attr_name == 'reactions':
            self.set_reactions(reactions = attr_value)
        else:
            self.__dict__[attr_name] = attr_value

        # Model
        if attr_name == 'model_id' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'model_id' for gene " + self.id + "! 'model_id' must be a string. A " + str(id) + "! 'name_aliases' must be a list of strings. A " + str(type(attr_value)) + " type object was entered instead")

    def set_reactions(self,reactions):
        """
        Sets reactions attribute
        """
        if reactions is not None and not isinstance(reactions,tuple) and not isinstance(reactions,list):
            raise TypeError("Invalid 'reactions' for compound " + str(self.id) + "! 'reactions'  must be a list of objects of type reaction. A " + str(type(reactions)) + " type object was entered instead")
        if len([n for n in reactions if not isinstance(n,reaction.reaction)]) > 0:
            raise TypeError("Invalid 'reactions' for compound " + str(self.id) + "! 'reactions'  must be a list of objects of type 'reaction'.  Objects that are not of type reaction found in the list:" + str([n for n in reactions if not isinstance(n,reaction.reaction)]))

        problem_rxns = [r.id for r in reactions if self not in r.genes]
        if len(problem_rxns) > 0:
            raise userError('The folloiwng reactions appear in reactions of gene {} but they do not appear in reaction.genes: {}'.format(self.id, problem_rxns))

        if isinstance(reactions,tuple):
            self.__dict__['reactions'] = reactions
        elif isinstance(reactions,list):
            self.__dict__['reactions'] = tuple(reactions)

   
