from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import userError
from copy import deepcopy

class gene(object):
    """
    A class holding the information about a gene 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 11-24-2014
    """

    def __init__(self, id, compartment = None, name = None, name_aliases = [], locus_pos = None, expression_level = None, notes = None): 
    
        # Gene id 
        self.id = id

        # Gene compartment (case insensitive string)
        self.compartment = compartment

        # Gene name (case insensitive string)
        self.name = name

        # Name name_aliases
        self.name_aliases = name_aliases

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

        self.__dict__[attr_name] = attr_value

   
