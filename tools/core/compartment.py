from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from copy import deepcopy

class compartment(object):
    """
    A class holding the information about a compartment (e.g., cytosol, mitochondria, etc) 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06-02-2016
    """

    def __init__(self, id, name = '', name_aliases = [], model_id = '', notes = ''): 
    
        # Gene id 
        self.id = id

        # na,e 
        if name == None or name == '':
            self.name = self.id
        else:
            self.name = name

        # name aliases
        self.name_aliases = name_aliases

        # The model in which this compartment is used
        self.model_id = model_id

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
            raise TypeError("Invalid 'id' for compartment " + str(attr_value) + "! 'id' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # Name
        if attr_name == 'name' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'name' for compartment " + self.id + "! 'name' must be a string. A " + str(type(attr_value)) + " type object was entered instead")
    
        # Model id
        if attr_name == 'model_id' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'model_id' for compartment " + self.id + "! 'model_id' must be a string. A " + str(type(attr_value)) + " type object was entered instead")
    
        self.__dict__[attr_name] = attr_value

