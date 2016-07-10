from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *

class organism(object):
    """
    A class holding the information about an organism

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 10-22-2014
    """

    def __init__(self, id, name = '', name_aliases = [], domain = '', genus = '', species = '', strain = '', model_id = None, ModelSEED_type = '', gDW_per_cell = None, gWW_per_cell = None, cells_per_ml = None, gDW_per_ml = None, mu = None, random_mortality_rate = None, notes = ''): 
    
        # Organism id 
        self.id = id

        # organism complete name (case insensitive string)
        self.name = name

        # Name aliases
        self.name_aliases = name_aliases

        # Domain (case insensitive string).Example of choices are:  
        # bacteria, archaea,Eukaryotes 
        self.domain = domain

        # Genus (case insensitive string). Example: Escherichia 
        self.genus = genus

        # Species ((case insensitive string)). Example: coli
        self.species = species

        # Strain ((case insensitive string)). Example: MG1655
        self.strain = strain

        # model id
        self.model_id = model_id

        # ModelSEED type (bacteria_GramPositive, bacteria_GramNegative, Human and Plant)
        self.ModelSEED_type = ModelSEED_type

        # Gram of dry weight per (one) cell
        self.gDW_per_cell = gDW_per_cell

        # Gram of wet weight per (one) cell
        self.gWW_per_cell = gWW_per_cell

        # Cell mass concnetration in number of cells per ml of the culture
        self.cells_per_ml = cells_per_ml

        # Cell mass concentration in gram of dry weight per ml of the culture
        self.gDW_per_ml = gDW_per_ml

        # Specific growth rate (in 1/h)
        self.mu = mu

        # Mortality rate
        if random_mortality_rate > 0:
            raise userError('Mortality rate for organism ' + self.id + ' must be non-positive')
        else:
            self.random_mortality_rate = random_mortality_rate

        # Notes and comments
        if isinstance(notes,str):
            self.notes = notes
        else:
            self.notes = ''

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
            raise TypeError("Invalid 'id' for organism " + str(attr_value) + "! 'id' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # Name
        if attr_name == 'name' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'name' for organism " + self.id + "! 'name' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # Name aliases
        if attr_name == 'name_aliases' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'name_aliases' for organism " + str(attr_value) + "! 'name_aliases' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # Name aliases
        if attr_name == 'name_aliases' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'name_aliases' for organism " + str(id) + "! 'name_aliases' must be a list of strings. A " + str(type(attr_value)) + " type object was entered instead")
        if attr_name == 'name_aliases' and len([n for n in attr_value if not isinstance(r,str)]) > 0:
            raise TypeError("Invalid 'name_aliases' for organism " + str(id) + "! 'name_aliases' must be a list of strings. Objects that are not string found in the list:" + str([n for n in attr_value if not isinstance(r,str)]))

        # Domain 
        if attr_name == 'domain' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'domain' for organism " + str(attr_value) + "! 'domain' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # ModelSEED type 
        if attr_name == 'ModelSEED_type' and attr_value != None and not isinstance(attr_value,str):
            raise TypeError("Invalid 'ModelSEED_type' for organism " + str(attr_value) + "! 'ModelSEED_type' must be a string. A " + str(type(attr_value)) + " type object was entered instead")
        elif attr_name == 'ModelSEED_type' and attr_value.lower() not in ['','bacteria_grampositive','bacteria_gramnegative','human','plant']:
            raise ValueError("Illegal value for 'ModelSEED_type' for organism " + str(attr_value) + "! Allowed values are bacteria_grampositive, bacteria_gramnegative, human and plant")

        # Random mortality rate
        if attr_name == 'random_mortality_rate' and attr_value > 0:
            raise ValueError('random_mortality_rate for organism ' + self.id + ' must be non-positive') 
   
        self.__dict__[attr_name] = attr_value

  
