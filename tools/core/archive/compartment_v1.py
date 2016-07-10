from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from copy import deepcopy

class compartment(object):
    """
    A class holding the information about a compartment (e.g., cytosol, mitochondria, etc) 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 11-24-2014
    """

    def __init__(self, id, name = None, synonyms = [], model = None, notes = None): 
    
        # Gene id 
        self.id = id

        # Gene name (case insensitive string)
        if name == None or name == '':
            self.name = self.id
        else:
            self.name = name

        # Name synonyms
        self.synonyms = synonyms

        # The model in which this compartment is used
        self.model = model

        # Notes and comments
        self.notes = notes
  
