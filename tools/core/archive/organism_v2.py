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

    def __init__(self, id, name = '', domain = '', genus = '', species = '', strain = '', gDW_per_cell = None, gWW_per_cell = None, cells_per_ml = None, gDW_per_ml = None, mu = None, mortality_rate = -1e-3, notes = ''): 
    
        # Organism id 
        self.id = id

        # organism complete name (case insensitive string)
        self.name = name

        # Domain (case insensitive string).Example of choices are:  
        # bacteria, archaea,Eukaryotes 
        self.domain = domain

        # Genus (case insensitive string). Example: Escherichia 
        self.genus = genus

        # Species ((case insensitive string)). Example: coli
        self.species = species

        # Strain ((case insensitive string)). Example: MG1655
        self.strain = strain

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
        if mortality_rate > 0:
            raise userError('Mortality rate for organism ' + self.id + ' must be non-positive')
        else:
            self.mortality_rate = mortality_rate

        # Notes and comments
        if isinstance(notes,str):
            self.notes = notes
        else:
            self.notes = ''

  
