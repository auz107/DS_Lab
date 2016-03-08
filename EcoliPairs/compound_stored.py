from __future__ import division
import sys
sys.path.append('../')
from tools.userError import *

class compound_stored(object):
    """
    A class holding the information about a shared compound.
    This is similar to the compound object but it does not have
    its methods in order to save only the required information 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 12-19-2014
    """

    def __init__(self, id, name = None, Kegg_id = None, ModelSeed_id = None, formula = None, compartment = None, reactions = [], reactant_reactions = [], product_reactions = [], concentration = None, deltaGf_range = []): 
    
        # Metabolite id or abbreviation in the model (string)
        self.id = id

        # Complete or expanded compound name (string)
        self.name = name

        # Name of the compartment for this compound (string) 
        self.compartment = compartment

        # Kegg id
        self.Kegg_id = Kegg_id

        # ModelSeed id
        self.ModelSeed_id = ModelSeed_id 

        # Chemical formula of this compound (string) 
        self.formula = formula

        # List of reaction objects in which this compound participates as a
        # reactant or product. 
        self.reactions = reactions

        # List of reaction objects in which this compound participates as a
        # reactant 
        self.reactant_reactions = reactant_reactions

        # List of reaction objects in which this compound participates as a
        # product 
        self.product_reactions = product_reactions

        # Concentration of this compound (real) 
        self.concentration = concentration

        # A tuple of the form (dGfmin,dGfmax) containing the min and max values of deltaG
        # (Gibbs free energy change) of formation for this reaction
        self.deltaGf_range = deltaGf_range 

