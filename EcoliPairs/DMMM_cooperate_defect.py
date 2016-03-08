from __future__ import division
import sys, time
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.userError import *
from tools.core.metabolite import metabolite
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from tools.fba.death_rate import death_rate

class DMMM_cooperate_defect(DMMM):
    """
    Performs dynamic multi-species metabolic modeling (DMMM) (PMID: 20668487)
    while solving a different fba problem in each time point. The fba problems
    simulate cooperate or defect.   

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: 02-13-2015
    """   

    def update_fba_model(self):
        """ 
        This is an abstract method that can be defined by subclasses. This is useful
        for cases one needs to solve a different fba model (i.e., different objective
        function or different constriants) at each time point

        The dynamic objective function coefficients are assigned to the model for each
        community member using an instance variable called strategy.
        It is a dictionary, where the keys are time points and values are either 'C'
        (cooperate)  or 'D' (defect)
        """ 
        for member in self.community_members:
            for rxn in member.reactions:
                rxn.objective_coefficient = 0
            member.biomass_reaction.objective_coefficient = 1
 
            if member.organism.strategy[self._t].upper() == 'D':
                # No need to produce the metabolite needed by the partner
                for exch_rxn in member.cooperative_rxns:
                    exch_rxn.flux_bounds[0] = 0

            elif member.organism.strategy[self._t].upper() == 'C':
                # Cooperation invovles the production of one unit of the
                # metabolite needed by the partner
                for exch_rxn in member.cooperative_rxns:
                    exch_rxn.flux_bounds[0] = 1

            else:
                raise userError('**ERROR! Invalid strategy (acceptable choices are C or D)')


