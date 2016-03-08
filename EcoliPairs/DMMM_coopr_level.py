from __future__ import division
import sys, time
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.userError import *
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from tools.fba.death_rate import death_rate

class DMMM_coopr_level(DMMM):
    """
    Performs DMMM for the cooperation level simulation 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: 01-14-2015 
    """   

    def update_coopr_exchrxn_bounds(self):
        """
        Udates the bounds on the flux of cooperative reactions
        """
        for member in self.community_members:

            for coopr_exchrxn in member.coopr_exchrxns:
                member.biomass_reaction.objective_coefficient = 0
                coopr_exchrxn.objective_coefficient = 1
                member.fba(build_new_optModel = False, reset_fluxes = False, store_opt_fluxes = False, flux_key = None, stdout_msgs = False)
                if member.fba_model.solution['exit_flag'] == 'globallyOptimal':
                    coopr_exchrxn.flux_bounds[0] = (member.cooperation_level/100)*member.fba_model.solution['objective_value']
                else:
                    coopr_exchrxn.flux_bounds[0] = 0 

                #print '\t\t',member.organism.id,':\tLB ',coopr_exchrxn.id,' = ',coopr_exchrxn.flux_bounds[0]
            # Reset the objective coefficients
            coopr_exchrxn.objective_coefficient = 0
            member.biomass_reaction.objective_coefficient = 1

    def uptake_rate_UB_calc(self):
        """ 
        Compute the upper bound on the uptake rates of the shared compounds 
        (LB on exchange fluxes) using kinetic expressions 
        """ 
        for shared_cmp in self.shared_compounds:
            for uptake_reaction in shared_cmp.reactant_reactions:         
                # Assume uptake reactions are exchange reactions with only one
                # reactant 
                uptake_reaction.reactants[0].concentration = shared_cmp.concentration[self._t] 
                uptake_reaction.kinetic_rate_calc(assignLB = True) 

        self.update_coopr_exchrxn_bounds()

