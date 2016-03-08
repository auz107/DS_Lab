from __future__ import division
import sys, time
sys.path.append('../../')
from tools.userError import userError 
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.death_rate import death_rate
from tools.fba.DMMM import DMMM

class DMMM_multipleCarbonSrc(DMMM):
    """
    Essentially the same as DMMM, except that (1) the upper bound on oxygen uptake flux is set
    to two times that for limiting carbon source uptake and (2) the specific grwoth rate for  
    each community member is set to the max biomass minus the max biomass values when no carbon
    source is taken up. This is because the kbase gap filled models can grow even if no 
    carbon source is supplied. The max biomass in the abscence of any carbon source in the 
    minimal medium should have already been computed and stored in a parameter called 
    noCarbonSrc_maxBiomass in the model object for each community member. 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: 07-23-2015
    """   

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
        
        # Set the LB on oxygen uptake two times that for the limiting carbon source. If there are more than 
        # one limiting resource set it to two times the minimum of that for all such limiting resources
        for member in self.community_members:
            member.oxygen_exchrxn.flux_bounds[0] = 2*min([exchrxn.flux_bounds[0] for exchrxn in member.limiting_carbonSrc_exchrxns])

    def mu_uptake_export_calc(self):
        """ 
        Compute the specific growth rate of community members and the uptake
        and export rates of shared compounds 
        """
        for member in self.community_members:
            member.fba(create_model = False, store_opt_fluxes = self.store_dynamic_fluxes, flux_key = self._t, screen_output = 'off')
            if member.fba_model.solution['exit_flag'] == 'globallyOptimal':
                member.organism.mu[self._t] = member.biomass_reaction.flux[self._t] + member.organism.mortality_rate - member.noCarbonSrc_maxBiomass
            else:
                if self.screen_output.lower() == 'on':
                    print '   Death rate ...'
                # Find the limiting nutrients that this member takes up 
                limiting_nutrients = [uptake_rxn.reactants[0] for shared_cmp in self.shared_compounds for uptake_rxn in shared_cmp.reactant_reactions if uptake_rxn in member.reactions]
                mu_death = death_rate(limiting_nutrients = limiting_nutrients,ms_biomassYield_calc = False, wildType_max_biomass = member.wildType_max_biomass, screen_output = 'off')
                member.organism.mu[self._t] = mu_death + member.organism.mortality_rate

                # Set the flux of all reactions (including export reactions for shared compounds)
                # to zero
                for r in member.reactions:
                    r.flux[self._t] = 0

