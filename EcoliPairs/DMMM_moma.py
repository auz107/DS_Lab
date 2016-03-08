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
from tools.fba.moma import moma
from tools.fba.DMMM import DMMM
from tools.fba.death_rate import death_rate

class DMMM_moma(DMMM):
    """
    Performs dynamic multi-species metabolic modeling (DMMM) (PMID: 20668487)
    while solving moma instead of fba at each time point. 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: 01-11-2015
    """   

    def mu_uptake_export_calc(self):
        """ 
        Compute the specific growth rate of community members and the uptake
        and export rates of shared compounds 
        """
        if self._t == self._t0:
            create_model = True
        else:
            create_model = False

        if self.store_dynamic_fluxes == False:
            flux_key = None 
        else:
            flux_key = self._t

        for member in self.community_members:
            if create_model == True:
                member.moma_model = moma(model = member, distance_type = '1-norm',create_model = create_model, flux_key = flux_key, store_opt_fluxes = True, screen_output = 'off')
            else:
                # Update the flux bounds
                for j in member.reactions:
                    member.moma_model.pyomo_momaModel.v[j].lb = j.flux_bounds[0]
                    member.moma_model.pyomo_momaModel.v[j].ub = j.flux_bounds[1]

            member.moma_model.run()

            if member.moma_model.solution['exit_flag'] == 'globallyOptimal':
                if self.store_dynamic_fluxes == False:
                    member.organism.mu[self._t] = member.biomass_reaction.flux
                else:
                    member.organism.mu[self._t] = member.biomass_reaction.flux[self._t]

            else:
                if self.screen_output.lower() == 'on':
                    print '   Death rate ...'
                # Find the limiting nutrients that this member takes up 
                limiting_nutrients = [uptake_rxn.reactants[0] for shared_cmp in self.shared_compounds for uptake_rxn in shared_metab.reactant_reactions if uptake_rxn in member.reactions]
                mu_death = death_rate(member,limiting_nutrients,screen_output = 'off')
                if mu_death is not None:
                    member.organism.mu[self._t] = mu_death
                else:
                    print '\n**WARNING! mu_death could not be computed successfully. 1% of max biomass flux for the wild-type under the same condition was assigned as the death rate'
                    member.organism.mu[self._t] = -0.01*member.wildType_max_biomass

                # Set the flux of all reactions (including export reactions for shared compounds)
                # to zero
                for r in member.reactions:
                    r.flux[self._t] = 0


