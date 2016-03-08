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
from tools.ancillary.cell_pool_conc import cell_pool_conc

class DMMM_mortality(DMMM):
    """
    Performs dynamic multi-species metabolic modeling (DMMM) (PMID: 20668487)
    while accounting for the amino acid pool added to the pool of the shared
    compounds due to cell death 

    Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
    Last updated: 12-10-2015
    """   

    def update_concentrations_batch(self):
        """ 
        Updates cell concentrationns and the concentration of extracellular compounds
        for a batch reactor 
        """  
        #--- Update the cell concentrations ---
        # dX_i/dt = mu_i*(1-rmp/100)*X_i*(1 - sum(i,(1-rmp/100)X(i))/carrying_capacity) or 
        # (X_i(t+dt) - X_i(t))/dt = mu*(1-rmp/100)*X_i(t)*(1 - sum(i,(1-rmp/100)*X_i(t))/carrying_capacity)
        # where rmp is the random mortality percentage
        # If concentration is negative set it to zero
        members_gDW_per_ml_total = sum([(1 - self.random_mortality_percentage/100)*member.organism.gDW_per_ml[self._t] for member in self.community_members])
        self._logisticG_factor_gDW_per_ml = 1 - members_gDW_per_ml_total/self.carrying_capacity['gDW_per_ml']
        if len([member for member in self.community_members if member.organism.cells_per_ml != None]) == len(self.community_members):
            members_cells_per_ml_total = sum([(1 - self.random_mortality_percentage/100)*member.organism.cells_per_ml[self._t] for member in self.community_members])
            self._logisticG_factor_cells_per_ml = 1 - members_cells_per_ml_total/self.carrying_capacity['cells_per_ml']

        for member in self.community_members:
            # We always need gDW_per_ml to update compound concentrations but
            # providing cells_per_ml is optional
            member.organism.gDW_per_ml[self._t + self._dt] = max(member.organism.mu[self._t]*(1-self.random_mortality_percentage/100)*member.organism.gDW_per_ml[self._t]*self._logisticG_factor_gDW_per_ml*self._dt + member.organism.gDW_per_ml[self._t],0)

            if member.organism.cells_per_ml is not None:
                member.organism.cells_per_ml[self._t + self._dt] = max(member.organism.mu[self._t]*(1 - self.random_mortality_percentage/100)*member.organism.cells_per_ml[self._t]*self._logisticG_factor_cells_per_ml*self._dt + member.organism.cells_per_ml[self._t],0)

            # Total death rate (** newly added for DMMM_mortality **)
            if member.organism.mu[self._t] < 0:
                # In thise case random_mortality_rate has already been incorporated into mu
                # (see DMMM.py)
                member.organism.total_death_rate[self._t] = member.organism.mu[self._t] 
            else:
                member.organism.total_death_rate[self._t] = member.organism.random_mortality_rate


        #--- Update shared compound concentrations ---
        # dC/dt = f where, f = sum(k,v_export_k*X_k) - sum(k,v_uptake_k*X_k) + dead_pool_rate
        # (C(t+dt) - c(t))/dt = sum(k,v_export_k*X_k) - sum(k,v_uptake_k*X_k) + dead_pool 
        # where, dead_pool_rate = sum(k,-self.cell_pool_factor*cell_pool_concentration_k*total_death_rate_k*X_k)
        # Here, cell_pool_concentration is the concentration of the compound pool 
        # per gDW or per cells, which should have already been assigned to each 
        # shared compound. The minus sign is because total_death_rate is negative
        # while dead_pool_rate must be non-negative. Here, self.cell_pool_factor is
        # the factor that should be multiplied by concentration of that compound 
        # (this is because sometimes we want to explore how the growth is affected if 
        # the cell pool is higher than the ones reported experimentally)
        total_cmps_conc = sum([cmp.concentration[self._t] for cmp in self.shared_compounds])
        self._logisticG_factor_cmps = 1 - total_cmps_conc/self.carrying_capacity['compounds_mM']
        self._f = dict([(cmp,None) for cmp in self.shared_compounds])

        if not hasattr(self,'cell_pool_factor'):
            self.cell_pool_factor = 1

        for shared_cmp in self.shared_compounds:
            dead_pool_rate = -sum([self.cell_pool_factor*shared_cmp.cell_pool_concentration[member.organism.id]*member.organism.total_death_rate[self._t]*1000*member.organism.gDW_per_ml[self._t] for member in self.community_members])
            if dead_pool_rate < 0:
                raise userError('dead_pool_rate is negative')

            f = sum([r.flux[self._t]*1000*r.model.organism.gDW_per_ml[self._t] for r in shared_cmp.reactions]) + dead_pool_rate
            self._f[shared_cmp] = f

            conc = f*self._logisticG_factor_cmps*self._dt + shared_cmp.concentration[self._t]

            if conc >= 0 or (conc < 0 and abs(conc) <= 1e-9):
                conc = max(conc,0)

            shared_cmp.concentration[self._t + self._dt] = conc


