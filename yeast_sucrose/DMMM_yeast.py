from __future__ import division
import sys, time
sys.path.append('../')
import numpy as np
from tools.fba.DMMM import DMMM
import copy
from death_rate import death_rate
from tools.fba.fba import fba

class DMMM_yeast(DMMM):
    """
    Same as DMMM with funcitons mu_uptake_export_calc and  dilute redefined 

    mu_uptake_export_calc: A new segment was added that examines whether cooperation is ever possible by cooperators when the FBA
                           problem to find mu becomes infeasible
                   dilute: Function was modified to replicated the serial dilutions in Jeff Gore's paper where instead of providing a 
                           diluiton factor is just says that "Serial dilutions were performed daily (23 h of growth) such that the starting 
                           optical density was 0.0025, corresponding to ~150,000 cells" 
    """   
    def mu_uptake_export_calc(self):
        """ 
        Compute the specific growth rate of community members and the uptake
        and export rates of shared compounds 
        """
        #-- Compute mu and uptake and export rate for any member who has not driven to extinsion --
        for member in [m for m in self.community_members if (m.organism.gDW_per_ml != None and m.organism.gDW_per_ml[self._t] > 0) or (m.organism.cells_per_ml != None and m.organism.cells_per_ml[self._t] > 0)]:

            # Lag phase
            if self._t < self.lag_phase_time:
                member.organism.mu[self._t] = member.organism.random_mortality_rate

                # Set the flux of exchange reactions for this member related to all shared compounds to zero
                for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]:
                        rxn.flux[self._t] = 0

            # After lag phae
            else:
                if self._t == self._t0:
                    member.fba(build_new_optModel = True, reset_fluxes = True, store_opt_fluxes= True, flux_key = self._t, stdout_msgs = False)
                else:
                    member.fba(build_new_optModel = False, reset_fluxes = False, store_opt_fluxes = True, flux_key = self._t, stdout_msgs = False)

                if member.fba_model.solution['exit_flag'] == 'globallyOptimal':
                    member.organism.mu[self._t] = member.biomass_reaction.flux[self._t] + member.organism.random_mortality_rate
                else:
                    if self.stdout_msgs:
                        print '   Death rate for {} ...'.format(member.id)
                    # Find the limiting nutrients that this member takes up 
                    limiting_nutrients = [uptake_rxn.reactants[0] for shared_cmp in self.shared_compounds for uptake_rxn in shared_cmp.reactant_reactions if uptake_rxn in member.reactions]
                    if hasattr(member,'wildType_max_biomass'):
                        mu_death = death_rate(limiting_nutrients = limiting_nutrients,ms_biomassYield_calc = False, wildType_max_biomass = member.wildType_max_biomass, warnings = False, stdout_msgs = False)
                    else:
                        mu_death = death_rate(limiting_nutrients = limiting_nutrients,ms_biomassYield_calc = False, wildType_max_biomass = None, warnings = False, stdout_msgs = False)
                    member.organism.mu[self._t] = mu_death + member.organism.random_mortality_rate

                    # Set the flux of exchange reactions for this member related to all shared compounds to zero
                    for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]:
                            rxn.flux[self._t] = 0
    
                    # If this is the cooperator we can still find out whether EX_glc_secreted and EX_fru_secreted 
                    # can carry a non-zero flux upon relaxing the constraint on ATPM (i.e., to determine whether the 
                    # cooperator strain is able to cooperate at all). 
                    if member.id == 'Cooperator':
                        if self.stdout_msgs:
                            print '\tChecking whether SUCRe can carry any flux if ATPM constraint is relaxed...', 
    
                        if not hasattr(self,'max_SUCRe_fbaModel'):
                            # Find reactions participating in the objective function
                            obj_rxns = [r for r in member.reactions if r.objective_coefficient != 0]
      
                            # Save the objective coefficients for these reactions in a dictionary
                            obj_coeff = dict([(r,r.objective_coefficient) for r in obj_rxns])
    
                            # Set the objective coefficient for all reactions to zero
                            for r in member.reactions:
                                r.objective_coefficient = 0
    
                            # A parameter showing that changes were made to the objective_coefficient and ATPM reaction
                            changes_made = True
    
                            member.reactions_by_id['SUCRe'].objective_coefficient = 1
    
                            atpm_rxn_bounds = member.reactions_by_id['ATPM'].flux_bounds
                            member.reactions_by_id['ATPM'].flux_bounds = [0,1000]
    
                            self.max_SUCRe_fbaModel = fba(model = member, build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False)
    
                        # If the fba model already created just update the flux bounds for the reactions of shared compounds
                        else:
                            changes_made = False
                            for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]:
                                self.max_SUCRe_fbaModel.optModel.v[rxn.id].setlb(rxn.flux_bounds[0])
                                self.max_SUCRe_fbaModel.optModel.v[rxn.id].setub(rxn.flux_bounds[1])
    
                        SUCRe_fba_solution = self.max_SUCRe_fbaModel.run()
                        if SUCRe_fba_solution['exit_flag'] == 'globallyOptimal':
                             member.reactions_by_id['EX_glc_secreted'].flux[self._t] = SUCRe_fba_solution['opt_rxnFluxes']['EX_glc_secreted']
                             member.reactions_by_id['EX_fru_secreted'].flux[self._t] = SUCRe_fba_solution['opt_rxnFluxes']['EX_fru_secreted']
                             if self.stdout_msgs:
                                 print '\t\tSUCRe = {}  ,  EX_glc_secreted = {}   , EX_fru_secreted = {}'.format(SUCRe_fba_solution['opt_rxnFluxes']['SUCRe'], SUCRe_fba_solution['opt_rxnFluxes']['EX_glc_secreted'],SUCRe_fba_solution['opt_rxnFluxes']['EX_fru_secreted'])
                        else:
                             if self.stdout_msgs:
                                 print ' Infeasible! SUCRe cannot carry any flux even if ATPM constraint is relaxed' 
    
                        if changes_made:
                            # Set the objective coefficients back to what they were before
                            member.reactions_by_id['SUCRe'].objective_coefficient = 0
                            for r in obj_rxns:
                                r.objective_coefficient = obj_coeff[r]
                            member.reactions_by_id['ATPM'].flux_bounds = atpm_rxn_bounds

        #-- Set mu and uptake and export fluxes to zero for any member that has been extincted --
        for member in [m for m in self.community_members if (m.organism.gDW_per_ml != None and m.organism.gDW_per_ml[self._t] == 0) or (m.organism.cells_per_ml != None and m.organism.cells_per_ml[self._t] == 0)]:
            member.organism.mu[self._t] = 0
            for rxn in [r for c in self.shared_compounds for r in c.reactions if r in member.reactions]:
                    rxn.flux[self._t] = 0

    def dilute(self):
        """ 
        Performs dilution (used when reactor_type is serial_dilution 
        
        According to Jeff Gore's paper: "Serial dilutions were performed daily (23 h of growth) such that the starting 
        optical density was 0.0025, corresponding to ~150,000 cells"
        """ 
        # Total cell concentration
        total_cell_conc = sum([member.organism.cells_per_ml[self._t] for member in self.community_members])

        # Find the fraciton of each cell
        member_fracs = dict([(member,member.organism.cells_per_ml[self._t]/total_cell_conc) for member in self.community_members])

        for member in self.community_members:
            member.organism.cells_per_ml[self._t] = int(member_fracs[member]*self.serial_dilution_params['total_cell_conc_beginCycle']) 
            member.organism.gDW_per_ml[self._t] = member.organism.cells_per_ml[self._t]*member.organism.gDW_per_cell

        #--- shared metabolite concentrations are set to those at the initial time point ---
        for shared_cmp in self.shared_compounds:
            shared_cmp.concentration[self._t] = shared_cmp.concentration[self._t0]

