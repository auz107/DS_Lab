from __future__ import division
import sys
sys.path.append('../')
from tools.globalVariables import *
from tools.fba.fba import fba
from tools.io.create_model import create_model
from tools.core.organism import organism
from tools.fba.set_specific_bounds import set_specific_bounds

# Max value for the uptake or secretion 
max_uptake_leakiness = 0.3*10

#-------------------- Sample implementation -----------------
def rankAAs_byGrowthCost(leakiness_level = max_uptake_leakiness, stdout_msgs = True, stdout_msgs_details = False):
    """
    Ranks AAs according to their growth cost, i.e., the reduction in growth upon forcing the model to produce one unit of each AA
    """
    # Import the list of exchange reactions for amino acids
    from AAexchrxns import AAexchrxns

    AAexchrxns = [ 'EX_ala_L(e)', 'EX_arg_L(e)', 'EX_asn_L(e)', 'EX_asp_L(e)', 'EX_cys_L(e)', 'EX_gln_L(e)', 'EX_glu_L(e)', 'EX_gly(e)', 'EX_his_L(e)', 'EX_ile_L(e)', 'EX_leu_L(e)', 'EX_lys_L(e)', 'EX_met_L(e)', 'EX_phe_L(e)', 'EX_pro_L(e)', 'EX_ser_L(e)', 'EX_thr_L(e)', 'EX_trp_L(e)', 'EX_tyr_L(e)', 'EX_val_L(e)']  

    exchrxns_AA_map = {'EX_gln_L(e)': 'gln', 'EX_glu_L(e)': 'glu', 'EX_asn_L(e)': 'asn', 'EX_asp_L(e)': 'asp', 'EX_thr_L(e)': 'thr', 'EX_ser_L(e)': 'ser', 'EX_arg_L(e)': 'arg', 'EX_ala_L(e)': 'ala', 'EX_gly(e)': 'gly', 'EX_val_L(e)': 'val', 'EX_lys_L(e)': 'lys', 'EX_leu_L(e)': 'leu', 'EX_phe_L(e)': 'phe', 'EX_met_L(e)': 'met', 'EX_ile_L(e)': 'ile', 'EX_cys_L(e)': 'cys', 'EX_tyr_L(e)': 'tyr', 'EX_his_L(e)': 'his', 'EX_pro_L(e)': 'pro', 'EX_trp_L(e)': 'trp'}

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds)

    # Max biomass for wild-type
    maxbm_WT = WT.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M']

    AAsCosts = []
    for AAexch in AAexchrxns:
        set_specific_bounds(model = WT, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)
        set_specific_bounds(model = WT, flux_bounds = {AAexch:[leakiness_level, None]}, reset_flux_bounds = False)

        WT.fba(stdout_msgs = False)
        if WT.fba_model.solution['exit_flag'] == 'globallyOptimal': 
            AAsCosts.append((exchrxns_AA_map[AAexch], maxbm_WT - WT.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M']))
            if stdout_msgs_details:
                print 'AA exch rxn = {} , AA cost = {}'.format(exchrxns_AA_map[AAexch], maxbm_WT - WT.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M'])
        else:
            AAsCosts.append((exchrxns_AA_map[AAexch], 100)) 
            if stdout_msgs_details:
                print 'AA exch rxn = {} , AA cost = {}'.format(exchrxns_AA_map[AAexch], None) 

    AAsCosts.sort(key=lambda t: t[1],reverse = False)
    if stdout_msgs:
        print '\nAmino acids growth costs = '
        for AAcost in AAsCosts:
            print '{}. {}'.format(AAsCosts.index(AAcost) + 1, AAcost) 

    return AAsCosts

def rankAAs_byGrowthGain(uptake_level = max_uptake_leakiness, stdout_msgs = True, stdout_msgs_details = False):
    """
    Ranks AAs according to how the growth rate improves by taking up one unit of each AA 
    """
    # Import the list of exchange reactions for amino acids
    from AAexchrxns import AAexchrxns

    AAexchrxns = [ 'EX_ala_L(e)', 'EX_arg_L(e)', 'EX_asn_L(e)', 'EX_asp_L(e)', 'EX_cys_L(e)', 'EX_gln_L(e)', 'EX_glu_L(e)', 'EX_gly(e)', 'EX_his_L(e)', 'EX_ile_L(e)', 'EX_leu_L(e)', 'EX_lys_L(e)', 'EX_met_L(e)', 'EX_phe_L(e)', 'EX_pro_L(e)', 'EX_ser_L(e)', 'EX_thr_L(e)', 'EX_trp_L(e)', 'EX_tyr_L(e)', 'EX_val_L(e)']  

    exchrxns_AA_map = {'EX_gln_L(e)': 'gln', 'EX_glu_L(e)': 'glu', 'EX_asn_L(e)': 'asn', 'EX_asp_L(e)': 'asp', 'EX_thr_L(e)': 'thr', 'EX_ser_L(e)': 'ser', 'EX_arg_L(e)': 'arg', 'EX_ala_L(e)': 'ala', 'EX_gly(e)': 'gly', 'EX_val_L(e)': 'val', 'EX_lys_L(e)': 'lys', 'EX_leu_L(e)': 'leu', 'EX_phe_L(e)': 'phe', 'EX_met_L(e)': 'met', 'EX_ile_L(e)': 'ile', 'EX_cys_L(e)': 'cys', 'EX_tyr_L(e)': 'tyr', 'EX_his_L(e)': 'his', 'EX_pro_L(e)': 'pro', 'EX_trp_L(e)': 'trp'}

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds)

    # Max biomass for wild-type
    maxbm_WT = WT.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M']

    #--- Load the list of rxns that must be off in each mutant ---
    AAsGains = []
    for AAexch in AAexchrxns:
        set_specific_bounds(model = WT, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)
        set_specific_bounds(model = WT, flux_bounds = {AAexch:[-uptake_level,None]}, reset_flux_bounds = False)

        WT.fba(stdout_msgs = False)
        if WT.fba_model.solution['exit_flag'] == 'globallyOptimal': 
            AAsGains.append((exchrxns_AA_map[AAexch], WT.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M'] - maxbm_WT))
            if stdout_msgs_details:
                print 'AA exch rxn = {} , AA cost = {}'.format(exchrxns_AA_map[AAexch], WT.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M'] - maxbm_WT)
        else:
            AAsGains.append((exchrxns_AA_map[AAexch], -100)) 
            if stdout_msgs_details:
                print 'AA exch rxn = {} , AA cost = {}'.format(exchrxns_AA_map[AAexch], None) 

    AAsGains.sort(key=lambda t: t[1],reverse = True)

    if stdout_msgs:
        print '\nAmino acids growth gains = '
        for AAgain in AAsGains:
            print '{}. {}'.format(AAsGains.index(AAgain) + 1, AAgain) 

    return AAsGains

def rank_mutants_byGrowthGain(uptake_level = max_uptake_leakiness, stdout_msgs = True, stdout_msgs_details = False):
    """
    Ranks mustants according to their growth rate by obtaining one unit ot the AA they need
    """
    # Import the list of exchange reactions for amino acids
    from AAexchrxns import AAexchrxns

    AAs_map = { 'ala_L':'ala', 'arg_L':'arg', 'asn_L':'asn', 'asp_L':'asp', 'cys_L':'cys', 'gln_L':'gln', 'glu_L':'glu', 'gly':'gly', 'his_L':'his', 'ile_L':'ile', 'leu_L':'leu', 'lys_L':'lys', 'met_L':'met', 'phe_L':'phe', 'pro_L':'pro', 'ser_L':'ser', 'thr_L':'thr', 'trp_L':'trp', 'tyr_L':'tyr', 'val_L':'val'}  

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    model = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds)

    # Max biomass for wild-type
    maxbm_WT = model.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M']

    #--- Load the list of rxns that must be off in each mutant ---
    from mutants_rxn_info_iJO1366 import mutants_rxn_info_AAs, genes_AA_map
    from mutants_auxotrophy_AAs import mutants_auxotrophy

    # Mutants that are equivalent, i.e., they block the biosynthesis of the same amino acid. We consider
    # only from each group. These lists are extracted according to the contents of auxoMetabs.py 
    equivalentMutants = {}
    equivalentMutants['hisB'] = ['hisC', 'hisD', 'hisI']
    equivalentMutants['leuB'] = ['leuC', 'leuD']
    equivalentMutants['trpA'] = ['trpB']
    # Mutants that should not be considered
    not_consider_equivalent_mutatns = [m for k in equivalentMutants.keys() for m in equivalentMutants[k]]

    mutants_gains = []

    for gen in [g for g in mutants_rxn_info_AAs.keys() if g not in not_consider_equivalent_mutatns]:
        if stdout_msgs:
            print gen,':\t',[exchrxn for exchrxn in mutants_auxotrophy[gen]]

        # Growth condition
        set_specific_bounds(model = model, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)

        # Knockouts 
        set_specific_bounds(model = model, flux_bounds = dict([(rxn_id,[0,0]) for rxn_id in mutants_rxn_info_AAs[gen]]), reset_flux_bounds = False)   
        set_specific_bounds(model = model, flux_bounds = dict([(exchrxn,[-uptake_level, None]) for exchrxn_list in mutants_auxotrophy[gen] for exchrxn in exchrxn_list]), reset_flux_bounds = False)

        model.fba(stdout_msgs = False)
        if model.fba_model.solution['exit_flag'] == 'globallyOptimal': 
            mutants_gains.append((AAs_map[genes_AA_map[gen]], model.fba_model.solution['opt_rxnFluxes']['Ec_biomass_iJO1366_core_53p95M']))
        else:
            mutants_gains.append((AAs_map[genes_AA_map[gen]],-100)) 

    mutants_gains.sort(key=lambda t: t[1],reverse = True)

    if stdout_msgs:
        print '\nMutants growth gains = '
        for AAgain in mutants_gains:
            print '{}. {}   {}'.format(mutants_gains.index(AAgain) + 1, AAgain, AAgain[1] > maxbm_WT) 

    return mutants_gains

if __name__ == "__main__":
    AAsCosts = rankAAs_byGrowthCost(stdout_msgs = False)
    AAsGains = rankAAs_byGrowthGain(stdout_msgs = False)
    mutants_gains = rank_mutants_byGrowthGain(stdout_msgs = False)

    AAsCosts = dict(AAsCosts)
    AAsGains = dict(AAsGains)
    mutants_gains = dict(mutants_gains)
    print '\nAA growth cost for WT vs. AA growth gain for WT:'
    for AA in AAsCosts.keys():
        print '{}\t{}\t{}'.format(AA,AAsCosts[AA], AAsGains[AA])

    print '\nAA growth cost for WT vs. AA growth gain for mutants:'
    for AA in AAsCosts.keys():
        print '{}\t{}\t{}'.format(AA,AAsCosts[AA], mutants_gains[AA])

    print '\nCost (to producer) to benefit (to the mutant)  ratio:'
    CostToBenefit_ratio = sorted([(AA, AAsCosts[AA]+mutants_gains[AA]) for AA in AAsCosts.keys()], key=lambda t: t[1])
    for AAratio in CostToBenefit_ratio:
        print '{}. {}'.format(CostToBenefit_ratio.index(AAratio) + 1, AAratio) 
    for AAratio in CostToBenefit_ratio:
        print "'",AAratio[0],"_L',", 



