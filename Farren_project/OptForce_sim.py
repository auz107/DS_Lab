from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../')
sys.path.append('results/')
from copy import deepcopy
import itertools
from tools.userError import userError
from tools.globalVariables import *
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.create_model import create_model
from tools.fba.find_blocked_rxns import find_blocked_rxns
from tools.fba.find_essential_rxns import find_essential_rxns
from tools.fba.fcf import fcf
from tools.strain_design.OptForce.OptForce import OptForce
from ss_analysis import add_nsAA_pathways
from imp import load_source
from multiprocessing import Process, Manager, Value, Array
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

def create_modifed_model(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True):
    """
    Creates the updated model after adding external reactions for the production of nsAAs

    Ali R. Zomorrodi - Segre's Lab @ BU
    Last updated: 03-17-2016
    """
    # Model path
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'

    if media_type.lower() == 'minimal' and aeration.lower() == 'aerobic':
        flux_bounds_dict = {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}
        flux_bounds_filename = model_path + 'iJO1366_minimal_glucose_aerobic.py'
    elif media_type.lower() == 'minimal' and aeration.lower() == 'anaerobic':
        flux_bounds_dict = {'EX_glc(e)':[-100,1000]} 
        flux_bounds_filename = model_path + 'iJO1366_minimal_glucose_anaerobic.py'
    elif media_type.lower() == 'rich' and aeration.lower() == 'aerobic':
        flux_bounds_dict = {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]} 
        flux_bounds_filename = model_path + 'iJO1366_rich_glucose_aerobic.py' 
    elif media_type.lower() == 'rich' and aeration.lower() == 'anaerobic':
        flux_bounds_dict = {'EX_glc(e)':[-100,1000]} 
        flux_bounds_filename = model_path + 'iJO1366_rich_glucose_anaerobic.py' 
    else:
        raise ValueError('Unknown media_type and/or aeraiton type: {}, {}'.format(media_type,aeration))

    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    # Orignal iJo1266 model
    model = create_model(model_organism = model_organism, model_info = {'id':'iJO1366', 'input_file_type':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, stdout_msgs = True, warnings = True)

    # Load the list of transport-like (pseudo-transport) reacitons and add "pseudo-transport" to their subsystem 
    load_source('dataFile', model_path + 'iJO1366_transport_like_rxns.py')
    import dataFile
    transport_like_rxns = dataFile.transport_like_rxns
    for r_id in transport_like_rxns:
        model.reactions_by_id[r_id].subsystem = model.reactions_by_id[r_id].subsystem + ' (pseudo-transport)' 

    # Add the nsAA biosynthesis pathways to the model
    if nsAA.lower() == 'pyrrolysine':
        add_pyrrolysine = True
        add_pAF = False
        product_exchrxn_id = 'EX_pyrlys_L(e)'
    elif nsAA.lower() == 'paf':
        add_pyrrolysine = False
        add_pAF = True
        product_exchrxn_id = 'EX_paphe(e)'
    else:
        raise ValueError('Invalid value for nsAA')

    if stdout_msgs:
        print '\nproduct_exchrxn_id = {}'.format(product_exchrxn_id)

    Ecoli_nsAA_producing = add_nsAA_pathways(model = deepcopy(model), add_pyrrolysine = add_pyrrolysine, add_pAF = add_pAF, stdout_msgs = stdout_msgs, warnings = warnings) 

    # Find the max biomass flux for Ecoli_nsAA_producing 
    print '\n--- FBA after adding new pathways ---'
    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = flux_bounds_filename, flux_bounds = flux_bounds_dict)

    Ecoli_nsAA_producing.fba(stdout_msgs = False)
    if Ecoli_nsAA_producing.fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_biomass = Ecoli_nsAA_producing.fba_model.solution['objective_value']
        if stdout_msgs:
            if nsAA.lower() == 'pyrrolysine':
                print 'max biomass = {} , Pyrrolysine exch flux = {}'.format(max_biomass,Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'])
            elif nsAA.lower() == 'paf':
                print 'max biomass = {} , p-amino-phenylalanine exch. flux = {}'.format(max_biomass,Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
       
    else:
        raise userError('FBA to find max biomass flux for Ecoli_nsAA_producing was not solved to optimality')

    # Find the max pyrrolysine exchange reaction flux 
    if add_pyrrolysine:
        print '\n--- FBA to find max pyrrolysine exchange reaction flux ---'
        for rxn in Ecoli_nsAA_producing.reactions:
            rxn.objective_coefficient = 0
        Ecoli_nsAA_producing.reactions_by_id['EX_pyrlys_L(e)'].objective_coefficient = 1
        set_specific_bounds(model = Ecoli_nsAA_producing, file_name = flux_bounds_filename, flux_bounds = flux_bounds_dict)

        Ecoli_nsAA_producing.fba(stdout_msgs = False)
        if Ecoli_nsAA_producing.fba_model.solution['exit_flag'] == 'globallyOptimal':
            max_EX_pyrlys =  Ecoli_nsAA_producing.fba_model.solution['objective_value']
            if stdout_msgs:
                print 'biomass = {} , max Pyrrolysine exch flux = {}'.format(Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes'][Ecoli_nsAA_producing.biomass_reaction.id],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'])
        else:
            raise userError('FBA to find max exchange reaction flux for L-pyrrolysine was not solved to optimality')

    # Find the max p-amino-phenylalanine exchange reaction flux 
    if add_pAF:
        print '\n--- FBA to find max p-amino-phenylalanine exchange reaction flux ---\n'
        for rxn in Ecoli_nsAA_producing.reactions:
            rxn.objective_coefficient = 0
        Ecoli_nsAA_producing.reactions_by_id['EX_paphe(e)'].objective_coefficient = 1
        set_specific_bounds(model = Ecoli_nsAA_producing, file_name = flux_bounds_filename, flux_bounds = flux_bounds_dict)
    
        Ecoli_nsAA_producing.fba(stdout_msgs = stdout_msgs)
        if Ecoli_nsAA_producing.fba_model.solution['exit_flag'] == 'globallyOptimal':
            max_EX_paphe =  Ecoli_nsAA_producing.fba_model.solution['objective_value']
            if stdout_msgs:
                print 'biomass = {} , max p-amino-phenylalanine exch. flux = {}'.format(Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes'][Ecoli_nsAA_producing.biomass_reaction.id],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
        else:
            raise userError('FBA to find max exchange reaction flux for p-amino-phenylalanine was not solved to optimality')

    return (Ecoli_nsAA_producing, product_exchrxn_id, flux_bounds_dict, flux_bounds_filename)


def find_inSilico_essential_rxns(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True): 
    """
    Finds the set of blocked reactions under the examined condition in the model
    """
    Ecoli_nsAA_producing = create_modifed_model(nsAA = nsAA, aeration = aeration, media_type = media_type, stdout_msgs = stdout_msgs, warnings = warnings)

    # Find blocked reactions under the examined condition
    results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_essential_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    essential_rxns = find_essential_rxns(model = Ecoli_nsAA_producing, viability_thr = 1, results_filename = results_filename,  stdout_msgs = True, warnings = True)

def find_all_blocked_rxns(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True): 
    """
    Finds the set of blocked reactions under the examined condition in the model
    """
    Ecoli_nsAA_producing = create_modifed_model(nsAA = nsAA, aeration = aeration, media_type = media_type, stdout_msgs = stdout_msgs, warnings = warnings)

    # Find blocked reactions under the examined condition
    if stdout_msgs:
        print '--- Finding blocked reactions under the given uptake and aeration conditions ---'
    results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_blocked_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    blocked_rxns = find_blocked_rxns(model = Ecoli_nsAA_producing, always_blocked_only = False, results_filename = results_filename, stdout_msgs = True)

    # Find always blocked reactions 
    if stdout_msgs:
        print '--- Finding always blocked reactions ---'
    results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_always_blocked_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    always_blocked_rxns = find_blocked_rxns(model = Ecoli_nsAA_producing, always_blocked_only = True, results_filename = results_filename, stdout_msgs = True, warnings = True)

def find_coupled_rxns(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True): 
    """
    Finds the set of fully coupled reactions 
    """
    (Ecoli_nsAA_producing, product_exchrxn_id, flux_bounds_dict, flux_bounds_filename) = create_modifed_model(nsAA = nsAA, aeration = aeration, media_type = media_type, stdout_msgs = stdout_msgs, warnings = warnings)

    # Load blocked rxns list 
    blocked_results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_blocked_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    load_source('dataFile',blocked_results_filename)
    import dataFile
    blocked_rxns = dataFile.blocked_rxns

    results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_coupled_sets_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'

    fcf_inst = fcf(model = Ecoli_nsAA_producing, blocked_rxns = blocked_rxns, results_filename = results_filename, warnings = True, stdout_msgs = True)
    fcf_inst.run()


def OptForce_sim(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True):
    """
    Performs OptForce simulations for Ecoli strains producing pyrrolysine and pAF 

    INPUTS:
    -------
          nsAA: non-standard amino acid to consider (pyrrolysine or pAF)
    media_type: Type of growth mediam: M9 (minimal(, LB (rich)
      aeration: Type of aeration 'aerobic' or 'anaerobic'

    Ali R. Zomorrodi - Segre's Lab @ BU
    Last updated: 03-17-2016
    """
    if not isinstance(nsAA,str):
        raise TypeError('nsAA must be a string')
    elif nsAA.lower() not in ['pyrrolysine','paf']:
        raise ValueError('Invalid value for nsAA. Allowed choices are Pyrrolysine and pAF')

    if not isinstance(aeration,str):
        raise TypeError('aeration must a a string')
    elif aeration.lower() not in ['aerobic','anaerobic']:
        raise ValueError('Invalid aeration value: {}'.format(aeration))

    if not isinstance(media_type,str):
        raise TypeError('media_type must be a string')
    elif media_type.lower() not in ['minimal','rich']:
        raise ValueError('Invalid media_type value! Allowed choices are: [M9, LB]')

    print '\n------------- {} , {} ---------------\n'.format(media_type, aeration)
    Ecoli_nsAA_producing, product_exchrxn_id, flux_bounds_dict, flux_bounds_filename = create_modifed_model(nsAA = nsAA, aeration = aeration, media_type = media_type, stdout_msgs = stdout_msgs, warnings = warnings)
    
    # Set all objective function coefficients to zero
    for rxn in Ecoli_nsAA_producing.reactions:
        rxn.objective_coefficient = 0

    # Create an instance of OptForce
    results_filename_base = home_dir + 'work/Farren_project/results/optforce_' + aeration + '_' + media_type + '_' + nsAA.lower()

    # Data files 
    read_exp_flux_bounds_ref_fromThisFile = '/usr2/postdoc/alizom/work/tools/strain_design/flux_data/Escherichia_coli/iJO1366_minimal_aerobic_glc_PMID_23036703.py'
    read_blocked_rxns_fromThisFile = home_dir + 'work/Farren_project/results/iJO1366nsAA_blocked_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    read_inSilico_essential_rxns_fromThisFile = home_dir + 'work/Farren_project/results/iJO1366nsAA_essential_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    read_inVivo_essential_rxns_fromThisFile = home_dir + 'work/models/Escherichia_coli/iJO1366/iJO1366_inVivo_essential_rxns_aerobic_glucose_minimal.py'

    MUST_singles_diff_thr = 1
    MUST_doubles_params = {'objective_thr':1, 'validate_results':True, 'results_filename':results_filename_base + '_MUST_doubles.py', 'stdout_msgs': False}

    OptForce_nsAA = OptForce(model = Ecoli_nsAA_producing, product_exchrxn_id = product_exchrxn_id, product_targetYield_percent = 80, min_biomass_percent = 20, growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, read_exp_flux_bounds_ref_fromThisFile = read_exp_flux_bounds_ref_fromThisFile, read_blocked_rxns_fromThisFile = read_blocked_rxns_fromThisFile, read_inSilico_essential_rxns_fromThisFile = read_inSilico_essential_rxns_fromThisFile, read_inVivo_essential_rxns_fromThisFile = read_inVivo_essential_rxns_fromThisFile, MUST_singles_diff_thr = MUST_singles_diff_thr, MUST_doubles_params = MUST_doubles_params, stdout_msgs = True)

    if nsAA.lower() == 'pyrrolysine':
        #OptForce_nsAA.save_flux_bounds_ref_toThisFile = results_filename_base + '_flux_bounds_ref.py' 
        OptForce_nsAA.read_flux_bounds_ref_fromThisFile = results_filename_base + '_flux_bounds_ref.py' 

        #OptForce_nsAA.save_flux_bounds_overprod_toThisFile = results_filename_base + '_flux_bounds_overprod.py'
        OptForce_nsAA.read_flux_bounds_overprod_fromThisFile = results_filename_base + '_flux_bounds_overprod.py'

        #OptForce_nsAA.save_MUST_singles_toThisFile = results_filename_base + '_MUST_singles.py'
        OptForce_nsAA.read_MUST_singles_fromThisFile = results_filename_base + '_MUST_singles.py'

        OptForce_nsAA.read_MUST_doubles_fromThisFile = results_filename_base + '_MUST_doubles.py'

        # First round: PylB or PylC or PylD1 or PylD2 --> product yield = 31.96 (80.0% of theoretical maximum = 39.95) 
        #fixed_U_rxns = []
        fixed_U_rxns = ['DAPDC']
        #ignored_U_rxns = []
        ignored_U_rxns = ['PylB','PylC','PylD1','PylD2']
        OptForce_nsAA.FORCE_params = {'total_interven_num':10, 'notMUST_total_interven_num':5, 'notXrxns_interven_num':5, 'notLrxns_interven_num':5, 'notUrxns_interven_num':5, 'fixed_X_rxns':[], 'fixed_L_rxns':[], 'fixed_U_rxns':fixed_U_rxns, 'ignored_X_rxns':[], 'ignored_L_rxns':[], 'ignored_U_rxns':ignored_U_rxns, 'build_new_optModel':True, 'dual_formulation_type': 'simplified', 'stopWith_product_yield_percent': 95, 'validate_results': True, 'results_filename':results_filename_base + '_FORCE_sets.py'}

    elif nsAA.lower() == 'paf':
        #OptForce_nsAA.save_flux_bounds_ref_toThisFile = results_filename_base + '_flux_bounds_ref.py' 
        OptForce_nsAA.read_flux_bounds_ref_fromThisFile = results_filename_base + '_flux_bounds_ref.py' 

        #OptForce_nsAA.save_flux_bounds_overprod_toThisFile = results_filename_base + '_flux_bounds_overprod.py'
        OptForce_nsAA.read_flux_bounds_overprod_fromThisFile = results_filename_base + '_flux_bounds_overprod.py'

        #OptForce_nsAA.save_MUST_singles_toThisFile = results_filename_base + '_MUST_singles.py'
        OptForce_nsAA.read_MUST_singles_fromThisFile = results_filename_base + '_MUST_singles.py'

        OptForce_nsAA.read_MUST_doubles_fromThisFile = results_filename_base + '_MUST_doubles.py'

        # First round: PapB or PapC or PAPHES or ADCS --> yield = 44.21 (80.0% of theoretical maximum = 55.27)
        fixed_U_rxns = []
        #fixed_U_rxns = ['PapB']
        #ignored_U_rxns = []
        ignored_U_rxns = ['PapB','PapC','PAPHES','ADCS']
        OptForce_nsAA.FORCE_params = {'total_interven_num':10, 'notMUST_total_interven_num':5, 'notXrxns_interven_num':5, 'notLrxns_interven_num':5, 'notUrxns_interven_num':5, 'fixed_X_rxns':[], 'fixed_L_rxns':[], 'fixed_U_rxns':fixed_U_rxns, 'ignored_X_rxns':[], 'ignored_L_rxns':[], 'ignored_U_rxns':ignored_U_rxns, 'build_new_optModel':True, 'dual_formulation_type': 'simplified', 'stopWith_product_yield_percent': 95, 'validate_results': True, 'results_filename':results_filename_base + '_FORCE_sets.py'}

    else:
        raise userError('Invalid nsAA: {}'.format(nsAA))

    # Find the Must single sets
    OptForce_nsAA.run()

def print_must_subsystems():
    """
    Prints how many reaction participates in the subsystem for each type of manipualtion
    """
    # Model path
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'

    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    # Orignal iJo1266 model
    model = create_model(model_organism = model_organism, model_info = {'id':'iJO1366', 'input_file_type':'sbml','model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, perform_fba = False, stdout_msgs = True, warnings = True)

    Ecoli_pyr_producing = add_nsAA_pathways(model = deepcopy(model), add_pyrrolysine = True, add_pAF = False)
    Ecoli_paf_producing = add_nsAA_pathways(model = deepcopy(model), add_pyrrolysine = False, add_pAF = True)

    must_sets = {}
    must_sets['X_pyr'] = ['ICL', 'MALS', 'COBALT2tex', 'S2FE2ST']
    must_sets['L_pyr'] = ['PTAr']
    must_sets['U_pyr'] = ['PylB', 'PylC', 'PylD1', 'PylD2', 'PyrLYStex', 'ASPK', 'DAPDC', 'SDPTA', 'DHDPS', 'GLUSy', 'GLUDy', 'SDPDS', 'PPC', 'THDPS', 'DHDPRy', 'DAPE', 'ASAD', 'ASPTA', 'ACKr', 'ACtex', 'CO2tex', 'CO2tpp', 'NH4tex', 'NH4tpp']
    must_sets['all_pyr'] = list(set(must_sets['X_pyr'] + must_sets['L_pyr'] + must_sets['U_pyr']))

    must_sets['X_paf'] = ['ICL', 'MALS', 'COBALT2tex', 'S2FE2ST']
    must_sets['L_paf'] = ['ENO','GAPD', 'PGI', 'TPI', 'ACONTa', 'ACONTb', 'ICDHyr', 'AKGDH', 'PTAr']
    must_sets['U_paf'] = ['PapB', 'PapC', 'PAPHEtex', 'PAPHEt2rpp', 'ADCS', 'CHORS', 'DDPA', 'DHQS', 'DHQTi', 'DRPA', 'RPI', 'RPE', 'GLUDy', 'GND', 'PGK', 'PGM', 'PPM2', 'PSCVT', 'SHKK', 'SHK3Dr', 'TKT1', 'TKT2', 'ACKr', 'ACtex', 'CO2tex', 'CO2tpp', 'NH4tex', 'NH4tpp']
    must_sets['all_paf'] = list(set(must_sets['X_paf'] + must_sets['L_paf'] + must_sets['U_paf']))

    for must_set in ['X_paf', 'L_paf', 'U_paf', 'all_paf', 'X_pyr', 'L_pyr', 'U_pyr', 'all_pyr']: 
        print '\n----- {} --------'.format(must_set)
        subsys_dict = {}
        for r in must_sets[must_set]:
            if 'pyr' in must_set.lower():
                subsys = Ecoli_pyr_producing.reactions_by_id[r].subsystem
            elif 'paf' in must_set.lower():
                subsys = Ecoli_paf_producing.reactions_by_id[r].subsystem
            else:
                raise ValueError('Unknown must)set: {}'.format(must_set))

            if subsys != '':
                if subsys in subsys_dict.keys():
                    subsys_dict[subsys].append(r) 
                else:
                    subsys_dict[subsys] = [r]    
        for k in sorted(subsys_dict.keys(), key = lambda x: len(subsys_dict[x]), reverse = False):
            print '{}\t{}'.format(k,len(subsys_dict[k]))
  
#------------------------------------
if __name__ == '__main__':
    OptForce_sim(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True)



