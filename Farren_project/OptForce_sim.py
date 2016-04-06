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
    model = create_model(model_organism = model_organism, model_info = {'id':'iJO1366', 'sbml_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, stdout_msgs = True, warnings = True)

    # Load reaction subsystems 
    load_source('dataFile', model_path + 'iJO1366_rxns_subsystems.py')
    import dataFile
    rxns_subsystems = dataFile.rxns_subsystems
    #for r_id in rxns_subsystems.keys():
    #    model.reactions_by_id[r_id].subsystem = rxns_subsystems[r_id]


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
    Ecoli_nsAA_producing = create_modifed_model(nsAA = nsAA, aeration = aeration, media_type = media_type, stdout_msgs = stdout_msgs, warnings = warnings)

    # Load blocked rxns list 
    blocked_results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_blocked_rxns_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'
    load_source('dataFile',blocked_results_filename)
    import dataFile
    blocked_rxns = dataFile.blocked_rxns

    results_filename = home_dir + 'work/Farren_project/results/iJO1366nsAA_fullyCoupled_sets_' + aeration + '_' + media_type + '_' + nsAA.lower() + '.py'

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
    read_inVivo_essential_rxns_fromThisFile = home_dir + 'work/models/Escherichia_coli/iJO1366/' + 'inVivo_essential_rxns_aerobic_glucose_minimal.py'

    MUST_doubles_params = {'objective_thr':1, 'validate_results':True, 'results_filename':results_filename_base + '_MUST_doubles.py'}

    FORCE_params = {'total_interven_num':10, 'notMUST_total_interven_num':0, 'notXrxns_interven_num':0, 'notLrxns_interven_num':0, 'notUrxns_interven_num':0, 'fixed_X_rxns':[], 'fixed_L_rxns':[], 'fixed_U_rxns':[], 'build_new_optModel':True, 'dual_formulation_type': 'simplified', 'stopWith_product_yield_percent': 95, 'validate_results': True, 'results_filename':results_filename_base + '_FORCE_sets.py'}

    OptForce_nsAA = OptForce(model = Ecoli_nsAA_producing, product_exchrxn_id = product_exchrxn_id, product_targetYield_percent = 80, min_biomass_percent = 10,  growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, read_exp_flux_bounds_ref_fromThisFile = read_exp_flux_bounds_ref_fromThisFile, read_blocked_rxns_fromThisFile = read_blocked_rxns_fromThisFile, read_inSilico_essential_rxns_fromThisFile = read_inSilico_essential_rxns_fromThisFile, read_inVivo_essential_rxns_fromThisFile = read_inVivo_essential_rxns_fromThisFile, MUST_doubles_params = MUST_doubles_params, FORCE_params = FORCE_params, stdout_msgs = True)

    #OptForce_nsAA.save_flux_bounds_ref_toThisFile = home_dir + 'work/Farren_project/results/optforce_' + aeration + '_' + media_type + '_flux_bounds_ref.py'
    OptForce_nsAA.read_flux_bounds_ref_fromThisFile = home_dir + 'work/Farren_project/results/optforce_' + aeration + '_' + media_type + '_flux_bounds_ref.py'

    OptForce_nsAA.save_flux_bounds_overprod_toThisFile = results_filename_base + '_flux_bounds_overprod.py'
    #OptForce_nsAA.read_flux_bounds_overprod_fromThisFile = results_filename_base + '_flux_bounds_overprod.py'

    OptForce_nsAA.save_MUST_singles_toThisFile = results_filename_base + '_MUST_singles.py'
    #OptForce_nsAA.read_MUST_singles_fromThisFile = results_filename_base + '_MUST_singles.py'

    #OptForce_nsAA.read_MUST_doubles_fromThisFile = results_filename_base + '_MUST_doubles.py'

    # Find the Must single sets
    OptForce_nsAA.run()

#------------------------------------
if __name__ == '__main__':
    OptForce_sim(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True)



