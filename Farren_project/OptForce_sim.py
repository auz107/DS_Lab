from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../')
sys.path.append('results/')
from copy import deepcopy
import itertools
from tools.userError import *
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.create_model import create_model
from tools.strain_design.OptForce.OptForce import OptForce
from ss_analysis import add_nsAA_pathways
from imp import load_source
from multiprocessing import Process, Manager, Value, Array
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

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

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

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
    
    # Set all objective function coefficients to zero
    for rxn in Ecoli_nsAA_producing.reactions:
        rxn.objective_coefficient = 0

    # Load experimental flux data
    load_source('dataFile','/usr2/postdoc/alizom/work/tools/strain_design/flux_data/Escherichia_coli/iJO1366_minimal_aerobic_glc_PMID_23036703.py')
    import dataFile
    exp_flux_bounds = dataFile.exp_flux_bounds

    # Create an instance of OptForce
    OptForce_nsAA = OptForce(model = Ecoli_nsAA_producing, product_exchrxn_id = product_exchrxn_id, exp_flux_bounds = exp_flux_bounds, product_targetYield_percent = 90, min_biomass_percent = 10,  growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, results_filename_base = 'results/optforce_' + aeration + '_' + media_type + '_' + nsAA.lower(), stdout_msgs = True, warnings = True)

    #OptForce_nsAA.flux_bounds_filename = '/usr2/postdoc/alizom/work/Farren_project/results/' + 'optforce_' + aeration + '_' + media_type + '_' + nsAA.lower() + '_flux_bounds.py'

    # Find the Must single sets
    OptForce_nsAA.run()

#------------------------------------
if __name__ == '__main__':
    OptForce_sim(nsAA = 'Pyrrolysine', aeration = 'aerobic', media_type = 'minimal', stdout_msgs = True, warnings = True)



