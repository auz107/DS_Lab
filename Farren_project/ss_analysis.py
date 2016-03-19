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
from multiprocessing import Process, Manager, Value, Array
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

"""
Thsi module contains the following functions:
        add_nsAA_pathways: Adds the pathways for non-standard amino acids to the model
            findTradeOffs: Performs various FBA simulations to identify trade-off between biomass  and 
                           nsAA production
              master_func: Creates the model for E. coli harboring the pathways producing the non-standard 
                           amino acids and assesses its impact on the cell growth 
      run_master_Tradeoff: Runs master func to find out trade-offs under different aeration and media types
    plot_tradeoff_results: Analyzes and plots the results 
run_plot_tradeoff_results: Plots the mortality results

Ali R. Zomorrodi - Segre's Lab @ BU
Last updated: 03-16-2016
"""

def add_nsAA_pathways(model, add_pyrrolysine = True, add_pAF = True, warnings = True, stdout_msgs = True):
    """
    Adds the pathways for non-standard amino acids to the model

    INPUTS:
    -------
    add_pyrrolysine: If True, reactions and compounds related to pyrrolysine biosynthesis are 
                     added to the model
            add_pAF: If True, reactions and compounds related to p-amino-phenylalanine
                     biosynthesis are added to the model
                 
    """
    if not isinstance(add_pyrrolysine,bool):
        raise TypeError('add_pyrrolysine must be either True or False')
    if not isinstance(add_pAF,bool):
        raise TypeError('add_pAF must be either True or False')

    model.id = model.id + '_nsAA_producing'

    cytosol = model.compartments_by_id['c']
    periplasm = model.compartments_by_id['p']
    extracellular = model.compartments_by_id['e']
    
    # Some required compounds
    atp_c = model.compounds_by_id['atp_c']
    adp_c = model.compounds_by_id['adp_c']
    h_c = model.compounds_by_id['h_c']
    h_p = model.compounds_by_id['h_p']
    pi_c = model.compounds_by_id['pi_c']
    h2o_c = model.compounds_by_id['h2o_c']
    nh4_c = model.compounds_by_id['nh4_c']
    nad_c = model.compounds_by_id['nad_c']
    nadh_c = model.compounds_by_id['nadh_c']

    #------ Pyrrolysine biosynthesis pathway -----
    if add_pyrrolysine:    
        #-- New compounds --
        metorn_c = compound(id = 'metorn_c', name = '(3R)-3-Methyl-D-ornithine', compartment = cytosol, ModelSEED_id = 'cpd24167', KEGG_id = 'C20277')
        metornlys_c = compound(id = 'metornlys_c', name = '(3R)-3-Methyl-D-ornithyl-Nepsilon-L-lysine', compartment = cytosol, ModelSEED_id = 'cpd24174')
        metglusaldlys_c = compound(id = 'metglusaldlys_c', name = '(3R)-3-methyl-D-glutamyl-semialdehyde-Nepsilon-L-Lysine', compartment = cytosol, ModelSEED_id = 'cpd24177')
        pyrlys_L_c = compound(id = 'pyrlys_L_c', name = 'L-pyrrolysine', compartment = cytosol, ModelSEED_id = 'cpd14859')
        pyrlys_L_p = compound(id = 'pyrlys_L_p', name = 'L-pyrrolysine', compartment = periplasm, ModelSEED_id = 'cpd14859')
        pyrlys_L_e = compound(id = 'pyrlys_L_e', name = 'L-pyrrolysine', compartment = extracellular, ModelSEED_id = 'cpd14859')
    
        model.add_compounds([metorn_c,metornlys_c,metglusaldlys_c,pyrlys_L_c,pyrlys_L_p,pyrlys_L_e])
    
        #-- New reactions --
        lys_L_c = model.compounds_by_id['lys_L_c']
    
        # PylB: L-lysine --> (3R)-3-methyl-D-ornithine
        PylB = reaction(id = 'PylB', name = 'methylornithine synthase', ModelSEED_id = 'rxn19159', stoichiometry = {lys_L_c:-1, metorn_c:1}, reversibility =  'irreversible',objective_coefficient = 0)
    
        # PylC: L-lysine  + (3R)-3-methyl-D-ornithine + ATP --> (3R)-3-Methyl-D-ornithyl-Neps-L-lysine + ADP + Pi + H+
        PylC = reaction(id = 'PylC', name = '(2R,3R)-3-methylornithyl-N6-lysine synthase', ModelSEED_id = 'rxn23366', stoichiometry = {lys_L_c:-1, metorn_c:-1, atp_c:-1, metornlys_c:1, adp_c:1, pi_c:1, h_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
    
        # PylD1: (3R)-3-Methyl-D-ornithyl-Neps-L-lysine + NAD+ + H2O --> (3R)-3-methyl-D-glutamyl-semialdehyde-Neps-L-Lysine + NADH +  NH3 + 2 H+ 
        # Given that ammonium in the model is in the form of NH4+, instead of NH3 + 2 H+, we write it as NH4+ + H+
        PylD1 = reaction(id = 'PylD1', name = 'pyrrolysine synthase', ModelSEED_id = 'rxn23370', KEGG_id = 'R10012', stoichiometry = {metornlys_c:-1, nad_c:-1, h2o_c:-1, metglusaldlys_c:1, nadh_c:1, nh4_c:1, h_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
    
        # PylD2: (3R)-3-methyl-D-glutamyl-semialdehyde-Neps-L-Lysine --> L-pyrrolysine + H2O + H+
        PylD2 = reaction(id = 'PylD2', name = 'pyrrolysine synthase', ModelSEED_id = 'rxn23371', stoichiometry = {metglusaldlys_c:-1, pyrlys_L_c:1, h2o_c:1, h_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
    
        # PyrLYSabcpp: atp_c + h2o_c + pyrlys_L_p --> adp_c + h_c + pyrlys_L_c + pi_c
        PyrLYSabcpp = reaction(id = 'PyrLYSabcpp', name = 'L-pyrrolysine transport via ABC system (periplasm)', stoichiometry = {atp_c:-1, h2o_c:-1, pyrlys_L_p:-1, adp_c:1, h_c:1, pyrlys_L_c:1, pi_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
    
        # PyrLYSt2pp: h_p + pyrlys_L_p --> h_c + pyrlys_L_c
        PyrLYSt2pp = reaction(id = 'PyrLYSt2pp', name = 'L-pyrrolysine transport in via proton symport (periplasm)', stoichiometry = {h_p:-1, pyrlys_L_p:-1, h_c:1, pyrlys_L_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
     
        # PyrLYSt3pp: h_p + pyrlys_L_c --> h_c + pyrlys_L_p
        PyrLYSt3pp = reaction(id = 'PyrLYSt3pp', name = 'L-pyrrolysine transport out via proton antiport (cytoplasm to periplasm)', stoichiometry = {h_p:-1, pyrlys_L_c:-1, h_c:1, pyrlys_L_p:1}, reversibility =  'irreversible', objective_coefficient = 0)
       
        # PyrLYStex:  pyrlys_L_e <==> pyrlys_L_p
        PyrLYStex = reaction(id = 'PyrLYStex', name = 'L-pyrrolysine transport via diffusion (extracellular to periplasm)', stoichiometry = {pyrlys_L_e:-1, pyrlys_L_p:1}, reversibility =  'reversible', objective_coefficient = 0)
    
        # EX_pyrlys_L_e: pyrlys_L_e <==>
        EX_pyrlys_L_e = reaction(id = 'EX_pyrlys_L(e)', name = 'L-pyrrolysine exchange', stoichiometry = {pyrlys_L_e:-1}, reversibility =  'exchange', objective_coefficient = 0)
    
        model.add_reactions([PylB,PylC,PylD1,PylD2,PyrLYSabcpp,PyrLYSt2pp,PyrLYSt3pp,PyrLYStex,EX_pyrlys_L_e])

        if stdout_msgs:
            print '\n---- Pyrrolysine biosynthesis reactions ---\n'
            for rxn in [PylB,PylC,PylD1,PylD2,PyrLYSabcpp,PyrLYSt2pp,PyrLYSt3pp,PyrLYStex,EX_pyrlys_L_e]:
                print rxn.id,':\t',rxn.get_equation()   
            print '\n'
     
    #----- P-amino-phenylalanine (pAF) biosynthesis pathway -----
    if add_pAF:
        #-- New compounds --
        adpphn_c = compound(id = 'adpphn_c', name = '4-amino-4-deoxyprephenate', compartment = cytosol)   
        aphpyr_c = compound(id = 'aphpyr_c', name = 'p-aminophenylpyruvate', compartment = cytosol)
        paphe_c = compound(id = 'paphe_c', name = 'p-amino-phenylalanine', compartment = cytosol)
        paphe_p = compound(id = 'paphe_p', name = 'p-amino-phenylalanine', compartment = periplasm)
        paphe_e = compound(id = 'paphe_e', name = 'p-amino-phenylalanine', compartment = extracellular)
    
        model.add_compounds([adpphn_c,aphpyr_c,paphe_c,paphe_p,paphe_e])
    
        #-- New reactions --
        co2_c = model.compounds_by_id['co2_c']   
        # L-glutamine  
        gln_L_c = model.compounds_by_id['gln_L_c']  
        # 2-Oxoglutarate
        akg_c = model.compounds_by_id['akg_c']       
        # 4-amino-4-deoxychorismate 
        cpd_4adcho_c = model.compounds_by_id['4adcho_c'] 
    
        # PapB: 4-amino-4-deoxychorismate --> 4-amino-4-deoxyprephenate
        PapB = reaction(id = 'PapB', name = 'chorismate mutase', stoichiometry = {cpd_4adcho_c:-1, adpphn_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
    
        # PapC: 4-amino-4-deoxyprephenate  + NAD+ --> p-aminophenylpyruvate + CO2 + NADH
        PapC = reaction(id = 'PapC', name = 'prephenate dehydrogenase', stoichiometry = {adpphn_c:-1, nad_c:-1, aphpyr_c:1, nadh_c:1, co2_c:1}, reversibility =  'irreversible', objective_coefficient = 0)
    
        # PAPHES: p-aminophenylpyruvate + L-glutamine <==>  p-amino-phenylalanine + 2-Oxoglutarate
        PAPHES = reaction(id = 'PAPHES', name = 'p-amino-phenylalanine synthase', stoichiometry = {aphpyr_c:-1, gln_L_c:-1, paphe_c:1, akg_c:1}, reversibility =  'reversible', objective_coefficient = 0)
     
        # PAPHEt2rpp: h_p + paphe_p <==> h_c + paphe_c 
        PAPHEt2rpp = reaction(id = 'PAPHEt2rpp', name = 'p-amino-phenylalanine reversible transport via proton symport (periplasm)', stoichiometry = {h_p:-1, paphe_p:-1, h_c:1, paphe_c:1}, reversibility =  'reversible', objective_coefficient = 0)
    
        # PAPHEtex: paphe_e <==> paphe_p 
        PAPHEtex = reaction(id = 'PAPHEtex', name = 'p-amino-phenylalanine transport via diffusion (extracellular to periplasm)', stoichiometry = {paphe_e:-1, paphe_p:1}, reversibility =  'reversible', objective_coefficient = 0)
    
        # EX_paphe(e): paphe_e <==>
        EX_paphe_e = reaction(id = 'EX_paphe(e)', name = 'p-amino-phenylalanine exchange', stoichiometry = {paphe_e:-1}, reversibility =  'exchange', objective_coefficient = 0)
    
        model.add_reactions([PapB,PapC,PAPHES,PAPHEt2rpp,PAPHEtex,EX_paphe_e])

        if stdout_msgs:
            print '\n---- p-amino-phenylalanine biosynthesis reactions ---\n'
            for rxn in [PapB,PapC,PAPHES,PAPHEt2rpp,PAPHEtex,EX_paphe_e]: 
                print rxn.id,':\t',rxn.get_equation()   
            print '\n'
    
    model.validate()

    return model

def findTradeOffs(input_data):
    """
    Performs various FBA simulations to identify trade-off between biomass  and nsAA production
    """
    model_path = input_data['model_path'] 
    model = input_data['model'] 
    media_type = input_data['media_type']
    aeration = input_data['aeration']
    fluxBounds_dict = input_data['fluxBounds_dict'] 
    fluxBounds_filename = input_data['fluxBounds_filename']
    max_biomass = input_data['max_biomass'] 
    max_biomass_percent = input_data['max_biomass_percent'] 
    max_pyrlys_percent = input_data['max_pyrlys_percent']
    results_filename = input_data['results_filename']
    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings']

    # First find max flux of EX_pyrlys_L
    if stdout_msgs:
        print '\n--- FBA to find max EX_pyrlys_L exchane reaction flux ---\n'
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.reactions_by_id['EX_pyrlys_L(e)'].objective_coefficient = 1

    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = fluxBounds_filename, flux_bounds = fluxBounds_dict, reset_flux_bounds = True)
    set_specific_bounds(model = Ecoli_nsAA_producing, flux_bounds = {biomass_reaction.id:[(max_biomass_percent/100)*max_biomass,1000]}, reset_flux_bounds = False)

    model.fba(stdout_msgs = stdout_msgs)

    if model.fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_EX_pyrlys =  model.fba_model.solution['objective_value']
        if stdout_msgs:
            print 'biomass = {} , max Pyrrolysine exch flux = {}  , p-amino-phenylalanine exch. flux = {}'.format(model.fba_model.solution['opt_rxnFluxes'][model.biomass_reaction.id],model.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'],model.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
    else:
        max_EX_pyrlys = 0
        if warnings:
             print '\tInfeasible FBA for finding max_EX_pyrlys at max_biomass_percent = {}'.format(max_biomass_percent) 

    # Find the max p-amino-phenylalanine exchange reaction flux 
    if stdout_msgs:
        print '\n--- FBA to find max p-amino-phenylalanine exchange reaction flux ---\n'
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.reactions_by_id['EX_paphe(e)'].objective_coefficient = 1

    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = fluxBounds_filename, flux_bounds = fluxBounds_dict, reset_flux_bounds = True)
    set_specific_bounds(model = Ecoli_nsAA_producing, flux_bounds = {biomass_reaction.id:[(max_biomass_percent/100)*max_biomass,1000]} , reset_flux_bounds = False)

    model.fba(stdout_msgs = stdout_msgs)

    if model.fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_EX_paphe =  model.fba_model.solution['objective_value']
        if stdout_msgs:
            print 'biomass = {} , Pyrrolysine exch flux = {}  , max p-amino-phenylalanine exch. flux = {}'.format(model.fba_model.solution['opt_rxnFluxes'][model.biomass_reaction.id],model.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'],model.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
    else:
        max_EX_paphe = 0
        if warnings:
             print '\tInfeasible FBA for finding max_EX_paphe at max_biomass_percent = {}'.format(max_biomass_percent) 
        
    # Find the max p-amino-phenylalanine exchange reaction flux at fixed max_biomass_percent and max_pyrlys_percent 
    if stdout_msgs:
        print '\n--- FBA to find max p-amino-phenylalanine exchange reaction flux at fixed max_biomass_percent and max_pyrlys_percent---\n'

    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = fluxBounds_filename, flux_bounds = fluxBounds_dict, reset_flux_bounds = True)
    set_specific_bounds(model = Ecoli_nsAA_producing, flux_bounds = {biomass_reaction.id:[(max_biomass_percent/100)*max_biomass,1000], 'EX_pyrlys_L(e)':[(max_pyrlys_percent/100)*max_EX_pyrlys,1000]} , reset_flux_bounds = False)

    model.fba(stdout_msgs = stdout_msgs)
    if model.fba_model.solution['exit_flag'] == 'globallyOptimal':
        EX_paphe =  model.fba_model.solution['objective_value']
        if stdout_msgs:
            print 'biomass = {} , Pyrrolysine exch flux = {}  , max p-amino-phenylalanine exch. flux = {}'.format(model.fba_model.solution['opt_rxnFluxes'][model.biomass_reaction.id],model.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'],model.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
    else:
        EX_paphe = 0
        if warnings:
             print '\tInfeasible FBA for finding max_EX_paphe at max_biomass_percent = {}'.format(max_biomass_percent) 

    # Store the results        
    if results_filename != '':
        res_key = (('max_biomass_percent',max_biomass_percent),('max_pyrlys_percent',max_pyrlys_percent))
        res_val = (('max_biomass',max_biomass),('max_EX_pyrlys',max_EX_pyrlys),('max_EX_paphe',max_EX_paphe),('biomass',(max_biomass_percent/100)*max_biomass), ('EX_pyrlys',(max_pyrlys_percent/100)*max_EX_pyrlys),('EX_paphe',EX_paphe))
        with open(results_filename,'a') as f:
            f.write("results[" + str(res_key) + "] = " + str(res_val) + "\n")    
 
def master_func(start_pos = None, end_pos = None, max_biomass_percentages = range(0,101,5), max_pyrrolysine_percentages = range(0,101,5), aeration = 'aerobic', media_type = 'minimal',results_filename = '', stdout_msgs = True, warnings = True):
    """
    Creates the model for E. coli harboring the pathways producing the non-standard amino acids
    and assesses its impact on the cell growth 

    INPUTS:
    -------
                      start_pos: Start position of the array containing all possble cases to consider (see all_cases variable)
                        end_pos: End position of the array containing all possble cases to consider (see all_cases variable)
        max_biomass_percentages: A vector containing the percentages of max biomass flux to consider
    max_pyrrolysine_percentages: A vector containing the percentage of max EX_pyrlys flux to consider
                     media_type: Type of growth mediam: M9 (minimal(, LB (rich)
                       aeration: Type of aeration 'aerobic' or 'anaerobic'
    """
    print '\n------------- {} , {} ---------------\n'.format(media_type, aeration)

    if not isinstance(aeration,str):
        raise TypeError('aeration must a a string')
    elif aeration.lower() not in ['aerobic','anaerobic']:
        raise ValueError('Invalid aeration value: {}'.format(aeration))

    if not isinstance(media_type,str):
        raise TypeError('media_type must be a string')
    elif media_type.lower() not in ['minimal','rich']:
        raise ValueError('Invalid media_type value! Allowed choices are: [M9, LB]')

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

    if media_type.lower() == 'minimal' and aeration.lower() == 'aerobic':
        fluxBounds_dict = {'EX_glc(e)':[-10,1000], 'EX_o2(e)':[-20,1000]}
        fluxBounds_filename = model_path + 'iJO1366_minimal_glucose_aerobic.py'
    elif media_type.lower() == 'minimal' and aeration.lower() == 'anaerobic':
        fluxBounds_dict = {'EX_glc(e)':[-10,1000]} 
        fluxBounds_filename = model_path + 'iJO1366_minimal_glucose_anaerobic.py'
    elif media_type.lower() == 'rich' and aeration.lower() == 'aerobic':
        fluxBounds_dict = {'EX_glc(e)':[-10,1000], 'EX_o2(e)':[-20,1000]} 
        fluxBounds_filename = model_path + 'iJO1366_rich_glucose_aerobic.py' 
    elif media_type.lower() == 'rich' and aeration.lower() == 'anaerobic':
        fluxBounds_dict = {'EX_glc(e)':[-10,1000]} 
        fluxBounds_filename = model_path + 'iJO1366_rich_glucose_anaerobic.py' 
    else:
        raise ValueError('Unknown media_type and/or aeraiton type: {}, {}'.format(media_type,aeration))

    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    # Orignal iJo1266 model
    model = create_model(model_organism = model_organism, model_info = {'id':'iJO1366', 'sbml_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_fluxBounds = {'fluxBounds_filename':model_path + fluxBounds_filename, 'fluxBounds_dict': fluxBounds_dict}, stdout_msgs = True, warnings = True)

    # Add the nsAA biosynthesis pathways to the model
    Ecoli_nsAA_producing = add_nsAA_pathways(model = deepcopy(model), stdout_msgs = stdout_msgs, warnings = warnings) 

    # Find the max biomass flux for Ecoli_nsAA_producing 
    print '\n--- FBA after adding new pathways ---\n'
    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = fluxBounds_filename, flux_bounds = fluxBounds_dict)

    Ecoli_nsAA_producing.fba(stdout_msgs = stdout_msgs)
    if Ecoli_nsAA_producing.fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_biomass = Ecoli_nsAA_producing.fba_model.solution['objective_value']
        print 'max biomass = {} , Pyrrolysine exch flux = {}  , p-amino-phenylalanine exch. flux = {}'.format(max_biomass,Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
    else:
        raise userError('FBA to find max biomass flux for Ecoli_nsAA_producing was not solved to optimality')

    # Find the max pyrrolysine exchange reaction flux 
    print '\n--- FBA to find max pyrrolysine exchange reaction flux ---\n'
    for rxn in Ecoli_nsAA_producing.reactions:
        rxn.objective_coefficient = 0
    Ecoli_nsAA_producing.reactions_by_id['EX_pyrlys_L(e)'].objective_coefficient = 1
    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = fluxBounds_filename, flux_bounds = fluxBounds_dict)

    Ecoli_nsAA_producing.fba(stdout_msgs = stdout_msgs)
    if Ecoli_nsAA_producing.fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_EX_pyrlys =  Ecoli_nsAA_producing.fba_model.solution['objective_value']
        print 'biomass = {} , max Pyrrolysine exch flux = {}  , p-amino-phenylalanine exch. flux = {}'.format(Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes'][Ecoli_nsAA_producing.biomass_reaction.id],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
    else:
        raise userError('FBA to find max exchange reaction flux for L-pyrrolysine was not solved to optimality')

    # Find the max p-amino-phenylalanine exchange reaction flux 
    print '\n--- FBA to find max p-amino-phenylalanine exchange reaction flux ---\n'
    for rxn in Ecoli_nsAA_producing.reactions:
        rxn.objective_coefficient = 0
    Ecoli_nsAA_producing.reactions_by_id['EX_paphe(e)'].objective_coefficient = 1
    set_specific_bounds(model = Ecoli_nsAA_producing, file_name = fluxBounds_filename, flux_bounds = fluxBounds_dict)

    Ecoli_nsAA_producing.fba(stdout_msgs = stdout_msgs)
    if Ecoli_nsAA_producing.fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_EX_paphe =  Ecoli_nsAA_producing.fba_model.solution['objective_value']
        print 'biomass = {} , Pyrrolysine exch flux = {}  , max p-amino-phenylalanine exch. flux = {}'.format(Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes'][Ecoli_nsAA_producing.biomass_reaction.id],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_pyrlys_L(e)'],Ecoli_nsAA_producing.fba_model.solution['opt_rxnFluxes']['EX_paphe(e)'])
    else:
        raise userError('FBA to find max exchange reaction flux for p-amino-phenylalanine was not solved to optimality')

    # Set all objective function coefficients to zero
    for rxn in Ecoli_nsAA_producing.reactions:
        rxn.objective_coefficient = 0

    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('results = {}\n')

    all_cases = [(max_biomass_percent, max_pyrlys_percent) for max_biomass_percent in max_biomass_percentages for max_pyrlys_percent in max_pyrrolysine_percentages] 
    print '\nThe total # of cases to consider = {}\n'.format(len(all_cases))

    if start_pos != None and end_pos != None:
        cases_to_consider = all_cases[start_pos - 1:end_pos]
        counter = start_pos - 1
        print 'Simulating slice {}\n'.format((start_pos,end_pos))
    else:
        cases_to_consider = all_cases
        counter = 0
    
    counter = 0

    for (max_biomass_percent, max_pyrlys_percent) in cases_to_consider:

        counter += 1
      
        print '{}. (max_biomass_percent, max_pyrlys_percent) = {}'.format(counter,(max_biomass_percent, max_pyrlys_percent))

        # Creating a shared memory using the manager
        input_data = {} 
        input_data['model_path'] = model_path
        input_data['model'] = Ecoli_nsAA_producing
        input_data['max_biomass'] = max_biomass
        input_data['max_biomass_percent'] = max_biomass_percent
        input_data['max_pyrlys_percent'] = max_pyrlys_percent
        input_data['media_type'] = media_type
        input_data['aeration'] = aeration
        input_data['fluxBounds_dict'] = fluxBounds_dict
        input_data['fluxBounds_filename'] = fluxBounds_filename
        input_data['results_filename'] = results_filename
        input_data['stdout_msgs'] = stdout_msgs
        input_data['warnings'] = warnings

        p = Process(target = findTradeOffs, args = (input_data,))
        p.start()
        p.join() 
        if p.exitcode > 0:
            raise userError('Error in python subprocess. Please check performDMMM\n')

def plot_tradeoff_results(results_filename, title = '', xaxis_label = 'biomass production flux (1/h)', yaxis_label = 'ns Amino Acid produciton flux (mmmol/gDW.h)', media_type = 'minimal',aeration = 'aerobic', figsize = None,output_filename_base = ''):
    """
    Plots the results of tradeoff analysis
    """
    import os
    from imp import load_source
    from tools.ancillary.plot import plot, axis, color_bar
    import matplotlib

    # Load the data from the files
    load_source('dataFile',results_filename)
    import dataFile
    results = dataFile.results

    #----- First make 2D tradeoff plots (production of each nsAA vs. biomass) ---
    # Create a dictionary where keys are biomass fluxes and values are max_EX_pyrlys and max_EX_paphe 
    # at that biomass flux
    max_EX_pyrlys_vs_biomass = {}
    max_EX_paphe_vs_biomass = {}

    # A dictionary whose keys are a tuple containing the biomass and EX_pyrlys reaction fluxes and values are
    # EX_paphe produciotn flux
    EX_paphe_vs_biomassAndEX_pyrlys = {}

    # An examples of results is as follows:
    # results[(('max_biomass_percent', 0), ('max_pyrlys_percent', 0))] = (('max_biomass', 0.950153566914), ('max_EX_pyrlys', 3.88179775281), ('max_EX_paphe', 5.36168757127), ('biomass', 0.0), ('EX_pyrlys', 0.0), ('EX_paphe', 5.36168757127))
    for k in results.keys():
        biomass_flux = dict(results[k])['biomass']

        max_EX_pyrlys_flux = dict(results[k])['max_EX_pyrlys']
        max_EX_pyrlys_vs_biomass[biomass_flux] = max_EX_pyrlys_flux

        max_EX_paphe_flux = dict(results[k])['max_EX_paphe']
        max_EX_paphe_vs_biomass[biomass_flux] = max_EX_paphe_flux

    # Plot the results
    #xaxis_label = r'$\boldsymbol{Biomass \,producton \, flux \, (h^-1)}$'
    xaxis_label = r'$Biomass \,producton \, flux \, (h^{-1})$'
    yaxis_label = r'$Non-standard \, Amino \, Acid$' + '\n' + r'$production \, flux \, (\frac{mmol}{gDW.h})$'
    nsAA2Dplot = plot(title = title,  xaxis = axis(label = xaxis_label, limits = (0,None), spines_format = {'bottom':{'linewidth':2},'top':{'linewidth':2}}), yaxis = axis(label = yaxis_label, limits = (0,None), spines_format = {'left':{'linewidth':2}, 'right':{'linewidth':2}}), plot_gridlines = False, show_legend = True, fig_format = {'use_LaTex':True}, output_filename = output_filename_base + '_2D.pdf')

    # Pyrrolysine
    nsAA2Dplot.plot2D(x = sorted(max_EX_pyrlys_vs_biomass.keys()), y = [max_EX_pyrlys_vs_biomass[k] for k in sorted(max_EX_pyrlys_vs_biomass.keys())], label = 'Pyrrolysine', create_new_figure = True, save_current = False)

    # p-Amino-phenylalanine
    nsAA2Dplot.plot2D(x = sorted(max_EX_paphe_vs_biomass.keys()), y = [max_EX_paphe_vs_biomass[k] for k in sorted(max_EX_paphe_vs_biomass.keys())], label = 'pAF', create_new_figure = False, save_current = False)

    nsAA2Dplot.customize_and_save()


    #---- 3D plot of EX_paphe as a function of biomass and EX_pyrlys production fluxes------
    # x = biomass,  y = EX_pyrlys,  z = EX_paphe
    x_biomass = []
    y_EX_pyrlys = []
    z_EX_paphe = []
    # An examples of results is as follows:
    # results[(('max_biomass_percent', 0), ('max_pyrlys_percent', 0))] = (('max_biomass', 0.950153566914), ('max_EX_pyrlys', 3.88179775281), ('max_EX_paphe', 5.36168757127), ('biomass', 0.0), ('EX_pyrlys', 0.0), ('EX_paphe', 5.36168757127))
    for k in sorted(results.keys() ,key=lambda x:(dict(results[x])['biomass'],dict(results[x])['EX_pyrlys'])):
        biomass_flux = dict(results[k])['biomass']
        EX_pyrlys_flux = dict(results[k])['EX_pyrlys']
        EX_paphe_flux = dict(results[k])['EX_paphe']

        x_biomass.append(biomass_flux)
        y_EX_pyrlys.append(EX_pyrlys_flux)
        z_EX_paphe.append(EX_paphe_flux)
 
        EX_paphe_vs_biomassAndEX_pyrlys[(('biomass',biomass_flux),('EX_pyrlys',EX_pyrlys_flux))] = EX_paphe_flux

    x_biomass = np.array(x_biomass)
    y_EX_pyrlys = np.array(y_EX_pyrlys)
    z_EX_paphe = np.array(z_EX_paphe)

    #-- Plot the results --
    # Specify the ticklabels format
    # x axis
    #if aeration.lower() == 'aerobic':
    #    xticklabels_format = {'rotation':0, 'verticalalignment':'baseline'} 
    #elif aeration.lower() == 'anaerobic':
    #    xticklabels_format = {'rotation':20, 'horizontalalignment':'right','verticalalignment':'baseline'} 
    #xticklabels_format = {'distance_from_ticklabels':20}
    xticklabels_format = {}
    #xlabel_3d = 'Biomass production flux \n(1/h)'    
    xlabel_3d = '' 

    # y axis
    yticklabels_format = {'rotation':-30, 'verticalalignment':'baseline'}
    #yticklabels_format = {'distance_from_ticklabels':20}
    #yticklabels_format = {}
    ylabel_3d = 'Pyrrolysine produciton flux \n(mmol/gDW.h)'
    ylabel_3d = '' 

    # z axis
    #if aeration.lower() == 'aerobic':
    #    zticklabels_format = {'rotation':0, 'horizontalalignment':'right', 'verticalalignment':'baseline'}
    #elif aeration.lower() == 'anaerobic':
    #    zticsklabels_format = {'rotation':0, 'horizontalalignment':'center', 'verticalalignment':'baseline'}
    #zticklabels_format = {'distance_from_ticklabels':20}
    zticklabels_format = {}
    #zlabel_3d = 'pAF production flux \n(mmol/gDW.h)' 
    zlabel_3d = '' 

    nsAA3Dplot = plot(title = title,  xaxis = axis(label = xlabel_3d, limits = (0,None), ticklabels_format = xticklabels_format, invert = True), yaxis = axis(label = ylabel_3d, limits = (0,None), ticklabels_format = yticklabels_format, invert = True), zaxis = axis(label = zlabel_3d, limits = (0,None), ticklabels_format = zticklabels_format), plot_gridlines = True, gridlines_format = {'linewidth':0.5}, output_filename = output_filename_base + '_3D.pdf', fig_format = {'figsize':(8,8)})

    nsAA3Dplot.plot3D(x = x_biomass, y = y_EX_pyrlys, z = z_EX_paphe, line_format = {'width':0}, clrbar = color_bar(colormap = matplotlib.cm.coolwarm))

    nsAA3Dplot.customize_and_save()


def run_plot_tradeoff_results():
    """
    Plots the mortality results
    """
    #--- Minimal ---
    # Aerobic 
    plot_tradeoff_results(results_filename = 'results/tradeoff_minimal_aerobic.py', title = 'Minimal Aerobic', xaxis_label = 'Biomass production flux (1/h)', yaxis_label = 'ns Amino Acid produciton\nflux (mmmol/gDW.h)', aeration = 'Aerobic', output_filename_base = 'results/tradeoff_minimal_aerobic')

    # Aerobic 
    plot_tradeoff_results(results_filename = 'results/tradeoff_minimal_anaerobic.py', title = 'Minimal Anaerobic', xaxis_label = 'Biomass production flux (1/h)', yaxis_label = 'ns Amino Acid produciton\nflux (mmmol/gDW.h)', aeration = 'Anaerobic', output_filename_base = 'results/tradeoff_minimal_anaerobic')

    #--- Rich ---
    # Aerobic 
    plot_tradeoff_results(results_filename = 'results/tradeoff_rich_aerobic.py', title = 'Rich Aerobic', xaxis_label = 'Biomass production flux (1/h)', yaxis_label = 'ns Amino Acid produciton\nflux (mmmol/gDW.h)', aeration = 'Aerobic', output_filename_base = 'results/tradeoff_rich_aerobic')

    # Aerobic 
    plot_tradeoff_results(results_filename = 'results/tradeoff_rich_anaerobic.py', title = 'Rich Anaerobic', xaxis_label = 'Biomass production flux (1/h)', yaxis_label = 'ns Amino Acid produciton\nflux (mmmol/gDW.h)', aeration = 'Anaerobic', output_filename_base = 'results/tradeoff_rich_anaerobic')

 
#-------------------------
if __name__ == '__main__':
    pass 
