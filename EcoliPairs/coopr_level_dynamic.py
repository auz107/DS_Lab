from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../')
sys.path.append('results/')
from copy import deepcopy
import cPickle as pk
import shelve
import itertools
from tools.userError import *
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from DMMM_coopr_level import DMMM_coopr_level
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.ancillary.importData import importData
import cobra
from read_exp_data import read_exp_data
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)
from multiprocessing import Process, Manager
import re

# Last updated: 01-13-2015


def performDMMM(input_data):

    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings']
    t0 = input_data['t0']
    tf = input_data['tf']
    delta_t = input_data['delta_t']
    time_points = input_data['time_points'][:]
    mutant1 = deepcopy(input_data['mutant1'])
    mutant2 = deepcopy(input_data['mutant2'])
    model_path = deepcopy(input_data['model_path'])
    shared_cmp_ids  = input_data['shared_cmp_ids'][:]
    shared_cmp_names  = dict(input_data['shared_cmp_names'][:])
    auxoMetabsMutants_m1  = input_data['auxoMetabsMutants_m1'][:]
    auxoMetabsMutants_m2  = input_data['auxoMetabsMutants_m2'][:]
    coopr_exchrxns_fluxRanges = input_data['coopr_exchrxns_fluxRanges'] 
    coopr_level_m1 = input_data['coopr_level_m1']
    coopr_level_m2 = input_data['coopr_level_m2']
    glc_conc = input_data['glc_conc']
    glc_excess = input_data['glc_excess']
    use_DMMM_coopr = input_data['use_DMMM_coopr']
    save_with_key = input_data['save_with_key']
    save_details = input_data['save_details'] 
    results_filename = input_data['results_filename']

    # Make a list of reacitons instead of a list of lists
    auxoMetabsMutants_m1 = [exch_rxn for exch_rxn_list in auxoMetabsMutants_m1 for exch_rxn in exch_rxn_list]
    auxoMetabsMutants_m2 = [exch_rxn for exch_rxn_list in auxoMetabsMutants_m2 for exch_rxn in exch_rxn_list]

    # List of exchange reaction ids for all shared compounds
    shared_cmps_exch_rxn_ids = []
    for shared_cmp_id in shared_cmp_ids:
        # FInd the exchange reaction for this shared compound
        exch_rxn_ids = [r.id for r in mutant1.reactions if r.type.lower() == 'exchange' and r.reactants[0].id == shared_cmp_id]
        if len(exch_rxn_ids) == 0:
            raise userError('No exchange reaction was found for ' + shared_cmp_id)
        elif len(exch_rxn_ids) > 1:
            raise userError('More than exchange reaction was found for ' + shared_cmp_id + ': ' + repr(exch_rxn_ids))
        else:
            shared_cmps_exch_rxn_ids.append(exch_rxn_ids[0])
    
    # Exchange rxns for all shared compounds each mutant takes up 
    uptake_rxns_m1 = []
    uptake_rxns_m2 = []
    # Exchange rxns for compounds each mutant should produce to help its partner 
    export_rxns_m1 = []
    export_rxns_m2 = []

    for exch_rxn_id in shared_cmps_exch_rxn_ids:
    # Consider each set of compounds rescuing a mutant separately
        if not glc_excess:
            if exch_rxn_id in auxoMetabsMutants_m1 or exch_rxn_id == 'EX_glc(e)':
                uptake_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
            if exch_rxn_id in auxoMetabsMutants_m2 or exch_rxn_id == 'EX_glc(e)':
                uptake_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 
        else:
            if exch_rxn_id in auxoMetabsMutants_m1:
                uptake_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
            if exch_rxn_id in auxoMetabsMutants_m2:
                uptake_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 

        if exch_rxn_id not in auxoMetabsMutants_m1 and exch_rxn_id != 'EX_glc(e)':
            export_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
        if exch_rxn_id not in auxoMetabsMutants_m2 and exch_rxn_id != 'EX_glc(e)':
            export_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 

    if stdout_msgs:
        print mutant1.id,':'
        print '\t\tupdate_rxns = ',[r.id for r in uptake_rxns_m1]
        print '\t\texport = ',[r.id for r in export_rxns_m1]
        print '\n',mutant2.id,':'
        print '\t\tupdate_rxns = ',[r.id for r in uptake_rxns_m2]
        print '\t\texport = ',[r.id for r in export_rxns_m2]

    mutant1.coopr_exchrxns = export_rxns_m1
    mutant2.coopr_exchrxns = export_rxns_m2
    mutant1.cooperation_level = coopr_level_m1
    mutant2.cooperation_level = coopr_level_m2

    # Create the list of shared compounds
    shared_cmps = []

    # Glucose as a shared compound (concentration in mM)
    if not glc_excess: 
        glucose = compound(id = 'glc_D_e', name = 'D-Glucose', KEGG_id = 'C00031', reactant_reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],concentration = {0:glc_conc})    

        shared_cmps.append(glucose)

    for shared_cmp_id in [id for id in shared_cmp_ids if id != 'glc_D_e']:
        reactant_rxns = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
        product_rxns = [r for r in export_rxns_m1 + export_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
        shared_cmp = compound (id = shared_cmp_id, name = shared_cmp_names[shared_cmp_id],reactant_reactions = reactant_rxns, product_reactions = product_rxns, concentration = {0:0}) 
        shared_cmps.append(shared_cmp)
        if stdout_msgs:
            print shared_cmp.id,': '
            print '\treactant reactions = ',[(r.model.id,r.id) for r in shared_cmp.reactant_reactions]
            print '\tproduct reactions = ',[(r.model.id,r.id) for r in shared_cmp.product_reactions]
            print '\treactions = ',[(r.model.id,r.id) for r in shared_cmp.reactions]

    # Set store_flux to True for all exchange reacitons related to shared metabolites
    for rxn in [r for c in shared_cmps for r in c.reactions]:
        rxn.store_flux = True
 
    for rxn in [m.biomass_reaction for m in [mutant1,mutant2]]: 
        rxn.store_flux = True

    # Growth medium
    if not glc_excess:
        if use_DMMM_coopr:
            set_specific_bounds(model = mutant1, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py', flux_bounds = dict([(rid,[0,0]) for rid in  mutant1.knockedout_rxn_ids]))
            set_specific_bounds(model = mutant2, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py', flux_bounds = dict([(rid,[0,0]) for rid in  mutant2.knockedout_rxn_ids]))
        else:
            set_specific_bounds(model = mutant1, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',flux_bounds = dict([(rid,[0,0]) for rid in  mutant1.knockedout_rxn_ids] + [(export_rxn.id,[(coopr_level_m1/100)*coopr_exchrxns_fluxRanges[export_rxn.id][1],1000]) for export_rxn in export_rxns_m1]))
            set_specific_bounds(model = mutant2, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',flux_bounds = dict([(rid,[0,0]) for rid in  mutant2.knockedout_rxn_ids] + [(export_rxn.id,[(coopr_level_m2/100)*coopr_exchrxns_fluxRanges[export_rxn.id][1],1000]) for export_rxn in export_rxns_m2]))
    else:
        if use_DMMM_coopr:
            set_specific_bounds(model = mutant1, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py', flux_bounds = dict([(rid,[0,0]) for rid in  mutant1.knockedout_rxn_ids] + [('EX_glc(e)',[-10,1000])]))
            set_specific_bounds(model = mutant2, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py', flux_bounds = dict([(rid,[0,0]) for rid in  mutant2.knockedout_rxn_ids] + [('EX_glc(e)',[-10,1000])]))
        else:
            set_specific_bounds(model = mutant1, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',flux_bounds = dict([(rid,[0,0]) for rid in  mutant1.knockedout_rxn_ids] + [('EX_glc(e)',[-10,1000])] + [(export_rxn.id,[(coopr_level_m1/100)*coopr_exchrxns_fluxRanges[export_rxn.id][1],1000]) for export_rxn in export_rxns_m1]))
            set_specific_bounds(model = mutant2, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',flux_bounds = dict([(rid,[0,0]) for rid in  mutant2.knockedout_rxn_ids] + [('EX_glc(e)',[-10,1000])] + [(export_rxn.id,[(coopr_level_m2/100)*coopr_exchrxns_fluxRanges[export_rxn.id][1],1000]) for export_rxn in export_rxns_m2]))

   
    if stdout_msgs:
        print '\n--- FBA for mutant 1 after resetting the bounds ---'
        mutant1.compounds_by_id['glc_D_e'].concentration = 111.01
        mutant1.reactions_by_id['EX_glc(e)'].kinetic_rate_calc(assignLB = True)
        print 'm1: LB EX_glc(e) = ',mutant1.reactions_by_id['EX_glc(e)'].flux_bounds[0]
        mutant1.fba()

        print '\n--- FBA for mutant 2 after resetting the bounds ---'
        mutant2.compounds_by_id['glc_D_e'].concentration = 111.01
        mutant2.reactions_by_id['EX_glc(e)'].kinetic_rate_calc(assignLB = True)
        print 'm2: LB EX_glc(e) = ',mutant2.reactions_by_id['EX_glc(e)'].flux_bounds[0]
        mutant2.fba()
        print '\nshared metabs ids = ',shared_cmp_ids,'\n'

    if use_DMMM_coopr:
        DMMM_m1m2 = DMMM_coopr_level(community_members = [mutant1,mutant2],shared_compounds = shared_cmps, time_range = [t0,delta_t,tf], store_dynamic_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)
    else:
        DMMM_m1m2 = DMMM(community_members = [mutant1,mutant2],shared_compounds = shared_cmps, time_range = [t0,delta_t,tf], store_dynamic_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)
    DMMM_m1m2.run()

    #--- Store the results ---
    if results_filename != '':
        results = {}
        # Store only the required information for each shared compound
        mutants_names = (mutant1.organism.id,mutant2.organism.id)
        cell_concs_stored = {mutant1.organism.id:mutant1.organism.cells_per_ml,mutant2.organism.id:mutant2.organism.cells_per_ml}
        results['cell_concs'] = cell_concs_stored

        if save_details:
            shared_cmps_concs_stored = {}
            for shared_cmp in [c for c in shared_cmps if c.id != 'glc_D_e' and max(c.concentration.values()) > 1e-9]:
                shared_cmps_concs_stored[shared_cmp.id] = shared_cmp.concentration
            if stdout_msgs:
                print '\nshared compounds with non-zero concentration = ',shared_cmps_concs_stored.keys()
            results['shared_cmps_concs'] = shared_cmps_concs_stored
    
            # Shared the flux of exchange (uptake) reactions for any shared metabolite, which is non-zero
            exch_rxns_stored = {}
            for rxn in list(set([r for c in shared_cmps for r in c.reactions if max([abs(v) for v in r.flux.values()]) > 1e-9])):
                exch_rxns_stored[(rxn.model.organism.id,rxn.id)] = rxn.flux
            results['exch_rxns'] = exch_rxns_stored
            
        # The following outputs are saved: shared_cmps and organism object 
    # for each community member 
        with open(results_filename ,'a') as f:
            if save_with_key:
                save_key = ((mutant1.organism.id + '_coopr_level',coopr_level_m1), (mutant2.organism.id + '_coopr_level', coopr_level_m2))
                f.write('results[' + str(mutants_names) + '][' + str(save_key) + '] = ' + str(results) + '\n')
            else:
                f.write('results[' + str(mutants_names) + '] = ' + str(results) + '\n')

def create_model(warnings = True, stdout_msgs = True):
    """
    Create the metabolic model
    """

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655',gDW_per_cell = 2.8e-13)

    model = read_sbml_model(file_name = model_path + 'iJO1366_updated.xml', model_id = 'iJO1366',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    model.biomass_reaction = model.reactions_by_id['Ec_biomass_iJO1366_core_53p95M']
    model.all_biomass_reactions = {'core':model.reactions_by_id['Ec_biomass_iJO1366_core_53p95M'],'WT':model.reactions_by_id['Ec_biomass_iJO1366_WT_53p95M']}

    # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
    # Example: EX_glc(E): glc_D_e <==>    Vmax*C['glc_D_e']/(Km + C['glc_D_e']) 
    # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
    for reaction in [r for r in model.reactions if r.type.lower() == 'exchange']:
        # The id of compound participating in the exchange reaction
        metab_id = [m.id for m in reaction.compounds][0]
        reaction.kinetics = "10*C['" + metab_id + "']/(10 + C['" + metab_id + "'])"

    # Glucose uptake kinetics 
    exch_rxns = model.get_reactions({'EX_glc(e)':'id','EX_lys_L(e)':'id','EX_ile_L(e)':'id'})
    exch_rxns['EX_glc(e)'].kinetics = "10*C['glc_D_e']/(0.15 + C['glc_D_e'])"
    exch_rxns['EX_lys_L(e)'].kinetics = "0.1964*C['lys_L_e']/(5e-4 + C['lys_L_e']) + 0.3055*C['lys_L_e']/(1e-2 + C['lys_L_e'])"
    exch_rxns['EX_ile_L(e)'].kinetics = "0.0346*C['ile_L_e']/(1.22e-3 + C['ile_L_e'])"
   
    # Assign the objective function coefficients
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.biomass_reaction.objective_coefficient = 1

    # Growth medium
    set_specific_bounds(model = model, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',limiting_nutrients = {'EX_glc(e)':[-10,1000]})

    # Perform FBA for the wild-type
    model.fba(assign_wildType_max_biomass = True, stdout_msgs = stdout_msgs)

    return model

def find_coopr_exchrxns_fluxRanges(results_filename = 'coopr_exchrxns_fluxRanges.py', warnings = True, stdout_msgs = True):
    """
    Finds the min and max value of each compound that can be produced by the wild-type strain
    """
    from auxoMetabs import auxoMetabsMutants

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

    # Exchange reactions needed for cooperaiton in all pairs
    cooperative_exchrxn_ids = list(set([r for m in auxoMetabsMutants.keys() for r_list in auxoMetabsMutants[m] for r in r_list]))
  
    WT = create_model(warnings = warnings, stdout_msgs = stdout_msgs)

    # Growth medium
    set_specific_bounds(model = WT, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',limiting_nutrients = {'EX_glc(e)':[-10,1000]})

    with open(results_filename,'w') as f:
        f.write('coopr_exchrxns_fluxRanges = {}\n')

    for exchrxn_id in cooperative_exchrxn_ids:

        if stdout_msgs:
            print '\n------- ',exchrxn_id,' ---------\n'

        exchrxn = WT.reactions_by_id[exchrxn_id]

        # Find maximum
        for rxn in WT.reactions:
            rxn.objective_coefficient = 0
        exchrxn.objective_coefficient = 1
 
        WT.fba(store_opt_fluxes = False, warnings = warnings, stdout_msgs = stdout_msgs)
        if WT.fba_model.solution['exit_flag'] == 'globallyOptimal':
             max_exch = WT.fba_model.solution['objective_value']
        else:
             max_exch = 0

        # Find minimum
        for rxn in WT.reactions:
            rxn.objective_coefficient = 0
        exchrxn.objective_coefficient = -1
 
        WT.fba(store_opt_fluxes = False, warnings = warnings, stdout_msgs = stdout_msgs)
        if WT.fba_model.solution['exit_flag'] == 'globallyOptimal':
             min_exch = -WT.fba_model.solution['objective_value']
        else:
             min_exch = 0

        if min_exch > max_exch:
            raise userError('min_exch = {} > max_exch = {}'.format(min_exch,max_exch))

        with open(results_filename,'a') as f:
            f.write("coopr_exchrxns_fluxRanges['" + exchrxn_id + "'] = (" + str(min_exch) + ',' + str(max_exch) + ')\n') 

    print '\nResults were written into coopr_exchrxns_fluxRanges.py ...\n'

def get_cmp_name_from_SBML(model,cmp_id_list):
    """
    INPUTS:
    ------
    cmp_id_list: A list of compound ids

    OUTPUTS:
    --------
    cmp_names: A dictionary where keys are compound ids and values are their names
               extracted from the SBML file
    """
    cmp_names = {}
    counter = 0
    for cmp_id in cmp_id_list:
        sbml_cmp = [c for c in model.compounds if remove_non_alphanumeric(cmp_id).lower() == remove_non_alphanumeric(c.id).lower()]
        if len(sbml_cmp) == 1:
            cmp_names[cmp_id] = sbml_cmp[0].name
            counter += 1
        # otherwise assign the id as its name
        else:
            cmp_names[cmp_id] = cmp_id

    return cmp_names


def master_func(t0,delta_t,tf,start_pos,end_pos,start_pos_pair = None, end_pos_pair = None, cooperation_level_m1 = 0, cooperation_level_m2 = 0, glc_conc = 111.01, glc_excess = False, per_capita_cell_conc_init = 7.5e6, save_details = True, save_with_key = False, use_DMMM_coopr = True, results_filename_base = '', stdout_msgs = True, warnings = True):
    """
     INPUTS:
     -------
                             t0: Initial simulaiton time
                        delta_t: Time step
                             tf: Total simulation time
                      strat_pos: Start position of the array containing all possible
                                 mutant pairs to examine                                 
                        end_pos: End position of the array containing all possible
                                 mutant pairs to examine                                 
                 strat_pos_pair: Start position of the array containing all possible
                                 cases to examine for a given mutant pair 
                   end_pos_pair: End position of the array containing all possible
                                 cases to examine for a given mutant pair 
           cooperation_level_m1: Level of cooperation (percentage of maximum produciton flux 
                                 of compounds needed by their partner) for the first mutant 
           cooperation_level_m1: Level of cooperation (percentage of maximum produciton flux 
                                 of compounds needed by their partner) for the second mutant
                       glc_conc: Glucose concentraiton
                     glc_excess: If true glc is not added to the set of shared compounds and its uptake flux is always set to -10
      per_capita_cell_conc_init: Initial per captial cell concentration
                   save_details: A parameter showing whether the details of DMMM simulations (such as exchange flux of
                                 shared compounds) should be saved (True) or not (False)
                  save_with_key: A parameter showing whether the results have to be saved with keys 
                                 (('mortality_percentage',A),('pool_conc_factor',B)), where A and B are some numbers showing
                                 the specified mortality percetage and cell_pool_conc_factors respectively.
                 use_DMMM_coopr: If True DMMM_coopr_level is used, otherwise the original DMMM is used. The former is suitable
                                 for simulations with excess glucose while the latter is better for those with limitted glucose 
                                 (i.e., when glucose concentration is provided)
          results_filename_base: The base for the file name storing the results
                                 Example: results_filename_base = 'results/emc_results'. The code
                                 will add the start and end positions to the file name.
                                 Example: 'results/emc_results_1_500.txt'

    NOTE: The user can enter both the start and end  positions numbers 
          assuming that they start at one. The code will take care of
          array indexing convention of python (indexing starts at zero) 
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs mut be either True to False')

    if not isinstance(warnings,bool):
        raise TypeError('warnings mut be either True to False')

    if not isinstance(save_details,bool):
        raise TypeError('save_details mut be either True to False')

    if not isinstance(save_with_key,bool):
        raise TypeError('save_with_key mut be either True to False')
    
    if not isinstance(glc_excess,bool):
        raise TypeError('glc_excess mut be either True to False')
    
    if not isinstance(use_DMMM_coopr,bool):
        raise TypeError('use_DMMM_coopr mut be either True to False')
    
    if isinstance(cooperation_level_m1,int) or isinstance(cooperation_level_m1,float):
        cooperation_level_m1 = [cooperation_level_m1]
    elif not isinstance(cooperation_level_m1,list):
        raise TypeError('cooperation_level_m1 must be an integer, float, or a list of integers and floats') 
    if len([cl for cl in cooperation_level_m1 if cl > 100]):
        raise ValueError('All entries of cooperation_level_m1 must be less than or equal to 100 as they are percentages.')

    if isinstance(cooperation_level_m2,int) or isinstance(cooperation_level_m2,float):
        cooperation_level_m2 = [cooperation_level_m2]
    elif not isinstance(cooperation_level_m2,list):
        raise TypeError('cooperation_level_m2 must be an integer, float, or a list of integers and floats') 
    if len([cl for cl in cooperation_level_m2 if cl > 100]):
        raise ValueError('All entries of cooperation_level_m2 must be less than or equal to 100 as they are percentages.')

    if not isinstance(per_capita_cell_conc_init,int) and not isinstance(per_capita_cell_conc_init,float):
        raise TypeError('per_capita_cell_conc_init must be an integer, float') 

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]    
    print '\ntime_points = ',time_points

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'
    WT = create_model(stdout_msgs = stdout_msgs, warnings = warnings)

    # Compute ms and biomass yield for glucose
    glc_D = WT.get_compounds({'glc_D_e':'id'})
    glc_D.ms_calc() 
    glc_D.biomass_yield_calc() 

    #--- Load the list of rxns that must be off in each mutant ---
    from mutants_rxn_info import mutants_rxn_info

    #--- Load the list of exchange rxns for compounds each mutant needs to survive ---
    from auxoMetabs import auxoMetabsMutants
    # Mutants that cannot be rescued at all
    not_rescued_mutants = [m for m in auxoMetabsMutants.keys() if len(auxoMetabsMutants[m]) == 0]

    # Load the min and max flux of exchange reactions used for cooperation
    from coopr_exchrxns_fluxRanges import coopr_exchrxns_fluxRanges
    
    # All possible pair combinations (do not consider the ones that cannot be rescued)
    mutant_pairs = list(itertools.combinations([m for m in mutants_rxn_info.keys() if m not in not_rescued_mutants],r=2))   
    print '\nThe total # of mutant pairs to examine = %i' % len(mutant_pairs)

    # IDs of the shared compounds
    if not glc_excess:
        shared_cmp_ids = ['glc_D_e']
    else:
        shared_cmp_ids = []

    # All possible cases to consider for each pair
    all_cases_m1m2 = [(coopr_level_m1,coopr_level_m2) for coopr_level_m1 in cooperation_level_m1 for coopr_level_m2 in cooperation_level_m2] 

    # Name of the output file storing the results
    if results_filename_base != '':
        if start_pos_pair != None and end_pos_pair != None:
            results_filename = results_filename_base + '_' + str(start_pos) + '_' + str(end_pos) + '_' + str(start_pos_pair) + '_' + str(end_pos_pair) + '.py'
        else:
            results_filename = results_filename_base + '_' + str(start_pos) + '_' + str(end_pos) + '.py'
    else:
        results_filename = '' 

    # Initialize the file
    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('results = {}\n')

    print 'Simulating slice {}\n'.format((start_pos,end_pos))

    #--- DMMM for the co-culture of mutant1 and mutant2 mutants ---
    #for (m1,m2) in [('lysA','ilvE')]:
    #for (m1,m2) in [('glnA','trpC'),('lysA','ilvE')]:
    for (m1,m2) in mutant_pairs[start_pos - 1:end_pos]:
      
        print '\n**** %i. (%s,%s) ****\n' % (mutant_pairs.index((m1,m2)) + 1,m1,m2)
    
        #-- Mutant 1 and the related reactions whose flux should be set to zero --
        if stdout_msgs:
            print '\n-- ' + m1 + '_Ecoli mutant --'
        mutant1 = deepcopy(WT)
        mutant1.id = WT.id + '_' + m1
        mutant1.organism.id = m1 + '_Ecoli'
        mutant1.organism.cells_per_ml = {0:per_capita_cell_conc_init}
        mutant1.organism.gDW_per_ml = {0:per_capita_cell_conc_init*mutant1.organism.gDW_per_cell}
        mutant1.knockedout_rxn_ids =  mutants_rxn_info[m1]
        for rxn_id in mutants_rxn_info[m1]:
            rxn = mutant1.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]

        mutant1.fba(build_new_optModel = False, store_opt_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)
 
        print '** Test how much ile_L this mutant can produce ...'
        EX_ile_L = mutant1.reactions_by_id['EX_ile_L(e)'] 
        for rxn in mutant1.reactions:
            rxn.objective_coefficient = 0
        EX_ile_L.objective_coefficient = 1
        mutant1.fba(store_opt_fluxes = False, warnings = warnings, stdout_msgs = stdout_msgs)
        for rxn in mutant1.reactions:
            rxn.objective_coefficient = 0
        mutant1.biomass_reaction.objective_coefficient = 1

        # Compute ms and biomass yield for compounds needed to rescue the mutant 
        for exch_rxn_id in list(set([r for rList in auxoMetabsMutants[m1] for r in rList])): 
            # The compound participating in the exchange reaction
            exch_rxn = mutant1.get_reactions({exch_rxn_id:'id'})
            cmp = exch_rxn.reactants[0] 
            shared_cmp_ids.append(cmp.id)
            cmp.ms_calc()    
            cmp.biomass_yield_calc()    
    
        # Exchange reactions for cooperation (i.e., to export compounds rescuing the partner)
        mutant1.cooperative_rxns = []

        #-- Mutant 2 and related reaction whose flux should be set to zero --
        if stdout_msgs:
            print '\n-- ' + m2 + '_Ecoli mutant --'
        mutant2 = deepcopy(WT)
        mutant2.id = WT.id + '_' + m2
        mutant2.organism.id = m2 + '_Ecoli'
        mutant2.organism.cells_per_ml = {0:per_capita_cell_conc_init}
        mutant2.organism.gDW_per_ml = {0:per_capita_cell_conc_init*mutant2.organism.gDW_per_cell}
        mutant2.knockedout_rxn_ids =  mutants_rxn_info[m2]
        for rxn_id in mutants_rxn_info[m2]:
            rxn = mutant2.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]
            
        mutant2.fba(build_new_optModel = False, store_opt_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)
    
        print '** Test how much lys_L this mutant can produce ...'
        EX_lys_L = mutant2.reactions_by_id['EX_lys_L(e)'] 
        for rxn in mutant2.reactions:
            rxn.objective_coefficient = 0
        EX_lys_L.objective_coefficient = 1
        mutant2.fba(store_opt_fluxes = False, warnings = warnings, stdout_msgs = stdout_msgs)
        for rxn in mutant2.reactions:
            rxn.objective_coefficient = 0
        mutant2.biomass_reaction.objective_coefficient = 1

        # Compute ms and biomass yield for compounds needed to rescue the mutant 
        for exch_rxn_id in list(set([r for rList in auxoMetabsMutants[m2] for r in rList])): 
            exch_rxn = mutant2.get_reactions({exch_rxn_id:'id'})
            # The compound participating in the exchange reaction
            cmp = exch_rxn.reactants[0] 
            shared_cmp_ids.append(cmp.id)
            cmp.ms_calc()    
            cmp.biomass_yield_calc()    
    
        # Exchange reactions for cooperation (i.e., to export compounds rescuing the partner)
        mutant2.cooperative_rxns = []

        if save_with_key:
            with open(results_filename,'a') as f:
                f.write('results[' + str((mutant1.organism.id,mutant2.organism.id)) + '] = {}\n')
    
        # --- Define the compounds available in the extracellular medium ---
        # Consider only unique elements as, in general, it is quite possible that
        # two mutants need the same compounds to survive
        shared_cmp_ids = sorted(list(set(shared_cmp_ids)))

        # Get the names of the shared compounds from the SBML file
        shared_cmp_names = get_cmp_name_from_SBML(model = WT, cmp_id_list = shared_cmp_ids)

        # Simulation scenarios for each mutant
        print '\nThe total # of cases to consider for mutant pair ({},{}) = {}\n'.format(m1,m2,len(all_cases_m1m2))
        
        if start_pos_pair != None and end_pos_pair != None:
            cases_to_consider = all_cases_m1m2[start_pos_pair - 1:end_pos_pair]
            counter = start_pos_pair - 1
            print 'Simulating slice {}\n'.format((start_pos_pair,end_pos_pair))
        else:
            cases_to_consider = all_cases_m1m2
            counter = 0

        for (coopr_level_m1, coopr_level_m2) in cases_to_consider:

            counter += 1
          
            print '{}. (coopr_leve_m1, coopr_level_m2) = {}'.format(counter,(coopr_level_m1, coopr_level_m2))

            # Creating a shared memory using the manager
            input_data = {} 
            input_data['stdout_msgs'] = stdout_msgs
            input_data['warnings'] = warnings
            input_data['t0'] = t0
            input_data['tf'] = tf
            input_data['delta_t'] = delta_t
            input_data['time_points'] = time_points
            input_data['model_path'] = model_path
            input_data['mutant1'] = mutant1
            input_data['mutant2'] = mutant2
            input_data['shared_cmp_ids'] = shared_cmp_ids 
            input_data['shared_cmp_names'] = shared_cmp_names.items() 
            input_data['auxoMetabsMutants_m1'] = auxoMetabsMutants[m1] 
            input_data['auxoMetabsMutants_m2'] = auxoMetabsMutants[m2] 
            input_data['coopr_exchrxns_fluxRanges'] = coopr_exchrxns_fluxRanges
            input_data['coopr_level_m1'] = coopr_level_m1 
            input_data['coopr_level_m2'] = coopr_level_m2 
            input_data['glc_conc'] = glc_conc 
            input_data['glc_excess'] = glc_excess 
            input_data['use_DMMM_coopr'] = use_DMMM_coopr 
            input_data['save_with_key'] = save_with_key 
            input_data['save_details'] = save_details 
            input_data['results_filename'] = results_filename
    
            p = Process(target = performDMMM, args = (input_data,))
            p.start()
            p.join() 
            if p.exitcode > 0:
                raise userError('Error in python subprocess. Please check performDMMM\n')

def run_master_func(start_pos_pair = None,end_pos_pair = None, coopr_level_m1 = [1], coopr_level_m2 = [1], glc_conc = 111.01, glc_excess = False, use_DMMM_coopr = True, run_main = False, run_test = True, results_filename_base = ''):
    """
    Runs the master_func
    
    INPUTS:
    -----
    run_main: If True the main job runs
    run_test: If True the test runs 
    """
 
    if run_main: 
        master_func(t0 = 0,delta_t = 0.5,tf = 96,start_pos = 652,end_pos = 652, start_pos_pair = start_pos_pair, end_pos_pair = end_pos_pair, cooperation_level_m1 = coopr_level_m1, cooperation_level_m2 = coopr_level_m2, glc_conc = glc_conc, glc_excess = glc_excess, use_DMMM_coopr = use_DMMM_coopr, save_details = False, save_with_key = True, results_filename_base = results_filename_base, stdout_msgs = False)

    if run_test:
        #coopr_level_m1 = [10,20] 
        #coopr_level_m2 = [10,20] 
        coopr_level_m1 = [0] 
        coopr_level_m2 = [0] 
        master_func(t0 = 0,delta_t = 0.5,tf = 96,start_pos = 652,end_pos = 652, start_pos_pair = 1, end_pos_pair = 1, cooperation_level_m1 = coopr_level_m1, cooperation_level_m2 = coopr_level_m2, glc_conc = 111.01, glc_excess = False, use_DMMM_coopr = True, save_details = True, save_with_key = True, results_filename_base = 'results/test_coopr_level', stdout_msgs = True)

def create_intervals(total_comb_num,interval_size):
    """
      total_comb_num: Length of the variable 'combinations'
    interval_size: Desired size of the iteration intervals

    Results are printed in the output
    """
    slices = []

    # Start positions
    start_positions = range(1,total_comb_num + 1,interval_size)
    for k in start_positions:
            if k + (interval_size - 1) <= total_comb_num:
                slices.append((k,k + (interval_size - 1)))
            else:
                slices.append((k,total_comb_num))

    print '\n'
    print '**Total # of intervals = ',len(slices),'\n'

    for slice in slices:
       print slice

    return slices

def create_job_files(total_comb_num,interval_size, job_filename_base, joboutput_filename_base, max_walltime, coopr_level_m1_str, coopr_level_m2_str, glc_conc, glc_excess, use_DMMM_coopr):
    """
    Creates the job files given:
             total_comb_num: Total number of combinations 
              interval_size: The desired of iteration intervals
          outfile_base_name: The base name of the output files
                       task: The task that should be performed ('analyze_games' or 'analyze_cost')
          job_filename_base: Base name of the output job file
    joboutput_filename_base: The name of the file storing the dynamic output of the job. It is essentially the same as job_filename_base
                             with the file path deleted (in order to create .out files in the same directory as the job file)
    """
    if not glc_excess:
        if glc_conc == 111.01:
            results_filename_base = 'results/coopr_level_expGlc'
        elif glc_conc == 1000:
            results_filename_base = 'results/coopr_level_fixedGlc'
        else:
            raise userError('Unknonw glc_conc:' + str(glc_conc))
    else:
        results_filename_base = 'results/coopr_level_excessGlc'

    coopr_level_m1 = range(0,101,5)

    slices = create_intervals(total_comb_num = total_comb_num,interval_size = interval_size)

    for slice in slices:
        job_filename = job_filename_base + '_' + str(slice[0]) + '_' + str(slice[1]) + '.sh'
        joboutput_filename = joboutput_filename_base + '_' + str(slice[0]) + '_' + str(slice[1]) + '.out'
        print 'creating file ',job_filename,'...'
        with open(job_filename,'w') as outfile:
            outfile.write('#!/bin/bash -l\n')
            outfile.write('#$-l h_rt=' + str(int(max_walltime))+ ':00:00\n\n')
            outfile.write('cd /usr2/postdoc/alizom/work/EcoliPairs\n\n')

            # This is to merge the .e[jobid] and .o[jobid] files
            outfile.write('#$ -j y\n')

            # Set the output file
            outfile.write('\n#$ -o ' + joboutput_filename + '\n')

            outfile.write('\nsource activate /projectnb/bioinfor/SEGRE/alizom\n\n')

            # python -c "import time;print '\n**Job started at ',time.strftime('%c'),'\n'" > job_dynamic_yeast_stsrt_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\"\n")

            start_pos = slice[0]
            end_pos = slice[1]

            outfile.write('\npython -c \"from coopr_level import run_master_func;run_master_func(start_pos_pair = ' + str(start_pos) + ', end_pos_pair = ' + str(end_pos) + ', coopr_level_m1 = ' + coopr_level_m1_str + ', coopr_level_m2 = ' + coopr_level_m2_str + ', glc_conc = ' + str(glc_conc) + ", glc_excess = " + str(glc_excess) + ", use_DMMM_coopr = " + str(use_DMMM_coopr) + ", run_main = True, run_test = False, results_filename_base = '" + results_filename_base + "')\"\n\n")

            # python -c "import time;print '\n**Job ended at ',time.strftime('%c'),'\n'" >> job_dynamic_yeast_start_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n")


def run_create_job_files():
    """
    Runs create_job_file()
    """

    coopr_level_m1 = range(0,101,5)
    coopr_level_m1_str = 'range(0,101,5)'
    coopr_level_m2 = range(0,101,5)
    coopr_level_m2_str = 'range(0,101,5)'

    total_comb_num = len(coopr_level_m1)*len(coopr_level_m1)
    print '\ntotal_comb_num = ',total_comb_num

    # Experimental glucose simulations
    glc_conc = 111.01
    create_job_files(total_comb_num = total_comb_num, interval_size = 25, job_filename_base = 'jobs/job_coopr_level_expGlc', joboutput_filename_base = 'job_coopr_level_expGlc', max_walltime = 72, coopr_level_m1_str = coopr_level_m1_str, coopr_level_m2_str = coopr_level_m2_str, glc_conc = glc_conc, glc_excess = False, use_DMMM_coopr = True)

    # Excess glcuose simulations
    create_job_files(total_comb_num = total_comb_num, interval_size = 25, job_filename_base = 'jobs/job_coopr_level_excessGlc', joboutput_filename_base = 'job_coopr_level_excessGlc', max_walltime = 72, coopr_level_m1_str = coopr_level_m1_str, coopr_level_m2_str = coopr_level_m2_str, glc_conc = glc_conc, glc_excess = True, use_DMMM_coopr = False)
    
def integrate_results_files(results_filenames, results_var_name, mutant_pair, output_file_name):
    """
    Integrates the results of all scripts (when splitting one job to multiple smaller jobs) inot one sinble file 
 
    INPUTS:
    -------
    results_file_names: A list of strings containing the names of the files containing the results
      results_var_name: Name of the variable storing the results
      output_file_name: Name of the output file name containing the integration of all results
           mutant_pair: A tuple containiing the names of the mutant pairs
    """
    from imp import load_source

    # sum of the entries in all result files
    entries_sum = 0

    results = {}
    results[mutant_pair] = {}
    # Import the data in the module stored in file_name
    for file_name in results_filenames:
        if type(file_name) == str:
            if not os.path.isfile(file_name):
                raise IOError("No such file was found :'" + file_name + "'")
            else:
                # First delete the model dataFile if it already exists. If it is not deleted
                # the new module data is merged with the previous ones
                try:
                    del sys.modules['dataFile']
                except:
                    pass
                load_source('dataFile',file_name)
                import dataFile
        exec('results_curr = dataFile.' + results_var_name)
        
        entries_sum += len(results_curr[mutant_pair].keys())
        print 'Total # of entries in ', file_name,' = ',len(results_curr[mutant_pair].keys())

        results_keys = results[mutant_pair].keys()
        for k in results_curr[mutant_pair].keys():
            if k in results_keys:
                raise userError(str(k) + ' in ' + file_name + ' already in results_keys\n')
            else:
                results[mutant_pair][k] = results_curr[mutant_pair][k]

    print '\nThe total # of entries in results = {},  entries_sum = {} '.format(len(results[mutant_pair].keys()),entries_sum)

    # Write the results inot the specified output file
    with open(output_file_name,'w') as f:
        f.write(results_var_name + ' = {}\n')
        f.write(results_var_name + '[' + str(mutant_pair) + '] = {}\n')
        for k in results[mutant_pair].keys():
            f.write(results_var_name + '[' + str(mutant_pair) + '][' + str(k) + '] = ' + str(results[mutant_pair][k]) + '\n')
    print '\nResults were integrated and written into ',output_file_name,'\n'

def del_results_files(results_file_names):
    """
    Deletes the results files after verification

    INPUTS:
    -------
    results_file_names: A list of strings containing the names of the files containing the results
    """
    for file_name in results_file_names:
        print 'deleteing ',file_name, ' ...'
        os.system('rm -rf ' + file_name)

def run_integrate_results_files():
    """
    Rns integrate_results_files
    """
    # Test results
    run_test =False
    if run_test:
        results_filenames = ['results/test_coopr_level_652_652_1_2.py','results/test_coopr_level_652_652_3_4.py'] 
        integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), output_file_name = 'results/test_all.py') 

    # Experimental glucose
    run_exp = True
    if run_exp:
        results_filenames = ['results/coopr_level_expGlc_652_652_1_25.py', 'results/coopr_level_expGlc_652_652_26_50.py', 'results/coopr_level_expGlc_652_652_51_75.py', 'results/coopr_level_expGlc_652_652_76_100.py', 'results/coopr_level_expGlc_652_652_101_125.py', 'results/coopr_level_expGlc_652_652_126_150.py', 'results/coopr_level_expGlc_652_652_151_175.py', 'results/coopr_level_expGlc_652_652_176_200.py', 'results/coopr_level_expGlc_652_652_201_225.py', 'results/coopr_level_expGlc_652_652_226_250.py', 'results/coopr_level_expGlc_652_652_251_275.py', 'results/coopr_level_expGlc_652_652_276_300.py', 'results/coopr_level_expGlc_652_652_301_325.py', 'results/coopr_level_expGlc_652_652_326_350.py', 'results/coopr_level_expGlc_652_652_351_375.py', 'results/coopr_level_expGlc_652_652_376_400.py', 'results/coopr_level_expGlc_652_652_401_425.py', 'results/coopr_level_expGlc_652_652_426_441.py']
        integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), output_file_name = 'results/coopr_level_expGlc_652_652.py') 


    # Excess glucose
    run_excess = False
    if run_excess:
        results_filenames = ['results/coopr_level_excessGlc_652_652_1_25.py', 'results/coopr_level_excessGlc_652_652_26_50.py', 'results/coopr_level_excessGlc_652_652_51_75.py', 'results/coopr_level_excessGlc_652_652_76_100.py', 'results/coopr_level_excessGlc_652_652_101_125.py', 'results/coopr_level_excessGlc_652_652_126_150.py', 'results/coopr_level_excessGlc_652_652_151_175.py', 'results/coopr_level_excessGlc_652_652_176_200.py', 'results/coopr_level_excessGlc_652_652_201_225.py', 'results/coopr_level_excessGlc_652_652_226_250.py', 'results/coopr_level_excessGlc_652_652_251_275.py', 'results/coopr_level_excessGlc_652_652_276_300.py', 'results/coopr_level_excessGlc_652_652_301_325.py', 'results/coopr_level_excessGlc_652_652_326_350.py', 'results/coopr_level_excessGlc_652_652_351_375.py', 'results/coopr_level_excessGlc_652_652_376_400.py', 'results/coopr_level_excessGlc_652_652_401_425.py', 'results/coopr_level_excessGlc_652_652_426_441.py' ]
        integrate_results_files(results_filenames = results_filenames, results_var_name = 'results', mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), output_file_name = 'results/coopr_level_excessGlc_652_652.py') 


def plot_results(results_filename, x, y, plot_func = 'pcolor', mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), title = '', xaxis_label = '', xaxis_label_format = None, yaxis_label = '', yaxis_label_format = None, x_minorticks_spacing = 5, y_minorticks_spacing = 5, set_minor_xticks = True, colorbar_label = '', colorbar_label_format = None, output_filename = ''):
    """
    Plots a heatmap of the anaylsis results.
 
    INPUTS:
    -------
       results_filenames: A list of strings containing the names of the files containing the results
                       x: A list containing the elements of horizontal axis 
                       x: A list containing the elements of vertical axis 
             mutant_pair: The mutant pair for which the results are going to be plotted
             xaxis_label: x-axis label
             yaxis_label: y-axis label
        set_minor_xticks: Whethter to set minor xtick_labels
    x_minorticks_spacing: Spacing between minor xtick_labels
        output_fiilename: Name of the output file containing the plot
    """
    import os
    from imp import load_source
    import numpy as np
    from matplotlib import colors
    from tools.ancillary.plot_heatmap import plot_heatmap

    #--- Load the data ---
    # Import the data in the module stored in results_filename
    if type(results_filename) == str:
        if not os.path.isfile(results_filename):
            raise IOError("No such file was found :'" + results_filename + "'")
        else:
            # First delete the model dataFile if it already exists. If it is not deleted
            # the new module data is merged with the previous ones
            try:
                del sys.modules['dataFile']
            except:
                pass
            load_source('dataFile',results_filename)
            import dataFile

    results = dataFile.results

    data = np.zeros((len(y),len(x)))

    results_summary = {}

    for i, yy in enumerate(y):
        for j, xx in enumerate(x):
            #res_key = (('coopr_level_m1',xx),('coopr_level_m2',yy))
            res_key = ((mutant_pair[0] + '_coopr_level',xx),(mutant_pair[1] + '_coopr_level',yy))
            if res_key not in results[mutant_pair].keys():
                raise userError('{} not a key in the results[{}]!'.format(res_key,mutant_pair))
            else:
                # Find the fold change in total community concentration at the final tiime point (tf) 
                # compared to that at the begining of th eexperiment (t = 0)
                time_points = sorted(results[mutant_pair][res_key]['cell_concs'][mutant_pair[0]].keys())
                tf = max(time_points)
                conc_total_t0 = results[mutant_pair][res_key]['cell_concs'][mutant_pair[0]][0] + results[mutant_pair][res_key]['cell_concs'][mutant_pair[1]][0]
                conc_total_tf = results[mutant_pair][res_key]['cell_concs'][mutant_pair[0]][tf] + results[mutant_pair][res_key]['cell_concs'][mutant_pair[1]][tf]
                fold_change = conc_total_tf/conc_total_t0
               
                results_summary[res_key] = fold_change

            data[i,j] = fold_change 

    print '\nmin data = {}   ,  max data = {}\n'.format(data.min(),data.max())

    import operator
    print '\n'
    for k in sorted(results_summary.items(), key=operator.itemgetter(1)): 
        print k
    print '\n'

    if plot_func in ['matshow','imshow']:
        interpolate = True
        invert_yaxis = True
        xticklabels_format = {'rotation':90}
    else:
        interpolate = False
        invert_yaxis = False
        xticklabels_format = {'rotation':0}

    
    plot_heatmap(x = np.array(x),y = np.array(y),data = data, plot_func = plot_func, interpolate = interpolate, title = title, xaxis_label = xaxis_label, xaxis_label_format = xaxis_label_format, yaxis_label = yaxis_label, yaxis_label_format = yaxis_label_format, set_minor_xticks = set_minor_xticks, set_minor_yticks = True, x_majorticks_spacing = None, x_minorticks_spacing = x_minorticks_spacing,y_majorticks_spacing = None, y_minorticks_spacing = y_minorticks_spacing, xticklabels_format = xticklabels_format, invert_xaxis = False, invert_yaxis = invert_yaxis, colorbar = None, colorbar_ticklabels = None, colorbar_label = colorbar_label, colorbar_label_format = colorbar_label_format,grid = False, figsize = (25,15), dpi = None, output_filename = output_filename)


def run_plot_results():
    """
    Plots the results
    """
    run_test = False
    run_exp = False
    run_excess = True

    if run_test:
        coopr_level_m1 = [10,20] 
        coopr_level_m2 = [10,20] 
        plot_func = 'pcolor'
        plot_results(results_filename = 'results/test_all.py', x = coopr_level_m1, y = coopr_level_m2, plot_func = plot_func, mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), title = 'Experimentql glucose uptake', xaxis_label = '% of max Isoleucine production level by lysA mutant', yaxis_label = '% of max lysine production level by ilvE mutant', x_minorticks_spacing = 10, y_minorticks_spacing = 10, set_minor_xticks = False, colorbar_label = 'Fold increase in total cell concentration', colorbar_label_format = {'distance_from_ticklabels':50}, output_filename = 'results/test.pdf')

    if run_exp:
        coopr_level_m1 = range(0,101,5)
        coopr_level_m2 = range(0,101,5)
        plot_func = 'pcolor'
        plot_results(results_filename = 'results/coopr_level_expGlc_652_652.py', x = coopr_level_m1, y = coopr_level_m2, plot_func = plot_func, mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), title = 'Experimentql glucose uptake', xaxis_label = '% of max Isoleucine production level by lysA mutant', xaxis_label_format = {'size':35}, yaxis_label = '% of max lysine production level by ilvE mutant', yaxis_label_format = {'size':35}, x_minorticks_spacing = 10, y_minorticks_spacing = 10, set_minor_xticks = False, colorbar_label = 'Fold increase in total cell concentration', colorbar_label_format = {'size':35,'distance_from_ticklabels':50}, output_filename = 'results/coopr_level_expGlc.pdf')

    if run_excess:
        coopr_level_m1 = range(0,101,5)
        coopr_level_m2 = range(0,101,5)
        plot_func = 'pcolor'
        plot_results(results_filename = 'results/coopr_level_excessGlc_652_652.py', x = coopr_level_m1, y = coopr_level_m2, plot_func = plot_func, mutant_pair = ('lysA_Ecoli', 'ilvE_Ecoli'), title = '', xaxis_label = '% of max Isoleucine production level by lysA mutant', xaxis_label_format = {'font_size':30}, yaxis_label = '% of max lysine production level by ilvE mutant', yaxis_label_format = {'font_size':30},x_minorticks_spacing = 10, y_minorticks_spacing = 10, set_minor_xticks = False, colorbar_label = 'Fold increase in total cell concentration', colorbar_label_format = {'font_size':30,'distance_from_ticklabels':30}, output_filename = 'results/coopr_level_excessGlc.pdf')

#-------------------------
if __name__ == '__main__':
    pass


