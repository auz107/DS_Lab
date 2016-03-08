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
from tools.ancillary.get_ModelSEED_ids import get_cmp_ModelSEED_id
from tools.ancillary.cell_pool_conc import cell_pool_conc
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
from DMMM_mortality import DMMM_mortality
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.ancillary.importData import importData
import cobra
from read_exp_data import read_exp_data
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)
from multiprocessing import Process, Manager
import re

# Last updated: 01-14-2015

def performDMMM(input_data):

    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings']
    t0 = input_data['t0']
    tf = input_data['tf']
    delta_t = input_data['delta_t']
    time_points = input_data['time_points'][:]
    model_path = input_data['model_path'] 
    mutant1 = deepcopy(input_data['mutant1'])
    mutant2 = deepcopy(input_data['mutant2'])
    shared_cmp_ids  = input_data['shared_cmp_ids'][:]
    shared_cmp_names  = dict(input_data['shared_cmp_names'][:])
    auxoMetabsMutants_m1  = input_data['auxoMetabsMutants_m1'][:]
    auxoMetabsMutants_m2  = input_data['auxoMetabsMutants_m2'][:]
    pc_cell_conc_init = input_data['pc_cell_conc_init'] 
    rand_mort_percent = input_data['rand_mort_percent'] 
    pool_conc_fac = input_data['pool_conc_fac'] 
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
        uptake_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
        uptake_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 
        if exch_rxn_id not in auxoMetabsMutants_m1 and exch_rxn_id != 'EX_glc(e)':
            export_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
        if exch_rxn_id not in auxoMetabsMutants_m2 and exch_rxn_id != 'EX_glc(e)':
            export_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 

    # Create the list of shared compounds
    shared_cmps = []

    # Glucose as a shared compound (concentration in mM)
    glucose = compound(id = 'glc_D_e', name = 'D-Glucose', KEGG_id = 'C00031', reactant_reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],concentration = {0:111.01})    

    shared_cmps.append(glucose)

    for shared_cmp_id in [id for id in shared_cmp_ids if id != 'glc_D_e']:
        reactant_rxns = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
        product_rxns = [r for r in export_rxns_m1 + export_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
        shared_cmp = compound (id = shared_cmp_id, name = shared_cmp_names[shared_cmp_id],reactant_reactions = reactant_rxns, product_reactions = product_rxns, concentration = {0:0}) 
        shared_cmps.append(shared_cmp)

    # Set store_flux to True for all exchange reacitons related to shared metabolites
    for rxn in [r for c in shared_cmps for r in c.reactions]:
        rxn.store_flux = True
 
    for rxn in [m.biomass_reaction for m in [mutant1,mutant2]]: 
        rxn.store_flux = True

    # Get the ModelSEED ids for shared compounds
    get_cmp_ModelSEED_id(cmp_list = shared_cmps,compart_list = [c.id for c in mutant1.compartments],stdout_msgs = False)

    # Get the concentration of the shared compounds in a cell pool
    cell_pool_conc(cmp_org = dict([(c,[mutant1.organism,mutant2.organism]) for c in shared_cmps]) ,stdout_msgs = False)

    set_specific_bounds(model = mutant1, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py', flux_bounds = dict([(rid,[0,0]) for rid in  mutant1.knockedout_rxn_ids]))
    set_specific_bounds(model = mutant2, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py', flux_bounds = dict([(rid,[0,0]) for rid in  mutant2.knockedout_rxn_ids]))

    if stdout_msgs:
        print '\nshared metabs ids = ',shared_cmp_ids,'\n'

    DMMM_m1m2 = DMMM_mortality(community_members = [mutant1,mutant2],shared_compounds = shared_cmps, time_range = [t0,delta_t,tf], random_mortality_percentage = rand_mort_percent, cell_pool_factor = pool_conc_fac, store_dynamic_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)
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
                save_key = (('pc_cell_conc_init',pc_cell_conc_init), ('rand_mort_percent', rand_mort_percent),('pool_conc_fac', pool_conc_fac))
                f.write('results[' + str(mutants_names) + '][' + str(save_key) + '] = ' + str(results) + '\n')
            else:
                f.write('results[' + str(mutants_names) + '] = ' + str(results) + '\n')

#--------- Extract the compound names from the SBML file -----------------
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

#---------- DMMM with mortality for all mutant combinations ---------------
def master_func(t0,delta_t,tf,start_pos,end_pos, per_capita_cell_conc_init = 7.5e6, random_mortality_percentage = 1, cell_pool_conc_factor = 1, save_details = True, save_with_key = False, results_filename_base = '', stdout_msgs = True, warnings = True):
    """
     INPUTS:
     -------
                             t0: Initial simulaiton time
                        delta_t: Time step
                             tf: Total simulation time
                      strat_pos: Start position of the array containing all possible
                                 strategy combinations (see variable 'combinations').
                        end_pos: End position of the array containing all possible
                                 strategy combinations (see variable 'combinations').
          results_filename_base: The base for the file name storing the results
                                 Example: results_filename_base = 'results/emc_results'. The code
                                 will add the start and end positions to the file name.
                                 Example: 'results/emc_results_1_500.txt'
      per_capita_cell_conc_init: Initial per captial cell concentration
    random_mortality_percentage: The percentage of the cell dying  
          cell_pool_conc_factor: The factor that should be multiplied by the total pool of the compound concentration
                                 becoming available due to the cell death
                   save_details: A parameter showing whether the details of DMMM simulations (such as exchange flux of
                                 shared compounds) should be saved (True) or not (False)
                  save_with_key: A parameter showing whether the results have to be saved with keys 
                                 (('mortality_percentage',A),('pool_conc_factor',B)), where A and B are some numbers showing
                                 the specified mortality percetage and cell_pool_conc_factors respectively.

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
    
    if isinstance(per_capita_cell_conc_init,int) or isinstance(per_capita_cell_conc_init,float):
        per_capita_cell_conc_init = [per_capita_cell_conc_init]
    elif not isinstance(per_capita_cell_conc_init,list):
        raise TypeError('per_capita_cell_conc_init must be an integer, float, or a list of integers and floats') 

    if isinstance(random_mortality_percentage,int) or isinstance(random_mortality_percentage,float):
        random_mortality_percentage = [random_mortality_percentage]
    elif not isinstance(random_mortality_percentage,list):
        raise TypeError('random_mortality_percentage must be an integer, float, or a list of integers and floats') 

    if isinstance(cell_pool_conc_factor,int) or isinstance(cell_pool_conc_factor,float):
        cell_pool_conc_factor = [cell_pool_conc_factor]
    elif not isinstance(cell_pool_conc_factor,list):
        raise TypeError('cell_pool_conc_factor must be an integer, float, or a list of integers and floats') 

    if results_filename_base != '':
        results_filename = results_filename_base + '_' + str(start_pos) + '_' + str(end_pos) + '.py'
    else:
        results_filename = '' 

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]    
    print '\ntime_points = ',time_points

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655',gDW_per_cell = 2.8e-13)

    WT = read_sbml_model(file_name = model_path + 'iJO1366_updated.xml', model_id = 'iJO1366',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    WT.biomass_reaction = WT.reactions_by_id['Ec_biomass_iJO1366_core_53p95M']
    WT.all_biomass_reactions = {'core':WT.reactions_by_id['Ec_biomass_iJO1366_core_53p95M'],'WT':WT.reactions_by_id['Ec_biomass_iJO1366_WT_53p95M']}
    WT.organism.total_death_rate = {0:0}

    # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
    # Example: EX_glc(E): glc_D_e <==>    Vmax*C['glc_D_e']/(Km + C['glc_D_e']) 
    # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
    for reaction in [r for r in WT.reactions if r.type.lower() == 'exchange']:
        # The id of compound participating in the exchange reaction
        metab_id = [m.id for m in reaction.compounds][0]
        reaction.kinetics = "10*C['" + metab_id + "']/(10 + C['" + metab_id + "'])"

    # Glucose uptake kinetics 
    exch_rxns = WT.get_reactions({'EX_glc(e)':'id','EX_lys_L(e)':'id','EX_ile_L(e)':'id'})
    exch_rxns['EX_glc(e)'].kinetics = "10*C['glc_D_e']/(0.15 + C['glc_D_e'])"
    exch_rxns['EX_lys_L(e)'].kinetics = "0.1964*C['lys_L_e']/(5e-4 + C['lys_L_e']) + 0.3055*C['lys_L_e']/(1e-2 + C['lys_L_e'])"
    exch_rxns['EX_ile_L(e)'].kinetics = "0.0346*C['ile_L_e']/(1.22e-3 + C['ile_L_e'])"
   
    # Growth medium
    set_specific_bounds(model = WT, file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',flux_bounds = {'EX_glc(e)':[-10,1000]})

    # Assign the objective function coefficients
    for rxn in WT.reactions:
        rxn.objective_coefficient = 0
    WT.biomass_reaction.objective_coefficient = 1

    # Perform FBA for the wild-type
    WT.fba(assign_wildType_max_biomass = True, stdout_msgs = stdout_msgs)

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

    # All possible pair combinations (do not consider the ones that cannot be rescued)
    mutant_pairs = list(itertools.combinations([m for m in mutants_rxn_info.keys() if m not in not_rescued_mutants],r=2))   
    print '\nThe total # of mutant pairs to examine = %i' % len(mutant_pairs)

    # IDs of the shared compounds
    shared_cmp_ids = ['glc_D_e']

    # The id of all compounds in the cell pool
    cell_pool_cmp_ids = ['ala_L_e','val_L_e','gly_e','ile_L_e','thr_L_e','leu_L_e','ser_L_e','pro_L_e','asp_L_e','cys_L_e','met_L_e','glu_L_e','phe_L_e','tyr_L_e','orn_e','lys_L_e','trp_L_e','arg_L_e']

    shared_cmp_ids += cell_pool_cmp_ids

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
        mutant1.knockedout_rxn_ids = mutants_rxn_info[m1]
        for rxn_id in mutants_rxn_info[m1]:
            rxn = mutant1.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]

        mutant1.fba(build_new_optModel = False, store_opt_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)

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
        mutant2.knockedout_rxn_ids = mutants_rxn_info[m2]
        for rxn_id in mutants_rxn_info[m2]:
            rxn = mutant2.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]
            
        mutant2.fba(build_new_optModel = False, store_opt_fluxes = False, stdout_msgs = stdout_msgs, warnings = warnings)
    
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
        shared_cmp_names['tlc-D_e'] = 'D-glucose'
        shared_cmp_names['ala_L_e'] = 'L-Alanine'
        shared_cmp_names['val_L_e'] = 'L-valine'
        shared_cmp_names['gly_e'] = 'glycine'
        shared_cmp_names['ile_L_e'] = 'L-isoleucine'
        shared_cmp_names['thr_L_e'] = 'L-threonine'
        shared_cmp_names['leu_L_e'] = 'L-leucine'
        shared_cmp_names['ser_L_e'] = 'L-serine'
        shared_cmp_names['pro_L_e'] = 'L-proline'
        shared_cmp_names['asp_L_e'] = 'L-aspartate'
        shared_cmp_names['cys_L_e'] = 'L-cysteine'
        shared_cmp_names['met_L_e'] = 'L-methionine'
        shared_cmp_names['glu_L_e'] = 'L-glutamate'
        shared_cmp_names['phe_L_e'] = 'L-phenylalanine'
        shared_cmp_names['tyr_L_e'] = 'L-tyrosine'
        shared_cmp_names['orn_e'] = 'Ornithine'
        shared_cmp_names['lys_L_e'] = 'L-lysine'
        shared_cmp_names['trp_L_e'] = 'L-tryptophan'
        shared_cmp_names['arg_L_e'] = 'L-arginine'

        # Simulation scenarios for each mutant
        print '\nrandom_mortality_percentage = ',random_mortality_percentage
        print '\ncell_pool_conc_factor = ',cell_pool_conc_factor

        all_cases = [(pc_cell_conc_init, rand_mort_percent,pool_conc_fac) for pc_cell_conc_init in per_capita_cell_conc_init for rand_mort_percent in random_mortality_percentage for pool_conc_fac in cell_pool_conc_factor]
        print '\nThe total # of cases to consider for mutant pair ({},{}) = {}\n'.format(m1,m2,len(all_cases))
        
        counter = 0

        for (pc_cell_conc_init, rand_mort_percent,pool_conc_fac) in all_cases:

            counter += 1
          
            print '{}. (pc_cell_conc_init, rand_mort_percent,pool_conc_fac) = {}'.format(counter,(pc_cell_conc_init, rand_mort_percent,pool_conc_fac))

            mutant1.organism.cells_per_ml = {0:pc_cell_conc_init}
            mutant1.organism.gDW_per_ml = {0:pc_cell_conc_init*mutant1.organism.gDW_per_cell}
            mutant2.organism.cells_per_ml = {0:pc_cell_conc_init}
            mutant2.organism.gDW_per_ml = {0:pc_cell_conc_init*mutant2.organism.gDW_per_cell}

            # Creating a shared memory using the manager
            input_data = {} 
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
            input_data['pc_cell_conc_init'] = pc_cell_conc_init 
            input_data['rand_mort_percent'] = rand_mort_percent 
            input_data['pool_conc_fac'] = pool_conc_fac 
            input_data['save_with_key'] = save_with_key 
            input_data['save_details'] = save_details 
            input_data['results_filename'] = results_filename
            input_data['stdout_msgs'] = stdout_msgs
            input_data['warnings'] = warnings
    
            p = Process(target = performDMMM, args = (input_data,))
            p.start()
            p.join() 
            if p.exitcode > 0:
                raise userError('Error in python subprocess. Please check performDMMM\n')

def inspect_init_conc():
    """
    This function investigates the effect of initial cell concentration on
    whether the free pool of amino acids released due to the cell death can 
    serve as a driver of growth and dominates the cell death itself
    The simulations are performed for ('lysA', 'ilvE') only
    """
    # Create a vector of initial concentrations
    per_capita_cell_conc_init_vec = [7.5e6,7.5e7,7.5e8,7.5e9,7.5e10]
    master_func(t0 = 0,delta_t = 0.5,tf = 96,start_pos = 652,end_pos = 652, per_capita_cell_conc_init = per_capita_cell_conc_init_vec, save_details = False, save_with_key = True, results_filename_base = 'results/mortality_init_conc_effect', stdout_msgs = False)

def process_init_conc_results():
    """
    Processes the results of inspect_init_conc    
    """
    from DMMM_mortality_init_conc_all import results
    
    # Time points 
    time_points = sorted(results[0]['cell_concs'][results[0]['cell_concs'].keys()[0]].keys())
    tf = time_points[-1] 

    # List of mutants
    mutants_list = results[0]['cell_concs'].keys()

    total_cell_concs = {}

    for res in results:
        # Create a MATLAB matrix where each row represents a 
        # Create a dicitonary where keys are the total initial cell concentrations 
        # and values are the time course total cell concentrations 

        # Total initial cell concentration
        init_cell_conc = sum([res['cell_concs'][m][0] for m in mutants_list])           
        total_cell_concs[init_cell_conc] = [sum([res['cell_concs'][m][t] for m in mutants_list]) for t in time_points]

    # Write the results into a MATLAB file
    with open('results/mortality_init_conc_results_toPlot_withMatlab.m','w') as f:
        # Write the time points
        f.write('t = [\n')
        for t in time_points:
            f.write(str(t) + '\n')
        f.write('];\n')

        # Write the total concentrations for each case as a row of a matrix
        f.write('\ntotal_concs_init_effect = [\n')
        for init_conc in total_cell_concs.keys():
            for conc in total_cell_concs[init_conc]:
                f.write(str(conc) + '\t')
            f.write('\n')
        f.write('];\n')

def inspect_mortality_percent():
    """
    This function investigates increasing the mortality rate (intrinsic death rate) on
    can serve as a driver of growth through the release of amino acid pool
    The simulations are performed for ('lysA', 'ilvE') only
    """
    # Create a vector of initial concentrations
    random_mortality_percentage_vec = [1] + range(5,100,5) 
    master_func(t0 = 0,delta_t = 0.5,tf = 96,start_pos = 652,end_pos = 652,random_mortality_percentage = random_mortality_percentage_vec, save_details = False, save_with_key = True, results_filename_base = 'results/mortality_percentage_effect', stdout_msgs = False)


def plot_mortalityResults(results_filename, results_key, start_pos = None, end_pos = None, title = '', xaxis_label_dynamic = 'time (h)', plot_labels_guide_dynamic = '', xaxis_label_tf = '% of cells dying', yaxis_label = 'log10(total cell conc) [cells/ml]', output_filename_base = ''):
    """
    Analyzes and plots the results of 
   
    INPUTS:
    -------
                  results_key: Name of the key in the results variable containing the parameter of interest (rand_mort_percent or ) 
    plot_labels_guide_dynamic: The explanation of the plot labels in the dynamic graph
    """
    import os
    from imp import load_source
    import matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize = (25,15))
    #ax.set_aspect('auto')

    # Update the matplotlib configuration parameters:
    matplotlib.rcParams.update({'font.size': 30,'font.weight':'bold','font.family': 'STIXG'})

    # Set the thcknes of axes 
    ax.spines['top'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    ax.spines['right'].set_linewidth(4)
    ax.spines['left'].set_linewidth(4)

    # Load the data from the files
    load_source('dataFile',results_filename)
    import dataFile
    results = dataFile.results

    # A dictionary storing the total cell concetration at the final time points 
    # for varying values of mortality percentage. 
    total_cell_conc_tf = {}

    # Specifies the positino of the lable
    t_pos = 4

    counter = 0

    if start_pos != None and end_pos != None:
        cases_to_consider = sorted(results[('lysA_Ecoli', 'ilvE_Ecoli')].keys(),key=lambda x:dict(list(x))[results_key])[start_pos - 1:end_pos]
    else:
        cases_to_consider = sorted(results[('lysA_Ecoli', 'ilvE_Ecoli')].keys(),key=lambda x:dict(list(x))[results_key])

    for k in cases_to_consider: 

        counter += 1

        # Key parameter in results (Mortality rate or fold increase in the release compounds into the medium due to death)
        results_key_param = dict(list(k))[results_key] 

        # organisms names
        org_names = ('lysA_Ecoli', 'ilvE_Ecoli') 

        # Time points
        time_points = sorted(results[('lysA_Ecoli', 'ilvE_Ecoli')][k]['cell_concs'][org_names[0]].keys()) 

        # Final time point
        tf = max(time_points)

        # Total cell concentration (conc mutant 1 + conc mutant 2)
        total_cell_conc = dict([(t,np.log10(sum([results[('lysA_Ecoli', 'ilvE_Ecoli')][k]['cell_concs'][org_name][t] for org_name in org_names]))) for t in time_points])   

        total_cell_conc_tf[results_key_param] = total_cell_conc[tf]
        total_cell_conc_init = total_cell_conc[0]
        print '{}. {}   {}'.format(counter,results_key_param,total_cell_conc_tf[results_key_param])

        # Plot the results
        ax.plot(sorted(total_cell_conc.keys()), [total_cell_conc[kk] for kk in sorted(total_cell_conc.keys())], linewidth = 4)
        if results_key == 'rand_mort_percent':
            t_pos += 4
            if t_pos > 96:
               t_pos = 4
            ax.text(t_pos,1e-8 + total_cell_conc[t_pos],str(results_key_param) + '%', fontsize = 25,fontweight = 'bold')
            if counter == 1:
                ax.text(t_pos - 5, 1.2 + total_cell_conc[t_pos],str(plot_labels_guide_dynamic), fontsize = 30,fontweight = 'bold')
        else:
            ax.text(96-5,1e-8 + total_cell_conc[96-5],str(results_key_param), fontsize = 25,fontweight = 'bold')
            if counter == 1:
                ax.text(96 - 27, -1.1 + total_cell_conc[0],str(plot_labels_guide_dynamic), fontsize = 30,fontweight = 'bold')


    ax.set_title(title,{'weight':'bold','size':35})
    ax.set_xlabel(xaxis_label_dynamic,{'weight':'bold','size':30})
    ax.set_ylabel(yaxis_label,{'weight':'bold','size':30})

    #ax.legend(loc=0)

    for ticklabel in ax.get_xmajorticklabels():
        ticklabel.set_fontsize(30)
        ticklabel.set_fontweight('bold')
    for ticklabel in ax.get_ymajorticklabels():
        ticklabel.set_fontsize(30)
        ticklabel.set_fontweight('bold')

    ax.set_xlim([0,96]) 
    ax.set_xticks([i for i in range(0,96+1,8)])

    ax.grid(True, color = 'k', linestyle='dashed', linewidth = 2)

    fig.savefig(output_filename_base + '_dynamic.pdf')

    plt.tight_layout()
    plt.show()  

    #---- Now plot the final concentrations ----
    fig, ax = plt.subplots(figsize = (25,15))

    # Update the matplotlib configuration parameters:
    matplotlib.rcParams.update({'font.size': 30,'font.weight':'bold','font.family': 'STIXG'})

    # Set the thcknes of axes 
    ax.spines['top'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    ax.spines['right'].set_linewidth(4)
    ax.spines['left'].set_linewidth(4)

    param_vec = sorted(total_cell_conc_tf.keys())
    ax.plot(param_vec, [total_cell_conc_tf[p] for p in param_vec], linewidth = 4, marker = 'o', markersize = 15, markerfacecolor = 'white', markeredgewidth = 4)
    ax.plot(param_vec, [total_cell_conc_init for p in param_vec], linewidth = 6, linestyle = '--')
    if results_key == 'rand_mort_percent':
        ax.text(100 - 45,0.5 + total_cell_conc_init,'Total initial cell concentration', fontsize = 30,fontweight = 'bold')
    else:
        ax.text(1000 - 400,0.1 + total_cell_conc_init,'Total initial cell concentration', fontsize = 30,fontweight = 'bold')

    print 'total_cell_conc_init = ', total_cell_conc_init
 
    ax.set_title(title,{'weight':'bold','size':35})
    ax.set_xlabel(xaxis_label_tf,{'weight':'bold','size':30})
    ax.set_ylabel(yaxis_label,{'weight':'bold','size':30})

    for ticklabel in ax.get_xmajorticklabels():
        ticklabel.set_fontsize(30)
        ticklabel.set_fontweight('bold')
    for ticklabel in ax.get_ymajorticklabels():
        ticklabel.set_fontsize(30)
        ticklabel.set_fontweight('bold')
 
    ax.set_xlim([min(param_vec),max(param_vec)]) 
    if results_key == 'rand_mort_percent':
        ax.set_xticks([i for i in range(5,max(param_vec)+1,5)])
    else:
        ax.set_xticks([1] + range(100,1001,100))

    ax.grid(True, color = 'k', linestyle='dashed', linewidth = 2)

    fig.savefig(output_filename_base + '_tf.pdf')

    plt.tight_layout()
    plt.show()  

def run_plot_mortalityResults():
    """
    Plots the mortality results
    """
    # Plot the results of inspect_mortality_percent()
    plot_mortalityResults(results_filename = 'results/mortality_percentage_effect_652_652_v2.py', results_key = 'rand_mort_percent', title = 'Impact of the cell death level', xaxis_label_dynamic = 'time (h)', xaxis_label_tf = '% of cells dying', yaxis_label = 'log10(total cell conc) [cells/ml]', plot_labels_guide_dynamic = '% of cells dying', output_filename_base = 'results/mortality_percentage')

    # Plot the results of inspect_cell_pool_conc() 
    #plot_mortalityResults(results_filename = 'results/cell_pool_conc_effect_652_652_v2.py', start_pos = 1, end_pos = 14, results_key = 'pool_conc_fac', title = 'Impact of the level of spilled amino acids due to cell death', xaxis_label_dynamic = 'time (h)', plot_labels_guide_dynamic = 'Fold increase in the level \nof released compounds', xaxis_label_tf = 'Fold increase in the level of released compounds', yaxis_label = 'log10(total cell conc) [cells/ml]', output_filename_base = 'results/cell_pool_conc_effect')
 
def inspect_cell_pool_conc():
    """
    This function investigates how increasing the concentration of the compounds 
    released to th environment due to cell death affects the community growth 
    The simulations are performed for ('lysA', 'ilvE') only
    """
    # Create a vector of the factor multiplied by the cell pool conc multiplied by each compound 
    factor_vec = [1] + range(50,1001,50) 
    master_func(t0 = 0,delta_t = 0.5,tf = 96,start_pos = 652,end_pos = 652, cell_pool_conc_factor = factor_vec, save_details = False, save_with_key = True, results_filename_base = 'results/cell_pool_conc_effect', stdout_msgs = False)

def process_results(metab_thr = 1e-10,total_cell_conc_frac = 0.9,exp_fold_growth_thr = 50):
    """
    Process the results

    INPUTS:
    ------
    metab_thr: Threshold above which metabolite concentration is cosidered to be significant
    total_cell_conc_thr: Percentabge of max fold change in total concetration of cells above which
                         is considered as cooperation 
    exp_fold_growth_thr: Threshold for cooperation accroding to experimental data

    """
    from DMMM_mortality_all import results
    
    print '\nThe total # of pairs = %s\n' % len(results)

    # List of all mutants in results. This excludes three mutants, which are not essential in silico
    # and four mutants which cannot be rescued in silico
    mutants_in_results = list(set([m1m2 for res in results for m1m2 in res['name']]))

    #---- Report based on compound concentration ---    
    # Final time point 
    tf = sorted(results[0]['cell_concs'][results[0]['cell_concs'].keys()[0]].keys())[-1]

    pairs_with_nonzero_shared_cmp = [p for p in results if len([m for m in p['shared_cmp_concs'].keys() if max(p['shared_cmp_concs'][m].values()) > metab_thr and m != 'glc_D_e']) > 0]
    print '\nThe total # of pairs with non-zero shared compounds = %i\n'% len(pairs_with_nonzero_shared_cmp)

    for pair in pairs_with_nonzero_shared_cmp:
        print '\n',pair['name'],'     ',
        for shared_cmp in [m for m in pair['shared_cmp_concs'].keys() if max(pair['shared_cmp_concs'][m].values()) > metab_thr and m != 'glc_D_e']: 
            print '(%s: %.9f) ' %(shared_cmp,max(pair['shared_cmp_concs'][shared_cmp].values())),

    print '\nThe total # of pairs with non-zero shared compounds = %i\n'% len(pairs_with_nonzero_shared_cmp)
    print 

    #---- Report based on single cell concentration ---   
    # Initial cell concentration
    per_capita_cell_conc_init = 7.5e6

    # Create a dictionary where keys and values as follows
    #   Keys: A tuple where the first element is the name of the mutant and the
    #         second element is another tuple containing the name of the mutant pair
    # Values: The cell conc at the final time point
    cell_concs_tf = dict([((m,tuple(res['cell_concs'].keys())),res['cell_concs'][m][tf]) for res in results for m in res['cell_concs'].keys()])

    # Find the mutant with maximum increase in the cell concentration
    max_conc_mutant = max(cell_concs_tf.iteritems(),key = lambda x: x[-1])
    max_fold_change_mutant = max_conc_mutant[1]/per_capita_cell_conc_init
    print '\nmutant %s in mutant pair %s has the maximum fold increase %.6f\n'%(max_conc_mutant[0][0],max_conc_mutant[0][1],max_fold_change_mutant)

    # Find the mutant with minimum increase in the cell concentration
    min_conc_mutant = min(cell_concs_tf.iteritems(),key = lambda x: x[-1])
    min_fold_change_mutant = min_conc_mutant[1]/per_capita_cell_conc_init
    print '\nmutant %s in mutant pair %s has the minimum fold increase %.6f\n'%(min_conc_mutant[0][0],min_conc_mutant[0][1],min_fold_change_mutant)

    #---- Report based on total cell concentration ---   
    # Create a dictionary where keys and values as follows
    #   Keys: A tuple containing the name of the mutant pair 
    # Values: The total conc of mutants at the final time point
    pair_concs_tf = dict([(tuple(res['cell_concs'].keys()),sum([res['cell_concs'][m][tf] for m in res['cell_concs'].keys()])) for res in results])

    # Write into a fiile to plot with matlab 
    with open('results/mortality_results_toPlot_withMatlab.m','w') as outfile:
        outfile.write('mortality_fold = [\n')
        for v in pair_concs_tf.values():
            outfile.write(str(v/(2*per_capita_cell_conc_init)) + '\n')
        outfile.write('];\n')

    # Find the mutant pairs with maximum and minimum total concentration
    max_total_conc_pair = max(pair_concs_tf.iteritems(),key = lambda x:x[1])
    max_fold_change_pair_conc = max_total_conc_pair[1]/(2*per_capita_cell_conc_init)
    print '\nMutant pair %s has the maximum fold increase in total concentration of  %.6f\n'%(max_total_conc_pair[0],max_fold_change_pair_conc)

    min_total_conc_pair = min(pair_concs_tf.iteritems(),key = lambda x:x[1])
    min_fold_change_pair_conc = min_total_conc_pair[1]/(2*per_capita_cell_conc_init)
    print '\nMutant pair %s has the minimum fold increase in total concentration of  %.6f\n'%(min_total_conc_pair[0],min_fold_change_pair_conc)

    total_cell_conc_thr = total_cell_conc_frac*max_fold_change_pair_conc

    # Find all mutant pairs whose total concentrations is above the threshold
    above_thr_pairs = dict([p for p in pair_concs_tf.iteritems() if p[1]/(2*per_capita_cell_conc_init) > total_cell_conc_thr])

    print '\nThe following pairs have a total concentration above the threshould:\n'
    for p in sorted(above_thr_pairs.iteritems(),key = lambda x:x[1], reverse = True):
        print '%s\t%s: %.6f\t%s: %.6f\ttotal: %.6f'%(p[0],p[0][0],cell_concs_tf[(p[0][0],p[0])]/per_capita_cell_conc_init,p[0][1],cell_concs_tf[(p[0][1],p[0])]/per_capita_cell_conc_init,p[1]/(2*per_capita_cell_conc_init))

    print '\nThe total # of pairs whose total concentration is above the threshould (%.6f) = %i\n'%(total_cell_conc_thr,len(above_thr_pairs))

    #--- Compare with experimental data -----
    # Create a dictionary matching the pair names with and without '_Ecoli'
    no_Ecoli_name_map = []
    for p in pair_concs_tf.keys():
        no_Ecoli_name_map.append((p,(re.sub('_Ecoli','',p[0]),re.sub('_Ecoli','',p[1]))))
    no_Ecoli_name_map = dict(no_Ecoli_name_map)

    # Load the experimental growth data
    [day1Ave,day2Ave,day3Ave,day4Ave] = read_exp_data() 

    # Experimental fold growth for the community. Consider only mutant pairs analyzed
    expFoldGrowth = dict([((k[0]+'_Ecoli',k[1]+'_Ecoli'),day4Ave[k]/(2*per_capita_cell_conc_init)) for k in day4Ave.keys() if k[0]+'_Ecoli' in mutants_in_results and k[1]+'_Ecoli' in mutants_in_results])

    # Write the results into a text fiel be plotted with matlab
    with open('results/mortality_results_toPlot_withMatlab.m','a') as outfile:
        outfile.write('\nexp_fold = [\n')
        for v in expFoldGrowth.values():
            outfile.write(str(v) + '\n')
        outfile.write('];\n')

    # Max fold growth for experimental data
    max_fold_growth_exp = max(expFoldGrowth.iteritems(),key = lambda x:x[1])
    print '\nThe max fold growth based experimental data is for %s and is equal to %.6f\n'%(max_fold_growth_exp[0],max_fold_growth_exp[1])

    # Mutant pairs whose fold growth (in total cell concentration) is above eight
    above_thr_pair_exp = dict([p for p in expFoldGrowth.iteritems() if p[1] >= exp_fold_growth_thr])   

    print '\nTotal # of pair with a fold growth above the threshold for experimental data: %i\n'%(len(above_thr_pair_exp.keys()))
    
    # Find how many of the predicted cooperative mutants cooperated accroding to experimental data
    mortality_exp_matches = list(set.intersection(set(above_thr_pairs.keys()),set(above_thr_pair_exp.keys())))

    print '\nMatches between cooperative phenotypes both in experimental data and mortality simulaitons:\n'
    for p in mortality_exp_matches:
        print '%s\tmortality %s: %.6f\tmortality %s: %.6f\tmortality total: %.6f\texp total:%.6f'%(p,p[0],cell_concs_tf[(p[0],p)]/per_capita_cell_conc_init,p[1],cell_concs_tf[(p[1],p)]/per_capita_cell_conc_init,above_thr_pairs[p]/(2*per_capita_cell_conc_init),above_thr_pair_exp[p])
    print '\nThe total # of cooeprative mataches (thr = %.6f) = %.i\n'%(exp_fold_growth_thr,len(mortality_exp_matches))

    #---- Data to plot fold growth in experimental data with cell_pool_concs ----
    # Load exchange reactions corresponding to compounds rescuing each mutant pair
    from auxoMetabs import auxoMetabsMutants
    # Mutants that cannot be rescued at all

    # A dicitonary where keys are exchange reactions corresponding to compounds present in the cell_pool_concs
    # and values are the respective concnetrations
    cell_pool_concs = {'EX_ala_L(e)':6.7e-3,'EX_val_L(e)':9.2e-3,'EX_gly(e)':1.7e-3,'EX_ile_L(e)':1.6e-3,'EX_thr_L(e)':0.2e-3,'EX_leu_L(e)':0.4e-3,'EX_ser_L(e)':0,'EX_pro_L(e)':0.4e-3,'EX_asp_L(e)':0,'EX_cys_L(e)':1.8e-3,'EX_met_L(e)':0.3e-3,'EX_glu_L(e)':25.5e-3,'EX_ph_L(e)':6.8e-3,'EX_tyr_L(e)':0.1e-3,'EX_orn_L(e)':0,'EX_lys_L(e)':0.2e-3,'EX_trp_L(e)':0,'EX_arg_L(e)':0}

    # For each mutant pair check what exchange reactions exist in cell_pool_concs
    # and add them together
    cell_pool_pair = dict([(k,0) for k in expFoldGrowth.keys()])
    
    for (m1,m2) in cell_pool_pair.keys():
        print (m1,m2),'\t',[exch_rxn for m in [m1,m2] for exch_rxn_list in auxoMetabsMutants[re.sub('_Ecoli','',m)] for exch_rxn in exch_rxn_list if exch_rxn in cell_pool_concs.keys()]

        # The cell_pool_pair for each mutant pair is sum of the values for excahange reactions 
        # in cell_pool_concs
        cell_pool_pair[(m1,m2)] = sum([cell_pool_concs[exch_rxn] for m in [m1,m2] for exch_rxn_list in auxoMetabsMutants[re.sub('_Ecoli','',m)] for exch_rxn in exch_rxn_list if exch_rxn in cell_pool_concs.keys()]) 

    # Create a MATLAB matrix where the first column are the values of cell_pool_pair and
    # the second column is the fold growth for each mutant
    with open('results/mortality_results_toPlot_withMatlab.m','a') as outfile:
        outfile.write('\ndeathPool_expFold = [\n')
        for (m1,m2) in expFoldGrowth.keys():
            outfile.write(str(cell_pool_pair[(m1,m2)]) + '\t' + str(expFoldGrowth[(m1,m2)]) + '\n')
        outfile.write('];\n')

def rankAAs():
    """
    This function ranks the amino acids based on their cost
    """
    #--- E. coli iJO1366 model ---
    WT = read_gams_model(gams_model_file = '../models/Ecoli/iJO1366/iJO1366ModelData.py',model_name = 'iJO1366',organism_name = 'E. coli',model_type = 'metabolic')
    WT.biomass_reaction = WT.get_reactions({'Ec_biomass_iJO1366_core_59p81M':'id'})
    WT.all_biomass_reactions = {'core':WT.get_reactions({'Ec_biomass_iJO1366_core_59p81M':'id'}),'WT':WT.get_reactions({'Ec_biomass_iJO1366_WT_59p81M':'id'})}
 
    # Growth medium
    set_specific_bounds(WT,specific_bounds_file = '../models/Ecoli/iJO1366/iJO1366_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign the objective function coefficients
    for rxn in WT.reactions:
        rxn.objective_coefficient = 0
    WT.biomass_reaction.objective_coefficient = 1

    # Perform FBA for the wild-type
    WT.fba(assign_wildType_max_biomass = True)

    from AAList import AAList

    AAcosts = {}

    for AA_exch in AAList:
        # Force to produce one unit of the amino acid
        WT.get_reactions({AA_exch:'id'}).flux_bounds[0] = 1
        WT.fba(build_new_optModel = False, store_opt_fluxes = False,stdout_msgs = False)
        AAcosts[AA_exch] = WT.wildType_max_biomass - WT.fba_model.solution['objective_value'] 
        WT.get_reactions({AA_exch:'id'}).flux_bounds[0] = 0
        
    counter = 0
    for AA in [k[0] for k in sorted(AAcosts.iteritems(),key = lambda x:x[1],reverse = True)]:
        counter += 1
        print '%i\t%s\t%.6f'%(counter,AA,AAcosts[AA]) 

    # Write the results into a file and in the output 
    with open('results/mortality_results_toPlot_withMatlab.m','a') as outfile:
        outfile.write('\nAAcosts = [\n')
        for AA in [k[0] for k in sorted(AAcosts.iteritems(),key = lambda x:x[1],reverse = True)]:
            outfile.write(str(AAcosts[AA]) + '\n')
        outfile.write('];\n')

        outfile.write('\nAAnames = {\n')
        for AA in [k[0] for k in sorted(AAcosts.iteritems(),key = lambda x:x[1],reverse = True)]:
            outfile.write("'" + AA + "'\n")
        outfile.write('};\n')

#---- Create intervals for iterations ------
def create_intervals(total_comb_num,interval_size, X = None):
    """
      total_comb_num: Length of the variable 'combinations'
    interval_size: Desired size of the iteration intervals
                X: The input array. The user can provide the actual input array
                   to test how the code works             

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
       if X == None:
           print slice
       else:
           print slice,'   slice = ',X[slice[0]-1:slice[1]],'\n'
    print '\n'

    return slices


def create_pbs_files(t0,delta_t,tf,total_comb_num,interval_size, outfile_base_name,results_filename_base,with_output = 0):
    """
    Creates the pbs files given:
        total_comb_num: Total number of combinations 
         interval_size: The desired of iteration intervals
     outfile_base_name: The base name of the output files
    results_filename_base: The base name for the files storing the results
                        (see the inputs of enumerate_caes for more details).
           with_output: If 1 the output of the pbs file is dynamically saved at
                        a file called outfile_base_name_start_end.out. If 0 no 
                        output file is stored.      
     The file names will be in the form of outfile_base_name_start_end.pbs

    Short guide on how to read/write files in python:
    http://www.pythonforbeginners.com/files/reading-and-writing-files-in-python
    """
    slices = create_intervals(total_comb_num = total_comb_num,interval_size = interval_size)

    for slice in slices:
        file_name = outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + '.pbs'
        print 'creaitng file ',file_name,'...'
        with open(file_name,'w') as outfile:
            outfile.write('#!/bin/bash\n')
            outfile.write('#PBS -l nodes=1:ppn=1\n')
            outfile.write('#PBS -l walltime=240:00:00\n\n')
            outfile.write('PUTDIR=/data/alizom/EcoliPairs\n')
            outfile.write('cd $PUTDIR\n\n')

            if with_output == 0:
                outfile.write("python -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\"\n\n")
    
                # python -c "from master_func import master_func;master_func(t0 = 0, delta_t = 0.5, tf = 96, start_pos = 1, end_pos = 100,results_filename_base = 'results/DMMM_mortality')" 
                outfile.write('python -c "from master_func import master_func;master_func(t0 = ' + str(t0) + ', delta_t = ' + str(delta_t) + ', tf = ' +str(tf) + ', start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", results_filename_base = '" + results_filename_base + "')\"\n\n")

                outfile.write("python -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n\n")

            elif with_output == 1:
                outfile.write("python -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\" > " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + ".out 2>&1\n\n")

                # python -c "from master_func import master_func;master_func(t0 = 0, delta_t = 0.5, tf = 96, start_pos = 1, end_pos = 100,results_filename_base = 'results/DMMM_mortality')" >> job_DMMM_mortality_1_100.out 2>&1 
                outfile.write('python -c "from master_func import master_func;master_func(t0 = ' + str(t0) + ', delta_t = ' + str(delta_t) + ', tf = ' +str(tf) + ', start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", results_filename_base = '" + results_filename_base + "')\" >> " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + '.out 2>&1\n\n')

                outfile.write("python -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\" >> " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + ".out 2>&1\n\n")
    
            else:
                raise userError('Invalid with_output value.')

            outfile.close()
          
            # make it executable
            os.system('chmod u+x ' + file_name)

#-------------------------
if __name__ == '__main__':
    process_results(metab_thr = 1e-10,total_cell_conc_frac = 0.999998, exp_fold_growth_thr = 50)
