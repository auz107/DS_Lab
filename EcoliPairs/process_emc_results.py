from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../')
from copy import deepcopy
import cPickle as pk
import shelve
import itertools
from tools.userError import *
from tools.io.read_gams_model import read_gams_model
from tools.core.metabolite import metabolite
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from DMMM_cooperate_defect import DMMM_cooperate_defect
from tools.fba.set_specific_bounds import set_specific_bounds
from metabolite_stored import metabolite_stored
from osize import asizeof
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)
from multiprocessing import Process, Manager
import json as js
import matplotlib.pyplot as plt
from scipy.cluster.vq import *
import numpy as np
import re
from copy import deepcopy

# Last updated: 01-26-2015


def process_results(results_file_names,coopr_thr):
    """
    results_file_names: An array containing all results file names
             coopr_thr: Threshold above which the pairs can be considered
                        as cooperative
    """

    #------- Find the max and min fold change in results -----    
    # Initial cell concentrations
    cell_conc_init = 7.5e6
    
    # The current maximum fold change
    max_foldChange_lysA = 1
    max_foldChange_ilvE = 1
    min_foldChange_lysA = 1
    min_foldChange_ilvE = 1

    number_of_results = 0

    for result_file in results_file_names:
        with open(result_file,'r') as f:
            line_number = 0
            for line in f:
                line_number += 1
                if line_number == 1:
                    exec(line)
                else:
                    number_of_results += 1
                    exec('result = ' + line) 

                    if number_of_results == 1:
                        # Time vector
                        t = result['strategies']['lysA_Ecoli'].keys()

                    # Identify the cases where both cooperate at all time points
                    if 'D' not in result['strategies']['lysA_Ecoli'].values() + result['strategies']['ilvE_Ecoli'].values():
                        ALLC_result = result

                    # Fold change for each organism
                    lysA_foldChange = result['organisms']['lysA_Ecoli']/cell_conc_init
                    if lysA_foldChange > max_foldChange_lysA:
                        max_foldChange_lysA = lysA_foldChange
                        max_lysA_result = result
                    min_foldChange_lysA = min(min_foldChange_lysA,lysA_foldChange)
                    ilvE_foldChange = result['organisms']['ilvE_Ecoli']/cell_conc_init
                    if ilvE_foldChange > max_foldChange_ilvE:
                        max_foldChange_ilvE = ilvE_foldChange
                        max_ilvE_result = result
                    min_foldChange_ilvE = min(min_foldChange_ilvE,ilvE_foldChange)
    
    # Find the maximum fold change in cell concetrationsn when both 
    # cooperate in all time points
    ALLC_lysA_foldChange = ALLC_result['organisms']['lysA_Ecoli']/cell_conc_init
    ALLC_ilvE_foldChange = ALLC_result['organisms']['ilvE_Ecoli']/cell_conc_init
    
    print '\nTotal # of results = %i\n' % number_of_results
    print '\nmax_foldChange_lysA = %.4f   ,    max_foldChange_ilvE = %.4f\n' % (max_foldChange_lysA,max_foldChange_ilvE)
    print '\nALLC_lysA_foldChange = %.4f   ,    ALLC_ilvE_foldChange = %.4f\n' % (ALLC_lysA_foldChange,ALLC_ilvE_foldChange)
    print '\nmax_lysA_result = ',max_lysA_result,'\n'
    print '\nmax_ilvE_result = ',max_ilvE_result,'\n'
    print '\nTotal # of results = %i ,  min_foldChange_lysA = %.4f   ,    min_foldChange_ilvE = %.4f\n' % (number_of_results,min_foldChange_lysA,min_foldChange_ilvE)

    #------- Find which strategies lead to cooperation -----    
    # Consider 90% of max_fold_change as cooperation
    cooper_fold_change_lysA = (ALLC_lysA_foldChange - 1)*coopr_thr + 1
    cooper_fold_change_ilvE = (ALLC_ilvE_foldChange - 1)*coopr_thr + 1

    print '\ncooper_fold_change_lysA = ',cooper_fold_change_lysA,'\n'
    print 'cooper_fold_change_ilvE = ',cooper_fold_change_ilvE,'\n\n'

    # Pairs of cooperative strategies defined as those whose fold increase is 90% of the maximum
    # fold increase for both species 
    successful_strategies = []
    y_vec = []
    coopr_counter = 0

    for result_file in results_file_names:
        with open(result_file,'r') as f:
            line_number = 0
            results = {}
            for line in f:
                line_number += 1
                if line_number == 1:
                    exec(line)
                else:
                    exec('result = ' + line) 
                    lysA_foldChange = result['organisms']['lysA_Ecoli']/cell_conc_init
                    ilvE_foldChange = result['organisms']['ilvE_Ecoli']/cell_conc_init
                    if lysA_foldChange >= cooper_fold_change_lysA and ilvE_foldChange >= cooper_fold_change_ilvE:
                        coopr_counter += 1

                        successful_strategies.append(result['strategies'])

                        print '(%i) ' % coopr_counter,'\n'
                        print 'lysA Ecoli:',result['strategies']['lysA_Ecoli'],'   fc = ',cooper_fold_change_lysA,'\n'
                        print 'ilvE Ecoli:',result['strategies']['ilvE_Ecoli'],'   fc = ',cooper_fold_change_ilvE,'\n\n'

                        # Vector to plot. At each time point assign 1 to CC, 2 to CD, 3 to DC and 4 to DD
                        # with first letter representing the strategy of lysA and the second that of ilvE
                        y = []
                        for tp in t:
                            if result['strategies']['lysA_Ecoli'][tp] == 'C' and result['strategies']['ilvE_Ecoli'][tp] == 'C':
                                y.append(1) 
                            elif result['strategies']['lysA_Ecoli'][tp] == 'C' and result['strategies']['ilvE_Ecoli'][tp] == 'D':
                                y.append(2) 
                            elif result['strategies']['lysA_Ecoli'][tp] == 'D' and result['strategies']['ilvE_Ecoli'][tp] == 'C':
                                y.append(3) 
                            elif result['strategies']['lysA_Ecoli'][tp] == 'D' and result['strategies']['ilvE_Ecoli'][tp] == 'D':
                                y.append(4) 
                            else:
                                raise userError('\n**ERROR! Unknown case. Not any of four cases CC, CD, DC or DD\n')

                        y_vec.append(y)

    with open('cooperative_strategies.txt','w') as f:
        for strategy in y_vec:
            for v in strategy:
                f.write(str(v) + '\t') 
            f.write('\n')

    print '\nTotal # of cases with cooperation  = %i\n' % len(y_vec)

    # Find sucessfull strategies that are unique. To this end, consider the 
    # strategies which are different only at their strategy in the final time
    # point as identifcal. This is becuase the strategy in the final time point
    # does not affect the concentrations
    successful_strategies_unique = deepcopy(successful_strategies)

    cs = 0

    for succ_strategy in successful_strategies_unique:
         # Final time point
         tf = succ_strategy['lysA_Ecoli'].keys()[-1]    
         for other_succ_strategy in [strat for strat in successful_strategies_unique if strat != succ_strategy]:
             s1 = deepcopy(succ_strategy)
             s2 = deepcopy(other_succ_strategy)
  
             del s1['lysA_Ecoli'][tf]             
             del s1['ilvE_Ecoli'][tf]             
             del s2['lysA_Ecoli'][tf]             
             del s2['ilvE_Ecoli'][tf]             
             if s1 == s2:
                 cs += 1
                 # index of other_succ_strategy
                 idx = successful_strategies_unique.index(other_succ_strategy)
                 del successful_strategies_unique[idx]

    # Add result_ALLC to the all outputs for the reference
    if ALLC_result['strategies'] not in successful_strategies:
        successful_strategies.append(ALLC_result['strategies'])
    if ALLC_result['strategies'] not in successful_strategies_unique:
        successful_strategies_unique.append(ALLC_result['strategies'])
        y_vec.append([1 for i in range(len(t))])

    return (y_vec,successful_strategies,successful_strategies_unique)

#----------------------------------------
def plot_results(v_vec):

    # Set the ytick labels
    plt.yticks([0,1,2,3,4])

    # The following sets the font of all items including xlable, ylabel, title, xtick and ytick labels
    plt.rc('font',**{'family':'normal','size':20,'weight':'bold'})

    plt.xlabel('Time',{'weight':'bold','size':36})
    plt.ylabel('1 = CC, 2 = CD, 3 == DC, 4 = DD',{'weight':'bold','size':30})

    plt.show()
    plt.savefig('graph.pdf')


#------------------------------------------------------
def performDMMM(input_data):

    t0 = input_data['t0']
    tf = input_data['tf']
    delta_t = input_data['delta_t']
    mutant1 = deepcopy(input_data['mutant1'])
    mutant2 = deepcopy(input_data['mutant2'])
    shared_metab_ids  = input_data['shared_metab_ids'][:]
    auxoMetabsMutants_m1  = input_data['auxoMetabsMutants_m1'][:]
    auxoMetabsMutants_m2  = input_data['auxoMetabsMutants_m2'][:]
    results_file_name_base = input_data['results_file_name_base']
    combination_time = deepcopy(input_data['combination_time'])

    # Consider each set of metabolites rescuing a mutant separately
    for exch_rxns_list_m1 in auxoMetabsMutants_m1:    
        for exch_rxns_list_m2 in auxoMetabsMutants_m2:   
            shared_metabs = []
            # Exchange rxns for all shared metabolites each mutant takes up 
            uptake_rxns_m1 = []
            uptake_rxns_m2 = []
            # Exchange rxns for metabolites each mutant should produce to help its partner 
            export_rxns_m1 = []
            export_rxns_m2 = []

            # First include glucose 
            uptake_rxns_m1.append(mutant1.get_reactions({'EX_glc(e)':'id'})) 
            uptake_rxns_m2.append(mutant2.get_reactions({'EX_glc(e)':'id'}))
            # Then the rest of shared metabolties
            for exch_rxn_id in exch_rxns_list_m1:
                uptake_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
                export_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 
                mutant2.cooperative_rxns.append(mutant2.get_reactions({exch_rxn_id:'id'})) 
            for exch_rxn_id in exch_rxns_list_m2:
                export_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
                uptake_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 
                mutant1.cooperative_rxns.append(mutant1.get_reactions({exch_rxn_id:'id'})) 

            # Glucose as a shared metabolite (concentration in mM)
            glucose = metabolite(id = 'glc-D[e]', name = 'D-Glucose', Kegg_id = 'C00031', reactant_reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],concentration = {0:111.01})    

            shared_metabs.append(glucose)

            for shared_metab_id in [id for id in shared_metab_ids if id != 'glc-D[e]']:
                reactant_rxns = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.reactants[0].id.lower() == shared_metab_id.lower()]
                product_rxns = [r for r in export_rxns_m1 + export_rxns_m2 if r.reactants[0].id.lower() == shared_metab_id.lower()]
                shared_metab = metabolite (id = shared_metab_id,reactant_reactions = reactant_rxns, product_reactions = product_rxns, concentration = {0:0}) 
                shared_metabs.append(shared_metab)
    

    # Perfomr some initializations
    #for shared_metab in shared_metabs:
    #    shared_metab.concentration = {0:shared_metab.concentration[0]}
    #mutant1.organism.cells_per_mL = {0:mutant1.organism.cells_per_mL[0]} 
    #mutant2.organism.cells_per_mL = {0:mutant2.organism.cells_per_mL[0]} 

    # Strategies for mutants 1 and 2
    mutant1.organism.strategy = dict([(t,dict(combination_time[t])['m1']) for t in combination_time.keys()])
    mutant2.organism.strategy = dict([(t,dict(combination_time[t])['m2']) for t in combination_time.keys()])

    DMMM_m1m2 = DMMM_cooperate_defect(community_members = [mutant1,mutant2],shared_metabolites = shared_metabs, time_range = [t0,delta_t,tf],screen_output = 'off')
    DMMM_m1m2.run()

    #--- Store the results ---
    # Store only the required information for each shared metabolite
    shared_metabs_stored = {}
    for shared_metab in shared_metabs:
        shared_metabs_stored[shared_metab.id] = shared_metab.concentration

    # The following outputs are saved: shared_metabs and organism object 
    # for each community member 
    with open(results_file_name_base + '.py','a') as f:
        f.write('full_results.append(' + repr({'name':(mutant1.organism.id,mutant2.organism.id),'strategies':{mutant1.organism.id:mutant1.organism.strategy,mutant2.organism.id:mutant2.organism.strategy},'organisms':{mutant1.organism.id:mutant1.organism.cells_per_mL,mutant2.organism.id:mutant2.organism.cells_per_mL},'shared_metabs':shared_metabs_stored}) + ')\n')


    time_points = mutant1.organism.strategy.keys()
    with open(results_file_name_base + '.m','a') as f:
        f.write('\ncurrent_index = length(results) + 1;\n')
        f.write("results(current_index).organism_names = {'" + mutant1.organism.id + "','" + mutant2.organism.id + "'};\n")
        f.write('results(current_index).time = ' + str(time_points) + ';\n')

        # mutant1_strategy and mutant2_strategy take a value of 0 (D) or 1 (C) for mutant 1 
        # at each time point
        mutant1_strategy = []
        for t in time_points:
            if mutant1.organism.strategy[t] == 'C':
                mutant1_strategy.append(1)
            elif mutant1.organism.strategy[t] == 'D':
                mutant1_strategy.append(0)
            else: 
                raise userError('**ERROR! Undefined strategy pair (not in CC, CD, DC, DD)\n')
        f.write('results(current_index).strategies.' + mutant1.organism.id + ' = ' + str(mutant1_strategy) + ';\n')

        mutant2_strategy = []
        for t in time_points:
            if mutant2.organism.strategy[t] == 'C':
                mutant2_strategy.append(1)
            elif mutant2.organism.strategy[t] == 'D':
                mutant2_strategy.append(0)
            else: 
                raise userError('**ERROR! Undefined strategy pair (not in CC, CD, DC, DD)\n')
        f.write('results(current_index).strategies.' + mutant2.organism.id + ' = ' + str(mutant2_strategy) + ';\n')

        # Cell concentrations
        mutant1_cell_concs = []
        mutant2_cell_concs = []
        for t in time_points:
            mutant1_cell_concs.append(mutant1.organism.cells_per_mL[t]) 
            mutant2_cell_concs.append(mutant2.organism.cells_per_mL[t]) 
        f.write('results(current_index).cell_concs.' + mutant1.organism.id + ' = ' + str(mutant1_cell_concs) + ';\n')
        f.write('results(current_index).cell_concs.' + mutant2.organism.id + ' = ' + str(mutant2_cell_concs) + ';\n')

        # Shared metab concnetrations
        for shared_metab_name in shared_metabs_stored.keys():
            shared_metab_name_forMATLAB = re.sub('\[','_',shared_metab_name)
            shared_metab_name_forMATLAB = re.sub('\]','',shared_metab_name_forMATLAB)
            shared_metab_name_forMATLAB = re.sub('-','_',shared_metab_name_forMATLAB)
            shared_metab_concs = []
            for t in time_points:
                shared_metab_concs.append(shared_metabs_stored[shared_metab_name][t])
            f.write('results(current_index).shared_metab_concs.' + shared_metab_name_forMATLAB + ' = ' + str(shared_metab_concs) + ';\n')
            

#-------------------------------------------------------------------
def  analyze_succ_strat(t0,delta_t,tf,successful_strategies,results_file_name_base):
    """
     Analyzes successful strategies further by redoing DMMM and storing all results

     INPUTS:
     -------
             t0: Initial simulaiton time
        delta_t: Time step
             tf: Total simulation time
      strat_pos: Start position of the array containing all possible
                 strategy combinations (see variable 'combinations').
        end_pos: End position of the array containing all possible
                 strategy combinations (see variable 'combinations').
     results_file_name_base: The base for the file name storing the results
                 Example: results_file_name_base = 'results/emc_results'. The code
                 will add the start and end positions to the file name.
                 Example: 'results/emc_results_1_500.txt'

    NOTE: The user can enter both the start and end  positions numbers 
          assuming that they start at one. The code will take care of
          array indexing convention of python (indexing starts at zero) 

    Shelve use: http://stackoverflow.com/questions/13851438/how-to-store-a-big-dictionary

    Example of how to open the shelve file:
    from contextlib import closing
    import shelve
    with closing(shelve.open('results/enumeration_results_1_2.shelf')) as d:
        for s in d.keys():
            print d[s],'\n'
    """

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]    
    print 'time_points = ',time_points

    #--- E. coli iAF1260 model ---
    print '\n--- Wild-type E.coli (iAF1260 model) ----'
    WT = read_gams_model(gams_model_file = '../models/Ecoli/iAF1260/iAF1260ModelData.py',model_name = 'iAF1260',organism_name = 'E. coli',model_type = 'metabolic')
    WT.biomass_reaction = WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'})
    WT.all_biomass_reactions = {'core':WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'}),'WT':WT.get_reactions({'Ec_biomass_iAF1260_WT_59p81M':'id'})}
    WT.organism.gDW_per_cell = 2.8e-13
 
    # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
    # Example: EX_glc(E): glc-D[e] <==>    Vmax*C['glc-D[e]']/(Km + C['glc-D[e]']) 
    # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
    for reaction in [r for r in WT.reactions if r.type.lower() == 'exchange']:
        # The id of metabolite participating in the exchange reaction
        metab_id = [m.id for m in reaction.metabolites][0]
        reaction.kinetics = "10*C['" + metab_id + "']/(10 + C['" + metab_id + "'])"

    # Glucose uptake kinetics 
    exch_rxns = WT.get_reactions({'EX_glc(e)':'id','EX_lys-L(e)':'id','EX_ile-L(e)':'id'})
    exch_rxns['EX_glc(e)'].kinetics = "10*C['glc-D[e]']/(0.15 + C['glc-D[e]'])"
    exch_rxns['EX_lys-L(e)'].kinetics = "0.1964*C['lys-L[e]']/(5e-4 + C['lys-L[e]']) + 0.3055*C['lys-L[e]']/(1e-2 + C['lys-L[e]'])"
    exch_rxns['EX_ile-L(e)'].kinetics = "0.0346*C['ile-L[e]']/(1.22e-3 + C['ile-L[e]'])"
   
    # Growth medium
    set_specific_bounds(WT,specific_bounds_file = '../models/Ecoli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign the objective function coefficients
    for rxn in WT.reactions:
        rxn.objective_coefficient = 0
    WT.biomass_reaction.objective_coefficient = 1

    # Perform FBA for the wild-type
    WT.fba(assign_wildType_max_biomass = True)

    # Compute ms and biomass yield for glucose
    glc_D = WT.get_metabolites({'glc-D[e]':'id'})
    glc_D.ms_calc() 
    glc_D.biomass_yield_calc() 

    #--- Load the list of rxns that must be off in each mutant ---
    from mutants_rxn_info import mutants_rxn_info

    # All possible pair combinations
    mutant_pairs = list(itertools.combinations(mutants_rxn_info.keys(),r=2))   

    #--- Load the list of exchange rxns for metabolites each mutant needs to survive ---
    with open('results/auxoMetabs.pk','rb') as inputFile:
        termCond_mutant,auxoMetabsMutants = pk.load(inputFile)

    # IDs of the shared metabolites
    shared_metab_ids = ['glc-D[e]']

    #--- DMMM for the co-culture of mutant1 and mutant2 mutants ---
    (m1,m2) = ('lysA','ilvE')

    #-- Mutant 1 mutant and the related reactions whose flux should be set to zero --
    mutant1 = deepcopy(WT)
    mutant1.id = WT.id + '_' + m1
    mutant1.organism.id = m1 + '_Ecoli'
    mutant1.organism.cells_per_mL = {0:7.5e6}
    mutant1.organism.gDW_per_mL = {0:7.5e6*mutant1.organism.gDW_per_cell}
    for rxn_id in mutants_rxn_info[m1]:
        rxn = mutant1.get_reactions({rxn_id:'id'})
        rxn.flux_bounds = [0,0]
        
    mutant1.fba(create_model = False)

    # Compute ms and biomass yield for metabolites needed to rescue the mutant 
    for exch_rxn_id in list(set([r for rList in auxoMetabsMutants[m1] for r in rList])): 
        # The metabolite participating in the exchange reaction
        exch_rxn = mutant1.get_reactions({exch_rxn_id:'id'})
        metab = exch_rxn.reactants[0] 
        shared_metab_ids.append(metab.id)
        metab.ms_calc()    
        metab.biomass_yield_calc()    

    # Exchange reactions for cooperation (i.e., to export metabolites rescuing the partner)
    mutant1.cooperative_rxns = []

    #-- Mutant 2 and related reaction whose flux should be set to zero --
    print '\n--- ' + m2 + '_Ecoli mutant ----'
    mutant2 = deepcopy(WT)
    mutant2.id = WT.id + '_' + m2
    mutant2.organism.id = m2 + '_Ecoli'
    mutant2.organism.cells_per_mL = {0:7.5e6}
    mutant2.organism.gDW_per_mL = {0:7.5e6*mutant2.organism.gDW_per_cell}
    for rxn_id in mutants_rxn_info[m2]:
        rxn = mutant2.get_reactions({rxn_id:'id'})
        rxn.flux_bounds = [0,0]
        
    mutant2.fba(create_model = False)

    # Compute ms and biomass yield for metabolites needed to rescue the mutant 
    for exch_rxn_id in list(set([r for rList in auxoMetabsMutants[m2] for r in rList])): 
        exch_rxn = mutant2.get_reactions({exch_rxn_id:'id'})
        # The metabolite participating in the exchange reaction
        metab = exch_rxn.reactants[0] 
        shared_metab_ids.append(metab.id)
        metab.ms_calc()    
        metab.biomass_yield_calc()    

    # Exchange reactions for cooperation (i.e., to export metabolites rescuing the partner)
    mutant2.cooperative_rxns = []

    # --- Define the metabolites available in the extracellular medium ---
    # Consider only unique elements as, in general, it is quite possible that
    # two mutants need the same compounds to survive
    shared_metab_ids = sorted(list(set(shared_metab_ids)))

    counter = 0

    # Initialize the file
    with open(results_file_name_base + '.py','w') as f:
        f.write('full_results = []\n')
    with open(results_file_name_base + '.m','w') as f:
        f.write('results = struct([]);\n')

    # Creating a shared memory using the manager
    input_data = {} 
    input_data['t0'] = t0
    input_data['tf'] = tf
    input_data['delta_t'] = delta_t
    input_data['mutant1'] = mutant1
    input_data['mutant2'] = mutant2
    input_data['shared_metab_ids'] = shared_metab_ids 
    input_data['auxoMetabsMutants_m1'] = auxoMetabsMutants[m1] 
    input_data['auxoMetabsMutants_m2'] = auxoMetabsMutants[m2] 
    input_data['results_file_name_base'] = results_file_name_base

    print 'Iterations started at',time.strftime('%c'),' ...\n'
    start_time = time.clock()

    # Convert successful_strategies to the format needed by performDMMM 
    combination_times = []
    for strategy_pair in successful_strategies:
        combination_time = {}
        for t in time_points:
            combination_time[t] = (('m1',strategy_pair['lysA_Ecoli'][t]),('m2',strategy_pair['ilvE_Ecoli'][t]))
        combination_times.append(combination_time)

    for combination_time in combination_times:

        counter += 1

        # Add time points to them
        input_data['combination_time'] = combination_time

        p = Process(target = performDMMM, args = (input_data,))
        p.start()
        p.join() 
        if p.exitcode > 0:
            raise userError('**ERROR! Error in python subprocess. Please check performDMMM\n')

        if counter/500 == int(counter/500):
            print '\ncounter = ',counter,'  current time is: ',time.strftime('%c')
            print 'elapsed time (h) = %.5f  , time per iteration (s) = %.5f' %((time.clock() - start_time)/3600,(time.clock() - start_time)/counter)


#-------------------------
if __name__ == '__main__':
    from emc_results_file_names import emc_results_file_names
    (y_vec,successful_strategies,successful_strategies_unique) = process_results(results_file_names = emc_results_file_names, coopr_thr = 0.90)
    print '# of successful_strategies = %i\n'%(len(successful_strategies)-1)
    print '# of successful_strategies_unique = %i\n'%(len(successful_strategies_unique)-1)
    for succ_strat in successful_strategies_unique:
        print succ_strat  

    if len(successful_strategies) > 1:
        # plot_results(y_vect)
        t0, delta_t, tf = 0, 1, 10
        #analyze_succ_strat(t0 = t0,delta_t = delta_t, tf = tf, successful_strategies = successful_strategies_unique,results_file_name_base = 'results/results_emc1_succ_strats_' + str(t0) + '_' + str(tf))
        analyze_succ_strat(t0 = t0,delta_t = delta_t, tf = tf, successful_strategies = successful_strategies_unique,results_file_name_base = 'results/succ_strats_' + str(t0) + '_' + str(tf))



