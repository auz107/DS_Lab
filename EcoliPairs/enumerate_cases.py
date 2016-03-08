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
import jsonpickle as jspk

# Last updated: 02-13-2015

#---------------------------------
# This funciton converts the output of json, which is in unicode format to byte stirng
# Source: http://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-ones-from-json-in-python
def byteify(input):
    if isinstance(input, dict):
        return {byteify(key):byteify(value) for key,value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

def performDMMM(input_data):

    t0 = input_data['t0']
    tf = input_data['tf']
    delta_t = input_data['delta_t']
    start_pos = input_data['start_pos']
    end_pos = input_data['end_pos']
    mutant1 = deepcopy(input_data['mutant1'])
    mutant2 = deepcopy(input_data['mutant2'])
    shared_metab_ids  = input_data['shared_metab_ids'][:]
    auxoMetabsMutants_m1  = input_data['auxoMetabsMutants_m1'][:]
    auxoMetabsMutants_m2  = input_data['auxoMetabsMutants_m2'][:]
    results_file_name_base = input_data['results_file_name_base']
    combination_time = byteify(deepcopy(input_data['combination_time']))


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

    DMMM_m1m2 = DMMM_cooperate_defect(community_members = [mutant1,mutant2],shared_metabolites = shared_metabs, time_range = [t0,delta_t,tf],store_dynamic_fluxes = False, screen_output = 'off')
    DMMM_m1m2.run()

    #--- Store the results ---
    # Store only the required information for each shared metabolite
    shared_metabs_stored = {}
    for shared_metab in shared_metabs:
        shared_metabs_stored[shared_metab.id] = shared_metab.concentration[tf]

    # The following outputs are saved: shared_metabs and organism object 
    # for each community member 
    with open(results_file_name_base + '_' + str(start_pos) + '_' + str(end_pos) + '.txt','a') as f:
        f.write(repr({'name':(mutant1.organism.id,mutant2.organism.id),'strategies':{mutant1.organism.id:mutant1.organism.strategy,mutant2.organism.id:mutant2.organism.strategy},'organisms':{mutant1.organism.id:mutant1.organism.cells_per_mL[tf],mutant2.organism.id:mutant2.organism.cells_per_mL[t]},'shared_metabs':shared_metabs_stored}) + '\n')

#--- Creating all possible combinations for strategies ---
def create_strategy_comb(input_data):

    start_pos = input_data['start_pos']
    end_pos = input_data['end_pos']
    time_points = input_data['time_points'][:]
    
    # Determine the number of time points
    # All possible cases at each time point
    possible_strategies = [(('m1','C'),('m2','C')),(('m1','C'),('m2','D')),(('m1','D'),('m2','C')),(('m1','D'),('m2','D'))] 

    # All possible combinations for the total number of time points
    # NOte that here we use len(time_points)-1 because the strategy ihe last time point
    # does not matter according to the time update scheme of DMMM using finite difference
    strategy_combinations = sorted(list(itertools.product(possible_strategies,repeat = len(time_points)-1)))

    # Assign cooperate-cooperate to the last time step and consider only
    # strategy combinations specified by start_pos and end_pos
    strategy_combinations = [c + ((('m1', 'C'), ('m2', 'C')),) for c in strategy_combinations][start_pos - 1:end_pos]

    with open('../tmp/strat_comb_' + str(start_pos) + '_' + str(end_pos) + '.js','w') as f:
        js.dump(strategy_combinations,f)    

#---------- Examine all possible cooperate/defect combinations  ---------------
def  enumerate_cases(t0,delta_t,tf,start_pos,end_pos,results_file_name_base):
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

    #--- Creating all possible combinations for strategies ---
    input_data = {}
    input_data['start_pos'] = start_pos
    input_data['end_pos'] = end_pos
    input_data['time_points'] = time_points
    p = Process(target = create_strategy_comb, args = (input_data,))
    p.start()
    p.join() 
    if p.exitcode > 0:
        raise userError('**ERROR! Error in python subprocess. Please check create_strategy_comb\n')

    # Now load the data from file. The advantage of using a subprocess here is that the occupied space 
    # in memory is only as larage as the slice specified by start_pos and end_pos 
    with open('../tmp/strat_comb_' + str(start_pos) + '_' + str(end_pos) + '.js','r') as f:
        strategy_combinations = js.load(f)    

    #--- DMMM for the co-culture of mutant1 and mutant2 mutants ---
    for (m1,m2) in [('lysA','ilvE')]:
    #for (m1,m2) in mutant_pairs:
      
        print '(',m1,',',m2,')'
    
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
        with open(results_file_name_base + '_' + str(start_pos) + '_' + str(end_pos) + '.txt','w') as f:
            f.write('start_pos = ' + str(start_pos) + '\n')


        # Creating a shared memory using the manager
        input_data = {} 
        input_data['t0'] = t0
        input_data['tf'] = tf
        input_data['delta_t'] = delta_t
        input_data['start_pos'] = start_pos
        input_data['end_pos'] = end_pos
        input_data['mutant1'] = mutant1
        input_data['mutant2'] = mutant2
        input_data['shared_metab_ids'] = shared_metab_ids 
        input_data['auxoMetabsMutants_m1'] = auxoMetabsMutants[m1] 
        input_data['auxoMetabsMutants_m2'] = auxoMetabsMutants[m2] 
        input_data['results_file_name_base'] = results_file_name_base

        print 'Iterations started at',time.strftime('%c'),' ...\n'
        start_time = time.clock()

        for combination_no_time in strategy_combinations:

            counter += 1

            # Add time points to them
            combination_time = dict([(time_points[k],combination_no_time[k]) for k in range(len(time_points))])    
            input_data['combination_time'] = combination_time

            p = Process(target = performDMMM, args = (input_data,))
            p.start()
            p.join() 
            if p.exitcode > 0:
                raise userError('**ERROR! Error in python subprocess. Please check performDMMM\n')

            if counter/500 == int(counter/500):
                print '\ncounter = ',counter,'  current time is: ',time.strftime('%c')
                print 'elapsed time (h) = %.5f  , time per iteration (s) = %.5f' %((time.clock() - start_time)/3600,(time.clock() - start_time)/counter)

        # Remove the temporary file --> Don't remove the file because other codes may still use it
        #os.system('rm -rf ../tmp/strat_comb_' + str(start_pos) + '_' + str(end_pos) + '.js')

#-------- Print results ---------
def print_results(file_name):
    """
    prints the results in the output
    """
    with open(file_name,'r') as f:
        counter = 0
        for line in f:
            counter += 1
            if counter == 1:
                exec(line)
            else:
                exec('res' + str(start_pos + counter - 2) + ' = ' + line)
                exec("print 'res" + str(start_pos + counter - 2) + " = '," + 'res' + str(start_pos + counter - 2) + ",'\\n\\n'")
    
    
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

#-------- Create pbs files ---------
def create_pbs_files(t0,delta_t,tf,interval_size, outfile_base_name,results_file_name_base,with_output = 0):
    """
    Creates the pbs files given:
        total_comb_num: Total number of combinations, and
         interval_size: The desired of iteration intervals
     outfile_base_name: The base name of the output files
    results_file_name_base: The base name for the files storing the results
                        (see the inputs of enumerate_caes for more details).
           with_output: If 1 the output of the pbs file is dynamically saved at
                        a file called outfile_base_name_start_end.out. If 0 no 
                        output file is stored.      
     The file names will be in the form of outfile_base_name_start_end.pbs

    Short guide on how to read/write files in python:
    http://www.pythonforbeginners.com/files/reading-and-writing-files-in-python
    """
    # Total number of combinations. Subtract delta_t from tf because the last time point
    # does not have any effect on simulations
    total_comb_num = np.power(4,int(((tf-delta_t) - t0)/delta_t) + 1)

    print 'total_comb_num = ',total_comb_num

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
    
                # python -c 'from enumerate_cases import enumerate_cases;enumerate_cases(t0=0,delta_t=1,tf=10,start_pos=1,end_pos=2)'
                outfile.write('python -c "from enumerate_cases import enumerate_cases;enumerate_cases(t0 = ' + str(t0) + ', delta_t = ' + str(delta_t) + ', tf = ' +str(tf) + ', start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", results_file_name_base = '" + results_file_name_base + "')\"\n\n")

                outfile.write("python -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n\n")

            elif with_output == 1:
                outfile.write("python -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\" > " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + ".out 2>&1\n\n")

                # python -c 'from enumerate_cases import enumerate_cases;enumerate_cases(t0=0,delta_t=1,tf=10,start_pos=1,end_pos=2)' > job_emc_startPos_endPos.out 2>&1
                outfile.write('python -c "from enumerate_cases import enumerate_cases;enumerate_cases(t0 = ' + str(t0) + ', delta_t = ' + str(delta_t) + ', tf = ' +str(tf) + ', start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", results_file_name_base = '" + results_file_name_base + "')\" >> " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + '.out 2>&1\n\n')

                outfile.write("python -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\" >> " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + ".out 2>&1\n\n")
    
            else:
                raise userError('Invalid with_output value.')

            outfile.close()
          
            # make it executable
            os.system('chmod u+x ' + file_name)
