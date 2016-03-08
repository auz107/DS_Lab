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
from tools.io.read_gams_model import read_gams_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from DMMM_moma import DMMM_moma
from tools.fba.set_specific_bounds import set_specific_bounds
from compound_stored import compound_stored
from tools.importData import importData
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)
from multiprocessing import Process, Manager
import re

# Last updated: 02-12-2015

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
    time_points = input_data['time_points'][:]
    mutant1 = deepcopy(input_data['mutant1'])
    mutant2 = deepcopy(input_data['mutant2'])
    shared_cmp_ids  = input_data['shared_cmp_ids'][:]
    auxoMetabsMutants_m1  = input_data['auxoMetabsMutants_m1'][:]
    auxoMetabsMutants_m2  = input_data['auxoMetabsMutants_m2'][:]
    results_file_name_base = input_data['results_file_name_base']

    print '\nauxoMetabsMutants_m1 = ',auxoMetabsMutants_m1
    print '\nauxoMetabsMutants_m2 = ',auxoMetabsMutants_m2,'\n'

    # Consider each set of compounds rescuing a mutant separately
    for exch_rxns_list_m1 in auxoMetabsMutants_m1:    
        for exch_rxns_list_m2 in auxoMetabsMutants_m2:   
            shared_cmps = []
            # Exchange rxns for all shared compounds each mutant takes up 
            uptake_rxns_m1 = []
            uptake_rxns_m2 = []
            # Exchange rxns for compounds each mutant should produce to help its partner 
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

            # Glucose as a shared compound (concentration in mM)
            glucose = compound(id = 'glc-D[e]', name = 'D-Glucose', Kegg_id = 'C00031', reactant_reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],concentration = {0:111.01})    

            shared_cmps.append(glucose)

            for shared_cmp_id in [id for id in shared_cmp_ids if id != 'glc-D[e]']:
                reactant_rxns = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
                product_rxns = [r for r in export_rxns_m1 + export_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
                shared_cmp = compound (id = shared_cmp_id,reactant_reactions = reactant_rxns, product_reactions = product_rxns, concentration = {0:0}) 
                shared_cmps.append(shared_cmp)
    

    print '\nshared metabs ids = ',shared_cmp_ids,'\n'

    DMMM_m1m2 = DMMM_moma(community_members = [mutant1,mutant2],shared_compounds = shared_cmps, time_range = [t0,delta_t,tf],store_dynamic_fluxes = False, screen_output = 'off')
    DMMM_m1m2.run()

    #--- Store the results ---
    # Store only the required information for each shared compound
    shared_cmps_stored = {}
    for shared_cmp in shared_cmps:
        shared_cmps_stored[shared_cmp.id] = shared_cmp.concentration

    # The following outputs are saved: shared_cmps and organism object 
    # for each community member 
    with open(results_file_name_base + '.py','a') as f:
        f.write('results.append(' + repr({'name':(mutant1.organism.id,mutant2.organism.id),'cell_concs':{mutant1.organism.id:mutant1.organism.cells_per_mL,mutant2.organism.id:mutant2.organism.cells_per_mL},'shared_cmp_concs':shared_cmps_stored}) + ')\n')

    # Write the results into a MATLAB file to plot later
    with open(results_file_name_base + '.m','a') as f:
        f.write('\ncurrent_index = length(results) + 1;\n')
        f.write("results(current_index).organism_names = {'" + mutant1.organism.id + "','" + mutant2.organism.id + "'};\n")
        f.write('results(current_index).time = ' + str(time_points) + ';\n')

        # Cell concentrations
        mutant1_cell_concs = []
        mutant2_cell_concs = []
        for t in time_points:
            mutant1_cell_concs.append(mutant1.organism.cells_per_mL[t])
            mutant2_cell_concs.append(mutant2.organism.cells_per_mL[t])
        f.write('results(current_index).cell_concs.' + mutant1.organism.id + ' = ' + str(mutant1_cell_concs) + ';\n')
        f.write('results(current_index).cell_concs.' + mutant2.organism.id + ' = ' + str(mutant2_cell_concs) + ';\n')

        # Shared cmp concnetrations
        for shared_cmp_name in shared_cmps_stored.keys():
            shared_cmp_name_forMATLAB = re.sub('\[','_',shared_cmp_name)
            shared_cmp_name_forMATLAB = re.sub('\]','',shared_cmp_name_forMATLAB)
            shared_cmp_name_forMATLAB = re.sub('-','_',shared_cmp_name_forMATLAB)
            shared_cmp_concs = []
            for t in time_points:
                shared_cmp_concs.append(shared_cmps_stored[shared_cmp_name][t])
            f.write('results(current_index).shared_cmp_concs.' + shared_cmp_name_forMATLAB + ' = ' + str(shared_cmp_concs) + ';\n')


#---------- DMMM with moma for all mutant combinations ---------------
def  mutants_moma(t0,delta_t,tf,start_pos,end_pos,results_file_name_base):
    """
     INPUTS:
     -------
             t0: Initial simulaiton time
        delta_t: Time step
             tf: Total simulation time
     results_file_name_base: The base for the file name storing the results
                 Example: results_file_name_base = 'results/emc_results'. The code
                 will add the start and end positions to the file name.
                 Example: 'results/emc_results_1_500.txt'
    """

    results_file_name_base = results_file_name_base + '_' + str(start_pos) + '_' + str(end_pos)

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]    
    print '\ntime_points = ',time_points

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
        # The id of compound participating in the exchange reaction
        cmp_id = [m.id for m in reaction.compounds][0]
        reaction.kinetics = "10*C['" + cmp_id + "']/(10 + C['" + cmp_id + "'])"

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

    # Store wild-type fluxes
    for rxn in WT.reactions:
        rxn.wildtype_flux = rxn.flux 

    # Compute ms and biomass yield for glucose
    glc_D = WT.get_compounds({'glc-D[e]':'id'})
    glc_D.ms_calc() 
    glc_D.biomass_yield_calc() 

    #--- Load the list of rxns that must be off in each mutant ---
    from mutants_rxn_info import mutants_rxn_info

    #--- Load the list of exchange rxns for compounds each mutant needs to survive ---
    with open('results/auxoMetabs.pk','rb') as inputFile:
        termCond_mutant,auxoMetabsMutants = pk.load(inputFile)
    # Mutants that cannot be rescued at all
    not_rescued_mutants = [m for m in auxoMetabsMutants.keys() if len(auxoMetabsMutants[m]) == 0]

    # All possible pair combinations (do not consider the ones that cannot be rescued)
    mutant_pairs = list(itertools.combinations([m for m in mutants_rxn_info.keys() if m not in not_rescued_mutants],r=2))   
    print '\nThe total # of mutant pairs to examine = %i' % len(mutant_pairs)

    # IDs of the shared compounds
    shared_cmp_ids = ['glc-D[e]']

    # Initialize the file
    with open(results_file_name_base + '.py','w') as f:
        f.write('results = []\n')
    with open(results_file_name_base + '.m','w') as f:
        f.write('results = struct([]);\n')

    #--- DMMM for the co-culture of mutant1 and mutant2 mutants ---
    #for (m1,m2) in [('glnA','trpC'),('lysA','ilvE')]:
    #for (m1,m2) in [('ppc','hisD')]:
    for (m1,m2) in mutant_pairs[start_pos - 1:end_pos]:
      
        print '\n**** %i. (%s,%s) ****\n' % (mutant_pairs.index((m1,m2)) + 1,m1,m2)
    
        #-- Mutant 1 and the related reactions whose flux should be set to zero --
        print '\n-- ' + m1 + '_Ecoli mutant --'
        mutant1 = deepcopy(WT)
        mutant1.id = WT.id + '_' + m1
        mutant1.organism.id = m1 + '_Ecoli'
        mutant1.organism.cells_per_mL = {0:7.5e6}
        mutant1.organism.gDW_per_mL = {0:7.5e6*mutant1.organism.gDW_per_cell}
        for rxn_id in mutants_rxn_info[m1]:
            rxn = mutant1.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]

        mutant1.fba(create_model = False, store_opt_fluxes = False)

        # Compute ms and biomass yield for compounds needed to rescue the mutant 
        for exch_rxn_id in list(set([r for rList in auxoMetabsMutants[m1] for r in rList])): 
            # The compound participating in the exchange reaction
            exch_rxn = mutant1.get_reactions({exch_rxn_id:'id'})
            cmp = exch_rxn.reactants[0] 
            shared_cmp_ids.append(metab.id)
            metab.ms_calc()    
            metab.biomass_yield_calc()    
    
        # Set all objective coefficients to 1 for the MOMA
        for rxn in mutant1.reactions:
            rxn.objective_coefficient = 1

        # Exchange reactions for cooperation (i.e., to export compounds rescuing the partner)
        mutant1.cooperative_rxns = []

        #-- Mutant 2 and related reaction whose flux should be set to zero --
        print '\n-- ' + m2 + '_Ecoli mutant --'
        mutant2 = deepcopy(WT)
        mutant2.id = WT.id + '_' + m2
        mutant2.organism.id = m2 + '_Ecoli'
        mutant2.organism.cells_per_mL = {0:7.5e6}
        mutant2.organism.gDW_per_mL = {0:7.5e6*mutant2.organism.gDW_per_cell}
        for rxn_id in mutants_rxn_info[m2]:
            rxn = mutant2.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]
            
        mutant2.fba(create_model = False, store_opt_fluxes = False)
    
        # Compute ms and biomass yield for compounds needed to rescue the mutant 
        for exch_rxn_id in list(set([r for rList in auxoMetabsMutants[m2] for r in rList])): 
            exch_rxn = mutant2.get_reactions({exch_rxn_id:'id'})
            # The compound participating in the exchange reaction
            cmp = exch_rxn.reactants[0] 
            shared_cmp_ids.append(metab.id)
            metab.ms_calc()    
            metab.biomass_yield_calc()    
    
        # Set all objective coefficients to 1 for the MOMA
        for rxn in mutant2.reactions:
            rxn.objective_coefficient = 1

        # Exchange reactions for cooperation (i.e., to export compounds rescuing the partner)
        mutant2.cooperative_rxns = []
    
        # --- Define the compounds available in the extracellular medium ---
        # Consider only unique elements as, in general, it is quite possible that
        # two mutants need the same compounds to survive
        shared_cmp_ids = sorted(list(set(shared_cmp_ids)))
    
        # Creating a shared memory using the manager
        input_data = {} 
        input_data['t0'] = t0
        input_data['tf'] = tf
        input_data['delta_t'] = delta_t
        input_data['time_points'] = time_points
        input_data['mutant1'] = mutant1
        input_data['mutant2'] = mutant2
        input_data['shared_cmp_ids'] = shared_cmp_ids 
        input_data['auxoMetabsMutants_m1'] = auxoMetabsMutants[m1] 
        input_data['auxoMetabsMutants_m2'] = auxoMetabsMutants[m2] 
        input_data['results_file_name_base'] = results_file_name_base

        p = Process(target = performDMMM, args = (input_data,))
        p.start()
        p.join() 
        if p.exitcode > 0:
            raise userError('**ERROR! Error in python subprocess. Please check performDMMM\n')

def process_results():
    """
    Process the results
    """
    from DMMM_moma_all import results
    
    print '\nThe total # of pairs = %s\n' % len(results)

    #---- Report based on compound concentration ---    
    # thereshold for compound concentration
    cmp_thr = 3e-5

    # Final time point 
    tf = sorted(results[0]['cell_concs'][results[0]['cell_concs'].keys()[0]].keys())[-1] 

    pairs_with_nonzero_shared_cmp = [p for p in results if len([m for m in p['shared_cmp_concs'].keys() if max(p['shared_cmp_concs'][m].values()) > cmp_thr and m != 'glc-D[e]']) > 0]
    print '\nThe total # of pairs with non-zero shared compounds = %i\n'% len(pairs_with_nonzero_shared_cmp)

    for pair in pairs_with_nonzero_shared_cmp:
        print '\n',pair['name'],'     ',
        for shared_cmp in [m for m in pair['shared_cmp_concs'].keys() if max(pair['shared_cmp_concs'][m].values()) > cmp_thr and m != 'glc-D[e]']: 
            print '(%s: %.9f) ' %(shared_cmp,max(pair['shared_cmp_concs'][shared_cmp].values())), 

    print '\nThe total # of pairs with non-zero shared compounds = %i\n'% len(pairs_with_nonzero_shared_cmp)
    print 

    #---- Report based on single cell concentration ---   
    # Initial cell concentration
    init_cell_conc = 7.5e6

    # Create a dictionary where keys and values as follows
    #   Keys: A tuple where the first element is the name of the mutant and the
    #         second element is another tuple containing the name of the mutant pair
    # Values: The cell conc at the final time point
    cell_concs_tf = dict([((m,tuple(res['cell_concs'].keys())),res['cell_concs'][m][tf]) for res in results for m in res['cell_concs'].keys()])

    # Find the mutant with maximum increase in the cell concentration
    max_conc_mutant = max(cell_concs_tf.iteritems(),key = lambda x: x[-1])
  
    max_fold_change_mutant = max_conc_mutant[1]/init_cell_conc
 
    print '\nmutant %s in mutant pair %s has the maximum fold increase %.6f\n'%(max_conc_mutant[0][0],max_conc_mutant[0][1],max_fold_change_mutant)

    #---- Report based on total cell concentration ---   
    # Create a dictionary where keys and values as follows
    #   Keys: A tuple containing the name of the mutant pair 
    # Values: The total conc of mutants at the final time point
    pair_concs_tf = dict([(tuple(res['cell_concs'].keys()),sum([res['cell_concs'][m][tf] for m in res['cell_concs'].keys()])) for res in results])

    # Wrfite into a fiile to plot with matlab 
    with open('results/moma_results_toPlot_withMatlab.m','w') as outfile:
        outfile.write('moma_fold = [\n')
        for v in pair_concs_tf.values():
            outfile.write(str(v/(2*init_cell_conc)) + '\n')
        outfile.write('];\n')

    # Find the mutant pair with maximum total concentration
    max_total_conc_pair = max(pair_concs_tf.iteritems(),key = lambda x:x[1])

    max_fold_change_pair_conc = max_total_conc_pair[1]/(2*init_cell_conc)
    print '\nMutant pair %s has the maximum fold increase in total concentration of  %.6f\n'%(max_total_conc_pair[0],max_fold_change_pair_conc)

    # Threshould for fold change in cell concentrations
    #total_cell_conc_thr = 0.4*max_fold_change_pair_conc
    total_cell_conc_thr = 1.6

    # Find all mutant pairs whose total concentrations is above the threshold
    above_thr_pairs = dict([p for p in pair_concs_tf.iteritems() if p[1]/(2*init_cell_conc) > total_cell_conc_thr])

    print '\nThe following pairs have a total concentration above the threshould:\n'
    for p in sorted(above_thr_pairs.iteritems(),key = lambda x:x[1], reverse = True):
        print '%s\t%s: %.6f\t%s: %.6f\ttotal: %.6f'%(p[0],p[0][0],cell_concs_tf[(p[0][0],p[0])]/init_cell_conc,p[0][1],cell_concs_tf[(p[0][1],p[0])]/init_cell_conc,p[1]/(2*init_cell_conc))

    print '\nThe total # of pairs whose total concentration is above the threshould (%.6f) = %i\n'%(total_cell_conc_thr,len(above_thr_pairs))

    #--- Compare with experimental data -----
    # Create a dictionary matching the pair names with and without '_Ecoli'
    no_Ecoli_name_map = []
    for p in pair_concs_tf.keys():
        no_Ecoli_name_map.append((p,(re.sub('_Ecoli','',p[0]),re.sub('_Ecoli','',p[1]))))
    no_Ecoli_name_map = dict(no_Ecoli_name_map)

    # Load the experimental growth data
    day1Rep1Inst = importData(inputFile = 'expData/day1Rep1.txt',delType = 'tab',dataType = 'float')
    day1Rep1 = day1Rep1Inst.run()

    day1Rep2Inst = importData(inputFile = 'expData/day1Rep2.txt',delType = 'tab',dataType = 'float')
    day1Rep2 = day1Rep2Inst.run()

    day4Rep1Inst = importData(inputFile = 'expData/day4Rep1.txt',delType = 'tab',dataType = 'float')
    day4Rep1 = day4Rep1Inst.run()

    day4Rep2Inst = importData(inputFile = 'expData/day4Rep2.txt',delType = 'tab',dataType = 'float')
    day4Rep2 = day4Rep2Inst.run()

    # Compute the measured fold change in growth of the pairs 
    # First find the average over replicates
    # Note that the data matrixes should actually be symmetric as (mutant1,mutant2) is the
    # same as (mutant2,mutant1), however, this is not always the case due to the experiiemntal
    # errors. As such we shouls also take average on the values of (mutant1,mutant2) and
    # (mutant2,mutant1). 
    day4Ave = dict([((m1+'_Ecoli',m2+'_Ecoli'),(day4Rep1[(m1,m2)] + day4Rep1[(m2,m1)] + day4Rep2[(m1,m2)] + day4Rep2[(m2,m1)])/4) for (m1,m2) in no_Ecoli_name_map.values()])

    # Experimental fold growth for the community. 
    expFoldGrowth = dict([(k,day4Ave[k]/(2*init_cell_conc)) for k in day4Ave.keys()])

    # Write the results into a text fiel be plotted with matlab
    with open('results/moma_results_toPlot_withMatlab.m','a') as outfile:
        outfile.write('\nexp_fold = [\n')
        for v in expFoldGrowth.values():
            outfile.write(str(v) + '\n')
        outfile.write('];\n')

    # Max fold growth for experimental data
    max_fold_growth_exp = max(expFoldGrowth.iteritems(),key = lambda x:x[1])
    print '\nThe max fold growth based experimental data is for %s and is equal to %.6f\n'%(max_fold_growth_exp[0],max_fold_growth_exp[1])

    # Threshold for cooperation accroding to experimental data
    exp_fold_growth_thr = 50

    # Mutant pairs whose fold growth (in total cell concentration) is above eight
    above_thr_pair_exp = dict([p for p in expFoldGrowth.iteritems() if p[1] >= exp_fold_growth_thr])   

    print '\nTotal # of pair with a fold growth above the threshold for experimental data: %i\n'%(len(above_thr_pair_exp.keys()))
    
    # Find how many of the preidcted cooperative mutants cooperated accroding to experimental data
    moma_exp_matches = list(set.intersection(set(above_thr_pairs.keys()),set(above_thr_pair_exp.keys())))

    print '\nMatches between cooperative phenotypes both in experimental data and moma simulaitons:\n'
    for p in moma_exp_matches:
        print '%s\tmoma %s: %.6f\tmoma %s: %.6f\tmoma total: %.6f\texp total:%.6f'%(p,p[0],cell_concs_tf[(p[0],p)]/init_cell_conc,p[1],cell_concs_tf[(p[1],p)]/init_cell_conc,above_thr_pairs[p]/(2*init_cell_conc),above_thr_pair_exp[p])
    print '\nThe total # of cooeprative mataches (thr = %.6f) = %.i\n'%(exp_fold_growth_thr,len(moma_exp_matches))


def rankAAs():
    """
    This function ranks the amino acids based on their cost
    """
    #--- E. coli iAF1260 model ---
    WT = read_gams_model(gams_model_file = '../models/Ecoli/iAF1260/iAF1260ModelData.py',model_name = 'iAF1260',organism_name = 'E. coli',model_type = 'metabolic')
    WT.biomass_reaction = WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'})
    WT.all_biomass_reactions = {'core':WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'}),'WT':WT.get_reactions({'Ec_biomass_iAF1260_WT_59p81M':'id'})}
 
    # Growth medium
    set_specific_bounds(WT,specific_bounds_file = '../models/Ecoli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

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
        WT.fba(create_model = False, store_opt_fluxes = False,screen_output = 'off')
        AAcosts[AA_exch] = WT.wildType_max_biomass - WT.fba_model.solution['objective_value'] 
        WT.get_reactions({AA_exch:'id'}).flux_bounds[0] = 0
        
    counter = 0
    for AA in [k[0] for k in sorted(AAcosts.iteritems(),key = lambda x:x[1],reverse = True)]:
        counter += 1
        print '%i\t%s\t%.6f'%(counter,AA,AAcosts[AA]) 

    # Write the results into a file and in the output 
    with open('results/moma_results_toPlot_withMatlab.m','a') as outfile:
        outfile.write('\nAAcosts = [\n')
        for AA in [k[0] for k in sorted(AAcosts.iteritems(),key = lambda x:x[1],reverse = True)]:
            outfile.write(str(AAcosts[AA]) + '\n')
        outfile.write('];\n')

        outfile.write('\nAAnames = {\n')
        for AA in [k[0] for k in sorted(AAcosts.iteritems(),key = lambda x:x[1],reverse = True)]:
            outfile.write("'" + AA + "'\n")
        outfile.write('};\n')

#------------------------------------------
if __name__ == '__main__':
    process_results()
