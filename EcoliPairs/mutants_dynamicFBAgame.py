from __future__ import division
import sys,os, time
sys.path.append('../')
sys.path.append('results/')
from copy import deepcopy
import itertools
from tools.io.read_gams_model import read_gams_model
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
from tools.fba.set_specific_bounds import set_specific_bounds
from read_exp_data import read_exp_data
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)
from multiprocessing import Process, Manager
import re

# Last updated: 07-08-2015

def performDMMM(input_data):

    t0 = input_data['t0']
    tf = input_data['tf']
    delta_t = input_data['delta_t']
    time_points = input_data['time_points'][:]
    mutant1 = deepcopy(input_data['mutant1'])
    mutant2 = deepcopy(input_data['mutant2'])
    shared_cmp_ids  = input_data['shared_cmp_ids'][:]
    shared_cmp_names  = dict(input_data['shared_cmp_names'][:])
    auxoMetabsMutants_m1  = input_data['auxoMetabsMutants_m1'][:]
    auxoMetabsMutants_m2  = input_data['auxoMetabsMutants_m2'][:]
    results_file_name_base = input_data['results_file_name_base']

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
        if exch_rxn_id not in auxoMetabsMutants_m1:
            export_rxns_m1.append(mutant1.get_reactions({exch_rxn_id:'id'})) 
        if exch_rxn_id not in auxoMetabsMutants_m2:
            export_rxns_m2.append(mutant2.get_reactions({exch_rxn_id:'id'})) 

    # Create the list of shared compounds
    shared_cmps = []

    # Glucose as a shared compound (concentration in mM)
    glucose = compound(id = 'glc-D[e]', name = 'D-Glucose', Kegg_id = 'C00031', reactant_reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],reactions = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.id == 'EX_glc(e)'],concentration = {0:111.01})    

    shared_cmps.append(glucose)

    for shared_cmp_id in [id for id in shared_cmp_ids if id != 'glc-D[e]']:
        reactant_rxns = [r for r in uptake_rxns_m1 + uptake_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
        product_rxns = [r for r in export_rxns_m1 + export_rxns_m2 if r.reactants[0].id.lower() == shared_cmp_id.lower()]
        shared_cmp = compound (id = shared_cmp_id, name = shared_cmp_names[shared_cmp_id],reactant_reactions = reactant_rxns, product_reactions = product_rxns, concentration = {0:0}) 
        shared_cmps.append(shared_cmp)

    # Set store_flux to True for all exchange reacitons related to shared metabolites
    for rxn in [r for c in shared_cmps for r in c.reactions]:
        rxn.store_flux = True
 
    for rxn in [m.biomass_reaction for m in [mutant1,mutant2]]: 
        rxn.store_flux = True

    # Get the ModelSeed ids for shared compounds
    get_cmp_ModelSeed_id(cmp_list = shared_cmps,compart_list = [c.id for c in mutant1.compartments],screen_output = 'off')

    # Get the concentration of the shared compounds in a cell pool
    cell_pool_conc(cmp_org = dict([(c,[mutant1.organism,mutant2.organism]) for c in shared_cmps]) ,screen_output = 'on')

    print '\nshared metabs ids = ',shared_cmp_ids,'\n'

    DMMM_m1m2 = DMMM_mortality(community_members = [mutant1,mutant2],shared_compounds = shared_cmps, time_range = [t0,delta_t,tf],store_dynamic_fluxes = False, screen_output = 'on')
    DMMM_m1m2.run()

    #--- Store the results ---
    # Store only the required information for each shared compound
    shared_cmps_concs_stored = {}
    for shared_cmp in [c for c in shared_cmps if c.id != 'glc-D[e]' and max(c.concentration.values()) > 1e-9]:
        shared_cmps_concs_stored[shared_cmp.id] = shared_cmp.concentration
    print '\nshared compounds with non-zero concentration = ',shared_cmps_concs_stored.keys()

    # Shared the flux of exchange (uptake) reactions for any shared metabolite, which is non-zero
    exch_rxns_stored = {}
    for rxn in list(set([r for c in shared_cmps for r in c.reactions if max([abs(v) for v in r.flux.values()]) > 1e-9])):
        exch_rxns_stored[(rxn.model.organism.id,rxn.id)] = rxn.flux

    # The following outputs are saved: shared_cmps and organism object 
    # for each community member 
    with open(results_file_name_base + '.py','a') as f:
        f.write('results.append(' + repr({'name':(mutant1.organism.id,mutant2.organism.id),'cell_concs':{mutant1.organism.id:mutant1.organism.cells_per_mL,mutant2.organism.id:mutant2.organism.cells_per_mL},'shared_cmp_concs':shared_cmps_concs_stored,'exch_rxns':exch_rxns_stored}) + ')\n')

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
        for shared_cmp_name in shared_cmps_concs_stored.keys():
            shared_cmp_concs = []
            for t in time_points:
                shared_cmp_concs.append(shared_cmps_concs_stored[shared_cmp_name][t])
            f.write('results(current_index).shared_cmp_concs.' + remove_non_alphanumeric(shared_cmp_name) + ' = ' + str(shared_cmp_concs) + ';\n')

        # Exchange reactions for shared cmps 
        for (org_id,exch_rxn_id) in exch_rxns_stored.keys():
            exch_rxn_fluxes = []
            for t in time_points:
                exch_rxn_fluxes.append(exch_rxns_stored[(org_id,exch_rxn_id)][t])
            f.write('results(current_index).exch_rxn_fluxes.' + re.sub('_LParen_e_RParen','',remove_non_alphanumeric(exch_rxn_id)) + '_' + remove_non_alphanumeric(org_id) + ' = ' + str(exch_rxn_fluxes) + ';\n')



def mutants_dynamicFBAgame(t0,delta_t,tf,start_pos,end_pos,results_file_name_base,cell_conc_init = 7.5e6):
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
    """

    results_file_name_base = results_file_name_base + '_' + str(start_pos) + '_' + str(end_pos)

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]    
    print '\ntime_points = ',time_points

    #--- E. coli iAF1260 model ---
    print '\n--- Wild-type E.coli (iAF1260 model) ----'
    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655',gDW_per_cell = 2.8e-13)

    WT = read_gams_model(file_name = '../models/Escherichia_coli/iAF1260/iAF1260ModelData.py',model_name = 'iAF1260',model_organism = model_organism,model_type = 'metabolic')

    WT.biomass_reaction = WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'})
    WT.all_biomass_reactions = {'core':WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'}),'WT':WT.get_reactions({'Ec_biomass_iAF1260_WT_59p81M':'id'})}
    WT.organism.total_death_rate = {0:0}

    # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
    # Example: EX_glc(E): glc-D[e] <==>    Vmax*C['glc-D[e]']/(Km + C['glc-D[e]']) 
    # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
    for reaction in [r for r in WT.reactions if r.type.lower() == 'exchange']:
        # The id of compound participating in the exchange reaction
        metab_id = [m.id for m in reaction.compounds][0]
        reaction.kinetics = "10*C['" + metab_id + "']/(10 + C['" + metab_id + "'])"

    # Glucose uptake kinetics 
    exch_rxns = WT.get_reactions({'EX_glc(e)':'id','EX_lys-L(e)':'id','EX_ile-L(e)':'id'})
    exch_rxns['EX_glc(e)'].kinetics = "10*C['glc-D[e]']/(0.15 + C['glc-D[e]'])"
    exch_rxns['EX_lys-L(e)'].kinetics = "0.1964*C['lys-L[e]']/(5e-4 + C['lys-L[e]']) + 0.3055*C['lys-L[e]']/(1e-2 + C['lys-L[e]'])"
    exch_rxns['EX_ile-L(e)'].kinetics = "0.0346*C['ile-L[e]']/(1.22e-3 + C['ile-L[e]'])"
   
    # Growth medium
    set_specific_bounds(WT,file_name = '../models/Escherichia_coli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign the objective function coefficients
    for rxn in WT.reactions:
        rxn.objective_coefficient = 0
    WT.biomass_reaction.objective_coefficient = 1

    # Perform FBA for the wild-type
    WT.fba(assign_wildType_max_biomass = True)

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

    # The id of all compounds in the cell pool
    cell_pool_cmp_ids = ['ala-L[e]','val-L[e]','gly[e]','ile-L[e]','thr-L[e]','leu-L[e]','ser-L[e]','pro-L[e]','asp-L[e]','cys-L[e]','met-L[e]','glu-L[e]','phe-L[e]','tyr-L[e]','orn[e]','lys-L[e]','trp-L[e]','arg-L[e]']

    shared_cmp_ids += cell_pool_cmp_ids

    # Initialize the file
    with open(results_file_name_base + '.py','w') as f:
        f.write('results = []\n')
    with open(results_file_name_base + '.m','w') as f:
        f.write('results = struct([]);\n')

    #--- DMMM for the co-culture of mutant1 and mutant2 mutants ---
    #for (m1,m2) in [('lysA','ilvE')]:
    #for (m1,m2) in [('glnA','trpC'),('lysA','ilvE')]:
    for (m1,m2) in mutant_pairs[start_pos - 1:end_pos]:
      
        print '\n**** %i. (%s,%s) ****\n' % (mutant_pairs.index((m1,m2)) + 1,m1,m2)
    
        #-- Mutant 1 and the related reactions whose flux should be set to zero --
        print '\n-- ' + m1 + '_Ecoli mutant --'
        mutant1 = deepcopy(WT)
        mutant1.id = WT.id + '_' + m1
        mutant1.organism.id = m1 + '_Ecoli'
        mutant1.organism.cells_per_mL = {0:per_capita_cell_conc_init}
        mutant1.organism.gDW_per_mL = {0:per_capita_cell_conc_init*mutant1.organism.gDW_per_cell}
        for rxn_id in mutants_rxn_info[m1]:
            rxn = mutant1.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]

        mutant1.fba(create_model = False, store_opt_fluxes = False)

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
        print '\n-- ' + m2 + '_Ecoli mutant --'
        mutant2 = deepcopy(WT)
        mutant2.id = WT.id + '_' + m2
        mutant2.organism.id = m2 + '_Ecoli'
        mutant2.organism.cells_per_mL = {0:per_capita_cell_conc_init}
        mutant2.organism.gDW_per_mL = {0:per_capita_cell_conc_init*mutant2.organism.gDW_per_cell}
        for rxn_id in mutants_rxn_info[m2]:
            rxn = mutant2.get_reactions({rxn_id:'id'})
            rxn.flux_bounds = [0,0]
            
        mutant2.fba(create_model = False, store_opt_fluxes = False)
    
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
    
        # --- Define the compounds available in the extracellular medium ---
        # Consider only unique elements as, in general, it is quite possible that
        # two mutants need the same compounds to survive
        shared_cmp_ids = sorted(list(set(shared_cmp_ids)))

        # Get the names of the shared compounds from the SBML file
        shared_cmp_names = get_cmp_name_from_SBML(cmp_id_list = shared_cmp_ids)
        shared_cmp_names['tlc-D[e]'] = 'D-glucose'
        shared_cmp_names['ala-L[e]'] = 'L-Alanine'
        shared_cmp_names['val-L[e]'] = 'L-valine'
        shared_cmp_names['gly[e]'] = 'glycine'
        shared_cmp_names['ile-L[e]'] = 'L-isoleucine'
        shared_cmp_names['thr-L[e]'] = 'L-threonine'
        shared_cmp_names['leu-L[e]'] = 'L-leucine'
        shared_cmp_names['ser-L[e]'] = 'L-serine'
        shared_cmp_names['pro-L[e]'] = 'L-proline'
        shared_cmp_names['asp-L[e]'] = 'L-aspartate'
        shared_cmp_names['cys-L[e]'] = 'L-cysteine'
        shared_cmp_names['met-L[e]'] = 'L-methionine'
        shared_cmp_names['glu-L[e]'] = 'L-glutamate'
        shared_cmp_names['phe-L[e]'] = 'L-phenylalanine'
        shared_cmp_names['tyr-L[e]'] = 'L-tyrosine'
        shared_cmp_names['orn[e]'] = 'Ornithine'
        shared_cmp_names['lys-L[e]'] = 'L-lysine'
        shared_cmp_names['trp-L[e]'] = 'L-tryptophan'
        shared_cmp_names['arg-L[e]'] = 'L-arginine'
 
        # Creating a shared memory using the manager
        input_data = {} 
        input_data['t0'] = t0
        input_data['tf'] = tf
        input_data['delta_t'] = delta_t
        input_data['time_points'] = time_points
        input_data['mutant1'] = mutant1
        input_data['mutant2'] = mutant2
        input_data['shared_cmp_ids'] = shared_cmp_ids 
        input_data['shared_cmp_names'] = shared_cmp_names.items() 
        input_data['auxoMetabsMutants_m1'] = auxoMetabsMutants[m1] 
        input_data['auxoMetabsMutants_m2'] = auxoMetabsMutants[m2] 
        input_data['results_file_name_base'] = results_file_name_base

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
    """
    # Create a vector of initial concentrations
    per_capita_cell_conc_init_vec = [7.5e6,7.5e7,7.5e8,7.5e9,7.5e10]
    for per_capita_cell_conc_init in per_capita_cell_conc_init_vec:
        mutants_mortality(t0 = 0,delta_t = 0.5,tf = 96,start_pos = 652,end_pos = 652,results_file_name_base = 'results/init_conc_effect_' + str(per_capita_cell_conc_init/7.5),per_capita_cell_conc_init = per_capita_cell_conc_init)

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

    pairs_with_nonzero_shared_cmp = [p for p in results if len([m for m in p['shared_cmp_concs'].keys() if max(p['shared_cmp_concs'][m].values()) > metab_thr and m != 'glc-D[e]']) > 0]
    print '\nThe total # of pairs with non-zero shared compounds = %i\n'% len(pairs_with_nonzero_shared_cmp)

    for pair in pairs_with_nonzero_shared_cmp:
        print '\n',pair['name'],'     ',
        for shared_cmp in [m for m in pair['shared_cmp_concs'].keys() if max(pair['shared_cmp_concs'][m].values()) > metab_thr and m != 'glc-D[e]']: 
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
    cell_pool_concs = {'EX_ala-L(e)':6.7e-3,'EX_val-L(e)':9.2e-3,'EX_gly(e)':1.7e-3,'EX_ile-L(e)':1.6e-3,'EX_thr-L(e)':0.2e-3,'EX_leu-L(e)':0.4e-3,'EX_ser-L(e)':0,'EX_pro-L(e)':0.4e-3,'EX_asp-L(e)':0,'EX_cys-L(e)':1.8e-3,'EX_met-L(e)':0.3e-3,'EX_glu-L(e)':25.5e-3,'EX_ph-L(e)':6.8e-3,'EX_tyr-L(e)':0.1e-3,'EX_orn-L(e)':0,'EX_lys-L(e)':0.2e-3,'EX_trp-L(e)':0,'EX_arg-L(e)':0}

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


def create_pbs_files(t0,delta_t,tf,total_comb_num,interval_size, outfile_base_name,results_file_name_base,with_output = 0):
    """
    Creates the pbs files given:
        total_comb_num: Total number of combinations 
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
    
                # python -c "from mutants_mortality import mutants_mortality;mutants_mortality(t0 = 0, delta_t = 0.5, tf = 96, start_pos = 1, end_pos = 100,results_file_name_base = 'results/DMMM_mortality')" 
                outfile.write('python -c "from mutants_mortality import mutants_mortality;mutants_mortality(t0 = ' + str(t0) + ', delta_t = ' + str(delta_t) + ', tf = ' +str(tf) + ', start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", results_file_name_base = '" + results_file_name_base + "')\"\n\n")

                outfile.write("python -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n\n")

            elif with_output == 1:
                outfile.write("python -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\" > " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + ".out 2>&1\n\n")

                # python -c "from mutants_mortality import mutants_mortality;mutants_mortality(t0 = 0, delta_t = 0.5, tf = 96, start_pos = 1, end_pos = 100,results_file_name_base = 'results/DMMM_mortality')" >> job_DMMM_mortality_1_100.out 2>&1 
                outfile.write('python -c "from mutants_mortality import mutants_mortality;mutants_mortality(t0 = ' + str(t0) + ', delta_t = ' + str(delta_t) + ', tf = ' +str(tf) + ', start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", results_file_name_base = '" + results_file_name_base + "')\" >> " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + '.out 2>&1\n\n')

                outfile.write("python -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\" >> " + outfile_base_name + '_' + str(slice[0]) + '_' + str(slice[1]) + ".out 2>&1\n\n")
    
            else:
                raise userError('Invalid with_output value.')

            outfile.close()
          
            # make it executable
            os.system('chmod u+x ' + file_name)

#-------------------------
if __name__ == '__main__':
    process_results(metab_thr = 1e-10,total_cell_conc_frac = 0.999998, exp_fold_growth_thr = 50)
