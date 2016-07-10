from __future__ import division
import sys,os, time
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
sys.path.append('../')
sys.path.append('results/')
from copy import deepcopy
import itertools
import numpy as np
from tools.globalVariables import *
from tools.userError import *
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.io.create_model import create_model
from tools.GAMETES.game import game
from tools.GAMETES.replicator_dynamics import replicator_dynamics
from tools.utilities.load_data_fromFile import load_data_from_python_file
import cobra
from read_exp_data import read_exp_data
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)
from multiprocessing import Process, Manager
import re

"""
This script investigates the type of games played for different values of the amino acids leakiness (capture efficiencies) 
"""

def gameCreator_BQH(input_data, output_data):
    """
    Creates the game and finds its Nash equilibrium for a single case

    Here, we have four strains:
    1. Mutant m1, which cannot produce amino acid 1 but synthesizes
       the leaky amino acid 2 required by mutant 2 (Cooperator m1)
    2. Mutant m1, which cannot produce amino acid 1 and does not 
       synthesize the leaky amino acid 2 required by mutant 2 (Cheater
       m1)
    4. Mutant m2, which cannot produce amino acid 2 but synthesizes
       the leaky amino acid 1 required by mutant 1 (Cooperator m1)
    5. Mutant m2, which cannot produce amino acid 2 and does not 
       synthesize the leaky amino acid 1 required by mutant 1 (Cheater
       m2)
    """
    strains = input_data['strains'] 
    growthMedium_flux_bounds = input_data['growthMedium_flux_bounds']
    death_rate = input_data['death_rate']
    coopr_exchrxns_fluxRanges = input_data['coopr_exchrxns_fluxRanges']
    max_game_players_num = input_data['max_game_players_num']
    WT_can_produce_all_leakyTraits = input_data['WT_can_produce_all_leakyTraits']
    t0 = input_data['t0']
    tf = input_data['tf']
    dt = input_data['dt']
    exclude_infeas_fba = input_data['exclude_infeas_fba']
    simulate_rep_dynamics = input_data['simulate_rep_dynamics']
    results_filename = input_data['results_filename'] 
    save_details = input_data['save_details'] 
    warnings = input_data['warnings']
    stdout_msgs = input_data['stdout_msgs'] 
    stdout_msgs_details = input_data['stdout_msgs_details'] 
    stdout_msgs_fba = input_data['stdout_msgs_fba'] 

    #--- Check2 whether the wild-type can satisfy the imposed leakienes levels --
    set_specific_bounds(model = strains['wild_type'], flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    
    set_specific_bounds(model = strains['wild_type'], flux_bounds = dict([(exchrxn_id,[strains['wild_type'].leaky_traits[rxn_set]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for rxn_set in strains['wild_type'].leaky_traits.keys() for exchrxn_id in rxn_set]), reset_flux_bounds = False)    
    if stdout_msgs_details:
        print 'leakiness levels = ', dict([(exchrxn_id,[strains['wild_type'].leaky_traits[rxn_set]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for rxn_set in strains['wild_type'].leaky_traits.keys() for exchrxn_id in rxn_set])
    strains['wild_type'].fba(stdout_msgs = stdout_msgs_fba)

    # A dictionary where keys are the number of players considered (number_of_players) and values are another dictionary 
    # with keys and values as follows:
    # Key: 'Nash eq' --> lists of lists containing the Nash equilibria for each given number of players
    # Key: 'payoff matrix' --> A dictionary containing the payoff matrix fo the game 
    games_info = {}

    # If max_game_players_num set it to the maximum possible, i.e., the total number of strians
    if max_game_players_num == None:
        max_game_players_num = len(strains.keys())

    for number_of_players in range(2,max_game_players_num + 1):
        # Here, instead of players, strategies. But since class "game" requires
        # player names we just use player1 and player2 as the player names
        # The strategies for each player is either to produce everything (wild-
        # type) or to not produce one or more leaks traits
        players_names = ['player' + str(i) for i in range(1,number_of_players + 1)]
        payoff_matrix = {}
        players_strategies = {}
        for player_name in players_names:
            players_strategies[player_name] = strains.keys()
    
        # ------ All possible combinations of strategies --------
        for strain_names_tuple in itertools.combinations_with_replacement(strains.keys(), r = number_of_players):
            if stdout_msgs:
                print '\n--- Compute payoffs for {}'.format(strain_names_tuple)
    
            """
            First, we need to find out what leaky traits are common among two or
            more species and then find out how much benefit is provided to each
            strain in the community. Consider a given leaky trait. Let N be the
            total number of strains in the community, Nt be the number of strains
            with this leaky trait in the community and x be the leakiness level 
            for this trait and A be the leaky compound produced. Then,
            Net export of A for a strain producing and leaking A =
                x - (Nt - 1)*x/(N - 1)
            Net uptake of A for a strian lacking this leaky trait (does not
            produce and leak A) = 
                Nt*x/(N - 1)
    
            Examples: Let the second trait in the following cases be the one 
                      producing compound A
                      (a) 1 and 1 --> N = 2, Nt = 2
                          vExport(A,1) = x - (2 - 1)*x/(2 - 1) = x - x = 0
                      (b) 11 and 11 --> N = 2 ,  Nt = 2
                          vExport(A,11) = x - (2 - 1)*x/(2-1) = x - x = 0
                      (c) 11 and 10 --> N = 2 , Nt = 1
                          vExport(A,11) = x - (1 - 1)*x/(2-1) = x
                          vUptake(A,10) = 1*x/(2-1) = x
                      (d) 11 and 11 and 11 --> N = 3, Nt = 3
                          vExport(111,A) = x - (3 - 1)*x/(3 - 1) = 0
                      (e) 11 and 11 and 10 --> N = 3, Nt = 2
                          vExport(11,A) = x - (2 - 1)*x/(3-1) = x/2
                          vUptake(10,A) = 2*x(3-1) = x 
                      (f) 11 and 10 and 10 -- N = 3, Nt = 1
                          vExport(11,A) = x - (1 - 1)*x/(3-1) = 0
                          vUptake(10,AA) = 1*x/(3-1) = x/2 
            """
    
            # A dictionary with keys being the strain names and values being 
            # their payoff
            strains_payoffs = dict([(strain_name, None) for strain_name in strains.keys()]) 
    
            # Strains names whose FBA problem turned out to be infeasible in any run of the
            # while loop
            infeasible_strain_names = []

            done = False    
            # This outer while loop is to redo the computations of payoff if the FBA problem 
            # for one or more community member turns out to be infeasible 
            while not done: 
     
                #--- Find out what strains produce each leaky trait ---
                # This is a dictionary with keys being reaction ids and values 
                # being a list of strain names having that leaky strait
                leaky_traits_producers = dict([(exchrxn_id,[]) for exchrxn_id in set([rxn for strain_name in strain_names_tuple if strain_name not in infeasible_strain_names for rxn_set in strains[strain_name].leaky_traits.keys() for rxn in rxn_set])])
                for strain_name in [s for s in strain_names_tuple if s not in infeasible_strain_names]:
                    for  exchrxn_id in list(set([r for rset in strains[strain_name].leaky_traits.keys() for r in rset])):
                        leaky_traits_producers[exchrxn_id].append(strain_name)
                for exchrxn_id in leaky_traits_producers.keys():
                   leaky_traits_producers[exchrxn_id] = list(leaky_traits_producers[exchrxn_id]) 
                if stdout_msgs_details:
                    print 'leaky_traits_producers',leaky_traits_producers 
        
                # Strains names whose FBA problem turned out to be infeasible for the current
                # run of the while loop
                infeasible_strain_names_curr = []
        
                # Perform FBA for any strain whose FBA was NOT infeasible
                # We use set(strain_names_tuple) here in order to avoid redoing 
                # FBA for a strain that may appear more than once in a 
                # particular combinations, e.g., 11 vs 11
                for strain_name in [s for s in set(strain_names_tuple) if s not in infeasible_strain_names]:
                    if stdout_msgs:
                        print '- FBA for {}'.format(strain_name)

                    # Growth condition
                    set_specific_bounds(model = strains[strain_name], flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    
        
                    # Knockouts 
                    set_specific_bounds(model = strains[strain_name], flux_bounds = dict([(rxn_id,[0,0]) for rxn_id in list(itertools.chain(*strains[strain_name].knockedout_rxn_ids))]), reset_flux_bounds = False)    
        
                    # Leaky traits. If leaky_traits = {('r1','r2'):0.1, 
                    # ('r3','r4'):0.2, ('r5',):0.3}, then flux_bounds = 
                    # [('r1', [0.1*max_coopr_level, None]), 
                    # ('r2', [0.1*max_coopr_level, None]), 
                    # ('r3', [0.2*max_coopr_level, None]), 
                    # ('r4', [0.2*max_coopr_level, None]), 
                    # ('r5', [0.3*max_coopr_level, None])] 
                    # Net leakiness (export) = x - (Nt - 1)*x/(N - 1)
                    net_leakiness = dict([(exchrxn_id, strains[strain_name].leaky_traits[rset] - (len(leaky_traits_producers[exchrxn_id]) - 1)*strains[strain_name].leaky_traits[rset]/(len(strain_names_tuple) - 1)) for rset in strains[strain_name].leaky_traits.keys() for exchrxn_id in rset])
                    set_specific_bounds(model = strains[strain_name], flux_bounds = dict([(exchrxn_id,[net_leakiness[exchrxn_id]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for exchrxn_id in net_leakiness.keys()]), reset_flux_bounds = False)    
                    if stdout_msgs_details:
                        print '{} net leaky traits (export) levels = {}\n'.format(strain_name, net_leakiness)
                
                    # Benefit from the leaky traits of the partner for traits this 
                    # strain is not leaky (compounds it is not synthesizing)
                    # 
                    # All exchrxn ids for compounds this strain is not leaky for
                    # (i.e., those it does not produce)
                    not_leaky_exchrxn_ids = [exchrxn_id for exchrxn_id in leaky_traits_producers.keys() if exchrxn_id not in [r for rset in strains[strain_name].leaky_traits.keys() for r in rset]]
                    # The leakiness level for each exchrxn in not_leaky_exchrxn_ids 
                    # is the same for all strains that leaky for this trait.
                    # it's enough just to pick one and find its leakiness level
                    partners_leakiness_levels = dict([(exchrxn_id, strains[leaky_traits_producers[exchrxn_id][0]].leaky_traits[rset]) for exchrxn_id in not_leaky_exchrxn_ids for rset in strains[leaky_traits_producers[exchrxn_id][0]].leaky_traits.keys() if exchrxn_id in rset])
                    if stdout_msgs_details:
                        print 'partners_leakiness_levels = ',partners_leakiness_levels
                    # Net uptake = Nt*x/(N - 1)
                    net_uptakes = dict([(exchrxn_id, len(leaky_traits_producers[exchrxn_id])*partners_leakiness_levels[exchrxn_id]/(len(strain_names_tuple) - 1)) for exchrxn_id in not_leaky_exchrxn_ids]) 
                    if stdout_msgs_details:
                        print 'Net uptakes % for {} = {}'.format(strain_name, net_uptakes)
                        print 'Net uptakes for {} = {}'.format(strain_name, [(exchrxn_id,[-net_uptakes[exchrxn_id]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for exchrxn_id in net_uptakes.keys()])
                    set_specific_bounds(model = strains[strain_name], flux_bounds = dict([(exchrxn_id,[-net_uptakes[exchrxn_id]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for exchrxn_id in net_uptakes.keys()]), reset_flux_bounds = False)    
               
                    # fba
                    strains[strain_name].fba(build_new_optModel = False, stdout_msgs = stdout_msgs_fba)
                    if strains[strain_name].fba_model.solution['exit_flag'] == 'globallyOptimal':
                        strains_payoffs[strain_name] = strains[strain_name].fba_model.solution['objective_value']
                    else: 
                        strains_payoffs[strain_name] = death_rate 
                        infeasible_strain_names.append(strain_name) 
                        infeasible_strain_names_curr.append(strain_name) 
                        # Uf infeasible, remove it from leaky_traits_producers
                        #for leaky_trait in [lt for lt in leaky_traits_producers.keys() if strain_name in leaky_traits_producers[lt]]:
                        #    # Index of strain name in leaky_traits_producers[leak_trait]
                        #    ind = leaky_traits_producers[leaky_trait].index(strain_name)    
                        #    del leaky_traits_producers[leaky_trait][ind]                   
                        #if stdout_msgs_details:
                        #    print '{} was removed form leaky_traits_producers was updated. leaky_traits_producers = {}'.format(strain_name, leaky_traits_producers)

                # Stop if there is no infeasible FBA problem or if all FBA problems were
                # infeasible    
                if exclude_infeas_fba:
                    if len(infeasible_strain_names_curr) == 0 or len([s for s in set(strain_names_tuple) if s not in infeasible_strain_names]) == 0:
                        done = True    
                else:
                    done = True    

            # Assign values of the payoff matirx
            for strain_permutation in itertools.permutations(strain_names_tuple , r = len(strain_names_tuple)):
                payoff_matrix_key = tuple([(p,strain_permutation[players_names.index(p)]) for p in players_names])
                payoff_matrix_value = dict([(p,strains_payoffs[strain_permutation[players_names.index(p)]]) for p in players_names])
                payoff_matrix[payoff_matrix_key] = payoff_matrix_value 
            #if stdout_msgs_details:
            #    print '\nstrain_names_tuple = {}, payoff_matrx_key = {}, payoff_matrix_value = {}'.format(strain_names_tuple, payoff_matrix_key, payoff_matrix_value)
    
        # Find Nash equilibria
        BQH_game = game(game_name = 'BQH Ecoli game', players_names = players_names, players_strategies = players_strategies, payoff_matrix = payoff_matrix)
        BQH_game.find_NashEq(stdout_msgs = False)
   
        # Find the payoff matrix of the symmetric game (to save in file)
        BQH_game.create_symmetric_payoff_matrix()
 
        output_data['BQH_game'] = BQH_game
    
        games_info[number_of_players] = {'NashEq':list(set([tuple(sorted(neq)) for neq in [dict(e).values() for e in BQH_game.pureNash_equilibria]])), 'payoff matrix':BQH_game.symmetric_payoff_matrix}

        if stdout_msgs_details:
            for number_of_players in [n for n in games_info.keys() if n < 4]:
                print '\n--- number of players = {} ---\n'.format(number_of_players)
                print '\nPayoff matrix of the games: '
                for k in games_info[number_of_players]['payoff matrix'].keys():
                    print k,' --> ',games_info[number_of_players]['payoff matrix'][k]

    #--- Simulate replicator's dynamics ---
    if simulate_rep_dynamics:
        input_data_repdyn = {}
        input_data_repdyn['strains_names'] = strains.keys()
        input_data_repdyn['all_leaky_traits'] = tuple(set([leaky_trait for strain in strains.values() for leaky_trait in strain.leaky_traits.items()])) 
        input_data_repdyn['games_info'] = games_info 
        input_data_repdyn['WT_can_produce_all_leakyTraits'] = WT_can_produce_all_leakyTraits 
        input_data_repdyn['t0'] = t0
        input_data_repdyn['tf'] = tf
        input_data_repdyn['dt'] = dt
        input_data_repdyn['save_details'] = save_details
        input_data_repdyn['results_filename'] = results_filename 
        input_data_repdyn['stdout_msgs'] = stdout_msgs 

        output_data_repdyn = Manager().dict()
        output_data_repdyn['replicator_dynamics_info'] = None
        
        p = Process(target = do_replicator_dynamics, args = (input_data_repdyn, output_data_repdyn))
        p.start()
        p.join() 
        if p.exitcode > 0:
            raise userError('Error in python subprocess. Please check do_replicator_dynamics\n')
        else:
            replicator_dynamics_info = output_data_repdyn['replicator_dynamics_info']

    else:
        if not save_details:
            for gk in games_info.keys():
                del games_info[gk]['payoff matrix']

        replicator_dynamics_info = {}
        results = {'WT_can_produce_all_leakyTraits': WT_can_produce_all_leakyTraits, 'games_info':games_info, 'replicator_dynamics_info':replicator_dynamics_info}
        if results_filename != '':
            all_leaky_traits = tuple(set([leaky_trait for strain in strains.values() for leaky_trait in strain.leaky_traits.items()])) 
            with open(results_filename,'a') as f:
                f.write('results[' + str(all_leaky_traits) + '] = ' + str(results) + '\n')


def do_replicator_dynamics(input_data_repdyn, output_data_repdyn):
    """
    Performs replicator's dynamic simulations
    """
    strains_names = input_data_repdyn['strains_names'] 
    games_info = input_data_repdyn['games_info'] 
    all_leaky_traits = input_data_repdyn['all_leaky_traits']
    WT_can_produce_all_leakyTraits = input_data_repdyn['WT_can_produce_all_leakyTraits'] 
    t0 = input_data_repdyn['t0'] 
    tf = input_data_repdyn['tf'] 
    dt = input_data_repdyn['dt'] 
    stdout_msgs = input_data_repdyn['stdout_msgs'] 
    results_filename = input_data_repdyn['results_filename'] 
    save_details = input_data_repdyn['save_details'] 

    replicator_dynamics_info = {}

    x_init = dict([(strain_name,0) for strain_name in strains_names])

    # Check if any of the mutatns can invade wild-type
    x_init['wild_type'] = 0.99
    for strain_name in [k for k in x_init.keys() if k != 'wild_type']:
        x_init[strain_name] = 0.01/(len(strains_names) - 1)
    strains_fracs = replicator_dynamics(payoff_matrix = dict([(k,v) for number_of_players in games_info.keys() for (k, v) in games_info[number_of_players]['payoff matrix'].iteritems() for number_of_players in games_info.keys()]), x_init = x_init, time_range = [t0,dt,tf])   
    # Strain fractions in the last time point
    replicator_dynamics_info[tuple(x_init.items())] = dict([(strain_name, strains_fracs[strain_name][tf]) for strain_name in x_init.keys()])

    # Check if wild-type can invade the population of mutants 
    x_init['wild_type'] = 0.01
    for strain_name in [k for k in x_init.keys() if k != 'wild_type']:
        x_init[strain_name] = 0.99/(len(strains_names) - 1)
    strains_fracs = replicator_dynamics(payoff_matrix = dict([(k,v) for number_of_players in games_info.keys() for (k, v) in games_info[number_of_players]['payoff matrix'].iteritems() for number_of_players in games_info.keys()]), x_init = x_init, time_range = [t0,dt,tf])   
    # Strain fractions in the last time point
    replicator_dynamics_info[tuple(x_init.items())] = dict([(strain_name, strains_fracs[strain_name][tf]) for strain_name in x_init.keys()])

    output_data_repdyn['replicator_dynamics_info'] = replicator_dynamics_info

    if stdout_msgs:
        print "\n-- Strain fractions according to Replicator's dynamics ---\n"
        for x0 in replicator_dynamics_info.keys():
            print '\nx_init = {} --> x_tf = {}\n'.format(x0,replicator_dynamics_info[x0])

        for number_of_players in games_info.keys():
            print '\n--- number of players = {} ---\n'.format(number_of_players)

            print '\nNash equilibria:'
            for neq in games_info[number_of_players]['NashEq']:
                print neq

    if results_filename != '':
        if not save_details:
            for gk in games_info.keys():
                del games_info[gk]['payoff matrix']

        results = {'WT_can_produce_all_leakyTraits': WT_can_produce_all_leakyTraits, 'games_info':games_info, 'replicator_dynamics_info':replicator_dynamics_info}

        with open(results_filename,'a') as f:
            f.write('results[' + str(all_leaky_traits) + '] = ' + str(results) + '\n')

def find_coopr_exchrxns_fluxRanges(analysis_type = 'AAs', results_filename = 'coopr_exchrxns_fluxRanges_AAs.py', warnings = True, stdout_msgs = True):
    """
    Finds the min and max value of each compound that can be produced by the wild-type strain

    NOTE:
    Since the max flux value for EX_ala_L(e) and EX_met_L(e) were zero, we added the following two
    reactions to the model:
    METt2rpp        1.0 met_L_p + 1.0 h_p <==> 1.0 met_L_c + 1.0 h_c        (rxn manually added to the model)
    GLNt2rpp        1.0 gln_L_p + 1.0 h_p <==> 1.0 gln_L_c + 1.0 h_c        (rxn manually added to the model)
    """
    if analysis_type == 'AAs':
        #mutants_auxotrophy = load_data_from_python_file(file_name = 'mutants_auxotrophy_AAs_all.py', var_names = ['mutants_auxotrophy'])['mutants_auxotrophy']
        mutants_auxotrophy = load_data_from_python_file(file_name = 'mutants_auxotrophy_AAs.py', var_names = ['mutants_auxotrophy'])['mutants_auxotrophy']

        # Exchange reactions needed for cooperaiton in all pairs
        cooperative_exchrxn_ids = list(set([r for m in mutants_auxotrophy.keys() for r_list in mutants_auxotrophy[m] for r in r_list]))

    elif analysis_type == 'PamSilver':
        from mutants_auxotrophy_PamSilver import auxoMetabsMutants

        # Exchange reactions needed for cooperaiton in all pairs
        cooperative_exchrxn_ids = list(set([r for m in auxoMetabsMutants.keys() for r_list in auxoMetabsMutants[m] for r in r_list]))
  

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds) 

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

def simulate_rep_dynamics(game_results_filename, start_pos = None, end_pos = None,  t0 = 0, tf = 1000, dt = 0.5, save_details = False, results_filename = '', stdout_msgs = True):
    """ 
    This function gets the game results and simualtes the replicator's dynamics
    """ 
    game_results = load_data_from_python_file(file_name = 'results/leakiness_games_BQH_2LT_all_old.py', var_names = ['results'])['results']

    # Initialize the file
    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('results = {}\n')

    if start_pos != None and end_pos != None:
        cases_to_consider = game_results.keys()[start_pos - 1:end_pos]
        counter = start_pos - 1
    else: 
        cases_to_consider = game_results.keys()
        counter = 0 

    for res_key in cases_to_consider:
    #for res_key in game_results.keys()[:2]:
    #for res_key in [((('EX_val_L(e)',), 0.1), (('EX_ala_L(e)',), 0.1))]:
        counter += 1
        if counter/500 == int(counter/500) or counter == len(game_results.keys()):
            print '{} cases considered'.format(counter)

        input_data_repdyn = {}
        input_data_repdyn['strains_names'] = list(set([n for gk in game_results[res_key]['games_info'][2]['payoff matrix'].keys() for n in gk]))
        input_data_repdyn['all_leaky_traits'] = res_key 
        input_data_repdyn['games_info'] = game_results[res_key]['games_info'] 
        input_data_repdyn['WT_can_produce_all_leakyTraits'] = game_results[res_key]['WT_can_produce_all_leakyTraits'] 
        input_data_repdyn['t0'] = t0
        input_data_repdyn['tf'] = tf
        input_data_repdyn['dt'] = dt
        input_data_repdyn['save_details'] = save_details
        input_data_repdyn['results_filename'] = results_filename 
        input_data_repdyn['stdout_msgs'] = stdout_msgs 
    
        output_data_repdyn = Manager().dict()
        output_data_repdyn['replicator_dynamics_info'] = None
        
        p = Process(target = do_replicator_dynamics, args = (input_data_repdyn, output_data_repdyn))
        p.start()
        p.join() 
        if p.exitcode > 0:
            raise userError('Error in python subprocess. Please check do_replicator_dynamics\n')
        else:
            replicator_dynamics_info = output_data_repdyn['replicator_dynamics_info']
    

def analyze_AAs(results_filename = '', stdout_msgs = True, warnings = True):
    """
    Performs the following tasks:
    (i) whether removing the reactions corresponding to gene knockotus given in mutants_rxn_info_AAs 
        (see mutants_rxn_info_iJO1366.py) lead to zero growth rate)
    (ii) Finds the compounds each mutant is auxotroph for

    INPUTS:
    -------
    results_filename: 
    File name storing the results
    """
    from tools.fba.auxotrophy_finder import auxotrophy_finder

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds) 

    max_biomass = WT.fba_model.solution['objective_value']

    # Load the list of rxns that must be off in each mutant
    from mutants_rxn_info_iJO1366 import mutants_rxn_info_AAs, genes_AA_map

    mutants_auxotrophy = {}
    mutants_exitflags = {}

    forbidden_cpds_exchrxn_ids = ['EX_o2(e)']

    #mutants_rxn_info_AAs['ilvE'] = ['KARA1', 'KARA2']
    #mutants_rxn_info_AAs['alaA_alaC_avtA__ilvE'] = ['ALATA_L', 'VPAMTr', 'ILETA', 'VALTA']
    #genes_AA_map['alaA_alaC_avtA__ilvE'] = 'ala_L_ile_L' 
    #mutants_rxn_info_AAs['ilvE_avtA__ilvE'] = ['VPAMTr', 'ILETA', 'VALTA']
    #genes_AA_map['ilvE_avtA__ilvE'] = 'val_L_ile_L' 
    
    for gen in sorted(mutants_rxn_info_AAs.keys()):
    #for gen in ['alaA_alaC_avtA', 'ilvE','alaA_alaC_avtA__ilvE']: 
    #for gen in ['ilvE_avtA', 'ilvE','ilvE_avtA__ilvE']: 
        if stdout_msgs:
            print '\n{}, {}'.format(gen, genes_AA_map[gen])

        set_specific_bounds(model = WT, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)

        # Impose rxn knockouts
        set_specific_bounds(model = WT, flux_bounds = dict([(r,[0, 0]) for r in mutants_rxn_info_AAs[gen]]), reset_flux_bounds = False)
        WT.fba(stdout_msgs = False)
        fba_soln = WT.fba_model.solution['objective_value']

        # Test
        #set_specific_bounds(model = WT, flux_bounds = {'EX_gly(e)':[-10,None]}, reset_flux_bounds = False)
        #WT.fba(stdout_msgs = True)
        #set_specific_bounds(model = WT, flux_bounds = {'EX_gly(e)':[0,None]}, reset_flux_bounds = False)
        #cpd = WT.compounds_by_id['gly_c']
        #for rxn in [r for r in cpd.reactions if (cpd in r.products and r.reversibility.lower() == 'irreversible') or (cpd in r.compounds and r.reversibility.lower() == 'reversible')]:
        #    print '{}:\t{}\t{}'.format(rxn.id,rxn.get_equation(),WT.fba_model.solution['opt_rxnFluxes'][rxn.id])
        #print '\n'

        # Find out for what compounds this mutant is auxotroph for
        auxotrophy_finder_inst = auxotrophy_finder(model = WT, mutant_name = gen, max_biomass = max_biomass, viability_thr = 0.01, max_soln_size = 0, max_soln_num = 20, max_cpds_uptake_flux = 1000, forbidden_cpds_exchrxn_ids = forbidden_cpds_exchrxn_ids, results_filename = '', warnings = True, stdout_msgs = False)
        (exit_flag, auxotroph_cpds) = auxotrophy_finder_inst.run() 
        mutants_auxotrophy[gen] = auxotroph_cpds
        mutants_exitflags[gen] = exit_flag
          
        if stdout_msgs:
            print 'fba_soln = {}\nauxotroph_cpds = {}\nexit flag = {}\n'.format(fba_soln, auxotroph_cpds, exit_flag)

    """
    #--- Find auxotrohps for serC ---
    mutants_auxotrophy[gen] = []
    set_specific_bounds(model = WT, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)
    # Impose rxn knockouts
    set_specific_bounds(model = WT, flux_bounds = dict([(r,[0, 0]) for r in mutants_rxn_info_AAs[gen]]), reset_flux_bounds = False)
    WT.fba(stdout_msgs = False)
    fba_soln = WT.fba_model.solution['objective_value']
    for exchrxn in [r.id for r in WT.reactions if r.is_exchange and r.flux_bounds[0] >= 0 and r.id not in forbidden_cpds_exchrxn_ids]:
        set_specific_bounds(model = WT, flux_bounds = {exchrxn:[-10,None]}, reset_flux_bounds = False)
        WT.fba(stdout_msgs = False)
        if WT.fba_model.solution['exit_flag'] == 'globallyOptimal' and WT.fba_model.solution['objective_value'] >= 0.01*max_biomass:
             mutants_auxotrophy['serC'].append(exchrxn)
        # Reset the flux bound 
        set_specific_bounds(model = WT, flux_bounds = {exchrxn:[0,None]}, reset_flux_bounds = False)
    # ['EX_pydx(e)', 'EX_pydxn(e)']
    print 'mutants_auxotrophy[serC] =', mutants_auxotrophy['serC']  
    """

    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('mutants_auxotrophy = {}\n')
            for gen in sorted(mutants_rxn_info_AAs.keys()):
                f.write("mutants_auxotrophy['{}'] = {}\n".format(gen, mutants_auxotrophy[gen]))
            print '\n'
            f.write('\nmutants_exitflags = {}\n')
            for gen in sorted(mutants_rxn_info_AAs.keys()):
                f.write("mutants_exitflags['{}'] = '{}'\n".format(gen, mutants_exitflags[gen]))


def find_prob_mutants_combs(leaky_traits_num = 2, results_filename = '', stdout_msgs = True, warnings = True):
    """
    Finds problematic mutant pairs
    There are cases where mutant 1 is auxotroph for AA1, mutant2 is auxotroph for AA2 but a strain having both
    mutation1 and mutations 2 is auxotroph for not only AA1 and AA2 but also for additional compounds. This 
    function finds these cases and remove them from further considerations. An example is (alaA_alaC_avtA, ilvE)
    where, alaA_alaC_avtA is auxotroph for ala_L, ilvE is auxotroph for ile_L but (alaA_alaC_avtA, ilvE) is 
    auxotorph not only for ala_L and ile_L but also for val_L

    INPUTS:
    -------
    leaky_traits_num: 
    Number of leaky traits. For example, if the value of this 
    parameter is two, we consider pairs of mutants. If three, we
    consider triples, etc.

    results_filename: 
    The base for the file name storing the results
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True to False')

    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True to False')

    if not isinstance(leaky_traits_num,int):
        raise TypeError('leaky_traits_num must be an integer')
    elif leaky_traits_num < 0:
        raise TypeError('leaky_traits_num must be a non-negative integer')

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds) 

    # Calculate death rate when not enough glucose is around
    glc_D_e = WT.compounds_by_id['glc_D_e']
    glc_D_e.ms_calc(model = WT, stdout_msgs = False)
    glc_D_e.biomass_yield_calc(model = WT, stdout_msgs = False)
    death_rate = - glc_D_e.biomass_yield*glc_D_e.ms

    # Load the list of rxns that must be off in each mutant
    from mutants_rxn_info_iJO1366 import mutants_rxn_info_AAs, genes_AA_map

    # Load the list of exchange rxns for compounds each mutant needs to survive
    from mutants_auxotrophy_AAs import mutants_auxotrophy 
    # Mutants that cannot be rescued at all
    not_rescued_mutants = [m for m in mutants_auxotrophy.keys() if len(mutants_auxotrophy[m]) == 0]

    # Mutants that are equivalent, i.e., they block the biosynthesis of the same amino acid. We consider
    # only from each group. These lists are extracted according to the contents of auxoMetabs.py 
    equivalentMutants = {}
    equivalentMutants['hisB'] = ['hisC', 'hisD', 'hisI']
    equivalentMutants['leuB'] = ['leuC', 'leuD']
    equivalentMutants['trpA'] = ['trpB']
    # Mutants that should not be considered
    not_consider_equivalent_mutatns = [m for k in equivalentMutants.keys() for m in equivalentMutants[k]]
    
    # All possible pair/triple/... combinations. The number of 
    # combinations is determined by leaky_traits_num. Here, we do not 
    # consider the ones that cannot be rescued. Also consider only one representative from
    # each group in equivalentMutants
    traits_combs = list(itertools.combinations([m for m in mutants_rxn_info_AAs.keys() if m not in not_rescued_mutants and m not in not_consider_equivalent_mutatns],r = leaky_traits_num))   
    print '\nThe total # of mutant combinations to examine = %i' % len(traits_combs)

    # Initialize the file
    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('results = {}\n')

    #---- Loop over mutants ----
    # Mutant names are actually trait names
    for trait_names in traits_combs: 
    #for trait_names in [('alaA_alaC_avtA', 'ilvE')]: 

        if stdout_msgs:
            print '{}. {}'.format(traits_combs.index(trait_names),trait_names)      

        # List of exchange reactions modeling the leaky behavior of each
        # trait (mutation), which is equivalent to the list of reactions
        # that can rescue a mutant.  
        # The following is a list of lists where the first element is the set
        # of exchange reactions needed by mutant 1 and the second element is 
        # the list fo exxchange reactions needed by mutant 2 to be rescued
        # Example: if mutants_auxotrophy[m1] = [['r1','r2'],['r3']] and 
        #             mutants_auxotrophy[m2 = [['r4','r5'],['r6']]
        # x = [[('m1', ['r1', 'r2']), ('m1', ['r3'])]
        #      [('m3', ['r7']), ('m3', ['r8', 'r9'])]
        #      [('m2', ['r4', 'r5']), ('m2', ['r6'])]]
        x = [[(m,v) for v in mutants_auxotrophy[m]] for m in trait_names]
        # and traits_leaky_exchrxns_combs = 
        #                          [{'m1':['r1', 'r2'], 'm2':['r4', 'r5']], 
        #                          ['m1':['r1', 'r2'], 'm2':['r6']], 
        #                          ['m1':['r3'], 'm2':['r4', 'r5']], 
        #                          ['m1':['r3'], 'm2':['r6']]]
        # This works for any number of combinations
        traits_leaky_exchrxns_combs = [dict(a) for a in itertools.product(*x)]
        traits_leaky_exchrxns_combs_considered = traits_leaky_exchrxns_combs

        #---- Identify all possible cooperations of leaky traits ----
        for traits_leaky_exchrxns in traits_leaky_exchrxns_combs_considered:

            #-- Create the strains serving as game players --
            # A dictionary with keys being names of strains and values
            # being the corresponding model objects
            strains = {}
 
            # The first strain is wild-type strain, which can synthesize
            # all amino acids (possesses all leaky traits)
            strains['wild_type'] = WT
            strains['wild_type'].knockedout_traits = ()
            strains['wild_type'].knockedout_rxn_ids = []

            # Find all possible combinations of trait knockouts including 
            # those with one, two, etc up to leaky_traits_num. For example,
            # if we have three traits m1, m2 and m3, then 
            # traits_ko_combs = [(m1,), (m2,), (m3,),       ---> r = 1 
            #                    (m1,m2), (m1,m3), (m3,m3), ---> r = 2
            #                    (m1, m2, m3)]              ---> r = 3
            traits_ko_combs = [a for r in range(1,leaky_traits_num+1) for a in itertools.combinations(trait_names, r = r)]

            for traits_ko in [tko for tko in traits_ko_combs if len(tko) == leaky_traits_num]:
                strain = deepcopy(WT)
                # Knockout traits (a tuple)
                strain.knockedout_traits = traits_ko
                # knockout_rxn_ids is a list of lists. We convert it to a 
                # list of strings later on
                strain.knockedout_rxn_ids = [mutants_rxn_info_AAs[trait_name] for trait_name in traits_ko]
                strain.id = '_'.join([genes_AA_map[tko] for tko in traits_ko])
                strains[strain.id] = strain

                set_specific_bounds(model = strain, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    
                set_specific_bounds(model = strain, flux_bounds = dict([(rxn_id,[0,0]) for rxn_id in list(itertools.chain(*strain.knockedout_rxn_ids))]), reset_flux_bounds = False)
                for r in strain.reactions:
                    r.objective_coefficient = 0
                strain.biomass_reaction.objective_coefficient = 1
                strain.fba(build_new_optModel = False, store_opt_fluxes = False, stdout_msgs = False)

                if strain.fba_model.solution['exit_flag'] == 'globallyOptimal' and strain.fba_model.solution['objective_value'] > 0.001:
                    raise userError('FBA for strain {} resultsed in a biomass flux value of greater than 0.001: {}'.format(traits_ko, strain.fba_model.solution['objective_value']))

                # Now check if supplying the required nutrients can rescue the mutant
                set_specific_bounds(model = strain, flux_bounds = dict([(exchrxn,[-1, None]) for ko in traits_ko for exchrxn in mutants_auxotrophy[ko][0]]), reset_flux_bounds = False)
                strain.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False)
                if strain.fba_model.solution['exit_flag'] != 'globallyOptimal':
                    with open(results_filename,'a') as f:
                        f.write("results[{}] = 'infeasible_fba'\n".format(traits_ko))
                    if stdout_msgs:
                        print '\t**infeasible fba'
                if strain.fba_model.solution['exit_flag'] == 'globallyOptimal' and strain.fba_model.solution['objective_value'] <= 0.001:
                    with open(results_filename,'a') as f:
                        f.write("results[{}] = 'zero_biomass'\n".format(traits_ko))
                    if stdout_msgs:
                        print '\t**zero bioass flux'


def master_func_BQH(start_pos = None, end_pos = None, start_pos_leak = None, end_pos_leak = None, selected_indices = None, selected_indices_leak = None, leakiness_levels = [0.5], leaky_traits_num = 2, max_game_players_num = None, simulate_rep_dynamics = True, t0 = 0, dt = 1, tf = 2000, exclude_infeas_fba = False, save_details = False, results_filename_base = '', results_filename = '', stdout_msgs = True, stdout_msgs_details = False, stdout_msgs_fba = False, warnings = True):
    """
    Simulates the game between wild-type strain synthesizing all amino acids 
    and mutant strains each lacking one ore more genes for amino acid 
    synthesis

    INPUTS:
    -------
    strat_pos & end_pos: 
    Start and end positions of the array containing all possible
    mutant pairs to examine                                 

    start_pos_leak & end_pos_leak: 
    Same as start_pos and end_pos but for leakiness_levels

    selected_indices:
    A list of selected indices corresponding to particular trait combinations, if one wants to do the simulations
    only for those selected pairs

    NOTE: The user can enter both the start and end  positions numbers as well as  
          selected_indices assuming that they start at one. The code will take care of
          array indexing convention of python (indexing starts at zero) 

    leakiness_levels: 
    A list containing the level of leakiness for each 
    strain. The elements of this must be floats between 0 and 1
    denoting the fraction of maximum possible leakiness

    leaky_traits_num: 
    Number of leaky traits. For example, if the value of this 
    parameter is two, we consider pairs of mutants. If three, we
    consider triples, etc.

    max_game_players_num:
    Max of number of games pleyers to consider when consutrcuting games. A good practice to set this 
    to max_game_players_num in order to allow the cross-feeding of all possible traits. This minimum value
    for this arguemnt is two. If No value is provided (None) then it is set to the maximum possible number, 
    which is equal to the total number of strains. For example, for a two-trait system, the possible strains
    are 11, 10, 01 and 00 so max_game_players_num is set to 4 (this is set in gameCreator_BQH function) 

    simulate_rep_dynamics:
    If True, simulates the replicator's dynamics
   
    t0, tf, dt:
    Start and final simulation times and delta for discretizing the ODE 

    exclude_infeas_fba:
    If True any community member with an infeasible FBA problem is excluded from the list of 
    community members and we redo the computations of payoffs (Default False)

    save_details:
    If True, the payoff matrix is saved into files. Note that this may increase the file size
    significantly when specially when we deal with a large number of players

    results_filename_base: 
    The base for the file name storing the results
    Example: results_filename_base = 'results/emc_results'. 
    The code will add the start and end positions to the file name.
    Example: 'results/emc_results_1_500.txt'

    results_filename:
    The same as results_filename_base except that if this parameter is provided
    the start and end positions are not added to the file name. Provide only 
    results_filename_base or results_filename 
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True to False')
    if not isinstance(stdout_msgs_details,bool):
        raise TypeError('stdout_msgs_details must be either True to False')
    if not isinstance(stdout_msgs_fba,bool):
        raise TypeError('stdout_msgs_fba must be either True to False')

    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True to False')

    if not isinstance(leaky_traits_num,int):
        raise TypeError('leaky_traits_num must be an integer')
    elif leaky_traits_num < 0:
        raise TypeError('leaky_traits_num must be a non-negative integer')

    if not isinstance(max_game_players_num,int) and max_game_players_num != None:
        raise TypeError('max_game_players_num must be either None or an integer')
    elif isinstance(max_game_players_num,int) and max_game_players_num < 2:
        raise TypeError('max_game_players_num must be an integer greater than or equal to two')

    if not isinstance(leakiness_levels,list):
        raise TypeError('A list expected for leakiness_levels. A {} was provided instead'.format(type(leakiness_level))) 
    elif len([v for v in leakiness_levels if (not isinstance(v,float) and not isinstance(v,int)) or (v < 0 or v > 1)]) > 0: 
       raise ValueError('Values of leakiness_levels must be integers or floats between zero and one')      

    if not isinstance(simulate_rep_dynamics,bool):
        raise TypeError('simulate_rep_dynamics must be either True to False')

    if not isinstance(save_details,bool):
        raise TypeError('save_details must be either True to False')

    if not isinstance(t0,int) and not isinstance(t0, float): 
        raise TypeError('t0 must be either an integer or a float')   
    if not isinstance(tf,int) and not isinstance(tf, float): 
        raise TypeError('tf must be either an integer or a float')   
    if not isinstance(dt,int) and not isinstance(dt, float): 
        raise TypeError('dt must be either an integer or a float')   

    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds) 

    # Calculate death rate when not enough glucose is around
    glc_D_e = WT.compounds_by_id['glc_D_e']
    glc_D_e.ms_calc(model = WT, stdout_msgs = False)
    glc_D_e.biomass_yield_calc(model = WT, stdout_msgs = False)
    death_rate = - glc_D_e.biomass_yield*glc_D_e.ms

    # Load the list of rxns that must be off in each mutant
    from mutants_rxn_info_iJO1366 import mutants_rxn_info_AAs, genes_AA_map

    # Load the list of exchange rxns for compounds each mutant needs to survive
    from mutants_auxotrophy_AAs import mutants_auxotrophy 
    # Mutants that cannot be rescued at all
    not_rescued_mutants = [m for m in mutants_auxotrophy.keys() if len(mutants_auxotrophy[m]) == 0]

    # Load the min and max flux of exchange reactions used for cooperation
    from coopr_exchrxns_fluxRanges_AAs import coopr_exchrxns_fluxRanges

    # Mutants that are equivalent, i.e., they block the biosynthesis of the same amino acid. We consider
    # only from each group. These lists are extracted according to the contents of auxoMetabs.py 
    equivalentMutants = {}
    equivalentMutants['hisB'] = ['hisC', 'hisD', 'hisI']
    equivalentMutants['leuB'] = ['leuC', 'leuD']
    equivalentMutants['trpA'] = ['trpB']
    # Mutants that should not be considered
    not_consider_equivalent_mutatns = [m for k in equivalentMutants.keys() for m in equivalentMutants[k]]
    
    # All possible pair/triple/... combinations. The number of 
    # combinations is determined by leaky_traits_num. Here, we do not 
    # consider the ones that cannot be rescued. Also consider only one representative from
    # each group in equivalentMutants
    traits_combs = list(itertools.combinations([m for m in mutants_rxn_info_AAs.keys() if m not in not_rescued_mutants and m not in not_consider_equivalent_mutatns],r = leaky_traits_num))   

    # Remove ('alaA_alaC_avtA', 'ilvE') from any combinations because any strain containing this double mutation is auxtorphic for
    # not only ala_L and ile_L but also for val_L
    traits_combs = [tc for tc in traits_combs if not ('alaA_alaC_avtA' in tc and 'ilvE' in tc)]
    print '\nThe total # of mutant combinations to examine = %i' % len(traits_combs)

    for tc in traits_combs:
        print traits_combs.index(tc) + 1,'.\t',tc
    print '\n'

    if start_pos == None and end_pos == None and selected_indices == None:
        traits_combs_to_consider = traits_combs
    elif start_pos != None and end_pos != None and selected_indices == None:
        traits_combs_to_consider = traits_combs[start_pos - 1:end_pos]
    elif start_pos == None and end_pos == None and selected_indices != None:
        traits_combs_to_consider = [traits_combs[i - 1] for i in selected_indices]
    else:
        raise ValueError('start_pos, end_pos and selected_indices cannot be provided at the same time')

    # Name of the output file storing the results
    if results_filename_base != '' and results_filename == '':
        results_filename = results_filename_base + '_' + str(leaky_traits_num) + 'LT_' + str(start_pos) + '_' + str(end_pos) + '.py'

    # Initialize the file
    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('results = {}\n')

    print 'Simulating slice {}\n'.format((start_pos,end_pos))

    BQH_games = []

    #---- Loop over mutants ----
    # Mutant names are actually trait names
    #for trait_names in [('ilvE',)]:
    #for trait_names in [('lysA','ilvE')]:
    for trait_names in traits_combs_to_consider: 
      
        print '\n**** {}. {} ****\n'.format(traits_combs.index(trait_names) + 1,trait_names)

        # List of exchange reactions modeling the leaky behavior of each
        # trait (mutation), which is equivalent to the list of reactions
        # that can rescue a mutant.  
        # The following is a list of lists where the first element is the set
        # of exchange reactions needed by mutant 1 and the second element is 
        # the list fo exxchange reactions needed by mutant 2 to be rescued
        # Example: if mutants_auxotrophy[m1] = [['r1','r2'],['r3']] and 
        #             mutants_auxotrophy[m2 = [['r4','r5'],['r6']]
        # x = [[('m1', ['r1', 'r2']), ('m1', ['r3'])]
        #      [('m3', ['r7']), ('m3', ['r8', 'r9'])]
        #      [('m2', ['r4', 'r5']), ('m2', ['r6'])]]
        x = [[(m,v) for v in mutants_auxotrophy[m]] for m in trait_names]
        # and traits_leaky_exchrxns_combs = 
        #                          [{'m1':['r1', 'r2'], 'm2':['r4', 'r5']], 
        #                          ['m1':['r1', 'r2'], 'm2':['r6']], 
        #                          ['m1':['r3'], 'm2':['r4', 'r5']], 
        #                          ['m1':['r3'], 'm2':['r6']]]
        # This works for any number of combinations
        traits_leaky_exchrxns_combs = [dict(a) for a in itertools.product(*x)]
        print 'Total possible combinations of leaky_exchrxn_sets (resucing reaction sets): {} '.format(traits_leaky_exchrxns_combs) 
        print 'The total possible combinations of leaky_exchrxn_sets (resucing reaction sets): {} '.format(len(traits_leaky_exchrxns_combs)) 
        print 'The total possible combinations of leakiness levels: {}'.format(len(list(itertools.product(*[leakiness_levels for trait_name in trait_names])))) 
        print 'The total # of cases to consider for {} considering all possible combinations of resucing reactions and leakiness levels = {}\n'.format(trait_names,len(traits_leaky_exchrxns_combs)*len(list(itertools.product(*[leakiness_levels for trait_name in trait_names]))))

        #-- Consider only combinations containing an amino acid exchange reaction to reduce runtime --
        # If having exchange reactions for all compounds in the file mutants_auxotrophy_AA.py
        if False:
            traits_leaky_exchrxns_combs_AAOnly = []
            # First find out the combination with the most number of AA exchange rxns
            leaky_comb_exchAA_num = {}
            for d in traits_leaky_exchrxns_combs:
                for k in d.keys():
                   d[k] = tuple(d[k])
                leaky_comb_exchAA_num[tuple(d.items())] = len(set([exchrxn for exchrxn_list in d.values() for exchrxn in exchrxn_list if '_L(e)' in exchrxn]))
            sorted_leaky_comb_exchAA_num = sorted(leaky_comb_exchAA_num.items(), key = lambda x: x[1], reverse = True)
            traits_leaky_exchrxns_combs_considered = [dict(list(sorted_leaky_comb_exchAA_num[0][0]))] 
        # If having only realted AAs in the file mutants_auxotrophy_AA.py
        else:
            traits_leaky_exchrxns_combs_considered = traits_leaky_exchrxns_combs
        print 'The leakiness level combination considered: {}'.format(traits_leaky_exchrxns_combs_considered)
        print 'The total # of cases to consider for {} considering all possible combinations of resucing reactions and leakiness levels = {}\n'.format(trait_names,len(traits_leaky_exchrxns_combs)*len(traits_leaky_exchrxns_combs_considered))

        #---- Identify all possible cooperations of leaky traits ----
        for traits_leaky_exchrxns in traits_leaky_exchrxns_combs_considered:

            print '\ntraits_leaky_exchrxns = ', traits_leaky_exchrxns
            print '\n--- cooperation strategies = {} ---\n'.format(traits_leaky_exchrxns)
            # Compute ms and biomass yield for compounds needed to rescue 
            # the mutant (or those in traits_leaky_exchrxns)
            set_specific_bounds(model = WT, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    
            for exch_rxn_id in [r for trait_name in trait_names for r in traits_leaky_exchrxns[trait_name]]: 
                exch_rxn = WT.reactions_by_id[exch_rxn_id]
                # The compound participating in the exchange reaction
                cpd = exch_rxn.reactants[0] 
                cpd.ms_calc(model = WT)    
                cpd.biomass_yield_calc(model = WT)    
    
            #-- Create the strains serving as game players --
            # A dictionary with keys being names of strains and values
            # being the corresponding model objects
            strains = {}
 
            # The first strain is wild-type strain, which can synthesize
            # all amino acids (possesses all leaky traits)
            strains['wild_type'] = WT
            strains['wild_type'].knockedout_traits = ()
            strains['wild_type'].knockedout_rxn_ids = []

            # Find all possible combinations of trait knockouts including 
            # those with one, two, etc up to leaky_traits_num. For example,
            # if we have three traits m1, m2 and m3, then 
            # traits_ko_combs = [(m1,), (m2,), (m3,),       ---> r = 1 
            #                    (m1,m2), (m1,m3), (m3,m3), ---> r = 2
            #                    (m1, m2, m3)]              ---> r = 3
            traits_ko_combs = [a for r in range(1,leaky_traits_num+1) for a in itertools.combinations(trait_names, r = r)]
            #print '\ntraits_ko_combs = ',traits_ko_combs
            for traits_ko in traits_ko_combs:
                strain = deepcopy(WT)
                # Knockout traits (a tuple)
                strain.knockedout_traits = traits_ko
                # knockout_rxn_ids is a list of lists. We convert it to a 
                # list of strings later on
                strain.knockedout_rxn_ids = [mutants_rxn_info_AAs[trait_name] for trait_name in traits_ko]
                strain.id = '_'.join([genes_AA_map[tko] for tko in traits_ko])
                strains[strain.id] = strain

                print '\n-- FBA for ' + str(traits_ko) + ' Ecoli mutant --'
                if stdout_msgs:
                    print 'knockedout rxn ids = ',strain.knockedout_rxn_ids
                set_specific_bounds(model = strain, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    
                set_specific_bounds(model = strain, flux_bounds = dict([(rxn_id,[0,0]) for rxn_id in list(itertools.chain(*strain.knockedout_rxn_ids))]), reset_flux_bounds = False)
                for r in strain.reactions:
                    r.objective_coefficient = 0
                strain.biomass_reaction.objective_coefficient = 1
                strain.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = True)
                if strain.fba_model.solution['exit_flag'] == 'globallyOptimal' and strain.fba_model.solution['objective_value'] > 0.001:
                    raise userError('FBA for strain {} resultsed in a biomass flux value of greater than 0.001: {}'.format(traits_ko, strain.fba_model.solution['objective_value']))
                # Now check if supplying the required nutrients can rescue the mutant
                print '\n-- FBA for checking the rescue of ' + str(traits_ko) + ' Ecoli mutant --'
                if stdout_msgs_details:
                   print 'rescuing rxns for this mutant = ', [exchrxn for ko in traits_ko for exchrxn in mutants_auxotrophy[ko][0]]
                set_specific_bounds(model = strain, flux_bounds = dict([(exchrxn,[-1, None]) for ko in traits_ko for exchrxn in mutants_auxotrophy[ko][0]]), reset_flux_bounds = False)
                strain.fba(build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = True)
                if strain.fba_model.solution['exit_flag'] != 'globallyOptimal':
                    raise userError('FBA problem to test whether the mutant {} can be rescued by the provided exchange reactions in mutants_auxotrophy was not solved to optimality'.format(traits_ko))
                if strain.fba_model.solution['exit_flag'] == 'globallyOptimal' and strain.fba_model.solution['objective_value'] <= 0.001:
                    raise userError('Strain {} cannot be rescued by the supply of nutrients specified by mutants_auxotrophy as the biomass flux from using FBA is less than 0.001: {}'.format(traits_ko, strain.fba_model.solution['objective_value']))

            # All possible combinations for the leakiness levels for all
            # leaky traits. Example: If traits_leaky_exchrxns_combs = 
            #       {'m1':['r1','r2'], 'm2':['r3']}
            # Then x = [ [(('r1','r2'), 0), (('r1','r2'), 0.1), (('r1','r2'),
            #              0.2), (('r1','r2'), 0.3)],
            #            [('r3', 0), ('r3', 0.1), ('r3', 0.2), ('r3', 0.3)] ]
            x = [[(tuple(traits_leaky_exchrxns[trait_name]),v) for v in leakiness_levels] for trait_name in traits_leaky_exchrxns.keys()]
            # and leakiness_levels_combs = [{('r1','r2'): 0, ('r2'): 0},
            #                               {('r1','r2'): 0, ('rr'): 0.1},...]
            leakiness_levels_combs = [dict(a) for a in itertools.product(*x)] 

            if start_pos_leak == None and end_pos_leak == None and selected_indices_leak == None:
                leakiness_levels_combs_considered = leakiness_levels_combs
            elif start_pos_leak != None and end_pos_leak != None:
                leakiness_levels_combs_considered = leakiness_levels_combs[start_pos_leak - 1:end_pos_leak]
            elif start_pos_leak == None and end_pos_leak == None and selected_indices_leak != None:
                leakiness_levels_combs_considered = [leakiness_levels_combs[i-1] for i in selected_indices]
            else:
                raise ValueError('start_pos_leak, end_pos_leak and selected_indices_leak cannot be provided at the same time')

            if stdout_msgs_details:
                print '\n***leakiness_levels_combs considered: {}'.format(leakiness_levels_combs_considered)
            print '\nThe total # of leakiness_levels_combs considered: {}'.format(len(leakiness_levels_combs_considered))

            counter = 0

            for leak_level in leakiness_levels_combs_considered:

                counter += 1

                # The following is to find out the required runtime for one case
                start_time_wt = time.time()
          
                print '{}. leakiness levels = {}'.format(counter,leak_level)

                # Assign the level of leakiness for mutants
                for strain_name in strains.keys():
                    # A dicitonay with keys being tuples containing leaky
                    # rxn sets and values being their leakiness levels
                    strains[strain_name].leaky_traits = {}
                    for leaky_trait in [lt for lt in leak_level.keys() if list(lt) not in [rlist for trait_name in strains[strain_name].knockedout_traits for rlist in mutants_auxotrophy[trait_name]]]:
                        strains[strain_name].leaky_traits[leaky_trait] = leak_level[leaky_trait] 
                    if stdout_msgs_details:
                        print 'Strain = {}, knocedout_traits = {} , Kockedout_rxn_ids = {} , leaky_traits = {}\n '.format(strain_name, strains[strain_name].knockedout_traits, strains[strain_name].knockedout_rxn_ids, strains[strain_name].leaky_traits)

                #-- First check if wild-type can produce all leaky traits simultaneously --
                # with this given leakiness combinations growth condition
                if stdout_msgs:
                    print '--- Check whether the wild-type can satisfy the imposed leakienes levels --'
                set_specific_bounds(model = strains['wild_type'], flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    

                # Leaky traits. If leaky_traits = {('r1','r2'):0.1, 
                # ('r3','r4'):0.2, ('r5',):0.3}, then flux_bounds = 
                # [('r1', [0.1*max_coopr_level, None]), 
                # ('r2', [0.1*max_coopr_level,None]), 
                # ('r3', [0.2*max_coopr_level, None]), 
                # ('r4', [0.2*max_coopr_level,None]), 
                # ('r5', [0.3*max_coopr_level, None])] 
                set_specific_bounds(model = strains['wild_type'], flux_bounds = dict([(exchrxn_id,[strains['wild_type'].leaky_traits[rxn_set]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for rxn_set in strains['wild_type'].leaky_traits.keys() for exchrxn_id in rxn_set]), reset_flux_bounds = False)    
                #print 'leakiness levels = ', dict([(exchrxn_id,[strains['wild_type'].leaky_traits[rxn_set]*coopr_exchrxns_fluxRanges[exchrxn_id][1], None]) for rxn_set in strains['wild_type'].leaky_traits.keys() for exchrxn_id in rxn_set])
                strains['wild_type'].fba(stdout_msgs = stdout_msgs)
                if strains['wild_type'].fba_model.solution['exit_flag'] == 'globallyOptimal':
                    WT_can_produce_all_leakyTraits = True 
                else:
                    WT_can_produce_all_leakyTraits = False 
                    if warnings:
                        print '\t**WARNING! This leakiness combinations makes the wild-type FBA problem is infeasible'

                # Creating a shared memory using the manager
                input_data = {} 
                input_data['strains'] = strains
                input_data['growthMedium_flux_bounds'] = growthMedium_flux_bounds
                input_data['death_rate'] = death_rate
                input_data['coopr_exchrxns_fluxRanges'] = coopr_exchrxns_fluxRanges
                input_data['max_game_players_num'] = max_game_players_num
                input_data['WT_can_produce_all_leakyTraits'] = WT_can_produce_all_leakyTraits
                input_data['exclude_infeas_fba'] = exclude_infeas_fba
                input_data['simulate_rep_dynamics'] = simulate_rep_dynamics
                input_data['t0'] = t0
                input_data['tf'] = tf
                input_data['dt'] = dt
                input_data['save_details'] = save_details
                input_data['results_filename'] = results_filename
                input_data['warnings'] = warnings
                input_data['stdout_msgs'] = stdout_msgs
                input_data['stdout_msgs_details'] = stdout_msgs_details
                input_data['stdout_msgs_fba'] = stdout_msgs_fba

                output_data = Manager().dict()
                output_data['BQH_game'] = None
        
                p = Process(target = gameCreator_BQH, args = (input_data, output_data))
                p.start()
                p.join() 
                if p.exitcode > 0:
                    raise userError('Error in python subprocess. Please check performDMMM\n')
                else:
                    BQH_games.append(output_data['BQH_game'])

                # The following is to find out the required runtime for one case
                elapsed_time_wt = str(timedelta(seconds = time.time() - start_time_wt))
                if stdout_msgs:
                    print '\nRequired time for one case is: ',elapsed_time_wt 

    return BQH_games

def create_jobs_files_BQH(leaky_traits_num, interval_size, jobs_to_create = 'game_jobs', leakiness_levels = None, leakiness_levels_str = None):
    """
    Creates the job files

    leakiness_levels_str:
    str of the command generating leakiness_levels. For example, 'range(0,2,100)'.
    This will avoid writing long lists in the python command
    """
    from tools.utilities.create_job_files import create_intervals, create_job_files

    # Load the list of rxns that must be off in each mutant
    from mutants_rxn_info_iJO1366 import mutants_rxn_info_AAs, genes_AA_map

    if jobs_to_create == 'game_jobs':
        # Load the list of exchange rxns for compounds each mutant needs to survive
        from mutants_auxotrophy_AAs import mutants_auxotrophy 
        # Mutants that cannot be rescued at all
        not_rescued_mutants = [m for m in mutants_auxotrophy.keys() if len(mutants_auxotrophy[m]) == 0]
    
        # Load the min and max flux of exchange reactions used for cooperation
        from coopr_exchrxns_fluxRanges_AAs import coopr_exchrxns_fluxRanges
    
        # Mutants that are equivalent, i.e., they block the biosynthesis of the same amino acid. We consider
        # only from each group. These lists are extracted according to the contents of auxoMetabs.py 
        equivalentMutants = {}
        equivalentMutants['hisB'] = ['hisC', 'hisD', 'hisI']
        equivalentMutants['leuB'] = ['leuC', 'leuD']
        equivalentMutants['trpA'] = ['trpB']
        # Mutants that should not be considered
        not_consider_equivalent_mutatns = [m for k in equivalentMutants.keys() for m in equivalentMutants[k]]
        
        # All possible pair/triple/... combinations. The number of 
        # combinations is determined by leaky_traits_num. Here, we do not 
        # consider the ones that cannot be rescued. Also consider only one representative from
        # each group in equivalentMutants
        traits_combs = list(itertools.combinations([m for m in mutants_rxn_info_AAs.keys() if m not in not_rescued_mutants and m not in not_consider_equivalent_mutatns],r = leaky_traits_num))   
        # Remove ('alaA_alaC_avtA', 'ilvE') from any combinations because any strain containing this double mutation is auxtorphic for
        # not only ala_L and ile_L but also for val_L
        traits_combs = [tc for tc in traits_combs if not ('alaA_alaC_avtA' in tc and 'alaA_alaC_avtA' in tc)]
        print '\nThe total # of mutant combinations to examine = %i' % len(traits_combs)
    
        # The total possible combinations of leakiness levels
        leak_level_combs = list(itertools.product(*[leakiness_levels for trait_name in traits_combs[0]]))
    
        print '\nTotal # of mutant combinations is ',len(traits_combs)
        print 'Total # of leakiness level combinations is ',len(leak_level_combs)
        print 'Total # of cases: ',len(leak_level_combs)*len(traits_combs)
     
        print 'Total # of cases in each interval: ',len(traits_combs)*interval_size
    
        #create_intervals(total_cases_num = len(traits_combs), interval_size = interval_size, one_case_req_time = 12*len(leak_level_combs), stdout_msgs = True)

        if leaky_traits_num == 2:    
            max_game_players_num = 3
            results_filename_base = 'results/leakiness_games_BQH_2LT'
            job_filename_base = 'jobs/job_leakiness_games_BQH_2LT'
        
            # NOTE: Do NOT specify directory jobs/ for the output files as any path to the output file uses 
            # the directory where the jobs itself is located as the source not what is specified by 'cd [directory]'.
            # Therefore, if one species 'jobs/output_filename' the codes will have the status of Eqw (error in jobs) 
            # because it looks for a directory 'jobs' in jobs directory.
            joboutput_filename_base = 'job_leakiness_games_BQH_2LT'
        
            commands = ['python -c """from __future__ import division;from leakiness_ss import master_func_BQH;master_func_BQH(start_pos = None, end_pos = None, leaky_traits_num = ' + str(leaky_traits_num) + ', max_game_players_num = ' + str(max_game_players_num) + ', leakiness_levels = ' + leakiness_levels_str + ', simulate_rep_dynamics = True, dt = 1, tf = 5000, save_details = False, results_filename = None, stdout_msgs = False, warnings = True)"""']
        
            create_job_files(total_cases_num = len(traits_combs), interval_size = interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos = None', end_pos_input_format = 'end_pos = None', main_code_dir = '/usr2/postdoc/alizom/work/EcoliPairs', commands = commands, max_walltime = 144)

        if leaky_traits_num == 3:    
            # ('lysA', 'ilvE_avtA', 'ilvE') --> 862
            case_to_consider = 'lys_ile_val'
            start_pos, end_pos = 862, 862

            # ('lysA', 'thrC', 'ilvE') --> 902
            #case_to_consider = 'lys_ile_thr'
            #start_pos, end_pos = 902, 902

            # ('lysA', 'ilvE', 'trpA') --> 916
            #case_to_consider = 'lys_ile_trp'
            #start_pos, end_pos = 916, 916

            results_filename_base = 'results/leakiness_games_BQH_3LT_' + case_to_consider
            job_filename_base = 'jobs/job_leakiness_games_BQH_3LT_' + case_to_consider
        
            max_game_players_num = 4

            # NOTE: Do NOT specify directory jobs/ for the output files as any path to the output file uses 
            # the directory where the jobs itself is located as the source not what is specified by 'cd [directory]'.
            # Therefore, if one species 'jobs/output_filename' the codes will have the status of Eqw (error in jobs) 
            # because it looks for a directory 'jobs' in jobs directory.
            joboutput_filename_base = 'job_leakiness_games_BQH_3LT_' +  case_to_consider
        
            commands = ['python -c """from __future__ import division;from leakiness_ss import master_func_BQH;master_func_BQH(start_pos = ' + str(start_pos) + ', end_pos = ' + str(end_pos) + ', start_pos_leak = None, end_pos_leak = None, leaky_traits_num = ' + str(leaky_traits_num) + ', max_game_players_num = ' + str(max_game_players_num) + ', leakiness_levels = ' + leakiness_levels_str + ', simulate_rep_dynamics = True, dt = 1, tf = 5000, save_details = False, results_filename = None, stdout_msgs = False, warnings = True)"""']

            create_job_files(total_cases_num = len(leak_level_combs), interval_size = interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos_leak = None', end_pos_input_format = 'end_pos_leak = None', main_code_dir = '/usr2/postdoc/alizom/work/EcoliPairs', commands = commands, max_walltime = 168)


    elif jobs_to_create == 'repdyn_jobs':

        game_results = load_data_from_python_file(file_name = 'results/leakiness_games_BQH_2LT_all_old.py', var_names = ['results'])['results']

        print 'Total # of cases  = ',len(game_results.keys())

        results_filename_base = 'results/leakiness_games_BQH_2LT'
        job_filename_base = 'jobs/job_leakiness_BQH_repdyn_2LT'
    
        # NOTE: Do NOT specify directory jobs/ for the output files as any path to the output file uses 
        # the directory where the jobs itself is located as the source not what is specified by 'cd [directory]'.
        # Therefore, if one species 'jobs/output_filename' the codes will have the status of Eqw (error in jobs) 
        # because it looks for a directory 'jobs' in jobs directory.
        joboutput_filename_base = 'job_leakiness_BQH_repdyn_2LT'

        commands = ['python -c """from leakiness_ss import simulate_rep_dynamics;simulate_rep_dynamics(game_results_filename = ' +  "'results/leakiness_games_BQH_2LT_all_old.py'" + ', start_pos = None, end_pos = None, tf = 5000, dt = 1, save_details = False, results_filename = None, stdout_msgs = False)"""']
    
        create_job_files(total_cases_num = len(game_results.keys()), interval_size = interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos = None', end_pos_input_format = 'end_pos = None', main_code_dir = '/usr2/postdoc/alizom/work/EcoliPairs', commands = commands, max_walltime = 72)

    else:
        raise ValueError('Unkonwn jobs_to_create')

def run_create_jobs_files_BQH():
    """   
    Runs create_jobs_files_BQH
    """   
    #--------- leaky_traits_num = 2 ----
    # Create jobs for pairs for games
    #create_jobs_files_BQH(leaky_traits_num = 2, leakiness_levels = [x/100 for x in range(0,105,5)], leakiness_levels_str = '[x/100 for x in range(0,105,5)]', interval_size = 10)

    # Create jobs for pairs for replicator dynamics
    #create_jobs_files_BQH(leaky_traits_num = 2, interval_size = 8000, jobs_to_create = 'repdyn_jobs', leakiness_levels = None, leakiness_levels_str = None)

    #--------- leaky_traits_num = 3 ----
    create_jobs_files_BQH(leaky_traits_num = 3, leakiness_levels = [0,0.05, 0.1, 0.15, 0.2] + [i/100 for i in range(30,101,10)], leakiness_levels_str = '[0,0.05, 0.1, 0.15, 0.2] + [i/100 for i in range(30,101,10)]', interval_size = 100)


def plot_gameResults_1LT(results_filename, number_of_gamePlayers, x, y, title = '', xaxis_label = '', yaxis_label = '', mixed_NashEq_label = 'Mixed', set_minor_xticks = True, invert_xaxis = True, invert_yaxis = False, output_filename_base = '', output_filetype = 'pdf'):
    """
    Plots a heatmap of the anaylsis results.
 
    INPUTS:
    -------
    results_filename: 
    A string containing the name of the file containing the results
    (use integrate_results_files function to integrate the results
    in different files

    number_of_gamePlayers:
    Number of game players for which the results should be plotted

    x & y: 
    A list containing the elements of horizontal and vertical axes 

    xaxis_label & yaxis_label: 
    x-axis and y-axis labels

    mixed_NashEq_label:
    The label of mixed Nash equilibria in the color bar

    output_filename_base: 
    Base name of the output file containing the plot. This funciton then creates 
    two files output_filename_base + '_games.pdf' and output_filename_base + '_freq.pdf'
    The first stores the plot for the games and the second for species frequencies. 

    output_filetype:
    File type of the output file (e.g., pdf, png, etc). Do not include dot "." here
    """
    from tools.utilities.load_data_fromFile import load_data_from_python_file 
    from matplotlib import colors
    from tools.utilities.plot import plot, axis, color_bar

    if not isinstance(set_minor_xticks,bool):
        raise TypeError('set_minor_xticks must be boolean')

    results = load_data_from_python_file(file_name = results_filename, var_names = ['results'])['results']

    # List of amino acids
    AAs_list = ['ala_L', 'arg_L', 'asn_L', 'asp_L', 'cys_L', 'gln_L', 'glu_L', 'gly', 'ile_L', 'leu_L', 'lys_L', 'met_L', 'phe_L', 'phe_L', 'pro_L', 'ser_L', 'thr_L', 'trp_L', 'tyr_L', 'val_L'] 
   
    #--------------- GAMES ---------------------------
    #--- Find out how many Nash eq we have in the results ---
    # (mutant, mutant)      --> Prisonder's Dilemma
    # (mutant,wild_type)    --> Snowdirft
    # (wild_type,wild_type) --> Mutually beneficial
    # Since mutant_name varies from one case to the other (depending on the AA), we'll replace
    # all actual mutant names with 'mutant' in order to identify Nash eq easier
    Nash_equilibria = []

    # The following is a dictionary where keys are the same results keys and values are Nash 
    # equilibria with above name replacements 
    results_neq = dict([(k,None) for k in results.keys()])

    for res_key in results.keys():

        AAexchrxns = [s for k in dict(res_key).keys() for s in k]

        # Extract the AA name
        AA_names = [re.sub('EX_|\(e\)','',exchrxn) for exchrxn in AAexchrxns] 
        AA_name = AA_names[0]

        # Nash eq
        Nash_eq = results[res_key]['games_info'][number_of_gamePlayers]['NashEq'] 

        # Add only the ones where AA does not appear there
        Nash_eq_AAreplaced = [n for n in Nash_eq if AA_name not in n]

        # Replace AA name with 'mutant'
        for neq in [n for n in Nash_eq if AA_name in n]:
            neq = list(neq)
            # If the Nash eq is ('lis_L','lys_L') then AA appears more than once there,
            # so repeat replacing AAname with 'mutant' as many times as AAname has been repeated
            for k in range(neq.count(AA_name)):
                neq[neq.index(AA_name)] = 'mutant'
            Nash_eq_AAreplaced.append(tuple(sorted(neq))) 

        Nash_eq_AAreplaced = sorted(Nash_eq_AAreplaced)

        if len(Nash_eq_AAreplaced) == 1:
            results_neq[res_key] = Nash_eq_AAreplaced[0]
            Nash_equilibria.append(Nash_eq_AAreplaced[0])
        else:
            results_neq[res_key] = tuple(Nash_eq_AAreplaced)
            Nash_equilibria.append(tuple(Nash_eq_AAreplaced))
            print '**WARNING! More than one Nash eq for {}: {}'.format(res_key, Nash_eq_AAreplaced)

    Nash_equilibria = list(set(Nash_equilibria))
    Nash_eq = Nash_equilibria  # just ot use a shorter name 

    print '\nNash equilibria of the systems are:'
    for neq in Nash_equilibria:
        print neq

    # 0 = PD --> red , 1 = MB --> green , 2 = SD --> cyan , 3 = Mixed --> Blue
    # Source: http://stackoverflow.com/questions/30893483/make-a-heatmap-with-a-specified-discrete-color-mapping-with-matplotlib-in-python
    pd = ('mutant', 'mutant')
    mb = ('wild_type', 'wild_type')
    sd = ('mutant', 'wild_type')
    
    # Number of known Nash equilibria (PD, MB and SD)
    if pd not in Nash_eq and mb not in Nash_eq and sd not in Nash_eq:
        known_NE_num = 0
    elif (pd in Nash_eq and mb not in Nash_eq and sd not in Nash_eq) or (pd not in Nash_eq and mb in Nash_eq and sd not in Nash_eq) or (pd not in Nash_eq and mb not in Nash_eq and sd in Nash_eq):
        known_NE_num = 1
    elif (pd in Nash_eq and mb in Nash_eq and sd not in Nash_eq) or (pd in Nash_eq and mb not in Nash_eq and sd in Nash_eq) or (pd not in Nash_eq and mb in Nash_eq and sd in Nash_eq):
        known_NE_num = 2
    elif pd in Nash_eq and mb in Nash_eq and sd in Nash_eq:
        known_NE_num = 3

    if known_NE_num == 0:
        colormap = colors.ListedColormap(['blue'])
        colorbar_ticklabels = [mixed_NashEq_label]
        data_value_mixed = 0
    elif known_NE_num == 1:
        if pd in Nash_eq:
             if len(Nash_eq) == 1:
                colormap = colors.ListedColormap(['red'])
                colorbar_ticklabels = ["Prisoner's \nDilemma"]
                data_value_pd = 0
             else:
                colormap = colors.ListedColormap(['red','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma","Mixed"]
                data_value_pd = 0
                data_value_mixed = 1
        elif mb in Nash_eq:
             if len(Nash_eq) == 1:
                colormap = colors.ListedColormap(['green'])
                colorbar_ticklabels = ["Mutually \nbeneficial"]
                data_value_mb = 0
             else:
                colormap = colors.ListedColormap(['green','blue'])
                colorbar_ticklabels = ["Mutually \nbeneficial","Mixed"]
                data_value_mb = 0
                data_value_mb = 1
        elif sd in Nash_eq:
             if len(Nash_eq) == 1:
                colormap = colors.ListedColormap(['cyan'])
                colorbar_ticklabels = ["Snowdrift"]
                data_value_sd = 0
             else:
                colormap = colors.ListedColormap(['cyan','blue'])
                colorbar_ticklabels = ["Snowdrift","Mixed"]
                data_value_sd = 0
                data_value_mixed = 1
        else:
            raise userError('Unknown case for known_NE_num = 1')

    elif known_NE_num == 2:
        if pd in Nash_eq and mb in Nash_eq:
             if len(Nash_eq) == 2:
                colormap = colors.ListedColormap(['red', 'green'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial']
                data_value_pd = 0
                data_value_mb = 1
             else:
                colormap = colors.ListedColormap(['red', 'green','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial',mixed_NashEq_label]
                data_value_pd = 0
                data_value_mb = 1
                data_value_mixed = 2
        elif pd in Nash_eq and sd in Nash_eq:
             if len(Nash_eq) == 2:
                colormap = colors.ListedColormap(['red', 'cyan'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Snowdrift']
                data_value_pd = 0
                data_value_sd = 1
             else:
                colormap = colors.ListedColormap(['red', 'cyan','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Snowdrift',mixed_NashEq_label]
                data_value_pd = 0
                data_value_sd = 1
                data_value_mixed = 2
        elif mb in Nash_eq and sd in Nash_eq:
             if len(Nash_eq) == 2:
                colormap = colors.ListedColormap(['green', 'cyan'])
                colorbar_ticklabels = ['Mutually \nBeneficial','Snowdrift']
                data_value_mb = 0
                data_value_sd = 1
             else:
                colormap = colors.ListedColormap(['green', 'cyan','blue'])
                colorbar_ticklabels = ['Mutually \nBeneficial','Snowdrift',mixed_NashEq_label]
                data_value_mb = 0
                data_value_sd = 1
                data_value_mixed = 2
        else:
            raise userError('Unknown case for known_NE_num = 2')
    elif known_NE_num == 3:
        if pd in Nash_eq and mb in Nash_eq and sd in Nash_eq:
            if len(Nash_eq) == 3:
                colormap = colors.ListedColormap(['red', 'green','cyan'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial','Snowdrift']
                data_value_pd = 0
                data_value_mb = 1
                data_value_sd = 2
            else:
                colormap = colors.ListedColormap(['red', 'green','cyan','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial','Snowdrift',mixed_NashEq_label]
                data_value_pd = 0
                data_value_mb = 1
                data_value_sd = 2
                data_value_mixed = 3
        else:
            raise userError('Unknown case for known_NE_num = 3')

    # Games
    dataG = np.zeros((len(y),len(x)))

    # Here, x is the amino acid (or cost) and y is capture efficiency
    for i, yy in enumerate(y):
        for j, xx in enumerate(x):
            # The dictionary key in results (this works only when we have one leaky trait) 
            res_key = ((('EX_' + xx + '(e)',), yy),) 

            if results_neq[res_key] == pd:
                dataG[i,j] = data_value_pd
            elif results_neq[res_key] == mb:
                dataG[i,j] = data_value_mb
            elif results_neq[res_key] == sd:
                dataG[i,j] = data_value_sd
            else:
                dataG[i,j] = data_value_mixed

    print 'min(data) = {}  , max(data) = {}'.format(dataG.min(), dataG.max())

    #-- Create a color map --
    yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) else '' for yy in y]

    # Map the short names of AAs to their full name
    #AA_fullname_map = {'ala_L':'L-Alanine', 'val_L':'L-valine', 'gly':'Glycine', 'asp_L':'L-Aspartate', 'ser_L':'L-Serine', 'glu_L':'L-Glutamate', 'asn_L':'L-Asparagine', 'gln_L':'L-Glutamine', 'leu_L':'L-Leucine', 'pro_L':'L-Proline', 'thr_L':'L-Threonine', 'ile_L':'L-Isoleucine', 'lys_L':'L-Lysine', 'cys_L':'L-Cysteine', 'phe_L':'L-Phenylalanine', 'arg_L':'L-Arginine', 'tyr_L':'Tyrosine', 'his_L':'L-Histidine', 'met_L':'L-Methionine', 'trp_L':'L-Tryptophan'} 
    AA_fullname_map = {'ala_L':'ala', 'val_L':'val', 'gly':'gly', 'asp_L':'asp', 'ser_L':'ser', 'glu_L':'glu', 'asn_L':'asn', 'gln_L':'gln', 'leu_L':'leu', 'pro_L':'pro', 'thr_L':'thr', 'ile_L':'ile', 'lys_L':'lys', 'cys_L':'cys', 'phe_L':'phe', 'arg_L':'arg', 'tyr_L':'tyr', 'his_L':'his', 'met_L':'met', 'trp_L':'trp'} 

    output_filename = output_filename_base + '_NEs.' + output_filetype
    
    x_fullname = []
    for xx in x:
        x_fullname.append(AA_fullname_map[xx])

    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':90}, plot_gridlines = True, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = invert_yaxis), fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = x_fullname, y = 100*np.array(y), data = dataG, plot_func = 'matshow', clrbar = color_bar(colormap = colormap, colorlimits = (dataG.min(),dataG.max()), label = '',  custom_ticklabels = colorbar_ticklabels, ticklabels_format = {'ticks_position':'middle'}))

    print 'The figure was saved into {}\n'.format(output_filename_base + '_games.' + output_filetype)

    #--------------- Evolutionary dynamics ---------------------------
    # Species frequences when wild-type invades (dataX_WT) or mutants invade (dataX_MT)
    dataX_WT = np.zeros((len(y),len(x)))
    dataX_MT = np.zeros((len(y),len(x)))

    # Here, x is the amino acid (or cost) and y is capture efficiency
    for i, yy in enumerate(y):
        for j, xx in enumerate(x):
            # The dictionary key in results (this works only when we have one leaky trait) 
            res_key = ((('EX_' + xx + '(e)',), yy),) 

            # Wild-type invades
            if ((xx, 0.99), ('wild_type', 0.01)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_WT[i,j] = results[res_key]['replicator_dynamics_info'][((xx, 0.99), ('wild_type', 0.01))]['wild_type']
            elif (('wild_type', 0.01), (xx, 0.99)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_WT[i,j] = results[res_key]['replicator_dynamics_info'][(('wild_type', 0.01), (xx, 0.99))]['wild_type']
            else:
                raise userError("Neither {} nor {} is in results[res_key]['replicator_dynamics_info'].keys()".format(((xx, 0.01), ('wild_type', 0.99)), (('wild_type', 0.99), (xx, 0.01))))

            # Mutant invades
            if ((xx, 0.01), ('wild_type', 0.99)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_MT[i,j] = results[res_key]['replicator_dynamics_info'][((xx, 0.01), ('wild_type', 0.99))]['wild_type']
            elif (('wild_type', 0.99), (xx, 0.01)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_MT[i,j] = results[res_key]['replicator_dynamics_info'][(('wild_type', 0.99), (xx, 0.01))]['wild_type']
            else:
                raise userError("Neither {} nor {} is in results[res_key]['replicator_dynamics_info'].keys()".format(((xx, 0.01), ('wild_type', 0.99)), (('wild_type', 0.99), (xx, 0.01))))

    print 'min(dataX_WT) = {}  , max(dataX_WT) = {}'.format(dataX_WT.min(), dataX_WT.max())
    print 'min(dataX_MT) = {}  , max(dataX_MT) = {}'.format(dataX_MT.min(), dataX_MT.max())

    custom_ticks = np.arange(0,1+0.2,0.2)

    output_filename = output_filename_base + '_freq_WTInvades.' + output_filetype
    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':90}, plot_gridlines = True, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = invert_yaxis), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = x_fullname, y = 100*np.array(y), data = dataX_WT, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)


    output_filename = output_filename_base + '_freq_MTInvades.' + output_filetype
    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':90}, plot_gridlines = True, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = invert_yaxis), plot_gridlines = False,fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = x_fullname, y = 100*np.array(y), data = dataX_MT, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)

    print 'The figure was saved into {} and {}\n'.format(output_filename_base + '_freq_WT.' + output_filetype, output_filename_base + '_freq_MT.' + output_filetype)

def plot_gameResults_2LT(results_filename, leaky_AAs, max_gamePlayers_num, x, y, title = '', xaxis_label = '', yaxis_label = '', mixed_NashEq_label = 'Mixed', set_minor_xticks = True, invert_xaxis = True, integrate_unknown_NashEq = True, output_filename_base = '', output_filetype = 'pdf'):
    """
    Plots a heatmap of the anaylsis results.
 
    INPUTS:
    -------
    results_filename: 
    A string containing the name of the file containing the results
    (use integrate_results_files function to integrate the results
    in different files

    leaky_AAs:
    A tuple of striings containing the name of amino acids to consider
    e.g., ('lys_L','ile_L')

    max_gamePlayers_num:
    Maximum Number of game players for which the Nash eq should be made. Note that
    the fractions plots is independent of this parameters as it has already taken into
    account all possible game players.

    x & y: 
    A list containing the elements of horizontal and vertical axes 

    xaxis_label & yaxis_label: 
    x-axis and y-axis labels

    mixed_NashEq_label:
    The label of mixed Nash equilibria in the color bar

    integrate_unknown_NashEq:
    If True, all unknowed Nash equilibria are integrated into one case called 'mixed'

    output_filename_base: 
    Base name of the output file containing the plot. This funciton then creates 
    two files output_filename_base + '_games.pdf' and output_filename_base + '_freq.pdf'
    The first stores the plot for the games and the second for species frequencies. 

    output_filetype:
    File type of the output file (e.g., pdf, png, etc). Do not include dot "." here
    """
    import numpy as np
    from tools.utilities.load_data_fromFile import load_data_from_python_file 
    from matplotlib import colors
    from tools.utilities.plot import plot, axis, color_bar

    if not isinstance(set_minor_xticks,bool):
        raise TypeError('set_minor_xticks must be boolean')

    print 'Importing results ...'

    results = load_data_from_python_file(file_name = results_filename, var_names = ['results'])['results']

    # List of amino acids
    AAs_list = ['ala_L', 'arg_L', 'asn_L', 'asp_L', 'cys_L', 'gln_L', 'glu_L', 'gly', 'ile_L', 'leu_L', 'lys_L', 'met_L', 'phe_L', 'phe_L', 'pro_L', 'ser_L', 'thr_L', 'trp_L', 'tyr_L', 'val_L'] 

    print 'Extracting results for the specified pair ...'
    res_keys_to_consider = []   
    for res_key in results.keys():

        AAexchrxns = [s for k in dict(res_key).keys() for s in k]

        # Extract the AA name
        AA_names = sorted(tuple([re.sub('EX_|\(e\)','',exchrxn) for exchrxn in AAexchrxns]))

        if sorted(leaky_AAs) == AA_names:
            res_keys_to_consider.append(res_key)

    #--------------- GAMES ---------------------------
    #--- Find out how many Nash eq we have in the results ---
    for gamePlayers_num in range(2,max_gamePlayers_num+1):
        print '\n-------- Number of game players = ',gamePlayers_num, ' ------------\n'
        print 'Extracting Nash equilibria for the specified pair ...'

        # (mutant, mutant)      --> Prisonder's Dilemma
        # (mutant,wild_type)    --> Snowdirft
        # (wild_type,wild_type) --> Mutually beneficial
        # Since mutant_name varies from one case to the other (depending on the AA), we'll replace
        # all actual mutant names with 'mutant' in order to identify Nash eq easier
        Nash_equilibria = []
    
        # The following is a dictionary where keys are the same results keys and values are Nash 
        # equilibria with above name replacements 
        results_neq = dict([(k,None) for k in results.keys()])
    
        # The following dictionary maps each AA to a general mutant name (to appeat in Nash equilibria)
        AA_mutant_map = dict([(AA, 'mutant' + str(leaky_AAs.index(AA) + 1)) for AA in leaky_AAs])
        print 'AA_mutant_map = ',AA_mutant_map
    
        for res_key in res_keys_to_consider:
            # Nash eq
            Nash_eq = results[res_key]['games_info'][gamePlayers_num]['NashEq'] 
    
            # Only the ones where AA does not appear there
            Nash_eq_AAreplaced = [neq for neq in Nash_eq if len([n for n in neq for AA in leaky_AAs if AA in n]) == 0] 
    
            # Replace AA names with 'mutant1' and 'mutant2'
            for neq in [neq for neq in Nash_eq if len([n for n in neq for AA in leaky_AAs if AA in n]) > 0]: 
                neq = list(neq)
    
                for AA in leaky_AAs:
                    for n in neq: 
                        neq[neq.index(n)] = re.sub(AA,AA_mutant_map[AA],n) 
                Nash_eq_AAreplaced.append(tuple(sorted(neq))) 
    
            Nash_eq_AAreplaced = sorted(Nash_eq_AAreplaced)
    
            if len(Nash_eq_AAreplaced) == 1:
                results_neq[res_key] = Nash_eq_AAreplaced[0]
                Nash_equilibria.append(Nash_eq_AAreplaced[0])
            else:
                results_neq[res_key] = tuple(Nash_eq_AAreplaced)
                Nash_equilibria.append(tuple(Nash_eq_AAreplaced))
                #print '**WARNING! More than one Nash eq for {}: {}'.format(res_key, Nash_eq_AAreplaced)
    
        Nash_equilibria = list(set(Nash_equilibria))
        Nash_eq = Nash_equilibria  # just ot use a shorter name 
    
        # Known Nash equilibira. Note that not all known Nash equilibria may appear in Nash_eq
        # Source: http://stackoverflow.com/questions/30893483/make-a-heatmap-with-a-specified-discrete-color-mapping-with-matplotlib-in-python
        if gamePlayers_num == 2:
            m1m1 = ('mutant1', 'mutant1')       # Red
            m2m2 = ('mutant2', 'mutant2')       # Orange
            m1m2 = ('mutant1', 'mutant2')       # Cyan
            wtwt = ('wild_type', 'wild_type')   # Green
            m1wt = ('mutant1', 'wild_type')     # Pink
            m2wt = ('mutant2', 'wild_type')     # Yellow
            m1m2_wtwt = (('mutant1','mutant2'), ('mutant1_mutant2', 'wild_type'))  # Magenta

            # Any other equilibrium will be blue
            known_Nash_eq = [m1m1, m2m2, m1m2, wtwt, m1wt, m2wt, m1m2_wtwt]

        elif gamePlayers_num == 3:
            m1_m1_m1 = ('mutant1', 'mutant1', 'mutant1')                 # Red
            m2_m2_m2 = ('mutant2', 'mutant2', 'mutant2')                 # Orange
            m1_m2_m1m2 = ('mutant1', 'mutant1_mutant2', 'mutant2')       # Cyan
            wt_wt_wt = ('wild_type', 'wild_type', 'wild_type')           # Green
            m1_m1_wt = ('mutant1', 'mutant1', 'wild_type')               # Pink
            m2_m2_wt = ('mutant2', 'mutant2', 'wild_type')               # Yellow
            m1_m2_m1m2_AND_m1m2_m1m2_wt = (('mutant1', 'mutant1_mutant2', 'mutant2'), ('mutant1_mutant2', 'mutant1_mutant2', 'wild_type'))  # Magenta

            # Any other equilibrium will be blue
            known_Nash_eq = [m1_m1_m1, m2_m2_m2, m1_m2_m1m2, wt_wt_wt, m1_m1_wt, m2_m2_wt, m1_m2_m1m2_AND_m1m2_m1m2_wt]

        elif gamePlayers_num == 4:
            m1_m1_m1_m1 = ('mutant1', 'mutant1', 'mutant1', 'mutant1')                 # Red
            m2_m2_m2_m2 = ('mutant2', 'mutant2', 'mutant2', 'mutant2')                 # Orange
            m1_m2_m1m2_m1m2 = ('mutant1', 'mutant1_mutant2', 'mutant1_mutant2', 'mutant2')       # Cyan
            wt_wt_wt_wt = ('wild_type', 'wild_type', 'wild_type', 'wild_type')           # Green
            m1_m1_m1_wt = ('mutant1', 'mutant1', 'mutant1', 'wild_type')               # Pink
            m2_m2_m2_wt = ('mutant2', 'mutant2', 'mutant2', 'wild_type')               # Yellow
            m1_m2_m1m2_m1m2_AND_m1m2_m1m2_m1m2_wt = (('mutant1', 'mutant1_mutant2', 'mutant1_mutant2', 'mutant2'), ('mutant1_mutant2', 'mutant1_mutant2', 'mutant1_mutant2', 'wild_type'))  # Magenta

            # Any other equilibrium will be blue
            known_Nash_eq = [m1_m1_m1_m1, m2_m2_m2_m2, m1_m2_m1m2_m1m2, wt_wt_wt_wt, m1_m1_m1_wt, m2_m2_m2_wt, m1_m2_m1m2_m1m2_AND_m1m2_m1m2_m1m2_wt] 

        unknown_Nash_eq = [n for n in Nash_eq if n not in known_Nash_eq]
    
        Nash_eq = [n for n in Nash_eq if n in known_Nash_eq] + unknown_Nash_eq 
     
        # A dictionary mapping known Nash equilibria to a name and a color and a data value (data value is used for plotting)
        # RGBs for various gray colors can be found here: http://www.tayloredmktg.com/rgb/ 
        if integrate_unknown_NashEq:
            if gamePlayers_num == 2:
                neq_color_map = {m1m1:'Red', m2m2:'Orange', m1m2: 'Cyan', wtwt:'Green', m1wt: 'Pink', m2wt:'Yellow', m1m2_wtwt:'Magenta', 'mixed':np.array([49,79,79])/255}
                neq_name_map = {m1m1:'(m1,m1)', m2m2:'(m2,m2)', m1m2: '(m1,m2)', wtwt:'(wt,wt)', m1wt: '(m1,wt)', m2wt:'(m2,wt)', m1m2_wtwt:'((m1,m2)+\n(m1m2,wt))', 'mixed':mixed_NashEq_label}
            elif gamePlayers_num == 3:
                neq_color_map = {m1_m1_m1:'Red', m2_m2_m2:'Orange', m1_m2_m1m2: 'Cyan', wt_wt_wt:'Green', m1_m1_wt: 'Pink', m2_m2_wt:'Yellow', m1_m2_m1m2_AND_m1m2_m1m2_wt:'Magenta', 'mixed':np.array([49,79,79])/255}
                neq_name_map = {m1_m1_m1:'(m1,m1,m1)', m2_m2_m2:'(m2,m2,m2)', m1_m2_m1m2: '(m1,m2,m1m2)', wt_wt_wt:'(wt,wt,wt)', m1_m1_wt: '(m1,m1,wt)', m2_m2_wt:'(m2,m2,wt)', m1_m2_m1m2_AND_m1m2_m1m2_wt:'((m1,m2,m1m2)+\n(m1m2,m1m2,wt))', 'mixed':mixed_NashEq_label}
            elif gamePlayers_num == 4:
                neq_color_map = {m1_m1_m1_m1:'Red', m2_m2_m2_m2:'Orange', m1_m2_m1m2_m1m2: 'Cyan', wt_wt_wt_wt:'Green', m1_m1_m1_wt: 'Pink', m2_m2_m2_wt:'Yellow', m1_m2_m1m2_m1m2_AND_m1m2_m1m2_m1m2_wt:'Magenta', 'mixed':np.array([49,79,79])/255}
                neq_name_map = {m1_m1_m1_m1:'(m1,m1,m1,m1)', m2_m2_m2_m2:'(m2,m2,m2,m2)', m1_m2_m1m2_m1m2: '(m1,m2,m1m2,m1m2)', wt_wt_wt_wt:'(wt,wt,wt,wt)', m1_m1_m1_wt: '(m1,m1,m1,wt)', m2_m2_m2_wt:'(m2,m2,m2,wt)', m1_m2_m1m2_m1m2_AND_m1m2_m1m2_m1m2_wt:'((m1,m2,m1m2,m1m2)+\n(m1m2,m1m2,m1m2,wt))', 'mixed':mixed_NashEq_label}
            
        else:
            if gamePlayers_num == 2:
                neq_color_map = {m1m1:'Red', m2m2:'Orange', m1m2: 'Cyan', wtwt:'Green', m1wt: 'Pink', m2wt:'Yellow', m1m2_wtwt:'Magenta'}
                neq_name_map = {m1m1:'(m1,m1)', m2m2:'(m2,m2)', m1m2: '(m1,m2)', wtwt:'(wt, wt)', m1wt: '(m1, wt)', m2wt:'(m2, wt)', m1m2_wtwt:'((m1,m2)+(m1_m2,wt))'}
            elif gamePlayers_num == 3:
                neq_color_map = {m1_m1_m1:'Red', m2_m2_m2:'Orange', m1_m2_m1m2: 'Cyan', wt_wt_wt:'Green', m1_m1_wt: 'Pink', m2_m2_wt:'Yellow', m1_m2_m1m2_AND_m1m2_m1m2_wt:'Magenta'}
                neq_name_map = {m1_m1_m1:'(m1,m1,m1)', m2_m2_m2:'(m2,m2,m2)', m1_m2_m1m2: '(m1,m2,m1m2)', wt_wt_wt:'(wt,wt,wt)', m1_m1_wt: '(m1,m1,wt)', m2_m2_wt:'(m2,m2,wt)', m1_m2_m1m2_AND_m1m2_m1m2_wt:'((m1,m2,m1m2)+\n(m1m2,m1m2,wt))'}
            elif gamePlayers_num == 4:
                neq_color_map = {m1_m1_m1_m1:'Red', m2_m2_m2_m2:'Orange', m1_m2_m1m2_m1m2: 'Cyan', wt_wt_wt_wt:'Green', m1_m1_m1_wt: 'Pink', m2_m2_m2_wt:'Yellow', m1_m2_m1m2_m1m2_AND_m1m2_m1m2_m1m2_wt:'Magenta'}
                neq_name_map = {m1_m1_m1_m1:'(m1,m1,m1,m1)', m2_m2_m2_m2:'(m2,m2,m2,m2)', m1_m2_m1m2_m1m2: '(m1,m2,m1m2,m1m2)', wt_wt_wt_wt:'(wt,wt,wt,wt)', m1_m1_m1_wt: '(m1,m1,m1,wt)', m2_m2_m2_wt:'(m2,m2,m2,wt)', m1_m2_m1m2_m1m2_AND_m1m2_m1m2_m1m2_wt:'((m1,m2,m1m2,m1m2)+\n(m1m2,m1m2,m1m2,wt))'}
    
            # Color and name for unknonw Nash equilibria    
            counter = 0
            gray_rgbs = [np.array([49,79,79])/255, np.array([105,105,105])/255, np.array([112, 138, 144])/255, np.array([190, 190, 190])/255, np.array([211, 211, 211])/255, np.array([238, 233, 233])/255, np.array([205, 201, 201])/255, np.array([220,220,220])/255, np.array([255,222,173])/255, np.array([238,203,173])/255, np.array([193,205,193])/255, np.array([100,149,237])/255]
            for neq in unknown_Nash_eq:
                neq_color_map[neq] = gray_rgbs[counter]
                neq_name_map[neq] = mixed_NashEq_label + str(counter + 1) 
                counter += 1
                print neq_name_map[neq],' --> ',neq
                
        print '\nNash equilibria of the systems are:'
        print 'Known:'
        for neq in [n for n in Nash_eq if n in known_Nash_eq]:
            print neq
        print '\nUnknown:'
        for neq in unknown_Nash_eq:
            print neq
    
        data_value = {}
    
        if integrate_unknown_NashEq:
            # Number of known Nash equilibria (m1m1,m1m2, m1m2, wtwt, m1wt,m2wt)
            known_NE_num = len([kneq for kneq in known_Nash_eq if kneq in Nash_eq])
        
            if known_NE_num == 0:
                colormap = colors.ListedColormap([neq_color_map['mixed']])
                colorbar_ticklabels = [neq_name_map['mixed']]
                data_value['mixed'] = 0
            else:
                # Find which known Nash equilibria appear in Nash_eq
                known_neq = [kneq for kneq in known_Nash_eq if kneq in Nash_eq]
                if len(Nash_eq) == known_NE_num:
                    colormap = colors.ListedColormap([neq_color_map[kneq] for kneq in known_neq])
                    colorbar_ticklabels = [neq_name_map[kneq] for kneq in known_eq]
                    for kneq in known_eq:
                        data_value[neq_name_map[neq]] = known_neq.index(kneq)
    
                elif len(Nash_eq) > known_NE_num:
                    colormap = colors.ListedColormap([neq_color_map[kneq] for kneq in known_neq] + [neq_color_map['mixed']])
                    colorbar_ticklabels = [neq_name_map[kneq] for kneq in known_neq] + [neq_name_map['mixed']]
                    for kneq in known_neq:
                        data_value[neq_name_map[kneq]] = known_neq.index(kneq)
                    data_value['mixed'] = len(known_neq)
                else:
                     raise userError('len(Nash_eq) = {} is less than known_NE_num = {}'.format(len(Nash_eq), known_NE_num)) 
        else:
            colormap = colors.ListedColormap([neq_color_map[neq] for neq in Nash_eq])
            colorbar_ticklabels = [neq_name_map[neq] for neq in Nash_eq]
            for neq in Nash_eq:
               data_value[neq_name_map[neq]] = Nash_eq.index(neq)
    
        # Games
        dataG = np.zeros((len(y),len(x)))
    
        # Here, x and y are leakiness levels for AA1 and AA2 
        for i, yy in enumerate(y):
            for j, xx in enumerate(x):
                # The dictionary key in results (this works only when we have one leaky trait) 
                if ((('EX_' + leaky_AAs[0] + '(e)',), xx), (('EX_' + leaky_AAs[1] + '(e)',), yy)) in results_neq.keys(): 
                    res_key = ((('EX_' + leaky_AAs[0] + '(e)',), xx), (('EX_' + leaky_AAs[1] + '(e)',), yy))
                elif ((('EX_' + leaky_AAs[1] + '(e)',), yy), (('EX_' + leaky_AAs[0] + '(e)',), xx)) in results_neq.keys(): 
                    res_key = ((('EX_' + leaky_AAs[1] + '(e)',), yy), (('EX_' + leaky_AAs[0] + '(e)',), xx))
                else:
                    raise userError('leaky_AAs = {} , x = {},  y = {} was not found in results.keys()'.format(leaky_AAs, xx, yy))
    
                if integrate_unknown_NashEq:
                    # Find out whether results_neq[res_key] is a knonw Nash equilibirum
                    neqs = [kneq for kneq in known_neq if kneq == results_neq[res_key]]
                    if len(neqs) == 1:
                        dataG[i,j] = data_value[neq_name_map[neqs[0]]]
                    elif len(neqs) == 0:
                        dataG[i,j] = data_value['mixed']
                    else:
                        raise userError('len(neqs) = {} is greater than 1'.format(len(neqs)))
                else: 
                    neqs = [neq for neq in Nash_eq if neq == results_neq[res_key]]
                    if len(neqs) == 1:
                        dataG[i,j] = data_value[neq_name_map[neqs[0]]]
                    else:
                        raise userError('len(neqs) = {} is not 1'.format(len(neqs)))
    
        print '\nmin(data) = {}  , max(data) = {}'.format(dataG.min(), dataG.max())
    
        #-- Create a color map --
        xaxis_ticklabels = [str(int(100*xx)) if 100*xx/20 == int(100*xx/20) else '' for xx in x]
        yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) else '' for yy in y]

        if integrate_unknown_NashEq:
            output_filename = output_filename_base + '_NEs' + str(gamePlayers_num) + 'players.' + output_filetype
        else:
            output_filename = output_filename_base + '_NEs' + str(gamePlayers_num) + 'players_notMixed.' + output_filetype
    
        curr_plt = plot(title = 'Number of players = ' + str(gamePlayers_num), xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
        curr_plt.heatmap(x = np.array(x), y = 100*np.array(y), data = dataG, plot_func = 'matshow', clrbar = color_bar(colormap = colormap, colorlimits = (dataG.min() ,dataG.max()), label = '',  custom_ticklabels =  colorbar_ticklabels, ticklabels_format = {'ticks_position':'middle'}))
    
        print 'The figure was saved into {}\n'.format(output_filename)
    
    #--------------- Evolutionary dynamics ---------------------------
    # Species frequences when wild-type invades (dataX_WT) or mutants invade (dataX_MT)
    dataX_WT_wt = np.zeros((len(y),len(x)))   # Wild-type frequency when wild_type invades
    dataX_WT_m1 = np.zeros((len(y),len(x)))   # mutant1 frequency when wild-type invades
    dataX_WT_m2 = np.zeros((len(y),len(x)))   # mutant2 frequency when wild-type invades
    dataX_WT_m1m2 = np.zeros((len(y),len(x))) # mutant1_mutant2 frequency when wild-type invades
    dataX_MT_wt = np.zeros((len(y),len(x)))   # Wild-type when mutans invade 
    dataX_MT_m1 = np.zeros((len(y),len(x)))   # mutant1 frequency when mutants invade
    dataX_MT_m2 = np.zeros((len(y),len(x)))   # mutant2 frequency when mutants invade
    dataX_MT_m1m2 = np.zeros((len(y),len(x))) # mutant1_mutant2 frequency when mutants invade

    # Here, x is the amino acid (or cost) and y is capture efficiency
    for i, yy in enumerate(y):
        for j, xx in enumerate(x):
            # The dictionary key in results (this works only when we have one leaky trait) 
            if ((('EX_' + leaky_AAs[0] + '(e)',), xx), (('EX_' + leaky_AAs[1] + '(e)',), yy)) in results_neq.keys(): 
                res_key = ((('EX_' + leaky_AAs[0] + '(e)',), xx), (('EX_' + leaky_AAs[1] + '(e)',), yy))
            elif ((('EX_' + leaky_AAs[1] + '(e)',), yy), (('EX_' + leaky_AAs[0] + '(e)',), xx)) in results_neq.keys(): 
                res_key = ((('EX_' + leaky_AAs[1] + '(e)',), yy), (('EX_' + leaky_AAs[0] + '(e)',), xx))
            else:
                raise userError('leaky_AAs = {} , x = {},  y = {} cannot be found in results.keys()'.format(leaky_AAs, xx, yy))

            # Wild-type invades
            game_key = [gk for gk in results[res_key]['replicator_dynamics_info'].keys() if ('wild_type', 0.01) in gk]
            if len(game_key) == 1: 
                # Wild-type
                dataX_WT_wt[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]]['wild_type']

                # Mutant 1
                dataX_WT_m1[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]][leaky_AAs[0]]

                # Mutant 2
                dataX_WT_m2[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]][leaky_AAs[1]]

                # (Mutant1, Mutant 2)
                if '_'.join((leaky_AAs[0], leaky_AAs[1])) in results[res_key]['replicator_dynamics_info'][game_key[0]].keys():
                    dataX_WT_m1m2[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]]['_'.join([leaky_AAs[0],leaky_AAs[1]])]
                elif '_'.join((leaky_AAs[1], leaky_AAs[0])) in results[res_key]['replicator_dynamics_info'][game_key[0]].keys():
                    dataX_WT_m1m2[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]]['_'.join([leaky_AAs[1],leaky_AAs[0]])]
                else:
                    raise userError("Neither {} nor {} is in results[res_key]['replicator_dynamics_info'][game_key[0]].keys() = {}".format('_'.join([leaky_AAs[0],leaky_AAs[1]]), '_'.join([leaky_AAs[1],leaky_AAs[0]]), results[res_key]['replicator_dynamics_info'][game_key[0]].keys()))
            else:
                raise userError('len(game_key) > 1, where game_key = {}'.format(game_key))

            # Mutant invades
            game_key = [gk for gk in results[res_key]['replicator_dynamics_info'].keys() if ('wild_type', 0.99) in gk]
            if len(game_key) == 1: 
                # WIld-type
                dataX_MT_wt[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]]['wild_type']

                # Mutant 1
                dataX_MT_m1[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]][leaky_AAs[0]]

                # Mutant 2
                dataX_MT_m2[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]][leaky_AAs[1]]

                # (Mutant1, Mutant 2)
                if '_'.join((leaky_AAs[0], leaky_AAs[1])) in results[res_key]['replicator_dynamics_info'][game_key[0]].keys():
                    dataX_MT_m1m2[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]]['_'.join([leaky_AAs[0],leaky_AAs[1]])]
                elif '_'.join((leaky_AAs[1], leaky_AAs[2])) in results[res_key]['replicator_dynamics_info'][game_key[0]].keys():
                    dataX_MT_m1m2[i,j] = results[res_key]['replicator_dynamics_info'][game_key[0]]['_'.join([leaky_AAs[1],leaky_AAs[0]])]
                else:
                    raise userError("Neither {} nor {} is in results[res_key]['replicator_dynamics_info'][game_key[0]].keys() = {}".format('_'.join([leaky_AAs[0],leaky_AAs[1]]), '_'.join([leaky_AAs[1],leaky_AAs[0]]), results[res_key]['replicator_dynamics_info'][game_key[0]].keys()))
            else:
                raise userError('len(game_key) > 1, where game_key = {}'.format(game_key))

    print 'min(dataX_WT_wt) = {}  , max(dataX_WT_wt) = {}'.format(dataX_WT_wt.min(), dataX_WT_wt.max())
    print 'min(dataX_MT_wt) = {}  , max(dataX_MT_wt) = {}'.format(dataX_MT_wt.min(), dataX_MT_wt.max())

    print 'The figures were saved into the following files:\n'

    custom_ticks = np.arange(0,1+0.2,0.2)

    # WT invades, wt
    output_filename_base + '_freq_WTInvades_wt.' + output_filetype
    curr_plt = plot(title = title + ' Wild-type (WT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_WT_wt, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_WT_wt.' + output_filetype

    # WT invades, m1
    output_filename_base + '_freq_WTInvades_' + leaky_AAs[0] + '.' + output_filetype
    curr_plt = plot(title = title + ' ' + leaky_AAs[0] + ' (WT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_WT_m1, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_WT_' + leaky_AAs[0] + '.' + output_filetype

    # WT invades, m2
    output_filename_base + '_freq_WTInvades_' + leaky_AAs[1] + '.' + output_filetype
    curr_plt = plot(title = title + ' ' + leaky_AAs[1] + ' (WT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_WT_m2, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_WT_' + leaky_AAs[1] + '.' + output_filetype

    # WT invades, m1_m2
    output_filename = output_filename_base + '_freq_WTInvades_' + leaky_AAs[0] + '_' + leaky_AAs[1] + '.' + output_filetype
    curr_plt = plot(title = title + ' ' + '_'.join(leaky_AAs) + ' (WT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_WT_m1m2, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_WT_' + leaky_AAs[0] + '_' + leaky_AAs[1] + '.' + output_filetype

    # MT invades, wt
    output_filename = output_filename_base + '_freq_MTInvades_wt.' + output_filetype
    curr_plt = plot(title = title + ' WIld-type (MT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False,fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_MT_wt, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_MT_wt.' + output_filetype

    # MT invades, m1 
    output_filename = output_filename_base + '_freq_MT_' + leaky_AAs[0] + '.' + output_filetype
    curr_plt = plot(title = title + ' '+ leaky_AAs[0] + ' (MT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False,fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_MT_m1, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_MT_' + leaky_AAs[0] + '.' + output_filetype

    # MT invades, m2 
    output_filename = output_filename_base + '_freq_MTInvades_' + leaky_AAs[1] + '.' + output_filetype
    curr_plt = plot(title = title + ' ' + leaky_AAs[1] + ' (MT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False,fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_MT_m2, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_MT_' + leaky_AAs[1] + '.' + output_filetype

    # MT invades, m1_m2 
    output_filename = output_filename_base + '_freq_MTInvades_' + leaky_AAs[0] + '_' + leaky_AAs[1] + '.' + output_filetype
    curr_plt = plot(title = title + ' ' + '_'.join(leaky_AAs) + ' (MT invades)', xaxis = axis(label = xaxis_label, custom_ticklabels = xaxis_ticklabels, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, plot_gridlines = False, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = True), plot_gridlines = False,fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = 100*np.array(x), y = 100*np.array(y), data = dataX_MT_m1m2, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print output_filename_base + '_freq_MT_' + leaky_AAs[0] + '_' + leaky_AAs[1] + '.' + output_filetype

def plot_AACostsGains(results_filename, x, y, invert_xaxis = False, invert_yaxis = False, output_filename_base = '', output_filetype = 'pdf'):
    """
    Plots a heatmap for the amino acids costs and benefits 
 
    INPUTS:
    -------
    results_filename: 
    A string containing the name of the file containing the results
    (use integrate_results_files function to integrate the results
    in different files

    x & y: 
    A list containing the elements of horizontal and vertical axes 

    output_filename_base: 
    Base name of the output file containing the plot. This funciton then creates 
    two files output_filename_base + '_games.pdf' and output_filename_base + '_freq.pdf'
    The first stores the plot for the games and the second for species frequencies. 

    output_filetype:
    File type of the output file (e.g., pdf, png, etc). Do not include dot "." here
    """
    import numpy as np
    from tools.utilities.load_data_fromFile import load_data_from_python_file 
    from matplotlib import colors
    from tools.utilities.plot import plot, axis, color_bar

    print 'Importing results ...'

    results = load_data_from_python_file(file_name = results_filename, var_names = ['cost_results','gain_results'])
    cost_results = results['cost_results']
    gain_results = results['gain_results']

    # Min cost other than -100
    max_cost = max([cost_results[leak_level][aa] for leak_level in cost_results.keys() for aa in cost_results[leak_level].keys() if cost_results[leak_level][aa] != 100])
    min_gain = min([gain_results[uptake_level][aa] for uptake_level in gain_results.keys() for aa in gain_results[uptake_level].keys() if gain_results[uptake_level][aa] != -100])

    print '\nmin cost (without 100) = {} ,  max gain (without -100) = {}'.format(max_cost, min_gain)

    # Species frequences when wild-type invades (dataX_WT) or mutants invade (dataX_MT)
    data_costs = np.zeros((len(y),len(x)))   # Wild-type frequency when wild_type invades
    data_gains = np.zeros((len(y),len(x)))   # mutant1 frequency when wild-type invades

    # Here, x is the amino acid (or cost) and y is capture efficiency
    for i, yy in enumerate(y):
        for j, xx in enumerate(x):

            if cost_results[yy][xx] != 100:
                data_costs[i,j] = cost_results[yy][xx]
            else:
                data_costs[i,j] = max_cost + 0.01 

            if gain_results[yy][xx] != -100:
                data_gains[i,j] = gain_results[yy][xx]
            else:
                data_gains[i,j] = min_gain - 0.01
            print '{} , {} , {}'.format(xx,yy, data_gains[i,j])

        print '\n'

    print 'min(data_costs) = {}  , max(data_costs) = {}'.format(data_costs.min(), data_costs.max())
    print 'min(dataX_gains) = {}  , max(data_gains) = {}'.format(data_gains.min(), data_gains.max())

    yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) or yy == 0.01 else '' for yy in y]

    colorbar_custom_ticks = [data_costs.min()] + [(data_costs.min() + data_costs.max())/2]+ [data_costs.max()] 

    title = ''
    xaxis_label = 'Amino acids'
    yaxis_label = 'Leakiness level (%)'
    output_filename = output_filename_base + '_costs.' + output_filetype
    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':90}, plot_gridlines = True, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = invert_yaxis), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = x, y = 100*np.array(y), data = data_costs, plot_func = 'matshow', clrbar = color_bar(colorlimits = (data_costs.min(), data_costs.max()), custom_ticks = colorbar_custom_ticks, label = 'Growth cost', label_format = {'distance_from_ticklabels':30}), interpolate = True)

    colorbar_custom_ticks = [data_gains.min()] + [(data_gains.min() + data_gains.max())/2]+ [data_gains.max()] 
    title = ''
    xaxis_label = 'Amino acids'
    yaxis_label = 'Uptake level (%)'
    output_filename = output_filename_base + '_gains.' + output_filetype
    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':90}, plot_gridlines = True, gridlines_format = {'linestyle':'solid', 'linewidth':0.3}), yaxis = axis(label = yaxis_label, custom_ticklabels = yaxis_ticklabels, invert = invert_yaxis), plot_gridlines = False, fig_format = {'figsize':(8,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = x, y = 100*np.array(y), data = data_gains, plot_func = 'matshow', clrbar = color_bar(colorlimits = (data_gains.min(), data_gains.max()), custom_ticks = colorbar_custom_ticks, label = 'Growth advantage', label_format = {'distance_from_ticklabels':30}), interpolate = True)


def run_plot_gameResults(leaky_traits_num = 1):
    """
    Runs run_plot_gameResults
    """
    #---- Amino acids costs and gains ------
    if leaky_traits_num == 0:
        # Here, x is a list of amino acids ordered according to their growth cost
        x = ['ala', 'val', 'asp', 'gly', 'ser', 'asn', 'glu', 'leu', 'pro', 'thr', 'gln', 'ile', 'lys', 'cys', 'phe', 'arg', 'tyr', 'met', 'his', 'trp'] 

        # Leakiness or uptake level
        y = [0.01] + [i/100 for i in range(5,101,5)] 

        plot_AACostsGains(results_filename = 'results/amino_acids_costs_gains.py', x = x, y = y, invert_xaxis = False, invert_yaxis = True, output_filename_base = 'results/figures/AA', output_filetype = 'pdf')

    #---- One leaky trait ------
    if leaky_traits_num == 1:
        # Here, x is a list of amino acids ordered according to their growth cost
        x = ['ala_L', 'val_L', 'asp_L', 'gly', 'ser_L', 'asn_L', 'glu_L', 'leu_L', 'pro_L', 'thr_L', 'gln_L', 'ile_L', 'lys_L', 'cys_L', 'phe_L', 'arg_L', 'tyr_L', 'met_L', 'his_L', 'trp_L'] 

        # List of AAs sorted according to the growth rate of an AA auxotroph mutant strain taking up one unit of that AA
        # in the order of decreasing growth rate
        #x = ['gln_L', 'asn_L', 'asp_L', 'thr_L', 'arg_L', 'ser_L', 'gly', 'lys_L', 'leu_L', 'phe_L', 'met_L', 'ile_L', 'cys_L', 'tyr_L', 'his_L', 'pro_L', 'val_L', 'trp_L', 'ala_L', 'glu_L'] 

        # List of amino acids sorted according to cost to benefit ratio
        #x = ['ala_L', 'val_L', 'asp_L', 'ser_L', 'gly', 'asn_L', 'leu_L', 'thr_L', 'gln_L', 'pro_L', 'ile_L', 'lys_L', 'glu_L', 'cys_L', 'phe_L', 'arg_L', 'tyr_L', 'met_L', 'his_L', 'trp_L']

        # List of amino acids sorted according to cost + benefit
        #x = ['glu_L', 'ala_L', 'val_L', 'gly', 'ser_L', 'leu_L', 'pro_L', 'asp_L', 'ile_L', 'asn_L', 'thr_L', 'lys_L', 'cys_L', 'gln_L', 'phe_L', 'tyr_L', 'met_L', 'arg_L', 'his_L', 'trp_L'] 

        #results_filename = 'results/leakiness_games_BQH_1LT_v1.py'
        results_filename = 'results/leakiness_games_BQH_1LT.py'

        # Leakiness levels
        y = [i/100 for i in range(0,101)]
 
        plot_gameResults_1LT(results_filename = 'results/leakiness_games_BQH_1LT.py', number_of_gamePlayers = 2, x = x, y = y, title = '', xaxis_label = 'Amino acids (in increasing growth cost)', yaxis_label = 'Leakiness level (%)', mixed_NashEq_label = "Mixed\n(PD + SD)",set_minor_xticks = True, invert_xaxis = False, invert_yaxis = True, output_filename_base = 'results/figures/BQH/leakiness_games_BQH_1LT')

    #---- Two leaky trait ------
    if leaky_traits_num == 2:
        results_filename = 'results/leakiness_games_BQH_2LT_all.py'
        #results_filename = 'results/leakiness_games_BQH_2LT_selected_old.py' 

        max_gamePlayers_num = 3
       
        # Leakiness levels
        x = [i/100 for i in range(0,105,5)]
        y = [i/100 for i in range(0,105,5)]
 
        # (L-lys, L-ile)
        #plot_gameResults_2LT(results_filename = results_filename, leaky_AAs = ('lys_L','ile_L'), max_gamePlayers_num = max_gamePlayers_num, x = x, y = y, title = '', xaxis_label = 'L-lysine leakiness level (%)', yaxis_label = 'L-isoleucine leakiness \nlevel (%)', mixed_NashEq_label = 'Mixed', set_minor_xticks = True, invert_xaxis = False, integrate_unknown_NashEq = True, output_filename_base = 'results/figures/BQH/leakiness_games_BQH_2LT_lysIle', output_filetype = 'pdf')

        # (L-ala, L-val)
        #plot_gameResults_2LT(results_filename = results_filename, leaky_AAs = ('ala_L','val_L'), max_gamePlayers_num = max_gamePlayers_num, x = x, y = y, title = '', xaxis_label = 'L-alanine leakiness level(%)', yaxis_label = 'L-valine leakiness level (%)', mixed_NashEq_label = 'Mixed', set_minor_xticks = True, invert_xaxis = False, integrate_unknown_NashEq = True, output_filename_base = 'results/figures/BQH/leakiness_games_BQH_2LT_alaVal', output_filetype = 'pdf')

        # (L-met, L-trp)
        #plot_gameResults_2LT(results_filename = results_filename, leaky_AAs = ('met_L','trp_L'), max_gamePlayers_num = max_gamePlayers_num, x = x, y = y, title = '', xaxis_label = 'L-methionine leakiness level(%)', yaxis_label = 'L-tryptophan leakiness level (%)', mixed_NashEq_label = 'Mixed', set_minor_xticks = True, invert_xaxis = True, integrate_unknown_NashEq = True, output_filename_base = 'results/figures/BQH/leakiness_games_BQH_2LT_metTrp', output_filetype = 'pdf')

        # (L-val, L-trp)
        plot_gameResults_2LT(results_filename = results_filename, leaky_AAs = ('val_L','trp_L'), max_gamePlayers_num = max_gamePlayers_num, x = x, y = y, title = '', xaxis_label = 'L-methionine leakiness level(%)', yaxis_label = 'L-tryptophan leakiness level (%)', mixed_NashEq_label = 'Mixed', set_minor_xticks = True, invert_xaxis = True, integrate_unknown_NashEq = True, output_filename_base = 'results/figures/BQH/leakiness_games_BQH_2LT_valTrp', output_filetype = 'pdf')


def integrate_pair_results():
    from tools.utilities.integrate_results_files import integrate_results_files, del_results_files

    results_filenames = ['results/leakiness_games_BQH_2LT_1_10.py', 'results/leakiness_games_BQH_2LT_11_20.py', 'results/leakiness_games_BQH_2LT_21_30.py', 'results/leakiness_games_BQH_2LT_31_40.py', 'results/leakiness_games_BQH_2LT_41_50.py', 'results/leakiness_games_BQH_2LT_51_60.py', 'results/leakiness_games_BQH_2LT_61_70.py', 'results/leakiness_games_BQH_2LT_71_80.py', 'results/leakiness_games_BQH_2LT_81_90.py', 'results/leakiness_games_BQH_2LT_91_100.py', 'results/leakiness_games_BQH_2LT_101_110.py', 'results/leakiness_games_BQH_2LT_111_120.py', 'results/leakiness_games_BQH_2LT_121_130.py', 'results/leakiness_games_BQH_2LT_131_140.py', 'results/leakiness_games_BQH_2LT_141_150.py', 'results/leakiness_games_BQH_2LT_151_160.py', 'results/leakiness_games_BQH_2LT_161_170.py', 'results/leakiness_games_BQH_2LT_171_171.py']   

    integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/leakiness_games_BQH_2LT_all.py')

    #del_results_files(results_filenames)

def run_rankAAs():
    """
    Ranks amino acids according to (1) the growth cost of producing them and (ii) the growth advtangae of amino acid
    auxotroph mutants taking up the missing amino acid for various levels of amino acid leakiness and uptake levels
    """
    from rankAAs import rankAAs_byGrowthCost, rank_mutants_byGrowthGain

    max_uptake_leakiness = 10

    # Leakiness and uptake levels percents
    levels_percent = [0.01] + [i/100 for i in range(5,101,5)] 

    AAsCosts = {}
    AAsGains = {}

    with open('results/amino_acids_costs_gains.py','w') as f:
        f.write('cost_results = {}\n')
        f.write('gain_results = {}\n')
        for level_perc in levels_percent:
            AAsCosts_list = rankAAs_byGrowthCost(leakiness_level = level_perc*max_uptake_leakiness, stdout_msgs = False, stdout_msgs_details = False)
            AAsCosts[level_perc] = dict(AAsCosts_list)
            f.write('cost_results[' + str(level_perc) + '] = ' + str(AAsCosts[level_perc]) + '\n')

            AAsGains_list = rank_mutants_byGrowthGain(uptake_level = level_perc*max_uptake_leakiness, stdout_msgs = False, stdout_msgs_details = False)
            AAsGains[level_perc] = dict(AAsGains_list)
            f.write('gain_results[' + str(level_perc) + '] = ' + str(AAsGains[level_perc]) + '\n')
           
 
def test():
    #--- E. coli iJO1366 model ---
    print '\n--- Wild-type E.coli (iJO1366 model) ----'
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'
    growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iJO1366_minimal_glucose_anaerobic.py', 'flux_bounds_dict': {'EX_glc(e)':[-10,1000]}}
    WT = create_model(model_organism = organism(id = 'Ecoli', name = 'Escherichia coli'), model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = growthMedium_flux_bounds) 

    #--- Load the list of rxns that must be off in each mutant ---
    from mutants_rxn_info_iJO1366 import mutants_rxn_info_AAs, genes_AA_map

    # Growth condition
    set_specific_bounds(model = WT, flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'], file_name = growthMedium_flux_bounds['flux_bounds_filename'], reset_flux_bounds = True)    

    mutant = 'gltBD_gdhA' 
    exchrxn_name = 'EX_glu_L(e)'
    exchrxn_lb = 10*0.6
       
    # Knockouts 
    set_specific_bounds(model = WT, flux_bounds = dict([(rxn_id,[0,0]) for rxn_id in mutants_rxn_info_AAs[mutant]]), reset_flux_bounds = False)    
    set_specific_bounds(model = WT, flux_bounds = {exchrxn_name:[-exchrxn_lb,None]}, reset_flux_bounds = False)    
        
    WT.fba() 

#----------------------------------------------------
if __name__ == '__main__':

    #test()

    #------- Two leaky traits -------
    #master_func_BQH(start_pos = 15, end_pos = 15, leaky_traits_num = 1, max_game_players_num = 2, leakiness_levels = [0.3], simulate_rep_dynamics = True, tf = 5000, dt = 1, save_details = True, results_filename = 't.py', stdout_msgs = True, stdout_msgs_details = True, stdout_msgs_fba = True, warnings = True)

    #------- Two leaky traits -------
    #master_func_BQH(start_pos = 1, end_pos = 1, leaky_traits_num = 1, max_game_players_num = 2, leakiness_levels = [0.005], results_filename_base = '', stdout_msgs = True, warnings = True)
    # (lys_L, ile_L) --> 121   (ala_L, val_L) --> 8     ('metC_malY', 'trpA') --> 161 
    # ('hisB', 'lysA') --> 41,   ('alaA_alaC_avtA', 'trpA') --> 17
    master_func_BQH(start_pos = 8, end_pos = 8, leaky_traits_num = 2, max_game_players_num = 4, leakiness_levels = [0.1], simulate_rep_dynamics = True, tf = 5000, dt = 1, save_details = False, results_filename = '', stdout_msgs = True, stdout_msgs_details = False, warnings = True)

    #------- Three leaky traits -------
    # ('alaA_alaC_avtA' (ala_L), 'lysA', 'ilvE') --> 102 (**Do NOT consider. ile_L_ala_L auxotroph for val_L)
    # ('lysA', 'thrC', 'ilvE') --> 902 
    # ('lysA', 'ilvE', 'trpA') --> 916,  ('lysA', 'ilvE_avtA', 'ilvE') --> 862
    #master_func_BQH(start_pos = 862, end_pos = 862, leaky_traits_num = 3, max_game_players_num = 4, leakiness_levels = [0.1], simulate_rep_dynamics = True, tf = 5000, dt = 1, save_details = False, results_filename = '', stdout_msgs = True, stdout_msgs_details = False, stdout_msgs_fba = False, warnings = True)

