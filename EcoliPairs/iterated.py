from __future__ import division
import re, sys, math, copy, time, random
import numpy as np
sys.path.append('../')
import os.path
import cPickle as pk
from copy import deepcopy
from tools.showTime import *
from models import *
from tools.FBA import *
from tools.GAMETES.game import *
from tools.GAMETES.iteratedGame import iteratedGame
from tools.GAMETES.stochasticStrategy import stochasticStrategy
from tools.customError import *
from tools.importData import *

# ------------------------------------
# Compares the experimentally observed co-growth with the results from iterated games 
# ------------------------------------

if __name__ == "__main__":

    startT = time.clock()
    print '\n**Job started at ',time.strftime("%Y-%m-%d %H:%M:%S"),'\n'

    #--- Load the games --- 
    if os.path.isfile('results/auxoEcoliGames_all.pk') == False:
        # Concatanate the two files
        with open('results/auxoEcoliGames_0_20.pk','rb') as inputFile:
            games1 = pk.load(inputFile)             
        with open('results/auxoEcoliGames_21_end.pk','rb') as inputFile:
            games2 = pk.load(inputFile)             
        games = dict(games1.items() + games2.items())
        with open('results/auxoEcoliGames_all.pk','wb') as outputFile:
            pk.dump(games,outputFile,-1)
    else:
        with open('results/auxoEcoliGames_all.pk','rb') as inputFile:
            games = pk.load(inputFile)

    print 'Loading experimental data ...\n'
    #--- Load the experimental flux data ---
    day1Rep1Inst = importData(inputFile = 'expData/day1Rep1.txt',delType = 'tab',dataType = 'float')
    day1Rep1 = day1Rep1Inst.run()

    day1Rep2Inst = importData(inputFile = 'expData/day1Rep2.txt',delType = 'tab',dataType = 'float')
    day1Rep2 = day1Rep2Inst.run()

    day4Rep1Inst = importData(inputFile = 'expData/day4Rep1.txt',delType = 'tab',dataType = 'float')
    day4Rep1 = day4Rep1Inst.run()

    day4Rep2Inst = importData(inputFile = 'expData/day4Rep2.txt',delType = 'tab',dataType = 'float')
    day4Rep2 = day4Rep2Inst.run()

    #--- Compute the measured fold change in growth of the pairs ----
    print 'Computing the fold change in growth ...\n'
    # First find the average over replicates
    # Note that the data matrixes should actually be symmetric as (mutant1,mutant2) is the
    # same as (mutant2,mutant1), however, this is not always the case due to the experiiemntal
    # errors. As such we shouls also take average on the values of (mutant1,mutant2) and
    # (mutant2,mutant1). 
    day1Ave = dict([((m1,m2),(day1Rep1[(m1,m2)] + day1Rep1[(m2,m1)] + day1Rep2[(m1,m2)] + day1Rep2[(m2,m1)])/4) for (m1,m2) in games.keys()])
    day4Ave = dict([((m1,m2),(day4Rep1[(m1,m2)] + day4Rep1[(m2,m1)] + day4Rep2[(m1,m2)] + day4Rep2[(m2,m1)])/4) for (m1,m2) in games.keys()])

    # Experimental fold growth for the community. 1.5e7 was taken from the readme tab of the
    # Excel sheet of the Wintermute paper
    expFoldGrowth = dict([(k,day4Ave[k]/(1.5e7)) for k in games.keys()])

    # Number/experiment comparison (model comes first folloed by experimental result)
    modelExpNumbers = dict([(k,0) for k in ['CC_CC','DD_CC','CD_CC','DC_CC','CDCD_CC','CC_DD','DD_DD','CD_DD','DC_DD','CDCD_DD']])

    # This holds the results of the model-experiment comparisons (i.e., the name of the games is a list)
    modelExpComp = dict([(k,[]) for k in ['CC_CC','DD_CC','CD_CC','DC_CC','CDCD_CC','CC_DD','DD_DD','CD_DD','DC_DD','CDCD_DD']])

    # Simulation time
    sim_time = 96

    #--- Compute the Nash equilibria ---
    list50 = [k for k in games.keys() if expFoldGrowth[k] >= 50]
    for m1m2 in games.keys():
    #for m1m2 in list50:
    #for m1m2 in [('lysA','ilvE')]:
        print m1m2

        game = games[m1m2]
        m1, m2 = m1m2[0],m1m2[1]

        # Create the keys for strategies_prob in stochastic strtategy class 
        prob_keys = game.payoffMatrix.keys()
    
        # Initialize strategies_prob with z zero probability for all elements
        strategies_prob_m1 = dict([(k1,dict([(k2,0) for k2 in game.players_strategies[m1]])) for k1 in prob_keys])
        strategies_prob_m2 = dict([(k1,dict([(k2,0) for k2 in game.players_strategies[m2]])) for k1 in prob_keys])
    
        # Define iStrategies
        iStrategies = {}
    
        iStrategies[(m1,'ALLD')] = stochasticStrategy(player_name = m1, strategy_name = 'ALLD',strategies_prob = strategies_prob_m1)
        iStrategies[(m1,'ALLC')] = stochasticStrategy(player_name = m1, strategy_name = 'ALLC',strategies_prob = strategies_prob_m1)
        iStrategies[(m1,'TFT')] = stochasticStrategy(player_name = m1, strategy_name = 'TFT',strategies_prob = strategies_prob_m1)
    
        iStrategies[(m2,'ALLD')] = stochasticStrategy(player_name = m2, strategy_name = 'ALLD',strategies_prob = strategies_prob_m2)
        iStrategies[(m2,'ALLC')] = stochasticStrategy(player_name = m2, strategy_name = 'ALLC',strategies_prob = strategies_prob_m2)
        iStrategies[(m2,'TFT')] = stochasticStrategy(player_name = m2, strategy_name = 'TFT',strategies_prob = strategies_prob_m2)
    
        # Initial concentrations
        strat_num = 2
        C_init = {(m1,'ALLD'):(1/strat_num)*(1e7),(m1,'ALLC'):(0/strat_num)*(1e7),(m1,'TFT'):(1/strat_num)*(1e7),(m2,'ALLD'):(1/strat_num)*(1e7),(m2,'ALLC'):(0/strat_num)*(1e7),(m2,'TFT'):(1/strat_num)*(1e7)}
    
        # Total number of initial cells
        C_init_tot = sum([C_init[k] for k in C_init.keys()])
    
        # Initial fractions
        x_init = {(m1,'ALLD'):C_init[(m1,'ALLD')]/C_init_tot,(m1,'ALLC'):C_init[(m1,'ALLC')]/C_init_tot,(m1,'TFT'):C_init[(m1,'TFT')]/C_init_tot,(m2,'ALLD'):C_init[(m2,'ALLD')]/C_init_tot,(m2,'ALLC'):C_init[(m2,'ALLC')]/C_init_tot,(m2,'TFT'):C_init[(m2,'TFT')]/C_init_tot}
    
        if abs(sum([x_init[k] for k in x_init.keys()]) - 1) >= 1e-5:
            raise customError('sum of the initial fractions does not add to one\n')
    
        ig = iteratedGame(game = game, iStrategies = iStrategies, x_init = x_init , sim_time = 96, dt = 1)
        x = ig.run()
    
    
        # Outcome of the simulation
        print '      (',m1,'ALLD)  --> ',x[(m1,'ALLD')][sim_time]
        print '      (',m1,'ALLC)  --> ',x[(m1,'ALLC')][sim_time]
        print '      (',m1,'TFT)  --> ',x[(m1,'TFT')][sim_time]
        print '      (',m2,'ALLD)  --> ',x[(m2,'ALLD')][sim_time]
        print '      (',m2,'ALLC)  --> ',x[(m2,'ALLC')][sim_time]
        print '      (',m2,'TFT)  --> ',x[(m2,'TFT')][sim_time]


        sim_outcome = ''
        if x[(m1,'ALLD')][sim_time] <= 0.05 and x[(m2,'ALLD')][sim_time] <= 0.05:
            sim_outcome = 'CC'  # TFT wins
        elif x[(m1,'TFT')][sim_time] <= 0.05 and x[(m2,'TFT')][sim_time] <= 0.05:
            sim_outcome = 'DD'  # Defectros win 
        elif x[(m1,'ALLD')][sim_time] <= 0.05 and x[(m2,'TFT')][sim_time] <= 0.05:
            sim_outcome = 'CD'  # m1 TFTs and m2 defectors win
        elif x[(m1,'TFT')][sim_time] <= 0.05 and x[(m2,'ALLD')][sim_time] <= 0.05:
            sim_outcome = 'DC'  # m1 defectors and m2 TFTs win
        else:
            sim_outcome = 'CDCD'  # All can live together
        
        print '      sim_outcome = ',sim_outcome,'\n'

        # Match with experimental data. The experimental data shows that the self-same pairs cannot
        # grow more than eight fold after four days, while 17% of the dissimilar pairs grow
        # more than 50 fold. We used 50 fold as the threshold to claim a match or dismatch
        # Both only cooperate in a unique Nash equilibrium based on the model
        if expFoldGrowth[m1m2] >= 50:    # Cooperate based on experimental data
            print '      **50+',
            if sim_outcome == 'CC':
                modelExpNumbers['CC_CC'] += 1
                modelExpComp['CC_CC'] += m1m2
            elif sim_outcome == 'DD':
                modelExpNumbers['DD_CC'] += 1
                modelExpComp['DD_CC'] += m1m2
            elif sim_outcome == 'CD':
                modelExpNumbers['CD_CC'] += 1
                modelExpComp['CD_CC'] += m1m2
            elif sim_outcome == 'DC':
                modelExpNumbers['DC_CC'] += 1
                modelExpComp['DC_CC'] += m1m2
            elif sim_outcome == 'CDCD':
                modelExpNumbers['CDCD_CC'] += 1
                modelExpComp['CDCD_CC'] += m1m2
        else:
            if sim_outcome == 'CC':
                modelExpNumbers['CC_DD'] += 1
                modelExpComp['CC_DD'] += m1m2
            elif sim_outcome == 'DD':
                modelExpNumbers['DD_DD'] += 1
                modelExpComp['DD_DD'] += m1m2
            elif sim_outcome == 'CD':
                modelExpNumbers['CD_DD'] += 1
                modelExpComp['CD_DD'] += m1m2
            elif sim_outcome == 'DC':
                modelExpNumbers['DC_DD'] += 1
                modelExpComp['DC_DD'] += m1m2
            elif sim_outcome == 'CDCD':
                modelExpNumbers['CDCD_DD'] += 1
                modelExpComp['CDCD_DD'] += m1m2


    print '\nStatistics (the first elmeent is model and the second is experiment = '
    print 'CC_CC = ',modelExpNumbers['CC_CC']
    print 'DD_DD = ',modelExpNumbers['DD_DD'],'\n'
    for k in modelExpNumbers.keys():
        if k not in ['CC_CC','DD_DD']:
            print k,'        ',modelExpNumbers[k]

    # Save the results
    print 'Save the final results into a file ...\n'
    with open('results/iteratedResults.pk','wb') as outputFile:
         pk.dump(modelExpNumbers,outputFile,-1)

    print '\nResults were written to iteratedResults.pk'

    elapsedT = (time.clock() - startT)/60
    print '\n**Job ended at ',time.strftime("%Y-%m-%d %H:%M:%S"),' (it took a total of ',elapsedT,' min to complete)\n'

