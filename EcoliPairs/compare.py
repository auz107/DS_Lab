from __future__ import division
import re, sys, math, copy, time, random
import numpy as np
sys.path.append('../')
import os.path
import cPickle as pk
from copy import deepcopy
from models import *
from tools.FBA import *
from tools.GAMETES.game import *
from tools.GAMETES.NashEqFinder import *
from tools.customError import *
from tools.importData import *

# ------------------------------------
# Compares the experimentally observed co-growth with the identifed Nash equilibria of
# the game
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
    # Excel sheet of the Wintermute paper.
    expFoldGrowth = dict([(k,day4Ave[k]/(1.5e7)) for k in games.keys()])

    print 'The total of cooperative pairs according to experimental data = %i\n'%(len([m1m2 for m1m2 in expFoldGrowth.keys() if expFoldGrowth[m1m2] >= 8]))

    # Number of (mis)matches between the game theory results and experimental data
    modelExpNumbers = dict([(k,0) for k in ['CC_noEq','DD_noEq','CC_CC','DD_CC','CC_DD','DD_DD','CC_CCandDD','CC_CCandD','CC_DDandC','CC_CandD']])

    #--- Compute the Nash equilibria ---
    print 'Computing the Nash equilibria ...\n'
    for m1m2 in games.keys():
        print m1m2
        NashEqFinderInst = NashEqFinder(games[m1m2], screenOutput = 'off')
        [pureNashEq,terminationFlag] = NashEqFinderInst.runPure()
        games[m1m2].pureNashEq = pureNashEq
        games[m1m2].expFoldGrowth = expFoldGrowth[m1m2]

        print '\npureNashEq = ',pureNashEq 

        #-- Determine what the mutants do in the Nash equilibria --
        # Strategies, where at least one defects
        atLeastOneDefect = [k1 for k1 in pureNashEq if True in [k2.lower() == 'defect' for k2 in dict(k1).values()]]

        # Strategy where both defect
        bothDefect = [k1 for k1 in pureNashEq if all([k2.lower() == 'defect' for k2 in dict(k1).values()])]

        # Strategy where both cooperate
        bothCoop = [k1 for k1 in pureNashEq if True not in [k2.lower() == 'defect' for k2 in dict(k1).values()]]

        # Compute for what percentage of the pairs the game theory predictions are consistent 
        # with experimental data. The experimental data shows that the self-same pairs cannot
        # grow more than eight fold after four days, while 17% of the dissimilar pairs grow
        # more than 50 fold. We used 50 fold as the threshold to claim a match or dismatch
        # Both only cooperate in a unique Nash equilibrium based on the model
        if expFoldGrowth[m1m2] >= 8:    # Cooperate based on experimental data
            # If there are no Nash equilibria based on the model 
            if len(pureNashEq) == 0:
                modelExpNumbers['CC_noEq'] += 1

            # BOth cooperate based on model
            elif len(bothCoop) >= 1 and len(atLeastOneDefect) == 0 and len(bothDefect) == 0: 
                modelExpNumbers['CC_CC'] += 1

            # Both defect based on the model 
            elif len(bothCoop) == 0 and len(atLeastOneDefect) == 1 and len(bothDefect) == 1: 
                modelExpNumbers['CC_DD'] += 1

            # Both cooperate or both defect based on the model 
            elif len(bothCoop) == 1 and len(atLeastOneDefect) == 1 and len(bothDefect) == 1: 
                modelExpNumbers['CC_CCandDD'] += 1

            # Some cooperate but some defect based on the model
            elif len(bothCoop) >= 1 and len(atLeastOneDefect) > 1: 
                modelExpNumbers['CC_CCandD'] += 1

            # Not both cooeprate, but one may cooperate based on the model
            elif len(bothCoop) == 0 and len(atLeastOneDefect) >= 1 and len(bothDefect) == 1: 
                modelExpNumbers['CC_DDandC'] += 1

            # Not both defect  but one may defect based on the model
            elif len(bothCoop) == 0 and (len(atLeastOneDefect) >= 1 or len(bothDefect) == 0): 
                modelExpNumbers['CC_CandD'] += 1

            else:
                print '   **Unnkonwn case1'
        else:                            # Do not cooperate based on experimental data
            # If there are no Nash equilibria based on the model 
            if len(pureNashEq) == 0:
                modelExpNumbers['DD_noEq'] += 1

            # BOth cooperate based on model
            elif len(bothCoop) >= 1 and len(atLeastOneDefect) == 0 and len(bothDefect) == 0: 
                modelExpNumbers['RD_CC'] += 1

            # Both defect based on the model 
            elif len(bothCoop) == 0 and len(atLeastOneDefect) == 1 and len(bothDefect) == 1: 
                modelExpNumbers['DD_DD'] += 1

            # Both cooperate or both defect based on the model 
            elif len(bothCoop) == 1 and len(atLeastOneDefect) == 1 and len(bothDefect) == 1: 
                modelExpNumbers['DD_CCandDD'] += 1

            # Some cooperate but some defect based on the model
            elif len(bothCoop) >= 1 and len(atLeastOneDefect) > 1: 
                modelExpNumbers['DD_CCandD'] += 1

            # Not both cooeprate, but one may cooperate based on the model
            elif len(bothCoop) == 0 and len(atLeastOneDefect) >= 1 and len(bothDefect) == 1: 
                modelExpNumbers['DD_DDandC'] += 1

            # Not both defect  but one may defect based on the model
            elif len(bothCoop) == 0 and (len(atLeastOneDefect) >= 1 or len(bothDefect) == 0): 
                modelExpNumbers['DD_CandD'] += 1

            else:
                print '   **Unnkonwn case1'

    print '\nStatistics = ',modelExpNumbers
    # Save the results
    print 'Save the final results into a file ...\n'
    with open('results/comparisonResults.pk','wb') as outputFile:
         pk.dump([games,modelExpNumbers],outputFile,-1)

    print '\nResults were written to finalResults.pk'
    elapsedT = (time.clock() - startT)/60
    print '\n**Job ended at ',time.strftime("%Y-%m-%d %H:%M:%S"),' (it took a total of ',elapsedT,' min to complete)\n'
