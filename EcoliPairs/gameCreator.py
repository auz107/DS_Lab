from __future__ import division
import re, sys, math, copy, time, random
import numpy as np
import sys
sys.path.append('../')
import cPickle as pk
from copy import deepcopy
from tools.showTime import *
from models import *
from tools.FBA import *
from metabolicModel import *
from fba import *
from auxoMetabsFinder import *
from tools.GAMETES.game import *
from tools.customError import *

class gameCreator(object):
     """
     Creates the instances of the class game for different pairs of auxotrophic E. coli mutants
     """

     def __init__(self, mutants_rxn_info, auxoMetabsMutants, wildTypeMetabolicModel,infiniteCostMetabs = None, optSolverName = None, loopStart = None,loopEnd = None, screenOutput = None, outputFile = None):
        """
        INPUTS (required):
        ------
              mutants_rxn_info: A dictionary whose keys are the names of the mutants and values
                             are the list of reactions that should be removed in the presence
                             of that mutant
          auxoMetabsMutants: A dictionary whose keys are the mutant names and values are the
                             list of exchange rxns for metabolites it needs to survive
     wildTypeMetabolicModel: An instance of the class metabolicModel containing the metabolic
                             model of the wild-type strain
          infiniteCostMetab: List of the exchange rxns corresponding to metabolites that 
                             cannot be produced by wild-type, i.e., their cost of production
                             is infinity
              optSolverName: Name of the optimization solver to be used to solve the MILP. 
                             Current allowable choices are cplex and gurobi
                  loopStart: Start index of the mutants list for the loop (default = 0)
                    loopEnd: End index of the mutants list for the loop (defaulty = (# of mutants) - 1)
               screenOutput: By default (None) writes  a summary including the 
                             solve  status, optimality status (if not optimal),
                             objective function value and the elapsed time on 
                             the screen if takes a value of 'silent' no resuults 
                             are written on the screen, in which case The user 
                             can instead specifiy  an output file using the 
                             option outputFile, or store them in the variable 
                             runOutput (see the 'run' method for details)
                 outputFile: Optional input. It is a string containg the path to a 
                             file and its name (e.g., 'results/fbaResults.txt'), 
                             where the results should be written to. 
  
        """
        # mutantINfo
        self.mutants_rxn_info = deepcopy(mutants_rxn_info)

        # List of metabolites that can rescue a mutant 
        self.auxoMetabsMutants = deepcopy(auxoMetabsMutants)

        # Metabolic model of the wild-type strain
        self.wildTypeMetabolicModel = deepcopy(wildTypeMetabolicModel)

        # List of the exchange rxns corresponding to metabolites that cannot be produced
        if infiniteCostMetabs == None:
            self.infiniteCostMetabs = infiniteCostMetabs
        else:
            self.infiniteCostMetabs = infiniteCostMetabs[:]

        # Start index of the mutants list for the loop (default = 0)
        if loopStart == None:
            self.loopStart = 0
        else:
            self.loopStart = loopStart

        # End index of the mutants list for the loop 
        if loopEnd == None:
            self.loopEnd = len(mutants_rxn_info.keys()) - 1
        else:
            self.loopEnd = loopEnd

        # Solver name
        if optSolverName == None:
            self.optSolverName = 'cplex'
        else:
            if optSolverName.lower() in ['cplex','gurobi']:
                self.optSolverName = optSolverName
            else:
                raise customError('**Error! Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Output to the screen 
        if screenOutput != None and screenOutput.lower() != 'silent':
            raise customError("**Error! The only eligible value for screenOutput is 'silent'")
        else:
             self.screenOutput = screenOutput

        # Details of the output to the screen (None or 'silent')
        self.screenOutputDetails = None

        # Output file
        self.outputFile = outputFile

        # Choose a value to show a big cost in the case a straiin cannot product a metabolite
        self._infCost = 10

     def _isInfiniteCostMetab(self,strategy):
        """
        Checks if any of the exchange rxns appearing in a strategy is in infiniteCostMetab 
        """
        output = False
        done = 0
        counter = 0
        while done == 0 and counter <= len(strategy)-1:
            if strategy[counter] in self.infiniteCostMetabs: 
                done = 1
                output = True

            counter += 1
        return output 

     def createPayoff(self):
        """
        Creates the payoff matrix for a mutant pair
        """
 
        # Examine the combination of all strategies
        for strategy_m1 in self._players_strategies[self._mutant1]:
            for strategy_m2 in self._players_strategies[self._mutant2]:

                #--- If both defect ---
                if strategy_m1.lower() == 'defect' and strategy_m2.lower() == 'defect':
                    payoff_mutant1 = 0
                    payoff_mutant2 = 0
  
                #--- If only one defects ---
                # For the payoff of the cooperating strain compute the cost of producing
                # the metabolites for its partnet 'assumiing' that its partnet also 
                # cooperates and then assign the negative of this cost as its net payoff
                elif strategy_m1.lower() == 'defect' or strategy_m2.lower() == 'defect':
                    if strategy_m1.lower() == 'defect' and strategy_m2.lower() != 'defect':
                        defector = self._mutant1
                        cooperator = self._mutant2
                        strategyDetails_cooperator = self._players_strategiesDetails[self._mutant2][strategy_m2]
                    elif strategy_m1.lower() != 'defect' and strategy_m2.lower() == 'defect':
                        defector = self._mutant2
                        cooperator = self._mutant1
                        strategyDetails_cooperator = self._players_strategiesDetails[self._mutant1][strategy_m1]
                    else:
                        raise customError('**Error! Both partnets either defect or cooperate. Check out the code for error.') 

                    #-- Defector --
                    fba_defector = deepcopy(self.fba_wildType)
                    fba_defector.createModel = 0
                    for rxn in self.mutants_rxn_info[defector]:
                        fba_defector.fbaModel.v[rxn] = 0                    
                        fba_defector.fbaModel.v[rxn].fixed = True

                    # Take up one unit of the metabolites produced by cooperator 
                    for auxoMetabExch in strategyDetails_cooperator:
                        fba_defector.fbaModel.v[auxoMetabExch].setlb(-1)

                    [fbaExitFlag,biomassFlux] = fba_defector.run()
                    if fbaExitFlag == 'globallyOptimal':
                        payoff_defector = biomassFlux
                    else:
                        payoff_defector = None
                        if fbaExitFlag == 'solverError' and self.screenOutput == None:
                            print '   ',fbaExitFlag

                        if fbaExitFlag != 'solverError' and self.screenOutputDetails == None:
                            print '   ',fbaExitFlag

                    #-- Cooperator --
                    fba_cooperator = deepcopy(self.fba_wildType)
                    fba_cooperator.createModel = 0 

                    # Perform FBA only if infiniteCostMetabs is not provided or otherwise if 
                    # at least one exchange reaction in the current strategy is in 
                    # infiniteCostMetabs 
                    if self.infiniteCostMetabs == None or (self.infiniteCostMetabs != None and self._isInfiniteCostMetab(strategyDetails_cooperator) == False):

                        # Since this mutant cooperates but its partnet defects compute the cost 
                        # of cooperations for the wild-type and assign the negative of this cost
                        # as the payoff
                        # Cooperator produces the compounds needed by defector 
                        for auxoMetabExch in strategyDetails_cooperator:
                            fba_cooperator.fbaModel.v[auxoMetabExch].setlb(1)

                        [fbaExitFlag,biomassFlux] = fba_cooperator.run()
                        if fbaExitFlag == 'globallyOptimal':
                            payoff_cooperator = -(self.biomassFlux_wildType - biomassFlux)

                        # If the model is infeasible assign the negative of the bigCost as the  
                        # payoff. This means that the cost of producing that metabolite is 
                        # infinity
                        else:
                            payoff_cooperator = -self._infCost
                            if fbaExitFlag == 'solverError' and self.screenOutput == None:
                                print '   ',fbaExitFlag

                            if fbaExitFlag != 'solverError' and self.screenOutputDetails == None:
                                print '   ',fbaExitFlag
                    else: 
                        payoff_cooperator = -self._infCost

                    if strategy_m1.lower() == 'defect' and strategy_m2.lower() != 'defect':
                        payoff_mutant1 = payoff_defector
                        payoff_mutant2 = payoff_cooperator

                    elif strategy_m1.lower() != 'defect' and strategy_m2.lower() == 'defect':
                        payoff_mutant1 = payoff_cooperator
                        payoff_mutant2 = payoff_defector
                    else:
                        raise customError('**Error! Both partnets either defect or cooperate. Check out the code for error.') 

                #--- if both cooperate ---                    
                else:
                    strategyDetails_mutant1 = self._players_strategiesDetails[self._mutant1][strategy_m1]
                    strategyDetails_mutant2 = self._players_strategiesDetails[self._mutant2][strategy_m2] 

                    #-- Mutant 1 --
                    # Perform FBA only if infiniteCostMetabs is not provided or otherwise if 
                    # at least one exchange reaction in the current strategy is in 
                    # infiniteCostMetabs 
                    if self.infiniteCostMetabs == None or (self.infiniteCostMetabs != None and self._isInfiniteCostMetab(strategyDetails_mutant1) == False):

                        fba_mutant1 = deepcopy(self.fba_wildType)
                        fba_mutant1.createModel = 0 
                        for rxn in self.mutants_rxn_info[self._mutant1]:
                            fba_mutant1.fbaModel.v[rxn] = 0                    
                            fba_mutant1.fbaModel.v[rxn].fixed = True

                        # Take up one unit of the metabolites produced by mutant 2 
                        for auxoMetabExch in strategyDetails_mutant2: 
                            fba_mutant1.fbaModel.v[auxoMetabExch].setlb(-1)


                        # Produce one unit of the metabolites needed by mutant 2 to survive
                        for auxoMetabExch in strategyDetails_mutant1: 
                            fba_mutant1.fbaModel.v[auxoMetabExch].setlb(1)

                        [fbaExitFlag,biomassFlux] = fba_mutant1.run()
                        if fbaExitFlag == 'globallyOptimal':
                            payoff_mutant1 = biomassFlux

                        # If the model is infeasible it means that it cannot produce that 
                        # metabolite or the cost of producing that metabolite is infinity. 
                        # Instead of infinity assign the negative of bigCost as the payoff
                        else:
                            payoff_mutant1 = -self._infCost
                            if fbaExitFlag == 'solverError' and self.screenOutput == None:
                                print '   ',fbaExitFlag
       
                            if fbaExitFlag != 'solverError' and self.screenOutputDetails == None:
                                print '   ',fbaExitFlag

                    else: 
                        payoff_mutant1 = -self._infCost

                    #-- Mutant 2 --
                    # Perform FBA only if infiniteCostMetabs is not provided or otherwise if 
                    # at least one exchange reaction in the current strategy is in 
                    # infiniteCostMetabs 
                    if self.infiniteCostMetabs == None or (self.infiniteCostMetabs != None and self._isInfiniteCostMetab(strategyDetails_mutant2) == False):

                        fba_mutant2 = deepcopy(self.fba_wildType)
                        fba_mutant2.createModel = 0 
                        for rxn in self.mutants_rxn_info[self._mutant2]:
                            fba_mutant2.fbaModel.v[rxn] = 0                    
                            fba_mutant2.fbaModel.v[rxn].fixed = True

                        # Take up one unit of the metabolites produced by mutant 1 
                        for auxoMetabExch in strategyDetails_mutant1: 
                            fba_mutant2.fbaModel.v[auxoMetabExch].setlb(-1)

                        # Produce one unit of the metabolites needed by mutant 1 to survive
                        for auxoMetabExch in strategyDetails_mutant2: 
                            fba_mutant2.fbaModel.v[auxoMetabExch].setlb(1)

                        [fbaExitFlag,biomassFlux] = fba_mutant2.run()
                        if fbaExitFlag == 'globallyOptimal':
                            payoff_mutant2 = biomassFlux

                        # see the comment for mutant 1
                        else:
                            payoff_mutant2 = -self._infCost

                            if fbaExitFlag == 'solverError' and self.screenOutput == None:
                                print '   ',fbaExitFlag

                            if fbaExitFlag != 'solverError' and self.screenOutputDetails == None:
                                print '   ',fbaExitFlag

                    else: 
                        payoff_mutant2 = -self._infCost

                # Assign the payoffs
                if self.screenOutputDetails == None:
                    print '    ',((self._mutant1,strategy_m1),(self._mutant2,strategy_m2)),' = ',{self._mutant1:payoff_mutant1,self._mutant2:payoff_mutant2}

                self.payoffMatrix[((self._mutant1,strategy_m1),(self._mutant2,strategy_m2))] = {self._mutant1:payoff_mutant1,self._mutant2:payoff_mutant2} 

        if self.screenOutputDetails == None:
            print '\n'

     def run(self):
        """
        OUTPUTS:
        -------
        games:  A dictionary whose keys are the tuples showing the mutant pairs and values 
                are instances of the class game 
        """

        # Create a general FBA model using the model for wild-type
        self.fba_wildType = fba(metabolicModel = self.wildTypeMetabolicModel, optSolverName = self.optSolverName) 
        [fbaExitFlag,biomassFlux] = self.fba_wildType.run()
        if fbaExitFlag == 'globallyOptimal':
            self.biomassFlux_wildType = biomassFlux
        else:
            self.biomassFlux_wildType = 0

        self.fba_wildType.screenOutput = 'silent'

        # List of mutants (do not consider mutants that cannot be rescued even if all exchange
        # rxns are allowed to have an upper bound of less than zero)
        mutantsList = [mutant for mutant in self.auxoMetabsMutants.keys() if len(self.auxoMetabsMutants[mutant]) >= 1] 

        games = {}

        for self._mutant1 in mutantsList[self.loopStart:self.loopEnd+1]:
            mutant1_index = mutantsList.index(self._mutant1)
            for self._mutant2 in mutantsList[mutant1_index+1:]:
                    
                if self.screenOutput == None:
                    print '(',self._mutant1,',',self._mutant2,')'

                game_name = self._mutant1 + '_' + self._mutant2
                players_names = [self._mutant1,self._mutant2]  

                # Players' strategies
                self._players_strategies = {}
 
                # Players' strategies details
                self._players_strategiesDetails = dict((player_name,{}) for player_name in players_names)

                # Mutant 1: Possible strategies are to either defect or produce the 
                # metabolites needed by its partner mutant 2
                self._players_strategies[self._mutant1] = [] 
                self._players_strategies[self._mutant1].append('Defect')
                self._players_strategiesDetails[self._mutant1]['Defect'] = ['Defect']
                for strategy in self.auxoMetabsMutants[self._mutant2]:
                    strategyName = ' and '.join(strategy)
                    self._players_strategies[self._mutant1].append(strategyName)
                    self._players_strategiesDetails[self._mutant1][strategyName] = strategy 

                # Mutant 2: Possible strategies are either to defect or to produce the 
                # metabolites needed by its partnet mutant 1
                self._players_strategies[self._mutant2] = [] 
                self._players_strategies[self._mutant2].append('Defect')
                self._players_strategiesDetails[self._mutant2]['Defect'] = ['Defect']
                for strategy in self.auxoMetabsMutants[self._mutant1]:
                    strategyName = ' and '.join(strategy)
                    self._players_strategies[self._mutant2].append(strategyName)
                    self._players_strategiesDetails[self._mutant2][strategyName] = strategy 

                if self.screenOutputDetails == None:
                    print '   Strategies of',self._mutant1,' are: ',self._players_strategies[self._mutant1]
                    print '   Strategies of',self._mutant2,' are: ',self._players_strategies[self._mutant2]

                # Payoff matrix      
                self.payoffMatrix = {}
                self.createPayoff()

                # Create the game
                games[(self._mutant1,self._mutant2)] = game(game_name = game_name, players_names = players_names, players_strategies = self._players_strategies, players_strategiesDetails = self._players_strategiesDetails, payoffMatrix = self.payoffMatrix)

                # Store the results into a binary file
                with open(self.outputFile,'wb') as pkOutputFile:
                    pk.dump(games,pkOutputFile,-1)

        return games

#-------------------- Sample implementation -----------------
#if __name__ == "__main__":

def mainRun(loopStart = None,loopEnd = None):

    from mutants_rxn_info import * 
    from gameCreator import *

    startT = time.clock()
    print '\n**Job started at ',time.strftime("%Y-%m-%d %H:%M:%S"),'\n'

    #---  E. coli iAF1260 model ----
    # Import general model files
    from iAF1260ModelData import *

    # Create a general model 
    print 'Creating the metabolic model ...\n'
    iAF1260Model = metabolicModel(microorganism_name = 'E.coli', model_name = 'iAF1260', metab_names = metab_names,rxn_names = rxn_names, stoic_matrix = stoic_matrix, rxn_types = rxn_types, biomass_rxn_name = biomass_rxn_name)

    # Import the medium files (anaerobic minimal glucose medium)
    from iAF1260_minimalGlucose_anaerobic import *
    for rxn in mediumSpecific_fluxBounds.keys():
        iAF1260Model.flux_bounds[rxn] = mediumSpecific_fluxBounds[rxn]

    # Objective vector
    obj_vec = dict((rxn,0) for rxn in iAF1260Model.rxn_names)
    obj_vec[biomass_rxn_name] = 1
    iAF1260Model.obj_vec = obj_vec

    optSolverName = 'gurobi'

    #--- Load the list of metabolites each mutant needs to survive ---
    with open('results/auxoMetabs.pk','rb') as inputFile:
        termCond_mutant,auxoMetabsMutants = pk.load(inputFile)             

    #--- infiniteCostMetabs ---
    # Load the list of exchange rxns corresponding to metabolites that cannot be produced
    # by the wildtype
    with open('results/infiniteCostMetabs.pk','rb') as inputFile:
        infiniteCostMetabs = pk.load(inputFile)

    # Name of the binary file name to store the results 
    if loopStart != None and loopEnd != None:
        outputFile = 'results/auxoEcoliGames_' + str(loopStart) + '_' + str(loopEnd) +'.pk'
    elif loopStart != None and loopEnd == None:
        outputFile = 'results/auxoEcoliGames_' + str(loopStart) + '_end.pk'
    elif loopStart == None and loopEnd != None:
        outputFile = 'results/auxoEcoliGames_start_' + str(loopEnd) +'.pk'
    else:
        outputFile = 'results/auxoEcoliGames_all.pk'

    # The following is used to test the code
    auxoMetabsMutants2 = {} 
    auxoMetabsMutants2['lysA'] = auxoMetabsMutants['lysA']
    auxoMetabsMutants2['ilvE'] = auxoMetabsMutants['ilvE']

    #--- Create the game ---
    EcoliGames = gameCreator(mutants_rxn_info = mutants_rxn_info, auxoMetabsMutants = auxoMetabsMutants, wildTypeMetabolicModel = iAF1260Model, infiniteCostMetabs = infiniteCostMetabs, optSolverName = optSolverName, loopStart = loopStart, loopEnd = loopEnd, outputFile = outputFile)
    EcoliGames.screenOutputDetails = 'silent'
    games = EcoliGames.run()

    print '\nThe results were written into ',outputFile
    elapsedT = (time.clock() - startT)/60
    print '\n**Job ended at ',time.strftime("%Y-%m-%d %H:%M:%S"),' (it took a total of ',elapsedT,' min to complete)\n'


