from __future__ import division
import re, sys, math, copy, time, random
from scipy.linalg import eig
from copy import deepcopy
sys.path.append('../../')
from tools.userError import *
from stochasticStrategy import *
from game import *

class iteratedGame(object):
    """
    This is a general class simulating an iterated game with N players 
    Replicator dynamics is used to simualte the evolution of the community. 

    Ali R. Zomorrodi - Daniel Segre lab @ BU
    Last updated: 06-11-2015
    """

    def __init__(self, game, iStrategies, x_init, sim_time, dt = None):
        """
        INPUTS:     
        -------
                   game: An instance of the class game 
            iStrategies: A dictionary holding the information about strategies 
                         that should be played in the iterated games. Keys and 
                         values of this dictionary are as follows:
                           Keys: A tuple with two elements: 
                                 1st element: Name of the player playing this strategy
                                 2nd element: is the Name of the strategy (e.g., ('lysA','ALLD'))
                         Values: An instance of the class stochasticStrategy
                 x_init: A dictionary containing the initial composition of the community 
                           Keys: The same keys as iStrategies
                         Values: The fraction of each community member playing that strategy
               sim_time: Simulation time for  replicator dynamics. The units are arbitrary
                     dt: Time step for simulations using replicator dynamics (default = 1) 
        """
        # Game
        self.game = game

        # Payoff matrix of the original game
        self._payoffMatrix = self.game.payoffMatrix

        # Total possible states of the game
        self._gameStates = self._payoffMatrix.keys() 

        # Iterated strategies
        self.iStrategies = deepcopy(iStrategies)

        # Initial concentrations
        self.x_init = x_init

        # Simulation time
        self.sim_time = sim_time

        # Time step for discretized replicator dynamics
        if dt == None:
            self.dt = 1
        else:         
            self.dt = dt

    def createPayoff_iStrategies(self):
        """
        This function creates the payoff matrix of the iterated strategies.
        The output of this function is a dictionary called payoff_iStrategies
        with keys and values as follows:
          keys: A tuple with two elements, where each element is a key of
                self.iStrategies, which is a tuple with the first element being
                the name of the player and the second the name of the iterated
                strategy
        values: Another dictionary with keys and values as follows: 
                  keys: Name of the player 
                values: The payoff of the player
        """
        self.payoff_iStrategies = {}

        for iStrategy1 in self.iStrategies.keys():
            for iStrategy2 in self.iStrategies.keys():
                # Players of the same type get nothing when facing with each other 
                if iStrategy1[0].lower() == iStrategy2[0].lower():
                    self.payoff_iStrategies[(iStrategy1,iStrategy2)] = {iStrategy1[0]:0, iStrategy2[0]:0}

                elif self.iStrategies.keys().index(iStrategy2) > self.iStrategies.keys().index(iStrategy1):
                    self._iStrategies_curr = {}
                    # If we have TFT vs TFT it's similar to ALLC vs. ALLC
                    if iStrategy1[1].lower() == 'tft' and iStrategy2[1].lower() == 'tft':
                        print 'iStrategy1[1] = ',iStrategy1[1],'\tiStrategy 2 = ',iStrategy2[1] 
                        iStrategy1[1] == 'ALLC'
                        iStrategy2[1] == 'ALLC'
                    # If we have sTFT vs sTFT it's similar to ALLD vs. ALLD
                    elif iStrategy1[1].lower() == 'stft' and iStrategy2[1].lower() == 'stft':
                        print 'iStrategy1[1] = ',iStrategy1[1],'\tiStrategy 2 = ',iStrategy2[1] 
                        iStrategy1[1] == 'ALLD'
                        iStrategy2[1] == 'ALLD'
 
                    self._iStrategies_curr[iStrategy1] = self.iStrategies[iStrategy1]
                    self._iStrategies_curr[iStrategy2] = self.iStrategies[iStrategy2]

                    # Create the Markov chain transition matrix
                    self.createTransitionMatrix()

                    # Find the stationary probabilities
                    self.findStationaryProbs()

                    # Compute the expected payoff for each player
                    # Player playing strategies iStrategy1 and iStrategy2
                    p1,p2 = iStrategy1[0], iStrategy2[0]
                    p1Payoff, p2Payoff = [], []
                    for stationaryProbs in self._stationaryProbs:
                        p1Payoff.append(sum([stationaryProbs[k]*self._payoffMatrix[k][p1] for k in stationaryProbs.keys()]))
                        p2Payoff.append(sum([stationaryProbs[k]*self._payoffMatrix[k][p2] for k in stationaryProbs.keys()]))

                    # If there are more than one set of stationary probabilities, the payoffs are
                    # computed by taking the average over all
                    self.payoff_iStrategies[(iStrategy1,iStrategy2)] = {p1:sum(p1Payoff)/len(p1Payoff),p2:sum(p2Payoff)/len(p2Payoff)} 

                # Otherwise this strategy pair has already been considered
                elif self.iStrategies.keys().index(iStrategy2) < self.iStrategies.keys().index(iStrategy1):
                    self.payoff_iStrategies[(iStrategy1,iStrategy2)] = self.payoff_iStrategies[(iStrategy2,iStrategy1)] 

                else:
                    raise userError('**Error! Unknonw strategy pair!\n')
        
    def createTransitionMatrix(self):
        """
        This funciton creates the Markov chain transition matrix for a set of iterated 
        strategies (ALLD, ALLC, TFT, etc), one for each player, stored in 
        self._iStrategies_curr, which is a dictionary created by createPayoff_iStrategies
        with keys and values as follows: 
          Keys: A tuple containing the player name and its strategy 
        Values: An instance of the class stochasticStrategy extracted from 
                self.iStrategies. The number of keys is equal to the number of players.
        """
        # Transition matrix
        self.transitionMatrix = {}

        for gameState_prev in self._gameStates: 
            for gameState_curr in self._gameStates: 
                counter = 0
                for player in self.game.players_names: 
                    # Strategy taken by this player in gameState_curr
                    player_strategy = dict(list(gameState_curr))[player]

                    # Find the name of iterated strategy this player has taken
                    for k in self._iStrategies_curr.keys():
                        if player.lower() == k[0].lower():
                            player_iStrategy = k

                    # The probability of taking this strategy by this player accordiing to the
                    # definition of the stochastic strategies
                    prob = self._iStrategies_curr[player_iStrategy].strategies_prob[gameState_prev][player_strategy]
                    if counter == 0:
                        self.transitionMatrix[(gameState_prev,gameState_curr)] = prob
                    else: 
                        self.transitionMatrix[(gameState_prev,gameState_curr)] = prob*self.transitionMatrix[(gameState_prev,gameState_curr)]

                    counter += 1

        """
        # Initial probability distribution vector
        self.probDistrInit = {}
        # If prob_distr_init has been provided
        if len([iStrategy for iStrategy in self._iStrategies_curr.values() if iStrategy.prob_distr_init == None]) == 0:
            for gameState in self._gameStates: 
                counter = 0
                for player in self.game.players_names: 
                    counter += 1
    
                    # Strategy taken by this player in gameState
                    player_strategy = dict(list(gameState))[player]
    
                    # Find the name of iterated strategy this player has taken
                    for k in self._iStrategies_curr.keys():
                        if player.lower() == k[0].lower():
                            player_iStrategy = k
    
                    # The probability of taking this strategy by this player accordiing to the
                    # definition of prob_distr_init of the stochastic strategies
                    prob = self._iStrategies_curr[player_iStrategy].prob_distr_init[player_strategy]
                    if counter == 1:
                        self.probDistrInit[gameState] = prob
                    else: 
                        self.probDistrInit[gameState] = prob*self.probDistrInit[gameState]
        """

    def findStationaryProbs(self):
        """
        This function finds the stationary probabilities of being in each state of the
        game using the Markov chain transition matrix. p = p*M. The output is a 
        dictionary self._stationaryProbs with keys and values as follows:
          keys: States of the game (elements of self._gameStates) 
        values: The stationay probability of being in each state when a particular
                set of iterated strategies face each other
        """
        
        # Initialize the transition matrix in the form of a numpy array
        M = np.zeros((len(self._gameStates),len(self._gameStates)))

        # Fill in the elements of transition matrix M 
        for row in self._gameStates:
            for col in self._gameStates: 
                M[self._gameStates.index(row),self._gameStates.index(col)] = self.transitionMatrix[(row,col)]

        # Sum of all columns for each row should be equal to one
        sumM = M.sum(axis=1)
        for m in sumM:
            if abs(sumM[m] - 1) >= 1e-6:
                raise userError('Sum of the columns for each row is not one for the transition matrix of  ' + srt(self._iStrategies_curr.keys()) + ': sum = ' + str(sumM))

        # Now find the stationary probabilities (left eigenvectors asociated with an 
        # eigenvalue of one)
        eigVals,eigVecs=eig(M,left=True,right=False)

        # Find the indices of eigen valueis, which is one (closest to one)         
        oneInds = [k for k in range(len(eigVals)) if abs(eigVals[k] - 1) < 1e-6]

        # Eigeven value of one selected for TFT strategies (if applicable)
        tft_oneInds = []

        # There should be only one eigenvalue equal to one, if not issue a warning and use the
        # first index in oneEigValInds
        if len(oneInds) > 1:
            print '**Warning! More than one eigenvalues are close to one for ',self._iStrategies_curr.keys()
            print '\tM = \n',M,'\n'
            print '\teigVals = \n\t',eigVals,'\n'
            print '\teigVecs = '
            for oneInd in oneInds:
                print '\t',eigVecs[:,oneInd]/sum(eigVecs[:,oneInd])

            """
            # Check which one corresponds to the given initial probability distribution
            if self.probDistrInit != {}:
                # Initialize the initial probability distribution of a game in the form of a numpy array
                p_init = np.zeros((1,len(self._gameStates)))

                # Fill in the elements of p_init 
                for row in self._gameStates:
                    p_init[self._gameStates.index(row)] = self.ProbDistrInit[row]

                # Check if the sum is equal to one
                if abs(p_init.sum() - 1) >= 1e-6:
                    raise userError('Sum of the initial probabilites in p_init is not one:' + str(p_init.sum()))

                # Consider a large power of M
                Mn = np.linalg.matrix_power(M,100)

                # FInd p_init*(M^100) 
                p_init_Mn = np.dot(p_init,Mn)
            """

            iStrategies = [s[1].lower() for s in self._iStrategies_curr.keys()]
            # If one strategy is TFT and the other is sTFT choose the eigen value whose eigen vector entries
            # are not just zero and one
            if set(iStrategies) == set(['tft','stft']):
                for oneInd in oneInds:
                    # Normalize the eigenvector corresponding to the eigen value of one
                    # This gives the stationary probabilities 
                    eigVecOne_norm = eigVecs[:,oneInd]/sum(eigVecs[:,oneInd])

                    isZeroOne = False
                    continueLoop = True
                    i = 0
                    while not isZeroOne and continueLoop:
                        # Check if any element is not equal to one
                        if abs(eigVecOne_norm[i]) > 1e-6 and abs(1 - eigVecOne_norm[i]) > 1e-6:
                            isZeroOne = True 
                        i += 1
                        if i == len(eigVecOne_norm) - 1:
                            continueLoop = False 

                    if isZeroOne:
                        tft_oneInds.append(oneInd) 

        if len(tft_oneInds) > 0:
            oneInds = tft_oneInds

        stationaryProbArrays = []
        for oneInd in oneInds:
            # Normalize the eigenvector corresponding to the eigen value of one
            # This gives the stationary probabilities 
            eigVecOne_norm = eigVecs[:,oneInd]/sum(eigVecs[:,oneInd])

            if len(tft_oneInds) > 0:
                print '\tChosen eigen vector = ',eigVecOne_norm

            # Check if any elements are negative
            negElems = [eigVecOne_norm[k] for k in range(max(eigVecOne_norm.shape)) if eigVecOne_norm[k] < 0 and abs(eigVecOne_norm[k]) > 1e-6]
            if len(negElems) >= 1:
                raise userError('Negative elements in Statioanry probability: ' + str(eigVecOne_norm))
            else:
                stationaryProbArrays.append(eigVecOne_norm)

        # Convert the numpy array to list of a python dictionary
        self._stationaryProbs = [] 
        for stationaryProbArray in stationaryProbArrays:
            self._stationaryProbs.append(dict([(self._gameStates[k],stationaryProbArray[k]) for k in range(len(self._gameStates))]))

    def replicatorDynamics(self):
        """
        This function implements the replicatory dynamics using the payoff 
        matrix of the iterated strategies. The output is the frequencies of
        community members at time t + dt. This information is held in 
        a dictionary self.x_tpdt, which is a dictionary with keys and values
        as follows:
          Keys: The same as keys of self.iStrategies
        Values: The corresponding frequency
        """

        # A dictionary whose keys are the same as those of self.iStrategies and values are
        # and values are their fitnesses
        self.fitness = {}

        # Compute the fitness of each player playing a certain strategy
        # Loop over each player that plays a certain strategy
        for player_strategy in self.x.keys():
            player_name = player_strategy[0]
            self.fitness[player_strategy] = sum([self.payoff_iStrategies[k][player_name]*self.x[k[1]][self.t] for k in self.payoff_iStrategies.keys() if k[0] == player_strategy ])             
 
        # Compute the expected fitness of the population (phi(x))
        self.phi = sum([self.fitness[k]*self.x[k][self.t] for k in self.x.keys()])
 
        # Update the abundances
        for k in self.x.keys():
            self.x[k][self.t + self.dt] = max(0,self.dt*self.x[k][self.t]*(self.fitness[k] - self.phi) + self.x[k][self.t])

    def run(self):
        """
        This function runs the entire iterated game. The output is 
        a dictionary self.x with keys and values as follows:
          keys: The same as keys of self.x_t or self.iStrategies (a tuple where the first 
                element is the name of the player and the second element is the name of 
                iterated strategy it plays)
        Values: A dictionary with keys and values as follows:
                  Keys: Time points
                Values: Frequency at that time point
        """

        # Create the payoff matrix of iterated strategies
        self.createPayoff_iStrategies()


        #------ Implement in silico evolution --------

        # Current time 
        self.t = 0

        # Dictionary containing the frequences over time
        self.x = {}

        for player_strategy in self.x_init.keys():
            self.x[player_strategy] = {0:self.x_init[player_strategy]}        

        while self.t <= self.sim_time:
            # Update frequencies using the replicator dynamics  
            self.replicatorDynamics()          

            self.t += self.dt

        # Return self. as the output
        return self.x



#--------------- Sample implementation -----------------
if __name__ == "__main__":

    import cPickle as pk
    import matplotlib
    import matplotlib.pyplot as plt
    from stochasticStrategy import *
    from game import *

    def plotResults(m):
        m1 = m[0]
        m2 = m[1]

        #--- Load the games --- 
        print 'Loading the games ...\n'
        with open('../../EcoliPairs/results/auxoEcoliGames_all_old.pk','rb') as inputFile:
            games = pk.load(inputFile)
    
        # Choose a specific game
        sgame = games[(m1,m2)]
    
        print 'payoffMatrix = \n'
        for k in sgame.payoffMatrix.keys():
            print k,' = ',sgame.payoffMatrix[k]
        print '\n'
    
        # Create the keys for strategies_prob in stochastic strtategy class 
        prob_keys = sgame.payoffMatrix.keys() 
    
        # Initialize strategies_prob with z zero probability for all elements
        strategies_prob_m1 = dict([(k1,dict([(k2,0) for k2 in sgame.players_strategies[m1]])) for k1 in prob_keys])
        strategies_prob_m2 = dict([(k1,dict([(k2,0) for k2 in sgame.players_strategies[m2]])) for k1 in prob_keys])
    
        # Define iStrategies
        iStrategies = {}
    
        iStrategies[(m1,'ALLD')] = stochasticStrategy(player_name = m1, strategy_name = 'ALLD',strategies_prob = strategies_prob_m1)      
        iStrategies[(m1,'TFT')] = stochasticStrategy(player_name = m1, strategy_name = 'TFT',strategies_prob = strategies_prob_m1)      
    
        iStrategies[(m2,'ALLD')] = stochasticStrategy(player_name = m2, strategy_name = 'ALLD',strategies_prob = strategies_prob_m2)      
        iStrategies[(m2,'sTFT')] = stochasticStrategy(player_name = m2, strategy_name = 'sTFT',strategies_prob = strategies_prob_m2)      
    
        # Initial concentrations
        strat_num = 2
        C_init = {(m1,'ALLD'):(1/strat_num)*(1e7),(m1,'TFT'):(1/strat_num)*(1e7),(m2,'ALLD'):(1/strat_num)*(1e7),(m2,'sTFT'):(1/strat_num)*(1e7)}
    
        # Total number of initial cells
        C_init_tot = sum([C_init[k] for k in C_init.keys()])
    
        # Initial fractions
        x_init = {(m1,'ALLD'):C_init[(m1,'ALLD')]/C_init_tot,(m1,'TFT'):C_init[(m1,'TFT')]/C_init_tot,(m2,'ALLD'):C_init[(m2,'ALLD')]/C_init_tot,(m2,'sTFT'):C_init[(m2,'sTFT')]/C_init_tot}
    
        if abs(sum([x_init[k] for k in x_init.keys()]) - 1) >= 1e-5:
            raise userError('sum of the initial fractions does not add to one\n')
        
        ig = iteratedGame(game = sgame, iStrategies = iStrategies, x_init = x_init , sim_time = 96, dt = 1)
        x = ig.run()
    
        sim_time = 96
        print '  (','m1,ALLD)',x[(m1,'ALLD')][sim_time]
        print '  (','m1,TFT)',x[(m1,'TFT')][sim_time]
        print '  (','m2,ALLD)',x[(m2,'ALLD')][sim_time]
        print '  (','m2,sTFT)',x[(m2,'sTFT')][sim_time]
        print '\n'   
 
        fig, ax = plt.subplots()
        plt.rc('xtick', labelsize=32)
        plt.rc('ytick', labelsize=32)
        plotList = [k for k in x.keys() if k[1] != 'ALLC']
        #for player_strategy in x.keys():
        for player_strategy in plotList:
            if player_strategy == (m1,'ALLD'):
                ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '--',color = 'k')
            elif player_strategy == (m1,'TFT'):
                ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '--',color = 'r')
            elif player_strategy == (m2,'ALLD'):
                ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '-',color = 'k')
            elif player_strategy == (m2,'sTFT'):
                ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '-',color = 'r')
            else:
                print 'Error!'
    
            plt.xlabel('Time',{'weight':'bold','size':36})
            plt.ylabel('Frequency',{'weight':'bold','size':36})
    
            #plt.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy)
            #plt.title(player_strategy)
            #plt.show()
    
        plt.legend(prop={'size':24,'weight':'bold'},loc='best')
        plt.show()
    

    plotResults(('lysA','ilvE'))   
    # NOTE: ALLD for all species will always have the exact same fractions if the strategies are only 
    # ALLD, TFT and/or sTFT. This is because when ALLD faces with a TFT or sTFT the final outcome is 
    # like when ALLD faces ALLD. So, all their expected payoffs will be zero. This is not the case though
    # if ALLD faces ALLC 
