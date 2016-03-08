from __future__ import division
import re, sys, math, copy, time, random
from scipy.linalg import eig
from copy import deepcopy
sys.path.append('../../')
from tools.customError import *
from stochasticStrategy import *
from game import *

class iteratedGame_dFBA(object):
    """
    This is a general class simulating an iterated game with N players 
    Replicator dynamics is used to simualte the evolution of the community. 

    Ali R. Zomorrodi - Daniel Segre lab @ BU
    Last updated: October-08-2014
    """

    def __init__(self, game, iStrategies, x_init , sim_time, dt = None):
        """
        INPUTS:     
        -------
               game: An instance of the class game 
        iStrategies: A dictionary holding the information about the strategies 
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
                    raise customError('**Error! Unknonw strategy pair!\n')
        
        # If both strains play TFT, then the payoffs are like when they both play ALLC
        for iStrategy1 in self.iStrategies.keys():
            for iStrategy2 in self.iStrategies.keys():
                # Players of the same type get nothing when facing with each other 
                if iStrategy1[1].lower() == 'tft' and iStrategy2[1].lower() == 'tft':
                    self.payoff_iStrategies[(iStrategy1,iStrategy2)] = self.payoff_iStrategies[((iStrategy1[0],'ALLC'),(iStrategy2[0],'ALLC'))] 


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

        # Fill in the elements of M
        for row in self._gameStates:
            for col in self._gameStates: 
                M[self._gameStates.index(row),self._gameStates.index(col)] = self.transitionMatrix[(row,col)]

        # Sum of all columns for each row should be equal to one
        sumM = M.sum(axis=1)
        for m in sumM:
            if abs(sumM[m] - 1) >= 1e-6:
                print '**Error! Sum of the columns for each row is not one for the transition matrix of  ',self._iStrategies_curr.keys()
                print 'sum = ',sumM
                raise customError()           

        # Now find the stationary probabilities (left eigenvectors asociated with an 
        # eigenvalue of one)
        eigVals,eigVecs=eig(M,left=True,right=False)

        # Find the indices of eigen valueis, which is one (closest to one)         
        oneInds = [k for k in range(len(eigVals)) if abs(eigVals[k] - 1) < 1e-6]

        # There should be only one eigenvalue equal to one, if not issue a warning and use the
        # first index in oneEigValInds
        if len(oneInds) > 1:
            print '**Warning! More than one eigenvalues are close to one for ',self._iStrategies_curr.keys()

            print 'game states = ','\n'
            for row in self._gameStates:
                print row,'\n'
            print 'M = \n',M,'\n'
            print 'eigVals = \n',eigVals,'\n'
            print 'eigVecs = '
            for oneInd in oneInds:
                print eigVecs[:,oneInd]/sum(eigVecs[:,oneInd])
            print '-------------------------\n'

        stationaryProbArrays = []
        for oneInd in oneInds:
            # Normalize the eigenvector corresponding to the eigen value of one
            # This gives the stationary probabilities 
            eigVecOne_norm = eigVecs[:,oneInd]/sum(eigVecs[:,oneInd])

            # Check if any elements are negative
            negElems = [eigVecOne_norm[k] for k in range(max(eigVecOne_norm.shape)) if eigVecOne_norm[k] < 0]
            if len(negElems) >= 1:
                print '**Error! Negative elements in Statioanry probability: ',eigVecOne_norm
                raise customError('\n')
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

    #--- Load the games --- 
    print 'Loading the games ...\n'
    with open('../../EcoliPairs/results/auxoEcoliGames_all.pk','rb') as inputFile:
        games = pk.load(inputFile)

    # Choose a specific game
    sgame = games[('lysA','ilvE')]
    sgame.payoffMatrix[(('lysA', 'EX_ile-L(e)'), ('ilvE', 'EX_lys-L(e)')) ] = {'ilvE': 0.17617/5, 'lysA': 0.18026/5}

    # Create the keys for strategies_prob in stochastic strtategy class 
    prob_keys = sgame.payoffMatrix.keys() 

    # Initialize strategies_prob with z zero probability for all elements
    strategies_prob_lysA = dict([(k1,dict([(k2,0) for k2 in sgame.players_strategies['lysA']])) for k1 in prob_keys])
    strategies_prob_ilvE = dict([(k1,dict([(k2,0) for k2 in sgame.players_strategies['ilvE']])) for k1 in prob_keys])

    # Define iStrategies
    iStrategies = {}

    iStrategies[('lysA','ALLD')] = stochasticStrategy(player_name = 'lysA', strategy_name = 'ALLD',strategies_prob = strategies_prob_lysA)      
    iStrategies[('lysA','ALLC')] = stochasticStrategy(player_name = 'lysA', strategy_name = 'ALLC',strategies_prob = strategies_prob_lysA)      
    iStrategies[('lysA','TFT')] = stochasticStrategy(player_name = 'lysA', strategy_name = 'TFT',strategies_prob = strategies_prob_lysA)      

    iStrategies[('ilvE','ALLD')] = stochasticStrategy(player_name = 'ilvE', strategy_name = 'ALLD',strategies_prob = strategies_prob_ilvE)      
    iStrategies[('ilvE','ALLC')] = stochasticStrategy(player_name = 'ilvE', strategy_name = 'ALLC',strategies_prob = strategies_prob_ilvE)      
    iStrategies[('ilvE','TFT')] = stochasticStrategy(player_name = 'ilvE', strategy_name = 'TFT',strategies_prob = strategies_prob_ilvE)      

    # Initial concentrations
    C_init = {('lysA','ALLD'):(1/3)*(1e7),('lysA','ALLC'):(1/3)*(1e7),('lysA','TFT'):(1/3)*(1e7),('ilvE','ALLD'):(1/3)*(1e7),('ilvE','ALLC'):(1/3)*(1e7),('ilvE','TFT'):(1/3)*(1e7)}

    # Total number of initial cells
    C_init_tot = sum([C_init[k] for k in C_init.keys()])

    # Initial fractions
    x_init = {('lysA','ALLD'):C_init[('lysA','ALLD')]/C_init_tot,('lysA','ALLC'):C_init[('lysA','ALLC')]/C_init_tot,('lysA','TFT'):C_init[('lysA','TFT')]/C_init_tot,('ilvE','ALLD'):C_init[('ilvE','ALLD')]/C_init_tot,('ilvE','ALLC'):C_init[('ilvE','ALLC')]/C_init_tot,('ilvE','TFT'):C_init[('ilvE','TFT')]/C_init_tot}

    print 'x_init = ',x_init,'\n'

    if abs(sum([x_init[k] for k in x_init.keys()]) - 1) >= 1e-5:
        raise customError('sum of the initial fractions does not add to one\n')
    
    ig = iteratedGame(game = sgame, iStrategies = iStrategies, x_init = x_init , sim_time = 96, dt = 1)
    x = ig.run()

    fig, ax = plt.subplots()
    plt.rc('xtick', labelsize=32)
    plt.rc('ytick', labelsize=32)
    for player_strategy in x.keys():
        if player_strategy == ('lysA','ALLD'):
            ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '--',color = 'k')
        elif player_strategy == ('lysA','ALLC'):
            ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '--',color = 'b')
        elif player_strategy == ('lysA','TFT'):
            ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '--',color = 'r')
        elif player_strategy == ('ilvE','ALLD'):
            ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '-',color = 'k')
        elif player_strategy == ('ilvE','ALLC'):
            ax.plot(x[player_strategy].keys(),x[player_strategy].values(),label=player_strategy,linewidth = 8,linestyle = '-',color = 'b')
        elif player_strategy == ('ilvE','TFT'):
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


