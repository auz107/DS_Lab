from __future__ import division
import re, sys, math, copy, time, random
from scipy.linalg import eig
from copy import deepcopy
sys.path.append('../../')
from tools.userError import userError
from stochastic_strategy import stochastic_strategy
from game import game

class iteratedGame(object):
    """
    This is a general class simulating an iterated game with N players 
    Replicator dynamics is used to simulate the evolution of the community. 

    Ali R. Zomorrodi - Daniel Segre lab @ BU
    Last updated: 03-09-2016
    """
    def __init__(self, game, stoch_strategies, x_init, time_range = [0,1,10], ignore_same_player_encounter = False, warnings = True, stdout_msgs = True):
        """
        INPUTS:     
        -------
                    game: An instance of the class game 
        stoch_strategies: A list of stochastic strategies. For example, if have two players p1 and p2 in the classic game 
                          and in the iterated game we have two versions of p1 one playing ALLD and the other TFT and also 
                          two versions of p2 one playing ALLD and the other ALLC, then the set of stochastic strategies will
                          by a list of instances of four stochastic_strategy for (p1,ALLD), (p1,TFT), (p2, ALLD) and (p2,ALLC)
                  x_init: A dictionary containing the frequency of each stochastic strategy at the beginning of simulation time 
                            Keys: Elements of stoch_strategies
                          Values: Fraction (frequency) of each stochastic strategy 
              time_range: Range of the simulation time. This is a list with 
                          either two or three elements.
                          Two-element format: [startTime,endTime]
                          Three-element format: [startTime,timeStep,endTime]. If the 
                          time step is not given a time step of one is used
        ignore_same_player_encounter: Ignores whether we should ignore (assign a zero payoff) the encouter of the same players with 
                          any stochastic strategies (True) or not (False). This is sometimes relevant, for example, when we 
                          deal with two E. coli mutants cross feeding each other the amino acids they cannot produce. In this case, if
                          for example, a lysA playing ALLD faces with another lysA mustant (no matter if the play the same or a different
                          stochastic strategy), they will gain nothing and zero is assigned as the payoff. This is mainly because in this
                          case cooperate and defect strategies are essentially defined in connection to a different mutant. 
        """
        # Game
        self.game = game

        # Payoff matrix of the original game
        self._payoff_matrix = self.game.payoff_matrix

        # Total possible states of the game
        # Example: For (lysA, ilvE) E. coli mutants, we have:
        # [(('lysA', 'EX_ile-L(e)'), ('ilvE', 'Defect')), (('lysA', 'EX_ile-L(e)'), ('ilvE', 'EX_lys-L(e)')), 
        #  (('lysA', 'Defect'), ('ilvE', 'EX_lys-L(e)')), (('lysA', 'Defect'), ('ilvE', 'Defect'))]
        self._game_states = self._payoff_matrix.keys() 

        # Players' strategies for iterated games
        self.stoch_strategies = deepcopy(stoch_strategies)

        # Initial concentrations
        self.x_init = x_init

        if len(time_range) == 3:
            # Initial time
            self._t0 = time_range[0]

            # Original dt. This mught be adjusted during the simulations
            self._dt = time_range[1]

            # Final simulation time
            self._tf = time_range[2]
        elif len(time_range) == 2:
            # Initial time
            self._t0 = time_range[0]

            # Original dt. This mught be adjusted during the simulations
            self._dt = 1

            # Final simulation time
            self._tf = time_range[1]   # Final simulation time
        else:
            userError("**ERROR! Invalid time_range (allowable size of the vectors are two or three)")

        # Ignore the encounter of the same players
        self.ignore_same_player_encounter = ignore_same_player_encounter

        # warnings and stdout_msgs
        self.warnings = warnings
        self.stdout_msgs = stdout_msgs
    
        # Dictionaries where keys are time points and values are the average fitness
        # or frequency (relative abundance) of each stochastic strategy 
        for stoch_strategy in self.stoch_strategies:
            stoch_strategy.fitness = {}
            stoch_strategy.frequency = {}

        # A dictionary with keys and values being the time pionts and the average fitness
        # of the community
        self.phi = {}

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # game
        if attr_name.lower() == 'game' and not isinstance(attr_value,game):
            raise TypeError('game must be an instance of class game')

        # Stochastic strategies in the iterated game 
        if attr_name.lower() == 'stoch_strategies' and not isinstance(attr_value,list):
            raise TypeError('stoch_strategies must be a list')

        # Players' frequencies 
        if attr_name.lower() == 'x_init' and not isinstance(attr_value,dict):
            raise TypeError('x_init must be a dictionary')

        # time range
        if attr_name.lower() == 'time_range' and not isinstance(attr_value,list):
            raise TypeError('time_range must be a list')

        # ignore_same_player_encounter
        if attr_name.lower() == 'ignore_same_player_encounter' and not isinstance(attr_value,bool):
            raise TypeError('ignore_same_player_encounter must be either True or False')

        # warnings and stdout_msgs 
        if attr_name.lower() == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError('warnings must be either True or False')
        if attr_name.lower() == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError('stdout_msgs must be either True or False')

        self.__dict__[attr_name] = attr_value

    def create_iteratedGame_payoffMatrix(self):
        """
        This function creates the payoff matrix of the iterated game with stochastic strategies.
        The output of this function is a dictionary called iteratedGame_payoffMatrix
        with keys and values as follows:
          keys: A tuple with two elements, where each element is a key of
                self.stoch_strategies, which is a tuple with the first element being
                the name of the player and the second the name of the iterated
                strategy it takes
        values: Another dictionary with keys and values as follows: 
                  keys: Name of the player 
                values: The payoff of the player
        """
        self.iteratedGame_payoffMatrix = {}

        for stoch_strategy1 in self.stoch_strategies:
            for stoch_strategy2 in self.stoch_strategies:
                # Players of the same type get nothing when facing with each other 
                if self.ignore_same_player_encounter and stoch_strategy1.player.lower() == stoch_strategy2.player.lower():
                    self.iteratedGame_payoffMatrix[(stoch_strategy1,stoch_strategy2)] = {stoch_strategy1:0, stoch_strategy2:0}

                elif self.stoch_strategies.index(stoch_strategy2) > self.stoch_strategies.index(stoch_strategy1):
                    self._stoch_strategies_pairs = [stoch_strategy1,stoch_strategy2]
 
                    # Create the Markov chain transition matrix
                    self._create_transition_matrix()

                    # Find the stationary probabilities
                    self._find_stationary_probs()

                    # Compute the expected payoff for each stoch_strategy1 and stoch_strategy2 when these two 
                    # strategies facing each other
                    stoch_strategy1_payoff, stoch_strategy2_payoff = [], []
                    for stationary_probs in self._stationary_probs_list:
                        stoch_strategy1_payoff.append(sum([stationary_probs[game_state]*self._payoff_matrix[game_state][stoch_strategy1.player] for game_state in stationary_probs.keys()]))
                        stoch_strategy2_payoff.append(sum([stationary_probs[game_state]*self._payoff_matrix[gaem_state][stoch_strategy2.player] for game_state in stationary_probs.keys()]))

                    # If there are more than one set of stationary probabilities, the payoffs are
                    # computed by taking the average over all
                    self.iteratedGame_payoffMatrix[(stoch_strategy1,stoch_strategy2)] = {stoch_strategy1:sum(stoch_strategy1_payoff)/len(stoch_strategy1_payoff),stoch_strategy2:sum(stoch_strategy2_payoff)/len(stoch_strategy2_payoff)} 
                    if self.warnings and len(self._stationary_probs_list) > 1:
                        print 'WARNING! more than one stationary probability distribution was found for {}. The payoffs were averaged.'.format(((stoch_strategy1.player,stoch_strategy1.name),(stoch_strategy2.player,stoch_strategy2.name)))

                # Otherwise this strategy pair has already been considered
                elif self.stoch_strategies.index(stoch_strategy2) < self.stoch_strategies.index(stoch_strategy1):
                    self.iteratedGame_payoffMatrix[(stoch_strategy1,stoch_strategy2)] = self.iteratedGame_payoffMatrix[(stoch_strategy2,stoch_strategy1)] 

                else:
                    raise userError('Unknonw stochastic strategy pair!\n')
        
    def _create_transition_matrix(self):
        """
        This funciton creates the Markov chain transition matrix for a pairs of stochastic 
        strategies (ALLD, ALLC, TFT, etc) stored in self._stoch_strategies_pairs, which is a list created by 
        funciton create_iteratedGame_payoffMatrix. The elements of this list are instances of the class stochastic_strategy 
        and pairs are selected from self.stoch_strategies. 
        The output of this function is a dictionary with the following keys and values:
          Keys: A tuple in the form (gameState_prev,gameState_curr) containing the states of the game in the previous and current tiime point
        Values: The probability of being in the current state of the game given a previous state for the game
        """
        # Transition matrix
        self.transitionMatrix = {}

        # Loop over the game states in the previous time point
        # Example: For (lysA, ilvE) E. coli mutants, we have:
        # [(('lysA', 'EX_ile-L(e)'), ('ilvE', 'Defect')), (('lysA', 'EX_ile-L(e)'), ('ilvE', 'EX_lys-L(e)')), 
        #  (('lysA', 'Defect'), ('ilvE', 'EX_lys-L(e)')), (('lysA', 'Defect'), ('ilvE', 'Defect'))]
        # FOr Yeast we have: [(('yeast1','C'),('yeast2','C')),(('yeast1','C'),('yeast2','D')),
        # (('yeast1','D'),('yeast2','C')),(('yeast1','D'),('yeast2','D'))]
        # NOTE: Even if both players are of the same type (e.g., if both are yeast) different names should be given to each player. 
        #       In such cases different names should be also assigned to stochastic_strategy.player as well. These names should be such
        #       that they can be uniquely found in game.players_names. This is because, we need to find a way of what strategy this current
        #       players took in the previous state of the game. 
        for gameState_prev in self._game_states: 
            # Loop over the game states in the current time point
            for gameState_curr in self._game_states: 
                counter = 0
                for player in self.game.players_names: 
                    # Strategy taken by this player in gameState_curr
                    player_strategy = dict(list(gameState_curr))[player]

                    # Find the name of stochastic strategy this player has taken
                    player_stoch_strategy = [stoch_strat for stoch_strat in self._stoch_strategies_pairs if stoch_strat.player.lower() == player.lower()]
                    if len(player_stoch_strategy) == 1:
                        player_stoch_strategy = player_stoch_strategy[0]
                    else:
                        raise userError('A unique stochastic strategy was not found for player {}. Stochastic strategies found are: {}'.format(player,[s.name for s in player_stoch_strategy]))

                    # The probability of taking this strategy by this player accordiing to the
                    # definition of the stochastic strategy it takes
                    prob = player_stoch_strategy.probabilities[gameState_prev][player_strategy]
                    if counter == 0:
                        self.transitionMatrix[(gameState_prev,gameState_curr)] = prob
                    else: 
                        self.transitionMatrix[(gameState_prev,gameState_curr)] = prob*self.transitionMatrix[(gameState_prev,gameState_curr)]

                    counter += 1

    def _find_stationary_probs(self):
        """
        This function finds the stationary probabilities of being in each state of the
        game using the Markov chain transition matrix. p = p*M. The output is a list of 
        dictionaries self._stationary_probs_list, where the keys and values of each dictionary
        are as follow:
          keys: States of the game (elements of self._game_states) 
        values: The stationay probability of being in each state when a particular
                set of stochastic strategies face each other
        """
        # Initialize the transition matrix in the form of a numpy array
        M = np.zeros((len(self._game_states),len(self._game_states)))

        # Fill in the elements of transition matrix M 
        # Note that self.transitionMatrix is in the form:
        # self.transitionMatrix[(gameState_prev,gameState_curr)] = prob
        for row in self._game_states:
            for col in self._game_states: 
                M[self._game_states.index(row),self._game_states.index(col)] = self.transitionMatrix[(row,col)]

        # Sum of all columns for each row should be equal to one
        sumM = M.sum(axis=1)
        for m in sumM:
            if abs(sumM[m] - 1) >= 1e-6:
                raise userError('Sum of the columns for each row is not one for the transition matrix of  ' + srt(self._stoch_strategies_pairs.keys()) + ': sum = ' + str(sumM))

        # Now find the stationary probabilities (left eigenvectors asociated with an 
        # eigenvalue of one)
        eigVals,eigVecs = eig(M, left = True, right = False)

        # Find the indices of eigen valueis, which is one (closest to one)         
        oneInds = [k for k in range(len(eigVals)) if abs(eigVals[k] - 1) < 1e-6]

        # Eigeven value of one selected for TFT strategies (if applicable)
        tft_oneInds = []

        # There should be only one eigenvalue equal to one, if not issue a warning and use the
        # first index in oneEigValInds
        if len(oneInds) > 1:
            if self.warnings:
                print 'WARNING! More than one eigenvalues are close to one for '.format([s.name for s in self._stoch_strategies_pairs])
                print '\tM = \n{}\n'.format(M)
                print '\teigVals = \n\t{}\n'.format(eigVals)
                print '\teigVecs = '
                for oneInd in oneInds:
                    print '\t{}\n'.format(eigVecs[:,oneInd]/sum(eigVecs[:,oneInd]))

        stationary_prob_arrays = []
        for oneInd in oneInds:
            # Normalize the eigenvector corresponding to the eigen value of one
            # This gives the stationary probabilities 
            eigVecOne_norm = eigVecs[:,oneInd]/sum(eigVecs[:,oneInd])

            # Check if any elements are negative
            negElems = [eigVecOne_norm[k] for k in range(max(eigVecOne_norm.shape)) if eigVecOne_norm[k] < 0 and abs(eigVecOne_norm[k]) > 1e-6]
            if len(negElems) >= 1:
                raise userError('Negative elements in Statioanry probabilities: {}'.format(eigVecOne_norm)
            else:
                stationary_prob_arrays.append(eigVecOne_norm)

        # Convert the numpy array to a list of python dictionaries
        self._stationary_probs_list = [] 
        for stationary_prob_array in stationary_prob_arrays:
            self._stationary_probs_list.append(dict([(self._game_states[k],stationary_prob_array[k]) for k in range(len(self._game_states))]))

    def replicatorDynamics(self):
        """
        This function implements the replicator dynamics using the payoff 
        matrix of the stochastic strategies. The output is the frequencies of
        community members at time t + dt. This information is held in 
        a dictionary self.x_tpdt, with the following keys and values:
          Keys: The same as keys of self._stoch_strategies
        Values: Another dictionary with keys and values being time and frequency, respectively.
        """
        # Compute the average fitness of each stochastic strategy (f_k = sum(j,a_kj*x_j) for j in {stochastic strategies} 
        for stoch_strategy in self.stoch_strategies:
            stoch_strategy.fitness[self._t] = sum([self.iteratedGame_payoffMatrix[stoch_strategy_pair][stoch_strategy]*stoch_strategy_pair[1].frequency[self._t] for stoch_strategy_pair in self.iteratedGame_payoffMatrix.keys() if stoch_strategy_pairs[0] == stoch_strategy])             
 
        # Compute the expected fitness of the population (phi(x))
        self.phi{self._t} = sum([self.stoch_strategy.fitness[self._t]*self.stoch_strategy.frequency[self._t] for stoch_strategy in self.stoch_strategies])
 
        # Update the abundances
        # dx_k/dt = x_k*(f_x - phi) or (x_k(t + dt) - x_k(t))/dt = x_k*(f_k - phi)
        for stoch_strategy in self.stoch_strategies:
            stoch_strategy.frequency[self._t + self._dt] = max(0,self._dt*stoch_strategy.frequency[self._t]*(stoch_strategy.fitness[self._t] - self.phi[self._t]) + stoch_strategy.frequency[self._t])

    def run(self):
        """
        This function runs the entire iterated game. The output is 
        a dictionary self.x with keys and values as follows:
          keys: The same as keys of self.x_t or self.stoch_strategies (a tuple where the first 
                element is the name of the player and the second element is the name of 
                iterated strategy it plays)
        Values: A dictionary with keys and values as follows:
                  Keys: Time points
                Values: Frequency at that time point
        """
        # Create the payoff matrix of iterated strategies
        self.create_iteratedGame_payoffMatrix()

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
    from stochastic_strategy import *
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
    
        print 'payoff_matrix = \n'
        for k in sgame.payoff_matrix.keys():
            print k,' = ',sgame.payoff_matrix[k]
        print '\n'
    
        # Create the keys for probabilities in stochastic strtategy class 
        # Example: [(('lysA', 'EX_ile-L(e)'), ('ilvE', 'Defect')), (('lysA', 'EX_ile-L(e)'), ('ilvE', 'EX_lys-L(e)'))
        prob_keys = sgame.payoff_matrix.keys() 
    
        # Initialize probabilities with z zero probability for all elements
        # Example: {(('lysA', 'EX_ile-L(e)'), ('ilvE', 'Defect')): {'Defect': 0, 'EX_ile-L(e)': 0}, 
        #           (('lysA', 'EX_ile-L(e)'), ('ilvE', 'EX_lys-L(e)')): {'Defect': 0, 'EX_ile-L(e)': 0}, 
        #           (('lysA', 'Defect'), ('ilvE', 'EX_lys-L(e)')): {'Defect': 0, 'EX_ile-L(e)': 0}, 
        #           (('lysA', 'Defect'), ('ilvE', 'Defect')): {'Defect': 0, 'EX_ile-L(e)': 0}} 
        probabilities_m1 = dict([(k1,dict([(k2,0) for k2 in sgame.players_strategies[m1]])) for k1 in prob_keys])
        probabilities_m2 = dict([(k1,dict([(k2,0) for k2 in sgame.players_strategies[m2]])) for k1 in prob_keys])
    
        # Define stoch_strategies
        stoch_strategies = {}
    
        stoch_strategies[(m1,'ALLD')] = stochastic_strategy(player = m1, name = 'ALLD',probabilities = probabilities_m1)      
        stoch_strategies[(m1,'TFT')] = stochastic_strategy(player = m1, name = 'TFT',probabilities = probabilities_m1)      
    
        stoch_strategies[(m2,'ALLD')] = stochastic_strategy(player = m2, name = 'ALLD',probabilities = probabilities_m2)      
        stoch_strategies[(m2,'sTFT')] = stochastic_strategy(player = m2, name = 'sTFT',probabilities = probabilities_m2)      
    
        # Initial concentrations
        strat_num = 2
        C_init = {(m1,'ALLD'):(1/strat_num)*(1e7),(m1,'TFT'):(1/strat_num)*(1e7),(m2,'ALLD'):(1/strat_num)*(1e7),(m2,'sTFT'):(1/strat_num)*(1e7)}
    
        # Total number of initial cells
        C_init_tot = sum([C_init[k] for k in C_init.keys()])
    
        # Initial fractions
        x_init = {(m1,'ALLD'):C_init[(m1,'ALLD')]/C_init_tot,(m1,'TFT'):C_init[(m1,'TFT')]/C_init_tot,(m2,'ALLD'):C_init[(m2,'ALLD')]/C_init_tot,(m2,'sTFT'):C_init[(m2,'sTFT')]/C_init_tot}
    
        if abs(sum([x_init[k] for k in x_init.keys()]) - 1) >= 1e-5:
            raise userError('sum of the initial fractions does not add to one\n')
        
        ig = iteratedGame(game = sgame, stoch_strategies = stoch_strategies, x_init = x_init , sim_time = 96, dt = 1)
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
