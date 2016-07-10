from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import userError
from game import game
import numpy as np

def replicator_dynamics(payoff_matrix, x_init, time_range = [0,1,10], warnings = True, stdout_msgs = True):
    """
    This function implements the replicator dynamics using the payoff 
    matrix of a symmetric games. 

    INPUTS:
    ------
    payoff_matirx: 
    This is the payoff matrix of a symmetric game as follows, where keys are a tuple of strateiges encountering each other
    and values are dictionaries with keys and values being strategies and their respected payoff, respectively.
    Example: 
    [
     ('C','C'):{'C':10, 'C:10'},  --> When C faces C
     ('C','D'):{'C':0, 'D':15},   --> When C faces D
     ('D','D'):{'D':0, 'D':0}     --> When D faces D
    ]
    x_init: 
    A dictionary containing the frequency of each strategy at the beginning of simulation time 
    The keys and values of this dictionary are as follow:
          Keys: Elements of stoch_strategies
        Values: Fraction (frequency) of each stochastic strategy 

    time_range: 
    Range of the simulation time. This is a list with either two or three elements.
          Two-element format: [startTime,endTime]
        Three-element format: [startTime,timeStep,endTime]. If the 
                              time step is not given a time step of one is used
   
    OUTPUTS:
    ------- 
    x:
    The frequency of each strategy at different time points. This is s dictionary with keys and 
    values as follows: 
          Keys: Strategy names 
        Values: Another dictionary with keys and values being time and frequency, respectively.

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06/01/2016
    """
    # Check if the strategies given in x_init are the same as those in payoff_matrix 
    if set(x_init.keys()) != set([s for strategy_set in payoff_matrix.values() for s in strategy_set.keys()]):
        payoff_matrix_strats = set([s for strategy_set in payoff_matrix.values() for s in strategy_set.keys()])
        raise userError('strategy names in x_init do not match those in the payoff matrix:\nStrategies in x_init not in the payoff matrix: {}\nStrategies in the payoff matrix but not in x_init = {}'.format([s for s in x_init.keys() if s not in payoff_matrix_strats],[s for s in payoff_matrix_strats if s not in x_init.keys()]))

    # Check all initial x values add up to one
    if abs(sum(x_init.values()) - 1) > 1e-6:
        raise userError('Sum of the entries in x_init does not add up to one (sum = {})'.format(sum(x_init.values())))

    if len(time_range) == 3:
        # Initial time
        t0 = time_range[0]

        # Original dt. This mught be adjusted during the simulations
        dt = time_range[1]

        # Final simulation time
        tf = time_range[2]
    elif len(time_range) == 2:
        # Initial time
        t0 = time_range[0]

        # Original dt. This mught be adjusted during the simulations
        dt = 1

        # Final simulation time
        tf = time_range[1]   # Final simulation time
    else:
        userError("**ERROR! Invalid time_range (allowable size of the vectors are two or three)")

    # Strategy names
    strategy_names = x_init.keys()

    # Current time
    t = 0

    # A dictionary where keys are strategy names and values are another dictionary with keys and values 
    # being time points and fraction of of the strategy at that time point 
    x = dict([(s,{0:x_init[s]}) for s in strategy_names])

    while t <= tf:

        # Compute the average fitness of each strategy k: 
        # f_k = sum(j,a_kj*x_j) for j in {strategies} + 
        #       sum(j,a_kjl*x_j*x_l) for j,l in {strategies} + ... 
        strategies_ave_fitness = dict([(s,0) for s in strategy_names])

        for strategy in strategy_names:
            for payoffMatrix_key in [s for s in  payoff_matrix.keys() if strategy in s]:
                # Remove the first instance of strategy from payoffMatrix_key. For example, ('C','C') is 
                # converted to ('C') or ('D','C','D') is converted to ('C','D'). Note that if we have repeated elements
                # in a list, index gives the index of the first element.
                payoffMatrix_key_no_strategy = list(payoffMatrix_key)
                payoffMatrix_key_no_strategy.pop(payoffMatrix_key.index(strategy))
                strategies_ave_fitness[strategy] += payoff_matrix[payoffMatrix_key][strategy]*np.prod([x[s][t] for s in payoffMatrix_key_no_strategy])

        # Compute the expected fitness of the population (phi(x))
        phi = sum([strategies_ave_fitness[s]*x[s][t] for s in strategy_names]) 

        # Update the abundances
        # dx_k/dt = x_k*(f_x - phi) or (x_k(t + dt) - x_k(t))/dt = x_k(t)*(f_k - phi) or
        # x_k(t + dt) = dt*x_k(t)*(f_k - phi) + x_k(t)
        for strategy in strategy_names:
            x[strategy][t + dt] = max(0,dt*x[strategy][t]*(strategies_ave_fitness[strategy] - phi) + x[strategy][t])
        
        # Check if the sum is one
        sum_x = sum([x[strategy][t + dt] for strategy in strategy_names])
        if abs(sum_x - 1) > 1e-5:
             raise userError("Sum of x's at time {} is not one: sum =  {}".format(t+dt, sum_x))

        t += dt

    return x

