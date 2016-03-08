from __future__ import division
import re, sys, math, copy, time, random
import numpy as np
from copy import deepcopy
sys.path.append('../../')
from tools.userError import *


class stochasticStrategy(object):
    """
    This is a general class holding information for a stochastic strategy 

    Ali R. Zomorrodi - Daniel Segre lab @ BU
    Last updated: 06-11-2015
    """
    def __init__(self,player_name, strategy_name, strategies_prob, prob_distr_init = None, eps = 0):
        """
        INPUTS:     
        -------
            player_name: A string containing the name of the player playing this strategy 
          strategy_name: A string containing the name of the stochastic strategy 
                         (e.g., ALLD, ALLC, TFT, gTFT, WSLS).
        strategies_prob: A dictionary containing the probability of taking a particular strategy
                         by this player given the state of the game in the previous
                         stage. The keys and values of this dictionary are as follows:
                           Keys: A tuple of tuples containing the state of the game in the previous stage 
                                 The first element of each inner tuple is the name of game players and 
                                 the second element is the strategy it has taken in the previous stage. 
                         Values: Another dictionary with keys and values as follows:
                                   Keys: The name of a strategy that this player may take
                                 Values: The probability that this player takes the strategy given 
                                         the state of the game in the previous stage as specified 
                                         in the key of the main dictionary
                        Example: Let assume the original game has three players p1, p2 and p3 
                                 Each player can take the following strategies:
                                 p1: D1, C11, C12    p2: D2, C21, C22    p3: D3, C3
                                 Additionally assume that the strategy in the current stage for player 
                                 p1 is defined over all possible 18 states of the game in the previous
                                 stage. As an example if the state of the game was: p1 played
                                 D1, p2 played D2 and p3 played D3. strategies_prob for this state 
                                 is defined as follows:
                                 strategies_prob[(('p1','D1'),('p2','D2'),('p3','D3'))] = 
                                                 {'D1':0.1,'C11':0.6,'C12':0.3} 
                                 Note that sum of these probabilities should be one.
                                 Remark: The user has the option to just input this field with all 
                                         probabilities initialized at zero. In this case there is 
                                         a function named assingProb that assigns the probabilites
                                         for well-known strategies such as ALLD, ALLC, TFT, gTFT, WSLS 
         prob_distr_init: A dictionary containing the initial probability distribution for this strategy. 
                            Keys: The name of a strategy that this player may take
                          Values: The probability that this player takes the strategy given 
                          For example, if the possible strategies of this player is 'D', 'C1' and 'C2', 
                          prob_distr_init = {'D':0,'C1':1,'C2':0}. Note that sum of the probablities must 
                          be one.
                     eps: A threshold that may be used instead of zero for assigning probabilities. For example,
                          ALLD may cooperate with probablity eps and defect with probability 1 - eps
        """

        # Name of the strategy (e.g., ALLD, ALLD, GRIM, TFT, gTFT, WSLS) 
        self.strategy_name = strategy_name

        # Name of the player 
        self.player_name = player_name

        # Theshould 
        self.eps = eps

        # Probabilities of taking a certain strategy by this player given a previous
        # state of the game 
        # But first check if the sum of the probabilites in strategies_prob is zero or one
        # If the sum is zero for all the states then the probabilites are assigned by the assingProb function
        # If the sum is not zero or one the code returns an error
        sumProb = dict([(state,0) for state in strategies_prob.keys()])

        for state in strategies_prob.keys():
            if abs(sum(strategies_prob[state].values()) - 1) <= 1e-6: 
                sumProb[state] = 1 
            elif abs(sum(strategies_prob[state].values()) - 0) <= 1e-6: 
                sumProb[state] = 0 
            else:
                print '**Error! Sum of the probabilities for ',sum(strategies_prob[state],values()),' (this sum must be one)'
                raise userError()


        # If sum of the probabilites is one for all cases
        if sum(sumProb.values()) == len(sumProb.keys()):
            self.strategies_prob = deepcopy(strategies_prob)
      
        # If sum of the probabilities is zero for all cases
        elif sum(sumProb.values()) == 0:
            self.strategies_prob = deepcopy(strategies_prob)
            self.__assignProb()
        else:
            raise userError('**Error! Sum of the probabilites is one for some cases and zero for some others\n')

        # Check sum of the probabilites again to make sure they add up to one
        sumProb = dict([(state,0) for state in self.strategies_prob.keys()])
        for state in strategies_prob.keys():
            if abs(sum(self.strategies_prob[state].values()) - 1) > 1e-6: 
                print '**Error! Sum of the probabilites for (',self.player_name,',',self.strategy_name,') --> State = ',state,'  ',strategies_prob[state],' is not one'
                raise userError()

        # Possible states of the game (list of dictionaries)
        self.gameStates = [dict(k) for k in  strategies_prob.keys()] 

        # Names of the game players (pick any arbitrary game state, which is a dictionary  
        # and find its keys)
        self.playersNames = self.gameStates[0].keys()

        # Initial strategy
        if prob_distr_init is not None:
            if abs(1 - sum([v for v in prob_distr_init.values()])) > 1e-6:
                raise userError('Sum of the probablities in prob_distr_init is not one: ' + str(sum([v for v in prob_distr_init.values()])))
            else:
                self.prob_distr_init = prob_distr_init
        # Otherwise create it for some known strategies
        else:
            self.__create_prob_distr_init()

    def __create_prob_distr_init(self):
        """
        Creates __create_prob_distr_init for some known strategies
        """            
        # Possible strategies this player can take
        strategies = self.strategies_prob.values()[0].keys()
        self.prob_distr_init = dict([(s,None) for s in strategies])
        if self.strategy_name.lower() in ['alld','stft']:
            for strategy in strategies:
                if strategy.lower() == 'defect':
                    self.prob_distr_init[strategy] = 1 - self.eps
                else:                 
                    self.prob_distr_init[strategy] = self.eps/(len(strategies) - 1)
        elif self.strategy_name.lower() in ['allc','tft']:
            for strategy in strategies:
                if strategy.lower() == 'defect':
                    self.prob_distr_init[strategy] = self.eps
                else:                 
                    self.prob_distr_init[strategy] = (1 - self.eps)/(len(strategies) - 1)


    def __assignProb(self):
       """
       This funciton assign the probabilites for a number of defined strategies including
       ALLD (always defect), ALLC (always cooperate), TFT (Tit-for-tat), gTFT (generous TFT)
       and WSLS (Win-State Loose-Shipt) 
       """

       if self.strategy_name.lower() == 'alld':
           for state in self.strategies_prob.keys():
               for strategy in self.strategies_prob[state].keys():
                   if strategy.lower() == 'defect':
                       self.strategies_prob[state][strategy] = 1 - self.eps
                   else:
                       self.strategies_prob[state][strategy] = self.eps/(len(self.strategies_prob[state].keys()) - 1)

       elif self.strategy_name.lower() == 'allc':
           for state in self.strategies_prob.keys():
               for strategy in self.strategies_prob[state].keys():
                   if strategy.lower() == 'defect':
                       self.strategies_prob[state][strategy] = self.eps
                   else:
                       # If there are multiple cooperation strategies divide the
                       # probabilities equally
                       self.strategies_prob[state][strategy] = (1-self.eps)/(len(self.strategies_prob[state].keys()) - 1)

       elif self.strategy_name.lower() == 'tft' or self.strategy_name.lower() == 'stft':
           for state in self.strategies_prob.keys():
               stateDict = dict(state)
               # Cooperate if everybody cooperated
               if 'defect' not in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player_name.lower()]:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = self.eps
                       # Cooperate with equal probabilities if there are more than 
                       # cooperation strategy
                       else:
                           self.strategies_prob[state][strategy] = (1-self.eps)/(len(self.strategies_prob[state].keys()) - 1)

               # Else if anyone defected
               else:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = 1 - self.eps
                       else:
                           self.strategies_prob[state][strategy] = self.eps/(len(self.strategies_prob[state].keys()) - 1)


       elif self.strategy_name.lower() == 'gtft':
           # Probability to cooperate if any of the opponents defected in the previous stage
           self.gTFT_prob = 1/3

           for state in self.strategies_prob.keys():
               stateDict = dict(state)
               # Cooperate if everybody cooperated
               if 'defect' not in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player_name.lower()]:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = self.eps
                       # Cooperate with equal probabilities if there are more than 
                       # cooperation strategy
                       else:
                           self.strategies_prob[state][strategy] = (1-self.eps)/(len(self.strategies_prob[state].keys()) - 1)

               # Else if anyone defected
               else:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = self.gTFT_prob
                       else:
                           self.strategies_prob[state][strategy] = (1 - self.gTFT_prob)/(len(self.strategies_prob[state].keys()) - 1) 

       elif self.strategy_name.lower() == 'wsls':
           for state in self.strategies_prob.keys():
               stateDict = dict(state)

               # Cooperate if averybody (including yourself) cooperated in the previous stage (stay)
               if 'defect' not in [stateDict[player].lower() for player in stateDict.keys()]:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = self.eps
                       # Cooperate with equal probabilities if there are more than 
                       # cooperation strategy
                       else:
                           self.strategies_prob[state][strategy] = (1-self.eps)/(len(self.strategies_prob[state].keys()) - 1)

               # Cooperate if every body (including yourself) defected in the previous stage (shift) 
               elif [stateDict[player].lower() for player in stateDict.keys()].count('defect') == len([stateDict[player].lower() for player in stateDict.keys()]):
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = self.eps 
                       else:
                           self.strategies_prob[state][strategy] = (1-self.eps)/(len(self.strategies_prob[state].keys()) - 1)

               # Defect if you cooperated but at least one other player defected (shift)
               elif stateDict[self.player_name].lower() != 'defect' and 'defect' in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player_name.lower()]:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = 1 - self.eps 
                       else:
                           self.strategies_prob[state][strategy] = self.eps/(len(self.strategies_prob[state].keys()) - 1)

               # Defect if you defected but everybody else cooperated (stay)
               elif stateDict[self.player_name].lower() == 'defect' and 'defect' not in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player_name.lower()]:
                   for strategy in dict(self.strategies_prob[state]).keys():
                       if strategy.lower() == 'defect':
                           self.strategies_prob[state][strategy] = 1 - self.eps 
                       else:
                           self.strategies_prob[state][strategy] = self.eps/(len(self.strategies_prob[state].keys()) - 1)

               # Otherwise control for errors in the code
               else:
                   raise userError('**Error! Unknown game state in the previous stage!')

       # Uknown pre-defined strategy
       else:
           print '**Error! Unknown strategy name ',self.strategy_name,'\n'
           raise userError('**Error! Unknown strategy name ')
        


