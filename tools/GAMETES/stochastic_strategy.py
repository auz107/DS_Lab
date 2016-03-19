from __future__ import division
import re, sys, copy, time, random
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
    def __init__(self, name, player, probabilities, prob_distr_init = {}, eps = 0):
        """
        In the following, by "strategy" we mean a regular (not stochastic) strategy unless stated otherwise. 

        INPUTS:     
        -------
                player: A string containing the name of the player playing this strategy 
                  name: A string containing the name of the stochastic strategy 
                        (e.g., ALLD, ALLC, TFT, gTFT, WSLS).
         probabilities: A dictionary containing the probability of taking a particular strategy
                        by this player given the state of the game in the previous
                        stage. The keys and values of this dictionary are as follows:
                          Keys: A tuple of tuples containing the state of the game in the previous stage 
                                The first element of each inner tuple is the name of game players and 
                                the second element is the strategy it took in the previous stage. 
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
                                 D1, p2 played D2 and p3 played D3. probabilities for this state 
                                 is defined as follows:
                                 probabilities[(('p1','D1'),('p2','D2'),('p3','D3'))] = 
                                                {'D1':0.1,'C11':0.6,'C12':0.3} 
                                 Note that sum of these probabilities should be one.
                                 Remark: The user has the option to just input this field with all 
                                         probabilities initialized at zero. In this case there is 
                                         a function named assingProb that assigns the probabilites
                                         for well-known strategies such as ALLD, ALLC, TFT, gTFT, WSLS 
         prob_distr_init: A dictionary containing the initial probability distribution for this strategy. 
                            Keys: The name of a strategy that this player may take
                          Values: The probability that this player takes the given strategy 
                          For example, if the possible strategies of this player is 'D', 'C1' and 'C2', 
                          prob_distr_init = {'D':0,'C1':1,'C2':0}. Note that sum of the probablities must 
                          be one.
                     eps: A threshold that may be used instead of zero for assigning probabilities. For example,
                          ALLD may cooperate with probablity eps and defect with probability 1 - eps
        """
        # Name of the strategy (e.g., ALLD, ALLD, GRIM, TFT, gTFT, WSLS) 
        self.name = name

        # Name of the player 
        self.player = player

        # Theshould 
        self.eps = eps

        # Probabilities of taking a certain strategy by this player given a previous
        # state of the game 
        # But first check if the sum of the probabilites in probabilities is zero or one
        # If the sum is zero for all the states then the probabilites are assigned by the assingProb function
        # If the sum is not zero or one the code returns an error
        sumProb = dict([(state,0) for state in probabilities.keys()])

        for state in probabilities.keys():
            if abs(sum(probabilities[state].values()) - 1) <= 1e-6: 
                sumProb[state] = 1 
            elif abs(sum(probabilities[state].values()) - 0) <= 1e-6: 
                sumProb[state] = 0 
            else:
                raise userError('**Error! Sum of the probabilities for {}. This sum must be zero or one.'.format(sum(probabilities[state],values()))

        # If sum of the probabilites is one for all cases
        if sum(sumProb.values()) == len(sumProb.keys()):
            self.probabilities = deepcopy(probabilities)
      
        # If sum of the probabilities is zero for all cases
        elif sum(sumProb.values()) == 0:
            self.probabilities = deepcopy(probabilities)
            self.__assign_prob()
        else:
            raise userError('Sum of the probabilites is one for some cases and zero for some others\n')

        # Check sum of the probabilites again to make sure they add up to one
        sumProb = dict([(state,0) for state in self.probabilities.keys()])
        for state in probabilities.keys():
            if abs(sum(self.probabilities[state].values()) - 1) > 1e-6: 
                print raise userError('Sum of the probabilites for ({},{}) --> State = {}  {} is not one.'.format(self.player,self.name,state,probabilities[state])

        # Initial strategy
        self.prob_distr_init = prob_distr_init
        if prob_distr_init {}:
            if abs(1 - sum([v for v in prob_distr_init.values()])) > 1e-6:
                raise userError('Sum of the probablities in prob_distr_init is not one: {}'.format(sum([v for v in prob_distr_init.values()])))
        # Otherwise create it for some known strategies
        else:
            self.__create_prob_distr_init()

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # name
        if attr_name.lower() == 'name' and not isinstance(attr_value,str):
            raise TypeError('name must be a string')

        # eps phase time
        if attr_name.lower() == 'eps' and not isinstance(attr_value,int) and not isinstance(attr_value,float):
            raise TypeError('eps must be non-negative integer or float')
        elif attr_name.lower() == 'eps' and attr_value < 0:
            raise ValueError('eps must be a non-negative number')

        if attr_name.lower() == 'prob_distr_init' and not isinstance(attr_value,dict):
            raise TypeError('prob_distr_init must be a dictionary')

        self.__dict__[attr_name] = attr_value

    def __create_prob_distr_init(self):
        """
        Creates prob_distr_init for some known strategies
        """            
        # Possible strategies this player can take
        strategies = self.probabilities.values()[0].keys()
        self.prob_distr_init = dict([(s,None) for s in strategies])
        if self.name.lower() in ['alld','stft']:
            for strategy in strategies:
                if strategy.lower() == 'defect':
                    self.prob_distr_init[strategy] = 1 - self.eps
                else:                 
                    self.prob_distr_init[strategy] = self.eps/(len(strategies) - 1)
        elif self.name.lower() in ['allc','tft']:
            for strategy in strategies:
                if strategy.lower() == 'defect':
                    self.prob_distr_init[strategy] = self.eps
                else:                 
                    self.prob_distr_init[strategy] = (1 - self.eps)/(len(strategies) - 1)

    def __assign_prob(self):
       """
       This funciton assigns the probabilites for a number of defined strategies including
       ALLD (always defect), ALLC (always cooperate), TFT (Tit-for-tat), gTFT (generous TFT)
       and WSLS (Win-State Loose-Shipt) 
       """
       if self.name.lower() == 'alld':
           for state in self.probabilities.keys():
               for strategy in self.probabilities[state].keys():
                   if strategy.lower() == 'defect':
                       self.probabilities[state][strategy] = 1 - self.eps
                   else:
                       self.probabilities[state][strategy] = self.eps/(len(self.probabilities[state].keys()) - 1)

       elif self.name.lower() == 'allc':
           for state in self.probabilities.keys():
               for strategy in self.probabilities[state].keys():
                   if strategy.lower() == 'defect':
                       self.probabilities[state][strategy] = self.eps
                   else:
                       # If there are multiple cooperation strategies divide the
                       # probabilities equally
                       self.probabilities[state][strategy] = (1-self.eps)/(len(self.probabilities[state].keys()) - 1)

       elif self.name.lower() == 'tft' or self.name.lower() == 'stft':
           for state in self.probabilities.keys():
               stateDict = dict(state)
               # Cooperate if everybody cooperated
               if 'defect' not in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player.lower()]:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = self.eps
                       # Cooperate with equal probabilities if there are more than 
                       # cooperation strategy
                       else:
                           self.probabilities[state][strategy] = (1-self.eps)/(len(self.probabilities[state].keys()) - 1)

               # Else if anyone defected
               else:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = 1 - self.eps
                       else:
                           self.probabilities[state][strategy] = self.eps/(len(self.probabilities[state].keys()) - 1)


       elif self.name.lower() == 'gtft':
           # Probability to cooperate if any of the opponents defected in the previous stage
           self.gTFT_prob = 1/3

           for state in self.probabilities.keys():
               stateDict = dict(state)
               # Cooperate if everybody cooperated
               if 'defect' not in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player.lower()]:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = self.eps
                       # Cooperate with equal probabilities if there are more than 
                       # cooperation strategy
                       else:
                           self.probabilities[state][strategy] = (1-self.eps)/(len(self.probabilities[state].keys()) - 1)

               # Else if anyone defected
               else:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = self.gTFT_prob
                       else:
                           self.probabilities[state][strategy] = (1 - self.gTFT_prob)/(len(self.probabilities[state].keys()) - 1) 

       elif self.name.lower() == 'wsls':
           for state in self.probabilities.keys():
               stateDict = dict(state)

               # Cooperate if averybody (including yourself) cooperated in the previous stage (stay)
               if 'defect' not in [stateDict[player].lower() for player in stateDict.keys()]:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = self.eps
                       # Cooperate with equal probabilities if there are more than 
                       # cooperation strategy
                       else:
                           self.probabilities[state][strategy] = (1-self.eps)/(len(self.probabilities[state].keys()) - 1)

               # Cooperate if every body (including yourself) defected in the previous stage (shift) 
               elif [stateDict[player].lower() for player in stateDict.keys()].count('defect') == len([stateDict[player].lower() for player in stateDict.keys()]):
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = self.eps 
                       else:
                           self.probabilities[state][strategy] = (1-self.eps)/(len(self.probabilities[state].keys()) - 1)

               # Defect if you cooperated but at least one other player defected (shift)
               elif stateDict[self.player].lower() != 'defect' and 'defect' in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player.lower()]:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = 1 - self.eps 
                       else:
                           self.probabilities[state][strategy] = self.eps/(len(self.probabilities[state].keys()) - 1)

               # Defect if you defected but everybody else cooperated (stay)
               elif stateDict[self.player].lower() == 'defect' and 'defect' not in [stateDict[player].lower() for player in stateDict.keys() if player.lower() != self.player.lower()]:
                   for strategy in dict(self.probabilities[state]).keys():
                       if strategy.lower() == 'defect':
                           self.probabilities[state][strategy] = 1 - self.eps 
                       else:
                           self.probabilities[state][strategy] = self.eps/(len(self.probabilities[state].keys()) - 1)

               # Otherwise control for errors in the code
               else:
                   raise userError('Unknown game state in the previous stage!')

       # Uknown pre-defined strategy
       else:
           raise userError('Unknown strategy name {}'.format(self.name))
        


