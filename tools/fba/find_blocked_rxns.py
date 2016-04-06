from __future__ import division
import sys, time
sys.path.append('../../')
from tools.globalVariables import *
from fva import fva
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from coopr import pyomo

def find_blocked_rxns(model, always_blocked_only = True, optimization_solver = default_optim_solver, save_to_model = False, results_filename = '', simulation_conditions = '', stdout_msgs = True, warnings = True):
    """
    Finds blocked reactions in the model under a given growth conditions specified by reactions
    flux bounds or reacions that are always blocked, i.e., when one allows unlimtted uptake 
    of all compounds with an exchange reaction in the model 

    INPUTS:
    -------
                    model: An instance of class model containing the information
                           about the metabolic model
      always_blocked_only: If True finds nnly alwasy blocked reactionss
      optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                           allowable choices are cplex and gurobi
             save_to_mode: If True, any reactions turned out to be essential is flagged by setting an attribute named "always_blocked"
                           to True, i.e., rxn.always_blocked = True where rxn is a reactions object in the model. If not always blocked 
                           rxn.always_blocked = False. 
         results_filename: A string containing the name of the file to save the results in. If an empty string is provided
    simulation_conditions: A string describing simulation conditions
                 warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                           screen or not. The default is True  
              stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                           Eligible values are True and False and the default is True 

    OUTPUT:
    -------
      always_blocked_rxns: A list of reaction ids containing the reactions that are essential under the examined conditions

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 04-01-2016 
    """   
    # always_blocked_only
    if not isinstance(always_blocked_only,bool):
        raise TypeError('always_blocked_only must be either True or False')

    # store_fva_flux_bounds
    if not isinstance(save_to_model,bool):
        raise TypeError('save_to_model must be either True or False')

    # optimization_solver
    if not isinstance(optimization_solver,str):
        raise TypeError('optimization_solver must be a string')
    elif optimization_solver.lower() not in ['gurobi','cplex']:
        raise ValueError('Invalid value for optimization_solver. Allowed choices are gurobi and cplex')

    # simulation_conditions 
    if not isinstance(simulation_conditions,str):
        raise TypeError('simulation_conditions must be a string')

    # warnings and stdout_msgs 
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True or False')
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True or False')

    # Store the original flux bounds 
    if always_blocked_only:
        rxns_flux_bounds = {}
        for rxn in model.reactions: 
            rxns_flux_bounds[rxn.id] = rxn.flux_bounds
        model.reset_flux_bounds()
        for exchrxn in [r for r in model.reactions if r.reversibility == 'exchange']:
            exchrxn.flux_bounds[0] = -1000

    # Perform fva
    fva_flux_bounds = fva(model, optimization_solver = optimization_solver, save_to_model = False, results_filename = '', simulation_conditions = simulation_conditions, warnings = True, stdout_msgs = False) 

    # List of blocked reactions 
    blocked_rxns = [] 
    for rxn in fva_flux_bounds.keys():
        if abs(fva_flux_bounds[rxn][0]) < 1e-6 and abs(fva_flux_bounds[rxn][1]) < 1e-6:
            blocked_rxns.append(rxn)   

    # Save to model
    if save_to_model:
        for rxn in model.reactions:
            if rxn.id in alwyas_blocked_rxns:
                if always_blocked_only:
                    rxn.always_blocked = True
                else:
                    rxn.blocked = True
            else: 
                if always_blocked_only:
                    rxn.always_blocked = False
                else:
                    rxn.blocked = False
  
    # Save results into a file 
    if results_filename != '':
        with open(results_filename,'w') as f:
            if always_blocked_only:
                f.write('alwyas_blocked_rxns = [\n')
            else: 
                f.write('blocked_rxns = [\n')
            for rxn in blocked_rxns: 
                f.write("'{}',\n".format(rxn))
            f.write(']')

    # Reset the original flux bounds
    if always_blocked_only:
        for rxn in model.reactions: 
            rxn.flux_bounds = rxns_flux_bounds[rxn.id]

    return blocked_rxns        

