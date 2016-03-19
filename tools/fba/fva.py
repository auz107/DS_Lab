from __future__ import division
import sys, time
sys.path.append('../../')
from tools.fba.fbaTools import fbaTools
from tools.fba.fba import fba
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from coopr import pyomo

def fva(model, optimization_solver = 'gurobi', store_fva_flux_bounds = False, simulation_conditions = '', stdout_msgs = True, warnings = True):
    """
    Performs flux variability analysis

    INPUTS:
    -------
                    model: An instance of class model containing the information
                           about the metabolic model
      optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                           allowable choices are cplex and gurobi
       build_new_optModel: A parameter indicating whether a new pyomo optimizaiton model should be 
                           created (True) or or an existing one should be used (False). The options is useful
                           for the cases a model is already created and one just 
                           wants to change some model attributes (e.g., flux bounds)
                           and rerun FBA. Setting this parameter to False will save 
                           some runtime as the model need not to be created again.
                 maximize: If True, the objective is maximized. If False, the objective funciton is minimized
         store_fva_flux_bounds: If True, it stores the identified bounds on reaciton fluxes in fva_flux_bounds. Otherwise
                           they are stored in a dictionary whose keys are ids and values are a list of two elements 
                           in the form [fva_LB,fva_UB], whith fva_LB and fva_UB being the FVA LB and UB on fluxes 
    simulation_conditions: A string describing simulation conditions
                 warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                           screen or not. The default is True  
              stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                           Eligible values are True and False and the default is True 

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 12-15-2016 
    """   
    # store_fva_flux_bounds
    if not isinstance(store_fva_flux_bounds,bool):
        raise TypeError('store_fva_flux_bounds must be either True or False')
 
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

    for rxn in model.reactions:
        rxn.objective_coefficient = 0

    # Make an FBA model
    model.biomass_reaction.objective_coefficient = 1
    fba_model = fba(model = model, optimization_solver = optimization_solver, build_new_optModel = True, maximize = True, store_fva_flux_bounds = False, simulation_conditions = simulation_conditions, stdout_msgs = stdout_msgs, warnings = warnings)
    fba_model.run()   # Run it once to create optModel and constraints
    fba_model.stdout_msgs = False

    # A dictionary holding the FVA flux bounds
    fva_flux_bounds = {} 

    for rxn in model.reactions:
        rxn.objective_coefficient = 1

        # Minimize flux
        fba_model.maximize = False
        fba_model.optModel.del_component('objectiveFunc')
        fba_model.optModel.objectiveFunc = pyomo.Objective(rule=fba_model.objectiveFunc_rule, sense = pyomo.minimize)

        fba_model.run()
        if fba_model.solution['exit_flag'] == 'globallyOptimal':
            LB = fba_model.solution['objective_value']
        else:
            LB = None

        # Maximize flux
        fba_model.maximize = True
        fba_model.optModel.del_component('objectiveFunc')
        fba_model.optModel.objectiveFunc = pyomo.Objective(rule = fba_model.objectiveFunc_rule, sense = pyomo.maximize)

        fba_model.run()
        if fba_model.solution['exit_flag'] == 'globallyOptimal':
            UB = fba_model.solution['objective_value']
        else:
            UB = None

        # Store the results
        if store_fva_flux_bounds:
            rxn.fva_flux_bounds = [LB,UB]
        else:
            fva_flux_bounds[rxn.id] = [LB,UB]

        rxn.objective_coefficient = 0

        if stdout_msgs:
            print '{}: {}'.format(rxn.id,fva_flux_bounds[rxn.id])

    # Return fva_flux_bounds if store_fva_flux_bounds is not True
    if not store_fva_flux_bounds:
        return fva_flux_bounds        

