from __future__ import division
import sys, time
sys.path.append('../../')
from tools.globalVariables import *
from fbaTools import fbaTools
from fba import fba
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from coopr import pyomo

def fva(model, optimization_solver = default_optim_solver, save_to_model = False, results_filename = '', simulation_conditions = '', warmstart = False, warnings = True, stdout_msgs = True):
    """
    Performs flux variability analysis

    INPUTS:
    -------
                    model: An instance of class model containing the information
                           about the metabolic model
      optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                           allowable choices are cplex and gurobi
            save_to_model: If True, it stores the identified bounds on reaciton fluxes in fva_flux_bounds. Otherwise
                           they are stored in a dictionary whose keys are ids and values are a list of two elements 
                           in the form [fva_LB,fva_UB], whith fva_LB and fva_UB being the FVA LB and UB on fluxes 
         results_filename: A string containing the name of the file to save the results in. If an empty string is provided
                           the results are not saved to a file
    simulation_conditions: A string describing simulation conditions
            warmstart: Uses warmstart if True (works only with gurobi_ampl or cplexamp)
                 warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                           screen or not. The default is True  
              stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                           Eligible values are True and False and the default is True 

    OUTPUTS:
    --------
          fva_flux_bounds: A dictionary with keys being reactions ids and values beiing a list ot two elements containing
                           the fva flux bounds in the form [LB, UB]

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 04-01-2016 
    """   
    # save_to_model
    if not isinstance(save_to_model,bool):
        raise TypeError('save_to_model must be either True or False')
 
    # optimization_solver
    if not isinstance(optimization_solver,str):
        raise TypeError('optimization_solver must be a string')
    elif optimization_solver.lower() not in ['gurobi','cplex','gurobi_ampl','cplexamp']:
        raise ValueError('Invalid value for optimization_solver. Allowed choices are gurobi and cplex')

    # simulation_conditions 
    if not isinstance(simulation_conditions,str):
        raise TypeError('simulation_conditions must be a string')

    # warmstart 
    if not isinstance(warmstart,bool):
        raise TypeError('use_warmsart must be either True or False')

    # warnings and stdout_msgs 
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True or False')
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True or False')

    # If warmstart is True use gurobi_ampl
    if warmstart:
        optimization_solver = 'gurobi_ampl'

    # A dictionary holding the FVA flux bounds
    fva_flux_bounds = dict([(r.id,[None,None]) for r in model.reactions]) 

    #--- Minimize rxn fluxes ---
    for rxn in model.reactions:
        rxn.objective_coefficient = 0

    counter = 0
    for rxn in model.reactions:
        counter += 1
       
        rxn.objective_coefficient = 1

        if counter == 1:
            fba_model = fba(model = model, optimization_solver = optimization_solver, build_new_optModel = True, maximize = False, save_to_model = False, simulation_conditions = simulation_conditions, warmstart = warmstart, warmings = warnings, stdout_msgs = False, show_solver_output = True)

        # From counter 2 on, turn off build_new_optModel and preprocessing and turn on warmstart 
        elif counter == 2:
            fba_model.build_new_optModel = False

        # Redefine the objective function if counter > 1
        if counter > 1:
            fba_model.optModel.del_component('objectiveFunc')
            fba_model.optModel.objectiveFunc = pyomo.Objective(rule = fba_model.objectiveFunc_rule, sense = pyomo.minimize)

            # Supply the current solution as the warm start
            for j in fba_model.optModel.J:
                fba_model.optModel.v[j] = fba_model.solution['opt_rxnFluxes'][j]

        fba_model.run()
        if fba_model.solution['exit_flag'] == 'globallyOptimal':
            LB = fba_model.solution['objective_value']
        else:
            raise userError('Infeasbie fba problem in fva')

        # Store the results
        if save_to_model:
            rxn.fva_flux_bounds[0] = LB
        else:
            fva_flux_bounds[rxn.id][0] = LB

        rxn.objective_coefficient = 0

    #--- Maximize rxn flux ---
    for rxn in model.reactions:
        rxn.objective_coefficient = 0

    counter = 0
    for rxn in model.reactions:
        counter += 1
       
        rxn.objective_coefficient = 1

        if counter == 1:
            fba_model = fba(model = model, optimization_solver = optimization_solver, build_new_optModel = True, maximize = True, save_to_model = False, simulation_conditions = simulation_conditions, warmstart = warmstart, warnings = warnings, stdout_msgs = False, show_solver_output = True)

        # From counter 2 on, turn off build_new_optModel and preprocessing and turn on warmstart 
        elif counter == 2:
            fba_model.build_new_optModel = False

        # Redefine the objective function if counter > 1
        if counter > 1:
            fba_model.optModel.del_component('objectiveFunc')
            fba_model.optModel.objectiveFunc = pyomo.Objective(rule = fba_model.objectiveFunc_rule, sense = pyomo.maximize)

            # Supply the current solution as the warm start
            for j in fba_model.optModel.J:
                fba_model.optModel.v[j] = fba_model.solution['opt_rxnFluxes'][j]

        fba_model.run()
        if fba_model.solution['exit_flag'] == 'globallyOptimal':
            UB = fba_model.solution['objective_value']
        else:
            raise userError('Infeasbie fba problem in fva')

        # Store the results
        if save_to_model:
            rxn.fva_flux_bounds[1] = UB
        else:
            fva_flux_bounds[rxn.id][1] = UB

        rxn.objective_coefficient = 0

        # Save results into a file 
        if results_filename != '':
            with open(results_filename,'w') as f:
                f.write('fva_flux_bounds = {\n')
                for rxn in fva_flux_bounds.keys():
                    f.write("'{}':{},\n".format(rxn, fva_flux_bounds[rxn]))
                f.write('}')

    # Return fva_flux_bounds if save_to_model is not True
    return fva_flux_bounds        

