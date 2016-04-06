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

def find_essential_rxns(model, viability_thr = 1, optimization_solver = default_optim_solver, save_to_model = False, results_filename = '', simulation_conditions = '', stdout_msgs = True, warnings = True):
    """
    Finds essential reactions in the model under a given uptake conditions 

    INPUTS:
    -------
                    model: An instance of class model containing the information
                           about the metabolic model
            viability_thr: Viability threshold that is the percentage of max biomass below which is considered as no growth
      optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                           allowable choices are cplex and gurobi
             save_to_mode: If True, any reactions turned out to be essential is flagged by setting an attribute named "essential"
                           to True, i.e., rxn.essential = True where rxn is a reactions object in the model. If not essential 
                           rxn.essential = False. 
         results_filename: A string containing the name of the file to save the results in. If an empty string is provided
    simulation_conditions: A string describing simulation conditions
                 warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                           screen or not. The default is True  
              stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                           Eligible values are True and False and the default is True 

    OUTPUT:
    -------
          essential_rxns: A list of reaction ids containing the reactions that are essential under the examined conditions

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 04-01-2016 
    """   
    # store_fva_flux_bounds
    if not isinstance(save_to_model,bool):
        raise TypeError('save_to_model must be either True or False')

    # Viability threshold
    if not isinstance(viability_thr,int) and not isinstance(viability_thr,float):
        raise TypeError('An integer or a float expected for viability_thr, a/an {} provided instead'.format(type(viability_thr))) 
    elif viability_thr < 0 or viability_thr > 100:
        raise ValueError('viability_thr must be between 0 and 100 as it is a percentage')

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

    # List of essential reactions 
    essential_rxns = [] 

    #--- Minimize rxn fluxes ---
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.biomass_reaction.objective_coefficient = 1

    # First find max biomass for the wild-type strain
    fba_model = fba(model = model, optimization_solver = optimization_solver, build_new_optModel = True, maximize = True, store_fva_flux_bounds = False, simulation_conditions = simulation_conditions, warnings = warnings, stdout_msgs = stdout_msgs)
    fba_model.run()
    if fba_model.solution['exit_flag'] == 'globallyOptimal':
        max_biomass = fba_model.solution['objective_value']
    else:
        raise userError('Infeasbie fba problem for the wild-tyeo strain')

    #---- Turn on the warm start ----
    fba_model.build_new_optModel = False

    fba_model.stdout_msgs = False   

    # Turn preprocessing off (this actually increases the run tiime)
    #if optimization_solver == 'gurobi':
    #    fba_model.optSolver.options['Presolve'] = 0  
    #elif optimization_solver == 'cplex':
    #    fba_model.optSolver.options['PreInd'] = 0  

    counter = 0
    for rxn in model.reactions:
        counter += 1
        rxn.orig_flux_bounds = rxn.flux_bounds

        # Eliminate reaction
        fba_model.optModel.v[rxn.id].setlb(0)
        fba_model.optModel.v[rxn.id].setub(0)

        fba_model.run()
        if (fba_model.solution['exit_flag'] == 'globallyOptimal' and fba_model.solution['objective_value'] < (viability_thr/100)*max_biomass) or fba_model.solution['exit_flag'].lower() == 'infeasible':
            essential_rxns.append(rxn.id)

        fba_model.optModel.v[rxn.id].setlb(rxn.flux_bounds[0])
        fba_model.optModel.v[rxn.id].setub(rxn.flux_bounds[1])


    # Store results into the model
    if save_to_model:
        for rxn in model.reactions:
            if rxn.id in essential_rxns:
                rxn.essential = True 
            else:
                rxn.essential = False

    # Save results into a file 
    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('essential_rxns = [\n')
            for rxn in essential_rxns: 
                f.write("'{}',\n".format(rxn))
            f.write(']\n')

    return essential_rxns        

