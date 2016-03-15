from __future__ import division
import sys, time
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from coopr.environ import *
from tools.pyomoSolverCreator import *
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction

# The following lines change the temporary directory for pyomo
from pyutilib.services import TempfileManager
TempfileManager.tempdir = '/tmp/'

class fbaTools(object):
    """
    A general class for performing various types of FBA methods (FBA, FVA, MOMA, etc) 

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 12-14-2015 
    """   

    def __init__(self,model, optimization_solver = 'gurobi', build_new_optModel = True, store_opt_fluxes = True, flux_key = None, simulation_conditions = None, stdout_msgs = True, warnings = True,  **additional_args): 
        """
        INPUTS (required):
        ------
                      model: An instance of class model containing the information
                             about the metabolic model

        INPUTS (optional):
        ------
        optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                             allowable choices are cplex and gurobi
                stdout_msgs: By default (on) writes  a summary including the solve 
                             status, optimality status (if not optimal), objective 
                             function value and the elapsed time on the screen.
                             if set to a value of 'off' no resuults are written on 
                             the screen, in which case The user can instead specifiy 
                             an output fiile using the option outputFile, or store 
                             them in the variable solution (see the 'run' method for
                             details)
         build_new_optModel: A parameter indicating whether a new pyomo optimizaiton model should be 
                             created (True) or or an existing one should be used (False). The options is useful
                             for the cases a model is already created and one just 
                             wants to change some model attributes (e.g., flux bounds)
                             and rerun FBA. Setting this parameter to False will save 
                             some runtime as the model need not to be created again.
           store_opt_fluxes: If True, it stores the optimal reaction flux value for any reaction  
                             with store_flux parameter set to True.  
                   flux_key: Optimal reaction fluxes after performing FBA are saved into
                             reaction.flux where reaction is the reaction object in the 
                             input model. If flux key is provided, then reaction fluxes 
                             are stored reaction.flux, but reaction.flux in this case is
                             a dictionary instead of a scalar and the current flux value
                             is stored in the dictionary with the key specified by flux_key.
                             This is useful, for example, when performing dynamic FBA,
                             where fluxes should be stored for each time, e.g.,
                             reaction.flux = {0:0.25,0.5:0.24,...}, where keys are tiime points
                             and values are fluxes
        simulation_conditions: A string describing simulation conditions
                   warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                             screen or not. The default is True  
                stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                             Eligible values are True and False and the default is True 
            additional_args: Additional arguments should be entered as normal but they are 
                             converted to a dictionary whose keys are the names of the arguments and 
                             values are the values of  those arguments

        OUTPUTS:
        ---------
        solution: A dictionary with the following keys:
                  exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                            or what is stored in OptSoln.solver.termination_condition
                  objValue: Optimal objective funtion value

        Optimal flux values are stored directly into the 'flux' field of the reaction
        objects. The value of the flux will be None if the problem is not solved to 
        optimality. 

        These are the outputs of the method 'run'
        """
       
        # Metabolic model
        self.model = model

        # Solver name
        self.optimization_solver = optimization_solver

        # Whether to create a pyomo model
        self.build_new_optModel = build_new_optModel
               
        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.warnings = warnings

        # flux_key
        self.flux_key = flux_key 

        # store_opt_fluxes
        self.store_opt_fluxes = store_opt_fluxes

        # Make sure that all reactions in the model have store_flux assigned
        if self.store_opt_fluxes:
            no_store_flux_rxns =  [r.id for r in self.model.reactions if not hasattr(r,'store_flux')]
            if len(no_store_flux_rxns) > 0: 
                if len(no_store_flux_rxns) > 10:
                    raise userError("'store_flux' has not been assigned for the following reactions: " + str(no_store_flux_rxns[:10]) + " and " + str(len(no_store_flux_rxns) - 10) + " other reactions")
                else:
                    raise userError("'store_flux' has not been assigned for the following reactions: " + str(no_store_flux_rxns[:10])) 

        # Simulation conditions
        self.simulation_conditions = simulation_conditions

        # Additoinal arguments. Additional arguments should be entered as normal but they are 
        # converted to a dictionary whose keys are the names of the arguments and values are 
        # the values of  those arguments
        argnames = additional_args.keys()
        argvals = additional_args.values()
        for argname in argnames:
           exec "self." + argname + " = " +"additional_args['" + argname + "']"

        # Check if compounds and reactions in the modle have unique ids
        rxn_ids_num = len(self.model.reactions_by_id.keys())
        uniq_rxn_ids_num = len(set(self.model.reactions_by_id.keys())) 
        if uniq_rxn_ids_num != rxn_ids_num: 
            raise userError('There are reactions with non unique ids in the model. (# of unique rxn ids = ' + str(uniq_rxn_ids_num) + ' , # of rxn ids = ' + str(rxn_ids_num))

        cpd_ids_num = len(self.model.compounds_by_id.keys())
        uniq_cpd_ids_num = len(set(self.model.compounds_by_id.keys())) 
        if uniq_cpd_ids_num != cpd_ids_num: 
            raise userError('There are compounds with non unique ids in the model. (# of unique cpd ids = ' + str(uniq_cpd_ids_num) + ' , # of cpd ids = ' + str(cpd_ids_num))

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
         attr_name: Attribute name
        attr_value: Attribute value
        """
        # model 
        if attr_name == 'model' and not isinstance(attr_value,model):
            raise TypeError('model must be instance of class model')

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Simulation conditions name
        if attr_name == 'simulation_conditions' and (attr_value is not None and not isinstance(attr_value,str)): 
            raise userError('Invalid simulation_conditions for fba model. simulation_conditions must be a striing')

        # build_new_optModel 
        if attr_name == 'build_new_optModel' and not isinstance(attr_value,bool):
            raise TypeError("'build_new_optModel' must be either True or False")

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("'stdout_msgs' must be either True or False")
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("'warnings' must be either True or False")

        self.__dict__[attr_name] = attr_value

    def objectiveFunc_rule(self,optModel):
        """
        Objective function
        """
        pass

    def assignFluxBounds(self,optModel,j):
        """
        Define the flux bounds
        """
        return self.model.reactions_by_id[j].flux_bounds 
        
    def massBalance_rule(self,optModel,i):
        """
        Mass balance 
        """
        return sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0 
        
    def build_optModel(self):
        """
        This optModel creates a pyomo model for FBA optModel
        """   
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of compounds 
        optModel.I = Set(initialize = self.model.compounds_by_id.keys())   

        # Set of rxns  
        optModel.J = Set(initialize = self.model.reactions_by_id.keys())     

        #--- Define the optModel variables --- 
        optModel.v = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)
        
        #--- Defiine the objective function and constraints ----
         # Objective function
        optModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = maximize)

        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule=self.massBalance_rule)

        self.optModel = optModel 
    
    def run(self):
        """ 
        This method runs the FBA tool. 

        OUTPUT:
        -------
        solution: A dictionary with the following keys:
                        exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                                   or what is stored in OptSoln.solver.termination_condition
                  objective_value: Optimal objective funtion value

        Optimal flux values are stored directly into the 'flux' field of the reaction
        objects. The value of the flux will be None if the problem is not solved to 
        optimality. 
        """

        # Total processing and wall time required to create the pyomo model, solve it and store the results
        start_total_pt = time.clock()
        start_total_wt = time.time()

        #---- Creating and instantiating the optModel ----
        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()

        # Create the pyomo model optModel only if self.build_new_optModel == 1        
        if self.build_new_optModel:
            self.build_optModel()

        # Instantiate the optModel
        self.optModel.preprocess()

        #---- Solve the model ----
        # Create a solver and set the options
        solverType = pyomoSolverCreator(self.optimization_solver)

        elapsed_pyomo_pt = (time.clock() - start_pyomo_pt)
        elapsed_pyomo_wt = (time.time() - start_pyomo_wt)

        #- Solve the optModel (tee=True shows the solver output) -
        try:
            # Time to solve the model
            start_solver_pt = time.clock()
            start_solver_wt = time.time()
            OptSoln = solverType.solve(self.optModel,tee=False)
            solverFlag = 'normal'

        # In the case of an error switch the solver
        except  Exception, e:
            if self.warnings:
                print '**WARNING (fba.py)! {} failed with the following error: \n{} \nAn alternative solver is tried'.format(self.optimization_solver,e)

            if self.optimization_solver.lower() == 'gurobi':
                self.optimization_solver = 'cplex'
            elif self.optimization_solver.lower() == 'cplex':
                self.optimization_solver = 'gurobi'

            # Try solving with the alternative solver
            solverType = pyomoSolverCreator(self.optimization_solver)
            try:
                start_solver_pt = time.clock()
                start_solver_wt = time.time()
                OptSoln = solverType.solve(self.optModel,tee=False)
                solverFlag = 'normal'
            except   Exception, e:
                solverFlag = 'solverError'
                if self.warnings:
                    print '**WARNING (fba.py)! {} failed with the following error: \n{} \nAn alternative solver is tried'.format(self.optimization_solver,e)

        elapsed_solver_pt = (time.clock() - start_solver_pt)
        elapsed_solver_wt = (time.time() - start_solver_wt)

        #----- Print the results in the output ------
        if solverFlag == 'normal' and str(OptSoln.solver.termination_condition).lower() == 'optimal':
        
            exit_flag = 'globallyOptimal'

            # Load the results
            self.optModel.load(OptSoln)
        
            # Value of the objective function
            opt_objValue = self.optModel.objectiveFunc()

            # Optimal value of reaction fluxes
            opt_rxnFluxes = {}
            for j in [r.id for r in self.model.reactions]:
                opt_rxnFluxes[j] = self.optModel.v[j].value

            # Print the results on the screen 
            if self.stdout_msgs:
                print "\nSolver.status = ",OptSoln.solver.termination_condition
                print "Optimality status = ",exit_flag
                print "Objective value = ",opt_objValue

            # Store the optimal flux values in the variable 'flux' of the reaction objects
            if self.store_opt_fluxes:
                for rxn in [r for r in self.model.reactions if r.store_flux]:
                    if self.flux_key is None:
                        rxn.flux = self.optModel.v[rxn.id].value
                    elif self.flux_key is not None and type(rxn.flux) is dict:
                        rxn.flux[self.flux_key] = self.optModel.v[rxn.id].value
                    else:
                        rxn.flux = {}
                        rxn.flux[self.flux_key] = self.optModel.v[rxn.id].value

        # If there was a solver error or if an optimal solution was not returned 
        else:
            if solverFlag == 'solverError':
                exit_flag = solverFlag
            else:
                exit_flag = str(OptSoln.solver.termination_condition)

            opt_objValue = None 

            # Optimal value of reaction fluxes
            opt_rxnFluxes = {}
            for j in [r.id for r in self.model.reactions]:
                opt_rxnFluxes[j] = None 

            if self.stdout_msgs:
                print "\n\n** No optimal solutions found (solution.solver.status = ",OptSoln.Solution.status,", solver.status =",OptSoln.solver.status,", solver.termination_condition = ",OptSoln.solver.termination_condition,")\n"

            if self.store_opt_fluxes:
                for rxn in [r for r in self.model.reactions if r.store_flux == True]:
                    if self.flux_key is None and type(rxn.flux) is not dict:
                        rxn.flux = None 
                    elif self.flux_key is not None and type(rxn.flux) is dict:
                        rxn.flux[self.flux_key] = None 
                    else:
                        rxn.flux = {}
                        rxn.flux[self.flux_key] = None 

        self.solution = {'exit_flag':exit_flag,'objective_value':opt_objValue,'opt_rxnFluxes':opt_rxnFluxes}

        # Time required to perform FBA
        elapsed_total_pt = (time.clock() - start_total_pt)
        elapsed_total_wt = (time.time() - start_total_wt)

        if self.stdout_msgs:
           print 'FBA elapsed (processing/wall) time (sec): pyomo = {:.3f}/{:.3f}  ,  solver = {:.3f}/{:.3f}  ,  Total = {:.3f}/{:.3f}\n'.format(elapsed_pyomo_pt,elapsed_pyomo_wt,elapsed_solver_pt,elapsed_solver_wt,elapsed_total_pt,elapsed_total_wt)

        return self.solution

#----------- Sample implementation -----------------
if __name__ == "__main__":
    pass 
