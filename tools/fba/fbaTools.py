from __future__ import division
import sys, time
sys.path.append('../../')
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from coopr.pyomo import *
from coopr.opt import *
from coopr.environ import *
import six
from tools.pyomoSolverCreator import pyomoSolverCreator
from tools.globalVariables import *
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
    Last updated: 11-15-2016 
    """   
    def __init__(self,model, optimization_solver = default_optim_solver, build_new_optModel = True, maximize = True, store_opt_fluxes = True, flux_key = None, simulation_conditions = '', show_solver_output = False, warmstart = False, stdout_msgs = True, warnings = True,  **additional_args): 
        """
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
         show_solver_output: If True, the solver output is shown
                  warmstart: If True, the solver uses a warm start based on the current values of variables
                             (see pyomo documentation for more details).
                stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                             Eligible values are True and False and the default is True 
            additional_args: Additional arguments should be entered as normal but they are 
                             converted to a dictionary whose keys are the names of the arguments and 
                             values are the values of  those arguments

        OUTPUTS:
        ---------
        solution: A dictionary with the following keys:
                  exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                            or what is stored in optSoln.solver.termination_condition
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
               
        # Whether to maximize the objective function 
        self.maximize = maximize

        # show_solver_output and warstart
        self.show_solver_output = show_solver_output
        self.warmstart = warmstart
               
        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.warnings = warnings

        # flux_key
        self.flux_key = flux_key 

        # store_opt_fluxes
        self.store_opt_fluxes = store_opt_fluxes

        # Make sure that all reactions in the model have store_flux assigned
        if self.store_opt_fluxes:
            no_store_flux_rxns = [r.id for r in self.model.reactions if not hasattr(r,'store_flux')]
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

        # Create a solver and set the options
        self.optSolver = pyomoSolverCreator(self.optimization_solver)

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
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi','gurobi_ampl','cplexamp']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Simulation conditions name
        if attr_name == 'simulation_conditions' and not isinstance(attr_value,str): 
            raise userError('simulation_conditions must be a striing')

        # build_new_optModel 
        if attr_name == 'build_new_optModel' and not isinstance(attr_value,bool):
            raise TypeError("'build_new_optModel' must be either True or False")

        # maximize 
        if attr_name == 'maximize' and not isinstance(attr_value,bool):
            raise TypeError("'maximize' must be either True or False")

        # show_solver_output and warmstart
        if attr_name == 'show_solver_output' and not isinstance(attr_value,bool):
            raise TypeError("'show_solver_output' must be either True or False")
        if attr_name == 'warmstart' and not isinstance(attr_value,bool):
            raise TypeError("'warmstart' must be either True or False")

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

    def build_optModel(self):
        """
        This optModel creates a pyomo model for FBA optModel
        """   
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])   

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])     

        #--- Define the optModel variables --- 
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds)
        
        #--- Defiine the objective function and constraints ----
        # Objective function
        if self.maximize:
            optModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = maximize)
        else:
            optModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = minimize)

        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)

        self.optModel = optModel 
    
    def run(self):
        """ 
        This method runs the FBA tool. 

        OUTPUT:
        -------
        solution: A dictionary with the following keys:
                        exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                                   or what is stored in optSoln.solver.termination_condition
                  objective_value: Optimal objective funtion value

        Optimal flux values are stored directly into the 'flux' field of the reaction
        objects. The value of the flux will be None if the problem is not solved to 
        optimality. 
        """
        # Create a solver and set the options
        self.optSolver = pyomoSolverCreator(self.optimization_solver)

        # Total processing and wall time required to create the pyomo model, solve it and store the results
        start_total_pt = time.clock()
        start_total_wt = time.time()

        #---- Creating and instantiating the optModel ----
        start_create_optModel_pt = time.clock()
        start_create_optModel_wt = time.time()

        # Create the pyomo model optModel only if self.build_new_optModel == 1        
        if self.build_new_optModel:
            self.build_optModel()
            elapsed_create_optModel_pt = str(timedelta(seconds = time.clock() - start_create_optModel_pt))
            elapsed_create_optModel_wt = str(timedelta(seconds = time.time() - start_create_optModel_wt))

            if self.warmstart:
                self.optimization_solver == 'gurobi_ampl'
                self.optSolver = pyomoSolverCreator(self.optimization_solver)

                #self.optSolver.options['presolve'] = 0  # Turn off presolve
                self.optSolver.options['basis'] = 3     # make sure gurobi_ampl returns and uses the sstatus suffix (solution basis info)
                self.optSolver.options['method'] = 0    # use primal simplex
                # According to the Gurobi documentation, sstatus suffix
                # values can be interpreted as:
                # 1: basic,   2: superbasic, 3: nonbasic <= (normally =) lower bound, 4: nonbasic >= (normally =) upper bound
                # 5: nonbasic at equal lower and upper bounds,  6: nonbasic between bounds
                # Ref: https://software.sandia.gov/trac/pyomo/browser/pyomo/trunk/examples/pyomo/suffixes/gurobi_ampl_basis.py?rev=11265&order=name
                self.optModel.sstatus = Suffix(direction=Suffix.IMPORT_EXPORT,datatype=Suffix.INT)
                self.optModel.dual = Suffix(direction=Suffix.IMPORT_EXPORT)

        else:
            elapsed_create_optModel_pt = 0 
            elapsed_create_optModel_wt = 0 

        # Processing and wall time for pyomo
        start_preproc_pyomo_pt = time.clock()
        start_preproc_pyomo_wt = time.time()

        # Instantiate the optModel
        self.optModel.preprocess()

        elapsed_preproc_pyomo_pt = str(timedelta(seconds = time.clock() - start_preproc_pyomo_pt))
        elapsed_preproc_pyomo_wt = str(timedelta(seconds = time.time() - start_preproc_pyomo_wt))

        #- Solve the optModel (tee=True shows the solver output) -
        start_solver_pt = time.clock()
        start_solver_wt = time.time()
        try:

            if self.optimization_solver not in ['gurobi_ampl','cplexamp']:
                optSoln = self.optSolver.solve(self.optModel,tee = self.show_solver_output, warmstart = self.warmstart)
            else:
                optSoln = self.optSolver.solve(self.optModel,tee = self.show_solver_output)
            solver_flag = 'normal'

        except  Exception, e:
            if self.warnings:
                print '**WARNING (fba.py)! {} failed with the following error: \n{}'.format(self.optimization_solver,e)
                solver_flag = 'solverError'
      
        elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
        elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))

        #----- Print the results in the output ------
        if solver_flag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':
        
            exit_flag = 'globallyOptimal'

            # Load the results
            self.optModel.load(optSoln)
        
            # Value of the objective function
            opt_objValue = self.optModel.objectiveFunc()

            # Optimal value of reaction fluxes
            opt_rxnFluxes = {}
            for j in [r.id for r in self.model.reactions]:
                opt_rxnFluxes[j] = self.optModel.v[j].value

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
            if solver_flag == 'solverError':
                exit_flag = solver_flag
            else:
                exit_flag = str(optSoln.solver.termination_condition)

            opt_objValue = None 

            # Optimal value of reaction fluxes
            opt_rxnFluxes = {}
            for j in [r.id for r in self.model.reactions]:
                opt_rxnFluxes[j] = None 

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
        elapsed_total_pt = str(timedelta(seconds = time.clock() - start_total_pt))
        elapsed_total_wt = str(timedelta(seconds = time.time() - start_total_wt))

        # Print the results on the screen 
        if self.stdout_msgs:
            if exit_flag == 'solverError':
                print '\nObjective value = None, Optimality status = None, Solution status = None, Solver run status = solverError'
            else:
                print '\nObjective value = {}, Optimality status = {}, Solution status = {}, Solver run status = {}'.format(opt_objValue, optSoln.solver.termination_condition, optSoln.Solution.status, solver_flag)
            print 'Took (hh:mm:ss) {}/{} of processing/walltime in total and {}/{} to solve the model\n'.format(elapsed_total_pt, elapsed_total_wt, elapsed_solver_pt, elapsed_solver_wt)

        return self.solution

#----------- Sample implementation -----------------
if __name__ == "__main__":
    pass 
