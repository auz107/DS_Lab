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

class auxotrophy_finder(object):
    """
    Identifies the set of exchange reactions an auxotrophic mutant strain needs to 
    survive 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06/07/2016 
    """   
    def __init__(self, model, mutant_name = '', max_biomass = None, viability_thr = 0.01, max_soln_size = 0, max_soln_num = 10, max_cpds_uptake_flux = 10, forbidden_cpds_exchrxn_ids = [], optimization_solver = 'gurobi', build_new_optModel = True, validate_solns = True, results_filename = '', warnings = True, stdout_msgs = True): 

        """
        INPUTS:
        ------
        model: 
        An instance of class model containing the information about the metabolic model

        mutant_name:
        Name of the mutant under consideration

        max_biomass: 
        Maximum theoretical biomass under the given growth condition

        viability_thr: 
        Viability threshold below which a cell assumed to be inviable (should be given as a fracton (i.e., a 
        number between zero and 1).

        max_soln_size: 
        Maximum size of the metabolite sets beyond which iterations for finding alternate sets of  
        metabolites that can rescue the mutants stop. The actual criterion that is used is 
        max(min_y,max_soln_size), where min_y is the minimal number of metabolites to rescue the mutant, 
        which is found in the first iteration. Therefore, the default value of zero for 
        max_soln_size means that min_y is used instead. 

        max_soln_num:
        Maximum number of solutions to find

        max_cpds_uptake_flux:
        Maximum allowed uptake flux for compounds in the medium that the mutant strain is 
        auxotrophic for. This typically should be set to the same uptake flux as the limiting
        carbon source. A positive number must be provided. 

        forbidden_cpds_exchrxn_ids:
        List of the exchange reaction ids for compounds that are not allowed to be taken up

        build_new_optModel:
        If True, a new pyomo optimization model is created. Otherwise, an existing one is used

        validate_solns:
        If True, each obtained solution is validated

        optimization_solver: 
        Name of the optimization solver to be used to solve the MILP. 
        Current allowable choices are cplex and gurobi

        results_filename:
        A string containing the filenmae to write the resutls in

        Ali R. Zomorrodi - Daniel Segre Lab @ Boston University
        Last updated: 06/07/2016 
        """
        self.model = model
        self.max_biomass = max_biomass
        self.viability_thr = viability_thr
        self.max_soln_size = max_soln_size
        self.max_soln_num = max_soln_num
        self.max_cpds_uptake_flux = max_cpds_uptake_flux
        self.forbidden_cpds_exchrxn_ids = forbidden_cpds_exchrxn_ids
        self.optimization_solver = optimization_solver 
        self.build_new_optModel = build_new_optModel
        self.validate_solns = validate_solns

        self.mutant_name = mutant_name
        if self.mutant_name == '':
            self.mutant_name = 'this mutant'

        self.warnings = warnings
        self.stdout_msgs = stdout_msgs

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

        # viability_thr 
        if attr_name == 'viability_thr' and not isinstance(attr_value,float) and not isinstance(attr_value,int):
            raise TypeError("'viability_thr' must be a float or integer")
        elif attr_name == 'viability_thr' and (attr_value < 0 or attr_value > 1):
            raise TypeError("'viability_thr' must be a float or integer between zero and one")

        # max_soln_size and max_soln_num 
        if attr_name in ['max_soln_size','max_soln_num'] and not isinstance(attr_value,int):
            raise TypeError("'{}' must be an integer".format(attr_name))
        elif attr_name in ['max_soln_size','max_soln_num'] and attr_value < 0:
            raise TypeError("'{}' must be a non-negative integer between zero and one".format(attr_name))

        # max_biomass 
        if attr_name == 'max_biomass' and not isinstance(attr_value,float) and not isinstance(attr_value,int):
            raise TypeError("'max_biomass' must be a float or integer")
        elif attr_name == 'max_biomass' and attr_value < 0:
            raise TypeError("'max_biomass' must be a non-negative float or integer")

        # max_cpds_uptake_flux 
        if attr_name == 'max_cpds_uptake_flux' and not isinstance(attr_value,int) and not isinstance(attr_value,float):
            raise TypeError("'max_cpds_uptake_flux' must be an integer")
        elif attr_name in 'max_cpds_uptake_flux' and attr_value <= 0:
            raise TypeError("'{}' must be a positive value".format(attr_name))

        # forbidden_cpds_exchrxn_ids 
        if attr_name == 'forbidden_cpds_exchrxn_ids' and not isinstance(attr_value,list):
            raise TypeError("'forbidden_cpds_exchrxn_ids' must be a list")
        elif attr_name == 'forbidden_cpds_exchrxn_ids' and len([s for s in attr_value if not isinstance(s,str)]) > 0:
            raise TypeError("'forbidden_cpds_exchrxn_ids' must be a strings. Objects of type {} were observed in the list.".format(list(set([type(s) for s in attr_value if not isinstance(s,str)]))))

        # build_new_optModel 
        if attr_name == 'build_new_optModel' and not isinstance(attr_value,bool):
            raise TypeError("'build_new_optModel' must be either True or False")

        # build_new_optModel, validate, warnings and messages in the standard output
        if attr_name in ['build_new_optModel','validate_solns', 'warnings', 'stdout_msgs'] and not isinstance(attr_value,bool):
            raise TypeError("'{}' must be either True or False".format(attr_name))

        self.__dict__[attr_name] = attr_value
               

    def build_optModel(self):
        """
        Creates a pyomo model 
        """   
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Define the optModel variables --- 
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds if j not in self._considered_exchrxns else [-self.max_cpds_uptake_flux,1000])

        optModel.y = Var([j for j in self._considered_exchrxns], domain=Boolean)
        
        #--- Defiine the objective function and constraints ----
        # Objective function
        # Minimize sum of the binary variables for exchange rxns corresponding
        # to metabolites that are not in the growth medium
        optModel.objectiveFunc = Objective(rule = lambda optModel: sum([optModel.y[j] for j in self._considered_exchrxns]), sense = minimize)
        
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)

        # Bound constraints for exchange reactions
        optModel.exchangeLB_const = Constraint([j for j in self._considered_exchrxns], rule = lambda optModel, j: optModel.v[j] >= (-self.max_cpds_uptake_flux)*optModel.y[j])

        # Constraint on minimum biomass yield
        optModel.minBiomass_const = Constraint(rule = lambda optModel: 10000*optModel.v[self.model.biomass_reaction.id] >= 10000*self.viability_thr*self.max_biomass)
           
        # Integer cuts
        optModel.integer_cuts = ConstraintList(noruleinit=True)

        # Assign the model to a global variable
        self.optModel = optModel 

    def test(self):
        print '\n------ test ----------\n'
        optModel = ConcreteModel()
        
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Define the optModel variables --- 
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds if j not in self._considered_exchrxns else [-self.max_cpds_uptake_flux,1000])
        #optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds) 
        optModel.y = Var([j for j in self._considered_exchrxns], domain=Boolean)

        #--- Defiine the objective function and constraints ----
        #optModel.objectiveFunc = Objective(rule = lambda optModel: optModel.v[self.model.biomass_reaction.id], sense = maximize)
        def objrule(optModel):
            return sum([optModel.y[j] for j in self._considered_exchrxns])
        #optModel.objectiveFunc = Objective(rule = lambda optModel: sum([optModel.y[j] for j in self._considered_exchrxns]), sense = minimize)
        optModel.objectiveFunc = Objective(rule = objrule, sense = minimize)
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)
        optModel.exchangeLB_const = Constraint([j for j in self._considered_exchrxns], rule = lambda optModel, j: optModel.v[j] >= -self.max_cpds_uptake_flux*optModel.y[j])
        optModel.minBiomass_const = Constraint(rule = lambda optModel: optModel.v[self.model.biomass_reaction.id] >= self.viability_thr*self.max_biomass)

        #for j in self._considered_exchrxns[:240]: 
        #   optModel.y[j] = 0
        #   optModel.y[j].fixed = True

        optSolver = pyomoSolverCreator('gurobi')
        optModel.preprocess()
        try:
            optSoln = optSolver.solve(optModel,tee = True)
            solver_flag = 'normal'
        except Exception, e:
            print '**WARNING (fba.py)! {} failed with the following error: \n{}'.format(self.optimization_solver,e)
            solver_flag = 'solverError'

        if solver_flag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':
            optModel.load(optSoln)
            opt_objValue = optModel.objectiveFunc()
            exit_flag = 'globallyOptimal'
            soln = [jj for jj in self._considered_exchrxns if abs(optModel.y[jj].value - 1) <= 1e-6]
        else:
            if solver_flag == 'solverError':
                exit_flag = solver_flag
            else:
                exit_flag = str(optSoln.solver.termination_condition)

            opt_objValue = None

        print 'solution = ',soln
        if exit_flag == 'solverError':
            print '\nObjective value = None, Optimality status = None, Solution status = None, Solver run status = solverError\n'
        else:
            print '\nObjective value = {}, Optimality status = {}, Solution status = {}, Solver run status = {}\n'.format(opt_objValue, optSoln.solver.termination_condition, optSoln.Solution.status, solver_flag)

    def validate(self, what_to_validate = 'optimal_soln'): 
        """
        Validates an obtained solution

        what_to_validate is a string determining what needs to be validated: 
            'optimal_soln': 
            an optimal solution (a set of exchange reactions for which
            the binary variables are one), or

            'already_viable_mutant': 
            An optimal solution with a zero objective funciton value meaning that the mutant
            is already viable and all binary variables are zero             
        """
        if not isinstance(what_to_validate,str):
            raise TypeError('what_to_validate must be a string')
        elif what_to_validate not in ['optimal_soln','already_viable_mutant']:
            raise ValueError('Invalid value for what_to_validate! Alowed values are optimal_soln and already_viable_mutant')

        if what_to_validate.lower() == 'optimal_soln':
            # First save the flux bounds for exchange reactions
            exchrxns_flux_bounds = dict([(r.id,r.flux_bounds) for r in self.model.reactions if r.is_exchange])
    
            # Set the flux bounds for compounds that can be taken up to -max_cpds_uptake_flux
            for exchrxn in [r for r in self.model.reactions if r.is_exchange and r.id in self._curr_soln]:
                exchrxn.flux_bounds[0] = - self.max_cpds_uptake_flux
            self.model.fba(stdout_msgs = False)
            if self.model.fba_model.solution['exit_flag'] != 'globallyOptimal':
                self._validated = False
                raise userError('Validation failed for solution {} for {} because of an infeasible fba problem'.format(self._curr_soln, self.mutant_name))
            elif self.model.fba_model.solution['exit_flag'] == 'globallyOptimal' and self.model.fba_model.solution['objective_value'] < self.viability_thr*self.max_biomass:
                self._validated = False
                raise userError('Validation failed for solution {} for {} because max biomass with the uptake of these compounds ({}) is than viability_thr*max_biomass ({})'.format(self._curr_soln, self.mutant_name, self.model.fba_model.solution['objective_value'], self.viability_thr*self.max_biomass))
            else:
                self._validated = True

                # Validated. Restore the flux bounds for exchange reactions
                for exchrxn_id in exchrxns_flux_bounds.keys():
                    self.model.reactions_by_id[exchrxn_id].flux_bounds = exchrxns_flux_bounds[exchrxn_id]           
      
                if self.stdout_msgs:
                    print '\tSolution validated ...'

        elif what_to_validate.lower() == 'already_viable_mutant':
            self.model.fba(stdout_msgs = False)
            if self.model.fba_model.solution['exit_flag'] != 'globallyOptimal':
                self._validated = False
                print '**WARNING! Validation of already_viable_mutant failed for {} because of an infeasible fba problem'.format(self.mutant_name)
            elif self.model.fba_model.solution['exit_flag'] == 'globallyOptimal' and self.model.fba_model.solution['objective_value'] < self.viability_thr*self.max_biomass:
                self._validated = False
                print '**WARNING! Validation of already_viable_mutant failed for {} because max biomass without the uptake of any compounds ({}) is than viability_thr*max_biomass ({})'.format(self.mutant_name, self.model.fba_model.solution['objective_value'], self.viability_thr*self.max_biomass)

            else:
                self._validated = True
                if self.warnings:
                    print '**WARNING! {} is viable in the current specified medium'.format(self.mutant_name)


    def run(self):
        """ 
        This method runs the optimization problem to find out the set of all exchange reactions for
        compounds that need to be taken up in order for he mutant strain to survive. 

        OUTPUS:
        ------
        exit_flag: 
        A striing showing which one of the termination conditions for the function run was satisfied. 
        It can take the following values:
            already_viable_mutant: 
            If an optimal solution is found for the optimization problem, however, its objective 
            function value turns out to be zero, implying that no metabolites are needed to be taken up

            obj_value_greater_than_max_soln_size: 
            If an optimal solution is found for the optimization problem, however, its objective 
            value is greater than max(min_y,max_soln_size) where min_y is the minimum objective 
            value found in the first run
           
            infeasible_opt_problem: 
            The optimization problem was not solved to optimality
 
        auxotroph_cpds: 
        A iist of lists, where each list contains the set of exchange reactions for compounds 
        needed for survivna 
        """
        start_total_pt = time.clock()
        start_total_wt = time.time()

        #---- Creating and instantiating the optModel ----
        # Set of exchange reacitons for which binary variables are defined
        self._considered_exchrxns = [j.id for j in self.model.reactions if j.is_exchange and j.flux_bounds[0] >= 0 and j.id not in self.forbidden_cpds_exchrxn_ids]

        #self.test()

        # Create the abstract optModel         
        # Create the pyomo model optModel only if self.build_new_optModel == 1        
        if self.build_new_optModel:
            self.build_optModel()

        # Instantiate the optModel
        self.optModel.preprocess()

        # Create a solver and set the options
        self.optSolver = pyomoSolverCreator(self.optimization_solver)

        #- Solve the optModel (tee=True shows the solver output) -
        # A parameter showing to tell python when to stop running the while loop
        done = False

        #-- Some initializations --
        auxotroph_cpds = []

        # Run counter
        counter = 0

        exit_flag = ''

        while not done:
            counter += 1

            if self.stdout_msgs:
                print '\n------- run #',counter+1,' ----------\n'

            # Add integer cuts 
            if counter > 1:
                # Add a new integer cut to exlucde the previously found solution 
                self.optModel.integer_cuts.add(sum([self.optModel.y[j] for j in self.prev_y.keys() if self.prev_y[j] == 1]) <= self.min_y - 1)

            # Instantiate the optModel with new fixed variables
            self.optModel.preprocess()

            #- Solve the optModel (tee=True shows the solver output) -
            start_solver_pt = time.clock()
            start_solver_wt = time.time()
    
            #self.optSolver.options['logfile']= 'my_gurobi.log'
            try:
                optSoln = self.optSolver.solve(self.optModel, tee = False)
                solver_flag = 'normal'
            except: 
                if self.warnings:
                    print '**WARNING (fba.py)! {} failed with the following error: \n{}'.format(self.optimization_solver,e)
                solver_flag = 'solverError'

            elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
            elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))
    
            if solver_flag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal': 
            
                self.exit_flag_curr = 'globallyOptimal'
    
                # Load the results
                self.optModel.load(optSoln)
            
                # Value of the objective function
                opt_objValue = self.optModel.objectiveFunc()

                # Minimum number of metabolites required to be taken up for survival
                # (identified in the first run)
                if counter == 1:
                    # Use Decimal for the equality testing to avoid problems with floadting points
                    if opt_objValue <= mip_integrality_tol:
                        self.validate(what_to_validate = 'already_viable_mutant')
                        done = True
                        if self._validated:
                            exit_flag = 'already_viable_mutant'
                        else:
                            exit_flag = 'Validation failed for already_viable_mutant'
                    else:
                        self.min_y = opt_objValue
  
                # Stop if the value of objective function in the current run is greater
                # than self.min_y
                if counter > 1 and opt_objValue > max(self.min_y,self.max_soln_size):
                    done = True
                    exit_flag = 'obj_value_greater_than_max_soln_size'

                # Write the results 
                if not done:
                    self._curr_soln = [j for j in self._considered_exchrxns if abs(self.optModel.y[j].value - 1) <= mip_integrality_tol]

                    # First validate the soluition before storing it
                    self.validate(what_to_validate = 'optimal_soln')

                    auxotroph_cpds.append(self._curr_soln)


                    if counter > 0:
                        # Previously found solution
                        self.prev_y = {}
                        for j in self._considered_exchrxns:
                            if j in self._curr_soln:
                                self.prev_y[j] = 1
                            else:
                                self.prev_y[j] = 0
           
                if counter == self.max_soln_num:
                    done = True
                    exit_flag = 'max allowed number of solutions reached'


            # If the optimization problem was not solved successfully
            else:
                # Stop running the while loop
                done = True
                exit_flag = 'infeasible_opt_problem'
                self.exit_flag_curr = str(optSoln.solver.termination_condition) 

        # Time required to run 
        elapsed_total_pt = str(timedelta(seconds = time.clock() - start_total_pt))
        elapsed_total_wt = str(timedelta(seconds = time.time() - start_total_wt))
    
        # Print the results on the screen 
        if self.stdout_msgs:
            if exit_flag == 'solverError':
                print '\nObjective value = None, Optimality status = None, Solution status = None, Solver run status = solverError'
            else:
                print '\nObjective value = {}, Optimality status = {}, Solution status = {}, Solver run status = {}'.format(opt_objValue, optSoln.
solver.termination_condition, optSoln.Solution.status, solver_flag)
            print 'Took (hh:mm:ss) {}/{} of processing/walltime in total and {}/{} to solve the model\n'.format(elapsed_total_pt, elapsed_total_wt, elapsed_solver_pt, elapsed_solver_wt)
 
        return (exit_flag, auxotroph_cpds)

#-------------------- Sample implementation -----------------
if __name__ == "__main__":
    pass
