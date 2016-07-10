from __future__ import division
import sys, time, copy
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.globalVariables import *
from tools.pyomoSolverCreator import *
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.utilities.get_ModelSEED_ids import get_ModelSEED_ids

class gapfill(object):
    """
    Performs gap filling for a metabolic model  

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 04-22-2016 
    """   
    def __init__(self, model, super_model = None, gapfilling_method = 'parsimony_based', penalties = {}, viability_thr = 0.01, max_soln_num = 5, 
                 growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, standard_to_model_compartID_map = {'c':'c','e':'e','p':'p'}, 
                 fixed_external_rxns = {}, validate_results = True, results_filename = '',
                 optimization_solver = default_optim_solver, build_new_optModel = True, obtain_ModelSEED_ids = False, warnings = True, stdout_msgs = True, stdout_msgs_details = True, **additional_args): 
        """
        INPUTS (required):
        ------
        model: 
        Original metabolic model to be gap filled. An instance of class model 

        super_model: 
        The super model integrating the original model and external reactions form ModelSEED 
        and  exchange and transport reactoins in the model. If no input is provided, super_model
        is created

        standard_to_model_compartID_map: 
        This input is needed only if super_model = None. It is a dictionary where keys are 
        one of the letters below (standard compartment ids and values are the corresponding 
        compartment ids in the model id in the original model. 
        c: Cytosol (cytoplasm), e: Extracellular, g: Golgi, m: Mitochondria, n: Nucleus
        p: Periplasm, r: Endoplasmic reticulum, x: Peroxisome
        For example, if a model has two compartments c and e, one can provide 
        {'c':'c id in the model', 'e':'e id in the model'}. One can also provide
        {'c':'', 'e': ''} in which the code searches for these two compartments 
        in the model

        gapfilling_method: 
        The method to perform gap filling. The default is parsimony-based gap filling, where the
        aim is to minimally perturb the model.

        penalties: 
        Penalties associated with performing each modification to the model

        viability_thr: 
        Viability threshold, i.e., the min biomass flux above which we assuem growth

        max_soln_num: 
        Maximum number of solution to obtain

        obtain_ModelSEED_ids: 
        If True ModelSEED ids of reactions and compounds are obtained

        fixed_external_rxns: 
        External reactions that have to be fixed (added) to the model. This is a dictionary with
        keys being reaction ids and values being either 'forward', 'backward' or 'reversible'
        meaning to add the reaction in the forward direction, backward direction or as a 
        reversible reaction (both directions)

        validate_results: 
        If True, each obtained results is validated

        results_filename: 
        A string containing the path and name of the file to store the results

        optimization_solver: 
        Name of the LP solver to be used to solve the LP. Current 
        allowable choices are cplex and gurobi

        build_new_optModel: 
        If True a new pyomo optimization model is constructed

        stdout_msgs: 
        If True writes a summary of the run in the output 


        OUTPUTS:
        ---------
        solution: 
        A dictionary with the following keys:
            exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                       or what is stored in OptSoln.solver.termination_condition
             objValue: Optimal objective funtion value

        These are the outputs of the method 'run'
        """
        # Original model
        self.model = model

        # Super model
        self.super_model = super_model

        # Gap filling method
        self.gapfilling_method = gapfilling_method

        # Penalties associated with adding each reaction type
        self.penalties = penalties

        # Obtain ModelSEED ids
        self.obtain_ModelSEED_ids = obtain_ModelSEED_ids

        # growthMedium_flux_bounds
        self.growthMedium_flux_bounds = growthMedium_flux_bounds

        # Viability threshold
        self.viability_thr = viability_thr

        # max number of solutions (gap filling strategies) to obtain
        self.max_soln_num = max_soln_num

        # External reactions that must be fixed (added to the model)
        self.fixed_external_rxns = fixed_external_rxns

        # Validate the obtained solutions
        self.validate_results = validate_results

        # Results file name
        self.results_filename = results_filename
               
        # Solver name
        self.optimization_solver = optimization_solver

        # Whether to create a pyomo model
        self.build_new_optModel = build_new_optModel

        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.stdout_msgs_details = stdout_msgs_details
        self.warnings = warnings

        # Additoinal arguments. Additional arguments should be entered as normal but they are 
        # converted to a dictionary whose keys are the names of the arguments and values are 
        # the values of  those arguments
        argnames = additional_args.keys()
        argvals = additional_args.values()
        for argname in argnames:
           exec "self." + argname + " = " +"additional_args['" + argname + "']"

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        if attr_name.lower() == 'gapfilling_method' and not isinstance(attr_value,str):
            raise TypeError('gapfilling_method must be a string')
        elif attr_name.lower() == 'gapfilling_method' and attr_value.lower() not in ['parsimony_based']: 
            raise ValueError('Invalid gapfilling_method value! Allowed choices are: parsimony_based')

        if attr_name.lower() == 'gapfilling_method' and not isinstance(attr_value,str):
            raise TypeError('gapfilling_method must be a string')
        elif attr_name.lower() == 'gapfilling_method' and attr_value.lower() not in ['parsimony_based']: 
            raise ValueError('Invalid gapfilling_method value! Allowed choices are: parsimony_based')

        if attr_name.lower() == 'penalties' and not isinstance(attr_value,dict):
            raise TypeError('penalties must be a dictionary')
        elif attr_name.lower() == 'penalties' and len(attr_value.keys()) > 0 and len([k for k in attr_value.keys() if v not in ['exchange_rxns']]) > 0:
            raise ValuesError('Invalid key for penalties: {}. Allowed keys are [exchange_rxns]'.format([k for k in attr_value.keys() if v not in ['exchange_rxns']]))

        if attr_name.lower() == 'viability_thr' and not isinstance(attr_value,int) and not isinstance(attr_value,float):
            raise TypeError('viability_thr must be an interger or float')
        elif attr_name.lower() == 'viability_thr' and attr_value < 0:
            raise ValueError('viability_thr must be a non-negative integer or float')

        if attr_name == 'growthMedium_flux_bounds' and not isinstance(attr_value,dict):
            raise TypeError('growthMedium_flux_bounds must be a dictionary')
        if attr_name == 'growthMedium_flux_bounds' and len([k for k in attr_value.keys() if k.lower() not in ['flux_bounds_filename','flux_bounds_dict']]) > 0:
            raise ValueError('Invalid key for growthMedium_flux_bounds. Allowed keys are flux_bounds_filename and flux_bounds_dict')

        if attr_name.lower() == 'max_soln_num' and not isinstance(attr_value,int):
            raise TypeError('max_soln_num must be an integer')
        elif attr_name.lower() == 'max_soln_num' and attr_value < 0: 
            raise TypeError('max_soln_num must be a non-negative integer')

        if attr_name.lower() == 'results_filename' and not isinstance(attr_value,str):
            raise TypeError('results_filename must be a string')

        if attr_name.lower() in ['validate_results','build_new_optModel','obtain_ModelSEED_ids','warnings', 'stdout_msgs', 'stdout_msgs_details'] and not isinstance(attr_value,bool):
            raise TypeError('{} must be either True or False'.format(attr_name))

        self.__dict__[attr_name] = attr_value

    #------------------------------------------------------------------------
    #--- Define parameters needed for the optimization problem ---
    #------------------------------------------------------------------------
    def set_model_flux_bounds(self):
        """
        Set the flux bounds for the growth medium and min biomass flux in the model
        This functions needs to be caleed before creating an optModel
        """
        self.super_model.reset_flux_bounds()

        # Set the lower bound of all irreversible reactions in the original model to -1000 
        # in order to allow them to go in the backward direction as well
        for rxn in [r for r in self.super_model.reactions if not r.external and r.reversibility.lower() == 'irreversible']:
            rxn.flux_bounds[0] = -1000

        # Set the flux bounds for all non-exchange exteranl reactions to [-1000,1000], 
        # irrespective of what their reversbility is to consider adding them both in the 
        # forward and backward directions
        for rxn in [r for r in self.super_model if r.external and (not r.is_exchange or 'exchange' not in r.reversibility)]:
            rxn.flux_bounds = [-1000,1000]        

        # Set the flux bounds for the model. NOTE: Do NOT reset flux bounds here!
        set_specific_bounds(model = self.super_model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = False)

        # Impose constraints on the min biomass formation flux
        set_specific_bounds(model = self.super_model, flux_bounds = {self.model.biomass_reaction.id:[self.viability_thr,None]}, reset_flux_bounds = False)


    def preproc(self):
        """
        Performs a number of preprocessing before the gap filling procedure. These include:
        - get ModelSEED ids for reactions and metabolites in the model
        - Create super_model if it is not provided in the input
        """
        # Get ModelSEED ids
        if self.obtain_ModelSEED_ids:
            get_ModelSEED_ids(model = self.model, msgs = True)

        # Create super_model
        if self.super_model == None:
            self.super_model = create_super_model(original_model = model, standard_to_model_compartID_map = standard_to_model_compartID_map, dd_c_cpds_to_otherComparts = False)

        # Set the flux bounds
        self.set_model_flux_bounds()

    def fix_known_variables(self):
        """
        Fixes a number of known variables
        """
        for rxn_id in self.fixed_external_rxns.keys():
            if self.fixed_external_rxns[rxn_id].lower() == 'forward':
                self.optModel.yf[rxn_id] = 1
                self.optModel.yf[rxn_id] = True
            elif self.fixed_external_rxns[rxn_id].lower() == 'backward':
                self.optModel.yb[rxn_id] = 1
                self.optModel.yb[rxn_id] = True
            elif self.fixed_external_rxns[rxn_id].lower() == 'reversible':
                self.optModel.yf[rxn_id] = 1
                self.optModel.yf[rxn_id] = True
                self.optModel.yb[rxn_id] = 1
                self.optModel.yb[rxn_id] = True
            else:
                raise ValueError('Unknonw value for fixed_external_rxns for rxn {}: {}'.format(rxn_id, fixed_external_rxns[rxn_id]))

    #--------------------------------------------------------
    #---------------- Create optimization models ------------        
    #--------------------------------------------------------
    def build_optModel(self):
        """
        Creates a pyomo optimization model for the primal problem 
        """
        # Create a pyomo model optimization model
        optModel = ConcreteModel()

        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Variables --- 
        # Reaction fluxes (NOTE: Make sure flux bounds are correctly assigned)
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.super_model.reactions_by_id[j].flux_bounds)

        # Binary two variables for reactions in the database, one for adding that reaction in the forward direciton and
        # the other for adding it in the backward direction. The binary variable in the backward direciton is also used to 
        # relax the irreversibility constraint on irreversible reactions in the original model 
        optModel.yf = Var(optModel.J, domain = Boolean)
        optModel.yb = Var(optModel.J, domain = Boolean)

        # A binary variable to check whether an external reaction must be included in the model 
        # (forward or backward or both). This is needed to take of the base cost. This binary variable
        # doesn't have to be decalred as binary as the related constriants would take care of it
        optModel.y = Var(optModel.J, domain = NonNegativeReals, bounds = [0,1])

        #--- Objective function ---
        optModel.objectiveFunc = Objective(rule = lambda optModel: 
                                 sum([self.super_model.reactions_by_id[j].base_cost*optModel.y[j] + self.super_model.reactions_by_id[j].forward_cost*optModel.yf[j] + 
                                      self.super_model.reactions_by_id[j].base_cost*optModel.y[j] + self.super_model.reactions_by_id[j].backward_cost*optModel.yb[j] 
                                      for j in optModel.J if self.super_model.reactions_by_id[j].external] + 
                                 sum([self.super_model.reactions_by_id[j].backward_cost*optModel.yb[j] 
                                      for j in optModel.J if not self.super_model.reactions_by_id[j].external and self.super_model.reactions_by_id[j].reversibility.lower() == 'reversible'])), 
                                 sense = minimize)

        #--- Constraints ----
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)

        # v(j) >= (-1000)*yb(j)). This constraint is written for both reactions in the external database and 
        # irreversible reactions int he original model to allow the relaxation of the irreversibility constraints
        optModel.LB_const = Constraint([j for j in optModel.J if self.super_model.reactions_by_id[j].external or (not self.super_model.reactions_by_id[j].external and self.super_model.reactions_by_id[j].reversibility.lower() == 'irreversible')], rule = lambda optModel, j: optModel.v[j] >= -1000*optModel.yb[j])

        # v(j) <= 1000*yf(j). This constraint is written only for reactions in the external database
        optModel.UB_const = Constraint([j for j in optModel.J if self.super_model.reactions_by_id[j].external], rule = lambda optModel, j: optModel.v[j] <= 1000*optModel.yf[j])

        # Constraints for y(j). In fact, y(j) = yf(j) OR yb(j)
        optModel.const1_y = Constraint(optModel.J, rule = lambda optModel, j: optModel.y[j] >= optModel.yf[j])
        optModel.const2_y = Constraint(optModel.J, rule = lambda optModel, j: optModel.y[j] >= optModel.yb[j])
        optModel.const3_y = Constraint(optModel.J, rule = lambda optModel, j: optModel.y[j] <= optModel.yf[j] + optModel.yb[j])

        # integer cut
        optModel.integer_cuts = ConstraintList(noruleinit = True)

        self.optModel = optModel

    def run_single(self):
        """
        Performs a single run of gap fill
        """
        # Processing and wall time for pyomo
        start_preproc_pyomo_pt = time.clock()
        start_preproc_pyomo_wt = time.time()

        # Instantiate the optModel
        self.optModel.preprocess()

        elapsed_preproc_pyomo_pt = str(timedelta(seconds = time.clock() - start_preproc_pyomo_pt))
        elapsed_preproc_pyomo_wt = str(timedelta(seconds = time.time() - start_preproc_pyomo_wt))

        #---- Solve the model ----
        #- Solve the optModel (tee=True shows the solver output) -
        try:
            # Processing and wall time for the solver
            start_solver_pt = time.clock()
            start_solver_wt = time.time()
            optSoln = self._optSolver.solve(self.optModel,tee=False)
            solver_flag = 'normal'

        # In the case of an error switch the solver
        except  Exception, e:
            solver_flag = 'solverError'
            if self.warnings:
                print '**WARNING (FORCE)! {} failed with the following error: \n{} \n'.format(self.optimization_solver,e)

        elapsed_solver_pt = str(timedelta(seconds = time.clock() - start_solver_pt))
        elapsed_solver_wt = str(timedelta(seconds = time.time() - start_solver_wt))

        if solver_flag == 'normal' and str(optSoln.solver.termination_condition).lower() == 'optimal':

            exit_flag = 'globallyOptimal'

            # Load the results
            self.optModel.load(optSoln)

            # Optimal value of the objective function
            opt_objValue = self.optModel.objectiveFunc()

            # Reactions added in forward or backward directions and reaction included in the model 
            yf_one_rxns = [j for j in self.optModel.J if abs(self.optModel.yf[j].value - 1) <= mip_integrality_tol]
            yb_one_rxns = [j for j in self.optModel.J if abs(self.optModel.yb[j].value - 1) <= mip_integrality_tol]
            y_one_rxns = [j for j in self.optModel.J if abs(self.optModel.y[j].value - 1) <= mip_integrality_tol]

        # If there was a solver error or if an optimal solution was not returned 
        else:
            opt_objValue = None
            if solver_flag == 'solverError':
                exit_flag = solver_flag
            else:
                exit_flag = str(optSoln.solver.termination_condition)

            # Reactions that must be knocked out, up-regulated or down-regulated
            yf_one_rxns = []
            yb_one_rxns = []
            y_one_rxns = []

        # Store the solution
        self._curr_soln = {'exit_flag':exit_flag,'objective_value':opt_objValue,'added_rxns_num': len(y_one_rxns),'forward_rxns':yf_one_rxns,'backward_rxns':yb_one_rxns}

        # Print the results on the screen 
        if self.stdout_msgs_details:
            print '\nObjective value = {}, Optimality status = {}, Solution status = {}, Solver run status = {}'.format(opt_objValue, optSoln.solver.termination_condition, optSoln.Solution.status, solver_flag)
            print 'Took (hh:mm:ss) {}/{} of processing/walltime to create a pyomo model, {}/{} to  preprcoess the model and {}/{} to solve the model\n'.format(self._elapsed_create_optModel_pt, self._elapsed_create_optModel_wt, elapsed_preproc_pyomo_pt,elapsed_preproc_pyomo_wt, elapsed_solver_pt,elapsed_solver_wt)


    def validate_soln(self):
        """
        Validates an obtained solution
        """
        updated_model = copy.deepcopy(super_model)

        forward_rxns = []
        forward_rxns =  [j for j in self.optModel.J if j in self._curr_soln['forward_rxns'] and j not in self._curr_soln['backward_rxns']]
        for rxn_f in forward_rxns: 
             rxn_obj = super_model.reactions_by_id[rxn_f]
             rxn_obj.reversibility = 'irreversible'
             rxn_obj.external = False
             for cpd in [c for c in rxn_obj.compounds if c.external]:
                 cpd.external = False
             forward_rxns.append(rxn_f)


        backward_rxns =  [j for j in self.optModel.J if j not in self._curr_soln['forward_rxns'] and j in self._curr_soln['backward_rxns']]
        for rxn_b in backward_rxns: 
             rxn_obj = super_model.reactions_by_id[rxn_b]
             # Multiply stoichiomteric coefficients by -1
             rxn_obj.set_stoichiometry(stoichiometry = dict([(cpd,-cdp_stoic) for (cpd,cpd_stoic) in rxn_obj.stoichiometry.items()]), replace = True)
             rxn_obj.reversibility = 'irreversible'
             rxn_obj.external = False
             for cpd in [c for c in rxn_obj.compounds if c.external]:
                 cpd.external = False

        reversible_rxns = [j for j in self.optModel.J if j in self._curr_soln['forward_rxns'] and j in self._curr_soln['backward_rxns']]
        for rxn in reversible_rxns: 
             rxn_obj = super_model.reactions_by_id[rxn]
             rxn_obj.reversibility = 'reversible'
             rxn_obj.external = False
             for cpd in [c for c in rxn_obj.compounds if c.external]:
                 cpd.external = False

        # Remove all compounds and reactions for which external = True
        updated_model.del_compounds([c for c in updated_model.compounds if c.external])
        updated_model.del_reactions([r for r in updated_model.reactions if r.external])

        # Now perform FBA and check biomass flux
        for rxn in updated_model.reactions:
            rxn.objective_coefficient = 0
        updated_model.biomass_reaction.objective_coefficient = 1

        # Set the flux bounds for the model. **Do not reset flux bounds here!
        set_specific_bounds(model = self.super_model, file_name = self.growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = self.growthMedium_flux_bounds['flux_bounds_dict'], reset_flux_bounds = True)

        updated_model.fba()
        if updated_model.fba_model.solution['exit_flag'] == 'globallyOptimal' and updated_model.fba_model.solution['objective_value'] >=  self.viability_thr:
            if self.stdout_msgs_details:
                print 'The current solution was successfully validated'
        elif updated_model.fba_model.solution['exit_flag'] == 'globallyOptimal' and updated_model.fba_model.solution['objective_value'] <  self.viability_thr:
            raise userError('The following solution was not validated because the biomass is equal to {}, which is less than the minimum viability threshould of {}:\nForward rxns = {}\nBackward rxns = {}, Reversible reaction s= {}'.format(updated_model.fba_model.solution['objective_value'], self.viability_thr,forward_rxns, backward_rxns, reversible_rxns)) 
        elif updated_model.fba_model.solution['exit_flag'] != 'globallyOptimal':
            raise userError('The following solution was not validated because the fba problem was not solved to optimality: exit_flag = {}:\nForward rxns = {}\nBackward rxns = {}, Reversible reaction s= {}'.format(updated_model.fba_model.solution['exit_flag'],forward_rxns, backward_rxns, reversible_rxns))


    def run(self):
        """
        Runs the gap filling procedure
        """
        # First perform some preprocessings
        self.preproc()

        # Number of solutions found so far
        found_solutions_num = 0

        if self.results_filename != '':
            with open(self.results_filename,'w') as f:
                f.write('gap_filling_solutions = [\n')

        # Creating the pyomo optModel
        if self.build_new_optModel:
            start_pyomo_pt = time.clock()
            start_pyomo_wt = time.time()

            self.build_optModel()

            self._elapsed_create_optModel_pt = str(timedelta(seconds = time.clock() - start_pyomo_pt))
            self._elapsed_create_optModel_wt = str(timedelta(seconds = time.time() - start_pyomo_wt))
        else:
            self._elapsed_create_optModel_pt = 0
            self._elapsed_create_optModel_wt = 0

        # Fix known variables
        self.fix_known_variables()

        # Create a solver and set the options
        self._optSolver = pyomoSolverCreator(self.optimization_solver)

        # A list of dictionaries containing the solutions obtained in different iterations
        self.solutions = []

        done = False

        while not done:

            # Solve the optimizaiton model 
            self.run_single()

            # Check whether the current run was successful
            if self._curr_soln['exit_flag'] == 'globallyOptimal':

                if validate_results:
                    self.validate_soln()

                self.solutions.append(self._curr_soln)

                found_solutions_num += 1

                # Add an integer cut excluding the current solution. 
                # Total number of yf and yb variables is 2*len(self.optModel.J)
                self.optModel.integer_cuts.add(
                     sum([self.optModel.yf[j] for j in self._curr_soln['forward_rxns']]) + \
                     sum([self.optModel.yb[j] for j in self._curr_soln['backward_rxns']]) + \
                     sum([self.optModel.yf[j] for j in self.fixed_to_one_yf_rxns]) + \
                     sum([self.optModel.yb[j] for j in self.fixed_to_one_yb_rxns]) + \
                     sum([1 - self.optModel.yf[j] for j in self.optModel.J if j not in self._curr_soln['forward_rxns']]) + \
                     sum([1 - self.optModel.yb[j] for j in self.optModel.J if j not in self._curr_soln['backward_rxns']]) <= \
                     2*len(self.optModel.J) - 1)

                # Save the results into file 
                if self.results_filename != '':
                    with open(self.results_filename,'a') as f:

                        f.write("\n\t{'objective_value': {}\n".format(self._curr_soln['objective_value']))

                        # Reactions added only in forward direction
                        f.write("\t{'forward_rxns': {\n")
                        for rxn_f in [j for j in self.optModel.J if j in self._curr_soln['forward_rxns'] and j not in self._curr_soln['backward_rxns']]:
                            rxn_obj = super_model.reactions_by_id[rxn_f]
                            f.write("\t\t{}:{{'name':{}, 'equation':{}, 'external_type':{}, 'template_rxn_type':{}, 'base_cost':{}, 'forward_cost':{}, 'backward_cost':{}}},\n".format(rxn_f, rxn_obj.name, rxn_obj.get_equation(show_cpds_by = 'name'), rxn_obj.external_type, rxn_obj.template_rxn_type, rxn_obj.base_cost, rxn_obj.forward_cost, rxn_obj.backward_cost))
                        f.write("\t}\n")

                    # Reactions added only in backward direction
                        f.write("\t{'backward_rxns': {\n")
                        for rxn_b in [j for j in self.optModel.J if j not in self._curr_soln['forward_rxns'] and j in self._curr_soln['backward_rxns']]:
                            rxn_obj = super_model.reactions_by_id[rxn_b]
                            f.write("\t\t{}:{{'name':{}, 'equation':{}, 'external_type':{}, 'template_rxn_type':{}, 'base_cost':{}, 'forward_cost':{}, 'backward_cost':{}}},\n".format(rxn_b, rxn_obj.name, rxn_obj.get_equation(show_cpds_by = 'name'), rxn_obj.external_type, rxn_obj.template_rxn_type, rxn_obj.base_cost, rxn_obj.forward_cost, rxn_obj.backward_cost))
                        f.write("\t}\n")

                        # Reactions added as reversible 
                        f.write("\t{'reversible_rxns': {\n")
                        for rxn in [j for j in self.optModel.J if j in self._curr_soln['forward_rxns'] and j in self._curr_soln['backward_rxns']]:
                            rxn_obj = super_model.reactions_by_id[rxn_f]
                            f.write("\t\t{}:{{'name':{}, 'equation':{}, 'external_type':{}, 'template_rxn_type':{}, 'base_cost':{}, 'forward_cost':{}, 'backward_cost':{}}},\n".format(rxn, rxn_obj.name, rxn_obj.get_equation(show_cpds_by = 'name'), rxn_obj.external_type, rxn_obj.template_rxn_type, rxn_obj.base_cost, rxn_obj.forward_cost, rxn_obj.backward_cost))
                        f.write("\t}\n")
                        f.write("\t}\n")

                # print a summary of results in the output
                if self.stdout_msgs:
                    print '\n{})\n'.format(found_solutions_num)
                    print 'Total # of added reactions = {}'.format(len(self._curr_soln['y_one_rxns']))
                    print 'Objective value: {}'.format(self._curr_soln['objective_value'])

                    # Reactions added only in forward direction
                    print '\nReactions added in the forward direction only:'
                    for rxn_f in [j for j in self.optModel.J if j in self._curr_soln['forward_rxns'] and j not in self._curr_soln['backward_rxns']]:
                        rxn_obj = super_model.reactions_by_id[rxn_f]
                        print '{}: name = {}  ,  equation = {} ,  external_type = {}  ,  template_rxn_type = {}  ,  base_cost = {} ,  forward cost = {}  , backward cost = {}\n'.format(rxn_f, rxn_obj.name, rxn_obj.get_equation(show_cpds_by = 'name'), rxn_obj.external_type, rxn_obj.template_rxn_type, rxn_obj.base_cost, rxn_obj.forward_cost, rxn_obj.backward_cost)

                    # Reactions added only in backward direction
                    print '\nReactions added in the backward direction only:'
                    for rxn_b in [j for j in self.optModel.J if j not in self._curr_soln['forward_rxns'] and j in self._curr_soln['backward_rxns']]:
                        rxn_obj = super_model.reactions_by_id[rxn_b]
                        print '{}: name = {}  ,  equation = {} ,  external_type = {}  ,  template_rxn_type = {}  ,  base_cost = {} ,  forward cost = {}  , backward cost = {}\n'.format(rxn_b, rxn_obj.name, rxn_obj.get_equation(show_cpds_by = 'name'), rxn_obj.external_type, rxn_obj.template_rxn_type, rxn_obj.base_cost, rxn_obj.forward_cost, rxn_obj.backward)

                    # Reactions added as reversible 
                    print '\nReactions added as reversible:'
                    for rxn in [j for j in self.optModel.J if j in self._curr_soln['forward_rxns'] and j in self._curr_soln['backward_rxns']]:
                        rxn_obj = super_model.reactions_by_id[rxn_f]
                        print '{}: name = {}  ,  equation = {} ,  external_type = {}  ,  template_rxn_type = {}  ,  base_cost = {} ,  forward cost = {}  , backward cost = {}\n'.format(rxn_f, rxn_obj.name, rxn_obj.get_equation(show_cpds_by = 'name'), rxn_obj.external_type, rxn_obj.template_rxn_type, rxn_obj.base_cost, rxn_obj.forward_cost, rxn_obj.backward)


            # If an optmal solution was not obtained 
            else: 
                done = True
 
                if self.stdout_msgs:
                    print '\nEnded because the optimization problem was not solved to optimality: exit_flag = {}\n'.format(self._curr_soln['exit_flag'])

        if self.results_filename != '':
            with open(self.results_filename,'w') as f:
                f.write(']\n')

        if found_solutions_num == max_soln_num:
            done = True
