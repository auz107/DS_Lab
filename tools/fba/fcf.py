from __future__ import division
import sys, time
sys.path.append('../../')
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from coopr.pyomo import *
from coopr.opt import *
from copy import deepcopy
from coopr.environ import *
from tools.pyomoSolverCreator import pyomoSolverCreator
from tools.globalVariables import *
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from copy import deepcopy
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

class fcf(object):
    """
    Performs flux coupling analysis 
    """

    def __init__(self,model, blocked_rxns = [], build_new_optModel = True, optimization_solver = default_optim_solver, results_filename = '', simulation_conditions = '', warmstart = False, warnings = True, stdout_msgs = True):
        """
        INPUTS:
        -------
                        model: An instance of class model containing the information
                               about the metabolic model
                 blocked_rxns: A list of reactions that are blocked under the given uptake and
                               aeration conditions
           build_new_optModel: If True, a new pyomo optimization model is created
          optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                               allowable choices are cplex and gurobi
             results_filename: A string containing the name of the file to save the results in
                             . If an empty string is provided the results are not saved to a 
                               file
        simulation_conditions: A string describing simulation conditions
                    warmstart: If True, warmstart is used
                     warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                               screen or not. The default is True  
                  stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                               Eligible values are True and False and the default is True 
    
        OUTPUTS:
        --------
    
        Ali R. Zomorrodi - Segre Lab @ Boston University
        Last updated: 04-04-2016 
        """   
        # model
        self.model = model

        # Always blocked reactions
        self.blocked_rxns = blocked_rxns

        # build_new_optModel
        self.build_new_optModel = build_new_optModel

        # optimization solver
        self.optimization_solver = optimization_solver

        # results_filename
        self.results_filename = results_filename

        # simulation_conditions
        self.simulation_conditions = simulation_conditions

        # warmstart
        self.warmstart = warmstart

        # warnings and stdout_msgs
        self.warnings = warnings
        self.stdout_msgs = stdout_msgs
        
        # Create the optimization solver
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

        # Always blocked reactions 
        if attr_value == 'blocked_rxns' and not isinstance(simulation_conditions,list):
            raise TypeError('blocked_rxns must be a list of strings')
        elif attr_value == 'blocked_rxns' and len([k for k in attr_value if not isinstance(k,str)]) > 0:
            raise TypeError('blocked_rxns must be a list of strings. Non-string objects observed in the list')
 
        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Simulation conditions name
        if attr_name == 'simulation_conditions' and not isinstance(attr_value,str):
            raise userError('simulation_conditions must be a striing')

        # build_new_optModel 
        if attr_name == 'build_new_optModel' and not isinstance(attr_value,bool):
            raise TypeError("'build_new_optModel' must be either True or False")

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("'stdout_msgs' must be either True or False")
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("'warnings' must be either True or False")

        self.__dict__[attr_name] = attr_value

    def split_rev_rxns(self):
        """
        Decomposes all reversible and exchange reactions in the model into irreversible forward and backward reactions
        """
        # First perform FBA for the original model
        for rxn in self.model.reactions:
            rxn.objective_coefficient = 0
        self.model.biomass_reaction.objective_coefficient = 1
        self.model.fba(stdout_msgs = False)
        max_biomass_model = self.model.fba_model.solution['objective_value']

        # map between reaction ids in model_split and those in the original model for reversible reactions
        # keys are ids of forward and backward reactions in model_split and values are their corresponding reaction id
        # in the original model
        self.model_split_model_map = {}

        # List of backward reactions
        backward_rxns = []
    
        self.model_split = deepcopy(self.model)

        for rxn in [r for r in self.model_split.reactions if r.reversibility.lower() in ['reversible','exchange']]:
            rxn_id = rxn.id
            # Forward reaction
            if rxn.reversibility.lower() == 'reversible':
                rxn.reversibility = 'reversible_forward'
                reversibility_backward = 'reversible_backward'
            elif rxn.reversibility.lower() == 'exchange':
                rxn.reversibility = 'exchange_forward'
                reversibility_backward = 'exchange_backward'

            rxn.id = rxn_id + '_f'
            self.model_split_model_map[rxn_id + '_f'] = rxn_id
 
            # backward reaction
            rxn_b = reaction(id = rxn_id + '_b', stoichiometry = dict([(c,-rxn.stoichiometry[c]) for c in rxn.compounds]), reversibility = reversibility_backward, name = rxn.name, name_aliases = rxn.name_aliases, KEGG_id = rxn.KEGG_id, ModelSEED_id = rxn.ModelSEED_id, BiGG_id = rxn.BiGG_id, EC_numbers = rxn.EC_numbers, subsystem = rxn.subsystem, pathways = rxn.pathways, compartment = rxn.compartment, genes = rxn.genes, gene_reaction_rule = rxn.gene_reaction_rule, objective_coefficient = rxn.objective_coefficient, flux = rxn.flux, store_flux = rxn.store_flux, flux_bounds = rxn.flux_bounds, deltaG = rxn.deltaG, deltaG_uncertainty = rxn.deltaG_uncertainty, deltaG_range = rxn.deltaG_range, kinetics = rxn.kinetics, kinetic_compounds = rxn.kinetic_compounds, confidence_level = rxn.confidence_level, notes = rxn.notes, warnings = rxn.warnings)
            rxn_b.id = rxn_id + '_b'
            self.model_split_model_map[rxn_id + '_b'] = rxn_id

            for cpd in rxn_b.compounds:
                cpd.product_reactions = list(set(cpd.product_reactions))
                if rxn_b.stoichiometry[cpd] < 0:
                    cpd.reactant_reactions.append(rxn_b)
                elif rxn_b.stoichiometry[cpd] > 0:
                    cpd.product_reactions.append(rxn_b)
                cpd.reactions.append(rxn_b)
                cpd.reactions = list(set(cpd.reactions))
                cpd.reactant_reactions = list(set(cpd.reactant_reactions))
    
            backward_rxns.append(rxn_b)

        self.model_split.add_reactions(backward_rxns)
        self.model_split.assign_props()
        self.model_split.validate()
     
        # Reset flux bounds
        self.model_split.reset_flux_bounds()

        # Reset the uptake conditions
        for rxn in [r for r in self.model.reactions if r.reversibility.lower() == 'exchange' and r.flux_bounds[0] < 0]:
            self.model_split.reactions_by_id[rxn.id + '_b'].flux_bounds[1] = -rxn.flux_bounds[0] 
        for rxn in [r for r in self.model.reactions if r.reversibility.lower() == 'exchange' and r.flux_bounds[0] > 0]:
            self.model_split.reactions_by_id[rxn.id + '_f'].flux_bounds[0] = rxn.flux_bounds[0] 
        for rxn in [r for r in self.model.reactions if r.reversibility.lower() == 'exchange' and r.flux_bounds[1] not in [0,1000]]:
            self.model_split.reactions_by_id[rxn.id + '_f'].flux_bounds[1] = rxn.flux_bounds[1] 

        for rxn in [r for r in self.model.reactions if (r.reversibility.lower() != 'exchange' and r.flux_bounds[0] not in [0,-1000])]:
            # Reversible reaction 
            if rxn.reversibility.lower() == 'reversible':
                if rxn.flux_bounds[0] > 0:
                    self.model_split.reactions_by_id[rxn.id + '_f'].flux_bounds[0] = rxn.flux_bounds[0]
                elif rxn.flux_bounds[0] < 0: # Here, LB in the reversible rxn corresponds to UB in the irreversible backward rxn
                    self.model_split.reactions_by_id[rxn.id + '_b'].flux_bounds[1] = -rxn.flux_bounds[0]
            # irreversible reaction 
            if rxn.reversibility.lower() == 'irreversible':
               self.model_split.reactions_by_id[rxn.id].flux_bounds[0] = rxn.flux_bounds[0]

        for rxn in [r for r in self.model.reactions if (r.reversibility.lower() != 'exchange' and r.flux_bounds[1] not in [0,1000])]:
            # Reversible reaction 
            if rxn.reversibility.lower() == 'reversible':
                if rxn.flux_bounds[1] > 0:
                    self.model_split.reactions_by_id[rxn.id + '_f'].flux_bounds[1] = rxn.flux_bounds[1]
                elif rxn.flux_bounds[1] < 0: # Here, UB in the reversible rxn corresponds to LB in the irreversible backward rxn
                    self.model_split.reactions_by_id[rxn.id + '_b'].flux_bounds[0] = -rxn.flux_bounds[1]
            # irreversible reaction 
            if rxn.reversibility.lower() == 'irreversible':
               self.model_split.reactions_by_id[rxn.id].flux_bounds[1] = rxn.flux_bounds[1]

        # Regulation reactions
        for rxn in [r for r in self.model.reactions if r.flux_bounds == [0,0]]:
            # Reversible reaction 
            if rxn.reversibility.lower() == 'reversible':
                    self.model_split.reactions_by_id[rxn.id + '_f'].flux_bounds = [0,0] 
                    self.model_split.reactions_by_id[rxn.id + '_b'].flux_bounds = [0,0] 
            else: 
               self.model_split.reactions_by_id[rxn.id].flux_bounds = [0,0] 

        #-- Make sure there are no reactions with a negative LB in irrev_model ---
        if len([r for r in self.model_split.reactions if r.flux_bounds[0] < 0 or r.flux_bounds[1] < 0]) > 0:
            raise userError('The following reactions have a negative LB or UB fpr the flux in the model with dcomposed reversible and exchange reactions: {}'.format([(r.id, r.flux_bounds) for r in self.model_split if r.flux_bounds[0] < 0 or r.flux_bounds[1] < 0])) 

        #--- pirform FBA for irrev_model --- 
        for rxn in self.model_split.reactions:
            rxn.objective_coefficient = 0
        self.model_split.biomass_reaction.objective_coefficient = 1
        self.model_split.fba(stdout_msgs = False)
        max_biomass_irrev_model = self.model_split.fba_model.solution['objective_value']

        if max_biomass_model != max_biomass_irrev_model:
            raise userError('max biomass for the original model ({}) is not equal to max biomass for the model with splitd reversible and exchange reactions ({})'.format(max_biomass_model, max_biomass_irrev_model))

    def create_ignored_list(self):
        """
        Creates the list of reactions that must be ignore in flux couplign analysis
        """
        self._ignore = self.blocked_rxns

        for rxn in self.model_split.reactions:
            if rxn.reversibility.lower() in ['exchange_forward', 'exchange_backward']:
                self._ignore.append(rxn.id)

    def build_optModel(self):
        """
        This optModel creates a pyomo model for FBA optModel
        """
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()

        #--- Define sets ---
        # Set of compounds 
        optModel.I = Set(initialize = self.model_split.compounds_by_id.keys())

        # Set of rxns  
        optModel.J = Set(initialize = self.model_split.reactions_by_id.keys())

        #--- Define the optModel variables --- 
        optModel.v = Var(optModel.J, domain = NonNegativeReals, bounds = lambda optModel, j: self.model_split.reactions_by_id[j].flux_bounds)

        optModel.t = Var(domain = NonNegativeReals)

        #--- Defiine the objective function and constraints ----
        # Objective function
        optModel.objectiveFunc = Objective(rule = lambda optModel: optModel.v[self._curr_j1], sense = maximize)

        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model_split.compounds_by_id[i]]*optModel.v[j.id] for j in self.model_split.compounds_by_id[i].reactions) == 0)

        # LB(j)*t <= v(j) <= UB(j)*t
        optModel.vLB_const = Constraint([j for j in optModel.J if self.model_split.reactions_by_id[j].flux_bounds[0] != 0], rule = lambda optModel, j: optModel.v[j] >= self.model_split.reactions_by_id[j].flux_bounds[0]*optModel.t)
        optModel.vUB_const = Constraint(optModel.J, rule = lambda optModel, j: optModel.v[j] <= self.model_split.reactions_by_id[j].flux_bounds[1]*optModel.t)

        # v(j2) = 1 for the current reaction j2
        optModel.denom_const = Constraint(rule = lambda optModel: optModel.v[self._curr_j2] == 1)

        self.optModel = optModel

    def run(self):
        """ 
        This method runs the FBA tool. 
        """ 
        # Total processing and wall time required to create the pyomo model, solve it and store the results
        start_total_pt = time.clock()
        start_total_wt = time.time()

        # Split reversible reactions into forward and backward
        print 'Split reversible rxns in the model ...'
        self.split_rev_rxns()
        print 'Done with spliting reversible rxns in the model ...'

        # Create the list of rections that must be ignored
        self.create_ignored_list()

        #---- Creating and instantiating the optModel ----
        start_create_optModel_pt = time.clock()
        start_create_optModel_wt = time.time()
        if self.build_new_optModel:
            # Assign two preliminary values to self._curr_j1 and self_curr_j2 so the model can be created.
            # These two will get updated inside the loops
            self._curr_j1 = self.model_split.reactions[0].id
            self._curr_j2 = self.model_split.reactions[1].id

            self.build_optModel()

            elapsed_create_optModel_pt = str(timedelta(seconds = time.clock() - start_create_optModel_pt))
            elapsed_create_optModel_wt = str(timedelta(seconds = time.time() - start_create_optModel_wt))

            if self.warmstart:
                self.optimization_solver == 'gurobi_ampl'
                self.optSolver = pyomoSolverCreator(self.optimization_solver)

                self.optSolver.options['presolve'] = 0  # Turn off presolve
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

        # list of reactions ids in the model
        rxn_ids = [r.id for r in self.model_split.reactions[:100]]

        # A parameter showing if a current reactions has already been coupled with another reaction
        already_coupled = dict([(r,False) for r in rxn_ids])

        # A list of dictionaries containing the details of directionally, partiually and fully coupled rxns
        self.directionally_coupled_split = []
        self.partially_coupled_split = {}
        self.fully_coupled_split = []

        # A dictioaary with keys being a rxns id representative of a fully coupled sets and values being a list of reactions 
        # ids fully coupled to it in the model_split and in the original model, e.g.,  {j1:[j5,j10], j20:[j40,j70,8100]}
        self.fully_coupled_reps_split = {}
        self.fully_coupled_reps = {}

        # By using rxn_ids[:len(rxn_ids)-1] we essential ask to not consider the last element of
        # of rxn_ids in the outer loop
        for j1 in [jj for jj in rxn_ids[:len(rxn_ids)-1] if (jj not in self._ignore and not already_coupled[jj])]:

            coupled_with_j1 = []

            for j2 in [jj for jj in rxn_ids if jj not in self._ignore and rxn_ids.index(jj) > rxn_ids.index(j1)]:
        
                self._curr_j1 = j1
                self._curr_j2 = j2

                #-- Maximize ---
                self.optModel.del_component('objectiveFunc')
                self.optModel.objectiveFunc = Objective(rule = lambda optModel: optModel.v[self._curr_j1], sense = maximize)
                self.optModel.preprocess()

                optSoln = self.optSolver.solve(self.optModel,tee = False, warmstart = True)

                if str(optSoln.solver.termination_condition).lower() in ['optimal','globallyoptimal','locallyoptimal']:
                    self.optModel.load(optSoln)
                    Rmax =  self.optModel.objectiveFunc()
                elif str(optSoln.solver.termination_condition).lower() == 'unbounded':
                    Rmax = float('inf') 
                else:
                    Rmax = None

                #-- Minimize ---
                self.maximize = False
                self.optModel.del_component('objectiveFunc')
                self.optModel.objectiveFunc = Objective(rule = lambda optModel: optModel.v[self._curr_j1], sense = minimize)
                self.optModel.preprocess()

                optSoln = self.optSolver.solve(self.optModel,tee = False, warmstart = True)

                if str(optSoln.solver.termination_condition).lower() in ['optimal','globallyoptimal','locallyoptimal']:
                    self.optModel.load(optSoln)
                    Rmin =  self.optModel.objectiveFunc()
                else:
                    Rmin = None 

                print '(j1,j2) = {}, (self_curr_j1, self_curr_j2) = {}, (Rmin, Rmax) = {}'.format((j1,j2), (self._curr_j1, self._curr_j2), (Rmin, Rmax))

                #--- Check for coupling relations ---
                if Rmin != None and Rmax != None and Rmax != float('inf') and Rmin == 0 and Rmax > 0:
                    self.directionally_coupled_split.append({'j1':j1,'j2':j2,'direction':'j1 --> j2', 'Rmin':Rmin, 'Rmax':Rmax})

                elif Rmin != None and Rmax != None and Rmin > 0 and Rmax == float('inf'):
                    self.directionally_coupled_split.append({'j1':j1,'j2':j2,'direction':'j2 --> j1', 'Rmin':Rmin, 'Rmax':Rmax})

                elif Rmin != None and Rmax != None and Rmax != float('inf') and Rmin > 0 and Rmax > 0 and Rmax - Rmin >= 1e-6:
                    self.partially_coupled_split.append({'j1':j1,'j2':j2,'relation':'j1 <--> j2', 'Rmin':Rmin, 'Rmax':Rmax})
                    already_coupled[j2] = True

                elif Rmin != None and Rmax != None and Rmax != float('inf') and (Rmin > 0 and Rmax > 0) and Rmax - Rmin < 1e-6:
                    self.fully_coupled_split.append({'j1':j1,'j2':j2,'relation':'j1 <==> j2', 'Rmin':Rmin, 'Rmax':Rmax})
                    already_coupled[j2] = True
                    coupled_with_j1.append(j2)     
    
            if len(coupled_with_j1) > 0:
                self.fully_coupled_reps_split[j1] = coupled_with_j1

        # Fully couples sets for the original model
        for j1 in self.fully_coupled_reps_split.keys():
            # Irreversible reactions fully coupled to this reaction
            irrev_rxns = [j2 for j2 in self.fully_coupled_reps_split[j1] if self.model_split.reactions_by_id[j1].reversibility.lower() == 'irreversible']

            # Coupled reversible forward and backward reactions
            rev_f_rxns = [j2 for j2 in self.fully_coupled_reps_split[j1] if self.model_split.reactions_by_id[j1].reversibility.lower() == 'reversible_forward']
            rev_b_rxns = [j2 for j2 in self.fully_coupled_reps_split[j1] if self.model_split.reactions_by_id[j1].reversibility.lower() == 'reversible_backward']

            # If the representative reaction is irreversible
            if self.model_split.reactions_by_id[j1].reversibility.lower() == 'irreversible':
                # All irreversible reactions are coupled to j1
                coupled_irrev_rxns = irrev_rxns

                # This reaction is fully coupled to a reversible reaction only only if it fully coupled to its forward and 
                # backward reaction (see
                coupled_rev_rxns = []
                for rf in rev_f_rxns:
                     if self.model_split_model_map[rf]  + '_b' in rev_b_rxns: # If its backward rxn is in rev_b_rxns too
                        coupled_rev_rxns.append(self.model_split_model_map[rf])

                if len(coupled_irrev_rxns + coupled_rev_rxns) > 0:
                    self.fully_coupled_reps[j1] = coupled_irrev_rxns + coupled_rev_rxns

            # If the representative reaction is reversible forward or reversible backward
            if self.model_split.reactions_by_id[j1].reversibility.lower() in ['reversible_forward','reversible_backward']:
                if (self.model_split.reactions_by_id[j1].reversibility.lower() in 'reversible_forward' and self.model_split_model_map[j1] + '_b' in rev_b_rxns) or (self.model_split.reactions_by_id[j1].reversibility.lower() in 'reversible_backward' and self.model_split_model_map[j1] + '_f' in rev_b_rxns): 
                    coupled_irrev_rxns = irrev_rxns

                    coupled_rev_rxns = []
                    for rf in rev_f_rxns:
                        if self.model_split_model_map[rf]  + '_b' in rev_b_rxns: # If its backward rxn is in rev_b_rxns too
                            coupled_rev_rxns.append(self.model_split_model_map[rf])

                        if len(coupled_irrev_rxns + coupled_rev_rxns) > 0:
                            self.fully_coupled_reps[j1] = coupled_irrev_rxns + coupled_rev_rxns

                else: # otherwise check if the rest of reversible reactions are coupled together 
                    coupled_rev_rxns = []
                    for rf in rev_f_rxns:
                        if self.model_split_model_map[rf]  + '_b' in rev_b_rxns: # If its backward rxn is in rev_b_rxns too
                            coupled_rev_rxns.append(self.model_split_model_map[rf])

                    if len(coupled_irrev_rxns + coupled_rev_rxns) > 1:
                        # Choose the first element as the representative
                        rep = coupled_irrev_rxns + coupled_rev_rxns[0]
                        self.fully_coupled_reps[rep] = (coupled_irrev_rxns + coupled_rev_rxns)[1:]


        if self.results_filename != '':
            # Save the results into a temporary file in case there are errors in save_results_toFile
            with open(results_filename + '.tmp','w') as f:
                f.write('directionarlly_coupled = {}\n'.format(self.directionally_coupled))
                f.write('\nparially_coupled = {}\n'.format(self.partially_coupled))
                f.write('\nfully_coupled = {}\n'.format(self.fully_coupled))
            self.save_results_toFile()

        # Time required to perform FBA
        elapsed_total_pt = str(timedelta(seconds = time.clock() - start_total_pt))
        elapsed_total_wt = str(timedelta(seconds = time.time() - start_total_wt))

        if self.stdout_msgs:
            print 'FCF took (hh:mm:ss) {}/{} of processing/walltime\n'.format(elapsed_total_pt, elapsed_total_wt)

        return self.directionally_coupled_split, self.partially_coupled_split, self.fully_coupled_reps_split, self.fully_coupled_reps_split, self.fully_coupled_reps

    def save_results_toFile(self):
        """
        Save the results into a file
        """
        with open(self.results_filename, 'w') as f:

            # Directionally coupled
            f.write('directionarlly_coupled_split = [\n')
            for dcoupled in self.directionarlly_coupled:
                f.write("{{'j1':'{}', 'j2':'{}', 'direction':'{}', 'Rmin':{}, 'Rmax':{}}},\n".format(dcoupled['j1'],dcoupled['j2'], dcoupled['direction'], dcoupled['Rmin'], doucpled['Rmax'])) 
            f.write(']\n\n')

            # Partially coupled
            f.write('partially_coupled_split = [\n')
            for pcoupled in self.partially_coupled_split:
                f.write("{{'j1':'{}', 'j2':'{}', 'Rmin':{}, 'Rmax':{}}},\n".format(pcoupled['j1'],pcoupled['j2'], pcoupled['Rmin'], doucpled['Rmax'])) 
            f.write(']\n\n')

            # Fully coupled
            f.write('fully_coupled_split = [\n')
            for fcoupled in self.fully_coupled_reps_split:
                f.write("{{'j1':'{}', 'j2':'{}', 'Rmin':{}, 'Rmax':{}}},\n".format(fcoupled['j1'],fcoupled['j2'], fcoupled['Rmin'], doucpled['Rmax'])) 
            f.write(']\n\n')

            # Fully coupled sets
            f.write('fully_coupled_reps_split = {\n')
            for k in self.fully_coupled_reps_split.keys():
                f.write("'{}':{},".format(k),fully_coupled_reps_split[j])
            f.write('}')

            # Fully coupled sets for the original model
            f.write('fully_coupled_reps = {\n')
            for k in self.fully_coupled_reps.keys():
                f.write("'{}':{},".format(k),fully_coupled_reps[j])
            f.write('}')



def find_coupling_with_revrxn(self, j1_f_j2_coupling, j1_b_j2_coupling):
    """
    Assume j1 and j2 are irreversible and reversible reactions in the original model. Then given the coupling relation
    between j1 and j2_f and j1 and j2_b, it identifies the coupling relation between (i,j)

    1. Given the coupling relations for (j1_f,j2_f) & (j1_f,j2_b) find the coupling relations between (j1_f,j2)
    2. Given the coupling relations for (j1_b,j2_f) & (j1_b,j2_b) find the coupling relations between (j1_b,j2)
    2. Given the coupling relations for (j1_f,j2) & (j1_b,j2) find the coupling relations between (j1,j2)

    According to the definition of this function an irreversible reaction is fully coupled to a reversible reaciton
    only if it is fully coupled to both its forward and backard reactions. Similarly, two reversible reactions are
    fully coupled only if the reverse and forward reactions of each reversible reaction are fully coupled to those of the
    other reactions 

    Reference: PMID: 21676263 (See Algorithm 2 [CouplingRelationRecompute] in Additionaly File 1)
    """
    if j1_f_j2_coupling == '<==>' and j1_b_j2_coupling == '<-->':
        coupling_relation = '<-->'
    elif j1_f_j2_coupling == '<==>' and j1_b_j2_coupling == '-->':
        coupling_relation = '-->'
    elif j1_f_j2_coupling == '<==>' and j1_b_j2_coupling == '<--':
        coupling_relation = '<--'
    elif j1_f_j2_coupling == '<==>' and j1_b_j2_coupling == 'uncoupled':
        coupling_relation = 'uncoupled'

    elif j1_f_j2_coupling == '<-->' and j1_b_j2_coupling == '-->':
        coupling_relation = '-->'
    elif j1_f_j2_coupling == '<-->' and j1_b_j2_coupling == '<--':
        coupling_relation = '<--'
    elif j1_f_j2_coupling == '<-->' and j1_b_j2_coupling == 'uncoupled':
        coupling_relation = 'uncoupled'

    elif j1_f_j2_coupling == '-->' and j1_b_j2_coupling == '<--':
        coupling_relation = 'uncoupled'
    elif j1_f_j2_coupling == '-->' and j1_b_j2_coupling == 'uncoupled':
        coupling_relation = 'uncoupled'

    elif j1_f_j2_coupling == '<--' and j1_b_j2_coupling == '-->':
        coupling_relation = 'uncoupled'
    elif j1_f_j2_coupling == '<--' and j1_b_j2_coupling == 'uncoupled':
        coupling_relation = 'uncoupled'

    else:
        coupling_relation = j1_f_j2_coupling

    return coupling_relation

