from __future__ import division
import sys, time
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.pyomoSolverCreator import *
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction

# The following lines change the temporary directory for pyomo
from pyutilib.services import TempfileManager
TempfileManager.tempdir = '/tmp/'

class break_cycles(object):
    """
    Breaks cycles in a metabolic network to avoid a non-zero biomass production in the 
    absence of a carbon source. 

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 10-23-2015 
    """   

    def __init__(self,model, max_modification_num = 1, max_solution_num = 1, max_biomass_thr = 1e-6, penalties = {'zeroLB_fwdProb':0,'zeroUB_fwdProb':2,'zeroLB_bwdProb':2,'zeroUB_bwdProb':0,'zeroLB_revProb':2,'zeroUB_revProb':2,'zeroLB_noProb':1,'zeroUB_noProb':1},optimization_solver = 'gurobi', build_new_optModel = True, stdout_msgs = True, warnings = True,  **additional_args): 
        """
        INPUTS (required):
        ------
                      model: An instance of class model containing the model to be gap filled.
                             Each reaciton object of the model must have an attribute named probMeas, 
                             which should be either a float between zero and one or None. 

        INPUTS (optional):
        ------
           gapfilling_method: The method to perform gap filling. The default is parsimony-based gap filling, where the
                              aim is to minimally perturb the model.
        max_modification_num: Maximum allowed number of modifications
            max_solution_num: Maximum number of alternative solutions that must be returned
             max_biomass_thr: Allowed threshold for maximum biomass formation when breaking cycles, i.e., break cycles such that
                              the max biomass formation in the network goes below this threshould
                   penalties: Penalties (cost) associated with performing each modification to the model
                                 zeroLB_fwdProb: Set LB to zero for reactions that must be 
                                                 irreversible in the forward direciton according to their 
                                                 probability measure.
                                 zeroUB_fwdProb: Set UB to zero for reactions that must be 
                                                 irreversible in the forward direciton according to their 
                                                 probability measure.
                                 zeroLB_bwdProb: Set LB to zero for reactions that must be 
                                                 irreversible in the backward direciton according to their 
                                                 probability measure.
                                 zeroUB_bwdProb: Set UB to zero for reactions that must be 
                                                 irreversible in the backward direciton according to their 
                                                 probability measure.
                                 zeroLB_revProb: Set LB to zero for reactions that must be 
                                                 reversible according to their probability measure.
                                 zeroUB_revProb: Set UB to zero for reactions that must be 
                                                 reversible according to their probability measure.
                                  zeroLB_noProb: Set LB to zero for reactions with no available probability measure
                                  zeroUB_noProb: Set UB to zero for reactions with no available probability measure
         optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                              allowable choices are cplex and gurobi
          build_new_optModel: A parameter indicating whether a new pyomo optimizaiton model should be 
                              created (True) or or an existing one should be used (False). The options is useful
                              for the cases a model is already created and one just 
                              wants to change some model attributes (e.g., flux bounds)
                              and rerun FBA. Setting this parameter to False will save 
                              some runtime as the model need not to be created again.
                 stdout_msgs: By default (on) writes  a summary including the solve 
                              status, optimality status (if not optimal), objective 
                              function value and the elapsed time on the screen.
                              if set to a value of 'off' no resuults are written on 
                              the screen, in which case The user can instead specifiy 
                              an output fiile using the option outputFile, or store 
                              them in the variable solution (see the 'run' method for
                              details)
                    warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                              screen or not. The default is True  
                 stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                              Eligible values are True and False and the default is True 
             additional_args: Additional arguments should be entered as normal but they are 
                              converted to a dictionary whose keys are the names of the arguments and 
                              values are the values of  those arguments

        These are the outputs of the method 'run'
        """
       
        # Metabolic model
        self.model = model

        # Maximum allowed number of modifications to the model
        self.max_modification_num = max_modification_num

        # Maximum number of alternative solutions to return
        self.max_solution_num = max_solution_num

        # Threshould on the maximum biomass formation when breaking cycles
        self.max_biomass_thr = max_biomass_thr

        # Penalties for incorporating each type of modification to the model
        self.penalties = penalties

        # Create new optimization model
        self.build_new_optModel = build_new_optModel

        # Solver name
        self.optimization_solver = optimization_solver

        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.warnings = warnings

        # Additoinal arguments. Additional arguments should be entered as normal but they are 
        # converted to a dictionary whose keys are the names of the arguments and values are 
        # the values of  those arguments
        argnames = additional_args.keys()
        argvals = additional_args.values()
        for argname in argnames:
           exec "self." + argname + " = " +"additional_args['" + argname + "']"

        # Define parameters required for the optimization model
        self.define_optModel_params()

    def check_attr(self,attr_name,attr_value):
        """
        Checks the conditions on the class attributes
 
        INPUTS:
        -------
         attr_name: Attribute name
        attr_value: Attribute vlaue
        """
        # model 
        if attr_name == 'model':
            # Reactions with no objective coefficient
            no_objCeoff_rxns = [rxn for rxn in attr_value.reactions if rxn.objective_coefficient == None]
            if len(no_objCeoff_rxns) > 0 and len(no_objCeoff_rxns) <= 10:
                raise ValueError('objective_coefficient not defined for ' + str(len(no_objCeoff_rxns)) + ' reactions including ' + str([r.id for r in no_objCeoff_rxns[:10]])) 
            elif len(no_objCeoff_rxns) > 0 and len(no_objCeoff_rxns) > 10:
                raise ValueError('objective_coefficient not defined for ' + str(len(no_objCeoff_rxns)) + ' reactions including ' + str([r.id for r in no_objCeoff_rxns[:10]]) + ' and ' + str(len(no_objCeoff_rxns) - 10) + 'more reactions.') 
            
            # Reactions with an invalid format for probMeas
            invalid_probMeasFormat_rxns = [rxn for rxn in attr_value.reactions if hasattr(rxn,'probMeas') and not isinstance(rxn.probMeas,float) and not isinstance(rxn.probMeas,int) and rxn.probMeas != Non]
            if len(invalid_probMeasFormat_rxns) > 0 and len(invalid_probMeasFormat_rxns) <= 10:
                raise ValueError('Invalid data type for probMeas for ' + str(len(invalid_probMeasFormat_rxns)) + ' reactions including ' + str([r.id for r in invalid_probMeasFormat_rxns[:10]])) 
            elif len(invalid_probMeasFormat_rxns) > 0 and len(invalid_probMeasFormat_rxns) > 10:
                raise ValueError('Invalid data type for probMeas for ' + str(len(invalid_probMeasFormat_rxns)) + ' reactions including ' + str([r.id for r in invalid_probMeasFormat_rxns[:10]]) + ' and ' + str(len(invalid_probMeasFormat_rxns) - 10) + 'more reactions.') 

            # Reactions with invalid values for probMeas (out of [0,1] range
            invalid_probMEasValue_rxns = [rxn for rxn in attr_value.reactions if hasattr(rxn,'probMeas') and rxn.probMeas != None and (rxn.probMeas < 0 or rxn.probMeas > 1)]
            if len(invalid_probMEasValue_rxns) > 0 and len(invalid_probMEasValue_rxns) <= 10:
                raise ValueError('Invalid data type for probMeas for ' + str(len(invalid_probMEasValue_rxns)) + ' reactions including ' + str([r.id for r in invalid_probMEasValue_rxns[:10]])) 
            elif len(invalid_probMEasValue_rxns) > 0 and len(invalid_probMEasValue_rxns) > 10:
                raise ValueError('Invalid data type for probMeas for ' + str(len(invalid_probMEasValue_rxns)) + ' reactions including ' + str([r.id for r in invalid_probMEasValue_rxns[:10]]) + ' and ' + str(len(invalid_probMEasValue_rxns) - 10) + 'more reactions.') 

        # Maximum allowed number of modifications 
        if attr_name == 'max_modification_num' and not isinstance(attr_value,int):
            raise TypeError('max_modification_num must be a non-negative integer') 

        # Maximum allowed number of modifications 
        if attr_name == 'max_solution_num' and not isinstance(attr_value,int):
            raise TypeError('max_solution_num must be a non-negative integer') 

        # Maximum allowed number of modifications 
        if attr_name == 'max_biomass_thr' and (not isinstance(attr_value,float) or (isinstance(attr_value,float) and attr_value < 0)):
            raise TypeError('max_biomass_thr must be a non-negative float') 

        # Penalties 
        if attr_name == 'penalties' and not isinstance(attr_value,dict):
            raise TypeError('penalties must be a dictionary') 

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise ValueError('Invalid solver name (eligible choices are cplex and gurobi)\n')          

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("stdout_msgs must be either True or False")

        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("warnings must be either True or False")

    def __setattr__(self,attr_name,attr_value):
         """
         Redefines funciton __setattr__
         INPUTS:
         -------
         attr_name: Attribute name
         attr_value: Attribute value
         """
         if attr_name in ['model','max_modification_num','penalties','optimization_solver','stdout_msgs','warnings']: 
             self.check_attr(attr_name,attr_value)
         self.__dict__[attr_name] = attr_value

    #------------------------------------------------------------------------
    #--- Define parameters needed for the optimizaiton problem ---
    #------------------------------------------------------------------------
    def define_optModel_params(self): 
        """
        This function
           - Assigns the type of reaction based on the probability measure (threshold from PMID: 26147299)
           - Assigns the cost of setting the lower or upper bound to zero
           - Reactions whose lower and/or upper bound must not be manipulated or there is no need to manipulate them
        """
        self._muLB_max = 1e4
        self._muUB_max = 1e4

        # Initialize ignoreLB and ignoreUB if it hasn't already been done
        for rxn in self.model.reactions:
            if not hasattr(rxn,'ignoreLB'):
                rxn.ignoreLB = False
            if not hasattr(rxn,'ignoreUB'):
                rxn.ignoreUB = False

            # If this ia blocked rxn under the examined condition
            if hasattr(rxn,'blocked') and rxn.blocked:
                rxn.ignoreLB = True
                rxn.ignoreUB = True

        for rxn in self.model.reactions:
            if not hasattr(rxn,'probMeas') or rxn.probMeas == None:
                rxn.probType = 'noProb' 
                rxn.costLB = self.penalties['zeroLB_noProb']
                rxn.costUB = self.penalties['zeroUB_noProb']
            elif rxn.probMeas > 0.7:
                rxn.probType = 'fwdProb'  # Irreversible in the forward direciton
                rxn.costLB = self.penalties['zeroLB_fwdProb']
                rxn.costUB = self.penalties['zeroUB_fwdProb']
            elif rxn.probMeas < 0.3:
                rxn.probType = 'bwdProb'  # Irreversible in the backward direciton
                rxn.costLB = self.penalties['zeroLB_bwdProb']
                rxn.costUB = self.penalties['zeroUB_bwdProb']
            else:
                rxn.probType = 'revProb'  # Reversible
                rxn.costLB = self.penalties['zeroLB_revProb']
                rxn.costUB = self.penalties['zeroUB_revProb']

            # Reactions whose lower and/or upper bound must not be manipulated or there is no need to manipulate them
            if rxn.type.lower() == 'exchange':
                rxn.ignoreLB = True
                rxn.ignoreUB = True

            if rxn.flux_bounds[0] == 0:
                rxn.ignoreLB = True
                rxn.ignoreUB = False

            if rxn.flux_bounds[1] == 0:
                rxn.ignoreLB = False
                rxn.ignoreUB = True
 
        # Do not manipulate biomass reaction and ATPM 
        self.model.biomass_reaction.ignoreLB = True
        self.model.biomass_reaction.ignoreUB = True

        ATPM = self.model.get_reactions({'ATPM':'id'}) 
        if ATPM != None:
            ATPM.ignoreLB = True
            ATPM.ignoreUB = True

    def assignFluxBounds(self,optModel,j):
        """
        Define the flux bounds
        """
        return j.flux_bounds 

    #------------------------------------------------------------------------
    #--- Define rules for the objective function and various constraints ----
    #------------------------------------------------------------------------

    #--- Constraints of the outer-level problem ---
    def outer_objectiveFunc_rule(self,optModel):
        """
        Objective function of the outer problem
        """
        return sum(j.costUB*(1 - optModel.yUB[j]) + j.costLB*(1 - optModel.yLB[j]) for j in optModel.J)

    def total_modifications_const_rule(self,optModel):
        """
        Restriction on the total number of modifications 
        """
        return sum((1 - optModel.yUB[j]) + (1 - optModel.yLB[j]) for j in optModel.J) <= self._max_modification_num_curr

    def max_biomass_const_rule(self,optModel):
        """
        Constraint on the maximum biomass formaiton flux 
        """
        return optModel.v[self.model.biomass_reaction] <= self.max_biomass_thr 

    def integer_cuts_rule(self,optModel):
        """
        Integer cuts 
        """
        return sum((1 - optModel.yUB[j]) for j in optModel.J if self._yUBopt_curr[j] == 0) + sum((1 - optModel.yLB[j]) for j in optModel.J if self._yLBopt_curr[j] == 0) <= self._max_modification_num_curr - 1

    #--- Constraints of the primal problem ---
    def primal_objectiveFunc_rule(self,optModel):
        """
        Objective function of the primal problem
        """
        return optModel.v[self.model.biomass_reaction] 

    def massBalance_const_rule(self,optModel,i):
        """
        Mass balance 
        """
        return sum(j.stoichiometry[i]*optModel.v[j] for j in i.reactions) == 0 

    def fluxLB_const_rule(self,optModel,j):
        """
        vj >= LB_j*yLB_j
        """
        return optModel.v[j] >= j.flux_bounds[0]*optModel.yLB[j] 

    def fluxUB_const_rule(self,optModel,j):
        """
        vj <= UB_j*yUB_j
        """
        return optModel.v[j] <= j.flux_bounds[1]*optModel.yUB[j] 

    #--- Constraints of the dual problem ---
    def dual_objectiveFunc_rule(self,optModel):
        """
        Objective function of the dual problem
        """
        return sum(j.flux_bounds[1]*optModel.muUB[j] for j in optModel.J if j.ignoreUB) + \
               sum(j.flux_bounds[1]*optModel.muUByUB[j] for j in optModel.J if not j.ignoreUB) + \
               sum(-j.flux_bounds[0]*optModel.muLB[j] for j in optModel.J if j.ignoreLB) + \
               sum(-j.flux_bounds[0]*optModel.muLByLB[j] for j in optModel.J if not j.ignoreLB)

    def dual_const_rule(self,optModel,j):
        """
        Constraints of the dual problem
        """
        return sum(j.stoichiometry[i]*optModel.Lambda[i] for i in j.compounds) + optModel.muUB[j] - optModel.muLB[j] == j.objective_coefficient 

    def strong_duality_const_rule(self,optModel):
        """
        Constraints of the dual problem
        """
        return optModel.v[self.model.biomass_reaction] == \
               sum(j.flux_bounds[1]*optModel.muUB[j] for j in optModel.J if j.ignoreUB) + \
               sum(j.flux_bounds[1]*optModel.muUByUB[j] for j in optModel.J if not j.ignoreUB) + \
               sum(-j.flux_bounds[0]*optModel.muLB[j] for j in optModel.J if j.ignoreLB) + \
               sum(-j.flux_bounds[0]*optModel.muLByLB[j] for j in optModel.J if not j.ignoreLB)

    def linearize_muUByUB_const1_rule(self,optModel,j):
        """
        muUByUB_j <= muUB_j_max*yUB_j 
        """
        return optModel.muUByUB[j] <= self._muUB_max*optModel.yUB[j] 

    def linearize_muUByUB_const2_rule(self,optModel,j):
        """
        muUByUB_j <= muUB_j
        """
        return optModel.muUByUB[j] <= optModel.muUB[j] 

    def linearize_muUByUB_const3_rule(self,optModel,j):
        """
        muUByUB_j >= muUB_j - muUB_j_max*(1 - yUB_j)
        """
        return optModel.muUByUB[j] >= optModel.muUB[j] - self._muUB_max*(1 - optModel.yUB[j]) 

    def linearize_muLByLB_const1_rule(self,optModel,j):
        """
        muLByLB_j <= muLB_j_max*yLB_j 
        """
        return optModel.muLByLB[j] <= self._muLB_max*optModel.yLB[j] 

    def linearize_muLByLB_const2_rule(self,optModel,j):
        """
        muLByLB_j <= muLB_j
        """
        return optModel.muLByLB[j] <= optModel.muLB[j] 

    def linearize_muLByLB_const3_rule(self,optModel,j):
        """
        muLByLB_j >= muLB_j - muLB_j_max*(1 - yLB_j)
        """
        return optModel.muLByLB[j] >= optModel.muLB[j] - self._muLB_max*(1 - optModel.yLB[j]) 

    #--------------------------------------------------------
    #---------------- Create Optimizaiton models ------------        
    #--------------------------------------------------------
    def build_primal_optModel(self):
        """
        Creates a pyomo optimization model for the primal problem 
        """   
        # Define parameters and scalars needed to define the optimizaiton problem
        self.define_optModel_params()

        # Create a pyomo model optimization model
        optModel = ConcreteModel()
        
        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = self.model.compounds)   

        # Set of rxns  
        optModel.J = Set(initialize = self.model.reactions)     

        #--- Variables --- 
        # Reaction fluxes
        optModel.v = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        # Binary variables in v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j 
        optModel.yLB = Var(optModel.J, domain=Boolean)
        optModel.yUB = Var(optModel.J, domain=Boolean)

        #--- Objective function ----
        # Objective function
        optModel.objectiveFunc = Objective(rule=self.primal_objectiveFunc_rule, sense = maximize)

        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule=self.massBalance_const_rule)

        # v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j
        optModel.fluxLB_const = Constraint(optModel.J, rule=self.fluxLB_const_rule)
        optModel.fluxUB_const = Constraint(optModel.J, rule=self.fluxUB_const_rule)

        self.optModel = optModel
    
    def build_dual_optModel(self):
        """
        Creates a pyomo optimization model for the dual problem 
        """   
        # Define parameters and scalars needed to define the optimizaiton problem
        self.define_optModel_params()

        # Create a pyomo model optimization model
        optModel = ConcreteModel()
        
        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = self.model.compounds)   

        # Set of rxns  
        optModel.J = Set(initialize = self.model.reactions)     

        #--- Variables --- 
        # Dual variables associated with steady-state mass balance constraints
        # It's better to not provide upper and lower bound on lambda otherwise the dual objective 
        # will be slightly different from the primal's
        optModel.Lambda = Var(optModel.I, domain=Reals)
        
        # Dual variables associated with v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j 
        optModel.muLB = Var(optModel.J, domain=Reals, bounds = (0,self._muLB_max))
        optModel.muUB = Var(optModel.J, domain=Reals, bounds = (0,self._muUB_max))

        # Product of muLB_j*yLB_j and muUB_j*yUB_j
        optModel.muLByLB = Var(optModel.J, domain=Reals, bounds = (0,self._muLB_max))
        optModel.muUByUB = Var(optModel.J, domain=Reals, bounds = (0,self._muUB_max))

        # Binary variables in v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j 
        optModel.yLB = Var(optModel.J, domain=Boolean)
        optModel.yUB = Var(optModel.J, domain=Boolean)

        #--- Objective function ----
        # Objective function
        optModel.objectiveFunc = Objective(rule=self.dual_objectiveFunc_rule, sense = minimize)

        #-- Constraints of the dual problem --
        # dual constraints 
        optModel.dual_const = Constraint(optModel.J, rule=self.dual_const_rule)
     
        # Constraints linearizing muLB_j*yLB_j and muUB_j*yUB_j
        optModel.linearize_muLByLB_const1 = Constraint(optModel.J, rule=self.linearize_muLByLB_const1_rule)
        optModel.linearize_muLByLB_const2 = Constraint(optModel.J, rule=self.linearize_muLByLB_const2_rule)
        optModel.linearize_muLByLB_const3 = Constraint(optModel.J, rule=self.linearize_muLByLB_const3_rule)
        optModel.linearize_muUByUB_const1 = Constraint(optModel.J, rule=self.linearize_muUByUB_const1_rule)
        optModel.linearize_muUByUB_const2 = Constraint(optModel.J, rule=self.linearize_muUByUB_const2_rule)
        optModel.linearize_muUByUB_const3 = Constraint(optModel.J, rule=self.linearize_muUByUB_const3_rule)

        self.optModel = optModel
    
    def build_bilevel_optModel(self):
        """
        Creates a pyomo optimization model for the bilevel problem 
        """   
        # Define parameters and scalars needed to define the optimizaiton problem
        self.define_optModel_params()

        # Create a pyomo model optimization model
        optModel = ConcreteModel()
        
        #--- Sets ---
        # Set of compounds 
        optModel.I = Set(initialize = self.model.compounds)   

        # Set of rxns  
        optModel.J = Set(initialize = self.model.reactions)     

        #--- Variables --- 
        # Reaction fluxes
        optModel.v = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        # Dual variables associated with steady-state mass balance constraints
        optModel.Lambda = Var(optModel.I, domain=Reals)
        
        # Dual variables associated with v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j 
        optModel.muLB = Var(optModel.J, domain=Reals, bounds = (0,self._muLB_max))
        optModel.muUB = Var(optModel.J, domain=Reals, bounds = (0,self._muUB_max))

        # Product of muLB_j*yLB_j and muUB_j*yUB_j
        optModel.muLByLB = Var(optModel.J, domain=Reals, bounds = (0,self._muLB_max))
        optModel.muUByUB = Var(optModel.J, domain=Reals, bounds = (0,self._muUB_max))

        # Binary variables in v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j 
        optModel.yLB = Var(optModel.J, domain=Boolean)
        optModel.yUB = Var(optModel.J, domain=Boolean)

        #--- Objective function (outer level)----
        # Objective function
        optModel.objectiveFunc = Objective(rule=self.outer_objectiveFunc_rule, sense = minimize)

        #--- Constraints of the outer-level problem ---
        # Constraint on the total number of modificaitons
        optModel.total_modifications_const = Constraint(rule=self.total_modifications_const_rule)

        # Constraint on the max biomass flux 
        optModel.biomass_const = Constraint(rule=self.max_biomass_const_rule)

        # Integer cuts
        optModel.integer_cuts = ConstraintList(noruleinit=True) 
       
        #-- Constraints of the primal problem --
        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule=self.massBalance_const_rule)

        # v_j >= LB_j*yLB_j and v_j <= UB_j*yUB_j
        optModel.fluxLB_const = Constraint(optModel.J, rule=self.fluxLB_const_rule)
        optModel.fluxUB_const = Constraint(optModel.J, rule=self.fluxUB_const_rule)

        #-- Constraints of the dual problem --
        # dual constraints 
        optModel.dual_const = Constraint(optModel.J, rule=self.dual_const_rule)
     
        # Strong duality 
        optModel.strong_duality_const = Constraint(rule=self.strong_duality_const_rule)

        # Constraints linearizing muLB_j*yLB_j and muUB_j*yUB_j
        optModel.linearize_muLByLB_const1 = Constraint(optModel.J, rule=self.linearize_muLByLB_const1_rule)
        optModel.linearize_muLByLB_const2 = Constraint(optModel.J, rule=self.linearize_muLByLB_const2_rule)
        optModel.linearize_muLByLB_const3 = Constraint(optModel.J, rule=self.linearize_muLByLB_const3_rule)
        optModel.linearize_muUByUB_const1 = Constraint(optModel.J, rule=self.linearize_muUByUB_const1_rule)
        optModel.linearize_muUByUB_const2 = Constraint(optModel.J, rule=self.linearize_muUByUB_const2_rule)
        optModel.linearize_muUByUB_const3 = Constraint(optModel.J, rule=self.linearize_muUByUB_const3_rule)
    
        self.optModel = optModel
     
    def fix_known_variables(self):
        """
        Fix all binary variables to one
        """
        for rxn in [r for r in self.model.reactions if r.ignoreLB]:
            self.optModel.yLB[rxn] = 1
            self.optModel.yLB[rxn].fixed = True
        for rxn in [r for r in self.model.reactions if r.ignoreUB]:
            self.optModel.yUB[rxn] = 1
            self.optModel.yUB[rxn].fixed = True

        # Fix all binary variables to one
        #for rxn in self.model.reactions:
        #    self.optModel.yLB[rxn] = 1
        #    self.optModel.yLB[rxn].fixed = True
        #    self.optModel.yUB[rxn] = 1
        #    self.optModel.yUB[rxn].fixed = True

        #self.optModel.yLB[self.model.get_reactions({'rxn00541_c0':'id'})] = 0
        #self.optModel.yLB[self.model.get_reactions({'rxn00541_c0':'id'})].fixed = True
        #self.optModel.yUB[self.model.get_reactions({'rxn00541_c0':'id'})] = 0
        #self.optModel.yUB[self.model.get_reactions({'rxn00541_c0':'id'})].fixed = True

    def validate_soln(self):
        """
        Validates a solution found by the code using fba. It returns an error if validation fails
        """
        if self.stdout_msgs:
            print 'Validating the solution of break_cycles ...',

        # Original bounds on the modified reactions
        modified_rxns_orig_bounds = {}
        for rxn in self.solution[-1]['zero_yLBopt_rxns']:
            modified_rxns_orig_bounds[rxn] = rxn.flux_bounds
            rxn.flux_bounds[0] = 0
        for rxn in self.solution[-1]['zero_yUBopt_rxns']:
            modified_rxns_orig_bounds[rxn] = rxn.flux_bounds
            rxn.flux_bounds[1] = 0

        self.model.fba(stdout_msgs = False,warnings = False)
        # If fba models is solved to optimality but max biomass is greater than max_biomass_thr
        if self.model.fba_model.solution['exit_flag'] == 'globallyOptimal' and self.model.fba_model.solution['objective_value'] > self.max_biomass_thr:
            raise userError('The max biomass obtained via FBA (' + str(self.model.fba_model.solution['objective_value']) + ') is greater than the threshould imposed in break_cycls code (' + str(self.max_biomass_thr) + ') for the following solution: zero_yLBopt_rxns = ' + str([r.id for r in self.solution['zero_yLBopt_rxns']]) + '    ,   zero_yUBopt_rxns = ' + str([r.id for r in self.solution['zero_yUBopt_rxns']]))

        elif self.model.fba_model.solution['exit_flag'] != 'globallyOptimal': 
            raise userError('FBA problem turned out to be infeasble while checking the validity of the following solution found by break_cycles: zero_yLBopt_rxns = ' + str([r.id for r in self.solution['zero_yLBopt_rxns']]) + '    ,   zero_yUBopt_rxns = ' + str([r.id for r in self.solution['zero_yUBopt_rxns']]))

        else: # The solution was validated successfully. Reset the original bound
            for rxn in self.solution[-1]['zero_yLBopt_rxns'] + self.solution[-1]['zero_yUBopt_rxns']:
                rxn.flux_bounds = modified_rxns_orig_bounds[rxn]  

            if self.stdout_msgs:
                print '  Validated'


    def run_single(self):
        """
        Performs a single run of the code
        """
        # Fix known variables
        self.fix_known_variables()

        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()

        # Instantiate the optModel
        self.optModel.preprocess()

        #---- Solve the model ----
        # Create a solver and set the options
        solverType = pyomoSolverCreator(self.optimization_solver)

        elapsed_preproc_pyomo_pt = time.clock() - start_pyomo_pt
        elapsed_preproc_pyomo_wt = time.time() - start_pyomo_wt

        #- Solve the optModel (tee=True shows the solver output) -
        try:
            start_solver_pt = time.clock()
            start_solver_wt = time.time()
            OptSoln = solverType.solve(self.optModel,tee=False)
            solverFlag = 'normal'

        # In the case of an error switch the solver
        except  Exception, e:
            if self.warnings:
                print '**WARNING (break_cycle.py)! {} failed with the following error: \n{} \nAn alternative solver is tried'.format(self.optimization_solver,e)

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
                    print '**WARNING (break_cycle.py)! {} failed with the following error: \n{} \nAn alternative solver is tried'.format(self.optimization_solver,e)

        elapsed_solver_pt = (time.clock() - start_solver_pt)
        elapsed_solver_wt = (time.time() - start_solver_wt)

        #----- Print the results in the output ------
        if solverFlag == 'normal' and str(OptSoln.solver.termination_condition).lower() == 'optimal':
        
            exit_flag = 'globallyOptimal'

            # Load the results
            self.optModel.load(OptSoln)
        
            # Optimal value of the objective function
            opt_objValue = self.optModel.objectiveFunc()

            # Print the results on the screen 
            if self.stdout_msgs:
                print "\nSolver.status = ",OptSoln.solver.termination_condition
                print "Optimality status = ",exit_flag
                print "Objective value = ",opt_objValue

            # List of reactions for which yLB or yUB is zero 
            zero_yLBopt_rxns = [rxn for rxn in self.model.reactions if self.optModel.yLB[rxn].value == 0]
            zero_yUBopt_rxns = [rxn for rxn in self.model.reactions if self.optModel.yUB[rxn].value == 0]

            # self._yLBopt_curr and self._yUBopt_curr hold the optimal value of the binary variables at the
            # current interation
            self._yLBopt_curr = {}
            self._yUBopt_curr = {}
            for rxn in self.model.reactions:
                self._yLBopt_curr[rxn] = 1
                self._yUBopt_curr[rxn] = 1
            for rxn in zero_yLBopt_rxns:
                self._yLBopt_curr[rxn] = 0
            for rxn in zero_yUBopt_rxns:
                self._yUBopt_curr[rxn] = 0

        # If there was a solver error or if an optimal solution was not returned 
        else:
            if solverFlag == 'solverError':
                exit_flag = solverFlag
            else:
                exit_flag = str(OptSoln.solver.termination_condition)

            opt_objValue = None 
            zero_yLBopt_rxns = [] 
            zero_yUBopt_rxns = [] 

            if self.stdout_msgs:
                print "\n\n** No optimal solutions found (solution.solver.status = ",OptSoln.Solution.status,", solver.status =",OptSoln.solver.status,", solver.termination_condition = ",OptSoln.solver.termination_condition,")\n"

        #---- Solution ----
        self.solution.append({'exit_flag':exit_flag,'objective_value':opt_objValue,'modifications_num':self._max_modification_num_curr,'zero_yLBopt_rxns':zero_yLBopt_rxns,'zero_yUBopt_rxns':zero_yUBopt_rxns})

        if self.stdout_msgs:
           print 'Time elapsed (processing/wall) time (sec): Preprocess pyomo model = {:.3f}/{:.3f}  ,  solver = {:.3f}/{:.3f}\n'.format(elapsed_preproc_pyomo_pt,elapsed_preproc_pyomo_wt,elapsed_solver_pt,elapsed_solver_wt)


    def run(self, optModel_name = 'bilevel', max_modification_num = None, max_solution_num = None, max_biomass_thr = None):
        """ 
        This method runs FBA. 

        INPUTS:
        -------
               optModel_name: Name of the optimizaiton model to solve
        max_modification_num: Same as that for the constructor
             max_biomass_thr: Same as that for the constructor
        OUTPUT:
        -------
        solution: A dictionary with the following keys:
                        exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                                   or what is stored in OptSoln.solver.termination_condition
                  objective_value: Optimal objective funtion value
                 zero_yLBopt_rxns: List of reactions for which the optimal value of yLB is zero
                 zero_yUBopt_rxns: List of reactions for which the optimal value of yUB is zero
        Optimal values of binary variables are stored in two field yLB and yUB for each reaction 
        objects. The value of these fields will be None if the problem is not solved to optimality. 
        """
        if max_modification_num != None:
            self.max_modification_num = max_modification_num
 
        if max_solution_num != None:
            self.max_solution_num = max_solution_num

        if max_biomass_thr != None:
            self.max_biomass_thr = max_biomass_thr
        
        # Total processing and wall time required to create the pyomo model, solve it and store the results 
        start_total_pt = time.clock()
        start_total_wt = time.time()

        self.solution = []

        # The current allowed number of modifications
        self._max_modification_num_curr = 1

        # Number of solutions found so far
        self._found_solutions_num = 0

        #---- Creating the pyomo optModel ----
        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()
        # Create the pyomo model optModel only if self.build_new_optModel == 1        
        if optModel_name.lower() == 'bilevel' and self.build_new_optModel:
            self.build_bilevel_optModel()
        elif optModel_name.lower() == 'primal' and self.build_new_optModel:
            self.build_primal_optModel()
        elif optModel_name.lower() == 'dual' and self.build_new_optModel:
            self.build_dual_optModel()
        else:
            raise ValueError("Invalid optModel_name value. Allowed choices are 'bilevel', 'primal' and 'dual'")
        # Time to create a pyomo model
        elapsed_create_pyomo_pt = time.clock() - start_pyomo_pt
        elapsed_create_pyomo_wt = time.time() - start_pyomo_wt

        # A list of dictionaries holding the optimal solution in different iterations
        while self._max_modification_num_curr <= self.max_modification_num and self._found_solutions_num < self.max_solution_num:
            if self.stdout_msgs:
                print '\n---- Total # of modifications = {}, # of solutions found so far = {}\n'.format(self._max_modification_num_curr,self._found_solutions_num)

            # Solve the optimizaiton model with the current maximum number of allowed modifications 
            self.run_single() 

            # Check whether the current run was successful
            if self.solution[-1]['exit_flag'] == 'globallyOptimal':
                # Validate the solution
                self.validate_soln()

                self._found_solutions_num += 1

                # Add integer cut
                self.optModel.integer_cuts.add(self.integer_cuts_rule(optModel = self.optModel)) 
            else:
                self._max_modification_num_curr += 1

        # Time required to perform FBA
        elapsed_total_pt = (time.clock() - start_total_pt)
        elapsed_total_wt = (time.time() - start_total_wt)

        if self.stdout_msgs:
           print 'FBA elapsed (processing/wall) time (sec): create pyomo model = {:.3f}/{:.3f}  ,  total = {:.3f}/{:.3f}\n'.format(elapsed_create_pyomo_pt,elapsed_create_pyomo_wt,elapsed_total_pt,elapsed_total_wt)

        # Sort the solution according to the total penalites associated with modifications (i.e., objective funciton value)
        self.solution = sorted(self.solution,key = lambda x:(x['objective_value'] is None,x['objective_value']))
        return self.solution

#----------- Sample implementation -----------------
if __name__ == "__main__":
    pass
