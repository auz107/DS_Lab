from __future__ import division
import sys, time
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.pyomoSolverCreator import *
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tool.ancillary.get_ModelSEED_ids import get_ModelSEED_ids

# The following lines change the temporary directory for pyomo
#from pyutilib.services import TempfileManager
#TempfileManager.tempdir = '/data/alizom/tmp/'

class gapfill(object):
    """
    Performs gap filling for a metabolic model  

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 02-18-2016 
    """   

    def __init__(self,model, gapfilling_method = 'parsimony_based', penalties = {},optimization_solver = 'gurobi', build_new_optModel = True, stdout_msgs = True, warnings = True,  **additional_args): 
        """
        INPUTS (required):
        ------
                      model: An instance of class model containing the model to be gap filled 

        INPUTS (optional):
        ------
          gapfilling_method: The method to perform gap filling. The default is parsimony-based gap filling, where the
                             aim is to minimally perturb the model.
                  penalties: Penalties associated with performing each modification to the model
        optimization_solver: Name of the LP solver to be used to solve the LP. Current 
                             allowable choices are cplex and gurobi
         build_new_optModel: If True a new pyomo optimization model is constructed
                stdout_msgs: If True writes a summary of the run in the output 
                   warnings: If True prints the warnings in the output 
            additional_args: Additional arguments should be entered as normal but they are 
                             converted to a dictionary whose keys are the names of the arguments and 
                             values are the values of  those arguments

        OUTPUTS:
        ---------
        solution: A dictionary with the following keys:
                  exit_flag: A string, which can be 'globallyOptimal', 'solverError'
                            or what is stored in OptSoln.solver.termination_condition
                  objValue: Optimal objective funtion value

        These are the outputs of the method 'run'
        """
       
        # Metabolic model
        self.model = model

        # Gap filling method
        self.gapfilling_method = gapfilling_method

        # Solver name
        self.optimization_solver = optimization_solver

        # Whether to create a pyomo model
        self.build_new_optModel = build_new_optModel

        # Penalties associated with adding each reaction type
        self.penalties = penalties
               
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

        # Get ModelSEED ids for compounds and reactions
        get_ModelSEED_ids(model = self.model, msgs = True)

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
       elif attr_name.lower() == 'gapfilling_method' and attr_value.lower() not in ['cplex','gurobi']: 
           raise ValueError('Invalid gapfilling_method value! Allowed choices are: [cplex,gurobi]')

       if attr_name.lower() == 'penalties' and not isinstance(attr_value,dict):
           raise TypeError('penalties must be a dictionary')
       elif attr_name.lower() == 'penalties' and len(attr_value.keys)) > 0 and len([k for k in attr_value.keys() if v not in ['exchange_rxns']]) > 0:
           raise ValuesError('Invalid key for penalties: {}. Allowed keys are [exchange_rxns]'.format([k for k in attr_value.keys() if v not in ['exchange_rxns']]))

       if attr_name.lower() == 'build_new_optModel' and not isinstance(attr_value,bool):
           raise TypeError('build_new_optModel must be either True or False')

       if attr_name.lower() == 'stdout_msgs' and not isinstance(attr_value,bool):
           raise TypeError('stdout_msgs must be either True or False')

       if attr_name.lower() == 'warnings' and not isinstance(attr_value,bool):
           raise TypeError('warnings must be either True or False')

       self.__dict__[attr_name] = attr_value

    def create_exchange_transport(self):
        """
        Creates exchange and transport reactions not existing in the model
        """

    def objectiveFunc_rule(self,pyomo_fbaModel):
        """
        Objective function
        """
        # Reactions for which the objective coefficient has not bee assigned
        non_obj_rxns = [j.id for j in pyomo_fbaModel.J if j.objective_coefficient == None]
        if len(non_obj_rxns) >= 1: 
            print("**ERROR! 'objective_coefficient' has not been defined for the following reactions:")
            print non_obj_rxns
            raise userError()
        return sum(j.objective_coefficient*pyomo_fbaModel.v[j] for j in pyomo_fbaModel.J)

    def assignFluxBounds(self,pyomo_fbaModel,j):
        """
        Define the flux bounds
        """
        return j.flux_bounds 
        
    def massBalance_rule(self,pyomo_fbaModel,i):
        """
        Mass balance 
        """
        return sum(j.stoichiometry[i]*pyomo_fbaModel.v[j] for j in i.reactions) == 0 
        
    def createPyomoModel(self):
        """
        This pyomo_fbaModel creates a pyomo model for FBA pyomo_fbaModel
        """   
        #--- Create a pyomo model pyomo_fbaModel ---
        pyomo_fbaModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of compounds 
        pyomo_fbaModel.I = Set(initialize = self.model.compounds)   

        # Set of rxns  
        pyomo_fbaModel.J = Set(initialize = self.model.reactions)     

        #--- Define the pyomo_fbaModel variables --- 
        pyomo_fbaModel.v = Var(pyomo_fbaModel.J, domain=Reals, bounds = self.assignFluxBounds)
        
        #--- Defiine the objective function and constraints ----
         # Objective function
        pyomo_fbaModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = maximize)

        # Mass balance 
        pyomo_fbaModel.massBalance_const = Constraint(pyomo_fbaModel.I, rule=self.massBalance_rule)

        self.pyomo_fbaModel = pyomo_fbaModel 
    
    
    def run(self):
        """ 
        This method runs FBA. 

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

        # Processing and walltime
        start_fba_pt = time.clock()
        start_fba_wt = time.time()

        #---- Creating and instantiating the pyomo_fbaModel ----
        start_pyomo_pt = time.clock()
        start_pyomo_wt = time.time()

        # Create the pyomo model pyomo_fbaModel only if self.build_new_optModel == 1        
        if self.build_new_optModel == True:
            self.createPyomoModel()

        # Instantiate the pyomo_fbaModel
        self.pyomo_fbaModel.preprocess()

        #---- Solve the model ----
        # Create a solver and set the options
        solverType = pyomoSolverCreator(self.optimization_solver)

        elapsed_pyomo_pt = (time.clock() - start_pyomo_pt)
        elapsed_pyomo_wt = (time.time() - start_pyomo_wt)

        #- Solve the pyomo_fbaModel (tee=True shows the solver output) -
        try:
            start_solver_pt = time.clock()
            start_solver_wt = time.time()
            OptSoln = solverType.solve(self.pyomo_fbaModel,tee=False)
            solverFlag = 'normal'

        # In the case of an error switch the solver
        except:
            if self.stdout_msgs:
                print "**Warning! ",self.optimization_solver," failed. An alternative solver is tried"        

            if self.optimization_solver.lower() == 'gurobi':
                self.optimization_solver = 'cplex'
            elif self.optimization_solver.lower() == 'cplex':
                self.optimization_solver = 'gurobi'

            # Try solving with the alternative solver
            solverType = pyomoSolverCreator(self.optimization_solver)
            try:
                start_solver_pt = time.clock()
                start_solver_wt = time.time()
                OptSoln = solverType.solve(self.pyomo_fbaModel,tee=False)
                solverFlag = 'normal'
            except:
                solverFlag = 'solverError'
                if self.stdout_msgs:
                    print '**Warning! The alternative solver failed. No solution was returned'

        elapsed_solver_pt = (time.clock() - start_solver_pt)
        elapsed_solver_wt = (time.time() - start_solver_wt)

        #----- Print the results in the output ------
        if solverFlag == 'normal' and str(OptSoln.solver.termination_condition).lower() == 'optimal':
        
            exit_flag = 'globallyOptimal'

            # Load the results
            self.pyomo_fbaModel.load(OptSoln)
        
            # Value of the objective function
            objValue = self.pyomo_fbaModel.objectiveFunc()

            # Optimal values of variables
            optVarValues = {}

            # Print the results on the screen 
            if self.stdout_msgs:
                print "\nSolver.status = ",OptSoln.solver.termination_condition
                print "Optimality status = ",exit_flag
                print "Objective value = ",objValue

            # Store the optimal flux values in the variable 'flux' of the reaction objects
            for rxn in [r for r in self.model.reactions if r.store_flux == True]:
                if self.flux_key is None:
                    rxn.flux = self.pyomo_fbaModel.v[rxn].value
                elif self.flux_key is not None and type(rxn.flux) is dict:
                    rxn.flux[self.flux_key] = self.pyomo_fbaModel.v[rxn].value
                else:
                    rxn.flux = {}
                    rxn.flux[self.flux_key] = self.pyomo_fbaModel.v[rxn].value

            self.solution = {'exit_flag':exit_flag,'objective_value':objValue}

        # If there was a solver error or if an optimal solution was not returned 
        else:
            if solverFlag == 'solverError':
                exit_flag = solverFlag
            else:
                exit_flag = str(OptSoln.solver.termination_condition)

            objValue = None 
            optVarValues = None 

            if self.stdout_msgs:
                print "\n\n** No optimal solutions found (solution.solver.status = ",OptSoln.Solution.status,", solver.status =",OptSoln.solver.status,", solver.termination_condition = ",OptSoln.solver.termination_condition,")\n"

            self.solution = {'exit_flag':exit_flag,'objective_value':objValue}

            for rxn in [r for r in self.model.reactions if r.store_flux == True]:
                if self.flux_key is None and type(rxn.flux) is not dict:
                    rxn.flux = None 
                elif self.flux_key is not None and type(rxn.flux) is dict:
                    rxn.flux[self.flux_key] = None 
                else:
                    if self.warnings:
                        print 'WARNING! The current reaction.flux is not in the formm of a dictionary for reaction ',rxn.id,' and will be overwritten'
                    rxn.flux = {}
                    rxn.flux[self.flux_key] = None 

        # Time required to perform FBA
        elapsed_fba_pt = (time.clock() - start_fba_pt)
        elapsed_fba_wt = (time.time() - start_fba_wt)

        if self.stdout_msgs:
           print 'FBA elapsed (processing/wall) time (sec): pyomo = {:.3f}/{:.3f}  ,  solver = {:.3f}/{:.3f}  ,  fba = {:.3f}/{:.3f}\n'.format(elapsed_pyomo_pt,elapsed_pyomo_wt,elapsed_solver_pt,elapsed_solver_wt,elapsed_fba_pt,elapsed_fba_wt)

        return self.solution

#----------- Sample implementation -----------------
if __name__ == "__main__":

    import time
    from tools.io.read_gams_model import read_gams_model
    from tools.io.read_sbml_model import read_sbml_model
    from set_specific_bounds import set_specific_bounds
    from cobra import test
 
    # Solver name
    optimization_solver = 'gurobi'

    #--- Test model ---
    print '\n--- Test model ---'
    testModel = read_gams_model(gams_model_file = '/fs/home06/alizom//models/test/testModelData.py',model_name = 'testModel',organism_name = 'testOrg',model_type = 'metabolic')

    # Growth medium
    #testModel = set_specific_bounds(testModel,file_name = '/fs/home06/alizom/models/test/testMedium.py')
    set_specific_bounds(testModel,file_name = '/fs/home06/alizom/models/test/testMedium.py')

    # Assign and objective function coefficients
    for rxn in testModel.reactions:
        rxn.objective_coefficient = 0

    for bm in testModel.biomass_reaction:
        bm.objective_coefficient = 1 
 
    fbaTest = fba(testModel, optimization_solver = optimization_solver) 
    solution = fbaTest.run()
    for r in testModel.reactions:
        print r.id,'   Vopt = ',r.flux 
      

    #--- E. coli iAF1260 model ---
    print '\n--- iAF1260 model ---'
    print '   Read the gams model ...'
    start = time.clock()
    iAF1260 = read_gams_model(gams_model_file = '/fs/home06/alizom//models/Ecoli/iAF1260/iAF1260ModelData.py',model_name = 'iAF1260',organism_name = 'E. coli',model_type = 'metabolic')
    print '        Reading the gams model took ',str(time.clock() - start)

    # Growth medium
    print '   Set the growth meidum ...'
    set_specific_bounds(iAF1260,file_name = '/data/alizom/models/Escherichia_coli/iAF1260/iAF1260_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')
    #set_specific_bounds(iAF1260,file_name = '/data/alizom/models/Escherichia_coli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign and objective function coefficients
    for rxn in iAF1260.reactions:
        rxn.objective_coefficient = 0

    for bm in iAF1260.biomass_reaction:
        if bm.id == 'Ec_biomass_iAF1260_core_59p81M':
            bm.objective_coefficient = 1 
            objective_rxn = bm 

    print '   create the pyomo model ...'
    fbaiAF1260 = fba(iAF1260, optimization_solver = optimization_solver) 
    fbaiAF1260.run()
    print 'flux value of biomass = ',objective_rxn.flux

    for rxn in [r for r in iAF1260.reactions if r.id in ['EX_glc(e)','EX_o2(e)']]:
        print rxn.id,'    ',rxn.flux_bounds,'    ',rxn.flux

    """
    #--- Salmonella ---
    print '\n--- Salmonella model ---'
    print '   Read the sbml model ...'
    salModel = read_sbml_model(file_name = test.salmonella_sbml,model_name = 'salmonella',model_organism = 'salmonella', model_type = 'metabolic',import_bounds = 1)
 
    # Assign and objective function coefficients
    for rxn in iAF1260.reactions:
        rxn.objective_coefficient = 0

    for bm in iAF1260.biomass_reaction:
        if bm.id == 'Ec_biomass_iAF1260_core_59p81M':
            bm.objective_coefficient = 1 
            objective_rxn = bm
 
    print '   create the pyomo model ...'
    fbaSal = fba(iAF1260, optimization_solver = optimization_solver) 
    fbaSal.run()
    print 'flux value of biomass = ',objective_rxn.flux

    """


