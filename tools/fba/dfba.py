from __future__ import division
import re, sys, math, copy, time, random
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.userError import *
from tools.pyomoSolverCreator import *
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model

class dfba(object):
    """
    Performs flux balance analysis

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 10-30-2014 
    """   

    def __init__(self,model, optSolverName = None, createModel = None, screenOutput = None): 

        """
        INPUTS (required):
        ------
                      model: An instance of class model containing the information
                             about the metabolic model

        INPUTS (optional):
        ------
              optSolverName: Name of the LP solver to be used to solve the LP. Current 
                             allowable choices are cplex and gurobi
               screenOutput: By default (on) writes  a summary including the solve 
                             status, optimality status (if not optimal), objective 
                             function value and the elapsed time on the screen.
                             if set to a value of 'off' no resuults are written on 
                             the screen, in which case The user can instead specifiy 
                             an output fiile using the option outputFile, or store 
                             them in the variable runOutput (see the 'run' method for
                             details)
                createModel: A parameter indicating whether a pyomo model should be 
                             created (1) or not. Default is 1. The options is useful
                             for the cases a model is already created and one just 
                             wants to change some model attributes (e.g., flux bounds)
                             and rerun FBA. Setting this parameter to zero will save 
                             soem runtime as the model need not to be created again.
        """
       
        # Metabolic model
        self.model = model

        # Solver name
        if optSolverName == None:
            self.optSolverName = 'cplex'
        else:
            if optSolverName.lower() in ['cplex','gurobi']:
                self.optSolverName = optSolverName
            else:
                raise userError('**Error! Invalid solver name (eligible choices are cplex and gurobi)\n')          

        # Whether to create a pyomo model
        if createModel == None: 
            self.createModel = 1   
        elif createModel not in [0,1]:    
            raise userError('**Error! Inalid value for createModel! The allowable values are 0 and 1 ')
        else:
            self.createModel = createModel
               
        # Output to the screen 
        if screenOutput == None:
            self.screenOutput = 'on'
        elif type(screenOutput) is not str:
            raise userError("**Error! screenOutput should be a string ('on' or 'off')")
        elif screenOutput.lower() not in ['on','off']:
            raise userError("**Error! The only eligible values for screenOutput are 'on' and 'off'")
        else:
             self.screenOutput = screenOutput

    # Objective function
    def objectiveFunc_rule(self,fbaModel):
        # Reactions for which the objective coefficient has not bee assigned
        non_obj_rxns = [j.id for j in fbaModel.J if j.objective_coefficient == None]
        if len(non_obj_rxns) >= 1: 
            print("**Error! 'objective_coefficient' has not been defined for the following reacitons:")
            print non_obj_rxns
            raise userError()
        return sum(j.objective_coefficient*fbaModel.v[j] for j in fbaModel.J)
        

    # Mass balance 
    def massBalance_rule(self,fbaModel, i):
        return sum(j.stoichiometry[i]*fbaModel.v[j] for j in i.reactions) == 0 
        
    def createPyomoModel(self):
        """
        This fbaModel creates a pyomo model for FBA fbaModel
        """   
        #--- Create a pyomo model fbaModel ---
        fbaModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of compounds 
        fbaModel.I = Set(initialize = self.model.compounds)   

        # Set of rxns  
        fbaModel.J = Set(initialize = self.model.reactions)     

        #--- Define the fbaModel variables --- 
        def assignFluxBounds(fbaModel,j):
            return j.flux_bounds 
        
        fbaModel.v = Var(fbaModel.J, domain=Reals, bounds = assignFluxBounds)
        
        #--- Defiine the objective function and constraints ----
         # Objective function
        fbaModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = maximize)
        

        # Mass balance 
        fbaModel.massBalance_const = Constraint(fbaModel.I, rule=self.massBalance_rule)

        self.fbaModel = fbaModel 
    
    
    def run(self):
        """ 
        This method runs FBA. 
        OUTPUT:
        -------
        runOutput: A list of three items:
                       exitflag: A string, which can be 'globallyOptimal', 'solverError'
                                 or what is stored in OptSoln.solver.termination_condition
                       objValue: Optimal objective funtion value

        NOTE: Optimal flux values are stored directly into the 'flux' field of the reaciton
              objects. The value of the flux will be None if the problem is not solved to 
              optimality. 
        """

        start_fba = time.clock()

        #---- Creating and instantiating the fbaModel ----
        start_pyomo = time.clock()

        # Create the pyomo model fbaModel only if self.createModel == 1        
        if self.createModel == 1:
            self.createPyomoModel()

        # Instantiate the fbaModel
        self.fbaModel.preprocess()

        #---- Solve the model ----
        # Create a solver and set the options
        solverType = pyomoSolverCreator(self.optSolverName)

        elapsed_pyomo = (time.clock() - start_pyomo)

        #- Solve the fbaModel (tee=True shows the solver output) -

        try:
            start_solver = time.clock()
            OptSoln = solverType.solve(self.fbaModel,tee=False)
            solverFlag = 'normal'

        # In the case of an error switch the solver
        except:
            if self.screenOutput.lower() == 'on':
                print "**Warning! ",self.optSolverName," failed. An alternative solver is tried"        

            if self.optSolverName.lower() == 'gurobi':
                self.optSolverName = 'cplex'
            elif self.optSolverName.lower() == 'cplex':
                self.optSolverName = 'gurobi'

            # Try solving with the alternative solver
            solverType = pyomoSolverCreator(self.optSolverName)
            try:
                start_solver = time.clock()
                OptSoln = solverType.solve(self.fbaModel,tee=False)
                solverFlag = 'normal'
            except:
                solverFlag = 'solverError'
                if self.screenOutput.lower() == 'on':
                    print '**Warning! The alternative solver failed. No solution was returned'

        elapsed_solver = (time.clock() - start_solver)

        #----- Print the results in the output ------
        if solverFlag == 'normal' and str(OptSoln.solver.termination_condition).lower() == 'optimal':
        
            exitflag = 'globallyOptimal'

            # Load the results
            self.fbaModel.load(OptSoln)
        
            # Value of the objective function
            objValue = self.fbaModel.objectiveFunc()

            # Optimal values of variables
            optVarValues = {}

            # Print the results on the screen 
            if self.screenOutput.lower() == 'on':
                print "\nSolver.status = ",OptSoln.solver.termination_condition
                print "Optimality status = ",exitflag
                print "Objective value = ",objValue
            elif self.screenOutput.lower() == 'off':
                pass

            # Store the optimal flux values in the field 'flux' of the reaction objects
            for rxn in self.model.reactions:
                rxn.flux = self.fbaModel.v[rxn].value

            runOutput = {'exitflag':exitflag,'objective_value':objValue,'model':self.model}


        # If there was a solver error or if an optimal solution was not returned 
        else:
            if solverFlag == 'solverError':
                exitflag = solverFlag
            else:
                exitflag = str(OptSoln.solver.termination_condition)

            objValue = None 
            optVarValues = None 

            if self.screenOutput.lower() == 'on':
                print "\n\n** No optimal solutions found (solution.solver.status = ",OptSoln.Solution.status,", solver.status =",OptSoln.solver.status,", solver.termination_condition = ",OptSoln.solver.termination_condition,")\n"

            runOutput = {'exitflag':exitflag,'objective_value':objValue,'model':self.model}

        # Time required to perform FBA
        elapsed_fba = (time.clock() - start_fba)

        if self.screenOutput.lower() == 'on':
           print 'elapsed time: (pyomo = ',elapsed_pyomo,'  ,  solver = ',elapsed_solver , '  ,  fba = ',elapsed_fba,')\n'

        return runOutput

#----------------------------
if __name__ == "__main__":

    import time
    from tools.io.read_gams_model import read_gams_model
    from tools.io.read_sbml_model import read_sbml_model
    from set_specific_bounds import set_specific_bounds
    from cobra import test
 
    # Solver name
    optSolverName = 'gurobi'

    #--- Test model ---
    print '\n--- Test model ---'
    testModel = read_gams_model(gams_model_file = '/fs/home06/alizom//models/test/testModelData.py',model_name = 'testModel',organism_name = 'testOrg',model_type = 'metabolic')

    # Growth medium
    #testModel = set_specific_bounds(testModel,specific_bounds_file = '/fs/home06/alizom/models/test/testMedium.py')
    set_specific_bounds(testModel,specific_bounds_file = '/fs/home06/alizom/models/test/testMedium.py')

    # Assign and objective function coefficients
    for rxn in testModel.reactions:
        rxn.objective_coefficient = 0

    for bm in testModel.biomass_reactions:
        bm.objective_coefficient = 1 
 
    fbaTest = fba(testModel, optSolverName = optSolverName) 
    runOutput = fbaTest.run()
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
    #iAF1260 = set_specific_bounds(iAF1260,specific_bounds_file = '/fs/home06/alizom/models/Ecoli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')
    set_specific_bounds(iAF1260,specific_bounds_file = '/fs/home06/alizom/models/Ecoli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign and objective function coefficients
    for rxn in iAF1260.reactions:
        rxn.objective_coefficient = 0

    for bm in iAF1260.biomass_reactions:
        if bm.id == 'Ec_biomass_iAF1260_core_59p81M':
            bm.objective_coefficient = 1 
            objective_rxn = bm 

    print '   create the pyomo model ...'
    fbaiAF1260 = fba(iAF1260, optSolverName = optSolverName) 
    fbaiAF1260.run()
    print 'flux value of biomass = ',objective_rxn.flux

    #--- Salmonella ---
    print '\n--- Salmonella model ---'
    print '   Read the sbml model ...'
    salModel = read_sbml_model(sbml_model_file = test.salmonella_sbml,model_name = 'salmonella',organism_name = 'salmonella', model_type = 'metabolic',import_bounds = 1)
 
    # Assign and objective function coefficients
    for rxn in iAF1260.reactions:
        rxn.objective_coefficient = 0

    for bm in iAF1260.biomass_reactions:
        if bm.id == 'Ec_biomass_iAF1260_core_59p81M':
            bm.objective_coefficient = 1 
            objective_rxn = bm
 
    print '   create the pyomo model ...'
    fbaSal = fba(iAF1260, optSolverName = optSolverName) 
    fbaSal.run()
    print 'flux value of biomass = ',objective_rxn.flux



