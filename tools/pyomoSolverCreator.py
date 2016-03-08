from coopr.pyomo import *
from coopr.opt import *
from coopr.environ import *
from userError import userError

def pyomoSolverCreator(optSolverName):
    """
    Creates a pyomo solver object and assigns the solver options

    INPUTS:
    ------
       optSolverName: A string containing the solver name

    OUTPUTS:
    -------
    pymoSolverObject: Pyomo solver object with solver options assigned
    """    

    # Pyomo solver object
    pymoSolverObject = SolverFactory(optSolverName)

    # - Set the solver options -        
    if optSolverName.lower() == 'cplex':
        # Memory
        pymoSolverObject.options["workmem"]=2500

        # Feasbility tolerance (eprhs). Defaul = 1e-6
        pymoSolverObject.options["simplex_tolerances_feasibility"]=1e-9

        # Optimality tolerance (epopt). Default = 1e-6
        pymoSolverObject.options["simplex_tolerances_optimality"]=1e-9

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["mip_tolerances_integrality"]=1e-9

        # MIP strategy variable select (varsel). Default = 0
        pymoSolverObject.options["mip_strategy_variableselect"]=3

        # Bound strengthening indicator (bndstrenind). Default = 0   
        pymoSolverObject.options["preprocessing_boundstrength"]=1

    elif optSolverName.lower() == 'gurobi':
        # Memory (in Gb). Default: Infinity
        pymoSolverObject.options["NodefileStart"]=1

        # Feasbility tolerance. Defaul = 1e-6
        pymoSolverObject.options["FeasibilityTol"]=1e-9

        # Optimality tolerance. Defaul = 1e-6
        pymoSolverObject.options["OptimalityTol"]=1e-9

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["IntFeasTol"]=1e-8

        # Branch variable selection strategy . Default = -1 (automatic)
        pymoSolverObject.options["VarBranch"]=3

    else:
        raise customError('**Error! Invalid optimization solver name ...')

    return pymoSolverObject 
