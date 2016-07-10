from coopr.pyomo import *
from coopr.opt import *
from coopr.environ import *
from userError import userError
from globalVariables import *

def pyomoSolverCreator(optSolverName, **solver_options):
    """
    Creates a pyomo solver object and assigns the solver options

    INPUTS:
    ------
     optSolverName: A string containing the solver name
    solver_options: Additional solver options. These are to a dictionary whose keys are 
                    are the names of the arguments and values are the values 

    OUTPUTS:
    -------
    pymoSolverObject: Pyomo solver object with solver options assigned
    """    
    # Pyomo solver object
    pymoSolverObject = SolverFactory(optSolverName)

    # - Set some default solver options -        
    if optSolverName.lower() == 'cplex':
        # Memory
        pymoSolverObject.options["workmem"]=2500

        # Feasbility tolerance (eprhs). Defaul = 1e-6
        pymoSolverObject.options["simplex_tolerances_feasibility"]=1e-9

        # Optimality tolerance (epopt). Default = 1e-6
        pymoSolverObject.options["simplex_tolerances_optimality"]=1e-9

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["mip_tolerances_integrality"] = mip_integrality_tol

        # Relative MIP optimality gap, Default = 1e-4
        pymoSolverObject.options["mip_tolerances_mipgap"] = 1e-9

        # Absolute MIP optimality gap, Default = 1e-6
        pymoSolverObject.options["mip_tolerances_absmipgap"] = 1e-10

        # MIP strategy variable select (varsel). Default = 0
        pymoSolverObject.options["mip_strategy_variableselect"]=3

        # Bound strengthening indicator (bndstrenind). Default = 0   
        pymoSolverObject.options["preprocessing_boundstrength"]=1

    elif optSolverName.lower() in 'gurobi':
        # Memory (in Gb). Default: Infinity
        pymoSolverObject.options["NodefileStart"]=1

        # Feasbility tolerance. Defaul = 1e-6
        pymoSolverObject.options["FeasibilityTol"]=1e-9

        # Optimality tolerance. Defaul = 1e-6
        pymoSolverObject.options["OptimalityTol"]=1e-9

        # Relative MIP optimality gap, Default = 1e-4
        pymoSolverObject.options["MIPGap"]=1e-9

        # Absolute MIP optimality gap, Default = 1e-10
        pymoSolverObject.options["MIPGapAbs"]=1e-10

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["IntFeasTol"]= mip_integrality_tol

        # Branch variable selection strategy . Default = -1 (automatic)
        pymoSolverObject.options["VarBranch"]=3

        # Number of cores to use for parallel computations (defaul 0 = all cores) 
        pymoSolverObject.options["Threads"] = 4

    elif optSolverName.lower() in 'gurobi_ampl':
        # Memory (in Gb). Default: Infinity
        pymoSolverObject.options["NodefileStart"]=1

        # Feasbility tolerance. Defaul = 1e-6
        pymoSolverObject.options["feastol"]=1e-9

        # Optimality tolerance. Defaul = 1e-6
        pymoSolverObject.options["opttol"]=1e-9

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["IntFeasTol"]= mip_integrality_tol

        # Branch variable selection strategy . Default = -1 (automatic)
        pymoSolverObject.options["VarBranch"]=3

    else:
        raise userError('**Error! Invalid optimization solver name ...')


    #--- Set provided solver options by the user (these overwrfite tho defaults)
    for option_name in solver_options.keys():
        exec 'pymoSolverObject.options["' + option_name+ '"] = ' + str(solver_options[option_name])


    return pymoSolverObject 
