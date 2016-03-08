from __future__ import division
import sys,time
sys.path.append('../')
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.customError import *
from copy import deepcopy
from models import *
from metabolicModel import *
from fba import *

# -------------------------------------------
# Tests the growth of each mutant in the presence/absence of selected metabolites
# -------------------------------------------

if __name__ == "__main__":

    
    #---  E. coli auxotrophic mutants using iAF1260 model ----
    # Import general model files
    import iAF1260ModelData
    
    # Create a general model 
    print 'Creating the metabolic model ...\n'
    iAF1260Model = metabolicModel(microorganism_name = 'E.coli', model_name = 'iAF1260', metab_names = iAF1260ModelData.metab_names,rxn_names = iAF1260ModelData.rxn_names, stoic_matrix = iAF1260ModelData.stoic_matrix, rxn_types = iAF1260ModelData.rxn_types, biomass_rxn_name = iAF1260ModelData.biomass_rxn_name)
    
    # Import the medium files
    import iAF1260_minimalGlucose_anaerobic
    for rxn in iAF1260_minimalGlucose_anaerobic.mediumSpecific_fluxBounds.keys():
        iAF1260Model.flux_bounds[rxn] = iAF1260_minimalGlucose_anaerobic.mediumSpecific_fluxBounds[rxn]
    
    # Objective vector
    obj_vec = dict((rxn,0) for rxn in iAF1260Model.rxn_names)
    obj_vec[iAF1260ModelData.biomass_rxn_name] = 1
    iAF1260Model.obj_vec = obj_vec

    optSolverName = 'cplex'

    #--- Test the wild-type network ---
    print '------- Testing the wildtype ------\n'
    iAF1260FBA = fba(iAF1260Model, optSolverName = optSolverName)
    [optimExitFlag,objVal] = iAF1260FBA.run() 

    maxBiomass = objVal

    # mutant 1
    m1 = 'hisD'

    print '******** ',m1,' ***********\n'
    # Create the mutant model instance 
    #m1 = deepcopy(iAF1260Model)
    #m1.flux_bounds['HISTD'] = [0,0] 
    #hisDFBA = fba(hisD,optSolverName = optSolverName)

    # Create a copy of the wild-type FBA model for mutant m1
    m1FBA = deepcopy(iAF1260FBA)
    m1FBA.createModel = 0

    #m1FBA.fbaModel.v['HISTD'].setlb(0)
    #m1FBA.fbaModel.v['HISTD'].setlb(0)
    m1FBA.fbaModel.v['HISTD'] = 0
    m1FBA.fbaModel.v['HISTD'].fixed = True
    [optimExitFlag,objVal] = m1FBA.run() 

    print '******** Allow for uptake by ',m1,' ***********\n'
    m1FBA.fbaModel.v['EX_his-L(e)'].setlb(-1)
    [optimExitFlag,objVal] = m1FBA.run() 

