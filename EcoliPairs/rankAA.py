from __future__ import division
import sys
sys.path.append('../')
from copy import deepcopy
from tools.customError import *
from tools.FBA.fba import *
from tools.FBA.moma import *
from models import *
from metabolicModel import *

#-------------------- Sample implementation -----------------
if __name__ == "__main__":

    # Import the list of exchange reactions for amino acids
    from AAList import *

    #---  E. coli auxotrophic mutants using iAF1260 model ----
    # Import general model files
    import iAF1260ModelData

    # Create a general model 
    print 'Creating the metabolic model ...\n'
    iAF1260Model = metabolicModel(microorganism_name = 'E.coli', model_name = 'iAF1260', metab_names = iAF1260ModelData.metab_names,rxn_names = iAF1260ModelData.rxn_names, stoic_matrix = iAF1260ModelData.stoic_matrix, rxn_types = iAF1260ModelData.rxn_types, biomass_rxn_name = iAF1260ModelData.biomass_rxn_name)

    # Import the medium files
    import iAF1260_minimalGlucose_aerobic
    for rxn in iAF1260_minimalGlucose_aerobic.mediumSpecific_fluxBounds.keys():
        iAF1260Model.flux_bounds[rxn] = iAF1260_minimalGlucose_aerobic.mediumSpecific_fluxBounds[rxn]

    # Objective vector
    obj_vec = dict((rxn,0) for rxn in iAF1260Model.rxn_names)
    obj_vec[iAF1260ModelData.biomass_rxn_name] = 1
    iAF1260Model.obj_vec = obj_vec

    optSolverName = 'gurobi'

    # FBA for the wildtype 
    iAF1260FBA = fba(iAF1260Model, optSolverName = optSolverName, screenOutput = 'off')
    [FBAexitFlag,bmWT] = iAF1260FBA.run()
    print 'WT      ',bmWT

    AACost = []
    if FBAexitFlag == 'globallyOptimal':
        for AAexch in AAList:
            AA_FBA = deepcopy(iAF1260FBA)
            AA_FBA.createModel = 0
            AA_FBA.fbaModel.v[AAexch].lb = 1
            [exitFlag,bmAA] = AA_FBA.run()
            if exitFlag == 'globallyOptimal':
                AACost.append((AAexch,bmWT - bmAA))
                print AAexch,'   bmAA = ',bmAA,'     (bmWT - bmAA) = ',bmWT - bmAA
            else:
                AACost.append((AAexch,None))
                print AAexch,'   bmAA = None (exitflag = ',exitFlag,')'
    else:
        raise customError('**Error! No optimal solution could be found for the wildtype strain ...') 

    AACost.sort(key=lambda t: t[1],reverse = True)
    print '\nAmino acids costs = '
    for AA in AACost:
        print AA 

    print '\n---- Test EX_met-L(e) ----'     
    metFBA = deepcopy(iAF1260FBA)
    metFBA.metabolicModel.obj_vec = dict((rxn,0) for rxn in iAF1260Model.rxn_names)
    metFBA.metabolicModel.obj_vec['EX_trp-L(e)'] = 1
    metFBA.fbaModel.del_component('objectiveFunc')
    metFBA.fbaModel.objectiveFunc = Objective(rule=metFBA.objectiveFunc_rule, sense = minimize)
    [exitFlag,metMax] = metFBA.run()
    if exitFlag == 'globallyOptimal':
        print 'EX_met-L(e) max = ',metMax
    else:
        print 'exitflag = ',exitflag

