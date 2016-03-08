if __name__ == "__main__":

    import sys
    sys.path.append('../')
    import pickle as pk
    from copy import deepcopy
    from models import *
    from tools.FBA import *
    from metabolicModel import *
    from fba import *
    from mutants_rxn_info import *
    from auxoMetabsFinder import *

    #---  E. coli auxotrophic mutants using iAF1260 model ----
    # Import general model files
    from iAF1260ModelData import *
    
    # Create a general model 
    print 'Creating the metabolic model ...\n'
    iAF1260Model = metabolicModel(microorganism_name = 'E.coli', model_name = 'iAF1260', metab_names = metab_names,rxn_names = rxn_names, stoic_matrix = stoic_matrix, rxn_types = rxn_types, biomass_rxn_name = biomass_rxn_name)
    
    # Import the medium files (anaerobic minimal glucose medium)
    from iAF1260_minimalGlucose_anaerobic import *
    for rxn in mediumSpecific_fluxBounds.keys():
        iAF1260Model.flux_bounds[rxn] = mediumSpecific_fluxBounds[rxn]
    
    # Objective vector
    obj_vec = dict((rxn,0) for rxn in iAF1260Model.rxn_names)
    obj_vec[biomass_rxn_name] = 1
    iAF1260Model.obj_vec = obj_vec

    optSolverName = 'gurobi'

    #--- Test the wild-type network ---
    print '------- Testing the wildtype ------\n'
    iAF1260FBA = fba(iAF1260Model, optSolverName = optSolverName)
    [exitflag,objVal] = iAF1260FBA.run() 

    maxBiomass = objVal

    #----------------------------------
    # First find out which metabolites can be produced by wild-type E.coli under the given
    # medium condition. If a metabolite cannot be produced, it means that the cost of production
    # of that metabolite is infinity.
    exchFBA = deepcopy(iAF1260FBA)
    exchFBA.screenOutput = 'silent'
    exchFBA.createModel = 0
    infiniteCostMetabs = []
    for exchRxn in iAF1260Model.rxn_names:
        # Do not cosider exchange rxns for metabolites that are already in the growth medium
        if iAF1260Model.rxn_types[exchRxn] == 3 and iAF1260Model.flux_bounds[exchRxn][0] >= 0:
            print exchRxn
            exchFBA.metabolicModel.obj_vec = dict((rxn,0) for rxn in iAF1260Model.rxn_names)
            exchFBA.metabolicModel.obj_vec[exchRxn] = 1
            exchFBA.fbaModel.del_component('objectiveFunc')
            exchFBA.fbaModel.objectiveFunc = Objective(rule=exchFBA.objectiveFunc_rule, sense = maximize)
            [exitflag,objVal] = exchFBA.run()
            if exitflag == 'globallyOptimal' and objVal == 0:
                infiniteCostMetabs.append(exchRxn)

    # Save the results into a binary file on disk
    with open('results/infiniteCostMetabs.pk','wb') as outputFile:
        pk.dump(infiniteCostMetabs,outputFile,-1)

    #----------------------------------
    # The results are saved in a dictionary whose keys are the mutant names and values
    # are lists of lists containing the set of reactions needed for the survivval of
    # that strain. For example:
    # {'pykA':[['EX_m1','EX_m2'],['EX_m1',EX_m5'],'hisD':[['EX_m2'],['Ex_m4']]]}
    auxoMetabs = {} 

    # Termination condition for each mutant is also reported
    termCond_mutant = {} 

    for mutantName in mutants_rxn_info.keys():
        print mutantName
        current_mutantModel = deepcopy(iAF1260Model)
        for mutant_rxns in mutants_rxn_info[mutantName]:
           current_mutantModel.flux_bounds[mutant_rxns] = [0,0]

        current_mutantAuxo = auxoMetabsFinder(current_mutantModel, maxBiomass = maxBiomass, viabilityThr = 10, maxSolnSize = 5, optSolverName = optSolverName, screenOutput ='silent' ,reportLevel = {'y':'all'})
        [runTerminationCond,runOutput] = current_mutantAuxo.run()   

        termCond_mutant[mutantName] = runTerminationCond

        auxoMetabs[mutantName] = []
        for soln in range(len(runOutput)):
            optVarValues = runOutput[soln][2]
            rxnsList = [rxn for rxn in optVarValues['y'].keys() if optVarValues['y'][rxn] == 1 and rxn not in infiniteCostMetabs]
            if len(rxnsList) > 0:
                auxoMetabs[mutantName].append(rxnsList)


    # Save the results into a binary file on disk
    with open('results/auxoMetabs.pk','wb') as outputFile:
        pk.dump([termCond_mutant,auxoMetabs],outputFile,-1)

