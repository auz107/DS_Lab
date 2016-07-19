from __future__ import division
import sys, time
sys.path.append('../../')
from tools.globalVariables import *
from tools.fba.fbaTools import fbaTools
from coopr.pyomo import *
from coopr.opt import *

class moma(fbaTools):
    """
    Performs flux balance analysis

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 10-27-2015 
    """   

    def __init__(self,model, distance_type = '2-norm',optimization_solver = default_optim_solver, build_new_optModel = True, flux_key = None, store_opt_fluxes = True, warmstart = False, show_solver_output = False, warnings = True, stdout_msgs = True): 

        """
        INPUTS:
        ------
        distance_type: 
        The way to compute the distance between the fluxes.
        Eligible cases are 'euclidean' or '2-norm' (default) 
        and '1-norm'

        The rest of inputs are as those in fbaTools
        """
       
        # Metabolic model
        self.model = model

        # Type of distance
        self.distance_type = distance_type

        # Solver name
        self.optimization_solver = optimization_solver

        # Whether to create a pyomo model
        self.build_new_optModel = build_new_optModel
    
        # warmstart and show_solver_output
        self.warmstart = warmstart
        self.show_solver_output = show_solver_output
           
        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.warnings = warnings

        # flux_key
        self.flux_key = flux_key 

        # store_opt_fluxes
        self.store_opt_fluxes = store_opt_fluxes
        if self.store_opt_fluxes == False:
            for rxn in model.reactions:
                rxn.flux = None

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
       -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # model
        if attr_name == 'model': 
            # Check if the all reactions in the model have a field called 
            # flux_wildtype containing the wild-type reaction fluxes
            no_wt_flux_rxns = [r for r in attr_value.reactions if 'flux_wildType' not in dir(r)]
            if len(no_wt_flux_rxns) > 0 and len(no_wt_flux_rxns) <= 10:
                raise AttributeError('No flux_wildType attribute for ' + str(len(no_wt_flux_rxns)) + ' reactions: ' + str(no_wt_flux_rxns))
            elif len(no_wt_flux_rxns) > 0 and len(no_wt_flux_rxns) > 10:
                raise AttributeError('No flux_wildType attribute for ' + str(len(no_wt_flux_rxns)) + ' reactions including: ' + str(no_wt_flux_rxns[:10]) + '\nand ' + str(len(no_wt_flux_rxns) - 10) + ' more reactions.')

        # distance_type 
        if attr_name == 'distance_type' and attr_value.lower() not in ['euclidean','2-norm','1-norm']:
            raise ValueError('Invalid value for distance_type. Eligible inputs are euclidean or 2-norm and 2-norm')

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Warnings and messages in the standard output
        if attr_name in ['build_new_optModel', 'warnings','stdout_msgs','warmstart', 'show_solver_output'] and not isinstance(attr_value,bool):
            raise TypeError("{} must be either True or False".format(attr_name))

        self.__dict__[attr_name] = attr_value

    # Objective function
    def objectiveFunc_rule(self,optModel):
        # Reactions for which the objective coefficient has not bee assigned
        non_obj_rxns = [j for j in optModel.J if self.model.reactions_by_id[j].objective_coefficient == None]
        if len(non_obj_rxns) >= 1: 
            raise userError("ERROR! 'objective_coefficient' has not been defined for the following reactions: " + str(non_obj_rxns))
        if self.distance_type.lower() in ['euclidean','2-norm']:
            return sum(self.model.reactions_by_id[j].objective_coefficient*(optModel.v[j] - self.model.reactions_by_id[j].flux_wildType)*(optModel.v[j] - self.model.reactions_by_id[j].flux_wildType) for j in optModel.J)
        else:
            return sum(self.model.reactions_by_id[j].objective_coefficient*optModel.abs_vdiff[j] for j in optModel.J)
        
    # Absolute value constraints
    def absolute_const1_rule(self,optModel,j):
        return optModel.abs_vdiff[j] >= optModel.v[j] - self.model.reactions_by_id[j].flux_wildType 
    def absolute_const2_rule(self,optModel,j):
        return optModel.abs_vdiff[j] >= -(optModel.v[j] - self.model.reactions_by_id[j].flux_wildType) 
        
    def build_optModel(self):
        """
        This optModel creates a pyomo model for FBA optModel
        """
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()

        #--- Define sets ---
        # Set of compounds 
        optModel.I = Set(initialize = [c.id for c in self.model.compounds])

        # Set of rxns  
        optModel.J = Set(initialize = [r.id for r in self.model.reactions])

        #--- Define the optModel variables --- 
        optModel.v = Var(optModel.J, domain=Reals, bounds = lambda optModel, j: self.model.reactions_by_id[j].flux_bounds)

        # abs(v - v_WT)
        optModel.abs_vdiff = Var(optModel.J, domain=Reals, bounds = (0,1000))
        
        #--- Defiine the objective function and constraints ----
         # Objective function
        optModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = minimize)

        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule = lambda optModel, i: sum(j.stoichiometry[self.model.compounds_by_id[i]]*optModel.v[j.id] for j in self.model.compounds_by_id[i].reactions) == 0)

        # Absolute value constraints
        if self.distance_type.lower() == '1-norm':
            optModel.absolute_const1 = Constraint(optModel.J, rule=self.absolute_const1_rule)
            optModel.absolute_const2 = Constraint(optModel.J, rule=self.absolute_const2_rule)

        self.optModel = optModel 

#----------------------------------------------------    
if __name__ == "__main__":

    import time
    from tools.io.read_gams_model import read_gams_model
    from tools.io.read_sbml_model import read_sbml_model
    from set_specific_bounds import set_specific_bounds
    from cobra import test
 
    # Solver name
    optimization_solver = 'gurobi'

    #--- E. coli iAF1260 model ---
    start = time.clock()
    WT = read_gams_model(gams_model_file = '/fs/home06/alizom//models/Ecoli/iAF1260/iAF1260ModelData.py',model_name = 'iAF1260',organism_name = 'E. coli',model_type = 'metabolic')
    print '        Reading the gams model took ',str(time.clock() - start)

    WT.biomass_reaction = WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'})

    # Growth medium
    set_specific_bounds(WT,specific_bounds_file = '/data/alizom/models/Ecoli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign and objective function coefficients
    for rxn in WT.reactions:
        rxn.objective_coefficient = 0
    WT.biomass_reaction.objective_coefficient = 1

    print '   Perfomring FBA ...'
    WT.fba()

    print '   Perfomring MOMA ...'
    for rxn in WT.reactions:
        rxn.flux_wildType = rxn.flux
        rxn.objective_coefficient = 1
    momaModel = moma(model = WT, distance_type = '1-norm', optimization_solver = 'cplex')
    momaModel.run()
    for rxn in WT.reactions:
        print rxn.id,'    fba = ',rxn.flux_wildType,'  moma = ',rxn.flux,'\n'
