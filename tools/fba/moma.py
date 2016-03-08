from __future__ import division
import sys, time
sys.path.append('../../')
from tools.fba.fbaTools import fbaTools

class moma(fbaTools):
    """
    Performs flux balance analysis

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 10-27-2015 
    """   

    # Class attributes that are set by __init__
    class_attr__init__ = ['penalties','distance_type','optimization_solver','create_optModel','flux_key','store_opt_fluxes','stdout_msgs','warnings']

    def __init__(self,model, distance_type = '2-norm',optimization_solver = 'gurobi', create_optModel = True, flux_key = None, store_opt_fluxes = True, warnings = True, stdout_msgs = True): 

        """
        INPUTS:
        ------
              distance_type: The way to compute the distance between the fluxes.
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
        self.create_optModel = create_optModel
               
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
            # Check if the all reactions in the model have a field called 
            # flux_wildtype containing the wild-type reaction fluxes
            no_wt_flux_rxns = [r for r in attr_value.reactions if 'wildtype_flux' not in dir(r)]
            if len(no_wt_flux_rxns) > 0 and len(no_wt_flux_rxns) <= 10:
                raise AttributeError('No wildtype_flux attribute for ' + str(len(no_wt_flux_rxns)) + ' reactions: ' + str(no_wt_flux_rxns))
            elif len(no_wt_flux_rxns) > 0 and len(no_wt_flux_rxns) > 10:
                raise AttributeError('No wildtype_flux attribute for ' + str(len(no_wt_flux_rxns)) + ' reactions including: ' + str(no_wt_flux_rxns[:10]) + '\nand ' + str(len(no_wt_flux_rxns) - 10) + ' more reactions.')

        # distance_type 
        if attr_name == 'distance_type' and attr_value.lower() not in ['euclidean','2-norm','1-norm']:
            raise ValueError('Invalid value for distance_type. Eligible inputs are euclidean or 2-norm and 2-norm')

        # Solver name
        if attr_name == 'optimization_solver' and attr_value.lower() not in ['cplex','gurobi']:
            raise userError('Invalid solver name (eligible choices are cplex and gurobi)\n')

        # Warnings and messages in the standard output
        if attr_name == 'create_optModel' and not isinstance(attr_value,bool):
            raise TypeError("create_optModel must be either True or False")

        # Warnings and messages in the standard output
        if attr_name == 'stdout_msgs' and not isinstance(attr_value,bool):
            raise TypeError("stdout_msgs must be either True or False")

        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("warnings must be either True or False")

    # Objective function
    def objectiveFunc_rule(self,optModel):
        # Reactions for which the objective coefficient has not bee assigned
        non_obj_rxns = [j.id for j in optModel.J if j.objective_coefficient == None]
        if len(non_obj_rxns) >= 1: 
            raise userError("ERROR! 'objective_coefficient' has not been defined for the following reactions: " + str(non_obj_rxns))
        if self.distance_type.lower() in ['euclidean','2-norm']:
            return sum(j.objective_coefficient*(optModel.v[j] - j.wildtype_flux)*(optModel.v[j] - j.wildtype_flux) for j in optModel.J)
        else:
            return sum(j.objective_coefficient*optModel.abs_vdiff[j] for j in optModel.J)
        
    # Absolute value constraints
    def absolute_const1_rule(self,optModel,j):
        return optModel.abs_vdiff[j] >= optModel.v[j] - j.wildtype_flux 
    def absolute_const2_rule(self,optModel,j):
        return optModel.abs_vdiff[j] >= -(optModel.v[j] - j.wildtype_flux) 
        
    def createPyomoModel(self):
        """
        This creates a pyomo model for moma 
        """   
        #--- Create a pyomo model optModel ---
        optModel = ConcreteModel()
        
        #--- Define sets ---
        # Set of compounds 
        optModel.I = Set(initialize = self.model.compounds)   

        # Set of rxns  
        optModel.J = Set(initialize = self.model.reactions)     

        #--- Define the optModel variables --- 
        # Reaction fluxes
        optModel.v = Var(optModel.J, domain=Reals, bounds = self.assignFluxBounds)

        # abs(v - v_WT)
        optModel.abs_vdiff = Var(optModel.J, domain=Reals, bounds = (0,1000))
        
        #--- Defiine the objective function and constraints ----
         # Objective function
        optModel.objectiveFunc = Objective(rule=self.objectiveFunc_rule, sense = minimize)

        # Mass balance 
        optModel.massBalance_const = Constraint(optModel.I, rule=self.massBalance_rule)

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
        rxn.wildtype_flux = rxn.flux
        rxn.objective_coefficient = 1
    momaModel = moma(model = WT, distance_type = '1-norm', optimization_solver = 'cplex')
    momaModel.run()
    for rxn in WT.reactions:
        print rxn.id,'    fba = ',rxn.wildtype_flux,'  moma = ',rxn.flux,'\n'
