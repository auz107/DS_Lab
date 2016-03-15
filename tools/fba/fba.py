from __future__ import division
import sys, time
sys.path.append('../../')
from tools.fba.fbaTools import fbaTools
from tools.userError import userError

class fba(fbaTools):
    """
    Performs flux balance analysis

    Ali R. Zomorrodi - Segre Lab @ Boston University
    Last updated: 10-27-2015 
    """   
    def objectiveFunc_rule(self,fba_optModel):
        """
        Objective function
        """
        # Reactions for which the objective coefficient has not bee assigned
        non_obj_rxns = [j.id for j in self.model.reactions if j.objective_coefficient == None]
        if len(non_obj_rxns) > 0: 
            if len(non_obj_rxns) <= 10:
                raise userError("'objective_coefficient' has not been defined for the following reactions: {}".format(non_obj_rxns[:10]))
            else:
                raise userError("'objective_coefficient' has not been defined for the following reactions: {} and {} more.".format(non_obj_rxns[:10],len(non_obj_rxns) - 10))
        return sum(j.objective_coefficient*fba_optModel.v[j.id] for j in self.model.reactions)

