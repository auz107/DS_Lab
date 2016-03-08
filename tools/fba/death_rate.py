from __future__ import division
import sys, time
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.userError import *
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model

def death_rate(limiting_nutrients, ms_biomassYield_calc = True, wildType_max_biomass = None, warnings = True, stdout_msgs = True):
    """
    Compute the death rate. If ms and biomass yield is not available for none of any of the limiting
    nutrients or if they could not be computed successfully, then the death rate is set to -1% of
    max biomass flux for the wild-type strain. If the max biomass flux for the wild-type is not 
    provided then the death rate is set to -0.001 (1/h). Note that the value of death rate is 
    non-positive

    INPUTS:
    -------
      limiting_nutrients: A list of objects of type compound containing 
                          the list of limiting nutrients
    ms_biomassYield_calc: If True ms and biomass yield is calculated for nutrients with no previously
                          assigned ms and biomass_yield. 
    wildType_max_biomass: In case ms and biomass_yield for limiting nutrients are not available
                          or cannot be computed successfully, the death rate is set to 1% of
                          max biomass flux for the wild-type. If wildType_max_biomass is not
                          provided either, then the death rate is set to 0.01 1/h
                warnings: Can be 'on' or 'off' shwoing whether the warnings should be writtten to the 
                          screen or not
             stdout_msgs: Can be 'on' or 'off' shwoing whether any messages should be written to the
                          screen

    OUTPUT:
    -------
              mu_death: death_rate (non_positive)

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 08-21-2015 
    """   
    # Warnings and messages in the standard output
    if not isinstance(stdout_msgs,bool):
        raise TypeError("stdout_msgs must be either True or False")
    if not isinstance(warnings,bool):
        raise TypeError("warnings must be either True or False")

    # Compute the death phase mu for each limiting nutrients
    mu_death = []

    for nutrient in limiting_nutrients:
        if not hasattr(nutrient,'biomass_yield') or (hasattr(nutrient,'biomass_yield') and nutrient.biomass_yield == None) and ms_biomassYield_calc:
            # Compute Yx/s (biomass yield)
            nutrient.biomass_yield_calc(warnings = warnings, stdout_msgs = stdout_msgs)  

        if not hasattr(nutrient,'ms') or (hasattr(nutrient,'ms') and nutrient.ms == None) and ms_biomassYield_calc:
            # Compute Yx/s (biomass yield)
            nutrient.ms_calc(warnings = warnings, stdout_msgs = stdout_msgs)  

        # If a value was computed for ms and biomass_yield
        if nutrient.biomass_yield != None and nutrient.ms != None:
 
            # Exchange (uptake) reaction for this compound
            exch_rxn = [r for r in nutrient.reactant_reactions if r.type.lower() == 'exchange']
            if len(exch_rxn) > 1:
                raise userError('compound ' + nutrient.id + ' has more than one exchange reaction')
            else:
                exch_rxn = exch_rxn[0]

            # Compute the lower bound on exchange reaction (upper bound on uptake)
            exch_rxn.kinetic_rate_calc()

            # mu_death = Yxs*(v_uptake - ms)
            if abs(exch_rxn.kinetic_rate) > abs(nutrient.ms):
                mu_death.append(0)
            else:
                mu_death.append((abs(exch_rxn.kinetic_rate) - abs(nutrient.ms))*nutrient.biomass_yield)

    # If more than one limitting nutrients have run out, set mu_death to the largest value
    # (the one with most negative value)
    if len(mu_death) > 0:
        mu_death = min(mu_death)
    elif wildType_max_biomass != None:
        mu_death = - 0.01*wildType_max_biomass
        if warnings:
            print 'WARNING! mu_death could not be computed successfully using ms and biomass yield for limiting nutrients. 1% of max biomass flux for the wild-type under the same condition was assigned as the death rate'

    else:
        mu_death = -0.001
        if warnings:
            print 'WARNING! mu_death could not be computed successfully using ms and biomass yield for limiting nutrients. Since max biomass flux for the wild-type was not also provided, a general value of -0.001 was assigned as the death rate'

    if mu_death > 0:
        raise userError('Positive death rate (death rate = ' + str(mu_death))

    if stdout_msgs:
        print 'mu_death = ',mu_death
                                 
    return mu_death

