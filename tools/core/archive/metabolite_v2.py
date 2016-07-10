from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from copy import deepcopy
from coopr.pyomo import *
from coopr.opt import *

class metabolite(object):
    """
    A class holding the information about a biochemical 

    METHODS:
    -------
             print_reactions_by: Prints in the out the list of reactions in which this 
                                 imetabolite partcipates as a reactnat or product 
    print_reactant_reactions_by: Prints in the output the list of reacitons in which this
                                 metabolite participates as a reactant
     print_product_reactions_by: Prints in the output the list of reactions in which this 
                                 metabolite participates as a product
             biomass_yield_calc: Calculates the biomass yield of this metabolite (if applicable)
                        ms_calc: Calculates the minimum required uptake of this metabolite to 
                                 satisfy the non-growth associated ATPM maintenance (if applicable)
 
    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 03-07-2015
    """

    def __init__(self, id, compartment = None, name = '', aliases = '', Kegg_id = '', ModelSeed_id = '', formula = '', charge = None, reactions = [], reactant_reactions = [], product_reactions = [], concentration = None, deltaGf_range = [], biomass_yield = None, ms = None, notes = ''): 
    
        # Metabolite id or abbreviation in the model (string)
        self.id = id

        # Complete or expanded metabolite name (string)
        self.name = name

        # Name aliases
        self.aliases = aliases

        # Name of the compartment for this metabolite (string) 
        self.compartment = compartment

        # Kegg id
        self.Kegg_id = Kegg_id

        # ModelSeed id
        self.ModelSeed_id = ModelSeed_id 

        # Chemical formula of this metabolite (string) 
        self.formula = formula

        # Charge of the compound
        self.charge = charge

        # List of reaction objects in which this metabolite participates as a
        # reactant or product. If the list is not provided, it can be constructed
        # from the list of reactant_reacitons and product_reactions (if any) 
        # NOTE: Construcing reactions from reactant_reactions or product_reacitons assume
        # that they contain the complete list of reactions this metabolite appears in as a
        # reactant or product
        self.reactions = reactions
        if self.reactions == [] and (len(reactant_reactions) >= 1 or len(product_reactions) >= 1):
            self.reactions = list(set(reactant_reactions + product_reactions))

        # List of reaction objects in which this metabolite participates as a
        # reactant 
        self.reactant_reactions = reactant_reactions
        # The following statement is int he try block because the self.reaciotns
        # may not be necessarily a list of reaciton objects. This is because sometimes
        # the user may just provide strings for the filed reacitons, instead of
        # reaciton objects
        try:
            if self.reactant_reactions == [] and len(self.reactions) >= 1 and isinstance(self.reactions[0],reaction):
                for rxn in [r for r in self.reactions if self in r.stoichiometry.keys()]:
                    if r.stoichiometry[self] < 0:
                        self.reactant_reactions.append(rxn)
        except:
            pass

        # List of reaction objects in which this metabolite participates as a
        # product 
        self.product_reactions = product_reactions
        try:
            if product_reactions == [] and len(self.reactions) >= 1:
                for rxn in [r for r in self.reactions if self in r.stoichiometry.keys()]:
                    if r.stoichiometry[self] > 0:
                        self.product_reactions.append(rxn)
        except:
            pass

        # Concentration of this metabolite (real) 
        self.concentration = concentration

        # A tuple of the form (dGfmin,dGfmax) containing the min and max values of deltaG
        # (Gibbs free energy change) of formation for this reaction
        self.deltaGf_range = deltaGf_range 

        # Biomass yield 
        self.biomass_yield = biomass_yield

        # ms (the minimum required uptake of this metabolite to satisfy ATP maintenance
        # requirements)
        self.ms = ms

        # Notes and comments (string)
        if isinstance(notes,str):
            self.notes = notes
        else:
            self.notes = ''

    def print_reactions_by(self,ref_type, metab_ref = None):
        """
        Prints in the output the list of all reactions this metabolite participates in 
        as a reactant or product in the model with the format specified by ref_type 
 
        INPUTS:
        -------
         ref_type: A string indicating  in what format the reactions should be 
                   printed in the output. Current eligible choices are 'id', 'name'
                   and formula. If name or id was not provided, id used instead.
        metab_ref: A string indicating in what format the metabolites should appear
                   in a reaciton equation if ref_type = 'equation'. If None, metabolite
                   ids are used
        """
        if ref_type.lower() not in ['id','name','equation']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")
        elif self.reactions == []:
            raise userError('**Error! List of the reactions is not defined for this metabolite ...')

        for reaction in self.reactions:
            if ref_type.lower() == 'id':
                rxn = reaction.id
            elif ref_type.lower() == 'name':
                if reaction.name is None:
                    rxn = reaction.id
                    userWARNING("Waarning! No name was not provided for '" + reaction.id + "'. id is used instead.")
                else:
                    rxn = reaction.name
            elif ref_type.lower() == 'equation':
                if metab_ref != None:
                    rxn = reaction.get_eqn_by(metab_ref)
                else:
                    rxn = reaction.get_eqn_by('id')

            print rxn 

    def print_reactant_reactions_by(self,ref_type):
        """
        Prints in the output the list of all reactions this metabolite participates in 
        as a reactant in the model with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the metabolites should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula. If name or id was not provided, id used instead.
        """
        if ref_type.lower() not in ['id','name','equation']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")
        elif self.reactions == []:
            raise userError('**Error! List of the reactions is not defined for this metabolite ...')

        for reaction in self.reactant_reactions:
            if ref_type.lower() == 'id':
                rxn = reaction.id
            elif ref_type.lower() == 'name':
                if reaction.name is None:
                    rxn = reaction.id
                    userWARNING("Waarning! No name was provided for '" + reaction.id + "'. id is used instead.")
                else:
                    rxn = reaction.name
            elif ref_type.lower() == 'equation':
                if reaction.equation is None:
                    rxn = reaction.id
                    userWARNING("Waarning! No equation was provided for '" + reaction.id + "'. id is used instead.")
                else:
                    rxn = reaction.equation

            print rxn 

    def print_product_reactions_by(self,ref_type):
        """
        Prints in the output the list of all reactions this metabolite participates in 
        as a product in the model with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the metabolites should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula. If name or id was not provided, id used instead.
        """
        if ref_type.lower() not in ['id','name','equation']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")
        elif self.reactions == []:
            raise userError('**Error! List of the reactions is not defined for this metabolite ...')

        for reaction in self.product_reactions:
            if ref_type.lower() == 'id':
                rxn = reaction.id
            elif ref_type.lower() == 'name':
                if reaction.name is None:
                    rxn = reaction.id
                    userWARNING("Waarning! No name was provided for '" + reaction.id + "'. id is used instead.")
                else:
                    rxn = reaction.name
            elif ref_type.lower() == 'equation':
                if reaction.equation is None:
                    rxn = reaction.id
                    userWARNING("Waarning! No equation was provided for '" + reaction.id + "'. id is used instead.")
                else:
                    rxn = reaction.equation

            print rxn


    def biomass_yield_calc(self,fba_model = None):
        """
        Calculates the biomass yield of this metabolite (if applicable).
        For the yield to be computed this metabolite should particiapte in
        in an exchange reaction (EX_m(e): m[e] <==>)

        INPUTS:
        -------
        fba_model: An instance of the class fba containing the fba model for the 
                   wild-type organism under the desired growth condition. If this is 
                   not provided, then the current fba_model of model is used 

        OUTPUS:
        -------
        The output is assigned to the global variable self.biomass_yield 
        """
        # Find the exchange reaction this metabolite participates in
        exch_rxn = [r for r in self.reactant_reactions if r.type.lower() == 'exchange']
        if len(exch_rxn) == 0:
            raise userError('**ERROR! Unable to compute the biomass yield for metabolite ' + self.id +'. The metabolite must have an excange reaction in reactant_reactions ...')
        elif len(exch_rxn) > 1:
             raise userError('**ERROR! metabolite ' + self.id + ' participates in more than one exchange reaciton ...')
        else:
            exch_rxn = exch_rxn[0]

            # Original LB on reaction flux 
            exch_rxn_LB = exch_rxn.flux_bounds[0]
  
            # Perform FBA for 10 mole uptake of this metabolite 
            exch_rxn.flux_bounds[0] = -10      

            if fba_model == None:
                # WARNING! If no fba_model is provided, the current model.fba_model is 
                # used to compute the biomass yield
                self.model.fba(create_model = False,store_opt_fluxes = True, screen_output = 'off')
                fba_solution = self.model.fba_model.solution
            else:
                fba_model.store_opt_fluxes = True
                fba_model.screen_output = 'off' 
                fba_model.run()
                fba_solution = fba_model.solution
 
            if fba_solution['exit_flag'] == 'globallyOptimal':
                if exch_rxn.flux < 0:
                    self.biomass_yield = fba_solution['objective_value']/(-exch_rxn.flux)
                else:
                    # if zero (no uptake) or positive (production)
                    self.biomass_yield = None 
                    print '\n**WARNING! The fba problem to compute the biomass yield for metabolite ' + self.id + ' resulted in a non-negative flux value for the corresponding exchange reaction. No value is assigned to biomass_yield'
            else:
                self.biomass_yield = None
                print '\n**WARNING! The fba problem to compute the biomass yield for metabolite ' + self.id + ' was not solved to optimality. No value is assigned to biomass_yield'

            # Set the lower bound on exchange reaction back to what it was before
            exch_rxn.flux_bounds[0] = exch_rxn_LB

    def ms_calc(self,fba_model = None):
        """
        Calculates the minimum required uptake rate of this metabolite to satisfy the 
        non-growth associated ATP maintenance. This can be computed by perfomring
        FBA simulations where the objective function is to minimize the flux of the 
        uptake reaction. Note that since we work with exchagne reactions (Where uptake
        is modeled by having a negative upper bound), minimizing the flux of the uptake
        reaction is equivalent to maximizing the flux of exchange reaction. If the 
        objective function value turns out to be negative, it represents the value for
        -ms. If it is zero or positive, it means that this metabolite is not a limiting 
        nutrient to maintain the required ATPM flux (ms = 0). 

        INPUTS:
        -------
        fba_model: An instance of the class fba containing the fba model for the 
                   ms simulations. If not fba_model is provided, then model.fba_model  
                   is used to construct the required model. 

        OUTPUS:
        -------
        The value of ms, which is assigned to the global variable self.ms 
        """
        # Find the exchange reaction this metabolite participates in
        exch_rxn = [r for r in self.reactant_reactions if r.type.lower() == 'exchange']
        if len(exch_rxn) == 0:
            raise userError('**ERROR! Unable to compute the biomass yield for metabolite ' + self.id +'. The metabolite must have an excange reaction in reactant_reactions ...')
        elif len(exch_rxn) > 1:
             raise userError('**ERROR! metabolite ' + self.id + ' participates in more than one exchange reaciton ...')
        else:
            exch_rxn = exch_rxn[0]

            if fba_model == None:
                # WARNING! If no fba_model is provided, the current model.fba_model 
                # is used to compute the biomass yield

                # Find reactions participating in the objective function
                obj_rxns = [r for r in self.model.reactions if r.objective_coefficient != 0]

                # Save the objective coefficients for these reactions in a dictionary
                obj_coeff = dict([(r,r.objective_coefficient) for r in obj_rxns])

                # Set the objective coefficient for all reactions to zero
                for r in self.model.reactions:
                    r.objective_coefficient = 0
                    
                # Set the coefficient for the exchange reaction of this metabolite to one
                exch_rxn.objective_coefficient = 1

                # current LB on the exchange reaction
                exch_rxn_LB = exch_rxn.flux_bounds[0]
                # Assign a large LB (allow for unlimitted uptake) 
                exch_rxn.flux_bounds[0] = -1000 

                self.model.fba(create_model = False,store_opt_fluxes = False, screen_output = 'off')
                fba_solution = self.model.fba_model.solution
            else:
                fba_model.store_opt_fluxes = False
                fba_model.screen_output = 'off'
                fba_model.run()
                fba_solution = fba_model.solution
 
            if fba_solution['exit_flag'] == 'globallyOptimal':
                self.ms = -min([0,fba_solution['objective_value']])
            else:
                self.ms = None
                print 'WARNING! The fba problem to compute ms for metabolite ' + self.id + ' was not solved to optimality. No value is assigned to biomass_yield'
        
            # Set the objective coefficients back to what they were before
            exch_rxn.objective_coefficient = 0
            for r in obj_rxns:
                r.objective_coefficient = obj_coeff[r]
 
            # Set the lower bound on exchange reaction to what it was before
            exch_rxn.flux_bounds[0] = exch_rxn_LB 

            # Re-define the objective function
            self.model.fba_model.pyomo_fbaModel.del_component('objectiveFunc')
            self.model.fba_model.pyomo_fbaModel.objectiveFunc = Objective(rule=self.model.fba_model.objectiveFunc_rule, sense = maximize)
