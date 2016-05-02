from __future__ import division
import sys
sys.path.append('../../')
import reaction 
from compartment import compartment
from tools.userError import *
from tools.custom_objects import customList, customDict
from copy import deepcopy

class compound(object):
    """
    A class holding the information about a compound 

    METHODS:
    -------
                   assign_props: Assign some properties to the compound related to other components of the model
                    reset_props: Resets the properties of the compound related to other components of the model
                  set_reactions: Sets the field reactions (direct modification of this attribute is not allowed)
          set_reactant_reactios: Sets attribute reactant_reactions (direct modification of this attribute is not allowed)
           set_product_reactios: Sets attribute product_reactions (direct modification of this attribute is not allowed)
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
    Last updated: 04-25-2016
    """
    def __init__(self, id, compartment = None, name = '', name_aliases = [], KEGG_id = [], ModelSEED_id = [], BiGG_id = [], formula = '', molecular_weight = None, charge = None, reactions = [], reactant_reactions = [], product_reactions = [], concentration = None, deltaG = None, deltaG_uncertainty = None, deltaG_range = [], notes = None, warnings = True): 

        # Warnings 
        self.warnings = warnings
    
        # Metabolite id or abbreviation in the model (string)
        self.id = id

        # Complete or expanded metabolite name (string)
        self.name = name

        # Name name_aliases
        self.name_aliases = name_aliases

        # Name of the compartment for this metabolite (string) 
        self.compartment = compartment

        # KEGG id
        self.KEGG_id = KEGG_id
        if isinstance(self.KEGG_id,str):
            self.KEGG_id = [self.KEGG_id]

        # ModelSEED id
        self.ModelSEED_id = ModelSEED_id 
        if isinstance(self.ModelSEED_id,str):
            self.ModelSEED_id = [self.ModelSEED_id]

        # BiGG id (BiGG is the database of metabolic models in Palsson's lab)
        self.BiGG_id = BiGG_id
        if isinstance(self.BiGG_id,str):
            self.BiGG_id = [self.BiGG_id]

        # Chemical formula of this metabolite (string) 
        self.formula = formula

        # Molecular weight
        self.molecular_weight = molecular_weight

        # Charge of the compound
        self.charge = charge

        # List of reaction objects in which this metabolite participates as a
        # reactant or product. 
        self.reactions = reactions

        # List of reaction objects in which this metabolite participates as a
        # reactant 
        self.reactant_reactions = reactant_reactions

        # List of reaction objects in which this metabolite participates as a
        # product 
        self.product_reactions = product_reactions

        # Concentration of this metabolite (real) 
        self.concentration = concentration

        # deltaG of reaction
        self.deltaG = deltaG

        # Uncertainty in the deltaG of reaction
        self.deltaG_uncertainty = deltaG_uncertainty

        # A tuple of the form (dGfmin,dGfmax) containing the min and max values of deltaG
        # (Gibbs free energy change) of formation for this reaction
        self.deltaG_range = deltaG_range 

        # Notes and comments (string)
        self.notes = notes

        # Assign other properties
        self.assign_props()

    def assign_props(self):
        """
        Assigns the properties of a compound related to the model they belong to
        """
        # List of reaction objects in which this metabolite participates as a
        # reactant or product. If the list is not provided, it can be constructed
        # from the list of reactant_reacitons and product_reactions (if any) 
        # NOTE: Construcing reactions from reactant_reactions or product_reacitons assume
        # that they contain the complete list of reactions this metabolite appears in as a
        # reactant or product
        if self.reactions == [] and (self.reactant_reactions != [] or self.product_reactions != []):
            self.reactions = list(set(self.reactant_reactions + self.product_reactions))

        # List of reaction objects in which this metabolite participates as a
        # reactant 
        if self.reactant_reactions == [] and len(self.reactions) >= 1:
            for rxn in [r for r in self.reactions if self in r.stoichiometry.keys() and r.stoichiometry[self] < 0]:
                self.reactant_reactions.append(rxn)

        # List of reaction objects in which this metabolite participates as a
        # product 
        if self.product_reactions == [] and len(self.reactions) >= 1:
            for rxn in [r for r in self.reactions if self in r.stoichiometry.keys() and r.stoichiometry[self] > 0]:
                self.product_reactions.append(rxn)

    def reset_props(self):
        """
        Resets the properties of a compound to default related to the model they belong to
        """
        self.model = None
        self.reactions = []
        self.reactant_reactions = []
        self.product_reactions = []
    
    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # Output messages and warnings 
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("warnings must be True or False")

        # id 
        if attr_name == 'id' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'id' for compound " + str(attr_value) + "! 'id' must be a string. A " + str(type(attr_value)) + ' was entered instead')

        # Name
        if attr_name == 'name' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'name' for compound " + self.id + "! 'name' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # Name aliases
        if attr_name == 'name_aliases' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'name_aliases' for compound " + str(self.id) + "! 'name_aliases' must be a list of strings. A " + str(type(attr_value)) + " type object was entered instead")
        if attr_name == 'name_aliases' and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError("Invalid 'name_aliases' for compound " + str(self.id) + "! 'name_aliases' must be a list of strings. Objects that are not string found in the list: " + str([n for n in attr_value if not isinstance(n,str)]))

        # Name of the compartment for this metabolite (string) 
        if attr_name == 'compartment' and (attr_value is not None and not isinstance(attr_value,compartment)):
            raise TypeError("Invalid 'compartment' for compound " + str(self.id) + "! 'compartment' must be an object of type compartment. A "  + str(type(attr_value)) + " type object was entered instead")

        # KEGG id
        if attr_name in ['ModelSEED_id', 'KEGG_id', 'BiGG_id'] and (not isinstance(attr_value,str) and not isinstance(attr_value,list)):
            raise TypeError('{} must be a list. A {} object were provided insteaded.'.format(attr_name, type(attr_value)))
        if attr_name in ['ModelSEED_id', 'KEGG_id', 'BiGG_id'] and isinstance(attr_value,list) and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError('Invalid {} for compound {}! {} must be a list of strings. Objects of type {} were observed in the list instead.'.format(attr_name, self.id, attr_name, list(set([type(n) for n in attr_value if not isinstance(n,str)])))) 

        # Formula 
        if attr_name == 'formula' and attr_value != None and not isinstance(attr_value,str):
            raise TypeError("Invalid 'formula' for compound " + str(self.id) + "! 'formula'  must be a string. A " + str(type(attr_value)) + " type object was entered instead")
        elif attr_name == 'formula' and attr_value == None: 
            self.__dict__[attr_name] = ''

        # Molecular weight
        if attr_name == 'molecular_weight' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invalid 'molecular_weight' for compound " + str(self.id) + "! 'molecular_weight'  must be either a float or an integer. A " + str(type(attr_value)) + " type object was entered instead")

        # Charge of the compound
        if attr_name == 'charge' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invalid 'id' for compound " + str(self.id) + "Invalid 'charge' for compound " + str(self.id) + "! 'charge'  must be either a float or an integer. A " + str(type(attr_value)) + " type object was entered instead")

        # deltaG
        if attr_name == 'deltaG' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invlaud 'deltaG' for compound " + self.id + "! 'deltaG'  must be either a float or an integer. A " + str(type(attr_value)) + " type object was entered instead")

        # deltaG_uncertainty
        if attr_name == 'deltaG_uncertainty' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invlaud 'deltaG_uncertainty' for compound " + self.id + "! deltaG_uncertainty  must be either a float or an integer. A " + str(type(attr_value)) + " type object was entered instead")

        # deltaG_range 
        if attr_name == 'deltaG_range' and not isinstance(attr_value,list):
            raise TypeError("Invlaud 'deltaG_range' for compound " + self.id + "! 'deltaG_range' must be a list. A " + str(type(attr_value)) + " type object was entered instead")


        # Concentration 
        if attr_name == 'concentration' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float) and not isinstance(attr_value,dict)):
            raise TypeError("Invalid 'concentration' for compound " + str(id) + "! 'concentration' must be either a float or an integer or a dictionary. A " + str(type(attr_value)) + " type object was entered instead")

        if attr_name in 'reactions':
            self.set_reactions(reactions = attr_value)
        elif attr_name in 'reactant_reactions':
            self.set_reactant_reactions(reactant_reactions = attr_value)
        elif attr_name in 'product_reactions':
            self.set_product_reactions(product_reactions = attr_value)
        else:
            self.__dict__[attr_name] = attr_value

    def set_reactions(self, reactions):
        """
        Makes changes to attribute reactions
        """
        if reactions is not None and not isinstance(reactions,list):
            raise TypeError("Invalid 'reactions' for compound " + str(self.id) + "! 'reactions'  must be a list of objects of type reaction. A " + str(type(reactions)) + " type object was entered instead")
        if len([n for n in reactions if not isinstance(n,reaction.reaction)]) > 0:
            raise TypeError("Invalid 'reactions' for compound " + str(self.id) + "! 'reactions'  must be a list of objects of type 'reaction'.  Objects that are not of type reaction found in the list:" + str([n for n in reactions if not isinstance(n,reaction.reaction)]))

        problem_rxns = [r.id for r in reactions if self not in r.stoichiometry.keys()]

        if len(problem_rxns) > 0: 
            raise userError('The following {} reactions appears in "reactions" of compound {} but do not appear in the stoichioemtry of this reaction: {}'.format(len(problem_rxns), self,id, problem_rxns))

        reactions = sorted(reactions,key=lambda x:x.id)

        self.__dict__['reactions'] = customList(reactions)

    def set_reactant_reactions(self, reactant_reactions):
        """
        Makes changes to attribute reactant_reactions
        """
        if (reactant_reactions is not None and not isinstance(reactant_reactions,list)):
            raise TypeError("Invalid 'reactant_reactions' for compound " + str(self.id) + "! 'reactant_reactions'  must be a list of objects of type reaction. A " + str(type(reactant_reactions)) + " type object was entered instead")
        if len([n for n in reactant_reactions if not isinstance(n,reaction.reaction)]) > 0:
            raise TypeError("Invalid 'reactant_reactions' for compound " + str(self.id) + "! 'reactant_reactions'  must be a list of objects of type reaction. Objects that are not of type reaction found in the list: " + str([n for n in reactant_reactions if not isinstance(n,reaction.reaction)]))

        problem_rxns = [r.id for r in reactant_reactions if self not in [c for c in r.stoichiometry.keys() if r.stoichiometry[c] < 0]]
        if len(problem_rxns) > 0: 
            raise userError('The following {} reactions appears in "reactant_reactions" of compound {} but it does not appear in the stoichioemtry of this reaction as a reactant {}'.format(len(problem_rxns), self,id, problem_rxns))

        reactant_reactions = sorted(reactant_reactions,key=lambda x:x.id)

        self.__dict__['reactant_reactions'] = customList(reactant_reactions)

    def set_product_reactions(self, product_reactions):
        """
        Makes changes to attribute product_reactions
        """
        if product_reactions is not None and not isinstance(product_reactions,list):
            raise TypeError("Invalid 'product_reactions' for compound " + str(self.id) + "! 'reactant_reactions'  must be a list of objects of type reaction. A " + str(type(product_reactions)) + " type object was entered instead")
        if len([n for n in product_reactions if not isinstance(n,reaction.reaction)]) > 0:
            raise TypeError("Invalid 'product_reactions' for compound " + str(self.id) + "! 'product_reactions'  must be a list of objects of type reaction. Objects that are not of type reaction found in the list: " + str([n for n in product_reactions if not isinstance(n,reaction.reaction)]))

        problem_rxns = [r.id for r in product_reactions if self not in [c for c in r.stoichiometry.keys() if r.stoichiometry[c] > 0]]
        if len(problem_rxns) > 0: 
            raise userError('The following {} reactions appears in "product_reactions" of compound {} but it does not appear in the stoichioemtry of this reaction as a product {}'.format(len(problem_rxns), self,id, problem_rxns))

        product_reactions = sorted(product_reactions,key=lambda x:x.id)

        self.__dict__['product_reactions'] = customList(product_reactions)

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
            raise userError("Invalid reference type (eligible choices are 'id' and 'name')")
        elif self.reactions == None:
            raise userError('List of the reactions is not defined for this metabolite ...')

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
            raise userError("Invalid reference type (eligible choices are 'id' and 'name')")
        elif self.reactions == None:
            raise userError('List of the reactions is not defined for this metabolite ...')

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
            raise userError("Invalid reference type (eligible choices are 'id' and 'name')")
        elif self.reactions == None:
            raise userError('List of the reactions is not defined for this metabolite ...')

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


    def biomass_yield_calc(self, warnings = None, stdout_msgs = False):
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
        from coopr import pyomo
        from tools.fba.fba import fba
        import copy

        if warnings == None:
            warnings = self.warnings
        elif not isinstance(warnings,bool):
            raise TypeError('warnings must be either True or False')

        # Find the exchange reaction this metabolite participates in
        exch_rxn = [r for r in self.reactant_reactions if r.reversibility.lower() == 'exchange']
        if len(exch_rxn) == 0:
            raise userError('Unable to compute the biomass yield for metabolite ' + self.id +'. The metabolite must have an excange reaction in reactant_reactions ...')
        elif len(exch_rxn) > 1:
             raise userError('metabolite ' + self.id + ' participates in more than one exchange reaciton ...')
        else:
            exch_rxn = exch_rxn[0]

        # Find reactions participating in the objective function
        obj_rxns = [r for r in self.model.reactions if r.objective_coefficient != 0]

        # Save the objective coefficients for these reactions in a dictionary
        obj_coeff = dict([(r,r.objective_coefficient) for r in obj_rxns])

        # Set the objective coefficient for all reactions to zero
        for r in self.model.reactions:
            r.objective_coefficient = 0
        if self.model.biomass_reaction == None:
            raise AttributeError("'biomass_reaction' is not defined for the model")
        else:
            self.model.biomass_reaction.objective_coefficient = 1

        # Assign a large LB (allow for unlimitted uptake) 
        # Original LB on reaction flux 
        exch_rxn_LB = exch_rxn.flux_bounds[0]
  
        # Perform FBA for 10 mole uptake of this metabolite 
        exch_rxn.flux_bounds[0] = -10      

        if hasattr(self.model,'fba_model'):
            biomass_yield_fba_model = copy.deepcopy(self.model.fba_model)
            biomass_yield_fba_model.optModel.del_component('objectiveFunc')
            biomass_yield_fba_model.optModel.objectiveFunc = pyomo.Objective(rule=self.model.fba_model.objectiveFunc_rule, sense = pyomo.maximize)
            biomass_yield_fba_model.build_new_optModel = False
            biomass_yield_fba_model.store_opt_fluxes = False
            biomass_yield_fba_model.stdout_msgs = False
            biomass_yield_fba_model.optModel.v[exch_rxn.id].setlb(exch_rxn.flux_bounds[0])
        else: 
            biomass_yield_fba_model = fba(build_new_optModel = True, reset_fluxes = False, store_opt_fluxes = False, stdout_msgs = False)

        # Solve the fba model
        fba_solution = biomass_yield_fba_model.run() 
 
        if fba_solution['exit_flag'] == 'globallyOptimal':
            if fba_solution['opt_rxnFluxes'][exch_rxn.id] < 0:
                self.biomass_yield = fba_solution['objective_value']/(-fba_solution['opt_rxnFluxes'][exch_rxn.id])
            else:
                self.biomass_yield = None 
                if  fba_solution['opt_rxnFluxes'][exch_rxn.id] == 0 and self.warnings:
                    print 'WARNING (compound.py)! The fba problem to compute the biomass yield for compound ' + self.id + ' resulted in a zero exchange flux value for the corresponding exchange reaction. No value is assigned to biomass_yield'
                     
                if  fba_solution['opt_rxnFluxes'][exch_rxn.id] > 0 and self.warnings:
                    print 'WARNING (compound.py)! The fba problem to compute the biomass yield for compound ' + self.id + ' resulted in a non-negative flux value for the corresponding exchange reaction. No value is assigned to biomass_yield'
        else:
            self.biomass_yield = None
            if  warnings:
                print 'WARNING (compound.py! The fba problem to compute the biomass yield for compound ' + self.id + ' was not solved to optimality. No value is assigned to biomass_yield'

        # Set the objective coefficients back to what they were before
        exch_rxn.objective_coefficient = 0
        for r in obj_rxns:
            r.objective_coefficient = obj_coeff[r]
 
        # Set the lower bound on exchange reaction back to what it was before
        exch_rxn.flux_bounds[0] = exch_rxn_LB

    def ms_calc(self, warnings = None, stdout_msgs = True):
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
        from coopr import pyomo
        from tools.fba.fba import fba
        import copy

        if warnings == None:
            warnings = self.warnings
        elif not isinstance(warnings,bool):
            raise TypeError('warnings must be either True or False')

        # Find the exchange reaction this metabolite participates in
        exch_rxn = [r for r in self.reactant_reactions if r.reversibility.lower() == 'exchange']
        if len(exch_rxn) == 0:
            raise userError('Unable to compute the biomass yield for metabolite ' + self.id +'. The metabolite must have an excange reaction in reactant_reactions ...')
        elif len(exch_rxn) > 1:
             raise userError('metabolite ' + self.id + ' participates in more than one exchange reaciton ...')
        else:
            exch_rxn = exch_rxn[0]

        # Find reactions participating in the objective function
        obj_rxns = [r for r in self.model.reactions if r.objective_coefficient != 0]

        # Save the objective coefficients for these reactions in a dictionary
        obj_coeff = dict([(r,r.objective_coefficient) for r in obj_rxns])

        # Set the objective coefficient for all reactions to zero
        for r in self.model.reactions:
            r.objective_coefficient = 0
                    
        # Set the coefficient for the exchange reaction of this compound to one
        exch_rxn.objective_coefficient = 1

        # Assign a large LB (allow for unlimitted uptake) 
        exch_rxn_LB = exch_rxn.flux_bounds[0]
        exch_rxn.flux_bounds[0] = -1000 

        if not hasattr(self.model,'fba_model'):
            ms_fba_model = copy.deepcopy(self.model.fba_model)
            ms_fba_model.optModel.del_component('objectiveFunc')
            ms_fba_model.optModel.objectiveFunc = pyomo.Objective(rule=self.model.fba_model.objectiveFunc_rule, sense = pyomo.maximize)
            ms_fba_model.build_new_optModel = False
            ms_fba_model.store_opt_fluxes = False
            ms_fba_model.stdout_msgs = False
            ms_fba_model.optModel.v[rxch_rxn.id].setlb(exch_rxn.flux_bounds[0])
        else: 
            ms_fba_model = fba(model = self.model, build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False)

        # Solve the fba model
        fba_solution = ms_fba_model.run()

        if fba_solution['exit_flag'] == 'globallyOptimal':
            self.ms = -min([0,fba_solution['objective_value']])
        else:
            self.ms = None
            if  warnings:
                print 'WARNING (compound.py)! The fba problem to compute ms for compound ' + self.id + ' was not solved to optimality. No value is assigned to biomass_yield'
        
        # Set the objective coefficients back to what they were before
        exch_rxn.objective_coefficient = 0
        for r in obj_rxns:
            r.objective_coefficient = obj_coeff[r]
 
        # Set the lower bound on exchange reaction to what it was before
        exch_rxn.flux_bounds[0] = exch_rxn_LB 

