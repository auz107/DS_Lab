from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from compound import compound
from gene import gene
from compartment import compartment

class reaction(object):
    """
    A class holding the information about a biochemical reaction 

    Methods:
         get_equation: Get the reaction equation
        get_compounds: Prints in the output the list of all compounds 
                       participating in this reaction as a reactant or product
        get_reactants: Prints in the output the list of all compounds 
                       participating in this reaction as a reactan 
         get_products: Prints in the output the list of all compounds 
                       participating in this reaction as a product
        add_compounds: Adds new compounds to the reaction
     remove_compounds: Removes compounds from the reaction
    kinetic_rate_calc: Computes the kinetic rate of the reaction

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 08-17-2015
    """
    def __init__(self, id, stoichiometry = {}, type = '', name = None, synonyms = [], Kegg_id = None, ModelSeed_id = None, EC_number = None, subsystem = None, compartment = [], genes = None, gene_reaction_rule = None, objective_coefficient = None, flux = None, store_flux = True, flux_bounds = None, deltaG_range = None, kinetics = None, kinetic_compounds = None, confidence_level =None, notes = None, warnings = True): 

        # Warnings and messages in the standard output
        if not isinstance(warnings,bool): 
            raise TypeError("warnings must be True or False")
        else:
            self.warnings = warnings
    
        # Reaction id or abbreviation in the model
        self.id = id

        # Reaction stoichiometry. A dictionary with keys and values as follows
        #   Keys: Instances of the object compound 
        # Values: Their stoichiometric coefficient in this reaction
        self.stoichiometry = stoichiometry
        
        # A string indicating the type of reaction. Allowable reaction types 
        # include (case insensitive): irreversible, reversible, reversible_forward,
        # reversible_backward, exchange, exchange_forward, exchange_backward 
        if not isinstance(type,str):
            raise TypeError('Invalid datatype for reaction type for reaction ' + self.id + '. String expected.')
        elif type.lower() not in ['irreversible','reversible','exchange','reversible_forward','reversible_backward','exchange_forward','exchange_backward']: 
            raise ValueError("Invalid reaction type for reaaction " + self.id+ "! Eligible choices are 'irreversible','reversible','exchange','reversible_forward','reversible_backward','exchange_forward','exchange_backward'") 
        else:
            self.type = type

        # A string containing the name of the reaction (complete or expanded)
        self.name = name

        # Name Synonyms
        self.synonyms = synonyms

        # Kegg id
        self.Kegg_id = Kegg_id

        # model Seed id
        self.ModelSeed_id = ModelSeed_id

        # EC number of the enzyme coding for the reaction
        self.EC_number = EC_number

        # A string containing the name of the reaction syb-system
        self.subsystem = subsystem

        # List of the gene objects coding for this reaction 
        self.genes = genes

        # A string containing the gene-reaction rules (genes are referred to by id)
        self.gene_reaction_rule = gene_reaction_rule 

        # The coefficinet of this reaction in the FBA objective function (real)
        self.objective_coefficient = objective_coefficient 

        # Flux of reaction under a specific condition (real)
        self.flux = flux

        # A list of the form [LB,UB] containing the flux bounds for this reaction 
        if flux_bounds == None: 
            self.assign_flux_bounds()
        else:
            self.flux_bounds = flux_bounds

        # A parameter showing whether or not to store the value of the flux after performing
        # FBA or other similar analysis 
        self.store_flux = store_flux

        # A tuple of the form (dGmin,dGmax) containing the min and max values of deltaG
        # (Gibbs free energy change) for this reaction
        self.deltaG_range = deltaG_range 

        # A stirng containing the reaction kinetics. See method kinetic_rate_calc
        # for more details (To be developed further in future!) 
        self.kinetics = kinetics

        # List of compounds appearing in the kinetic expression
        self.kinetic_compounds = kinetic_compounds

        # List of compound objects containing all compounds participating in this reaction
        self.compounds = sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] < 0 ],key=lambda x:x.id) + sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] > 0 ],key=lambda x:x.id) 

        # List of compound objects containing all reactants of this reaction
        self.reactants = sorted([m for m in self.stoichiometry.keys() if self.stoichiometry[m] < 0],key=lambda x:x.id)

        # List of compound objects containing all products of this reaction
        self.products = sorted([m for m in self.stoichiometry.keys() if self.stoichiometry[m] > 0],key=lambda x:x.id)

        # A list of compartments (each is an instance of the object compartment)
        # This is particularly hepful when a reaction occurs in more than one compartment
        # Here, we have a list instead of a single element because different compounds
        # in the reaction may participate in different compartment
        if compartment == []:
            self.compartment = list(set([c.compartment for c in self.compounds]))
        else:
            self.compartment = compartment

        # Confidence level for including this reaction in a model
        self.confidence_level = confidence_level

        # Notes and comments
        self.notes = notes

    def assign_flux_bounds(self):
        """
        Assigns general bounds to fluxes based on the reaction type 
        """
        if self.type.lower() == 'irreversible':
            self.flux_bounds = [0,1000]
        elif self.type.lower() == 'reversible':
            self.flux_bounds = [-1000,1000]
        elif self.type.lower() == 'reversible_forward':
            self.flux_bounds = [0,1000]
        elif self.type.lower() == 'reversible_backward':
            self.flux_bounds = [0,1000]
        elif self.type.lower() == 'exchange':
            self.flux_bounds = [0,1000]
        elif self.type.lower() == 'exchange_forward':
            self.flux_bounds = [0,1000]
        elif self.type.lower() == 'exchange_backward':
            self.flux_bounds = [0,0]
        else:
            raise userError('Unknown reaction type')

    def get_compounds(self,ref_type = 'id'):
        """
        Reports the list of compounds participating in this reaction in the output 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the compounds should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula. If name or id was not provided, id used instead.
        """
        if ref_type.lower() not in ['id','name','formula']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")

        for metab in self.compounds:
            if ref_type.lower() == 'id':
                mm = metab.id
            elif ref_type.lower() == 'name':
                if metab.name is None:
                    mm = metab.id
                else:
                    mm = metab.name
            elif ref_type.lower() == 'formula':
                if metab.formula is None:
                    mm = metab.id
                else:
                    mm = metab.formula

            return mm

    def get_reactants(self,ref_type = 'id'):
        """
        Reports the list of reactants of this reaction in the output 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the compounds should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula, If name or id was not provided, id used instead.

        """
        if ref_type.lower() not in ['id','name','formula']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")

        for metab in self.reactants: 
            if ref_type.lower() == 'id':
                mm = metab.id
            elif ref_type.lower() == 'name':
                if metab.name is None:
                    mm = metab.id
                else:
                    mm = metab.name
            elif ref_type.lower() == 'formula':
                if metab.formula is None:
                    mm = metab.id
                else:
                    mm = metab.formula

            return mm

    def get_products(self,ref_type = 'id'):
        """
        Reports the list of products of this reaction in the output 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the compounds should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula. If name or id was not provided, id used instead.
        """
        if ref_type.lower() not in ['id','name']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")

        for metab in self.products: 
            if ref_type.lower() == 'id':
                mm = metab.id
            elif ref_type.lower() == 'name':
                if metab.name is None:
                    mm = metab.id
                else:
                    mm = metab.name
            elif ref_type.lower() == 'formula':
                if metab.formula is None:
                    mm = metab.id
                else:
                    mm = metab.formula

            return mm

    def get_equation(self,ref_type = 'id'):
        """
        Prints in the output the reaction equation where compounds are referenced
        with the format specified in ref_type. If name or id was not provided, id used instead.
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the compounds should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula

        """
        if ref_type.lower() not in ['id','name','formula']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")
        
        # Type of the reaction arrow
        if self.type.lower() in ['reversible','exchange']:
            arrow_type = '<==>'
        else:
            arrow_type = '-->'

        # Loop over reactants
        reactants_eqn = ''
        for metab in self.reactants: 
            if ref_type.lower() == 'id':
                mm = metab.id
            elif ref_type.lower() == 'name':
                if metab.name is None:
                    mm = metab.id
                else:
                    mm = metab.name
            elif ref_type.lower() == 'formula':
                if metab.formula is None:
                    mm = metab.id
                else:
                    mm = metab.formula

            if self.reactants.index(metab) < len(self.reactants)-1: 
                reactants_eqn += str(-self.stoichiometry[metab]) + ' ' + mm + ' + ' 
            else:
                reactants_eqn += str(-self.stoichiometry[metab]) + ' ' + mm 
                
        # Loop over products
        products_eqn = ''
        for metab in self.products: 
            if ref_type.lower() == 'id':
                mm = metab.id
            elif ref_type.lower() == 'name':
                if metab.name is None:
                    mm = metab.id
                else:
                    mm = metab.name
            elif ref_type.lower() == 'formula':
                if metab.formula is None:
                    mm = metab.id
                else:
                    mm = metab.formula

            if self.products.index(metab) < len(self.products)-1: 
                products_eqn += str(self.stoichiometry[metab]) + ' ' + mm + ' + ' 
            else:
                products_eqn += str(self.stoichiometry[metab]) + ' ' + mm 

        return  reactants_eqn + ' ' + arrow_type + ' ' + products_eqn

    def add_compounds(self,compounds_list):
        """
        Adds new compounds to the reaction
      
        INPUTS:
        -------
        compounds_list: A dictionary whose keys are compounds objects and values are their
                        stoichiometric coefficinet in the reaction. The compound objects should
                        have already been defined and added to the model
        """
        if not isinstance(compounds_list,dict):
            raise TypeError('compounds_list must be a dictionary')
        
        # Return an error if there are new compound objects not already in the model 
        cmp_notInModel = [c for c in compounds_list.keys() if c not in self.model.compounds]
        if len(cmp_notInModel) > 0:
            raise ValueError('The following compounds not in the model ' + self.model.id + ': ' + str([c.id for c in cmp_notInModel])) 

        for cmp in compounds_list:
            if compounds_list[cmp] == 0 and self.warnings:
                print 'WARNING! Compound {} has a stoichiometric coefficinet of zero and was not added to reaction {}'.format(cmp.id,self.id)
            else:
                self.compounds.append(cmp)
                self.stoichiometry[cmp] = compounds_list[cmp]
                cmp.reactions.append(self)
                if compounds_list[cmp] < 0:
                    self.reactants.append(cmp)
                    cmp.reactant_reactions.append(self)
                elif compounds_list[cmp] > 0:
                    self.products.append(cmp)
                    cmp.product_reactions.append(self)

        self.compounds = sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] < 0 ],key=lambda x:x.id) + sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] > 0 ],key=lambda x:x.id) 
        self.reactants = sorted(self.reactants,key=lambda x:x.id) 
        self.products = sorted(self.products,key=lambda x:x.id) 

    def remove_compounds(self,compounds_list):
        """
        Removes compounds from the reaction
      
        INPUTS:
        -------
        compounds_list: A list of compound objects that must be removed from the reaction 
        """
        for cmp in compounds_list:
            if cmp not in set([c for c in self.compounds or c in self.reactants or c in self.products or c in self.stoichiometry.keys()]) and self.warnings:
                print 'WARNING! Compound {} not in the reaction {} compounds and cannot be removed from this reaction reaction'.format(cmp.id,self.id)
       
            if cmp in self.compounds:
                del self.compounds[self.compounds[cmp]]
                del cmp.reactions[cmp.reactions.index(self)]
            if cmp in self.reactants:
                del self.reactants[self.reactants[cmp]]
                del cmp.reactant_reactions[cmp.reactant_reactions.index(self)]
            if cmp in self.products:
                del self.products[self.products[cmp]]
                del cmp.product_reactions[cmp.product_reactions.index(self)]
            if cmp in self.stoichiometry.keys():
                del self.stoichiometry[cmp]

            if len(cmp.reactions) == 0:
                self.model.remove_compounds([cmp])

        self.compounds = sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] < 0 ],key=lambda x:x.id) + sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] > 0 ],key=lambda x:x.id) 
        self.reactants = sorted(self.reactants,key=lambda x:x.id) 
        self.products = sorted(self.products,key=lambda x:x.id) 

    def kinetic_rate_calc(self,assignLB = False, assignUB= False, assignFlux = False, concentration_key = None):
        """
        Computes the kinetic rate of the reaction. 
        
        INPUTS:
        ------
          assignUB: Assigns the computed kinetic rate to the upper bound
                    on reaction flux, if True, and does not do this if False
          assignLB: Assigns the "negative" of the computed kinetic rate (assuming that 
                    this rate is positive) to the lower bound on reaction flux
                    (e.g., for exchange reactions), if True, and does 
                    not do this if False
        assignFlux: Assigns the computed kinetic rate to the reaction flux
                    if True, and does not do this if False 
        concentration_key: The key to access the concetration in case compounds'
                    concentrations are stored in the form of a dictionary. For example, 
                    if the keys of this dictionary represent time points,  
                    then flux_key can refer to a speciic time point. 

        Note that concentrations must appear as 'C[compound_id]'
        in the kinetic expression to be able to compute the kinetic
        rate automatically (C will be in the form of a dictionary) 

        OUTPUT:
        ------
        kinetic_rate:  A global variable storing the kientic rate of the reaction
        """
        # A list of tuples, which is then converted to a dicitonary containing the list of
        # metabolites participating in the kinetic equation and their concentration
        C = []
        if self.kinetic_compounds == None:
            self.kinetic_compounds = self.compounds

        for m in self.kinetic_compounds:
            if m.concentration == None:
                raise userError("'concentration' for metaboltie " + m.id + ' is not defined')
            else:
                if concentration_key == None: 
                    if type(m.concentration) is int or type(m.concentration) is float:
                        C.append((m.id,m.concentration))
                    else:
                        raise TypeError('A number expected for concentraiton of compound ' + m.id + '. A ' + str(type(m.concentration)) + ' is given.')
                elif concentration_key != None:
                    if type(m.concentration) is dict:
                        C.append((m.id,m.concentration[concentration_key]))
                    else:
                        raise TypeError('A dictionary expected for concentraiton of compound ' + m.id + '. A ' + str(type(m.concentration)) + ' is given.')

        C = dict(C)
        exec 'self.kinetic_rate = ' +  self.kinetics        
        if assignLB == True:
            self.flux_bounds[0] = - self.kinetic_rate
        if assignUB == True:
            self.flux_bounds[1] = self.kinetic_rate
        if assignFlux == True:
            if concentration_key == None and type(self.flux) is not dict:
                self.flux = self.kinetic_rate            
            elif concentration_key == None and type(self.flux) is dict and self.warnings:
                print 'WARNING! Kineetic rate of the reaction will overwrite reaction.flux, which is the form of a dictionary.'
            elif concentration_key != None and type(self.flux) is dict:
                self.flux[concentration_key] = self.kinetic_rate            
            elif concentration_key != None and type(self.flux) is not dict:
                self.flux = self.kinetic_rate            
                 
