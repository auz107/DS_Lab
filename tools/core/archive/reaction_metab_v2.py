from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from metabolite import metabolite
from gene import gene
from compartment import compartment

class reaction(object):
    """
    A class holding the information about a biochemical reaciton 

    Methods:
         get_equation: Get the reaction equation
      get_metabolites: Prints in the output the list of all metabolites 
                       participating in this reaction as a reactant or product
        get_reactants: Prints in the output the list of all metabolites 
                       participating in this reaction as a reactan 
         get_products: Prints in the output the list of all metabolites 
                       participating in this reaction as a product
    kinetic_rate_calc: Computes the kinetic rate of the reaction

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 02-18-2015
    """
    def __init__(self, id, stoichiometry, type, name = '', Kegg_id = '', ModelSeed_id = '', EC_number = '', subsystem = '', compartments = None, genes = [], gene_reaction_rule = '', objective_coefficient = None, flux = None, flux_bounds = [], deltaG_range = [], kinetics = None, kinetic_metabolites = None, notes = ''): 
    
        # Reaction id or abbreviation in the model
        self.id = id

        # Reaction stoichiometry. A dictionary with keys and values as follows
        #   Keys: Instances of the object metabolite 
        # Values: Their stoichiometric coefficient in this reaction
        self.stoichiometry = stoichiometry
        
        # A string indicating the type of reaction. Allowable reaction types 
        # include (case insensitive): irreversible, reversible, reversible_forward,
        # reversible_backward, exchange, exchange_forward, exchange_backward 
        if type.lower() not in ['irreversible','reversible','exchange','reversible_forward','reversible_backward','exchange_forward','exchange_backward']: 
            raise userError("**Error! Invalid type! Eligible choices are 'irreversible','reversible','exchange','reversible_forward','reversible_backward','exchange_forward','exchange_backward'") 
        else:
            self.type = type

        # A string containing the name of the reaction (complete or expanded)
        self.name = name

        # Kegg id
        self.Kegg_id = Kegg_id

        # model Seed id
        self.ModelSeed_id = ModelSeed_id

        # EC number of the enzyme coding for the reaciton
        self.EC_number = EC_number

        # A string containing the name of the reaction syb-system
        self.subsystem = subsystem

        # A list of compartments (each is an instance of the object compartment)
        # Here, we have a list instead of a single element because different metabolites
        # in the reaction may participate in different compartments
        self.compartments = compartments

        # List of the gene objects coding for this reaction 
        self.genes = genes

        # A string containing the gene-reaciton rules (genes are referred to by id)
        self.gene_reaction_rule = gene_reaction_rule 

        # The coefficinet of this reaction in the FBA objective function (real)
        self.objective_coefficient = objective_coefficient 

        # Flux of reaction under a specific condition (real)
        self.flux = flux

        # A list of the form [LB,UB] containing the flux bounds for this reaction 
        if flux_bounds == []: 
            self._assignFluxBounds()
        else:
            self.flux_bounds = flux_bounds

        # A tuple of the form (dGmin,dGmax) containing the min and max values of deltaG
        # (Gibbs free energy change) for this reaction
        self.deltaG_range = deltaG_range 

        # A stirng containing the reaction kinetics. See method kinetic_rate_calc
        # for more details (To be developed further in future!) 
        self.kinetics = kinetics

        # List of metabolites appearing in the kinetic expression
        self.kinetic_metabolites = kinetic_metabolites

        # List of metabolite objects containing all metabolites participating in this reaction
        self.metabolites = self.stoichiometry.keys()

        # List of metabolite objects containing all reactants of this reaction
        self.reactants = [m for m in self.stoichiometry.keys() if self.stoichiometry[m] < 0]

        # List of metabolite objects containing all products of this reaction
        self.products = [m for m in self.stoichiometry.keys() if self.stoichiometry[m] > 0]

        # Notes and comments
        if isinstance(notes,str):
            self.notes = notes
        else:
            self.notes = ''

    def _assignFluxBounds(self):
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
            raise userError('Unknown reaciton type')

    def bolites(self,ref_type = 'id'):
        """
        Reports the list of metabolites participating in this reaciton in the output 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the metabolites should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula. If name or id was not provided, id used instead.
        """
        if ref_type.lower() not in ['id','name','formula']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")

        for metab in self.metabolites:
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
        Reports the list of reactants of this reaciton in the output 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the metabolites should be 
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
        Reports the list of products of this reaciton in the output 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the metabolites should be 
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
        Prints in the output the reaction equation where metabolites are referenced
        with the format specified in ref_type. If name or id was not provided, id used instead.
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the metabolites should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula

        """
        if ref_type.lower() not in ['id','name','formula']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")
        
        # Type of the reaciton arrow
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

    def kinetic_rate_calc(self,assignLB = False, assignUB= False, assignFlux = False, concentration_key = None):
        """
        Computes the kinetic rate of the reaction. 
        
        INPUTS:
        ------
          assignUB: Assigns the computed kinetic rate to the upper bound
                    on reaction flux, if True, and does not do this if False
          assignLB: Assigns the "negative" of the computed kinetic rate (assuming that 
                    this rate is positive) to the lower on reaction flux
                    (e.g., for exchange reactions), if True, and does 
                    not do this if False
        assignFlux: Assigns the computed kinetic rate to the reaction flux
                    if True, and does not do this if False 
        concentration_key: In case metabolite concentrations are stored in the form of a 
                    dictionary. For example, if the keys of this dictionary represent 
                    time points, then flux_key can refer to a speciic time point. 

        Note that concentrations must appear as 'C[metabolite_id]'
        in the kinetic expression to be able to compute the kinetic
        rate automatically (C will be in the form of a dictionary) 

        OUTPUT:
        ------
        kinetic_rate:  A global variable storing the kientic rate of the reaction
        """

        C = []
        if self.kinetic_metabolites == None:
            self.kinetic_metabolites = self.metabolites

        for m in self.kinetic_metabolites:
            if m.concentration == None:
                raise userError("**ERROR! 'concentration' for metaboltie " + m.id + ' is not defined')
            else:
                if concentration_key == None: 
                    if type(m.concentration) is int or type(m.concentration) is float:
                        C.append((m.id,m.concentration))
                    else:
                        raise userError('ERROR! A number expected for concentraiton of metabolite ' + m.id + '. A ' + str(type(m.concentration)) + ' is given.')
                elif concentration_key != None:
                    if type(m.concentration) is dict:
                        C.append((m.id,m.concentration[concentration_key]))
                    else:
                        raise userError('ERROR! A dictionary expected for concentraiton of metabolite ' + m.id + '. A ' + str(type(m.concentration)) + ' is given.')
                    

        C = dict(C)
        exec 'self.kinetic_rate = ' +  self.kinetics        
        if assignLB == True:
            self.flux_bounds[0] = - self.kinetic_rate
        if assignUB == True:
            self.flux_bounds[1] = self.kinetic_rate
        if assignFlux == True:
            if concentration_key == None and type(self.flux) is not dict:
                self.flux = self.kinetic_rate            
            elif concentration_key == None and type(self.flux) is dict:
                print 'WARNING! Kineetic rate of the reaciton will overwrite reaction.flux, which is the form of a dictionary.'
            elif concentration_key != None and type(self.flux) is dict:
                self.flux[concentration_key] = self.kinetic_rate            
            elif concentration_key != None and type(self.flux) is not dict:
                self.flux = self.kinetic_rate            
                 
