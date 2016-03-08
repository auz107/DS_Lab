from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
import compound 
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
        del_compounds: Removes compounds from the reaction
    kinetic_rate_calc: Computes the kinetic rate of the reaction

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 01-05-2016
    """
    def __init__(self, id, stoichiometry = {}, type = '', name = '', name_aliases = [], KEGG_id = [], ModelSEED_id = [], EC_numbers = [], subsystem = None, pathways = [], compartment = [], genes = [], gene_reaction_rule = '', objective_coefficient = None, flux = None, store_flux = True, flux_bounds = [], deltaG = None, deltaG_uncertainty = None, deltaG_range = [], kinetics = None, kinetic_compounds = None, confidence_level =None, notes = None, warnings = True): 

        # Warnings and messages in the standard output
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
        self.type = type

        # A string containing the name of the reaction (complete or expanded)
        self.name = name

        # Name Synonyms
        self.name_aliases = name_aliases

        # KEGG id
        self.KEGG_id = KEGG_id
        if isinstance(self.KEGG_id,str):
            self.KEGG_id = [self.KEGG_id]

        # ModelSEED id
        self.ModelSEED_id = ModelSEED_id
        if isinstance(self.ModelSEED_id,str):
            self.ModelSEED_id = [ModelSEED_id]
 
        # EC number of the enzyme coding for the reaction
        self.EC_numbers = EC_numbers

        # A list of strings containing the names of the reaction syb-system
        self.subsystem = subsystem

        # A list of strings containing the names of the pathways the reaciton is invovled in 
        self.pathways = pathways

        # List of the gene objects coding for this reaction 
        self.genes = genes

        # A string containing the gene-reaction rules (genes are referred to by id)
        self.gene_reaction_rule = gene_reaction_rule 

        # The coefficinet of this reaction in the FBA objective function (real)
        self.objective_coefficient = objective_coefficient 

        # Flux of reaction under a specific condition (real)
        self.flux = flux

        # A list of the form [LB,UB] containing the flux bounds for this reaction 
        if flux_bounds == []: 
            self.assign_flux_bounds()
        else:
            self.flux_bounds = flux_bounds

        # A parameter showing whether or not to store the value of the flux after performing
        # FBA or other similar analysis 
        self.store_flux = store_flux

        # deltaG of reaction
        self.deltaG = deltaG

        # Uncertainty in the deltaG of reaction
        self.deltaG_uncertainty = deltaG_uncertainty

        # A list of the form [dGmin,dGmax] containing the min and max values of deltaG
        # (Gibbs free energy change) for this reaction
        self.deltaG_range = deltaG_range 

        # Either a stirng containing the reaction kinetics or a python funciton defining the uptake kientics.
        # See method kinetic_rate_calc for more details. If it is a string the concentration of different compounds participating
        # in the kinetic equation must be in the form of a diciontary, e.g., '2*C[compound_id[/(25 + C[compound_id1])'. In case it is
        # a python function, the input to this funciton must be a dictionary with keys being compound ids and values their concentration
        self.kinetics = kinetics

        # List of compounds appearing in the kinetic expression
        self.kinetic_compounds = kinetic_compounds

        # This is particularly hepful when a reaction occurs in more than one compartment
        # Here, we have a list instead of a single element because different compounds
        # in the reaction may participate in different compartment
        self.compartment = compartment

        # Confidence level for including this reaction in a model
        self.confidence_level = confidence_level

        # Notes and comments
        self.notes = notes

        # Assign other properties
        self.assign_props() 
       
    def assign_props(self):
        """
        Assigns the properties of reaction
        """
        # A list of the form [dGmin,dGmax] containing the min and max values of deltaG
        # (Gibbs free energy change) for this reaction
        if self.deltaG_range == [] and self.deltaG != None and self.deltaG_uncertainty != None:
            self.deltaG_range = [min(self.deltaG - self.deltaG_uncertainty,self.deltaG + self.deltaG_uncertainty),max(self.deltaG - self.deltaG_uncertainty,self.deltaG + self.deltaG_uncertainty)]

        # List of compound objects containing all compounds participating in this reaction
        self.compounds = sorted(list(set([c for c in self.stoichiometry.keys() if self.stoichiometry[c] < 0] + [c for c in self.stoichiometry.keys() if self.stoichiometry[c] > 0])),key=lambda x:x.id) 

        # List of compound objects containing all reactants of this reaction
        self.reactants = sorted(list(set([m for m in self.stoichiometry.keys() if self.stoichiometry[m] < 0])),key=lambda x:x.id)

        # List of compound objects containing all products of this reaction
        self.products = sorted(list(set([m for m in self.stoichiometry.keys() if self.stoichiometry[m] > 0])),key=lambda x:x.id)

        # This is particularly hepful when a reaction occurs in more than one compartment
        # Here, we have a list instead of a single element because different compounds
        # in the reaction may participate in different compartment
        if self.compartment == []:
            self.compartment = list(set([c.compartment for c in self.compounds]))

        # model
        models = list(set([c.model for c in self.compounds if hasattr(c,'model') and c.model != None]))
        if len(models) == 1:
           self.model = models[0]
        elif len(models) > 1:
            raise userError('More than one model assigned to "compounds" of reaction ' + self.id + ': ' + str([m.id for m in models]))
        if len(models) == 0:
           self.model = None 

    def check_attr(self,attr_name,attr_value):
        """
        Checks the conditions on the class attributes
 
        INPUTS:
        -------
         attr_name: Attribute name
        attr_value: Attribute vlaue
        """
        # warnings 
        if attr_name == 'warnings' and not isinstance(attr_value,bool):
            raise TypeError("warnings must be True or False")

        # id 
        if attr_name == 'id' and not isinstance(attr_value,str):
            raise TypeError("Invalid reaction id " + str(attr_value) + "! reaction id must be a string. A " + str(attr_value) + " type object was entered instead")

        # Type 
        if attr_name == 'type' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'type' for reaaction " + self.id +"! reaction type must be a string. A " + str(attr_value) + " type object was entered instead")
        elif attr_name == 'type' and attr_value.lower() not in ['irreversible','reversible','exchange','reversible_forward','reversible_backward','exchange_forward','exchange_backward']: 
            raise ValueError("Invalid 'type' for reaaction " + self.id + "! Eligible choices are 'irreversible','reversible','exchange','reversible_forward','reversible_backward','exchange_forward','exchange_backward'") 

        # stoichiometry
        if attr_name == 'stoichiometry' and not isinstance(attr_value,dict):
            raise TypeError("Invalid 'stoichiometry' for reaaction " + self.id + "! 'stoichiometry' must be a dictionary. A " + str(attr_value) + " type object was entered instead")
        
        # Name
        if attr_name == 'name' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'name' for reaction " + self.id + "! 'name' must be a string. A " + str(attr_value) + " type object was entered instead")

        # Name aliases
        if attr_name == 'name_aliases' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'name_aliases' for reaction " + self.id + "! 'name_aliases'  must be a list of strings. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'name_aliases' and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError("Invalid 'name_aliases' for reaction " + self.id + "! 'name_aliases'  must be a list of strings. Objects that are not string found in the list:" + str([n for n in attr_value if not isinstance(n,str)]))

        # KEGG id
        if attr_name == 'KEGG_id' and (not isinstance(attr_value,str) and not isinstance(attr_value,list)):
            raise TypeError("Invalid 'KEGG_id' for reaction " + self.id + "! 'KEGG_id'  must be eitehr a string or a list of strings. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'KEGG_id' and isinstance(attr_value,list) and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError("Invalid 'KEGG_id' for reaction " + self.id + "! 'KEGG_id'  must be a list of strings. Objects that are not string found in the list: " + str([n for n in attr_value if not isinstance(n,str)]))

        # ModelSEED id
        if attr_name == 'ModelSEED_id' and (not isinstance(attr_value,str) and not isinstance(attr_value,list)):
            raise TypeError("Invalid 'ModelSEED_id' for reaction " + self.id + "! 'ModelSEED_id'  must be either a string or a list of strings. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'ModelSEED_id' and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError("Invalid 'ModelSEED_id' for reaction " + self.id + "! 'ModelSEED_id'  must be a list of strings. Objects that are not string found in the list:" + str([n for n in attr_value if not isinstance(n,str)]))

        # EC numbers of the enzyme coding for the reaction
        if attr_name == 'EC_numbers' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'EC_numbers' for reaction " + self.id + "! 'EC_numbers'  must be a list of strings. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'EC_numbers' and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError("Invalid 'EC_numbers' for reaction " + self.id + "! 'EC_numbers'  must be a list of strings. Objects that are not string found in the list:" + str([n for n in attr_value if not isinstance(n,str)]))

        # Subsystem
        if attr_name == 'subsystem' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'subsystem' for reaction " + self.id + "! 'subsystem'  must be a strings. A " + str(attr_value) + " type object was entered instead")

        # Pathways
        if attr_name == 'pathways' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'pathways' for reaction " + self.id + "! 'pathways'  must be a list of strings. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'pathways' and len([n for n in attr_value if not isinstance(n,str)]) > 0:
            raise TypeError("Invalid 'pathways' for reaction " + self.id + "! 'pathways' must be a list of strings. Objects that are not string found in the list: " + str([n for n in attr_value if not isinstance(n,str)]))

        # Compounds
        if attr_name == 'compounds' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'compounds' format for reaction " + self.id + "! 'compounds'  must be a list of objects of type compound. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'compounds' and len([n for n in attr_value if not isinstance(n,compound.compound)]) > 0:
            raise TypeError("Invalid 'compounds' format for reaction " + self.id + "! 'compounds'  must be a list of objects of type compound. Objects that are not of type compound found in the list: " + str([r for r in attr_value if not isinstance(n,compound.compound)]))

        # reactants
        if attr_name == 'reactants' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'reactants' format for reaction " + self.id + "! 'reactants'  must be a list of objects of type compound. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'reactants' and len([n for n in attr_value if not isinstance(n,compound.compound)]) > 0:
            raise TypeError("Invalid 'reactants' format for reaction " + self.id + "! 'reactants'  must be a list of objects of type compound. Objects that are not of type compound found in the list:" + str([r for r in attr_value if not isinstance(n,compound.compound)]))

        # products
        if attr_name == 'products' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'products' format for reaction " + self.id + "! 'products'  must be a list of objects of type compound. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'products' and len([n for n in attr_value if not isinstance(n,compound.compound)]) > 0:
            raise TypeError("Invalid 'products' format for reaction " + self.id + "! 'products'  must be a list of objects of type compound. Objects that are not of type compound found in the list: " + str([r for r in attr_value if not isinstance(n,compound.compound)]))

        # Compartments 
        if attr_name == 'compartments' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'compartments' format for reaction " + self.id + "! compartments for a model must be a list of objects of type compartment. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'compartments' and len([n for n in attr_value if not isinstance(n,compartment)]) > 0:
            raise TypeError("Invalid 'compartments' format for reaction " + self.id + "! compartments for a model must be a list of objects of type compartment. Objects that are not of type compartment found in the list: " + str([r for r in attr_value if not isinstance(n,compartment)]))

        # Genes 
        if attr_name == 'genes' and (attr_value is not None and not isinstance(attr_value,list)):
            raise TypeError("Invalid 'genes' format for reaction " + self.id + "! 'genes' for a model must be a list of objects of type gene. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'genes' and len([n for n in attr_value if not isinstance(n,gene)]) > 0:
            raise TypeError("Invalid 'genes' format for reaction " + self.id + "! 'genes' for a model must be a list of objects of type gene. Objects that are not of type gene found in the list: " + str([r for r in attr_value if not isinstance(n,gene)]))

        # gene_reaction_rule
        if attr_name == 'gene_reaction_rule' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'gene_reaction_rule' for reaction " + self.id + "! 'gene_reaction_rule' must be a string. A " + str(attr_value) + " type object was entered instead")

        # objective_coefficient
        if attr_name == 'objective_coefficient' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invlaud 'objective_coefficient' for reaction " + self.id + "! 'objective_coefficient'  must be either a float or an integer. A " + str(attr_value) + " type object was entered instead")

        # Flux bounds
        if attr_name == 'flux_bounds' and (attr_value is not None and not isinstance(attr_value,list)):
            raise TypeError("Invlaud 'flux_bounds' for reaction " + self.id + "! 'flux_bounds' must be a list of two integer or float elements. A " + str(type(attr_value)) + " type object was entered instead")
        if attr_name == 'flux_bounds' and (attr_value is not None and isinstance(attr_value,list) and attr_value[0] > attr_value[1]):
            raise ValueError("Lower bound is greater than upper bound for reaction " + self.id + "! flux_bounds = " + str(attr_value))

        # store_flux 
        if attr_name == 'store_flux' and not isinstance(attr_value,bool):
            raise TypeError("Invlaud 'store_flux' for reaction " + self.id + "! 'store_flux' must be True or False. A " + str(attr_value) + " type object was entered instead")

        # deltaG
        if attr_name == 'deltaG' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invlaud 'deltaG' for reaction " + self.id + "! 'deltaG' must be either a float or an integer. A " + str(attr_value) + " type object was entered instead")

        # deltaG_uncertainty
        if attr_name == 'deltaG_uncertainty' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invlaud 'deltaG_uncertainty' for reaction " + self.id + "! 'deltaG_uncertainty' must be either a float or an integer. A " + str(attr_value) + " type object was entered instead")

        # deltaG_range 
        if attr_name == 'deltaG_range' and not isinstance(attr_value,list):
            raise TypeError("Invlaud 'deltaG_range' for reaction " + self.id + "'deltaG_range' must be a listj. A " + str(attr_value) + " type object was entered instead")

        # Confidence level for including this reaction in a model
        if attr_name == 'confidence_level' and (attr_value is not None and not isinstance(attr_value,int) and not isinstance(attr_value,float)):
            raise TypeError("Invlaud 'confidence_level' for reaction " + self.id + "! 'confidence_level'  must be either a float or an integer. A " + str(attr_value) + " type object was entered instead")
        if attr_name == 'confidence_level' and attr_value is not None and  (attr_value < 0 or attr_value > 1):
            raise TypeError("Invlaud 'confidence_level' for reaction " + self.id + "! 'confidence_level' must be either a float or an integer between zero and one")

    def __setattr__(self,attr_name,attr_value):
       """
       Redefines funciton __setattr__
       INPUTS:
       -------
       attr_name: Attribute name
       attr_value: Attribute value
       """
       if attr_name in ['warnings','id','type','stoichiometry','type','name','name_aliases','KEGG_id','ModelSEED_id','EC_numbers','subsystem','pathways','reactants','reactants','products','compartmenets','genes','gene_reaction_rule','objective_coefficient','store_flux','deltaG','deltaG_uncertainty','deltaG_range','confidence_level']: 
           self.check_attr(attr_name,attr_value)
       self.__dict__[attr_name] = attr_value

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

    def add_compounds(self,new_compounds):
        """
        Adds new compounds to the reaction
      
        INPUTS:
        -------
        new_compounds: A dictionary whose keys are compounds objects and values are their
                       stoichiometric coefficinet in the reaction. The compound objects should
                       have already been defined and added to the model
        """
        if not isinstance(new_compounds,dict):
            raise TypeError('new_compounds must be a dictionary')

        # Values of the dictionary that are not integer or float
        non_int_float_values = [v for v in new_compounds.values() if not isinstance(v,int) and not isinstance(v,float)]
        if len(non_int_float_values) > 0:
            raise TypeError('The values of dictionary new_compounds must be an integer or float. The following non-integer and non-float values were observed: {}'.format(non_int_float_values))
       
        for cmp in new_compounds.keys():
            if new_compounds[cmp] == 0 and self.warnings:
                print 'WARNING (reaction.py)! Compound {} has a stoichiometric coefficinet of zero and was not added to reaction {}'.format(cmp.id,self.id)
            else:
                if cmp in self.compounds:
                    print 'WARNING (reaction.py)! Compound {} already participates in reaction {}. Its sstoichiometric coefficient in the reaction was updated.'.format(cmp.id,self.id)
                else:
                    self.compounds.append(cmp)
                self.stoichiometry[cmp] = new_compounds[cmp]
                if self not in cmp.reactions:
                    cmp.reactions.append(self)

                if new_compounds[cmp] < 0:
                    if cmp not in self.reactants:
                        self.reactants.append(cmp)
                    if self not in cmp.reactant_reactions:                    
                        cmp.reactant_reactions.append(self)

                if new_compounds[cmp] > 0:
                    if cmp not in self.products:
                        self.products.append(cmp)
                    if self not in cmp.product_reactions:
                        cmp.product_reactions.append(self)

        self.compounds = sorted(list(set([c for c in self.stoichiometry.keys() if self.stoichiometry[c] < 0])),key=lambda x:x.id) + sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] > 0 ],key=lambda x:x.id) 
        self.reactants = sorted(list(set(self.reactants)),key=lambda x:x.id) 
        self.products = sorted(list(set(self.products)),key=lambda x:x.id) 

        # Update the model if the reaction is associated with a model
        if self.model != None:
            # Compounds not in the model
            cmps_notInModel = [c for c in new_compounds.keys() if c not in self.model.compounds]
            if len(cmps_notInModel) > 0:
                self.model.add_compounds(cmps_notInModel)

        self.assign_props()

    def del_compounds(self,compounds_list):
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
                del self.compounds[self.compounds.index(cmp)]
                del cmp.reactions[cmp.reactions.index(self)]
            if cmp in self.reactants:
                del self.reactants[self.reactants.index(cmp)]
                del cmp.reactant_reactions[cmp.reactant_reactions.index(self)]
            if cmp in self.products:
                del self.products[self.products.index(cmp)]
                del cmp.product_reactions[cmp.product_reactions.index(self)]
            if cmp in self.stoichiometry.keys():
                del self.stoichiometry[cmp]

            if len(cmp.reactions) == 0 and self.model != None:
                self.model.del_compounds([cmp])

        self.compounds = sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] < 0 ],key=lambda x:x.id) + sorted([c for c in self.stoichiometry.keys() if self.stoichiometry[c] > 0 ],key=lambda x:x.id) 
        self.reactants = sorted(self.reactants,key=lambda x:x.id) 
        self.products = sorted(self.products,key=lambda x:x.id) 

    def kinetic_rate_calc(self,assignLB = False, assignUB= False, assignFlux = False, conc_key = None):
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
        conc_key: The key to access the concetration in case compounds'
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
                raise ValueError("'concentration' for compound " + m.id + ' is not defined')
            else:
                if conc_key == None: 
                    if type(m.concentration) is int or type(m.concentration) is float:
                        C.append((m.id,m.concentration))
                    else:
                        raise TypeError('A number expected for concentraiton of compound ' + m.id + '. A ' + str(type(m.concentration)) + ' is given.')
                else: # if conc_key is not known
                    if type(m.concentration) is dict:
                        C.append((m.id,m.concentration[conc_key]))
                    else:
                        raise TypeError('A dictionary expected for concentraiton of compound ' + m.id + '. A ' + str(type(m.concentration)) + ' is given.')

        C = dict(C)
        if isinstance(self.kinetics,str):
            exec 'self.kinetic_rate = ' +  self.kinetics        
        else: # if it is a function
            self.kinetic_rate = self.kinetics(C)

        if assignLB == True:
            self.flux_bounds[0] = - self.kinetic_rate
        if assignUB == True:
            self.flux_bounds[1] = self.kinetic_rate
        if assignFlux == True:
            if conc_key == None and type(self.flux) is not dict:
                self.flux = self.kinetic_rate            
            elif conc_key == None and type(self.flux) is dict and self.warnings:
                print 'WARNING! Kineetic rate of the reaction will overwrite reaction.flux, which is the form of a dictionary.'
            elif conc_key != None and type(self.flux) is dict:
                self.flux[conc_key] = self.kinetic_rate            
            elif conc_key != None and type(self.flux) is not dict:
                self.flux = self.kinetic_rate            
                 
