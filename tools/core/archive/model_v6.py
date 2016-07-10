from __future__ import division
import sys, time
from datetime import timedelta  # To convert elapsed time to hh:mm:ss format
from collections import Counter
sys.path.append('../../')
from tools.globalVariables import *
from tools.userError import userError
from organism import organism
from compound import compound
from compartment import compartment
from gene import gene
from reaction import reaction
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
from coopr.pyomo import *

class model(object):
    """
    This is a general class holding information for a model

    METHODS:
    --------
            assign_props: Assigns a number of properties to the model computed from existing inputs
           get_compounds: Returns a dictionary of the selected compound objects
           get_reactions: Returns a dictionary of the selected reaction objects
           add_compounds: Add new compound to the model
           add_reactions: Add new reactions to the model
           del_compounds: Remove compound from the model
           del_reactions: Add reactions to the model
                validate: Check for probable issues in the model and either fixes them
                          or issues a warning in the output 
       reset_flux_bounds: Restores the original bounds (bounds based on reaction types) for
                          all reactions in the model
         print_compounds: Prints in the output the list of all compounds in the model
                          with a desired format 
         print_reactions: Prints in the output the list of all reactions in the model with
                          a desired format 
        convert_to_cobra: Converts the model to COBRApy format 
            export_model: Export the model to various formats 
                     fba: Creates a fba model        


    Ali R. Zomorrodi, Segre Lab @ BU
    Last updated: 04-25-2016
    """
    def __init__(self, id, type = 'metabolic', organism = None, compounds = (), reactions = (), genes = (), compartments = (), name = None, biomass_reaction = None, atpm_reaction = None, notes = None, validate = True, stdout_msgs = True, warnings = True):

        # Warnings and messages in the standard output
        self.stdout_msgs = stdout_msgs
        self.warnings = warnings

        # Name of the model (string) 
        self.id = id 

        # Type of the model (string). Examples include cmpolic, cmpolic_expressions, genetic
        self.type = type

        # An organism object 
        self.organism = organism

        # A list of reaction objects 
        self.reactions = reactions

        # A list of compounds objects containing all compounds in the model 
        self.compounds = compounds

        # A list of gene objects containing all genes in the model
        self.genes = genes

        # A list of compartment objects containing compounds compartments 
        # for all compound in the model
        self.compartments = compartments

        # Model name
        self.name = name
      
        # Reaction object corresponding to biomass reaction. Note that some 
        # models may have more than one biomass reaction. For example, iAF1260 model 
        # of E. coli has a core and a wild-type biomass reaction. In that case the field
        # biomass_reaction must contain the reaction that is used as the objective function
        # of FBA. One can add a new instance variable to store all biomass reaction. For example,
        # model.all_biomass_reactions = {'core':model.get_reactions({'core_biomass_rxn_name':'id',
        # 'wildType_biomass_rxn_name':'id'})} (this must be done outside the class definition) 
        self.biomass_reaction = biomass_reaction

        # Reaciton object containing the ATP maintenance reaction
        self.atpm_reaction = atpm_reaction

        # Notes and comments
        self.notes = notes

        # Assign model properties using existing inputs
        self.assign_props()

        # Check for probable issues in the model
        if validate:
            self.validate(reassign_props = False)

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
       -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # Output messages and warnings 
        if attr_name in ['warnings','stdout_msgs'] and not isinstance(attr_value,bool):
            raise TypeError("{} must be True or False".format(attr_name))

        # id 
        if attr_name == 'id' and not isinstance(attr_value,str):
            raise TypeError("Invlaid id for model " + str(attr_value) + "! id must be a string. A " + str(attr_value) + " type object was entered instead")

        # Type 
        if attr_name == 'type' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'type' for model " + self.id + "! 'type' must be a string. A " + str(attr_value) + " type object was entered instead")

        # Organism 
        if attr_name == 'organism' and (attr_value is not None and not isinstance(attr_value,organism)):
            raise TypeError("Invalid 'organism' for model " + self.id + "! 'organism' must be an object of type organism. A " + str(attr_value) + " type object was entered instead. A " + str(attr_value) + " type object was entered instead")

        # Reactions by id
        if attr_name == 'reactions_by_id' and not isinstance(attr_value,dict):
            raise TypeError('A dictionary expected for reactions_by_id a {} provided instead'.format(type(attr_value)))
        elif attr_name == 'reactions_by_id' and len([r for r in attr_value.values() if not isinstance(r,reaction)]) > 0:
            invalid_rxn_objs = [(rid,attr_value[rid],type(attr_value[rid])) for rid in attr_value.keys() if not isinstance(attr_value[rid],reaction)]
            if len(invalid_rxn_objs) <= 10:
                raise TypeError("Invalid value for 'reactions_by_id' for model {}! Values of 'reactions_by_id' must be a list of objects of type 'reaction'. Objects of other types found instead in the following ten reactions: {}".format(self.id, invalid_rxn_objs))
            else:
                raise TypeError("Invalid value for 'reactions_by_id' for model {}! Values of 'reactions_by_id' must be a list of objects of type 'reaction'. Objects of other types found instead in the following ten reactions: {} and {} more".format(self.id, invalid_rxn_objs, len(invalid_rxn_objs) - 10))

        # Compartments 
        if attr_name == 'compartments' and not isinstance(attr_value,list):
            raise TypeError("Invalid 'compartments' for model " + self.id + "! 'compartments' for a model must be a list of objects of type compartment. A " + str(attr_value) + " type object was entered instead") 
        if attr_name == 'compartments' and len([r for r in attr_value if not isinstance(r,compartment)]) > 0: 
            raise TypeError("Invalid 'compartments' for model " + self.id + "! 'compartments' for a model must be a list of objects of type compartment. Objects that are not of type compartment found in the list: " + str([r for r in attr_value if not isinstance(r,compartment)]))

        # Model name
        if attr_name == 'name' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'name' for model " + self.id + "! 'name' must be a string. A " + str(attr_value) + " type object was entered instead")

        # Biomass reaction 
        if attr_name == 'biomass_reactions' and (attr_value is not None and not isinstance(attr_value,reaction)): 
            raise TypeError("Invalid 'biomass_reaction' for model " + self.id + "! 'biomass_reaction' must be of type reaction. A " + str(attr_value) + " type object was entered instead")

        # ATPM reaction 
        if attr_name == 'atpm_reactions' and (attr_value is not None and not isinstance(attr_value,reaction)): 
            raise TypeError("Invalid 'atpm_reactions' for model " + self.id + "! 'atpm_reaction' must be of type reaction. A " + str(attr_value) + " type object was entered instead")


        if attr_name == 'compounds':
            self.set_compounds(compounds = attr_value)
        elif attr_name == 'reactions':
            self.set_reactions(reactions = attr_value)
        elif attr_name == 'genes':
            self.set_genes(genes = attr_value)
        else: 
            self.__dict__[attr_name] = attr_value

    def assign_props(self):
        """
        Assigns a number of model properties from existing inputs 
        """
        # Model name
        if self.name == None:
            self.name = self.id

        #-- organism --
        self.organism.model_id = self.id

        #-- Compartments --
        if self.compartments == []:
            self.compartments = sorted(list(set([c.compartment for c in self.compounds if c.compartment != None])), key = lambda x: x.id)
        for cpt in self.compartments:
            cpt.model_id = self.id
        self.compartments_by_id = dict([(cpt.id,cpt) for cpt in self.compartments])    

        #-- Compounds --
        if self.compounds == []:
            self.compounds = sorted(tuple(set([c for r in self.reactions for c in r.stoichiometry.keys()], key = lambda x: x.id)))
        # A dictionary whose keys are compound ids and values are compound objects
        self.compounds_by_id = dict([(cpd.id,cpd) for cpd in self.compounds])    
        self.compounds_by_clean_id = dict([(remove_non_alphanumeric(cpd.id).lower(),cpd) for cpd in self.compounds])    

        #-- Reactions --
        if self.reactions == []:
            self.reactions = sorted(tuple(set([r for c in self.compounds for r in c.reactions])), key = lambda x: x.id)
        # A dictionary whose keys are reaciotn ids and values are reaction objects
        self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    
        self.reactions_by_clean_id = dict([(remove_non_alphanumeric(rxn.id).lower(),rxn) for rxn in self.reactions])    

        #-- Genes --
        if self.genes == []:
            self.genes = [g for r in self.reactions for g in r.genes] 

    def set_compounds(self, compounds): 
        """
        Set attribute compounds
        """
        if not isinstance(compounds,tuple) and not isinstance(compounds,list):
            raise TypeError("Invalid 'compounds' format for model {}! Compounds must be a tuple or a list but a {} object was provided instead".format(self.id, type(compounds)))
        if len([n for n in compounds if not isinstance(n,compound)]) > 0:
            raise TypeError("Invalid 'compounds' format for model {}! Compounds must be a list of 'compound' object but objects of {} were observed in the list instead. ".format(self.id, list(set([type(n) for n in compounds if not isinstance(n,compound.compound)]))))

        compounds = sorted(tuple(set(compounds)), key = lambda x: x.id)

        if isinstance(compounds,tuple):
            self.__dict__['compounds'] = compounds
        elif isinstance(compounds,list):
            self.__dict__['compounds'] = tuple(compounds)

        for m in self.compounds:
            m.model_id = self.id

        # A dictionary whose keys are compound ids and values are compound objects
        self.compounds_by_id = dict([(cpd.id,cpd) for cpd in self.compounds])    
        self.compounds_by_clean_id = dict([(remove_non_alphanumeric(cpd.id).lower(),cpd) for cpd in self.compounds])    

        cpds_reactions = dict((cpd,[]) for cpd in self.compounds)
        cpds_reactant_reactions = dict((cpd,[]) for cpd in self.compounds)
        cpds_product_reactions = dict((cpd,[]) for cpd in self.compounds)
        for rxn in self.reactions:
            for cpd in rxn.compounds: 
                cpds_reactions[cpd].append(rxn) 
            for cpd in rxn.reactants: 
                cpds_reactant_reactions[cpd].append(rxn) 
            for cpd in rxn.products: 
                cpds_product_reactions[cpd].append(rxn)
        for cpd in self.compounds:
            cpd.set_reactions(reactions = cpds_reactions[cpd])
            cpd.set_reactant_reactions(reactant_reactions = cpds_reactant_reactions[cpd])
            cpd.set_product_reactions(product_reactions = cpds_product_reactions[cpd])

    def set_reactions(self, reactions): 
        """
        Set attribute reactions
        """
        if reactions is not None and not isinstance(reactions,tuple) and not isinstance(reactions,list):
            raise TypeError("Invalid 'reactions' for compound " + str(self.id) + "! 'reactions'  must be a tuple or a list of objects of type reaction. A " + str(type(reactions)) + " type object was entered instead")
        if len([r for r in reactions if not isinstance(r,reaction)]) > 0:
            raise TypeError("Invalid 'reactions' for compound " + str(self.id) + "! 'reactions'  must be a list of objects of type 'reaction'.  Objects that are not of type reaction found in the list:" + str([n for n in reactions if not isinstance(n,reaction.reaction)]))

        reactions = sorted(tuple(set(reactions)), key = lambda x: x.id)

        if not isinstance(reactions,tuple):
            self.__dict__['reactions'] = reactions
        elif not isinstance(reactions,list):
            self.__dict__['reactions'] = list(reactions)

        for r in self.reactions:
            r.model_id = self.id

        # A dictionary whose keys are reaction ids and values are reaction objects
        self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    
        self.reactions_by_clean_id = dict([(remove_non_alphanumeric(rxn.id).lower(),rxn) for rxn in self.reactions])    

    def set_genes(self, genes): 
        """
        Set attribute genes
        """
        if not isinstance(genes, tuple) and  not isinstance(genes,list):
            raise TypeError("Invalid 'genes' for model " + self.id + "! 'genes' for a model must be a list of objects of type gene. A " + str(genes) + " type object was entered instead")
        if len([g for g in genes if not isinstance(g,gene)]) > 0:
            raise TypeError("Invalid 'compounds' for model " + self.id + "! 'genes' for a model must be a list of objects of type gene. Objects that are not of type gene found in the list: " + str([r for r in genes if not isinstance(r,gene)]))

        genes = sorted(tuple(set(genes)), key = lambda x: x.id)

        if not isinstance(genes, tuple):
            self.__dict__['genes'] = genes
        elif not isinstance(genes, list):
            self.__dict__['genes'] = tuple(genes)

        for gn in self.genes:
            gn.model_id = self.id

        # A dictionary whose keys are gene ids and values are gene objects
        self.genes_by_id = dict([(g.id,g) for g in self.genes])    

        if hasattr(self,'reactions'):
            for gen in [g for g in self.genes if g.reactions == []]:
                gen.reactions = [r for r in self.reactions if gen in r.genes]

    def get_reactions(self,reactions_ref,search_by_clean_ref = False):
        """
        Returns a dictionary of the selected reaction objects

        INPUTS:
        ------- 
              reactions_ref: A dictionary with keys and values as follows:
                               keys: A reaction id or name (string) 
                             values: A string indicating the format by which the reactions
                                     should be searched. Eligible choices include
                                     'id': Reaction id 
                                     'name': Reaction name
                                     'Kegg_id': Reaction Kegg id
                                     'ModelSeed_id: Reaction ModelSeed id 
                             In addition, since the same reaction (with the same name, Kegg_id,
                             ModelSeed_id and name)  may appear in more than one compartment
                             the following identifiers are defined by connecting name, Kegg_id ,
                             and ModelSeed_id with their compartment ids using an underline
                             'name_compart', 'Kegg_id_compart' and 'ModelSeed_id_compart'
        search_by_clean_ref: If True removes all non-alpphanumeric charcaterz from the names 
                             before searching

        Example: {'rxn1_id':'id','rxn2_name':'name'}
 
        OUTPUTS:
        -------
        If the input is a dictionary with a single key the the ouput is the corresponding
        reaction object, otherwise the output is a dictionary whose keys are the reaction 
        reference values (id or name) and values are reaction objects. None is returned if 
        a reaction is not found
        """
        reaction_objects = []

        if not isinstance(reactions_ref,dict):
            raise TypeError('reactions_ref must be a dictionary. A ' + str(type(reactions_ref)) + ' type object was entered instead')

        for rxn_ref in reactions_ref.keys(): 
            if reactions_ref[rxn_ref].lower() not in ['id','name','modelseed_id','kegg_id','modelseed_id_compart','kegg_id_compart','name_compart']:
                    raise ValueError("**Error! Invalid reaction reference '" + reactions_ref[rxn_ref] + "'")

            if search_by_clean_ref:
                rxn_object = []
                if reactions_ref[rxn_ref].lower() == 'id':
                    try:
                        rxn_object = [self.reactions_by_clean_id[remove_non_alphanumeric(rxn_ref).lower()]]
                    except:
                        rxn_object = []
                elif reactions_ref[rxn_ref].lower() == 'name':
                    rxn_object = [r for r in self.reactions if remove_non_alphanumeric(r.name).lower() == remove_non_alphanumeric(rxn_ref).lower() or remove_non_alphanumeric(rxn_ref).lower() in [remove_non_alphanumeric(n) for n in r.synonyms]]
                elif reactions_ref[rxn_ref].lower() == 'modelseed_id':
                    rxn_object = [r for r in self.reactions if r.ModelSeed_id != None and isinstance(r.ModelSeed_id,str) and r.ModelSeed_id.lower() == rxn_ref.lower()]
                elif reactions_ref[rxn_ref].lower() == 'kegg_id':
                    rxn_object = [r for r in self.reactions if r.Kegg_id != None and isinstance(r.Kegg_id,str) and r.Kegg_id.lower() == rxn_ref.lower()]
                elif reactions_ref[rxn_ref].lower() == 'name_compart':
                    rxn_object = [r for r in self.reactions if len(r.compartment) == 1 and remove_non_alphanumeric(r.name + '_' + r.compartment[o].id).lower() == remove_non_alphanumeric(rxn_ref).lower()]
                elif reactions_ref[rxn_ref].lower() == 'modelseed_id_compart':
                    rxn_object = [r for r in self.reactions if len(r.compartment) == 1 and remove_non_alphanumeric(r.ModelSeed_id + '_' + r.compartment[o].id).lower() == remove_non_alphanumeric(rxn_ref).lower()]
                elif reactions_ref[rxn_ref].lower() == 'kegg_id_compart':
                    rxn_object = [r for r in self.reactions if len(r.compartment) == 1 and remove_non_alphanumeric(r.Kegg_id + '_' + r.compartment[o].id).lower() == remove_non_alphanumeric(rxn_ref).lower()]

            # Search by actual (not clean) names
            else:
                rxn_object = []
                if reactions_ref[rxn_ref].lower() == 'id':
                    try:
                        rxn_object = [self.reactions_by_id[rxn_ref]]
                    except:
                        rxn_object = []
                elif reactions_ref[rxn_ref].lower() == 'name':
                    rxn_object = [r for r in self.reactions if r.name.lower() == rxn_ref.lower() or r.name.lower() in [n.lower() for n in r.synonyms]]
                elif reactions_ref[rxn_ref].lower() == 'modelseed_id':
                    rxn_object = [r for r in self.reactions if r.ModelSeed_id != None and isinstance(r.ModelSeed_id,str) and r.ModelSeed_id.lower() == rxn_ref.lower()]
                elif reactions_ref[rxn_ref].lower() == 'kegg_id':
                    rxn_object = [r for r in self.reactions if r.Kegg_id != None and isinstance(r.Kegg_id,str) and r.Kegg_id.lower() == rxn_ref.lower()]
                elif reactions_ref[rxn_ref].lower() == 'name_compart':
                    rxn_object = [r for r in self.reactions if len(r.compartment) == 1 and r.name + '_' + r.compartment[o].id.lower() == rxn_ref.lower()]
                elif reactions_ref[rxn_ref].lower() == 'modelseed_id_compart':
                    rxn_object = [r for r in self.reactions if len(r.compartment) == 1 and r.ModelSeed_id + '_' + r.compartment[o].id.lower() == rxn_ref.lower()]
                elif reactions_ref[rxn_ref].lower() == 'kegg_id_compart':
                    rxn_object = [r for r in self.reactions if len(r.compartment) == 1 and r.Kegg_id + '_' + r.compartment[o].id.lower() == rxn_ref.lower()]

            if len(rxn_object) > 1:
                reaction_objects.append((rxn_ref,rxn_object)) 
                if self.warnings:
                    print '\nWARNING! More than one reaction object found for ',rxn_ref,': ',[r.id for r in rxn_object]
            elif len(rxn_object) == 0:
                if self.warnings:
                    print '\nWARNING! No reaction object found for ' + rxn_ref
                reaction_objects.append((rxn_ref,None)) 
            elif len(rxn_object) == 1:
                reaction_objects.append((rxn_ref,rxn_object[0])) 

        if len(reactions_ref.keys()) == 1:
            output_reactions = reaction_objects[0][1]
        else:
            output_reactions = dict(reaction_objects)

        return output_reactions 


    def get_compounds(self,compounds_ref, search_by_clean_ref = True):
        """
        Returns a dictionary of the selected compound objects

        INPUTS:
        ------- 
              compounds_ref: A dictionary with keys and values as follows:
                              keys: A compound id or name (string) 
                             values: A string indicating the format by which the compound
                                     should be search. Eligible choices include 
                                     'id': Compound id 
                                     'name': Compound name
                                     'Kegg_id': Compound Kegg id
                                     'ModelSeed_id: Compound ModelSeed id 
                                     'formula': Compound formula
                             In addition, since the same metabolite (with the same name, Kegg_id,
                             ModelSeed_id and formula)  may appear in more than one compartment
                             the following identifiers are defined by connecting name, Kegg_id,
                             ModelSeed_id and formula with their compartment ids using an underline
                             'name_compart', 'Kegg_id_compart', 'ModelSeed_id_compart' and 'formula_compart'
        search_by_clean_ref: If True removes all non-alpphanumeric charcaterz from the names 
                             before searching

        Example: {'m1_id':'id','m2_name':'name','m3':'formula','h2o_c':'name_compart'}
        
        OUTPUTS:
        -------
        If the input is a dictionary with a single key the the ouput is the corresponding
        compound object, otherwise the output is a dictionary whose keys are the compound 
        reference values (id or name) and values are compound objects
        """
        compound_objects = []
        
        for cmp_ref in compounds_ref.keys(): 
            cmp_object = []
            # Search by clean reference
            if search_by_clean_ref:
                if compounds_ref[cmp_ref].lower() == 'id':
                    try:
                        cmp_object = [self.compounds_by_clean_id[remove_non_alphanumeric(cmp_ref).lower()]]
                    except:
                        cmp_object = []
                elif compounds_ref[cmp_ref].lower() == 'name':
                    cmp_object = [m for m in self.compounds if isinstance(m.name,str) and remove_non_alphanumeric(m.name).lower() == remove_non_alphanumeric(cmp_ref).lower() or isinstance(m.name,list) and remove_non_alphanumeric(cmp_ref).lower() in [remove_non_alphanumeric(n) for n in m.name]]
                elif compounds_ref[cmp_ref].lower() == 'name_compart':
                    cmp_object = [m for m in self.compounds if isinstance(m.name,str) and remove_non_alphanumeric(m.name + '_' + m.compartment.id).lower() == remove_non_alphanumeric(cmp_ref).lower() or isinstance(m.name,list) and remove_non_alphanumeric(cmp_ref).lower() in [remove_non_alphanumeric(n + '_' + m.compartment.id).lower() for n in m.name]]
                elif compounds_ref[cmp_ref].lower() == 'kegg_id':
                    cmp_object = [m for m in self.compounds if m.Kegg_id != None and m.Kegg_id.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'kegg_id_compart':
                    cmp_object = [m for m in self.compounds if m.Kegg_id != None and remove_non_alphanumeric(m.Kegg_id + '_' + m.compartment).id.lower() == remove_non_alphanumeric(cmp_ref).lower()]
                elif compounds_ref[cmp_ref].lower() == 'modelseed_id':
                    cmp_object = [m for m in self.compounds if m.ModelSeed_id != None and m.ModelSeed_id.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'modelseed_id_compart':
                    cmp_object = [m for m in self.compounds if m.ModelSeed_id != None and remove_non_alphanumeric(m.ModelSeed_id + '_' + m.compartment.id).lower() == remove_non_alphanumeric(cmp_ref).lower()]
                elif compounds_ref[cmp_ref].lower() == 'formula':
                    cmp_object = [m for m in self.compounds if m.formula != None and remove_non_alphanumeric(m.formula).lower() == remove_non_alphanumeric(cmp_ref).lower()]
                elif compounds_ref[cmp_ref].lower() == 'formula_compart':
                    cmp_object = [m for m in self.compounds if m.formula != None and remove_non_alphanumeric(m.formula + '_' + m.compartment.id).lower() == remove_non_alphanumeric(cmp_ref).lower()]
                else:
                    raise ValueError("Invalid compound reference (" + cmp_ref + ")")

            # Search with actual (not clean) reference
            else: 
                if compounds_ref[cmp_ref].lower() == 'id':
                    try:
                        cmp_object = [self.compounds_by_id[cmp_ref]]
                    except:
                        cmp_object = [] 
                elif compounds_ref[cmp_ref].lower() == 'name':
                    cmp_object = [m for m in self.compounds if isinstance(m.name,str) and m.name.lower() == cmp_ref.lower() or isinstance(m.name,lsit) and mp_ref.lower() in [n.lower() for n in m.name]]
                elif compounds_ref[cmp_ref].lower() == 'name_compart':
                    cmp_object = [m for m in self.compounds if isinstance(m.name,str) and m.name.lower() + '_' + m.compartment.id.lower() == cmp_ref.lower() or isinstance(m.name,list) and cmp_ref.lower() in [n.lower() + '_' + m.compartment.id.lower() for n in m.name]]
                elif compounds_ref[cmp_ref].lower() == 'kegg_id':
                    cmp_object = [m for m in self.compounds if m.Kegg_id != None and m.Kegg_id.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'kegg_id_compart':
                    cmp_object = [m for m in self.compounds if m.Kegg_id != None and m.Kegg_id.lower() + '_' + m.compartment.id.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'modelseed_id':
                    cmp_object = [m for m in self.compounds if m.ModelSeed_id != None and m.ModelSeed_id.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'modelseed_id_compart':
                    cmp_object = [m for m in self.compounds if m.ModelSeed_id != None and m.ModelSeed_id.lower() + '_' + m.compartment.id.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'formula':
                    cmp_object = [m for m in self.compounds if m.formula != None and m.formula.lower() == cmp_ref.lower()]
                elif compounds_ref[cmp_ref].lower() == 'formula_compart':
                    cmp_object = [m for m in self.compounds if m.formula != None and m.formula.lower() + '_' + m.compartment.id.lower() == cmp_ref.lower()]
                else:
                    raise ValueError("Invalid compound reference (" + cmp_ref + ")")

    
            if len(cmp_object) == 1:
                compound_objects.append((cmp_ref,cmp_object[0])) 
            elif len(cmp_object) > 1:
                if self.warnings:
                    print 'WARNING! More than one compound object found for ', cmp_ref,': ', [m.id for m in cmp_object]
                compound_objects.append((cmp_ref,cmp_object[0])) 
            elif len(cmp_object) == 0:
                if self.warnings:
                    print 'WARNING! No compound object found for ', cmp_ref,': ', [m.id for m in cmp_object]
                compound_objects.append((cmp_ref,None)) 

        if len(compounds_ref.keys()) == 1:
            if compound_objects[0] != None:
                output_compounds = compound_objects[0][1]
            else:
                output_compounds = None 
        else:
            output_compounds = dict(compound_objects)

        return output_compounds 

    def get_genes(self,genes_ref):
        """
        Returns a dictionary of the selected gene objects

        INPUTS:
        ------- 
        genes_ref: A dictionary with keys and values as follows:
                         keys: A gene id or name (string) 
                       values: A string indicating the format by which the genes
                               should be search. Eligible choices include 'id', 'name'
                               'Kegg_id' or 'ModelSeed_id'
                      Example: {'gn1_id':'id','gn2_name':'name'}
        
        OUTPUTS:
        -------
        If the input is a dictionary with a single key the the ouput is the corresponding
        gene object, otherwise the output is a dictionary whose keys are the gene 
        reference values (id or name) and values are gene objects
        """
        gene_objects = []
        
        for gn_ref in genes_ref.keys(): 
            gn_object = []
            if genes_ref[gn_ref].lower() == 'id':
                gn_object = [g for g in self.genes if g.id.lower() == gn_ref.lower()]
            elif genes_ref[gn_ref].lower() == 'name':
                gn_object = [g for g in self.genes if g.name.lower() == gn_ref.lower()]
            elif genes_ref[gn_ref].lower() == 'Kegg_id':
                gn_object = [g for g in self.genes if g.Kegg_id.lower() == gn_ref.lower()]
            elif genes_ref[gn_ref].lower() == 'ModelSeed_id':
                gn_object = [g for g in self.genes if g.ModelSeed_id.lower() == gn_ref.lower()]
            else:
                raise userError("**Error! Invalid gene reference '" + genes_ref[gn_ref] + "'")

            if len(gn_object) > 1:
                raise userError('More than one gene object found for ',gn)
            else:
                gene_objects.append((gn_ref,gn_object[0])) 

        if len(genes_ref.keys()) == 1:
            output_genes = gene_objects[0][1]
        else:
            output_genes = dict(gene_objects)

        return output_genes 

    def add_compounds(self,new_compounds):
        """
        Add new compounds to the model
        
        INPUTS:
        -------
        new_compounds: A list of new compound objects to be added to the model
        
        OUTPUTS:
        --------
        There is no actual output. The model is updated with the addition of
        new compounds 
        """
        # First check if the id of any new compounds is the same as that an existing compound
        # in the model
        if len(set(self.compounds_by_id.keys()).intersection(set([c.id for c in new_compounds]))) > 0:
            raise userError('The following compound ids already exist in the model: {}'.format(set(self.compounds_by_id.keys()).intersection(set([c.id for c in new_compounds]))))

        self.compounds += new_compounds

        # Sort according to id
        self.compounds = sorted(tuple(set(self.compounds)),key=lambda x:x.id)

        # Update compounds by id
        for new_cmp in  new_compounds:   
            self.compounds_by_id[new_cmp.id] = new_cmp
            self.compounds_by_clean_id[remove_non_alphanumeric(new_cmp.id).lower()] = new_cmp

        # Assign model 
        for new_cmp in  new_compounds:   
            new_cmp.model_id = self.id

    def add_reactions(self,new_reactions):
        """
        Add new reactions to the model
        
        INPUTS:
        -------
        new_reactions: A list of new reaction objects. If these new reactions contain
                       any new compounds, these compounds must be added to the model first
        
        OUTPUTS:
        --------
        There is no actual output. The model is updated with the addition of
        new reactions 
        """
        # First check if the id of any new reactions is the same as that an existing reaction
        # in the model
        if len(set(self.reactions_by_id.keys()).intersection(set([r.id for r in new_reactions]))) > 0:
            raise userError('The following reaction ids already exist in the model: {}'.format(set(self.reactions_by_id.keys()).intersection(set([r.id for r in new_reactions]))))

        self.reactions += new_reactions

        # Sort according to id
        self.reactions = sorted(tuple(set(self.reactions)),key=lambda x:x.id)

        # Update reactions_by_id
        for new_rxn in new_reactions:
            self.reactions_by_id[new_rxn.id] = new_rxn 
            self.reactions_by_clean_id[remove_non_alphanumeric(new_rxn.id).lower()] = new_rxn 

        # Assign model 
        for new_rxn in new_reactions:
            new_rxn.model_id = self.id

        # Replace compounds in new_reactions whose id is the same as the id of on existing 
        # compound in the model with the corresponding compound object in the model
        for rxn in new_reactions:
            rxn_stoic = {}
            for cpd in rxn.compounds:
                if cpd.id in self.compounds_by_id.keys():
                    rxn_stoic[self.compounds_by_id[cpd.id]] = rxn.stoichiometry[cpd]
                else:
                    rxn_stoic[cpd] = rxn.stoichiometry[cpd]
            rxn.set_stoichiometry(stoichiometry = rxn_stoic, replace = True)

        # Add any compounds in the reaction not in the model
        cpds_notInModel = [c for r in new_reactions for c in r.compounds if c not in self.compounds]
        if len(cpds_notInModel) > 0:
            self.add_compounds(cpds_notInModel) 

        # Update compound objects participating in these new reactions
        for rxn in new_reactions:
            for cpd in rxn.compounds:
                if rxn not in cpd.reactions:
                    cpd.reactions.append(rxn)
                if rxn.stoichiometry[cpd] < 0 and rxn not in cpd.reactant_reactions:
                    cpd.reactant_reactions.append(rxn) 
                if rxn.stoichiometry[cpd] > 0 and rxn not in cpd.product_reactions:
                    cpd.product_reactions.append(rxn) 

    def del_compounds(self,compounds_list):
        """
        Remove selected compounds from the model
        
        INPUTS:
        ------
        compounds_list: A list of compound objectc that must be removed 
        """
        if type(compounds_list) is not list:
            userError("**Error! the input to 'remove_compounds' should be a list of compound objects")

        # First remove this compound from all relevant reaction objects
        for cmp in compounds_list:
            for rxn in cmp.reactions:
                if rxn.stoichiometry[cmp] < 0:
                    if rxn.reactants != []:
                        del rxn.reactants[rxn.reactants.index(cmp)]
                elif rxn.stoichiometry[cmp] > 0:
                    if rxn.products != []:
                        del rxn.products[rxn.products.index(cmp)]

                del rxn.stoichiometry[cmp]

            # Check its presence in all other reactions
            cmp_rxns = list(set([r for r in self.reactions if cmp in r.compounds + r.reactants + r.products + r.stoichiometry.keys()])) 
            if len(cmp_rxns) > 0:
                for r in cmp_rxns:
                    try:
                        del r.compounds[r.compounds.index(cmp)]
                    except:
                        continue
                    try:
                        del r.reactants[r.reactant.index(cmp)]
                    except:
                        continue
                    try:
                        del r.products[r.products.index(cmp)]
                    except:
                        continue
                    try:
                        del r.stoichiometry.keys()[r.stoichiometry.keys().index(cmp)]
                    except:
                        continue

            # Remove it from the list of compounds
            del self.compounds[self.compounds.index(cmp)] 
            self.compounds = sorted(set(list(self.compounds)),key=lambda x:x.id)

            # Remove it from compounds_by_id
            del self.compounds_by_id[cmp.id]
            del self.compounds_by_clean_id[remove_non_alphanumeric(cmp.id).lower()]

    def del_reactions(self,reactions_list):
        """
        Remove selected reactions from the model
        
        INPUTS:
        ------
        reactions_list: A list of reaction objectc that must be removed 
        """
        if type(reactions_list) is not list:
            userError("**Error! the input to 'remove_reactions' should be a list of reaction objects")

        model_rxns = list(model.reactions)
        for rxn in reactions_list:
            # Remove it from the list of reactions
            del model_rxns[self.reactions.index(rxn)] 

        #--- Reset model properties ---
        model_cpds = sorted(tuple(set([c for r in model_rxns for c in r.stoichiometry.keys()], key = lambda x: x.id)))
        model_cpts = sorted(list(set([c.compartment for c in model_cpds if c.compartment != None])), key = lambda x: x.id)

        self.compartments_by_id = dict([(cpt.id,cpt) for cpt in self.compartments])    
        self.compounds_by_id = dict([(cpd.id,cpd) for cpd in self.compounds])    
        self.compounds_by_clean_id = dict([(remove_non_alphanumeric(cpd.id).lower(),cpd) for cpd in self.compounds])    
        self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    
        self.reactions_by_clean_id = dict([(remove_non_alphanumeric(rxn.id).lower(),rxn) for rxn in self.reactions])    


    def validate(self, reassign_props = True):
        """
        Checks the possible errors in the model and fixes them 
        It returns a warning in the output for issues that cannot be fixed.

        INPUTS:
        -------
        reassign_props: Must be True or False showng whether to ru self.assign_props()

        This function performs the following tasks:
        Fixs the following issues:
        ----------------------
          - Removes any duplicates in the list of compounds
          - Removes any duplicates in the list of reactions    
          - Removes any duplicates in fields reactions, reactant_reactions and product_reactions
            for each compound in the model
          - Removes any duplicates in fields reactants, products and compound for each reaction
            in the model
          - Checks if any compound participates in at least one reaction as a 
            reactant or product and fixes it otherwise. A warningn is issued, in case 
            this issue cannot be fixed 
          - Checks whether a compound appear in a reaction with a zero stoichiometry and fixes it
          - Checks if any reaction contains at least one participating compound and
            fixes it if otherwise. A warning is issued, if it cannot be fixed 
          - Checks and fixes any compound object that is referred to in reaction.compounds,
            but is not in the list of compounds in the model
          - Checks and fixes any reaction object that is referred to in compound.reactions
            but is not present in the list of reactions in the model.
          - Checks and fixes any compound that is referred to in a reaction.stoichiometry,
            reaction.compounds, or reaction.reactants but is not in the lsit of compound.reactions,
            compound.reactant_reactions or compound.product_reactions
          - Checks and fixes any reaction that is referred to in a compound.reactions, 
            compound.reactant_reactions, or compound.product_reactions, but is not in reaction.compounds,
            reaction.reactants or reaction.products
          - Check and fixes any repetitions in the compound.reactions, compound.reactant_reactions and 
            compound.product_reactions. This is very imporant to check otherwise it may cause serious
            issues in flux balance analysis or similar methods. 

        Issues a warning in the output for the folloiwng issues:
        ----------------------
          - Any replicates in compound ids
          - Any replicates in reaction ids
 
        """
        if self.stdout_msgs:
            print 'Validating the model ...'

        if not isinstance(reassign_props,bool):
            raise TypeError("'reassign_props' must be either True or False")

        if reassign_props:
            self.assign_props()

        errors_in_model = False

        #-- Check problems with compounds ---
        # Check for duplicates in the list of compound objects and fix them
        if len(set(self.compounds)) < len(self.compounds):
            errors_in_model = True
            self.compounds = sorted(tuple(set(self.compounds)),key=lambda x:x.id)
            if self.stdout_msgs:
                print 'Duplicates in the list of compounds were fixed'

        # Removes any duplicates in fields reactions, reactant_reactions and product_reactions
        # for each compound in the model
        for cpd in self.compounds:
            cpd.reactions = list(set(cpd.reactions))            
            cpd.reactant_reactions = list(set(cpd.reactant_reactions))            
            cpd.product_reactions = list(set(cpd.product_reactions))            
   
        # Check if there are any duplicates in compound ids and return an error if any
        cpd_ids = [m.id for m in self.compounds]
        if len(set(cpd_ids)) < len(cpd_ids):
            errors_in_model = True

            # Count how many times each compound id has been repated
            cpd_id_counts = dict((id,cpd_ids.count(id)) for id in  list(set(cpd_ids)))

            # compound ids repated more than once
            repeated_cpd_ids = [(id,cpd_id_counts[id]) for id in cpd_id_counts.keys() if cpd_id_counts[id] > 1]
            raise userError('The following compound ids are repeated: ' + str(repeated_cpd_ids)) 
            
        # Check if each compound participates in at least a reaction in the model
        cpds_rxns_dict = dict([(c.id,len(c.reactions)) for c in self.compounds])
        rxns_cpds_dict = dict([(l, []) for l in range(max(cpds_rxns_dict.values()) + 1)])
        for cpd in self.compounds:
            rxns_cpds_dict[cpds_rxns_dict[cpd.id]].append(cpd)
        if len(rxns_cpds_dict[0]) > 0:
            errors_in_model = True
            for cpd in rxns_cpds_dict[0]:
                cpd.reactions = [r for r in self.reactions if cpd in r.stoichiometry.keys()]
                cpd.reactant_reactions = [r for r in cpd.reactions if r.stoichiometry[cpd] < 0]
                cpd.product_reactions = [r for r in cpd.reactions if r.stoichiometry[cpd] > 0]

            no_rxns_cpds_afterFix = [c for c in rxns_cpds_dict[0] if len(c.reactions) == 0]
            if len(rxns_cpds_dict[0]) > len(no_rxns_cpds_afterFix) and self.stadout_msgs:
                print '\t{} compounds not participated in any reactions in the model were fixed'.format(len(rxns_cpds_dict[0]) - len(no_rxns_cpds_afterFix))
            if len(no_rxns_cpds_afterFix) > 0:
                if len(no_rxns_cpds_afterFix) > 100:
                    print '**WARNING (model.py): {} compounds do not participate in any reactions. These compounds include: {} and {} more'.format(len(rxns_cpds_dict[0]), [c.id for c in no_rxns_cpds_afterFix[:100]], len(no_rxns_cpds_afterFix) - 100)
                else:
                    print '**WARNING (model.py): {} compounds do not participate in any reactions. These compounds include: {}'.format(len(rxns_cpds_dict[0]), [c.id for c in no_rxns_cpds_afterFix])

        # compounds that are being referred to in a reaction object but that are not
        # present in the list of compounds in the model
        missing_cpds = list(set([m for r in self.reactions for m in r.stoichiometry.keys()]) - set(self.compounds)) 
        if len(missing_cpds) > 0:
            errors_in_model = True
            self.compounds += missing_cpds
            self.compounds = sorted(tuple(set(self.compounds)),key=lambda x:x.id)
            self.compounds_by_id = dict([(cpd.id,cpd) for cpd in self.compounds])    
            if self.stdout_msgs:
                print len(missing_cpds),' compounds were added to the list of compounds'

        # Checks and fixes any compounds that appear in a reaction stoichiomtery with a stoichiometric coefficient of zero
        cpd_counter = 0
        if 0 in list(set([s for r in self.reactions for s in r.stoichiometry.values()])):

            for cpd in [c for c in self.compounds if len([r for r in c.reactions if r.stoichiometry[c] == 0]) > 0]:
                cpd_counter += 1
                for rxn in [r for r in cpd.reactions if r.stoichiometry[cpd] == 0]:
                    # Remove this reaction from cpd.reactions, cpd.reactant_reactions or cpd.product_reactions
                    del rxn.stoichiometry[cpd]

                    try: 
                        del cpd.reactions[cpd.reactions.index(rxn)] 
                    except:
                        pass
                    try: 
                        del cpd.reactant_reactions[cpd.reactant_reactions.index(rxn)]       
                    except:
                        pass
                    try: 
                        del cpd.product_reactions[cpd.product_reactions.index(rxn)]       
                    except:
                        pass

                    # Remove this compound from rxn.compounds, rxn.reactants and rxn.products
                    try: 
                        del rxn.compounds[rxn.compounds.index(cpd)]  
                    except:
                        pass
                    try: 
                        del rxn.reactants[rxn.reactants.index(cpd)]              
                    except:
                        pass
                    try: 
                        del rxn.products[rxn.products.index(cpd)]              
                    except:
                        pass

        if cpd_counter > 0:
            errors_in_model = True
            if self.stdout_msgs:
                print '\t{} compounds appearing in reactions with zero stoichiometry were fixed ...'.format(cpd_counter)

        #--- Check problems with reactions ----
        self.reactions = sorted(tuple(set(self.reactions)),key=lambda x:x.id)
        for rxn in self.reactions:
            rxn.compounds = list(set(rxn.compounds))
            rxn.reactants = list(set(rxn.reactants))
            rxn.products = list(set(rxn.products))

        # Check if there are any replicates in reaction ids and return an error if any 
        rxn_ids = [r.id for r in self.reactions]
        if len(set(rxn_ids)) < len(rxn_ids):
            errors_in_model = True

            # Count how many times each reaction id has been repated
            rxn_id_counts = dict((id,rxn_ids.count(id)) for id in  rxn_ids)

            # reaction ids repated more than once
            repeated_rxn_ids = [(id,rxn_id_counts[id]) for id in rxn_id_counts.keys() if rxn_id_counts[id] > 1]
            raise userError('The following reaction ids are repeated: ' + str(repeated_rxn_ids)) 

        # Check if there are any reactions with no defined stoichiometry 
        if 0 in [len(r.stoichiometry) for r in self.reactions] and self.warnings:
            errors_in_model = True
            no_stoic_rxns = [r.id for r in self.reactions if len(r.stoichiometry) == 0]
            print "WARNING! 'stoichiometry' is not defined for these reactions: ",no_stoic_rxns,'\n'

        # reactions that are being referred to in a compound object but that are not
        # present in the list of reactions in the model
        missing_rxns = list(set([r for c in self.compounds for r in c.reactions]) - set(self.reactions)) 
        if len(missing_rxns) > 0:
            errors_in_model = True
            self.reactions += missing_rxns
            self.reactions = sorted(set(list(self.reactions)),key=lambda x:x.id)
            self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    
            if self.stdout_msgs:
                print len(missing_rxns),' missing reactions were added to the list of reactions'

        # Checks and fixes any reactions that is referred to in a compound.reactions, 
        # compound.reactant_reactions, or compound.product_reactions, but is not in reaction.compounds,
        # reaction.reactants or reaction.products
        rxn_counter = 0
        rxns_cpds = dict((rxn,[]) for rxn in self.reactions)
        rxns_reactants = dict((rxn,[]) for rxn in self.reactions)
        rxns_products = dict((rxn,[]) for rxn in self.reactions)
        for cpd in self.compounds:
            counter = 0
            for rxn in cpd.reactions:
                rxns_cpds[rxn].append(cpd)
            for rxn in cpd.reactant_reactions:
                rxns_reactants[rxn].append(cpd)
            for rxn in cpd.product_reactions:
                rxns_products[rxn].append(cpd)

        for rxn in self.reactions:
            if set(rxns_cpds[rxn]) != set(rxn.compounds):
                raise userError('The set of compounds in the model in which rxn ' + rxn.id + ' appears in their "reactions" field does not match the set of compounds in the stoichiometry of this reaction')

            if set(rxns_reactants[rxn]) != set(rxn.reactants):
                raise userError('The set of compounds in the model in which rxn ' + rxn.id + ' appears in their "reactant_reactions" field does not match the set of reactants of this reaction based on its stoichiometry')

            if set(rxns_products[rxn]) != set(rxn.products):
                raise userError('The set of compounds in the model in which rxn ' + rxn.id + ' appears in their "product_reactions" field (i.e.,' + str(set([c.id for c in self.compounds if rxn in c.product_reactions])) + ')does not match the set of products of this reaction based on its stoichiometry (i.e., ' + str(set([c.id for c in rxn.stoichiometry.keys() if rxn.stoichiometry[c] > 0])) + ')')

        if rxn_counter > 0:
            errors_in_model = True
            if self.stdout_msgs:
                print '\t{} reactions with with issues in their "compounds", "reactant" or "products" were fixed ...'.format(rxn_counter)

        if errors_in_model == False and self.stdout_msgs:
            print '\tChecked the model. All are OK'
 
    def print_reactions(self,ref_type = 'id', print_equation = True, cmp_ref = 'id'):
        """
        Prints in the output the list of all reactions and in the model 
        and their equations with the format specified by ref_type 
 
        INPUTS:
        -------
              ref_type: A string indicating the in what format the reactions should be 
                        printed in the output. Current eligible choices are 'id' or 'name'.
                        Default is 'id'.
        print_equation: If True reaction equation is also printed wherein compounds appear 
                        with the format specified by cmp_ref  
             cmp_ref: A string indicating in what format the compounds should appear
                        in a reaction equation if ref_type = 'equation'. If None, compound
                        ids are used
        """
        if ref_type.lower() not in ['id','name','equation']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")

        for reaction in self.reactions:
            if ref_type.lower() == 'id':
                rxn = reaction.id
            elif ref_type.lower() == 'name':
                if reaction.name is None:
                    rxn = reaction.id
                else:
                    rxn = reaction.name

            if print_equation == True:
                rxn = rxn + ':   ' + reaction.get_equation(show_cpds_by = cmp_ref)

            print rxn
        return ''

    def reset_flux_bounds(self):
        """
        Restores default flux bounds, i.e., flux bounds based on reaction types
        """
        for r in self.reactions:
            r.assign_flux_bounds()

    def print_compounds(self,ref_type = 'id'):
        """
        Prints in the output the list of all compounds in the model 
        with the format specified by ref_type 
 
        INPUTS:
        -------
        ref_type: A string indicating the in what format the compounds should be 
                  printed in the output. Current eligible choices are 'id', 'name'
                  and formula. If name or id was not provided, id used instead.
        """
        if ref_type.lower() not in ['id','name','formula']:
            raise userError("**Error! Invalid reference type (eligible choices are 'id' and 'name')")

        for cmp in self.compounds:
            if ref_type.lower() == 'id':
                mm = cmp.id
            elif ref_type.lower() == 'name':
                if cmp.name is None:
                    mm = cmp.id
                else:
                    mm = cmp.name
            elif ref_type.lower() == 'formula':
                if cmp.formula is None:
                    mm = cmp.id
                else:
                    mm = cmp.formula

            print mm
        return ''

    def convert_to_cobra(self):
        """
        Converts the model to COBRApy format 
        """
        import cobra

        # Create an empty COBRA model
        cobra_model = cobra.Model(self.id)

        # Compartments
        cobra_model.compartments = dict([(compart.id,compart.name) for compart in self.compartments])

        # Create compound objects
        for cmp in self.compounds:
            cobra_cmp = cobra.Metabolite(id = cmp.id, name = cmp.name, formula = cmp.formula)
            cobra_cmp.notes = cmp.notes
            if cmp.compartment != None:
                cobra_cmp.compartment = cmp.compartment.id
            cobra_model.add_metabolites(cobra_cmp)

        # Create reaction objects and add them to the model
        for rxn in self.reactions:
            cobra_rxn = cobra.Reaction(rxn.id)
            cobra_rxn.name = rxn.name
            cobra_rxn.lower_bound = rxn.flux_bounds[0]
            # Cobra does not allow to set the reaction reversibility
            if rxn.reversibility.lower() == 'exchange': 
                cobra_rxn.lower_bound = -1000
            else:
                cobra_rxn.subsyste = rxn.subsystem
            cobra_rxn.upper_bound = rxn.flux_bounds[1]

            if rxn.objective_coefficient != None:
                cobra_rxn.objective_coefficient = rxn.objective_coefficient
            stoic = {}
            for m in rxn.stoichiometry.keys():
                stoic[cobra_model.metabolites.get_by_id(m.id)] = rxn.stoichiometry[m]
            cobra_rxn.add_metabolites(stoic)
            if isinstance(rxn.gene_reaction_rule,str):
                cobra_rxn.gene_reaction_rule = rxn.gene_reaction_rule
            else:
                if rxn.genes != None:
                    cobra_rxn.gene_reaction_rule = ' or '.join([g.id for g in rxn.genes])
            cobra_model.add_reaction(cobra_rxn)

        return cobra_model

    def export(self, output_format = 'sbml', output_filename = ''):
        """
        Exports the model to various output formats.

        INPUTS:
        ------
        output_format: 
        Format of the output file. Eligible choices are:
          'sbml': SBML
        'pydict': Python dictionaries in .py files --> To be added
          'json': Python's json format --> To be added
        'pickle': Python's pickle format (.pk) --> To be added 

        output_filename: 
        Name of the output file
        """
        if not isinstance(output_format,str):
            raise TypeError('output_format must be a string.')
        elif output_format.lower() not in ['sbml','pydict','json','pickle']:
            raise ValueError('Invalid value for output_format (allowed choices are sbml, pydict, json and pickle)')

        if output_filename == '':
            output_filename = self.id

        #------- COBRA ------
        if output_format.lower() == 'sbml':
            cobra_model = self.convert_to_cobra(self)
            cobra.io.write_sbml_model(cobra_model,output_filename + '.xml',use_fbc_package=False)       
            if stdout_msgs:
                print 'The model was exported to {}'.format(output_filename + '.xml')

        #--------- pydict -----------
        elif output_format.lower() == 'pydict':
            import inspect

            model = {}
            
            #-- organism --
            model['organism'] = {}
            model['organism_global_attrs'] = {}

            if self.organism != None:
                # Attributes that are class methods
                method_attrs = [a for a in dir(self.organism) if inspect.ismethod(self.organism.__getattribute__(a))]
    
                # known attributes that are an instance of another class
                known_classType_attrs = ['model']
    
                # Attributes that are an instance of an unknown  class
                unknown_classType_attrs = [a for a in dir(self.organism) if '__' not in a and ('class' in str(type(self.organism.__getattribute__(a))) or 'object at' in repr(self.organism.__getattribute__(a))) and a not in known_classType_attrs + method_attrs]
                if len(unknown_classType_attrs) > 0 and self.warnings:
                    print '**WARNING (model.export_model)! The following unknonw attributes for "organism" were not exported: {}'.format(unknown_classType_attrs)
    
                # Global class attributes that are not an instance of a class
                global_attrs = [a for a in dir(self.organism) if '__' not in a  and a not in self.organism.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
    
                # Instance variables (attributes) that are not an instance of a class
                instance_attrs = [a for a in dir(self.organism) if a in self.organism.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
        
                for attr in global_attrs:
                    model['organism_global_attrs'][attr] = self.organism.__getattribute__(attr) 
                for attr in instance_attrs:
                    model['organism'][attr] = self.organism.__getattribute__(attr) 
                model['organism']['model_id'] = self.organism.model_id  # just save the id

            #-- Compartments --
            model['compartments'] = {}
            model['compartments_global_attrs'] = {}

            if self.compartments != ():
                for cpt in self.compartments:
                    attr_dict = {}
    
                    # Attributes that are class methods
                    method_attrs = [a for a in dir(cpt) if inspect.ismethod(cpt.__getattribute__(a))]
    
                    # known attributes that are an instance of another class
                    known_classType_attrs = ['model']
        
                    # Attributes that are an instance of an unknown class
                    unknown_classType_attrs = [a for a in dir(cpt) if '__' not in a and ('class' in str(type(cpt.__getattribute__(a))) or 'object at' in repr(cpt.__getattribute__(a))) and a not in known_classType_attrs + method_attrs]
                    if len(unknown_classType_attrs) > 0 and self.warnings:
                        print '**WARNING (model.export_model)! The following unknonw attributes for "compartment" were not exported: {}'.format(unknown_classType_attrs)
        
                    # Global class attributes that are not an instance of a class
                    global_attrs = [a for a in dir(cpt) if '__' not in a  and a not in cpt.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
        
                    # Instance variables (attributes) that are not an instance of a class
                    instance_attrs = [a for a in dir(cpt) if a in cpt.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
                    for attr in instance_attrs:
                        attr_dict[attr] = cpt.__getattribute__(attr) 
                    attr_dict['model_id'] = cpt.model_id
                    if cpt.id in model['compartments'].keys():
                       raise userError('duplicated compartment id: {}'.format(cpt.id))
                    else:
                        model['compartments'][cpt.id] = attr_dict
                # Global attributes are the same for all compartments, so use the one from the last cpt
                model['compartments']['global_attrs'] = {}
                for attr in global_attrs:
                    model['compartments_global_attrs'][attr] = cpt.__getattribute__(attr) 
    
            #-- Compounds --
            model['compounds'] = {} 
            model['compounds_global_attrs'] = {} 

            if self.compounds != ():
                for cpd in self.compounds:
                    attr_dict = {}
    
                    # Attributes that are class methods
                    method_attrs = [a for a in dir(cpd) if inspect.ismethod(cpd.__getattribute__(a))]
    
                    # known attributes that are an instance of another class
                    known_classType_attrs = ['model','compartment','reactions','reactant_reactions','product_reactions']
        
                    # Attributes that are an instance of an unknown class
                    unknown_classType_attrs = [a for a in dir(cpd) if '__' not in a and ('class' in str(type(cpd.__getattribute__(a))) or 'object at' in repr(cpd.__getattribute__(a))) and a not in known_classType_attrs + method_attrs]
                    if len(unknown_classType_attrs) > 0 and self.warnings:
                        print '**WARNING (model.export_model)! The following unknonw attributes for "compound" were not exported: {}'.format(unknown_classType_attrs)
        
                    # Global class attributes that are not an instance of a class
                    global_attrs = [a for a in dir(cpd) if '__' not in a  and a not in cpd.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
        
                    # Instance variables (attributes) that are not an instance of a class
                    instance_attrs = [a for a in dir(cpd) if a in cpd.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
            
                    for attr in instance_attrs:
                        attr_dict[attr] = cpd.__getattribute__(attr) 
                    attr_dict['model_id'] = cpd.model_id
                    attr_dict['compartment'] = cpd.compartment.id
                    attr_dict['reactions'] = [r.id for r in cpd.reactions]
                    attr_dict['reactant_reactions'] = [r.id for r in cpd.reactant_reactions]
                    attr_dict['product_reactions'] = [r.id for r in cpd.product_reactions]
                    if cpd.id in model['compounds'].keys():
                        raise userError('Duplicated compound id: {}'.format(cpd.id))
                    else:
                        model['compounds'][cpd.id] = attr_dict
                model['compounds']['global_attrs'] = {}
                for attr in global_attrs:
                    model['compounds_global_attrs'][attr] = cpd.__getattribute__(attr) 
                   
            #-- Reactions --
            model['reactions'] = {} 
            model['reactions_global_attrs'] = {} 

            if self.reactions != ():            
                for rxn in self.reactions:
                    attr_dict = {}
    
                    # Attributes that are class methods
                    method_attrs = [a for a in dir(rxn) if inspect.ismethod(rxn.__getattribute__(a))]
    
                    # known attributes that are an instance of another class
                    known_classType_attrs = ['model','compounds','reactants','products','stoichiometry','compartment','genes']
        
                    # Attributes that are an instance of an unknown class
                    unknown_classType_attrs = [a for a in dir(rxn) if '__' not in a and ('class' in str(type(rxn.__getattribute__(a)) or 'object at' in repr(rxn.__getattribute__(a)))) and a not in known_classType_attrs]
                    if len(unknown_classType_attrs) > 0 and self.warnings:
                        print '**WARNING (model.export_model)! The following unknonw attributes for "reaction" were not exported: {}'.format(unknown_classType_attrs)
        
                    # Global class attributes that are not an instance of a class
                    global_attrs = [a for a in dir(rxn) if '__' not in a  and a not in rxn.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
        
                    # Instance variables (attributes) that are not an instance of a class
                    instance_attrs = [a for a in dir(rxn) if a in rxn.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
            
                    for attr in instance_attrs:
                        attr_dict[attr] = rxn.__getattribute__(attr) 
                    attr_dict['model_id'] = rxn.model_id
                    attr_dict['compartments'] = [c.id for c in rxn.compartments]
                    attr_dict['compounds'] = [c.id for c in rxn.compounds]
                    attr_dict['reactants'] = [c.id for c in rxn.reactants]
                    attr_dict['products'] = [c.id for c in rxn.products]
                    attr_dict['genes'] = [g.id for g in rxn.genes]
                    attr_dict['stoichiometry'] = dict([(c.id,rxn.stoichiometry[c]) for c in rxn.compounds])
                    if rxn.id in model['reactions'].keys():
                        raise userError('Duplicated reaction id: {}'.format(rxn.id))
                    else:
                        model['reactions'][rxn.id] = attr_dict
                model['reactions']['global_attrs'] = {}
                for attr in global_attrs:
                    model['reactions_global_attrs'][attr] = rxn.__getattribute__(attr) 
    
            #-- Genes --
            model['genes'] = {} 
            model['genes_global_attrs'] = {} 

            if self.genes != ():
                for gen in self.genes:
                    attr_dict = {}
    
                    # Attributes that are class methods
                    method_attrs = [a for a in dir(gen) if inspect.ismethod(gen.__getattribute__(a))]
    
                    # known attributes that are an instance of another class
                    known_classType_attrs = ['model','compartment','reactions']
        
                    # Attributes that are an instance of an unknown class
                    unknown_classType_attrs = [a for a in dir(gen) if '__' not in a and ('class' in str(type(gen.__getattribute__(a))) or 'object at' in repr(gen.__getattribute__(a))) and a not in known_classType_attrs + method_attrs]
                    if len(unknown_classType_attrs) > 0 and self.warnings:
                        print '**WARNING (model.export_model)! The following unknonw attributes for "gene" were not exported: {}'.format(unknown_classType_attrs)
        
                    # Global class attributes that are not an instance of a class
                    global_attrs = [a for a in dir(gen) if '__' not in a  and a not in gen.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
        
                    # Instance variables (attributes) that are not an instance of a class
                    instance_attrs = [a for a in dir(gen) if a in gen.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + method_attrs]
            
                    for attr in instance_attrs:
                        attr_dict[attr] = gen.__getattribute__(attr) 
                    attr_dict['model_id'] = gen.model.id
                    attr_dict['compartment'] = [c.id for c in gen.compartment]
                    attr_dict['reactions'] = [r.id for r in gen.reactions]
                    if gen.id in model['genes'].keys():
                        raise userError('Duplicated gene id: {}'.formaT(gen.id))
                    else:
                        model['genes'][gen.id] = attr_dict
                model['genes']['global_attrs'] = {}
                for attr in global_attrs:
                    model['genes_global_attrs'][attr] = gen.__getattribute__(attr) 
    
            #-- Other model properties --
            model['model_global_attrs'] = {}
            model['model_instance_attrs'] = {}

            # Attributes already taken care of (see agove) or others that should not be saved (as they are created during
            # the instantiation of class 'model'
            ignore_attrs = ['organism','compartments','compounds','reactions','genes','reactions_by_id','reactions_by_clean_id','compounds_by_id','compounds_by_clean_id', 'genes_by_id','compartments_by_id']

            # Attributes that are class methods
            method_attrs = [a for a in dir(self) if inspect.ismethod(self.__getattribute__(a))]

            # known attributes that are an instance of another class
            known_classType_attrs = ['biomass_reaction','atpm_reaction']
    
            # Attributes that are an instance of an unknown class
            unknown_classType_attrs = [a for a in dir(self) if '__' not in a and ('class' in str(type(self.__getattribute__(a))) or 'object at' in repr(self.__getattribute__(a))) and a not in known_classType_attrs + method_attrs + ignore_attrs]
            if len(unknown_classType_attrs) > 0 and self.warnings:
                print '**WARNING (model.export_model)! The following unknonw attributes for "model" were not exported: {}'.format(unknown_classType_attrs)
    
            # Global class attributes that are not an instance of a class
            global_attrs = [a for a in dir(self) if '__' not in a  and a not in self.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + ignore_attrs + method_attrs]
    
            # Instance variables (attributes) that are not an instance of a class
            instance_attrs = [a for a in dir(self) if a in self.__dict__.keys() and a not in known_classType_attrs + unknown_classType_attrs + ignore_attrs + method_attrs]
        
            for attr in global_attrs:
                model['model_global_attrs'][attr] = self.__getattribute__(attr) 

            for attr in instance_attrs:
                model['model_instance_attrs'][attr] = self.__getattribute__(attr) 
            if self.biomass_reaction != None:
                model['model_instance_attrs']['biomass_reaction'] = self.biomass_reaction.id
            else:
                model['model_instance_attrs']['biomass_reaction'] = None 
            if self.atpm_reaction != None:
                model['model_instance_attrs']['atpm_reaction'] = self.atpm_reaction.id
            else:
                model['model_instance_attrs']['atpm_reaction'] = None 

            #-- Save model to file --
            with open(output_filename,'w') as f:
                f.write('model = {\n')

                # organism
                f.write("\t'organism_global_attrs': {")
                for k in sorted(model['organism_global_attrs'].keys()):
                    f.write("'{}': {}, ".format(k, repr(model['organism_global_attrs'][k])))        
                f.write("},\n")

                f.write("\t'organism': {")
                for k in sorted(model['organism'].keys()):
                    f.write("'{}': {},".format(k,repr(model['organism'][k])))        
                f.write("},\n")

                # compartments, compounds, reactions and genes
                for model_attr in ['compartments','compounds','reactions','genes']:
                    f.write("\n\t'{}_global_attrs':{{".format(model_attr))
                    for k in sorted(model[model_attr + '_global_attrs'].keys()):
                        f.write("'{}: {},\n'".format(k,repr(model[model_attr + '_global_attrs'][k])))
                    f.write("},\n")

                    f.write("\n\t'{}': {{\n".format(model_attr))
                    for entry in sorted([e for e in model[model_attr].keys() if e != 'global_attrs']):
                        f.write("\t\t'{}':{{".format(entry))
                        for k in sorted(model[model_attr][entry].keys()):
                            f.write("'{}':{}, ".format(k,repr(model[model_attr][entry][k])))
                        f.write("},\n")
                    f.write("\t},\n")

                # Instance and global variables of the model
                f.write("\n\t'model_global_attrs': {\n")
                for k in sorted(model['model_global_attrs'].keys()):
                    f.write("\t\t'{}': {},\n".format(k,repr(model['model_global_attrs'][k])))        
                f.write("\t},\n")
                f.write("\n\t'model_instance_attrs': {\n")
                for k in sorted(model['model_instance_attrs'].keys()):
                    f.write("\t\t'{}': {},\n".format(k,repr(model['model_instance_attrs'][k])))        
                f.write("\t}\n")

                f.write('}\n')

        #----- json ----
        elif output_format.lower() == 'json':
            import json as js
            js.dump(self)

        #----- Pickle or cPickle ----
        elif output_format.lower() == 'pickle':
            import cPickle as pk
            pk.dump(self) 

    def fba(self,optimization_solver = default_optim_solver, build_new_optModel = True, maximize = True, reset_fluxes = True, store_opt_fluxes = True, store_all_rxn_fluxes = False, flux_key = None, run = True, assign_wildType_max_biomass = False, simulation_conditions = '', stdout_msgs = True, warnings = True):
        """
        Creates a fba model for this model and runs fba 

        INPUTS:
        -------- 
               store_all_rxn_fluxes: If True store_flux is set to True for all reactions in the model
                       reset_fluxes: If True all reaction fluxes are set to None before performign any FBA-based
                                     simulation. This parameter is useful when ine has set store_flux to True only
                                     for a subset of rxns in the model. So setting all reaction fluxes to None avoids
                                     any confusions if a reaction flux with store_flux = False has already been assigned
                                     from another FBA problem
                                run: Set to False to create a fba model 
                                     but do not run fba 
        assign_wildType_max_biomass: If True the optimal objective function
                                     value of the fba is assigned to global 
                                     variable wildType_max_biomass

        The rest of inputs are the same as those in fba
        """
        from tools.fba.fba import fba

        # Reaction with no "store_flux" variable
        if store_all_rxn_fluxes:
            for rxn in self.reactions:
                rxn.store_flux = True
        else:
            for rxn in [r for r in self.reactions if not hasattr(r,'store_flux')]:
                rxn.store_flux = False

        # Reset reaction fluxes
        if not isinstance(reset_fluxes,bool):
            raise TypeError("'reset_fluxes' must be either True or False")
        if reset_fluxes:
            for rxn in [r for r in self.reactions if not r.store_flux]:
                r.flux = None

        if build_new_optModel == True:
            self.fba_model = fba(model = self,optimization_solver = optimization_solver, build_new_optModel = build_new_optModel, flux_key = flux_key, store_opt_fluxes = store_opt_fluxes, simulation_conditions = simulation_conditions, maximize = maximize, warnings = warnings, stdout_msgs = stdout_msgs)
        else:
            self.fba_model.build_new_optModel = build_new_optModel
            self.fba_model.flux_key = flux_key
            self.fba_model.optimization_solver = optimization_solver
            self.fba_model.store_opt_fluxes = store_opt_fluxes
            self.fba_model.maximize = maximize
            self.fba_model.simulation_conditions = simulation_conditions
            self.fba_model.warnings = warnings 
            self.fba_model.stdout_msgs = stdout_msgs 
            # Update the flux bounds
            for j in self.reactions:
                self.fba_model.optModel.v[j.id].setlb(j.flux_bounds[0])
                self.fba_model.optModel.v[j.id].setub(j.flux_bounds[1])

            # Update the objective function
            self.fba_model.optModel.del_component('objectiveFunc')
            self.fba_model.optModel.objectiveFunc = Objective(rule=self.fba_model.objectiveFunc_rule, sense = maximize)

        if run == True:
            self.fba_model.run()
            if assign_wildType_max_biomass == True:
                if self.fba_model.solution['exit_flag'] == 'globallyOptimal':
                    self.wildType_max_biomass = self.fba_model.solution['objective_value']


