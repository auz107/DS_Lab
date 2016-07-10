from __future__ import division
import sys
from collections import Counter
sys.path.append('../../')
from tools.userError import userError
from organism import organism
from compound import compound
from gene import gene
from reaction import reaction
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
from coopr.pyomo import *

class model(object):
    """
    This is a general class holding information for a model

    METHODS:
    --------
           get_compounds: Returns a dictionary of the selected compound objects
           get_reactions: Returns a dictionary of the selected reaction objects
           add_compounds: Add new compound to the model
           add_reactions: Add new reactions to the model
        remove_compounds: Remove compound from the model
        remove_reactions: Add reactions to the model
                validate: Check for probable issues in the model and either fixes them
                          or issues a warning in the output 
     restore_flux_bounds: Restores the original bounds (bounds based on reaction types) for
                          all reactions in the model
         print_compounds: Prints in the output the list of all compounds in the model
                          with a desired format 
         print_reactions: Prints in the output the list of all reactions in the model with
                          a desired format 
        convert_to_cobra: Converts the model to COBRApy format 
            export_model: Export the model to various formats 
                     fba: Creates a fba model        


    Ali R. Zomorrodi, Segre Lab @ BU
    Last updated: 08-17-2015
    """

    def __init__(self, id, type, organism = None, reactions = [], compounds = [], genes = [], compartments = [], name = None, biomass_reaction = None, atpm_reaction = None, simulation_condition = None, fba_model = None, wildType_max_biomass = None, notes = None,stdout_msgs = True, warnings = True):

        # Warnings and messages in the standard output
        if not isinstance(stdout_msgs,bool): 
            raise TypeError("stdout_msgs must be True or False")
        else:
            self.stdout_msgs = stdout_msgs

        if not isinstance(warnings,bool): 
            raise TypeError("warnings must be True or False")
        else:
            self.warnings = warnings

        # Name of the model (string) 
        self.id = id 

        # Type of the model (string). Examples include cmpolic, cmpolic_expressions, genetic
        self.type = type

        # An organism object 
        self.organism = organism

        # A list of reaction objects 
        if not isinstance(reactions,list):
            raise TypeError("'reactions' for a model must be a list")
        if reactions != []:
            for r in reactions:
                r.model = self
            self.reactions = sorted(reactions,key=lambda x:x.id)
        else:
            self.reactions = None

        # A list of compounds objects containing all compounds in the model 
        if not isinstance(compounds,list):
            raise TypeError("'compounds' for a model must be a list")
        if compounds == []:
            self._create_compounds_list()
        else:
            for m in compounds:
                m.model = self
            self.compounds = sorted(compounds,key=lambda x:x.id)

        # A list of gene objects containing all genes in the model
        if not isinstance(genes,list):
            raise TypeError("'genes' for a model must be a list")
        if genes == []:
            self._create_genes_list()
        else:
            for gn in genes:
                gn.model = self
            self.genes = genes

        # A list of compartment objects containing compounds compartments 
        # for all compound in the model
        if not isinstance(compartments,list):
            raise TypeError("'compartments' for a model must be a list")
        if compartments == []:
            self._create_compartments_list()
        else:
            for cm in compartments:
                cm.model = self
            self.compartments = compartments

        # Create the stoichiometric matrix
        self._create_stoic_matrix()

        # Model name
        if name == None:
            self.name = self.id
        else:
            self.name = anme

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

        # Simulaiton conditions (string) 
        self.simulation_condition = simulation_condition

        # fba model
        self.fba_model = fba_model

        # Maximum biomass flux for the wild-type
        self.wildType_max_biomass = wildType_max_biomass

        # Notes and comments
        self.notes = notes

        # A dictionary whose keys are reaciotn ids and values are reaction objects
        self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    
        self.reactions_by_clean_id = dict([(remove_non_alphanumeric(rxn.id).lower(),rxn) for rxn in self.reactions])    

        # A dictionary whose keys are compound ids and values are compound objects
        self.compounds_by_id = dict([(cmp.id,cmp) for cmp in self.compounds])    
        self.compounds_by_clean_id = dict([(remove_non_alphanumeric(cmp.id).lower(),cmp) for cmp in self.compounds])    

        # A dictionary whose keys are gene ids and values are gene objects
        if self.genes != None:
            self.genes_by_id = dict([(gene.id,gene) for gene in self.genes])    

        # Check for probable issues in the model
        self.validate()
           
    def _create_compounds_list(self):
        """
        Creates a list of compound objects for all compounds present in the model 
        """
        if self.reactions != None and len([r for r in self.reactions if r.compounds != None]):
            self.compounds = []

        for reaction in self.reactions:
            self.compounds += reaction.compounds

        self.compounds = sorted(list(set(self.compounds)),key=lambda x:x.id) 
        for m in self.compounds:
            m.model = self

    def _create_genes_list(self):
        """
        Creates a list of compound objects for all compounds present in the model 
        """
        if len([r for r in self.reactions if r.genes != None]) == 0:
            self.genes = None 
        else:
            for reaction in [r for r in self.reactions if r.genes != None]:
                self.genes += reaction.genes

            self.genes = list(set(self.genes)) 
            for gn in self.genes:
                gn.model = self

    def _create_compartments_list(self):
        """
        Creates a list of compartment objects for compartments of all 
        compounds in the model 
        """
        if len([m for m in self.compounds if m.compartment != None]) == 0:
            self.compartments = None 
        else:
            for compound in [m for m in self.compounds if m.compartment != None]:
                self.compartments += compound.compartment

            self.compartments = list(set(self.compartments)) 
            for cm in self.compartments:
                cm.model = self 

    def _create_stoic_matrix(self):
        """
        Creates the stoichiometric matrix of the network. This is in fact 
        not a matrix. It adds to each compound object the list of reaction
        objects in which that compounds participates as a reactnat or product
        The stoichiometric coefficients can then be found by referring to the 
        particuular reaction object
        """
        for compound in self.compounds:
            compound.reactions = [reaction for reaction in self.reactions if compound in reaction.compounds]
            compound.reactant_reactions = [reaction for reaction in self.reactions if compound in reaction.reactants]
            compound.product_reactions = [reaction for reaction in self.reactions if compound in reaction.products]

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
        new_compounds: A dictionary of new compound objects with
                         keys and values as follows                   
                           Keys: New compound objects
                         Values: Can be an empty dictionary '{}' if a new compound participates 
                                 in a new reaction that does not exist in the model yet, or if
                                 it has already been assigned to reactions using add_compounds
                                 function of reaction object.
                                 Note that new reactions are added to the model after
                                 creating new compound objects. 
                                 If this new compound participates in some existing 
                                 reactions in the model, then Values can be another 
                                 dictionary with keys and values as follows: 
                                   keys: Reaction objects, showing existing reactions in 
                                         the model in which the compound participates 
                                         as a reactant or product
                                 Values: The stoichiometric coefficient of the compound
                                         in existing reactions in the model
        
        OUTPUTS:
        --------
        There is no actual output. The model is updated with the addition of
        new compounds 
        """
        self.compounds += new_compounds.keys()

        # Sort according to id
        self.compounds = sorted(list(set(self.compounds)),key=lambda x:x.id)

        # Update compounds by id
        for new_cmp in  new_compounds.keys():   
            self.compounds_by_id[new_cmp.id] = new_cmp
            self.compounds_by_clean_id[remove_non_alphanumeric(new_cmp.id).lower()] = new_cmp

        # Update the stoichiometry of existing reactions in the model 
        for cmp in [c for c in new_compounds.keys() if new_compounds[c] != {}]:
            for rxn in new_compounds[cmp].keys():
                rxn.add_compounds({cmp:new_compounds[cmp][rxn]})

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
        self.reactions += new_reactions

        # Sort according to id
        self.reactions = sorted(list(set(self.reactions)),key=lambda x:x.id)

        # Update reactions_by_id
        for new_rxn in new_reactions:
            self.reactions_by_id[new_rxn.id] = new_rxn 
            self.reactions_by_clean_id[remove_non_alphanumeric(new_rxn.id).lower()] = new_rxn 

        # Check if there are compounds in the reactions not in the model
        cmp_notInMoDel = [c for r in new_reactions for c in r.compounds if c not in self.compounds]
        if len(cmp_notInMoDel) > 0:
            raise ValueError('The following compounds not in the model ' + self.id + ': ' + str([c.id for c in cmp_notInMoDel]))

        # Update compound objects participating in these new reactions
        for rxn in new_reactions:
            for cmp in rxn.compounds:
                if rxn not in cmp.reactions:
                    cmp.reactions.append(rxn)
                if rxn.stoichiometry[cmp] < 0 and rxn not in cmp.reactant_reactions:
                    cmp.reactant_reactions.append(rxn) 
                if rxn.stoichiometry[cmp] > 0 and rxn not in cmp.product_reactions:
                    cmp.product_reactions.append(rxn) 

    def remove_compounds(self,removed_compounds):
        """
        Remove selected compounds from the model
        
        INPUTS:
        ------
        removed_compounds: A list of compound objectc that must be removed 
        """
        if type(removed_compounds) is not list:
            userError("**Error! the input to 'remove_compounds' should be a list of compound objects")

        # First remove this compound from all relevant reaction objects
        for cmp in removed_compounds:
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

    def remove_reactions(self,removed_reactions):
        """
        Remove selected reactions from the model
        
        INPUTS:
        ------
        removed_reactions: A list of reaction objectc that must be removed 
        """
        if type(removed_reactions) is not list:
            userError("**Error! the input to 'remove_reactions' should be a list of reaction objects")

        # First remove this reaction from all relevant reaction objects
        for rxn in removed_reactions:
            for cmp in rxn.compounds:
                del cmp.reactions[cmp.reactions.index(rxn)]
                if rxn.stoichiometry[cmp] < 0 and cmp.reactant_reactions != []:
                    del cmp.reactant_reactions[cmp.reactant_reactions.index(rxn)]
                elif rxn.stoichiometry[cmp] > 0 and cmp.product_reactions != []:
                    del cmp.product_reactions[cmp.product_reactions.index(rxn)]

            # Check the presence in all other compounds
            rxn_cmps = list(set([m for m in self.compounds if rxn in m.reactions + m.reactant_reactions + m.product_reactions])) 
            if len(rxn_cmps) > 0:
                for m in rxn_cmps:
                    try: 
                        del m.reactions[m.reactions.index(rxn)]
                    except:
                        continue            
                    try: 
                        del m.reactant_reactions[m.reactant_reactions.index(rxn)]
                    except:
                        continue            
                    try: 
                        del m.product_reactions[m.product_reactions.index(rxn)]
                    except:
                        continue            

            # Remove it from the list of reactions
            del self.reactions[self.reactions.index(rxn)] 
            self.reactions = sorted(set(list(self.reactions)),key=lambda x:x.id)

            # Remove it from reactions_by_id
            del self.reactions_by_id[rxn.id]
            del self.reactions_by_clean_id[remove_non_alphanumeric(rxn.id).lower()]

    def validate(self):
        """
        Checks the possible errors in the model and fixes them 
        It returns a warning in the output for issues that cannot be fixed.

        This function performs the following tasks:
        Fixs the following issues:
        ----------------------
          - Removes any duplicates in the list of compounds
          - Removes any duplicates in the list of reactions    
          - Checks if any compound participates in at least one reaction as a 
            reactant or product and fixes it otherwise. A warningn is issued, in case 
            this issue cannot be fixed 
          - Checks if any reaction contains at least one participating compound and
            fixes it if otherwise. A warning is issued, if it cannot be fixed 
          - Checks and fixes any compound object that is referred to in reaction.compounds,
            but is not in the list of compounds in the model
          - Checks and fixes any reaction object that is referred to in compound.reactions
            but is not present in the list of reactions in the model.
          - Checks and fixes any compound that is referred to in a reaction.stoichiomtery,
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
        errors_in_model = False

        #-- Check problems with compounds ---
        # Check for duplicates in the list of compound objects and fix them
        if len(set(self.compounds)) < len(self.compounds):
            errors_in_model = True
            self.compounds = sorted(list(set(self.compounds)),key=lambda x:x.id)
            if self.warnings:
                print 'WARNINGS! Duplicates in the list of compounds were fixed'

        # Check if there are any duplicates in compound ids and issue a warning
        cmp_ids = [m.id for m in self.compounds]
        if len(set(cmp_ids)) < len(cmp_ids):
            errors_in_model = True

            # Count how many times each compound id has been repated
            cmp_id_counts = dict((id,cmp_ids.count(id)) for id in  cmp_ids)

            # compound ids repated more than once
            repeated_cmp_ids = [(id,cmp_id_counts[id]) for id in cmp_id_counts.keys() if cmp_id_counts[id] > 1]
            if self.warnings:
                print 'WARNING! The following compound ids are repeated: ',repeated_cmp_ids 
            
        # Check if each compound participates in at least a reaction in the model
        no_rxn_cmps = [m for m in self.compounds if len(m.reactions) == 0]
        fixed_cmp_num = 0
        if len(no_rxn_cmps) > 0:
            errors_in_model = True
            for cmp in no_rxn_cmps:
                # Try to fix 
                m_rxns = [r for r in self.reactions if cmp in r.stoichiometry.keys()]
                if len(m_rxns) == 0 and self.warnings: 
                    print "WARNING! compound '",cmp.id,"' does not participate in any reactions in the model"         
                else: # fix it
                    cmp.reactions = m_rxns
                    cmp.reactant_reactions = [r for r in self.reactions if r.stoichiometry[cmp] < 0]
                    cmp.product_reactions = [r for r in self.reactions if r.stoichiometry[cmp] > 0]
                    fixed_cmp_num += 1

        if fixed_cmp_num > 0 and self.stadout_msgs.lower() == 'on':
            print '\t',fixed_cmp_num,' compounds not participated in any reactions in the model were fixed'

        # compounds that are being referred to in a reaction object but that are not
        # present in the list of compounds in the model
        missing_cmps = list(set([m for r in self.reactions for m in r.compounds + r.reactants + r.products + r.stoichiometry.keys() if m not in self.compounds])) 
        if len(missing_cmps) > 0:
            errors_in_model = True
            self.compounds += missing_cmps
            self.compounds = sorted(list(set(self.compounds)),key=lambda x:x.id)
            self.compounds_by_id = dict([(cmp.id,cmp) for cmp in self.compounds])    
            if self.stdout_msgs:
                print len(missing_cmps),' compounds were added to the list of compounds'

        # Checks and fixes any compound that is referred to in a reaction.stoichiomtery,
        # reaction.compounds, or reaction.reactants but is not in the lsit of compound.reactions,
        # compound.reactant_reactions or compound.product_reactions
        cmp_counter = 0
        for cmp in self.compounds:
            # Rxns in which cmp appears as a reactant or product
            cmp_rxns = [r for r in self.reactions if cmp in r.compounds or cmp in r.reactants or cmp in r.products or cmp in r.stoichiometry.keys()]
            if set(cmp_rxns) != set(cmp.reactions):
                cmp_counter += 1
                for r in [r for r in cmp_rxns if r not in cmp.reactions]:
                    cmp.reactions.append(r)
                for r in [r for r in cmp_rxns if cmp in r.reactants or r.stoichiometry[cmp] < 0 and r not in cmp.reactant_reactions]:
                    cmp.reactant_reactions.append(r)
                for r in [r for r in cmp_rxns if cmp in r.products or r.stoichiometry[cmp] > 0 and r not in cmp.product_reactions]:
                    cmp.product_reactions.append(r)
            
        if cmp_counter > 0:
            errors_in_model = True
            if self.stdout_msgs:
                print '\t{} compounds with missing reactions were fixed ...'.format(cmp_counter)

        # Check and fixes any repetitions in the compound.reactions, compound.reactant_reactions and 
        # compound.product_reactions. This is very imporant to check otherwise it may cause serious
        # issues in flux balance analysis or similar methods. 
        # To add        

        #-- Check problems with reactions ---
        # Check for duplicates in the list of reaction objects and fix them
        if len(set(self.reactions)) < len(self.reactions):
            self.reactions = sorted(list(set(self.reactions)),key=lambda x:x.id)
            if self.stdout_msgs:
                print '\tDuplicates in the list of reactions were fixed'

        # Check if there are ny replicates in reaction ids and issue a warning
        rxn_ids = [r.id for r in self.reactions]
        if len(set(rxn_ids)) < len(rxn_ids):
            # Count how many times each reaction id has been repated
            rxn_id_counts = dict((id,rxn_ids.count(id)) for id in  rxn_ids)

            # reaction ids repated more than once
            repeated_rxn_ids = [(id,rxn_id_counts[id]) for id in rxn_id_counts.keys() if rxn_id_counts[id] > 1]
            if self.warnings:
                print 'WARNING! The following reaction ids are repeated: ',repeated_rxn_ids 

        # Check if there are any reactions with no defined stoichiometry 
        no_stoic_rxns = [r.id for r in self.reactions if len(r.stoichiometry) == 0]
        if len(no_stoic_rxns) > 0 and self.warnings:
            print "WARNING! 'stoichiometry' is not defined for these reactions: ",no_stoic_rxns,'\n'

        # Check if each reaction has at least one participating compound
        no_cmp_rxns = [r for r in self.reactions if len(r.compounds) == 0]
        fixed_rxns_num = 0
        if len(no_cmp_rxns) > 0:
            errors_in_model = True
            for rxn in no_cmp_rxns:
                # Try to fix it 
                rxn.compounds = rxn.stoichiometry.keys()
                fixed_rxns_num += 1

        if fixed_rxns_num > 0 and self.stdout_msgs:
            print '\t',fixed_rxns_num," reactions with empty field 'compounds' were fixed"

        # reactions that are being referred to in a compound object but that are not
        # present in the lis of reactions in the model
        missing_rxns = list(set([r for m in self.compounds for r in m.reactions + m.reactant_reactions + m.product_reactions if r not in self.reactions])) 
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
        for rxn in self.reactions:
            # Compounds in which rxn appears in compounds reactions, reactant reactions or product reactions 
            rxn_cmps = [c for c in self.compounds if rxn in c.reactions or rxn in c.reactant_reactions or rxn in c.product_reactions]
            if set(rxn_cmps) != set(rxn.compounds):
                rxn_counter += 1
                for c in [c for c in rxn_cmps if c not in rxn.compounds]:
                    rxn.compounds.append(c)

                for c in [c for c in rxn_cmps if rxn in c.reactant_reactions and c not in rxn.reactants]:
                    rxn.reactants.append(c)
                    if c not in rxn.stoichiometry.keys():
                        raise userError('Unknown stoichiometric coefficient for compound ' + c.id + ' in rxn ' + rxn.id)
                    elif rxn.stoichiometry[c] > 0 and self.warnings:
                        print 'WARNING! Compound ' + c.id + ' has a stoichiometric coefficinet of ' + str(rxn.stoichiometry[c]) + ' in reaction ' + rxn.id + ' but reaction ' + rxn.id + ' appears in reactnat_reactions of ' + c.id

                for c in [c for c in rxn_cmps if rxn in c.product_reactions and c not in rxn.products]:
                    rxn.products.append(c)
                    if c not in rxn.stoichiometry.keys():
                        raise userError('Uknonwn stoichiometric coefficient for compound ' + c.id + ' in reaction ' + rxn.id)
                    elif rxn.stoichiometry[c] < 0 and self.warnings:
                        print 'WARNING! Compound ' + c.id + ' has a stoichiometric coefficinet of ' + str(rxn.stoichiometry[c]) + ' in reaction ' + rxn.id + ' but reaction ' + rxn.id + ' appears in product_reactions of ' + c.id
            
        if rxn_counter > 0:
            errors_in_model = True
            if self.stdout_msgs:
                print '\t{} reactions with missing compounds were fixed ...'.format(rxn_counter)

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
                rxn = rxn + ':   ' + reaction.get_equation(ref_type = cmp_ref)

            print rxn
        return ''

    def restore_flux_bounds(self):
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
            if rxn.type.lower() == 'exchange': 
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

    def export(self,output_format = 'sbml', output_filename = ''):
        """
        Export the model to various output formats

        INPUTS:
        ------
          output_format: Format of the output file. Eligible choices are:
                           'sbml': SBML
                         'pydict': Python dictionaries in .py files --> To be added
                           'json': Python's json format --> To be added
                         'pickle': Python's pickle format (.pk) --> To be added 
        output_filename: Name of the output file
        """
        if not isinstance(output_format,str):
            raise TypeError('output_format must be a string.')
        elif output_format.lower() not in ['sbml','pydict','json','pickle']:
            raise ValueError('Invalid value for output_format (allowed choices are sbml, pydict, json and pickle)')

        if output_filename == '':
            output_filename = self.id

        if output_format.lower() == 'sbml':
            cobra_model = convert_to_cobra(self)
            cobra.io.write_sbml_model(cobra_model,output_filename + '.xml',use_fbc_package=False)       
            if stdout_msgs:
                print 'The model was exported to {}'.format(output_filename + '.xml')
        elif output_format.lower() == 'pydict':
            pass
        elif output_format.lower() == 'json':
            pass
        elif output_format.lower() == 'pickle':
            pass
   
    def fba(self,optimization_solver = 'gurobi', build_new_optModel = True, store_opt_fluxes = True, flux_key = None, run = True, assign_wildType_max_biomass = False, stdout_msgs = True, warnings = True):
        """
        Creates a fba model for this model and runs fba 
        All inputs, except run, are the same as fba
                                run: Set to False to create a fba model 
                                     but do not run fba 
        assign_wildType_max_biomass: If True the optimal objective function
                                     value of the fba is assigned to global 
                                     variable wildType_max_biomass
        """
        from tools.fba.fba import fba

        if build_new_optModel == True:
            self.fba_model = fba(model = self,optimization_solver = optimization_solver, build_new_optModel = build_new_optModel, flux_key = flux_key, store_opt_fluxes = store_opt_fluxes, warnings = warnings, stdout_msgs = stdout_msgs)
        else:
            self.fba_model.build_new_optModel = build_new_optModel
            self.fba_model.flux_key = flux_key
            self.fba_model.optimization_solver = optimization_solver
            self.fba_model.store_opt_fluxes = store_opt_fluxes
            self.fba_model.warnings = warnings 
            self.fba_model.stdout_msgs = stdout_msgs 
            # Update the flux bounds
            for j in self.reactions:
                self.fba_model.optModel.v[j].setlb(j.flux_bounds[0])
                self.fba_model.optModel.v[j].setub(j.flux_bounds[1])

            # Update the objective function
            self.fba_model.optModel.del_component('objectiveFunc')
            self.fba_model.optModel.objectiveFunc = Objective(rule=self.fba_model.objectiveFunc_rule, sense = maximize)

        if run == True:
            self.fba_model.run()
            if assign_wildType_max_biomass == True:
                if self.fba_model.solution['exit_flag'] == 'globallyOptimal':
                    self.wildType_max_biomass = self.fba_model.solution['objective_value']


