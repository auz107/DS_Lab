from __future__ import division
import sys
from collections import Counter
sys.path.append('../../')
from coopr.pyomo import *
from coopr.opt import *
from tools.userError import userError
from organism import organism
from metabolite import metabolite
from gene import gene
from reaction import reaction
from tools.fba.fba import fba
import cobra

class model(object):
    """
    This is a general class holding information for a model

    METHODS:
    --------
         get_metabolites: Returns a dictionary of the selected metabolite objects
           get_reactions: Returns a dictionary of the selected reaction objects
         add_metabolites: Add new metabolites to the model
           add_reactions: Add new reactions to the model
      remove_metabolites: Remove metabolites from the model
        remove_reactions: Add reactions to the model
                validate: Check for probable issues in the model and either fixes them
                          or issues a warning in the output 
       print_metabolites: Prints in the output the list of all metabolites in the model
                          with a desired format 
         print_reactions: Prints in the output the list of all reactions in the model with
                          a desired format 
         export_to_cobra: Export the model to a COBRApy model
                     fba: Creates a fba model        


    Ali R. Zomorrodi, Segre Lab @ BU
    Last updated: 03-09-2015
    """

    def __init__(self, id, type, organism, reactions = [], metabolites = [], genes = [], compartments = [], name = '', biomass_reaction = None, atpm_reaction = None, simulation_condition = '', fba_model = None, wildType_max_biomass = None, notes = ''):

        # Name of the model (string) 
        self.id = id 

        # Type of the model (string). Examples include metabolic, metabolic_expressions, genetic
        self.type = type

        # An organism object 
        self.organism = organism

        # A list of reaction objects 
        for r in reactions:
            r.model = self
        self.reactions = sorted(reactions,key=lambda x:x.id)

        # A list of metabolite objects containing all metabolites in the model 
        if metabolites == []:
            self._createMetabolitesList()
        else:
            for m in metabolites:
                m.model = self
            self.metabolites = sorted(metabolites,key=lambda x:x.id)

        # A list of gene objects containing all genes in the model
        if genes == []:
            self._createGenesList()
        else:
            for gn in genes:
                gn.model = self
            self.genes = genes

        # A list of compartment objects containing metabolite compartments 
        # for all metabolites in the model
        if compartments == []:
            self._createCompartmentsList()
        else:
            for cm in compartments:
                cm.model = self
            self.compartments = compartments

        # Create the stoichiometric matrix
        self._createStoicMatrix()

        # Model name
        if name == '':
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
        if isinstance(notes,str):
            self.notes = notes
        else:
            self.notes = ''

        # A dictionary whose keys are reaciotn ids and values are reaction objects
        self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    

        # A dictionary whose keys are metabolite ids and values are metabolite objects
        self.metabolites_by_id = dict([(metab.id,metab) for metab in self.metabolites])    

        # A dictionary whose keys are gene ids and values are gene objects
        self.genes_by_id = dict([(gene.id,gene) for gene in self.genes])    

        # Check for probable issues in the model
        self.validate()
           
    def _createMetabolitesList(self):
        """
        Creates a list of metabolite objects for all metabolites present in the model 
        """
        self.metabolites = []

        for reaction in self.reactions:
            self.metabolites += reaction.metabolites

        self.metabolites = sorted(list(set(self.metabolites)),key=lambda x:x.id) 
        for m in self.metabolites:
            m.model = self

    def _createGenesList(self):
        """
        Creates a list of metabolite objects for all metabolites present in the model 
        """
        self.genes = []

        for reaction in [r for r in self.reactions if r.genes != []]:
            self.genes += reaction.genes

        self.genes = list(set(self.genes)) 
        for gn in self.genes:
            gn.model = self

    def _createCompartmentsList(self):
        """
        Creates a list of compartment objects for compartments of all 
        metabolites in the model 
        """
        self.compartments = []

        for metabolite in [m for m in self.metabolites if m.compartment != None]:
            self.compartments += metabolite.compartment

        self.compartments = list(set(self.compartments)) 
        for cm in self.compartments:
            cm.model = self 

    def _createStoicMatrix(self):
        """
        Creates the stoichiometric matrix of the network. This is in fact 
        not a matrix. It adds to each metabolite object the list of reaction
        objects in which that metabolites participates as a reactnat or product
        The stoichiometric coefficients can then be found by referring to the 
        particuular reaction object
        """
        for metabolite in self.metabolites:
            metabolite.reactions = [reaction for reaction in self.reactions if metabolite in reaction.metabolites]
            metabolite.reactant_reactions = [reaction for reaction in self.reactions if metabolite in reaction.reactants]
            metabolite.product_reactions = [reaction for reaction in self.reactions if metabolite in reaction.products]

    def get_reactions(self,reactions_ref):
        """
        Returns a dictionary of the selected reaction objects

        INPUTS:
        ------- 
        reactions_ref: A dictionary with keys and values as follows:
                         keys: A reaction id or name (string) 
                       values: A string indicating the format by which the reactions
                               should be searched. Eligible choices include 'id' and 'name'
                      Example: {'rxn1_id':'id','rxn2_name':'name'}
        
        OUTPUTS:
        -------
        If the input is a dictionary with a single key the the ouput is the corresponding
        reaction object, otherwise the output is a dictionary whose keys are the reaction 
        reference values (id or name) and values are reaction objects
        """
        reaction_objects = []

        for rxn_ref in reactions_ref.keys(): 
            rxn_object = []
            if reactions_ref[rxn_ref].lower() == 'id':
                rxn_object = [r for r in self.reactions if r.id.lower() == rxn_ref.lower()]
            elif reactions_ref[rxn_ref].lower() == 'name':
                rxn_object = [r for r in self.reactions if r.name.lower() == rxn_ref.lower()]
            else:
                raise userError("**Error! Invalid reaction reference '" + reactions_ref[rxn] + "'")

            if len(rxn_object) > 1:
                raise userError('More than one reaction object found for ' + rxn_ref)
            elif len(rxn_object) == 0:
                raise userError('No reaction object found for ' + rxn_ref)
            else:
                reaction_objects.append((rxn_ref,rxn_object[0])) 

        if len(reactions_ref.keys()) == 1:
            output_reactions = reaction_objects[0][1]
        else:
            output_reactions = dict(reaction_objects)

        return output_reactions 


    def get_metabolites(self,metabolites_ref):
        """
        Returns a dictionary of the selected metabolite objects

        INPUTS:
        ------- 
        metabolites_ref: A dictionary with keys and values as follows:
                            keys: A metabolite id or name (string) 
                          values: A string indicating the format by which the metabolite
                                  should be search. Eligible choices include 'id' and 'name'
                                  or 'formula'
                         Example: {'m1_id':'id','m2_name':'name','m3':'formula'}
        
        OUTPUTS:
        -------
        If the input is a dictionary with a single key the the ouput is the corresponding
        metabolite object, otherwise the output is a dictionary whose keys are the metabolite 
        reference values (id or name) and values are metabolite objects
        """
        metabolite_objects = []
        
        for metab_ref in metabolites_ref.keys(): 
            metab_object = []
            if metabolites_ref[metab_ref].lower() == 'id':
                metab_object = [m for m in self.metabolites if m.id.lower() == metab_ref.lower()]
            elif metabolites_ref[metab_ref].lower() == 'name':
                metab_object = [m for m in self.metabolites if m.name.lower() == metab_ref.lower()]
            elif metabolites_ref[metab_ref].lower() == 'formula':
                metab_object = [m for m in self.metabolites if m.formula.lower() == metab_ref.lower()]
            else:
                raise userError("**Error! Invalid metabolite reference '" + metabolites_ref[metab] + "'")

            if len(metab_object) > 1:
                raise userError('More than one reaction object found for ',metab)
            else:
                metabolite_objects.append((metab_ref,metab_object[0])) 

        if len(metabolites_ref.keys()) == 1:
            output_metabolites = metabolite_objects[0][1]
        else:
            output_metabolites = dict(metabolite_objects)

        return output_metabolites 


    def get_genes(self,genes_ref):
        """
        Returns a dictionary of the selected gene objects

        INPUTS:
        ------- 
        genes_ref: A dictionary with keys and values as follows:
                         keys: A gene id or name (string) 
                       values: A string indicating the format by which the genes
                               should be search. Eligible choices include 'id' and 'name'
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

    def add_metabolites(self,new_metabolites):
        """
        Add new metabolites to the model
        
        INPUTS:
        -------
        new_metabolites: A dictionary of new metabolite objects with
                         keys and values as follows                   
                           Keys: New metabolite objects
                         Values: Can be an empty dictionary '{}' if a new metabolite participates 
                                 in a new reaction that does not exist in the model yet.
                                 Note that new reactions are added to the model after
                                 creating new metabolite objects. 
                                 If this new metabolite participates in some existing 
                                 reactions in the model, then Values can be another 
                                 dictionary with keys and values as follows: 
                                   keys: Reaction objects, showing existing reactions in 
                                         the model in which the metabolite participates 
                                         as a reactant or product
                                 Values: The stoichiometric coefficient of the metabolite
                                         in existing reactions in the model
        
        OUTPUTS:
        --------
        There is no actual output. The model is updated with the addition of
        new metabolites 
        """
        self.metabolites += new_metabolites.keys()

        # Sort according to id
        self.metabolites = sorted(list(set(self.metabolites)),key=lambda x:x.id)

        # Update metabolites by id
        self.metabolites_by_id = dict([(metab.id,metab) for metab in self.metabolites])    

        # Update the stoichiometry of existing reactions in the model 
        for metab in [m for m in new_metabolites.keys() if new_metabolites[m] != {}]:
            for rxn in new_metabolites[metab].keys():
                rxn.stoichiometry[metab] = new_metabolites[metab][rxn] 

    def add_reactions(self,new_reactions):
        """
        Add new reactions to the model
        
        INPUTS:
        -------
        new_reactions: A list of new reaction objects
        
        OUTPUTS:
        --------
        There is no actual output. The model is updated with the addition of
        new reactions 
        """
        self.reactions += new_reactions

        # Sort according to id
        self.reactions = sorted(list(set(self.reactions)),key=lambda x:x.id)

        # Update reactions_by_id
        self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    

        # Add to the model any metagbolite objects that serves as a reactant
        # or the product of this reaction, which is not present in the model
        new_metabolites = [m for r in new_reactions for m in r.metabolites if m not in self.metabolites]
        if len(new_metabolites) > 0:
            self.add_metabolites(new_metabolites)

        # Update metabolite objects participating in these new reactions
        for rxn in new_reactions:
            for metab in rxn.metabolites:
                if rxn not in metab.reactions:
                    metab.reactions.append(rxn)
                    if rxn.stoichiometry[metab] < 0:
                        metab.reactant_reactions.append(rxn) 
                    elif rxn.stoichiometry[metab] > 0:
                        metab.product_reactions.append(rxn) 

    def remove_metabolites(self,removed_metabolites):
        """
        Remove selected metabolites from the model
        
        INPUTS:
        ------
        removed_metabolites: A list of metabolite objectc that must be removed 
        """
        if type(removed_metabolites) is not list:
            userError("**Error! the input to 'remove_metabolites' should be a list of metabolite objects")

        # First remove this metabolite from all relevant reaction objects
        for metab in removed_metabolites:
            for rxn in metab.reactions:
                if rxn.stoichiometry[metab] < 0:
                    if rxn.reactants != []:
                        del rxn.reactants[rxn.reactants.index(metab)]
                elif rxn.stoichiometry[metab] > 0:
                    if rxn.products != []:
                        del rxn.products[rxn.products.index(metab)]

                del rxn.stoichiometry[metab]

            # Check its presence in all other reactions
            metab_rxns = list(set([r for r in self.reactions if metab in r.metabolites + r.reactants + r.products + r.stoichiometry.keys()])) 
            if len(metab_rxns) > 0:
                for r in metab_rxns:
                    try:
                        del r.metabolites[r.metabolites.index(metab)]
                    except:
                        continue
                    try:
                        del r.reactants[r.reactant.index(metab)]
                    except:
                        continue
                    try:
                        del r.products[r.products.index(metab)]
                    except:
                        continue
                    try:
                        del r.stoichiometry.keys()[r.stoichiometry.keys().index(metab)]
                    except:
                        continue

            # Remove it from the list of metabolites
            del self.metabolites[self.metabolites.index(metab)] 
            self.metabolites = sorted(set(list(self.metabolites)),key=lambda x:x.id)

            # Remove it from metabolites_by_id
            del self.metabolites_by_id[metab.id]

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
            for metab in rxn.metabolites:
                del metab.reactions[metab.reactions.index(rxn)]
                if rxn.stoichiometry[metab] < 0 and metab.reactant_reactions != []:
                    del metab.reactant_reactions[metab.reactant_reactions.index(rxn)]
                elif rxn.stoichiometry[metab] > 0 and metab.product_reactions != []:
                    del metab.product_reactions[metab.product_reactions.index(rxn)]

            # Check the presence in all other metabolites
            rxn_metabs = list(set([m for m in self.metabolites if rxn in m.reactions + m.reactant_reactions + m.product_reactions])) 
            if len(rxn_metabs) > 0:
                for m in rxn_metabs:
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

    def validate(self):
        """
        Checks the possible errors in the model and fixes them 
        It returns a warning in the output for issues that cannot be fixed.

        This function performs the following checks:
        Issues that are fixed:
        ----------------------
          - Any metabolite object replicates in the list of metabolties
          - Any reaction object replicates in the list of reactions    
          - Any metabolite participates in at elast one reaction as a 
            reactant or product. A warningn is issued, in case this issue 
            cannot be fixed 
          - Any reaction contains at least one participating metabolite.
            A warning is issued, in case this issue cannot be fixed 
          - Any metabolite objects that are referred to reaction
            objects, but are not in the list of metabolites in the model
          - Any metabolite objects that are referred to in reaction objects
            but that are not present in the list of metabolites in the model.
            The issue is fixed by adding the missing metabolites to the list 
            of metabolites in the model
          - Any reaction objects that are referred to in metabolite objects
            but that are not present in the list of reactions in the model.
            The issue is fixed by adding the missing reactions to the list 
            of reactions in the model

        Issues for whcih a warning is issued in the output:
        ----------------------
          - Any replicates in metabolite ids
          - Any replicates in reaction ids
 
        """
        print 'Validating the model ...'
        errors_in_model = False

        #-- Check problems with metabolites ---
        # Check for duplicates in the list of metabolite objects and fix them
        if len(set(self.metabolites)) < len(self.metabolites):
            errors_in_model = True
            self.metabolites = sorted(list(set(self.metabolites)),key=lambda x:x.id)
            print '   Duplicates in the list of metabolites were fixed'

        # Check if there are any replicates in metabolite ids and issue a warning
        metab_ids = [m.id for m in self.metabolites]
        if len(set(metab_ids)) < len(metab_ids):
            errors_in_model = True

            # Count how many times each metabolite id has been repated
            metab_id_counts = dict((id,metab_ids.count(id)) for id in  metab_ids)

            # metabolite ids repated more than once
            repeated_metab_ids = [(id,metab_id_counts[id]) for id in metab_id_counts.keys() if metab_id_counts[id] > 1]
            print '   The following metabolite ids are repeated: ',repeated_metab_ids 
            
        # Check if each metabolite participates in at least a reaction in the model
        no_rxn_metabs = [m for m in self.metabolites if len(m.reactions) == 0]
        fixed_metab = 0
        if len(no_rxn_metabs) > 0:
            errors_in_model = True
            for metab in no_rxn_metabs:
                # Try to fix 
                m_rxns = [r for r in self.reactions if metab in r.stoichiometry.keys()]
                if len(m_rxns) == 0: 
                    print "   WARNING! metabolite '",metab.id,"' does not participate in any reactions in the model"         
                else: # fix it
                    metab.reactions = m_rxns
                    metab.reactant_reactions = [r for r in self.reactions if r.stoichiometry[metab] < 0]
                    metab.product_reactions = [r for r in self.reactions if r.stoichiometry[metab] > 0]
                    fixed_metab += 1

        if fixed_metab > 0:
            print '   ',fixed_metab,' metabolites not participated in any reactions in the model were fixed'

        # metabolites that are being referred to in a reaction object but that are not
        # present in the lis of metabolites in the model
        missing_metabs = list(set([m for r in self.reactions for m in r.metabolites + r.reactants + r.products + r.stoichiometry.keys() if m not in self.metabolites])) 
        if len(missing_metabs) > 0:
            errors_in_model = True
            self.metabolites += missing_metabs
            self.metabolites = sorted(list(set(self.metabolites)),key=lambda x:x.id)
            self.metabolites_by_id = dict([(metab.id,metab) for metab in self.metabolites])    
            print '   ',len(missing_metabs),' metabolites were added to the list of metabolites'

        #-- Check problems with reactions ---
        # Check for duplicates in the list of reaction objects and fix them
        if len(set(self.reactions)) < len(self.reactions):
            self.reactions = sorted(list(set(self.reactions)),key=lambda x:x.id)
            print '   Duplicates in the list of reactions were fixed'

        # Check if there are ny replicates in reaction ids and issue a warning
        rxn_ids = [r.id for r in self.reactions]
        if len(set(rxn_ids)) < len(rxn_ids):
            # Count how many times each reaction id has been repated
            rxn_id_counts = dict((id,rxn_ids.count(id)) for id in  rxn_ids)

            # reaction ids repated more than once
            repeated_rxn_ids = [(id,rxn_id_counts[id]) for id in rxn_id_counts.keys() if rxn_id_counts[id] > 1]
            print '   The following reaction ids are repeated: ',repeated_rxn_ids 

        # Check if there are any reactions with no defined stoichiometry 
        no_stoic_rxns = [r.id for r in self.reactions if len(r.stoichiometry) == 0]
        if len(no_stoic_rxns) > 0:
            print "\n**WARNING! 'stoichiometry' is not defined for these reactions: ",no_stoic_rxns,'\n'

        # Check if each reaction has at least one participating metabolite
        no_metab_rxns = [r for r in self.reactions if len(r.metabolites) == 0]
        fixed_rxn = 0
        if len(no_metab_rxns) > 0:
            errors_in_model = True
            for rxn in no_metab_rxns:
                # Try to fix it 
                rxn.metabolites = rxn.stoichiometry.keys()
                fixed_rxn += 1

        if fixed_rxn > 0:
            print '   ',fixed_rxn," reactions with empty field 'metabolites' were fixed"

        # reactions that are being referred to in a metabolite object but that are not
        # present in the lis of reactions in the model
        missing_rxns = list(set([r for m in self.metabolites for r in m.reactions + m.reactant_reactions + m.product_reactions if r not in self.reactions])) 
        if len(missing_rxns) > 0:
            errors_in_model = True
            self.reactions += missing_rxns
            self.reactions = sorted(set(list(self.reactions)),key=lambda x:x.id)
            self.reactions_by_id = dict([(rxn.id,rxn) for rxn in self.reactions])    
            print '   ',len(missing_rxns),' missing reactions were added to the list of reactions'

        if errors_in_model == False:
            print '   Checked the model. All are OK'
 

    def print_reactions(self,ref_type = 'id', print_equation = True, metab_ref = 'id'):
        """
        Prints in the output the list of all reactions and in the model 
        and their equations with the format specified by ref_type 
 
        INPUTS:
        -------
              ref_type: A string indicating the in what format the reactions should be 
                        printed in the output. Current eligible choices are 'id' or 'name'.
                        Default is 'id'.
        print_equation: If True reaction equation is also printed wherein metabolites appear 
                        with the format specified by metab_ref  
             metab_ref: A string indicating in what format the metabolites should appear
                        in a reaction equation if ref_type = 'equation'. If None, metabolite
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
                rxn = rxn + ':   ' + reaction.get_equation(ref_type = metab_ref)

            print rxn
        return ''

    def print_metabolites(self,ref_type = 'id'):
        """
        Prints in the output the list of all metabolites in the model 
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

            print mm
        return ''

    def export_to_cobra(self):
        """
        Export the model to a COBRApy model 
        """
        # Create an empty COBRA model
        cobra_model = cobra.Model(self.id)

        # Create metabolite objects
        for metab in self.metabolites:
            cobra_metab = cobra.Metabolite(id = metab.id, name = metab.name, formula = metab.formula)
            cobra_metab.notes = metab.notes
            if metab.compartment != None:
                cobra_metab.compartment = metab.compartment.id
            cobra_model.add_metabolites(cobra_metab)

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
                if len(rxn.genes) > 0:
                    cobra_rxn.gene_reaction_rule = ' or '.join([g.id for g in rxn.genes])
            cobra_model.add_reaction(cobra_rxn)

        return cobra_model

    def fba(self,optimization_solver = 'gurobi', create_model = True, store_opt_fluxes = True, flux_key = None, run = True, assign_wildType_max_biomass = False, screen_output = 'on'):
        """
        Creates a fba model for this model and runs fba 
        All inputs, except run, are the same as fba
                                run: Set to False to create a fba model 
                                     but do not run fba 
        assign_wildType_max_biomass: If True the optimal objective function
                                     value of the fba is assigned to global 
                                     variable wildType_max_biomass
        """
        if create_model == True:
            self.fba_model = fba(model = self,optimization_solver = optimization_solver, create_model = create_model, flux_key = flux_key, store_opt_fluxes = store_opt_fluxes, screen_output = screen_output)
        else:
            self.fba_model.create_model = create_model
            self.fba_model.flux_key = flux_key
            self.fba_model.optimization_solver = optimization_solver
            self.fba_model.store_opt_fluxes = store_opt_fluxes
            self.fba_model.screen_output = screen_output 
            # Update the flux bounds
            for j in self.reactions:
                self.fba_model.pyomo_fbaModel.v[j].lb = j.flux_bounds[0]
                self.fba_model.pyomo_fbaModel.v[j].ub = j.flux_bounds[1]

            # Update the objective function
            self.fba_model.pyomo_fbaModel.del_component('objectiveFunc')
            self.fba_model.pyomo_fbaModel.objectiveFunc = Objective(rule=self.fba_model.objectiveFunc_rule, sense = maximize)

        if run == True:
            self.fba_model.run()
            if assign_wildType_max_biomass == True:
                if self.fba_model.solution['exit_flag'] == 'globallyOptimal':
                    self.wildType_max_biomass = self.fba_model.solution['objective_value']


