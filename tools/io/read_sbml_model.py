from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from tools.core.organism import organism
from tools.core.compartment import compartment
from tools.core.gene import gene
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
import cobra
from libsbml import readSBML
import re

def read_sbml_model(file_name,model_organism, model_id = '', model_type = 'cmpolic',model_name = '', import_params = False, use_cobra_parser = True, remove_M_R = True, validate = True, warnings = True, stdout_msgs = True): 
    """
    A class to read a sbml model into a model object. The sbml model is 
    first converted into a COBRApy model, which is used as an input to create
    the corresponding model object. 

    INPUTS: 
    ------
           file_name: The name of the sbml file containing the model (String)
                      This string can contain the full path to the file  
      model_organism: Can be either a string containing the name of the organism or it can be 
                      an instance of the object organism
            model_id: Name of the model (string) 
          model_type: Type of the model (string, e.g., 'cmpolic')
       import_params: Indicates whether the model parameters (flux bounds, objective coefficients)  
                      should be imported from the sbml file (True) or not (False). The default is zero
    use_cobra_parser: If true uses cobra SBML parser. If not, parses by itself
          remove_M_R: If true, 'M_' and 'R_' are removed from the begining of metabolite
                      and reaction ids. This input is used only if use_cobra_parser = False
            warnings: Can be 'on' or 'off' shwoing whether the warnings should be writtten to the 
                      screen or not
         stdout_msgs: Can be 'on' or 'off' shwoing whether any messages should be written to the
                      screen

    OUTPUT:
    -------
               model: A model object

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 12-15-2015
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError("stdout_msgs must be True or False")
    if not isinstance(warnings,bool):
        raise TypeError("warnings must be True or False")
    if not isinstance(validate,bool):
        raise TypeError("validate must be True or False")


    # organism
    if type(model_organism) is str:             # If name is provided
        model_organism = organism(id = model_organism)
    elif type(model_organism) is not organism:  # if an object of type organism is not provided
        raise userError('Model_organism must be either a string containing the name (id) of the organism or an instance of the class "organism"')

    compartments = []
    compartments_by_id = []
    genes = []
    genes_by_id = []
    compounds = []
    compounds_by_id = []
    reactions = []

    #--------- Use cobra parser ------
    if use_cobra_parser:
        # Import the COBRApy model
        cobra_model = cobra.io.read_sbml_model(file_name) 
    
        cobra_model.optimize()
        for rxn in [r for r in cobra_model.reactions if len(r.metabolites) == 0]:
            if warnings:
                print 'WARNING! ',rxn.id,' in the original sbml file has no associated compounds: ',rxn.reaction
    
        #-- Compartments --
        for compartment_id in sorted(cobra_model.compartments.keys()):
            comp = compartment(id = compartment_id, name = cobra_model.compartments[compartment_id])
            compartments.append(comp)
            compartments_by_id.append((compartment_id,comp))
        compartments_by_id = dict(compartments_by_id)
    
        #-- Genes --
        for gene_cobra in sorted(cobra_model.genes,key=lambda x:x.id):
            if gene_cobra.notes != {}:
                notes = gene_cobra.notes  
            else:
                notes = None
            g = gene(id = gene_cobra.id, name = gene_cobra.name, notes = notes)
            genes.append(g)
            genes_by_id.append((gene_cobra.id,g))
        genes_by_id = dict(genes_by_id)
    
        #-- Compounds --
        # Note that the field 'reactions' for each compound is added
        # after constructing the reaction objects
        for cmp in sorted(cobra_model.metabolites, key = lambda x:x.id):
            if cmp.notes != {}:
                notes = cmp.notes
            else:
                notes = None
            m = compound(id = cmp.id, compartment = compartments_by_id[cmp.compartment], name = cmp.name, formula = cmp.formula.id, molecular_weight = cmp.formula.weight,charge = cmp.charge, notes = notes)
            compounds.append(m)
            compounds_by_id.append((cmp.id,m))
        compounds_by_id = dict(compounds_by_id)
    
        #-- Reactions --
        for rxn in sorted(cobra_model.reactions, key = lambda x:x.id):
            r_is_exchange = False
            # Reaction reversibility
            #** Fix this after removing 'exchange' from reaction reversiblities
            if rxn.id.find('EX_') == 0 or 'exchange' in rxn.name.lower():
                rxn_rev = 'exchange'
                r_is_exchange = True
            elif rxn.reversibility == False: 
                rxn_rev = 'irreversible'
            elif rxn.reversibility == True: 
                rxn_rev = 'reversible'
    
            # Reaction stoichiometry
            stoichiometry = []
            for cmp in rxn.metabolites: 
                stoichiometry.append((compounds_by_id[cmp.id],rxn.get_coefficient(cmp)))   
    
            if len(stoichiometry) == 0 and warnings:
                print "WARNING! The field 'stoichiometry' was not assigned for reaction " + rxn.id
    
            # Genes for this reaction 
            r_genes = []
            for cobra_gene in rxn.genes:
                r_genes.append(genes_by_id[cobra_gene.id])
    
            if import_params == True:
                r_flux_bounds = [rxn.lower_bound,rxn.upper_bound]
                r_obj_coeff = rxn.objective_coefficient
                r_notes = rxn.notes
            else:
                r_flux_bounds = []
                r_obj_coeff = None
                r_notes = ''

            r = reaction(id = rxn.id, stoichiometry = dict(stoichiometry), reversibility = rxn_rev, name = rxn.name, subsystem = rxn.subsystem, genes = r_genes, gene_reaction_rule = rxn.gene_reaction_rule, is_exchange = r_is_exchange, objective_coefficient = r_obj_coeff, flux_bounds = r_flux_bounds, notes = r_notes)
            reactions.append(r)
    
    #--------- Do not use cobra parser ------
    else:
        sbml_doc = readSBML(file_name)
        if sbml_doc.getNumErrors() < 0 and warnings:
            print 'WARNING! ',sbml_doc.getNumErrors(),' errors in reading SBML file ...'
            # print sbml_doc.printErrors()

        if stdout_msgs:
            print 'SBML version ',sbml_doc.getVersion() 
        model_notes = sbml_doc.getNotesString()
        sbml_model = sbml_doc.getModel()

        if model_id == '' and sbml_model.getId() != '':
            model_id = sbml_model.getId()
        elif model_id == '' and sbml_model.getId() == '':
            raise userError('Please provide a model id') 
 
        if model_name == '':
            model_name = sbml_model.getName() 

        #-- Compartments --
        for sbml_compart in sbml_model.getListOfCompartments():
            comp = compartment(id = sbml_compart.getId(), name = sbml_compart.getName(),notes = sbml_compart.getNotesString())
            compartments.append(comp)
            compartments_by_id.append((sbml_compart.getId(),comp))
        compartments_by_id = dict(compartments_by_id)
    
        #-- Compounds --
        # Note that the field 'reactions' for each compound is added
        # after constructing the reaction objects
        for cmp in sorted([cmp for cmp in sbml_model.getListOfSpecies() if cmp.boundary_condition == False], key = lambda x:x.id):
            if remove_M_R:
                cmp_id = re.sub('M_','',cmp.getId())
            else:
                cmp_id = cmp.getId()

            c_notes_dict = parse_sbml_notes(cmp.getNotesString())
            # Charge (from sbml.py of COBRApy)
            if 'CHARGE' in c_notes_dict.keys():
                cmp_charge = c_notes_dict["CHARGE"][0]
                # Remove this element from the dict as we don't need it anymore
                del c_notes_dict['CHARGE']
                try:
                    cmp_charge = float(cmp_charge)
                    if cmp_charge == int(cmp_charge):
                       cmp_charge = int(cmp_charge)
                except:
                    if warnings:
                        print 'WARNING! Charge of ',cmp_id,' is not a number (',cmp_charge
                    cmp_charge = 0 
            else:
                cmp_charge = 0 

            # Formula
            formula_key = [k for k in c_notes_dict.keys() if k.lower() == 'formula']
            if len(formula_key) >= 1:
                if len(formula_key) > 1 and warnings:
                    print 'WARNING! More than two fields was detected for formula for compound ',cmp_id,': ',formula_key,'. The first element was used to extract the formula'
                formula_key = formula_key[0] 
                cmp_formula = c_notes_dict[formula_key]
                # Remove this element from the dict as we don't need it anymore
                del c_notes_dict[formula_key]
            else:
                cmp_formula = ''
                
            c = compound(id = cmp_id, compartment = compartments_by_id[cmp.getCompartment()], name = cmp.getName(), charge = cmp_charge, formula = cmp_formula, notes = str(c_notes_dict))
            compounds.append(c)
            compounds_by_id.append((c.id,c))

        compounds_by_id = dict(compounds_by_id)

        #-- Genes --
        # Genes are added using gene-reaction rules when parsing reactions
        genes_by_id = {}

        #-- Reactions --
        for rxn in sorted(sbml_model.getListOfReactions(), key = lambda x:x.id):
            if remove_M_R:
                rxn_id = re.sub('R_','',rxn.getId())
            else:
                rxn_id = rxn.getId()

            # Reaction type
            #** Fix this after removing 'exchange' from reaction reversiblities
            r_is_exchange = False
            if 'EX_' in rxn_id or 'EX_' in rxn.getName() or 'exchange' in rxn.getName().lower():
                r_rev = 'exchange'
                r_is_exchange = True
            elif rxn.reversible == False: 
                r_rev = 'irreversible'
            elif rxn.reversible == True: 
                r_rev = 'reversible'
    
            # Reaction stoichiometry and compartments
            r_stoichiometry = {}
            r_comparts = []
            for cmp in rxn.getListOfReactants(): 
                if remove_M_R:
                    cmp_id = re.sub('M_','',cmp.getSpecies())
                else:
                    cmp_id = cmp.getSpecies()
                # if cmp_id not in the list of compounds it must be a metabolite with a True boundary_condition
                if cmp_id in compounds_by_id.keys():
                    r_stoichiometry[compounds_by_id[cmp_id]] = -cmp.getStoichiometry()
                    r_comparts.append(compounds_by_id[cmp_id].compartment)
            for cmp in rxn.getListOfProducts(): 
                if remove_M_R:
                    cmp_id = re.sub('M_','',cmp.getSpecies())
                else:
                    cmp_id = cmp.getSpecies()
                # if cmp_id not in the list of compounds it must be a metabolite with a True boundary_condition
                if cmp_id in compounds_by_id.keys():
                    r_stoichiometry[compounds_by_id[cmp_id]] = cmp.getStoichiometry()
                    r_comparts.append(compounds_by_id[cmp_id].compartment)
 
            if len(r_stoichiometry) == 0 and warnings:
                print "WARNING! The field 'stoichiometry' was not assigned for reaction " + rxn_id
    
            # List of compartments for the reaction
            r_comparts = list(set(r_comparts)) 

            # Get parameters
            r_param_dict = {}
            if rxn.getKineticLaw():
                for param in rxn.getKineticLaw().getListOfParameters():
                    r_param_dict[param.getId().lower()] = param.getValue()

                # Bounds
                if 'lower_bound' in r_param_dict:
                    r_LB = r_param_dict['lower_bound']
                elif 'lower bound' in r_param_dict:
                    r_LB = r_param_dict['lower bound']
                else:
                    r_LB = None

                if 'upper_bound' in r_param_dict:
                    r_UB = r_param_dict['upper_bound']
                elif 'upper bound' in r_param_dict:
                    r_UB = r_param_dict['upper bound']
                else:
                    r_UB = None
                if r_LB == None or r_UB == None:
                    r_flux_bounds = []
 
                # Objective coefficient
                if 'objective_coefficient' in r_param_dict:
                    r_obj_coeff = r_param_dict['objective_coefficient']
                elif 'objective coefficient' in r_param_dict:
                    r_obj_coeff = r_param_dict['objective coefficient']
                else:
                    r_obj_coeff = None

            else:
                r_flux_bounds = []
                r_obj_coeff = None
        
            r_notes_dict = parse_sbml_notes(rxn.getNotesString())
            if 'SUBSYSTEM' in r_notes_dict.keys():
                r_subsystem = r_notes_dict['SUBSYSTEM']
                del r_notes_dict['SUBSYSTEM']
            else:
                r_subsystem = None

            if 'SYNONYMS' in r_notes_dict.keys():
                r_name = [rxn.getName()] + r_notes_dict['SYNONYMS']
                del r_notes_dict['SYNONYMS']
            else:
                r_name = rxn.getName()

            if 'EC NUMBER' in r_notes_dict.keys():
                r_EC_num = r_notes_dict['EC NUMBER'] 
                del r_notes_dict['EC NUMBER'] 
            else:
                r_EC_num = None

            if 'GENE ASSOCIATION' in r_notes_dict.keys():
                r_gene_rxn_rule = r_notes_dict['GENE ASSOCIATION'][0]
                # List of gene ids for this reaction
                r_genes_list = parse_grr(r_gene_rxn_rule)
                # List of gene objects for this reaction
                r_genes = []
                for g_id in r_genes_list:
                    # Add the gene to the model if it is not already there
                    if len(genes_by_id) == 0 or g_id not in genes_by_id.keys():
                        g = gene(id = g_id)
                        genes.append(g)
                        genes_by_id[g_id] = g
                        r_genes.append(g)
                    else:
                        r_genes.append(genes_by_id[g_id])

                del r_notes_dict['GENE ASSOCIATION']
            else:
                r_gene_rxn_rule = None

            if 'CONFIDENCE LEVEL' in r_notes_dict.keys():
                r_confidence = r_notes_dict['CONFIDENCE LEVEL'][0]
                del r_notes_dict['CONFIDENCE LEVEL']
            else:
                r_confidence = None

            r = reaction(id = rxn_id, stoichiometry = r_stoichiometry, reversibility = r_rev, is_exchange = r_is_exchange, name = r_name, EC_number = r_EC_num, subsystem = r_subsystem, compartment = r_comparts, genes = r_genes, gene_reaction_rule = r_gene_rxn_rule, objective_coefficient = r_obj_coeff, flux_bounds = [r_LB,r_UB])
            reactions.append(r)

    #--- Biomass and ATPM ---
    # Now add the field 'reactions' for each compound 
    for cmp in compounds:
       cmp.reactions = [r for r in reactions if cmp in r.stoichiometry.keys()]        
       cmp.reactant_reactions = [r for r in cmp.reactions if r.stoichiometry[cmp] < 0]        
       cmp.product_reactions = [r for r in cmp.reactions if r.stoichiometry[cmp] > 0]        
    
    # biomass reactions
    biomass_reactions = [r for r in reactions if r.id.lower().find('biomass') >= 0 or len([n for n in r.name if 'biomass' in n.lower() or 'growh' in n.lower()]) > 0]
    if len(biomass_reactions) == 0:
        biomass_reaction = None
        if warnings:
            print 'WARNING! No biomass reactions detected. Check the model manually\n'
    elif len(biomass_reactions) > 1:
        biomass_reaction = None
        if stdout_msgs:
            print '\nMore than one biomass reactions detected: ',[r.id for r in biomass_reactions],'\nPlease choose one manually.\n'
    elif len(biomass_reactions) == 1:
        biomass_reaction = biomass_reactions[0]
        if stdout_msgs:
            print '\nA biomass reaction detected: ',biomass_reaction.id,'\n'
    
    # ATP maintenance reaction 
    atpm_reaction = [r for r in reactions if 'atpm' in r.id.lower() or len([n for n in r.name if 'atpm' in n.lower() and 'maintenance' in n.lower()]) > 0 or len([n for n in r.name if 'ngam' in n.lower()]) > 0]
    if len(atpm_reaction) == 0 and warnings:
        print 'WARNING! No NGAM  reaction found. Check the model manually\n'
        atpm_reaction = None
    elif len(atpm_reaction) > 1 and warnings:
        print 'WARNING! More than one NGAM reactions found: ',[r.id for r in atpm_reaction],'\n'
    elif len(atpm_reaction) == 1:
        atpm_reaction = atpm_reaction[0] 

    #--- model ---
    return model(id = model_id, type = model_type, organism = model_organism, reactions = reactions, compounds = compounds, genes = genes, compartments = compartments, biomass_reaction = biomass_reaction,atpm_reaction = atpm_reaction, validate = validate, warnings = warnings, stdout_msgs = stdout_msgs)

def parse_sbml_notes(note_string, note_delimiter = ':'):
    """
    This is function parse_legacy_sbml_notes from sbml.py in COBRApy.
    It convers the notes string into a dictionary 
    """
    note_dict = {}
    start_tag = '<p>'
    end_tag = '</p>'
    if '<html:p>' in note_string:
        start_tag = '<html:p>'
        end_tag = '</html:p>'
    while start_tag in note_string and end_tag in note_string:
        note_start = note_string.index(start_tag)
        note_end = note_string.index(end_tag)
        the_note = note_string[(note_start + len(start_tag)):note_end].lstrip(' ').rstrip(' ')
        if note_delimiter in the_note:
            note_delimiter_index = the_note.index(note_delimiter)
            note_field = the_note[:note_delimiter_index].lstrip(' ').rstrip(' ').replace('_',' ').upper()
            note_value = the_note[(note_delimiter_index+1):].lstrip(' ').rstrip(' ')
            if note_field in note_dict:
                note_dict[note_field ].append(note_value)
            else:
                note_dict[note_field] = [note_value]
        note_string = note_string[(note_end+len(end_tag)): ]

    if 'CHARGE' in note_dict and note_dict['CHARGE'][0].lower() in ['none', 'na', 'nan']:
        note_dict.pop('CHARGE') #Remove non-numeric charges
        
    return(note_dict)

def parse_grr(grr):
    """
    Converts gene-reaction a rule to a list of genes
    grr:  A string containing the gene-reaction rule
    """
    # Remove ( and )
    grr = re.sub('\(|\)','',grr)

    # Remove and and or
    grr = re.sub(' and | AND | or | OR ',' ',grr)

    # Add a space to the start and beginning of grr
    grr = ' ' + grr + ' '
   
    # Replace two subsequent spaces with one space
    while '  ' in grr:
        grr = re.sub('  ',' ',grr)

    genes_list = []

    p = re.compile('(?<= ).*?(?= )')
    for g in p.finditer(grr):
        if g != None:
            genes_list.append(g.group())    

    return list(set(genes_list))

#------------ Sample run -------------
if __name__ == "__main__":
    pass
