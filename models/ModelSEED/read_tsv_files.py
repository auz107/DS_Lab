import sys, re
sys.path.append('/usr2/postdoc/alizom/work/')
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from imp import load_source
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
import copy
import scipy.special as sps
import numpy as np

def read_tsv_files_website(compounds_tsv_file,reactions_tsv_file):
    """
    Read the tsv files containing the model seed reactions and metabolites and creates a 
    python dictionary out of them 
  
    INPUTS:
    ------
      reactions_tsv_file: A tsv file downloaded from the ModelSEED database containing the 
                          list of all reacitons in this database.
    metabolites_tsv_file: A tsv file downloaded from the ModelSEED database containing the
                          list of all metabolites
    """

    
    # The following writes the metabolite and reaction information into a dicitonary and stores
    # the results into a file
    
    #---------------------------------------------------------
    #---------------------- Reactions ------------------------
    #---------------------------------------------------------
    print '\tParsing compound information ...'
    with open('ModelSEED_cpds_website.py','w') as f:
        f.write('ModelSEED_cpds = {}\n')
    
    compounds = {}
    
    with open(compounds_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                headers = line

                # Define new header ids
                new_headers = {}
                new_headers['Compound'] = 'ModelSEED_id'
                new_headers['Name'] = 'name'
                new_headers['Formula'] = 'formula'
                new_headers['Mass'] = 'mass'
                new_headers['KEGG maps'] = 'Kegg_maps'
                new_headers['KEGG CID'] = 'Kegg_id'

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of raw_metabolites.txt. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                if 'Compound' not in headers:
                    raise userError('**ERROR! "Compound" is not in the list of headers')
                else:
                    seed_mid_ind = headers.index('Compound')
                    compounds[line[seed_mid_ind]] = {}
                    for header in headers:
                        # index of header in headers
                        header_ind = headers.index(header)
                        if line[header_ind] == 'None' or line[header_ind] == ' ':
                            compounds[line[seed_mid_ind]][new_headers[header]] = None

                        # If the header is Name convert to a list, because sometimes there are multiple
                        # names for the smae compound separated by ,<br> 
                        elif header == 'Name':
                            p = re.compile('(?:(?<=^)|(?<=,\<br\>)).*?(?=$|,\<br\>)')
                            compounds[line[seed_mid_ind]][new_headers[header]] = re.findall(p,line[header_ind]) 

                        # Sometimes there are more than one kegg ids for each compound, separated by
                        # a comma. They are converted to a list
                        elif header == 'KEGG CID':
                            p = re.compile('(?:(?<=^)|(?<=,)).*?(?=$|,)')
                            compounds[line[seed_mid_ind]][new_headers[header]] = re.findall(p,line[header_ind]) 
                        else:
                            compounds[line[seed_mid_ind]][new_headers[header]] = line[header_ind] 

    # Write the results into a file   
    with open('ModelSEED_cpds_website.py','a') as f:
        for seed_id in compounds.keys():
            f.write("ModelSEED_cpds['" + seed_id + "'] = " + repr(compounds[seed_id]) + '\n')
    

    #--- Now create new dictionaries to access compounds by their name, Kegg id or formula --- 
    print '\tWriting compounds as dictionaries with names, clean names and Kegg ids as keys ...'
    from ModelSEED_cpds import ModelSEED_cpds

    # A dicitonary, which has compounds names as keys and returns the corresponding ModelSEED ids
    ModelSEED_cpds_by_name = {}
    for sid in ModelSEED_cpds.keys():
        for name in ModelSEED_cpds[sid]['name']:
            if name in ModelSEED_cpds_by_name.keys():
                ModelSEED_cpds_by_name[name] += [sid]
            else:
                ModelSEED_cpds_by_name[name] = [sid]

    with open('ModelSEED_cpds_website.py','a') as f:
        f.write('\nModelSEED_cpds_by_name = {}\n') 
        for name in ModelSEED_cpds_by_name.keys():
            if "'" in name:
                f.write('ModelSEED_cpds_by_name["' + name + '"] = ' + repr(ModelSEED_cpds_by_name[name]) + '\n') 
            else:
                f.write("ModelSEED_cpds_by_name['" + name + "'] = " + repr(ModelSEED_cpds_by_name[name]) + '\n') 

    # A dicitonary, which has compounds' clean names (names where all non-alphanumeric characters
    # are removed) as keys and returns the corresponding ModelSEED ids
    ModelSEED_cpds_by_clean_name = {}
    for sid in ModelSEED_cpds.keys():
        for name in ModelSEED_cpds[sid]['name']:
            clean_name = remove_non_alphanumeric(name).lower()
            if clean_name in ModelSEED_cpds_by_clean_name.keys():
                ModelSEED_cpds_by_clean_name[clean_name] += [sid]
            else:
                ModelSEED_cpds_by_clean_name[clean_name] = [sid]

    with open('ModelSEED_cpds_website.py','a') as f:
        f.write('\nModelSEED_cpds_by_clean_name = {}\n') 
        for clean_name in ModelSEED_cpds_by_clean_name.keys():
            # Use list(set()) to get only the unique clean names. This is because sometimes
            # the clean names for multiple names are identical and therefore the same ModelSEED id
            # is added more than once
            ModelSEED_cpds_by_clean_name[clean_name] = list(set(ModelSEED_cpds_by_clean_name[clean_name]))
            if "'" in name:
                f.write('ModelSEED_cpds_by_clean_name["' + clean_name + '"] = ' + repr(ModelSEED_cpds_by_clean_name[clean_name]) + '\n') 
            else:
                f.write("ModelSEED_cpds_by_clean_name['" + clean_name + "'] = " + repr(ModelSEED_cpds_by_clean_name[clean_name]) + '\n') 

    # A dicitonary, which has compounds Kegg ids as keys and returns the corresponding ModelSEED ids
    ModelSEED_cpds_by_KEGG_id = {}
    for sid in [sid for sid in ModelSEED_cpds.keys() if ModelSEED_cpds[sid]['Kegg_id'] != None]:
        for KEGG_id in ModelSEED_cpds[sid]['Kegg_id']:
            if KEGG_id in ModelSEED_cpds_by_KEGG_id.keys():
                ModelSEED_cpds_by_KEGG_id[KEGG_id] += [sid]
            else:
                ModelSEED_cpds_by_KEGG_id[KEGG_id] = [sid]

    with open('ModelSEED_cpds_website.py','a') as f:
        f.write('\nModelSEED_cpds_by_KEGG_id = {}\n') 
        for KEGG_id in ModelSEED_cpds_by_KEGG_id.keys():
            f.write("ModelSEED_cpds_by_KEGG_id['" + KEGG_id + "'] = " + repr(ModelSEED_cpds_by_KEGG_id[KEGG_id]) + '\n') 

    # A dicitonary, which has compounds formula as keys and returns the corresponding ModelSEED ids
    ModelSEED_cpds_by_formula = {}
    for sid in [sid for sid in ModelSEED_cpds.keys() if ModelSEED_cpds[sid]['formula'] != None]:
        formula = ModelSEED_cpds[sid]['formula']
        if formula in ModelSEED_cpds_by_formula.keys():
            ModelSEED_cpds_by_formula[formula] += [sid]
        else:
            ModelSEED_cpds_by_formula[formula] = [sid]

    with open('ModelSEED_cpds_website.py','a') as f:
        f.write('\nModelSEED_cpds_by_formula = {}\n') 
        for formula in ModelSEED_cpds_by_formula.keys():
            f.write("ModelSEED_cpds_by_formula['" + formula + "'] = " + repr(ModelSEED_cpds_by_formula[formula]) + '\n') 

    # A dicitonary, which has compounds clean formula (i.e., formula where all non-alpahnumeric 
    # characters are removed and all uppercase letters replaced with lower case) as keys and 
    # returns the corresponding ModelSEED ids
    ModelSEED_cpds_by_clean_formula = {}
    for sid in [sid for sid in ModelSEED_cpds.keys() if ModelSEED_cpds[sid]['formula'] != None]:
        clean_formula = remove_non_alphanumeric(ModelSEED_cpds[sid]['formula']).lower()
        if clean_formula in ModelSEED_cpds_by_clean_formula.keys():
            ModelSEED_cpds_by_clean_formula[clean_formula] += [sid]
        else:
            ModelSEED_cpds_by_clean_formula[clean_formula] = [sid]

    with open('ModelSEED_cpds_website.py','a') as f:
        f.write('\nModelSEED_cpds_by_clean_formula = {}\n') 
        for clean_formula in ModelSEED_cpds_by_clean_formula.keys():
            ModelSEED_cpds_by_clean_formula[clean_formula] = list(set(ModelSEED_cpds_by_clean_formula[clean_formula]))
            f.write("ModelSEED_cpds_by_clean_formula['" + clean_formula + "'] = " + repr(ModelSEED_cpds_by_clean_formula[clean_formula]) + '\n') 

    #------- Perform some diagnonsis tests ------------------
    print '\tPerforming some diagnostic tests on reactions ...'
    #-- Find compounds with identical name --
    non_unique_name_cpds = [n for n in ModelSEED_cpds_by_name.keys() if len(ModelSEED_cpds_by_name[n]) > 1]
    if len(non_unique_name_cpds) > 0:
        with open('non_unique_name_cpds_website.py','w') as f:
            f.write('# A dictionary where ModelSEED names and their corresponding ModelSEED ids serve as keys and values, respectively.\n\n') 
            f.write('\nnon_unique_name_cpds = {}\n') 
            for name in non_unique_name_cpds:
                if "'" in name:
                    f.write('non_unique_name_cpds["' + name + '"] = ' + repr(ModelSEED_cpds_by_name[name]) + '\n') 
                else:
                    f.write("non_unique_name_cpds['" + name + "'] = " + repr(ModelSEED_cpds_by_name[name]) + '\n') 

    #-- Find compounds with identical formula --
    non_unique_formula_cpds = [f for f in ModelSEED_cpds_by_formula.keys() if len(ModelSEED_cpds_by_formula[f]) > 1]
    if len(non_unique_formula_cpds):
        with open('non_unique_formula_cpds_website.py','w') as f:
            f.write('# A dictionary where ModelSEED formulas and their corresponding ModelSEED ids serve as keys and values, respectively.\n\n') 
            f.write('\nnon_unique_formula_cpds = {}\n') 
            for formula in non_unique_formula_cpds:
                f.write("non_unique_formula_cpds['" + formula + "'] = " + repr(ModelSEED_cpds_by_formula[formula]) + '\n') 

    #----- Write the list of all compounds with the same name or the same clean name ----
    # and the same formula in a list to show that they are actually the same compounds
    print '\tFinding equivalent compounds ...'
    # Compounds with the same clean names
    equivalent_cpds = []
    if len(non_unique_name_cpds) > 0:
        # A dictionary with ModelSEED ids as keys and clean formulas as values
        equiv_dict = {}
        for n in non_unique_name_cpds:
            for sid in ModelSEED_cpds_by_name[n]:
                if ModelSEED_cpds[sid]['formula'] != None:
                    equiv_dict[sid] = ModelSEED_cpds[sid]['formula'].lower()
        unique_formulas = list(set(equiv_dict.values())) 
        for uf in unique_formulas:
            eq_cmps = [sid for sid in equiv_dict.keys() if equiv_dict[sid] == uf]
            if len(eq_cmps) > 1:
                equivalent_cpds.append(eq_cmps)
   
    print 'equivalent_cpds = ',equivalent_cpds
    # Compounds with the same clean names
    non_unique_clean_name_cpds = [n for n in ModelSEED_cpds_by_clean_name.keys() if len(ModelSEED_cpds_by_clean_name[n]) > 1]
    print '# of non_unique_clean_name_cpds =',len(non_unique_clean_name_cpds) 
    if len(non_unique_clean_name_cpds) > 0:
        # A dictionary with ModelSEED ids as keys and clean formulas as values
        equiv_dict = {}
        for n in non_unique_clean_name_cpds:
            for sid in ModelSEED_cpds_by_clean_name[n]:
                if ModelSEED_cpds[sid]['formula'] != None:
                    equiv_dict[sid] = ModelSEED_cpds[sid]['formula'].lower()
        unique_formulas = list(set(equiv_dict.values())) 
        for uf in unique_formulas:
            eq_cmps = [sid for sid in equiv_dict.keys() if equiv_dict[sid] == uf]
            if len(eq_cmps) > 1 and eq_cmps not in equivalent_cpds:
                equivalent_cpds.append(eq_cmps)

    print 'equivalent_cpds = ',equivalent_cpds
    with open('ModelSEED_cpds_website.py','a') as f:
        f.write('\nequivalent_cpds = [\n') 
        if len(equivalent_cpds) > 0:
            for eq_cmps in equivalent_cpds:
                if equivalent_cpds.index(eq_cmps) == len(equivalent_cpds) - 1:
                    f.write(repr(eq_cmps) + '\n') 
                else:
                    f.write(repr(eq_cmps) + ',\n') 
        f.write(']\n') 


    print '\nMetabolites were written in a python dictionary stored in ModelSEED_cpds_website.py\n'   
 
    #---------------------------------------------------------
    #---------------------- Reactions ------------------------
    #---------------------------------------------------------
    print '\tParsing reaction information ...'
    with open('ModelSEED_rxns_website.py','w') as f:
        f.write('ModelSEED_rxns = {}\n')
    
    reactions = {}
    
    with open(reactions_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                headers = line

                # Define new header ids
                new_headers = {}
                new_headers['Reaction'] = 'ModelSEED_id'
                new_headers['Name'] = 'name'
                new_headers['Equation'] = 'equation'
                new_headers['Roles'] = 'roles'
                new_headers['Subsystems'] = 'subsystems'
                new_headers['KEGG maps'] = 'KEGG_maps'
                new_headers['Enzyme'] = 'EC_number'
                new_headers['KEGG RID'] = 'Kegg_id'

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of raw_metabolites.txt. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                if 'Reaction' not in headers:
                    raise userError('**ERROR! "Reaction" is not in the list of headers')
                else:
                    seed_mid_ind = headers.index('Reaction')
                    reactions[line[seed_mid_ind]] = {}
                    for header in headers:
                        # index of header in headers
                        header_ind = headers.index(header)
                        if line[header_ind] == 'None' or line[header_ind] == ' ':
                            reactions[line[seed_mid_ind]][new_headers[header]] = None

                        # Remove and from reaction names
                        elif header == 'Name' and ' and' in line[header_ind]:
                            reactions[line[seed_mid_ind]][new_headers[header]] = re.sub(' and$','',line[header_ind]) 

                        # Sometimes there are more than one Kegg id for each reaction
                        # spearated by a comma. They are converted to a list
                        elif header == 'KEGG RID':
                            p = re.compile('(?:(?<=^)|(?<=,)).*?(?=$|,)')
                            reactions[line[seed_mid_ind]][new_headers[header]] = re.findall(p,line[header_ind]) 

                        # Sometimes there are more than one subsystems for each reaction, separated by
                        # "|". They are converted to a list
                        elif header == 'Subsystems':
                            p = re.compile('(?:(?<=^)|(?<=\|)).*?(?=$|\|)')
                            reactions[line[seed_mid_ind]][new_headers[header]] = re.findall(p,line[header_ind]) 

                        # Sometimes there are more than one Enzymes (EC numbers) for each reaction
                        # spearated bya comma. They are converted to a list
                        elif header == 'Enzyme':
                            p = re.compile('(?:(?<=^)|(?<=,)).*?(?=$|,)')
                            reactions[line[seed_mid_ind]][new_headers[header]] = re.findall(p,line[header_ind]) 

                        else:
                            reactions[line[seed_mid_ind]][new_headers[header]] = line[header_ind] 

    # Write the results into a file   
    with open('ModelSEED_rxns_website.py','a') as f:
        for seed_id in reactions.keys():
            f.write("ModelSEED_rxns['" + seed_id + "'] = " + repr(reactions[seed_id]) + '\n')

    #--- Now create new dictionaries to access compounds by their name, Kegg id or formula --- 
    print '\tWriting reactions as dictionaries with names, clean names and Kegg ids as keys ...'
    from ModelSEED_rxns import ModelSEED_rxns

    # A dicitonary, which has reactions names as keys and returns the corresponding ModelSEED ids
    ModelSEED_rxns_by_name = {}
    for sid in ModelSEED_rxns.keys():
        # Note that reaction names are always a string (there is only one name) 
        name = ModelSEED_rxns[sid]['name']
        if name in ModelSEED_rxns_by_name.keys():
            ModelSEED_rxns_by_name[name] += [sid]
        else:
            ModelSEED_rxns_by_name[name] = [sid]

    with open('ModelSEED_rxns_website.py','a') as f:
        f.write('\nModelSEED_rxns_by_name = {}\n') 
        for name in ModelSEED_rxns_by_name.keys():
            if "'" in name:
                f.write('ModelSEED_rxns_by_name["' + name + '"] = ' + repr(ModelSEED_rxns_by_name[name]) + '\n') 
            else:
                f.write("ModelSEED_rxns_by_name['" + name + "'] = " + repr(ModelSEED_rxns_by_name[name]) + '\n') 

    # A dicitonary, which has reactions clean names (i.e., name where all non-alphanumeric characters
    # are removed)  as keys and returns the corresponding ModelSEED ids
    ModelSEED_rxns_by_clean_name = {}
    for sid in ModelSEED_rxns.keys():
        # Note that reaction names are always a string (there is only one name) 
        clean_name = remove_non_alphanumeric(ModelSEED_rxns[sid]['name']).lower()
        if clean_name in ModelSEED_rxns_by_clean_name.keys():
            ModelSEED_rxns_by_clean_name[clean_name] += [sid]
        else:
            ModelSEED_rxns_by_clean_name[clean_name] = [sid]

    with open('ModelSEED_rxns_website.py','a') as f:
        f.write('\nModelSEED_rxns_by_clean_name = {}\n') 
        for clean_name in ModelSEED_rxns_by_clean_name.keys():
            # Use list(set()) to get only the unique clean names. This is because sometimes
            # the clean names for multiple names are identical and therefore the same ModelSEED id
            # is added more than once
            ModelSEED_rxns_by_clean_name[clean_name] = list(set(ModelSEED_rxns_by_clean_name[clean_name]))
            if "'" in clean_name:
                f.write('ModelSEED_rxns_by_clean_name["' + clean_name + '"] = ' + repr(ModelSEED_rxns_by_clean_name[clean_name]) + '\n') 
            else:
                f.write("ModelSEED_rxns_by_clean_name['" + clean_name + "'] = " + repr(ModelSEED_rxns_by_clean_name[clean_name]) + '\n') 


    # A dicitonary, which has reactions Kegg ids as keys and returns the corresponding ModelSEED ids
    ModelSEED_rxns_by_KEGG_id = {}
    for sid in [sid for sid in ModelSEED_rxns.keys() if ModelSEED_rxns[sid]['Kegg_id'] != None]:
        for KEGG_id in ModelSEED_rxns[sid]['Kegg_id']:
            if KEGG_id in ModelSEED_rxns_by_KEGG_id.keys():
                ModelSEED_rxns_by_KEGG_id[KEGG_id] += [sid]
            else:
                ModelSEED_rxns_by_KEGG_id[KEGG_id] = [sid]

    with open('ModelSEED_rxns_website.py','a') as f:
        f.write('\nModelSEED_rxns_by_KEGG_id = {}\n') 
        for KEGG_id in ModelSEED_rxns_by_KEGG_id.keys():
            f.write("ModelSEED_rxns_by_KEGG_id['" + KEGG_id + "'] = " + repr(ModelSEED_rxns_by_KEGG_id[KEGG_id]) + '\n') 

    #------- Perform some diagnonsis tests ------------------
    print '\tPerforming some diagnostic tests on reactions ...'
    #-- Find reactions with identical name --
    non_unique_name_rxns = [n for n in ModelSEED_rxns_by_name.keys() if len(ModelSEED_rxns_by_name[n]) > 1]
    if len(non_unique_name_rxns) > 0:
        with open('non_unique_name_rxns_website.py','w') as f:
            f.write('# A dictionary where ModelSEED names and their corresponding ModelSEED ids serve as keys and values, respectively.\n\n') 
            f.write('\nnon_unique_name_rxns = {}\n') 
            for name in non_unique_name_rxns:
                if "'" in name:
                    f.write('non_unique_name_rxns["' + name + '"] = ' + repr(ModelSEED_rxns_by_name[name]) + '\n') 
                else:
                    f.write("non_unique_name_rxns['" + name + "'] = " + repr(ModelSEED_rxns_by_name[name]) + '\n') 

    #----- Write the list of all reactions with the same name or the same clean name ----
    # and the same equation in a list to show that they are actually the same reactions
    print '\tFinding equivalent reactions ...'
    # Compounds with the same clean names
    equivalent_rxns = []
    if len(non_unique_name_rxns) > 0:
        # A dictionary with ModelSEED ids as keys and clean equations as values
        equiv_dict = {}
        for n in non_unique_name_rxns:
            for sid in ModelSEED_rxns_by_name[n]:
                if ModelSEED_rxns[sid]['equation'] != None:
                    equiv_dict[sid] = ModelSEED_rxns[sid]['equation']
        unique_equations = list(set(equiv_dict.values())) 
        for ue in unique_equations:
            eq_rxns = [sid for sid in equiv_dict.keys() if equiv_dict[sid] == ue]
            if len(eq_rxns) >  1:
                equivalent_rxns.append(eq_rxns)
   
    # Compounds with the same clean names
    non_unique_clean_name_rxns = [n for n in ModelSEED_rxns_by_clean_name.keys() if len(ModelSEED_rxns_by_clean_name[n]) > 1]
    if len(non_unique_clean_name_rxns) > 0:
        # A dictionary with ModelSEED ids as keys and equations as values
        equiv_dict = {}
        for n in non_unique_clean_name_rxns:
            for sid in ModelSEED_rxns_by_clean_name[n]:
                if ModelSEED_rxns[sid]['equation'] != None:
                    equiv_dict[sid] = ModelSEED_rxns[sid]['equation']
        unique_equations = list(set(equiv_dict.values())) 
        for ue in unique_equations:
            eq_rxns = [sid for sid in equiv_dict.keys() if equiv_dict[sid] == ue]
            if len(eq_rxns) > 1 and eq_rxns not in equivalent_rxns:
                equivalent_rxns.append(eq_rxns)

    print 'equivalent_rxns = ',equivalent_rxns

    with open('ModelSEED_rxns_website.py','a') as f:
        f.write('\nequivalent_rxns = [\n') 
        if len(equivalent_rxns) > 0:
            for eq_rxns in equivalent_rxns:
                if equivalent_rxns.index(eq_rxns) == len(equivalent_rxns) - 1:
                    f.write(repr(eq_rxns) + '\n') 
                else:
                    f.write(repr(eq_rxns) + ',\n') 
        f.write(']\n') 


    print 'Reactions were written into a python dictionary stored in ModelSEED_rxns_website.py\n'   

def read_name_aliases(aliases_tsv_file, cpds_results_file,rxns_results_file):
    """
    Read the tsv file containing the name aliases and writes the results into a python dictionary 

      aliases_tsv_file: Input tsv file name
     cpds_results_file: File name containing the formatted results for compounds
     rxns_results_file: File name containing the formatted results for reactions
    """

    print '\nParsing name aliases information ...'
    sys.stdout.flush()

    with open(cpds_results_file,'w') as f:
        f.write('cpds_by_name = {}\n')
        f.write('cpds_by_clean_name = {}\n')
        f.write('cpds_by_clean_name_underline = {}\n')
    with open(rxns_results_file,'w') as f:
        f.write('rxns_by_name = {}\n')
        f.write('rxns_by_clean_name = {}\n')
        f.write('rxns_by_clean_name_underline = {}\n')
   
    # Dictionaries with kesy being Kegg ids and values being ModelSEED ids  
    cpds_by_name = {}
    cpds_by_clean_name = {}           # Uses remove_non_alphanumeric with replace_with_underline = False
    cpds_by_clean_name_underline = {} # Uses remove_non_alphanumeric with replace_with_underline = True
    rxns_by_name = {}                 
    rxns_by_clean_name = {}           # Uses remove_non_alphanumeric with replace_with_underline = False
    rxns_by_clean_name_underline = {} # Uses remove_non_alphanumeric with replace_with_underline = True

    # List of cpds and rxns (by ModelSEED ids)
    cpds = []
    rxns = []

    with open(aliases_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            elif line[1] != 'null':
                # If there are multiple compounds convert it to a list
                if '|' in line[1]:
                    p = re.compile('(?:(?<=\|)|(?<=^)).*?(?=$|\|)') 
                    ModelSEED_ids = re.findall(p,line[1]) 
                else:
                    ModelSEED_ids = [line[1]]

                #--- cpds and rxns by name ---
                if 'cpd' in line[1]:
                    # Sometimes a name is used for both a compound and a reaction. In such ases we need to remove the
                    # rxn ids from the list
                    cpd_ids = [sid for sid in ModelSEED_ids if 'cpd' in sid]
                    if len(cpd_ids) == 0:
                        raise userError('cpd_ids is empty for ' + line[1])
                    cpds_by_name[line[0]] = ModelSEED_ids 
                    clean_name = remove_non_alphanumeric(line[0], replace_with_underline = False).lower()
                    clean_name_underline = remove_non_alphanumeric(line[0], replace_with_underline = True).lower()
                    with open(cpds_results_file,'a') as f:
                        if "'" not in line[0]:
                            f.write("cpds_by_name['" + line[0] + "'] = " + str(ModelSEED_ids) + "\n")
                        elif '"' not in line[0]:
                            f.write('cpds_by_name["' + line[0] + '"] = ' + str(ModelSEED_ids) + "\n")
                        else: # If both ' and " are in the string
                            f.write('cpds_by_name["""' + line[0] + '"""] = ' + str(ModelSEED_ids) + "\n")
                        f.write("cpds_by_clean_name['" + clean_name + "'] = " + str(ModelSEED_ids) + "\n")
                        f.write("cpds_by_clean_name_underline['" + clean_name_underline + "'] = " + str(ModelSEED_ids) + "\n")
                        cpds += cpd_ids 
                
                
                # Don't use elif here becuae somes the same name is assigned to both a reaction and a compound        
                if 'rxn' in line[1]:
                    # Sometimes a name is used for both a compound and a reaction. In such ases we need to remove the
                    # cpd ids from the list
                    rxn_ids = [sid for sid in ModelSEED_ids if 'rxn' in sid]
                    if len(rxn_ids) == 0:
                        raise userError('rxn_ids is empty for ' + line[1] + ', ModelSEED_ids = ' + str(ModelSEED_ids))
                    rxns_by_name[line[0]] = ModelSEED_ids 
                    clean_name = remove_non_alphanumeric(line[0], replace_with_underline = False).lower() 
                    clean_name_underline = remove_non_alphanumeric(line[0], replace_with_underline = True).lower() 

                    with open(rxns_results_file,'a') as f:
                        if "'" not in line[0]:
                            f.write("rxns_by_name['" + line[0] + "'] = " + str(ModelSEED_ids) + "\n")
                        elif '"' not in line[0]:
                            f.write('rxns_by_name["' + line[0] + '"] = ' + str(ModelSEED_ids) + "\n")
                        else: # If both ' and " are in the string
                            f.write('rxns_by_name["""' + line[0] + '"""] = ' + str(ModelSEED_ids) + "\n")
                        f.write("rxns_by_clean_name['" + clean_name + "'] = " + str(ModelSEED_ids) + "\n")
                        f.write("rxns_by_clean_name_underline['" + clean_name_underline + "'] = " + str(ModelSEED_ids) + "\n")
                        rxns += rxn_ids

                if 'cpd' not in line[1] and 'rxn' not in line[1]:
                    raise userError('Unknown entry! Expected cpd or rxn but got ' + line[1] + ' in line ' + str(counter) + '\n')

    #-- The Kegg id of each compound --
    with open(cpds_results_file,'a') as f:
        f.write('\ncpds_name_aliases = {}\n')
        for cpd in cpds:
             name = [n for n in cpds_by_name.keys() if cpd in cpds_by_name[n]]
             if len(name) > 0:
                 f.write("cpds_name_aliases['" + str(cpd) + "'] = " + str(name) + '\n')
             elif len(name) == 0:
                 raise userError('No names found for ' + cpd)


    with open(rxns_results_file,'a') as f:
        f.write('\nrxns_name_aliases = {}\n')
        for rxn in rxns:
             name = [n for n in rxns_by_name.keys() if rxn in rxns_by_name[n]]
             if len(name) > 0:
                 f.write("rxns_name_aliases['" + str(rxn) + "'] = " + str(name) + '\n')
             elif len(name) == 0:
                 raise userError('No names found for ' + rxn)

    print '\tThe results were written into {} and {}\n'.format(cpds_results_file,rxns_results_file)

def read_kegg_aliases(aliases_tsv_file, cpds_results_file,rxns_results_file):
    """
    Read the tsv file containing the kegg aliases and writes the results into a python dictionary 

    aliases_tsv_file: Input tsv file name
     cpds_results_file: File name containing the formatted results for compounds
     rxns_results_file: File name containing the formatted results for reactions
    """

    print '\nParsing KEGG aliases information ...'
    sys.stdout.flush()

    with open(cpds_results_file,'w') as f:
        f.write('cpds_by_KEGG_id = {}\n')
    with open(rxns_results_file,'w') as f:
        f.write('rxns_by_KEGG_id = {}\n')
   
    # Dictionaries with kesy being Kegg ids and values being ModelSEED ids  
    cpds_by_KEGG_id = {}
    rxns_by_KEGG_id = {}

    # List of cpds and rxns (by ModelSEED ids)
    cpds = []
    rxns = []

    with open(aliases_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            elif line[1] != 'null':
                # If there are multiple compounds convert it to a list
                if '|' in line[1]:
                    p = re.compile('(?:(?<=\|)|(?<=^)).*?(?=$|\|)') 
                    ModelSEED_id = re.findall(p,line[1]) 
                else:
                    ModelSEED_id = line[1]

                if 'cpd' in line[1]:
                    cpds_by_KEGG_id[line[0]] = ModelSEED_id 
                    with open(cpds_results_file,'a') as f:
                        if type(ModelSEED_id) is list:
                            f.write("cpds_by_KEGG_id['" + line[0] + "'] = " + str(ModelSEED_id) + "\n")
                            cpds += ModelSEED_id 
                        elif type(ModelSEED_id) is str:
                            f.write("cpds_by_KEGG_id['" + line[0] + "'] = '" + str(ModelSEED_id) + "'\n")
                            cpds.append(ModelSEED_id)
                        
                elif 'rxn' in line[1]:
                    rxns_by_KEGG_id[line[0]] = ModelSEED_id 
                    with open(rxns_results_file,'a') as f:
                        if type(ModelSEED_id) is list:
                            f.write("rxns_by_KEGG_id['" + line[0] + "'] = " + str(ModelSEED_id) + "\n")
                            rxns += ModelSEED_id
                        if type(ModelSEED_id) is str:
                            f.write("rxns_by_KEGG_id['" + line[0] + "'] = '" + str(ModelSEED_id) + "'\n")
                            rxns.append(ModelSEED_id)
                else:
                    raise userError('Unknown entry! Expected cpd or rxn but got ' + line[1] + ' in line ' + str(counter) + '\n')

    #-- The Kegg id of each compound --
    with open(cpds_results_file,'a') as f:
        f.write('\ncpds_KEGG_ids = {}\n')
        for cpd in cpds:
             KEGG_ids = [c for c in cpds_by_KEGG_id.keys() if cpds_by_KEGG_id[c] == cpd or cpd in cpds_by_KEGG_id[c]]
             if len(KEGG_ids) > 0:
                 f.write("cpds_KEGG_ids['" + str(cpd) + "'] = " + str(KEGG_ids) + '\n')
             elif len(KEGG_ids) == 0:
                 raise userError('No KEGG id found for ' + cpd)


    with open(rxns_results_file,'a') as f:
        f.write('\nrxns_KEGG_ids = {}\n')
        for rxn in rxns:
             KEGG_ids = [r for r in rxns_by_KEGG_id.keys() if rxns_by_KEGG_id[r] == rxn or rxn in rxns_by_KEGG_id[r]]
             if len(KEGG_ids) > 0:
                 f.write("rxns_KEGG_ids['" + str(rxn) + "'] = " + str(KEGG_ids) + '\n')
             elif len(KEGG_ids) == 0:
                 raise userError('No KEGG id found for ' + rxn)

    print '\tThe results were written into {} and {}\n'.format(cpds_results_file,rxns_results_file)

def read_ECnumber_aliases(aliases_tsv_file, rxns_results_file):
    """
    Read the tsv file containing the Enzyme class aliases and writes the results into a python dictionary 

    aliases_tsv_file: Input tsv file name
     rxns_results_file: File name containing the formatted results for EC numbers of reactions
    """

    print '\nParsing EC number aliases information ...'
    sys.stdout.flush()

    with open(rxns_results_file,'w') as f:
        f.write('rxns_by_EC_number = {}\n')
   
    # Dictionaries with kesy being Kegg ids and values being ModelSEED ids  
    rxns_by_EC_number = {}

    # List of cpds and rxns (by ModelSEED ids)
    rxns = []

    with open(aliases_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            elif line[1] != 'null':
                # If there are multiple EC numbers convert it to a list
                if '|' in line[1]:
                    p = re.compile('(?:(?<=\|)|(?<=^)).*?(?=$|\|)') 
                    ModelSEED_id = re.findall(p,line[1]) 
                else:
                    ModelSEED_id = [line[1]]

                if 'rxn' in line[1]:
                    if line[0] in rxns_by_EC_number.keys():
                        rxns_by_EC_number[line[0]] += ModelSEED_id 
                    else:
                        rxns_by_EC_number[line[0]] = ModelSEED_id 
                    with open(rxns_results_file,'a') as f:
                        if type(ModelSEED_id) is list:
                            f.write("rxns_by_EC_number['" + line[0] + "'] = " + str(ModelSEED_id) + "\n")
                            rxns += ModelSEED_id
                        if type(ModelSEED_id) is str:
                            f.write("rxns_by_EC_number['" + line[0] + "'] = '" + str(ModelSEED_id) + "'\n")
                            rxns.append(ModelSEED_id)
                else:
                    raise userError('Unknown entry! Expected cpd or rxn but got ' + line[1] + ' in line ' + str(counter) + '\n')

    #-- EC number of the enzymes coding for each reaction --
    with open(rxns_results_file,'a') as f:
        f.write('\nrxns_EC_numbers = {}\n')
        for rxn in rxns:
             EC_number = [ec for ec in rxns_by_EC_number.keys() if rxns_by_EC_number[ec] == rxn or rxn in rxns_by_EC_number[ec]]
             if len(EC_number) > 0:
                 f.write("rxns_EC_numbers['" + str(rxn) + "'] = " + str(EC_number) + '\n')
             elif len(EC_number) == 0:
                 raise userError('No EC number found for ' + rxn)

    print '\tThe results were written into {}\n'.format(rxns_results_file)


def read_BiGG_aliases(aliases_tsv_file, cpds_results_file,rxns_results_file):
    """
    Read the tsv file containing the kegg aliases and writes the results into a python dictionary 

    aliases_tsv_file: Input tsv file name
     cpds_results_file: File name containing the formatted results for compounds
     rxns_results_file: File name containing the formatted results for reactions
    """

    print '\nParsing BiGG aliases information ...'
    sys.stdout.flush()

    with open(cpds_results_file,'w') as f:
        f.write('cpds_by_BiGG_id = {}\n')
        f.write('cpds_by_clean_BiGG_id = {}\n')
    with open(rxns_results_file,'w') as f:
        f.write('rxns_by_BiGG_id = {}\n')
        f.write('rxns_by_clean_BiGG_id = {}\n')
   
    # Dictionaries with kesy being Kegg ids and values being ModelSEED ids  
    cpds_by_BiGG_id = {}
    cpds_by_clean_BiGG_id = {}
    rxns_by_BiGG_id = {}
    rxns_by_clean_BiGG_id = {}

    # List of cpds and rxns (by ModelSEED ids)
    cpds = []
    rxns = []

    with open(aliases_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            elif line[1] != 'null':
                # If there are multiple compounds convert it to a list
                if '|' in line[1]:
                    p = re.compile('(?:(?<=\|)|(?<=^)).*?(?=$|\|)') 
                    ModelSEED_ids = re.findall(p,line[1]) 
                else:
                    ModelSEED_ids = [line[1]]

                if 'cpd' in line[1]:
                    cpds_by_BiGG_id[line[0]] = ModelSEED_ids 
                    cpds_by_clean_BiGG_id[remove_non_alphanumeric(line[0], replace_with_underline = True).lower()] = ModelSEED_ids 
                    with open(cpds_results_file,'a') as f:
                        if "'" not in line[0]:
                            f.write("cpds_by_BiGG_id['" + line[0] + "'] = " + str(ModelSEED_ids) + "\n")
                        elif '"' not in line[0]:
                            f.write('cpds_by_BiGG_id["' + line[0] + '"] = ' + str(ModelSEED_ids) + "\n")
                        else: # If both ' and " are in the string
                            f.write('cpds_by_BiGG_id["""' + line[0] + '"""] = ' + str(ModelSEED_ids) + "\n")
                        f.write("cpds_by_clean_BiGG_id['" + remove_non_alphanumeric(line[0], replace_with_underline = True).lower() + "'] = " + str(ModelSEED_ids) + "\n")
                        cpds += ModelSEED_ids 
                        
                if 'rxn' in line[1]:
                    rxns_by_BiGG_id[line[0]] = ModelSEED_ids 
                    rxns_by_clean_BiGG_id[remove_non_alphanumeric(line[0], replace_with_underline = True).lower()] = ModelSEED_ids 
                    with open(rxns_results_file,'a') as f:
                        if "'" not in line[0]:
                            f.write("rxns_by_BiGG_id['" + line[0] + "'] = " + str(ModelSEED_ids) + "\n")
                        elif '"' not in line[0]:
                            f.write('rxns_by_BiGG_id["' + line[0] + '"] = ' + str(ModelSEED_ids) + "\n")
                        else: # If both ' and " are in the string
                            f.write('rxns_by_BiGG_id["""' + line[0] + '"""] = ' + str(ModelSEED_ids) + "\n")
                        f.write("rxns_by_clean_BiGG_id['" + remove_non_alphanumeric(line[0], replace_with_underline = True).lower() + "'] = " + str(ModelSEED_ids) + "\n")
                        rxns += ModelSEED_ids 

                if 'cpd' not in line[1] and 'rxn' not in line[1]:
                    raise userError('Unknown entry! Expected cpd or rxn but got ' + line[1] + ' in line ' + str(counter) + '\n')

    #-- The BiGG id of each compound --
    with open(cpds_results_file,'a') as f:
        f.write('\ncpds_BiGG_ids = {}\n')
        for cpd in cpds:
             BiGG_ids = [b for b in cpds_by_BiGG_id.keys() if cpd in cpds_by_BiGG_id[b]]
             if len(BiGG_ids) > 0:
                 f.write("cpds_BiGG_ids['" + str(cpd) + "'] = " + str(BiGG_ids) + '\n')
             elif len(BiGG_ids) == 0:
                 raise userError('No BiGG id found for ' + cpd)

    with open(rxns_results_file,'a') as f:
        f.write('\nrxns_BiGG_ids = {}\n')
        for rxn in rxns:
             BiGG_ids = [b for b in rxns_by_BiGG_id.keys() if rxn in rxns_by_BiGG_id[b]]
             if len(BiGG_ids) > 0:
                 f.write("rxns_BiGG_ids['" + str(rxn) + "'] = " + str(BiGG_ids) + '\n')
             elif len(BiGG_ids) == 0:
                 raise userError('No BiGG id found for ' + rxn)

    print '\tThe results were written into {} and {}\n'.format(cpds_results_file,rxns_results_file)

def read_cpds_master(cpds_master_file, cpds_results_file):
    """
    Read the master files for compounds 

        cpds_master_file: master file for compounds 
    cpds_results_file: File name containing the formatted results for compounds
    """

    # Import name aliases, kEGG name(s) and BiGG names
    from ModelSEED_cpds_name_aliases import cpds_name_aliases 
    from ModelSEED_cpds_KEGG_aliases import cpds_KEGG_ids
    from ModelSEED_cpds_BiGG_aliases import cpds_BiGG_ids

    print '\nParsing compounds master file ...'
    sys.stdout.flush()

    with open(cpds_results_file,'w') as f:
        f.write('cpds_master = {}\n')
   
    # List of cpds and rxns (by ModelSEED ids)
    cpds_master = []

    with open(cpds_master_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                # headers = ['id','abbreviation','name','formula','mass','source','structure','charge','is_core','is_obsolete','linked_compound','is_cofactor','deltag','deltagerr','pka','pkb','abstract_compound','comprised_of','aliases'] 
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                cpd_info = {}

                # ModelSeed id
                ModelSEED_id = line[headers.index('id')]
                cpd_info['ModelSEED_id'] = ModelSEED_id 

                # Abbreviation
                cpd_info['abbreviation'] = line[headers.index('abbreviation')]

                # Primary name
                cpd_info['name'] = line[headers.index('name')] 

                # formula
                cpd_info['formula'] = line[headers.index('formula')]
                if cpd_info['formula'] == 'null':
                    cpd_info['formula'] = None

                # Name aliases
                if ModelSEED_id in cpds_name_aliases.keys():
                    cpd_info['name_aliases'] = cpds_name_aliases[ModelSEED_id]
                else:
                    cpd_info['name_aliases'] = [] 

                # KEEG id(s)
                if ModelSEED_id in cpds_KEGG_ids.keys():
                    cpd_info['KEGG_id'] = cpds_KEGG_ids[ModelSEED_id]
                else:
                    cpd_info['KEGG_id'] = []

                # BiGG name
                if ModelSEED_id in cpds_BiGG_ids.keys():
                    cpd_info['BiGG_id'] = cpds_BiGG_ids[ModelSEED_id]
                else:
                    cpd_info['BiGG_id'] = []

                # structure
                cpd_info['structure'] = line[headers.index('structure')]
                if cpd_info['structure'] == 'null':
                    cpd_info['structure'] = None

                # mass
                mass_str = line[headers.index('mass')]
                if mass_str != 'null':
                    cpd_info['mass'] = float(mass_str)
                else:
                    cpd_info['mass'] = None

                # charge
                charge_str = line[headers.index('charge')]
                if charge_str != 'null':
                    cpd_info['charge'] = float(charge_str)
                else:
                    cpd_info['charge'] = None

                # Linked compounds
                cpd_info['linked_compounds'] = line[headers.index('linked_compound')]
                if cpd_info['linked_compounds'] == 'null':
                    cpd_info['linked_compounds'] = None

                # is_cofactor
                is_cofactor_str = line[headers.index('is_cofactor')]
                if is_cofactor_str != 'null':
                    if is_cofactor_str == '1':
                        cpd_info['is_cofactor'] = True
                    elif is_cofactor_str == '0':
                        cpd_info['is_cofactor'] = False
                    else:
                        raise userError('Invalid is_cofactor value: ' + is_cofactor_str)
                else:
                    cpd_info['is_cofactor'] = None

                # deltaG
                deltaG_str = line[headers.index('deltag')]
                if deltaG_str != 'null' and deltaG_str != '10000000':
                    cpd_info['deltaG'] = float(deltaG_str)
                else:
                    cpd_info['deltaG'] = None

                # deltaG_uncertainty
                deltaG_uncertainty_str = line[headers.index('deltagerr')]
                if deltaG_uncertainty_str != 'null' and deltaG_uncertainty_str != '10000000':
                    cpd_info['deltaG_uncertainty'] = float(deltaG_uncertainty_str)
                else:
                    cpd_info['deltaG_uncertainty'] = None

                # pKa
                pKa_str = line[headers.index('pka')]
                if pKa_str != 'null':
                    cpd_info['pKa'] = pKa_str
                else:
                    cpd_info['pKa'] = None

                # pKb
                pKb_str = line[headers.index('pkb')]
                if pKb_str != 'null':
                    cpd_info['pKb'] = pKb_str
                else:
                    cpd_info['pKb'] = None

                with open(cpds_results_file,'a') as f:
                    f.write("cpds_master['" + ModelSEED_id + "'] = " + str(cpd_info) + '\n')
   
def read_comparts(comparts_tsv_file, comparts_results_file):
    """
    Read the the compartments information 
         coparts_tsv_file: name of the tsv file containing the compartments information
    comparts_results_file: name of the file containing the formatted results for compartments
    """
    print '\nGetting compartment indices info ...'
    sys.stdout.flush()

    compartments = {}
    compartments_by_index = {}

    with open(comparts_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                # index
                index = line[headers.index('index')]

                # id
                id = line[headers.index('id')]

                # name
                name = line[headers.index('name')]

                compartments[id] = name
                compartments_by_index[index] = {'id':id,'name':name}

    # Compartments by id 
    with open(comparts_results_file,'a') as f:
        f.write('ModelSEED_compartments = {}\n')
        for compart_id in compartments.keys():
            f.write("ModelSEED_compartments{'" + compart_id + "'] = '" + compartments[compart_id] +"'\n")

    # Compartments by iindex
    with open(comparts_results_file,'w') as f:
        f.write('\nModelSEED_compartments_by_index = {}\n')
        for index in compartments_by_index.keys():
            f.write("ModelSEED_compartments_by_index['" + index + "'] = {'id':'" + compartments_by_index[index]['id'] + "', 'name':'" + compartments_by_index[index]['name'] + "'}\n")


    print 'Results were written into {}\n'.format(comparts_results_file) 

def read_rxn_pathways(pathways_tsv_file, pathways_results_file):
    """
    Read the the pathways for KEGG reactions 
        pathways_tsv_file: name of the tsv file containing the pathways information
    pathways_results_file: name of the file containing the formatted results for pathways
    """
    print '\nParsing reaction pathways ...'
    sys.stdout.flush()

    KEGG_rxns_pathways = {}
    KEGG_rxns_by_pathways = {}

    with open(pathways_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                # Pathway
                pathway = line[headers.index('SCENARIO')]

                # KEGG id
                KEGG_id = line[headers.index('REACTION')]

                if KEGG_id in KEGG_rxns_pathways.keys():
                    KEGG_rxns_pathways[KEGG_id] += [pathway]
                else:
                    KEGG_rxns_pathways[KEGG_id] = [pathway]

                if pathway in KEGG_rxns_by_pathways.keys():
                    KEGG_rxns_by_pathways[pathway] += [KEGG_id]
                else:
                    KEGG_rxns_by_pathways[pathway] = [KEGG_id]

    # KEGG reactions pathways 
    with open(pathways_results_file,'w') as f:
        f.write('KEGG_rxns_pathways = {}\n')
        for kegg_id in KEGG_rxns_pathways.keys():
            f.write("KEGG_rxns_pathways['" + kegg_id + "'] = " + str(KEGG_rxns_pathways[kegg_id]) +"\n")

        f.write("\nKEGG_rxns_by_pathways = {}\n")
        for pathway in KEGG_rxns_by_pathways.keys():
            if "'" not in pathway:
                f.write("KEGG_rxns_by_pathways['" + pathway + "'] = " + str(KEGG_rxns_by_pathways[pathway]) + "\n")
            elif '"' not in pathway:
                f.write('KEGG_rxns_by_pathways["' + pathway + '"] = ' + str(KEGG_rxns_by_pathways[pathway]) + "\n")
            

    print 'Results were written into {}\n'.format(pathways_results_file) 


def read_rxns_master(rxns_master_file, rxns_default_file, rxns_results_file):
    """
    Read the master files for reactions 

        rxns_master_file: master file for reactions
       rxns_default_file: default file for reactions (needed to extract reaction abbreviations) 
            coparts_file: name of the file containing the compartments information
    rxns_results_file: File name containing the formatted results for reactions
    """

    # Import name aliases, kEGG name(s) and BiGG names
    from ModelSEED_rxns_name_aliases import rxns_name_aliases 
    from ModelSEED_rxns_KEGG_aliases import rxns_KEGG_ids
    from ModelSEED_rxns_BiGG_aliases import rxns_BiGG_ids 
    from ModelSEED_rxns_EC_number_aliases import rxns_EC_numbers 
    from KEGG_rxns_pathways import KEGG_rxns_pathways
    from ModelSEED_compartments import ModelSEED_compartments_by_index

    #--- Get reaction abbreviations from the default file ----
    print '\nGetting reaction abbreviations from reactions default file ...'
    sys.stdout.flush()

    # List of rxns and rxns (by ModelSEED ids)
    rxns_abbr = {}

    with open(rxns_default_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                # ModelSeed id
                ModelSEED_id = line[headers.index('id')]

                # name 
                name = line[headers.index('name')]

                # Abbreviation 
                abbr = line[headers.index('abbreviation')]
                # If the abbreviation is the same as name change it to the ModelSEED_id
                if abbr == name:
                     abbr = ModelSEED_id

                rxns_abbr[ModelSEED_id] = abbr

    #--- Parsing reactions master file from the default file ----
    print '\nParsing reactions master file ...'
    with open(rxns_results_file,'w') as f:
        f.write('rxns_master = {}\n')
   
    # List of rxns and rxns (by ModelSEED ids)
    rxns_master = []

    with open(rxns_master_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                # headers = ['id','abbreviation','name','code','stoichiometry','is_transport','equation','definition','reversibility','direction','abstract_reaction','pathways','aliases','ec_numbers','deltag','deltagerr','compound_ids','status','is_obsolete','linked_reaction'] 
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                rxn_info = {}

                # ModelSeed id
                ModelSEED_id = line[headers.index('id')]
                rxn_info['ModelSEED_id'] = ModelSEED_id

                # Abbreviation
                if ModelSEED_id in rxns_abbr.keys():
                    rxn_info['abbreviation'] = rxns_abbr[ModelSEED_id] 
                else:
                    rxn_info['abbreviation'] = ModelSEED_id

                # Primary name
                rxn_info['name'] = line[headers.index('name')] 

                # Name aliases
                if ModelSEED_id in rxns_name_aliases.keys():
                    rxn_info['name_aliases'] = rxns_name_aliases[ModelSEED_id]
                else:
                    rxn_info['name_aliases'] = []
    
                # KEEG id(s)
                if ModelSEED_id in rxns_KEGG_ids.keys():
                    rxn_info['KEGG_id'] = rxns_KEGG_ids[ModelSEED_id]
                else:
                    rxn_info['KEGG_id'] = []
    
                # BiGG name
                if ModelSEED_id in rxns_BiGG_ids.keys():
                    rxn_info['BiGG_id'] = rxns_BiGG_ids[ModelSEED_id]
                else:
                    rxn_info['BiGG_id'] = []

                # deltaG
                deltaG_ModelSEED_str = line[headers.index('deltag')]
                if deltaG_ModelSEED_str != 'null' and deltaG_ModelSEED_str != '10000000':
                    deltaG_ModelSEED = float(deltaG_ModelSEED_str)
                    rxn_info['deltaG_ModelSEED'] = deltaG_ModelSEED 
                else:
                    deltaG_ModelSEED = None
                    rxn_info['deltaG_ModelSEED'] = deltaG_ModelSEED

                # deltaG_uncertainty
                deltaG_uncertainty_ModelSEED_str = line[headers.index('deltagerr')]
                if deltaG_uncertainty_ModelSEED_str != 'null' and deltaG_uncertainty_ModelSEED_str != '10000000':
                    deltaG_uncertainty_ModelSEED = float(deltaG_uncertainty_ModelSEED_str)
                    rxn_info['deltaG_uncertainty_ModelSEED'] = deltaG_uncertainty_ModelSEED 
                else:
                    deltaG_uncertainty_ModelSEED = None
                    rxn_info['deltaG_uncertainty_ModelSEED'] = deltaG_uncertainty_ModelSEED 

                # deltaG_min and deltaG_max
                if deltaG_ModelSEED != None and deltaG_uncertainty_ModelSEED != None:
                    rxn_info['deltaG_min_ModelSEED'] = min(deltaG_ModelSEED - deltaG_uncertainty_ModelSEED, deltaG_ModelSEED + deltaG_uncertainty_ModelSEED) 
                    rxn_info['deltaG_max_ModelSEED'] = max(deltaG_ModelSEED - deltaG_uncertainty_ModelSEED, deltaG_ModelSEED + deltaG_uncertainty_ModelSEED) 

                    # Probability that dG is less than zero (Equation 3 of PMID: 19783351)
                    probability_dG_lessThanZero = 0.5*(1 + sps.erf(-deltaG_ModelSEED/(np.power(2,0.5)*deltaG_uncertainty_ModelSEED)))   
                    if probability_dG_lessThanZero >= 0 and probability_dG_lessThanZero <= 1:
                        rxn_info['probability_dG_lessThanZero'] = probability_dG_lessThanZero
                    else:
                        raise userError('probability_dG_lessThanZero for rxn ' + ModelSEED_id + ' not between zero and one: ' + str(probability_dG_lessThanZero))
                else:
                    rxn_info['deltaG_min_ModelSEED'] = None 
                    rxn_info['deltaG_max_ModelSEED'] = None 
                    rxn_info['probability_dG_lessThanZero'] = None 

                # Reaction reversability (this seems to be consistent with deltaG value.
                rev_str = line[headers.index('reversibility')]
                if rev_str == '=':
                    rxn_info['reversibility_ModelSEED_deltaG'] = 'reversible'
                elif rev_str == '>':
                    rxn_info['reversibility_ModelSEED_deltaG'] = 'irreversible_forward'
                elif rev_str == '<':
                    rxn_info['reversibility_ModelSEED_deltaG'] = 'irreversible_backward'
                elif rev_str == '?': # If unknown set it to reversible
                    rxn_info['reversibility_ModelSEED_deltaG'] = 'reversible'
                else:
                    raise userError('Unknown reversibility for reaction ' + ModelSEED_id + ': ' + rev_str)

                # Reaction direction (inferred from a combination of deltaG values and heuristic rules 
                # given in PMID: 19555510). It can be different from reversibility 
                rev_str = line[headers.index('direction')]
                if rev_str == '=':
                    rxn_info['reversibility_ModelSEED_curated_master'] = 'reversible'
                elif rev_str == '>':
                    rxn_info['reversibility_ModelSEED_curated_master'] = 'irreversible_forward'
                elif rev_str == '<':
                    rxn_info['reversibility_ModelSEED_curated_master'] = 'irreversible_backward'
                elif rev_str == '?': # If unknown set it to reversible
                    rxn_info['reversibility_ModelSEED_curated_master'] = 'reversible'
                else:
                    raise userError('Unknown direction_master for reaction ' + ModelSEED_id + ': ' + rev_str)

                # Reaction equation baed ModelSEED ids for compounds
                equation_by_cpd_ids = line[headers.index('equation')]
                if rxn_info['reversibility_ModelSEED_curated_master'] == 'reversible':
                    rxn_info['equation_by_cpd_ids'] = re.sub('<=>','<==>',equation_by_cpd_ids)
                elif rxn_info['reversibility_ModelSEED_curated_master'] == 'irreversible_forward':
                    rxn_info['equUation_by_cpd_ids'] = re.sub('<=>','-->',equation_by_cpd_ids)
                elif rxn_info['reversibility_ModelSEED_curated_master'] == 'irreversible_backward':
                    rxn_info['equation_by_cpd_ids'] = re.sub('<=>','<--',equation_by_cpd_ids)

                # Reaciton equation based names for compounds 
                equation_by_cpd_names = line[headers.index('definition')]
                if rxn_info['reversibility_ModelSEED_curated_master'] == 'reversible':
                    rxn_info['equation_by_cpd_names'] = re.sub('<=>','<==>',equation_by_cpd_names)
                elif rxn_info['reversibility_ModelSEED_curated_master'] == 'irreversible_forward':
                    rxn_info['equation_by_cpd_names'] = re.sub('<=>','-->',equation_by_cpd_names)
                elif rxn_info['reversibility_ModelSEED_curated_master'] == 'irreversible_backward':
                    rxn_info['equation_by_cpd_names'] = re.sub('<=>','<--',equation_by_cpd_names)

                # Stoichiometry
                stoic_str = line[headers.index('stoichiometry')]
                if stoic_str != 'null':
                    stoic = {}
                    # Parse stoic_str
                    p1 = re.compile('(?:(?<=\;)|(?<=^)).*?(?=$|\;)')
                    entries = re.findall(p1,stoic_str)
                    for entry in entries:
                        p2 = re.compile('(?:(?<=\:)|(?<=^)).*?(?=$|\:)')
                        parsed_entry = re.findall(p2,entry)
                        #print ModelSEED_id,"\t",entry,"\t",parsed_entry
                        stoic_coeff = float(parsed_entry[0])
                        cpd_id = parsed_entry[1] 
                        compart_index = parsed_entry[2] 
                        stoic[(cpd_id,ModelSEED_compartments_by_index[compart_index]['id'])] = stoic_coeff
                    rxn_info['stoichiometry'] = stoic
                else:
                    rxn_info['stoichiometry'] = None
                    
                # Status 
                status_str = line[headers.index('status')]
                # Status definitions
                status_def = {'OK':'balanced','MI':'mass_imbalance','CI':'charge_imbalance','HB':'balanced_by_H','EMPTY':'reactants_and_products_cancel_out','CPDFORMERROR':'no_formula_compounds'}
                # Parse those separated by |
                p1 = re.compile('(?:(?<=\|)|(?<=^)).*?(?=$|\|)')
                entries = re.findall(p1,status_str)
                # Extract the status at the begining before : e.g., parse MI:H:-3/O:-2 to MI and H:-3/O:-2
                status_details = {} # Keys are status abbreivations and values are their details e.g, MI and H:-3/O:-2 
                for entry in entries:
                    p2 = re.compile('(?<=^).*?(?=$|:)|(?<=:).*(?=$)')
                    stat_abbr_details = re.findall(p2,entry)
                    if len(stat_abbr_details) == 1:
                        stat_abbr = stat_abbr_details[0]
                        stat_det = ''
                    elif len(stat_abbr_details) == 2:
                        stat_abbr = stat_abbr_details[0]
                        stat_det = stat_abbr_details[1]
                    else: # If lenth is more than two
                        raise userError('length of the parsed status is more than two for rxn ' + ModelSEED_id + '  status = ' + status_srt)
                    status_details[stat_abbr] = stat_det
                if 'OK' in status_details.keys():
                    status = status_def['OK']
                else:
                    status = '|'.join([status_def[s] + ':(' + status_details[s] + ')' if s in ['MI','CI'] else status_def[s] for s in status_details.keys()])
                rxn_info['status'] = status
      
                # is_transport
                is_transport_str = line[headers.index('is_transport')]
                if is_transport_str != 'null':
                    if int(is_transport_str) == 1:
                        rxn_info['is_transport'] = True
                    elif int(is_transport_str) == 0:
                        rxn_info['is_transport'] = False
                    else:
                        raise userError('Invalid is_transport value: ' + is_transport_str)
                else:
                    rxn_info['is_transport'] = None

                # EC numbers 
                if ModelSEED_id in rxns_EC_numbers.keys():
                    rxn_info['EC_numbers'] = rxns_EC_numbers[ModelSEED_id]
                else:
                    rxn_info['EC_numbers'] = None
    
                # Reaciton pathways
                pathways = []
                if rxn_info['KEGG_id'] != None and rxn_info['KEGG_id'] in KEGG_rxns_pathways.keys():
                    for kid in rxn_info['KEGG_id']:
                        pathways += KEGG_rxns_pathways[kid]
                else:
                    pathways = None
                rxn_info['pathways'] = pathways
 
                with open(rxns_results_file,'a') as f:
                    f.write("rxns_master['" + ModelSEED_id + "'] = " + str(rxn_info) + "\n")

def read_template_rxns(rxns_tsv_file, rxns_results_file, results_var_name): 
    """
    Read the master files for reactions 

        rxns_tsv_file: tsv file containing reactions
    rxns_results_file: File name containing the formatted results for reactions
     results_var_name: Name of the variable containing the results
    """
    from ModelSEED_rxns_master import rxns_master as ModelSEED_rxns_master

    print '\nParsing template reactions file {} ...'.format(rxns_tsv_file)
    with open(rxns_results_file,'w') as f:
        f.write(results_var_name + ' = {}\n')
   
    with open(rxns_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                # headers = ['id','compartment','direction','gfdir','type','base_cost','forward_cost','reverse_cost','complexes'] 
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + aliases_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                rxn_info = {}

                # ModelSeed id
                ModelSEED_id = line[headers.index('id')]

                # Direction
                direction_str = line[headers.index('direction')]
                if direction_str == '>':
                    rxn_info['reversibility_ModelSEED_curated_template'] = 'irreversible_forward'
                elif direction_str == '<':
                    rxn_info['reversibility_ModelSEED_curated_template'] = 'irreversible_backward'
                elif direction_str == '=':
                    rxn_info['reversibility_ModelSEED_curated_template'] = 'reversible'
                elif direction_str == '?':
                    rxn_info['reversibility_ModelSEED_curated_template'] = 'unknown'
                else:
                    raise ValueError('Unknown direction for reaction ' + ModelSEED_id + ' (direction: ' + direction_str + ')\n')

                # Reaction type 
                rxn_info['template_rxn_type'] = line[headers.index('type')] 
                if rxn_info['template_rxn_type'] not in ['conditional','gapfilling','spontaneous','universal']:
                    raise ValueError('Invalid "type" for reaction ' + ModelSEED_id + ': type = ' + rxn_info['template_type'])

                # Base cost of adding this reaction to a model
                base_cost_str = line[headers.index('base_cost')] 
                rxn_info['base_cost'] = float(base_cost_str)  

                # The cost of adding this reaction to a model in the forward direction
                fwd_cost_str = line[headers.index('forward_cost')] 
                rxn_info['forward_cost'] = float(fwd_cost_str)  
 
                # The cost of adding this reaction to a model in the backward (reverse) direction
                bwd_cost_str = line[headers.index('reverse_cost')] 
                rxn_info['backward_cost'] = float(bwd_cost_str)  

                with open(rxns_results_file,'a') as f:
                    f.write(results_var_name + "['" + ModelSEED_id + "'] = " + str(rxn_info) + "\n")


def read_template_biomass(biomasses_tsv_file, biomass_cpds_tsv_file, biomass_results_file, results_var_name):
    """
    Read the biomass reaction from templates' biomass files

    INPUTS:
    -------
       biomasses_tsv_file: The input file for Biomasses.tsv
    biomass_cpds_tsv_file: The input file for BiomassCompounds.tsv
     biomass_results_file: Results file name   
         results_var_name: Name of the variable storing the results
    """
    #--- First extract information from Biomasses.tsv ---
    print '\nParsing template biomass reactions files {} ...'.format(biomasses_tsv_file)

    rxns_general_info = {}

    with open(biomasses_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                # headers = ['id','name','type','other','dna','rna','protein','lipid','cellwall','cofactor','energy'] 
                headers = line

            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of ' + biomasses_tsv_file +'. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                # id
                id = line[headers.index('id')]
                rxns_general_info[id] = {}

                # type
                rxns_general_info[id]['biomass_rxn_name'] = line[headers.index('name')]

                # type
                rxns_general_info[id]['biomass_rxn_type'] = line[headers.index('type')]

                # dna, rna, protein, lipid, cellwall, cofactor and energy  coefficient
                dna_coeff = float(line[headers.index('dna')])
                rna_coeff = float(line[headers.index('rna')])
                protein_coeff = float(line[headers.index('protein')])
                lipid_coeff = float(line[headers.index('lipid')])
                cellwall_coeff = float(line[headers.index('cellwall')])
                cofactor_coeff = float(line[headers.index('cofactor')])
                energy_coeff = float(line[headers.index('energy')])
                rxns_general_info[id]['macromolecules_stoic'] = {'dna':dna_coeff,'rna':rna_coeff, 'protein':protein_coeff, 'lipid':lipid_coeff, 'cellwall':cellwall_coeff, 'cofactor':cofactor_coeff, 'energy':energy_coeff}

    #--- Next, extract information from BiomassCompounds.tsv ---
    with open(biomass_results_file,'w') as f:
        f.write(results_var_name + ' = {}\n')

    # List of biomass reactions ids
    biomass_ids = [] 
   
    bio_rxns_info = {}

    with open(biomass_cpds_tsv_file,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read headers
                # headers = ['biomass_id','id','coefficient','coefficient_type','class','linked_compounds','compartment'] 
                headers = line

            elif len(line) != len(headers):
                raise userError('Line {} of {}. The number of items ({}) is not equal to the number of headers ({})\nline = {}'.format(counter,biomass_cpds_tsv_file,len(line),len(headers), line))
            else:
                # biomass id
                biomass_id = line[headers.index('biomass_id')]
                if biomass_id not in biomass_ids:
                    bio_rxns_info[biomass_id] = {} 
                    bio_rxns_info[biomass_id]['biomass_rxn_name'] = rxns_general_info[biomass_id]['biomass_rxn_name']
                    bio_rxns_info[biomass_id]['biomass_rxn_type'] = rxns_general_info[biomass_id]['biomass_rxn_type']
                    bio_rxns_info[biomass_id]['macromolecular_stoichiometry'] = rxns_general_info[biomass_id]['macromolecules_stoic'] 
                    bio_rxns_info[biomass_id]['stoichiometry'] = {}
                    bio_rxns_info[biomass_id]['compounds_info'] = {}
                    biomass_ids.append(biomass_id)

                # Compound id 
                cpd_id = line[headers.index('id')]

                # Stoichiometric coefficient
                cpd_cpeff = float(line[headers.index('coefficient')])
                bio_rxns_info[biomass_id]['stoichiometry'][cpd_id] = cpd_cpeff

                # Infomration about this compound
                bio_rxns_info[biomass_id]['compounds_info'][cpd_id] = {}

                # Stoichiometric coefficient type
                cpd_cpeff_type = line[headers.index('coefficient_type')]

                # Compound class
                cpd_class = line[headers.index('class')]

                # linked_compounds
                linked_cpds = line[headers.index('linked_compounds')]

                # compartment
                cpd_compart = line[headers.index('compartment')]

                # Store parameters
                bio_rxns_info[biomass_id]['compounds_info'][cpd_id]['stoic_coeff_tyype'] = cpd_cpeff_type
                bio_rxns_info[biomass_id]['compounds_info'][cpd_id]['cpd_class'] = cpd_class
                bio_rxns_info[biomass_id]['compounds_info'][cpd_id]['linked_cpds'] = linked_cpds
                bio_rxns_info[biomass_id]['compounds_info'][cpd_id]['compartment'] = cpd_compart

    with open(biomass_results_file,'a') as f:
        for biomass_id in bio_rxns_info.keys():
            f.write('\n' + results_var_name + "['" + bio_rxns_info[biomass_id]['biomass_rxn_name'] + "'] = {}\n")
            f.write('\n' + results_var_name + "['" + bio_rxns_info[biomass_id]['biomass_rxn_name'] + "']['biomass_rxn_type'] = " + bio_rxns_info[biomass_id]['biomass_rxn_type'] + "\n")
            f.write('\n' + results_var_name + "['" + bio_rxns_info[biomass_id]['biomass_rxn_name'] + "']['macromolecular_stoichiometry'] = " + str(bio_rxns_info[biomass_id]['macromolecular_stoichiometry']) + "\n")
            for cpd_id in bio_rxns_info[biomass_id]['stoichiometry'].keys():
                f.write('\n' + results_var_name + "['" + bio_rxns_info[biomass_id]['biomass_rxn_name'] + "']['stoichiometry']['" + cpd_id + "'] = " + str(bio_rxns_info[biomass_id]['stoichiometry'][cpd_id]) + "\n")
            for cpd_id in bio_rxns_info[biomass_id]['compounds_info'].keys():
                f.write('\n' + results_var_name + "['" + bio_rxns_info[biomass_id]['biomass_rxn_name'] + "']['compounds_info']['" + cpd_id + "'] = " + str(bio_rxns_info[biomass_id]['compounds_info'][cpd_id]) + "\n")


def write_cpds_formula_aliases(cpds_results_file):
    """
    Creates a dictionary to access compounds by their formula 

     cpds_results_file: File name containing the formatted results for compounds
    """
    print '\nCreating formula aliases for compounds ...'
    sys.stdout.flush()

    from ModelSEED_cpds_master import cpds_master

    with open(cpds_results_file,'w') as f:
        f.write('cpds_by_formula = {}\n')
   
    # Dictionaries with kesy being Kegg ids and values being ModelSEED ids  
    cpds_by_formula = {}

    for ModelSEED_id in [mid for mid in cpds_master.keys() if cpds_master[mid]['formula'] != None ]:
        formula = cpds_master[ModelSEED_id]['formula']
        if formula in cpds_by_formula.keys():
            cpds_by_formula[formula] += [ModelSEED_id]
        else:
            cpds_by_formula[formula] = [ModelSEED_id]


    with open(cpds_results_file,'a') as f:
        for formula in cpds_by_formula.keys(): 
            f.write("cpds_by_formula['" + formula + "'] = " + str(cpds_by_formula[formula]) + "\n")

#-------------------------------------------------------------------------
if __name__ == '__main__':
    #read_tsv_files_website(compounds_tsv_file = 'ModelSEED_cpds_website.tsv',reactions_tsv_file = 'ModelSEED_rxns_website.tsv')

    #read_name_aliases(aliases_tsv_file = 'github/Aliases/name.aliases', cpds_results_file = 'ModelSEED_cpds_name_aliases.py',rxns_results_file = 'ModelSEED_rxns_name_aliases.py')

    #read_kegg_aliases(aliases_tsv_file = 'github/Aliases/KEGG.aliases', cpds_results_file = 'ModelSEED_cpds_KEGG_aliases.py',rxns_results_file = 'ModelSEED_rxns_KEGG_aliases.py')

    #read_ECnumber_aliases(aliases_tsv_file = 'github/Aliases/Enzyme_Class.aliases', rxns_results_file = 'ModelSEED_rxns_EC_number_aliases.py')

    #read_BiGG_aliases(aliases_tsv_file = 'github/Aliases/BiGG.aliases', cpds_results_file = 'ModelSEED_cpds_BiGG_aliases.py',rxns_results_file = 'ModelSEED_rxns_BiGG_aliases.py')

    #read_cpds_master(cpds_master_file = 'github/Biochemistry/compounds.master.tsv', cpds_results_file = 'ModelSEED_cpds_master.py')

    #read_comparts(comparts_tsv_file = 'github/Templates/Human/Compartments.tsv', comparts_results_file = 'ModelSEED_compartments.py')

    #read_rxn_pathways(pathways_tsv_file = 'github/Pathways/HopeScenarios.txt', pathways_results_file = 'KEGG_rxns_pathways.py')

    #read_rxns_master(rxns_master_file = 'github/Biochemistry/reactions.master.tsv', rxns_default_file = 'github/Biochemistry/reactions.default.tsv', rxns_results_file = 'ModelSEED_rxns_master.py')

    # Template rxns
    read_template_rxns(rxns_tsv_file = 'github/Templates/GramPositive/Reactions.tsv', rxns_results_file = 'ModelSEED_rxns_GramPositive.py', results_var_name = 'rxns_GramPositive')

    read_template_rxns(rxns_tsv_file = 'github/Templates/GramNegative/Reactions.tsv', rxns_results_file = 'ModelSEED_rxns_GramNegative.py', results_var_name = 'rxns_GramNegative')

    read_template_rxns(rxns_tsv_file = 'github/Templates/Human/Reactions.tsv', rxns_results_file = 'ModelSEED_rxns_Human.py', results_var_name = 'rxns_Human')

    # Biomass
    #read_template_biomass(biomasses_tsv_file = 'github/Templates/GramPositive/Biomasses.tsv', biomass_cpds_tsv_file = 'github/Templates/GramPositive/BiomassCompounds.tsv',biomass_results_file = 'ModelSEED_rxns_GramPositive_Biomass.py', results_var_name = 'biomass_rxn_GramPositive')

    #read_template_biomass(biomasses_tsv_file = 'github/Templates/GramNegative/Biomasses.tsv', biomass_cpds_tsv_file = 'github/Templates/GramNegative/BiomassCompounds.tsv', biomass_results_file = 'ModelSEED_rxns_GramNegative_Biomass.py', results_var_name = 'biomass_rxn_GramNegative')

    #read_template_biomass(biomasses_tsv_file = 'github/Templates/Human/Biomasses.tsv', biomass_cpds_tsv_file = 'github/Templates/Human/BiomassCompounds.tsv', biomass_results_file = 'ModelSEED_rxns_Human_Biomass.py',results_var_name = 'biomass_rxn_Human')

    #write_cpds_formula_aliases(cpds_results_file = 'ModelSEED_cpds_formula_aliases.py')

