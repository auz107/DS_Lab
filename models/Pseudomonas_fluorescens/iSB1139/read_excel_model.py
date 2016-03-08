import sys, re
sys.path.append('/data/alizom')
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.compartment import compartment
from tools.core.model import model
from imp import load_source

def create_python_model(metabs_raw_data_file_name,rxns_raw_data_file_name,model_file_name_base):
    """
    Creates the python file containing the model in the Excel file
  
    INPUTS:
    ------
         metabs_raw_data: A string containing the path and file name containing the 
                          raw data from compounds extracted from the Excel file
           rxns_raw_data: A string containing the path and file name containing the 
                          raw data from reactions extracted from the Excel file
    model_file_name_base: Name of the output python file containing the file    
    """

    
    # The following writes the compound and reaction information into a dicitonary and stores
    # the results into a file
    
    #--- compounds ----
    with open(model_file_name_base + '_tmp.py','w') as f:
        f.write('metabs_tmp = {}\n')
    
    metabs_tmp = {}
    
    with open(metabs_raw_data_file_name,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read hears
                headers = line
            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of raw_compounds.txt. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                # id
                if 'replacement_id' not in headers:
                    raise userError('**ERROR! replacement_id is not in the list of headers')
                else:
                    rep_id = headers.index('replacement_id')
                    metabs_tmp[line[rep_id]] = {}
                    for header in headers:
                        # index of header in headers
                        header_ind = headers.index(header)
                        if line[header_ind] == 'NULL':
                            metabs_tmp[line[rep_id]][header] = None
                        else:
                            metabs_tmp[line[rep_id]][header] = line[header_ind] 
                    #with open(model_file_name_base + '_tmp.py','a') as f:
                    #    f.write("metabs_tmp['" + line[rep_id] + "'] = " + repr(metabs_tmp[line[rep_id]]) + '\n')
    
    #--- reactions ---
    with open(model_file_name_base + '_tmp.py','a') as f:
        f.write('\n\nrxns_tmp = {}\n')
    
    rxns_tmp = {}
    
    with open(rxns_raw_data_file_name,'r') as f:
        counter = 0
        for line in f:
            line = line.strip().split('\t')
            counter += 1
            if counter == 1:
                # Read hears
                headers = line
            elif len(line) != len(headers):
                print '**ERROR! in line ', counter, ' of raw_reactions.txt. The number of items (', len(line), ') is not equal to the number of headers (', len(headers),')\n'
                print 'line = ',line,'\n'
                raise userError()
            else:
                # id
                if 'id' not in headers:
                    raise userError('**ERROR! id is not in the list of headers')
                else:
                    rep_id = headers.index('id')
                    rxns_tmp[line[rep_id]] = {}
                    for header in headers:
                        # index of header in headers
                        header_ind = headers.index(header)
                        if line[header_ind] == 'NULL':
                            rxns_tmp[line[rep_id]][header] = None
                        else:
                            rxns_tmp[line[rep_id]][header] = line[header_ind] 
                    #with open(model_file_name_base + '_tmp.py','a') as f:
                    #    f.write("rxns_tmp['" + line[rep_id] + "'] = " + repr(rxns_tmp[line[rep_id]]) + '\n')
    
    #------------------------------------------------------------------
    #from iSB1139_tmp import metabs_tmp
    #from iSB1139_tmp import rxns_tmp
    
    # Replace => and <=> in reaction equations with --> and  <==>, respectively
    for rxn_id in rxns_tmp.keys():
        if ' <=> ' in rxns_tmp[rxn_id]['equation']:
            rxns_tmp[rxn_id]['equation'] = re.sub(' <=> ',' <==> ',rxns_tmp[rxn_id]['equation'])
        elif ' => ' in rxns_tmp[rxn_id]['equation']:
            rxns_tmp[rxn_id]['equation'] = re.sub(' => ',' --> ',rxns_tmp[rxn_id]['equation'])
        else:
            raise userError('Neither => nor <=> found in the equation of reaction ' + rxn_id)
    
    # Remove " from compound ids and reaction names 
    for metab_rid in [mid for mid in metabs_tmp.keys() if '"' in metabs_tmp[mid]['id']]:
        metabs_tmp[metab_rid]['id'] = re.sub('"','',metabs_tmp[metab_rid]['id'])
    for rxn_rid in [rid for rid in rxns_tmp.keys() if '"' in rxns_tmp[rid]['equation']]:
        rxns_tmp[rxn_rid]['equation'] = re.sub('"','',rxns_tmp[rxn_rid]['equation'])
    
    # Replace space in any compound ids with an underline both in the compound id and 
    # in the reaction equation 
    for metab_rid in [mid for mid in metabs_tmp.keys() if ' ' in metabs_tmp[mid]['id']]:
        # Metabolite id with no space
        no_space_id = re.sub(' ','_',metabs_tmp[metab_rid]['id'])
        if len([rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            raise userError('**ERROR! compound ' + metabs_tmp[metab_rid]['id'] + ' was not found in any reactions\n')
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace the old id in all reaction equations with new compound id 
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',no_space_id,rxns_tmp[rxn_id]['equation'])
        metabs_tmp[metab_rid]['id'] = no_space_id 

    # Replace [extracellular] in compound ids and reaction equations with [e]
    # Spaces should be removed before other corrections
    for metab_rid in [mid for mid in metabs_tmp.keys() if '[extracellular]' in metabs_tmp[mid]['id']]:
        # Metabolite id with no [extracellular]
        no_ex_id = re.sub('\[extracellular\]','[e]',metabs_tmp[metab_rid]['id'])
        if len([rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            raise userError('**ERROR! compound ' + metabs_tmp[metab_rid]['id'] + ' was not found in any reactions\n')

        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace the old id in all reaction equations with new compound id 
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',no_ex_id,rxns_tmp[rxn_id]['equation'])
        metabs_tmp[metab_rid]['id'] = no_ex_id 
    
    # Replace comma (,) in any compound ids with an underline both in the compound id and 
    # in the reaction equation 
    for metab_rid in [mid for mid in metabs_tmp.keys() if ',' in metabs_tmp[mid]['id']]:
        # Metabolite id with no comma
        no_comma_id = re.sub(',','_',metabs_tmp[metab_rid]['id'])
        if len([rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            raise userError('**ERROR! compound ' + metabs_tmp[metab_rid]['id'] + ' was not found in any reactions\n')
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace the old id in all reaction equations with new compound id 
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',no_comma_id,rxns_tmp[rxn_id]['equation'])
        metabs_tmp[metab_rid]['id'] = no_comma_id 
    
    # Replace __ in any compound ids with an underline both in the compound id and 
    # in the reaction equation 
    for metab_rid in [mid for mid in metabs_tmp.keys() if '__' in metabs_tmp[mid]['id']]:
        # Metabolite id with no space
        one_underline_id = re.sub('__','_',metabs_tmp[metab_rid]['id'])
        if len([rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            raise userError('**ERROR! compound ' + metabs_tmp[metab_rid]['id'] + ' was not found in any reactions\n')
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace the old id in all reaction equations with new compound id 
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',one_underline_id,rxns_tmp[rxn_id]['equation'])
        metabs_tmp[metab_rid]['id'] = one_underline_id 
    
    # Replace _ if it appears at the end of a compound id
    p = re.compile('_$')
    for metab_rid in [mid for mid in metabs_tmp.keys() if p.search(metabs_tmp[mid]['id']) != None]:
        # Metabolite id with no space
        no_underline_id = re.sub('_$','',metabs_tmp[metab_rid]['id'])
        if len([rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            raise userError('**ERROR! compound ' + metabs_tmp[metab_rid]['id'] + ' was not found in any reactions\n')
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace the old id in all reaction equations with new compound id 
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',no_underline_id,rxns_tmp[rxn_id]['equation'])
        metabs_tmp[metab_rid]['id'] = no_underline_id 
    
    
    # Replace '+' in any compound ids with 'plus' both in the compound id and 
    # in the reaction equation 
    for metab_rid in [mid for mid in metabs_tmp.keys() if '+' in metabs_tmp[mid]['id']]:
        # Metabolite id with no space
        no_plus_id = re.sub('\+','plus',metabs_tmp[metab_rid]['id'])
        if len([rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            raise userError('**ERROR! compound ' + metabs_tmp[metab_rid]['id'] + ' was not found in any reactions\n')
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace the old id in all reaction equations with new compound id 
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',no_plus_id,rxns_tmp[rxn_id]['equation'])
        metabs_tmp[metab_rid]['id'] = no_plus_id 
    
    # -- Create exchange and transport reactions ---
    # (1) Remove all compound containing _drain in their name from the list of compounds and
    # replace all reactions like A ==> A_draint to exchange rxns of the form A[e] <==> and add  
    # EX_ to the begining of their id. Also, change their compartment to e 
    for metab_rid in [mid for mid in metabs_tmp.keys() if '_drain' in metabs_tmp[mid]['id']]:
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + metabs_tmp[metab_rid]['id'] + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Remove metabName_drain from the reaction equation
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))','',rxns_tmp[rxn_id]['equation'])
                 
            # Replace A in A --> with A[e]
            no_drain_name = re.sub('_drain','',metabs_tmp[metab_rid]['id'])
            if no_drain_name + '[e]' in [metabs_tmp[mid]['id'] for mid in metabs_tmp.keys()]:
                rxns_tmp[rxn_id]['equation'] = re.sub(no_drain_name,no_drain_name + '[e]',rxns_tmp[rxn_id]['equation'])
            else:
                print metabs_tmp[metab_rid]['id'],'     ',no_drain_name + ' has no [e] version in the model'
    
            # Replace --> with <==>
            rxns_tmp[rxn_id]['equation'] = re.sub(' --> ',' <==> ',rxns_tmp[rxn_id]['equation'])
    
            # Add 'EX_' to the begining of the reaction id
            new_id = 'EX_' + rxns_tmp[rxn_id]['id']
            rxns_tmp[rxn_id]['id'] = new_id 
            rxns_tmp[new_id] = rxns_tmp.pop(rxn_id)
            rxns_tmp[new_id]['compartment'] = 'e' 
            rxns_tmp[new_id]['subsystem'] = 'exchange' 
    
        # Delete compound from the list of compounds
        del metabs_tmp[metab_rid]
    
    # Remove all compounds containing _input in their name from the list of compounds and
    # replace all reactions like A_input <==> A[e] to a transport reaction with the following 
    # equation: A <==> A[e]
    if 'Co2plus_input' not in [metabs_tmp[mid]['id'] for mid in metabs_tmp.keys() if '_input' in metabs_tmp[mid]['id']]:
        raise userError('**Ohhhhh')
    
    for metab_rid in [mid for mid in metabs_tmp.keys() if '_input' in metabs_tmp[mid]['id']]:
        for rxn_id in [rid for rid in rxns_tmp.keys() if re.compile('((?<=\s)|(?<=^))' + metabs_tmp[metab_rid]['id'] + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]:
            # Replace A_input in reaction equaiton with A
            no_input_name = re.sub('_input','',metabs_tmp[metab_rid]['id'])
    
            # Remove _input from the reaction equation
            rxns_tmp[rxn_id]['equation'] = re.sub('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))',no_input_name,rxns_tmp[rxn_id]['equation'])
    
            # Replace --> with <==>
            #rxns_tmp[rxn_id]['equation'] = re.sub(' --> ',' <==> ',rxns_tmp[rxn_id]['equation'])
                   
            # Add 'Transport_' to the begining of the reaction id
            new_id = 'Transport_' + rxns_tmp[rxn_id]['id']
            rxns_tmp[rxn_id]['id'] = new_id 
            rxns_tmp[new_id] = rxns_tmp.pop(rxn_id)
    
        # Delete compound from the list of compounds
        del metabs_tmp[metab_rid]
    
    #--- Extract the Kegg_ids of compounds from DB LINK, if they exist ---
    for metab_rid in metabs_tmp.keys():
        if metabs_tmp[metab_rid]['DB LINK'] != None:
            pattern = re.compile('""(C[0-9][0-9][0-9][0-9][0-9])""')
            s = pattern.search(metabs_tmp[metab_rid]['DB LINK'])
            if s != None:
                metabs_tmp[metab_rid]['Kegg_id'] = s.group(1) 
            else:
                metabs_tmp[metab_rid]['Kegg_id'] = None
        else:
            metabs_tmp[metab_rid]['Kegg_id'] = None

    #--- Add Kegg_id for compounds manually
    metabs_tmp['r809']['Kegg_id'] = 'C00076'  # Ca2plus 
    metabs_tmp['r810']['Kegg_id'] = 'C00076'  # Ca2plus[e]
    metabs_tmp['r851']['Kegg_id'] = 'C05776'  # coenzyme_B12 
    metabs_tmp['r280']['Kegg_id'] = 'C05776'  # coenzyme_B12[e]
    metabs_tmp['r821']['Kegg_id'] = 'C00115'  # chloride[e]. For somoe reaons kegg has a different id (C00698)
    metabs_tmp['r281']['Kegg_id'] = 'C00011'  # CO2[e].
    metabs_tmp['r282']['Kegg_id'] = 'C00175'  # Co2plus[e]
    metabs_tmp['r293']['Kegg_id'] = 'C00023'  # Fe2plus[e]
    metabs_tmp['r294']['Kegg_id'] = 'C14819'  # Fe3plus
    metabs_tmp['r295']['Kegg_id'] = 'C14819'  # Fe3plus[e]
    metabs_tmp['r298']['Kegg_id'] = 'C00001'  # H2O[e]
    metabs_tmp['r306']['Kegg_id'] = 'C00034'  # Mn2plus 
    metabs_tmp['r307']['Kegg_id'] = 'C00034'  # Mn2plus[e]
    metabs_tmp['r771']['Kegg_id'] = 'C00014'  # NH3 
    metabs_tmp['r772']['Kegg_id'] = 'C00014'  # NH3[e]
    metabs_tmp['r581']['Kegg_id'] = 'C00031'  # beta-D-glucose 
    metabs_tmp['r284']['Kegg_id'] = 'C00031'  # beta-D-glucse[e]
    metabs_tmp['r309']['Kegg_id'] = 'C00007'  # O2 

    #--- Find all [e] compounds with no exchange reaction in the model ----
    # and add a corresponding exchange reaction to the model 
    for metab_rid in [mid for mid in metabs_tmp.keys() if '[e]' in metabs_tmp[mid]['id']]:
        if len([rid for rid in rxns_tmp.keys() if 'EX_' in rxns_tmp[rid]['id'] and re.compile('((?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '((?=\s)|(?=$))').search(rxns_tmp[rid]['equation']) != None]) == 0:
            rxns_tmp['EX_' + metabs_tmp[metab_rid]['replacement_id']] = {'id':'EX_' + metabs_tmp[metab_rid]['replacement_id'],'name':metabs_tmp[metab_rid]['name'] + ' exchange','equation': metabs_tmp[metab_rid]['id'] + ' <==> ', 'subsystem':'exchange','genes': None,'EC. Number': None,'compartment':'e'}

    #---- replace compound names in reaction equations with their replacement ids -----
    # re.sub(r'(?:(?<=\s)|(?<=^))H\+(?=\s|$)','C0001','X_H+ + H+ --> H+[c]') 
    for metab_rid in metabs_tmp.keys():
        for rxn_id in rxns_tmp.keys():
            rxns_tmp[rxn_id]['equation'] = re.sub('(?:(?<=\s)|(?<=^))' + re.escape(metabs_tmp[metab_rid]['id']) + '(?=\s|$)',metabs_tmp[metab_rid]['replacement_id'],rxns_tmp[rxn_id]['equation'])
    
    #--- Delete Transport_MIRXN_88 ---
    # This reaciotn needs to be removed because the only tranposrt reaction for 
    # rehalose in the model is ENZRXN3EM_3828 the model does not contain trehalose
    # as a compound
    del rxns_tmp['Transport_MIRXN_88']
    
    #------------ Add stoichiometric coefficients equal to one -----
    for rxn_id in rxns_tmp.keys():
            rxns_tmp[rxn_id]['equation'] = re.sub('^r','1 r',rxns_tmp[rxn_id]['equation'])
            rxns_tmp[rxn_id]['equation'] = re.sub('\+ r','+ 1 r',rxns_tmp[rxn_id]['equation'])
            rxns_tmp[rxn_id]['equation'] = re.sub('--> r','--> 1 r',rxns_tmp[rxn_id]['equation'])
            rxns_tmp[rxn_id]['equation'] = re.sub('<==> r','<==> 1 r',rxns_tmp[rxn_id]['equation'])
    
    #------- Create the stoichiometric matrix --------------
    # First extract the reactnats and products
    for rxn_id in rxns_tmp.keys():
        rxns_tmp[rxn_id]['stoic'] = {}
        p = re.compile('(^.*(?= --> ))|(^.*(?= <==> ))')
        s = p.search(rxns_tmp[rxn_id]['equation'])
        if s != None:
            reactants = s.group() 
        else:
            raise userError('**ERROR! No reactants found')
    
        p = re.compile('((?<= --> ).*$)|((?<= <==> ).*$)')
        s = p.search(rxns_tmp[rxn_id]['equation'])
        if s != None:
            products = s.group() 
        else:
            products = None
    
        # Add '+ ' to the begining and ' +' to the end of the reactants and products
        reactants = '+ ' + reactants + ' +'
        products = '+ ' + products + ' +'
    
        #-- reactans ---
        p = re.compile('(?<=\+ ).*?(?= \+)')
        # Pairs of stoichiometric 'coefficients compound'
        for s in p.finditer(reactants):
            if s != None:
                sm = s.group()
    
                # Stoichioimetric coefficient
                ps = re.compile('^.*(?= )')
                s = ps.search(sm)
                if s != None:
                    scoeff = s.group()
                else:
                    raise userError('**ERROR! No stoichiometric coefficient for ' + sm + 'in reaction ' + rxn_id + ': ' + rxns_tmp[rxn_id]['equation'])
    
                # Metabolite name 
                pm = re.compile('(?<= ).*$')
                s = pm.search(sm)
                if s != None:
                    mname = s.group()
                else:
                    raise userError('**ERROR! No compound name for (' + sm + ') in reaction ' + rxn_id)
    
                try:
                    rxns_tmp[rxn_id]['stoic'][mname] = float('-' + scoeff)
                except:
                    raise userError('**ERROR! Invlid extracted stoichiometric coefficient for compound ' + mname + ' in reaction ' + + rxn_id + ': ' + rxns_tmp[rxn_id]['equation'] + '   s = ' + scoeff + '\n')
            else:
                raise userError('**ERROR! Not stoichiometric-compound pari for reaction ' + rxn_id)
    
        #-- products ---
        # Do not consider exchange reactions
        if 'EX_' not in rxns_tmp[rxn_id]['id']:
            p = re.compile('(?<=\+ ).*?(?= \+)')
            # Pairs of stoichiometric 'coefficients compound'
            for s in p.finditer(products):
                if s != None:
                    sm = s.group()
        
                    # Stoichioimetric coefficient
                    ps = re.compile('^.*(?= )')
                    s = ps.search(sm)
                    if s != None:
                        scoeff = s.group()
                    else:
                        raise userError('**ERROR! No stoichiometric coefficient for ' + sm + 'in reaction ' + rxn_id + ': ' + rxns_tmp[rxn_id]['equation'])
        
                    # Metabolite name 
                    ps = re.compile('(?<= ).*$')
                    s = pm.search(sm)
                    if s != None:
                        mname = s.group()
                    else:
                        raise userError('**ERROR! No compound name for (' + sm + ') in reaction ' + rxn_id)
        
                    try:
                        rxns_tmp[rxn_id]['stoic'][mname] = float(scoeff)
                    except:
                        raise userError('**ERROR! Invlid extracted stoichiometric coefficient for compound ' + mname + ' in reaction ' + rxn_id + ': ' + rxns_tmp[rxn_id]['equation'] + '   s = ' + scoeff + '\n')
                else:
                    raise userError('**ERROR! Not stoichiometric-compound pari for reaction ' + rxn_id)
    

    #----- Writing ghe fianl results into a file ---------------------
    model_file_name = 'iSB1139.py'
    
    compartments = ['c','e']
    metabs = {}
    rxns = {}
    
    with open(model_file_name,'w') as f:
        f.write('metabs = {}\n')
    
    with open(model_file_name,'a') as f:
        for metab_rid in metabs_tmp.keys():
            id = metabs_tmp[metab_rid]['replacement_id']     
            metabs[id] = {}
            metabs[id]['id'] = id 
            metabs[id]['name'] = metabs_tmp[metab_rid]['id']     
            metabs[id]['Kegg_id'] = metabs_tmp[metab_rid]['Kegg_id']     
            metabs[id]['compartment'] = metabs_tmp[metab_rid]['compartment']     
                
            f.write("metabs['" + id + "'] = " + repr(metabs[id]) + '\n')
    
        # Reactions (store only the fields that are needed)
        f.write('\n\nrxns = {}\n')
    
        for rxn_id in rxns_tmp.keys():
            id = rxns_tmp[rxn_id]['id']
            rxns[id] = {}
            rxns[id]['id'] = rxns_tmp[rxn_id]['id']
            rxns[id]['name'] = rxns_tmp[rxn_id]['name']
            rxns[id]['equation'] = rxns_tmp[rxn_id]['equation']
            if ' --> ' in rxns_tmp[rxn_id]['equation']:
                rxns[id]['type'] = 'irreversible'
            elif ' <==> ' in rxns_tmp[rxn_id]['equation'] and 'EX_' in rxns_tmp[rxn_id]['id']:
                rxns[id]['type'] = 'exchange'
            elif ' <==> ' in rxns_tmp[rxn_id]['equation'] and 'EX_' not in rxns_tmp[rxn_id]['id']:
                rxns[id]['type'] = 'reversible'
            else:
                raise userError('**ERROR! Unknown reaction type for reaction ' + id + '\n')
            rxns[id]['subsystem'] = rxns_tmp[rxn_id]['subsystem']
            rxns[id]['genes'] = rxns_tmp[rxn_id]['genes']
            rxns[id]['EC_number'] = rxns_tmp[rxn_id]['EC. Number']
            rxns[id]['compartment'] = rxns_tmp[rxn_id]['compartment']
            rxns[id]['stoic'] = rxns_tmp[rxn_id]['stoic']
    
            f.write("rxns['" + id + "'] = " + repr(rxns[id]) + '\n')
            #print rxns[rxn_id]['id'],':   ',rxns[rxn_id]['equation'],'\n'

    
def create_model(python_model_file_name,organism_name, model_type, model_name):
    """
    Creates the model file

    INPUTS:
    ------
    python_model_file_name: The python file name containing the model and 
                            created by create_python_model
             organism_name: Name of the organism
                model_type: Type of the model (string, e.g., 'metabolic')
                model_name: Name of the model
    """
    # Import the data in module stored in gams_model_file
    if type(python_model_file_name) == str:
        load_source('python_model',python_model_file_name)
    import python_model

    # organism
    org = organism(id = organism_name)

    # Compartments
    compartments = []
    get_compt_by_id = {} 
    compt_c = compartment(id = 'c', name = 'cytosol')
    compartments.append(compt_c)
    get_compt_by_id['c'] = compt_c
    compt_e = compartment(id = 'e', name = 'extracellular')
    compartments.append(compt_e)
    get_compt_by_id['e'] = compt_e

    # Metabolites
    get_metab_by_id = {}
    compounds = []
    for mid in sorted(python_model.metabs.keys()):
        m = compound(id = python_model.metabs[mid]['id'], name = python_model.metabs[mid]['name'], Kegg_id = python_model.metabs[mid]['Kegg_id'], compartment = get_compt_by_id[python_model.metabs[mid]['compartment']])
        compounds.append(m)
        get_metab_by_id[mid] = m

    # reactions
    get_rxn_by_id = {}
    reactions = []
    for rid in sorted(python_model.rxns.keys()):
        rxn_stoic = dict([(get_metab_by_id[ms[0]],ms[1]) for ms in  python_model.rxns[rid]['stoic'].items()])
        if len(rxn_stoic) == 0:
            raise userError('\n**ERROR! Empty reaction stoichiometry for reaction ' + rid,'\n')
        r = reaction(id = python_model.rxns[rid]['id'], stoichiometry = rxn_stoic, type = python_model.rxns[rid]['type'], name = python_model.rxns[rid]['name'], EC_number = python_model.rxns[rid]['EC_number'], subsystem = python_model.rxns[rid]['subsystem'], compartments = get_compt_by_id[python_model.rxns[rid]['compartment']])
        reactions.append(r)
        get_rxn_by_id[rid] = r

    # Create the model
    return model(id = model_name, type = model_type, organism = org, reactions = reactions, compounds = compounds, compartments = compartments, biomass_reaction = get_rxn_by_id['EX_GROWTH'])

#-------------------------------------------------------------------------
if __name__ == '__main__':
    create_python_model(metabs_raw_data_file_name = 'raw_compounds.txt',rxns_raw_data_file_name = 'raw_reactions.txt',model_file_name_base = 'iSB1139')
    #model = create_model(python_model_file_name = 'iSB1139.py',organism_name = 'Pseudomonas fluorescens', model_type = 'metabolic', model_name = 'iSB1139')
