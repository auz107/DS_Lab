from __future__ import division
import sys
sys.path.append('../')
from tools.ancillary.importData import importData
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
from models.model_seed_database.ModelSeed_compounds import ModelSeed_compounds, ModelSeed_compounds_by_clean_name
import re

if __name__ == '__main__':
    """
    Reads the experimental data and returns the results as python dicitonaries
    """
    # Load the experimental growth data on different carbon sources. This will be a dictionary with 
    # keys tuples of (species_name,carbon_soruce_name) and values being their correspondig ODs 
    growth_data = importData(inputFile = 'expData/growth_single.txt',delType = 'tab',dataType = 'float')

    # List of all organisms 
    organism_list = sorted(list(set([k[0] for k in growth_data.keys()])))

    # List of all carbon sources
    carbonSrc_list = sorted(list(set([k[1] for k in growth_data.keys()])))

    # Remove water from the list of carbon sources
    del carbonSrc_list[carbonSrc_list.index('water')]

    #--- Get ModelSeed ids for carbon sources ----
    # Create a dictionary with keys being the names of carbon sources and values their corresponding ModelSeed ids
    carbonSrc_dict = dict([(k,None) for k in carbonSrc_list])
    ModelSeed_cleanNames = ModelSeed_compounds_by_clean_name.keys() 
    for carbonSrc in carbonSrc_list: 
         clean_name = remove_non_alphanumeric(input_string = carbonSrc).lower()
         if clean_name in ModelSeed_cleanNames:
             ModelSeed_id = ModelSeed_compounds_by_clean_name[clean_name]
             if len(ModelSeed_id) > 1:
                 print '**More than one ModelSeed id for ',carbonSrc,' (',clean_name,')'
             ModelSeed_id = ModelSeed_id[0]
             carbonSrc_dict[carbonSrc] = ModelSeed_id
             print carbonSrc,'\t',carbonSrc_dict[carbonSrc],'\t',ModelSeed_compounds[ModelSeed_id]['name']

    carbonSrc_dict['4-hydroxy benzoic acid'] = 'cpd00136'
    carbonSrc_dict['alpha-D-lactose'] = 'cpd00208'
    carbonSrc_dict['D-cellobiose'] = 'cpd00158'
    carbonSrc_dict['2-hydroxy benzoic acid'] = 'cpd00599' 
    carbonSrc_dict['hydroxybutyric acid'] = 'cpd00728'    
    carbonSrc_dict['alpha-D-lactose'] = 'cpd00208'  
    carbonSrc_dict['D,L-alpha-glycerol phosphate'] = 'cpd00080'       
    carbonSrc_dict['i-erythritol'] = 'cpd00392'       
    carbonSrc_dict['galacturonic acid'] = ['cpd02143','cpd00280']  
    carbonSrc_dict['phenylalanine'] = ['cpd01400','cpd00066']  
    carbonSrc_dict['D-galactonic acid gamma lactone'] = 'cpd02143'  

    # Carbon sources with no ModelSeed id    
    noModelSeed_id = [c for c in carbonSrc_dict.keys() if carbonSrc_dict[c] == None]

    # Organisms with annotated genomes
    annotated_orgs = ['Enterobacter aerogenes', 'Pseudomonas chlororaphis', 'Pseudomonas fluorescens', 'Pseudomonas putida', 'Pseudomonas veronii', 'Serratia marcescens']

    #------ Identifying the list of carbon sources each organism can grow on --------
    # Any OD of greater than or equal to 0.025 is considered as grwoth
    growth_OD_thr = 0.025
 
    # Consider only organisms with annotated genomes that grow on a given carbon source with a ModelSeed id
    grow_annotated_ModelSeedId = [k for k in growth_data.keys() if growth_data[k] >= growth_OD_thr and carbonSrc_dict[k[1]] != None and k[0] in annotated_orgs]

    # List of all carbon sources in grow_annotated_ModelSeedId
    print '\nList of all carbon sources in the experiimental data for mutants that can grow:\n'
    carbonSrc_in_grow_annotated_ModelSeedId = sorted(set([k[1] for k in grow_annotated_ModelSeedId]))
    for carbonSrc in carbonSrc_in_grow_annotated_ModelSeedId:
        print carbonSrc_in_grow_annotated_ModelSeedId.index(carbonSrc) + 1,'\t',carbonSrc,'\t',carbonSrc_dict[carbonSrc]

    print '\n\n'
    for org in list(set([k[0] for k in grow_annotated_ModelSeedId])): 
        print '\n---- ',org,' ----'
        print ['minimal_aerobic_' + remove_non_alphanumeric(k[1]) for k in grow_annotated_ModelSeedId if k[0] == org]
        carbonSrc_counter = 1
        for carbonSrc in [k[1] for k in grow_annotated_ModelSeedId if k[0] == org]:
            print '\t%i\t%s\t%s\t\t%.4f' % (carbonSrc_counter,carbonSrc,carbonSrc_dict[carbonSrc],growth_data[k])
            carbonSrc_counter += 1
 
    print '\nThe total number of case = %i' % len(growth_data)
    print 'The total number of species that could grow = %i' % len(grow_annotated_ModelSeedId)
    print 'No ModelSeed id was found for the following ',len([c for c in carbonSrc_dict.keys() if carbonSrc_dict[c] == None]),' compounds'
    for c in noModelSeed_id:
        print '\t',c

    #---- Check which models can grow on which carbon sources after kbase gap filling ----
    inSilico_grwoth_data = {'Pseudomonas putida':['L-arginine','L-serine','L-asparagine','galacturonic acid','putrescine','phenylalanine','4-hydroxy benzoic acid']}




