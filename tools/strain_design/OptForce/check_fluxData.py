from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../../../')
from copy import deepcopy
from tools.userError import userError
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from imp import load_source
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

def check_fluxData(flux_data_filename, organism_info = {'id':'', 'name':''}, model_info = {'id':'', 'sbml_filename':None, 'biomassrxn_id':None}, media_info = {'media_filename':None, 'other_flux_bounds': {}}, stdout_msgs = True, warnings = True):
    """
    Check the feasibility of an FBA problem for a given flux data  
    Creates the metabolic model from an SBML file

    INPUTS:
    ------
     flux_data_filenam: Name and path of file containing experimental flux data
         organism_info: A dictionary containing the information about he organism (id and anme)
            model_info: A dictionary containing the metabolic model information with the following keys:
                                    'id': model id
                         'sbml_filename': Name of and path of the SBML file containing the model
                         'biomassrxn_id': The if of biomass reaction
            media_info: Information about growth media
                         'media_filename': Name of and path of the file containing the flux bounds
                                           file for exchange reactions corresponding to compounds in 
                                           in the growth media
                        other_flux_bounds: A dictionary containing the flux bounds for other reactions (such as carbon
                                           source uptake, oxygen uptake, etc), that have to be set for reactions not 
                                           in flux_data_filenam or media_filename

    Ali R. Zomorrodi - Segre's Lab @ BU
    Last updated: 03-11-2016
    """
    # Import the flux bounds in flux_data_filename
    load_source('dataFile',flux_data_filename)
    import dataFile
    exp_flux_bounds = dataFile.exp_flux_bounds

    # Combine the flux bounds in other_flux_bounds and exp_flux_bounds
    flux_intersect = set(media_info['other_flux_bounds'].keys()).intersection(exp_flux_bounds.keys())
    if len(flux_intersect) > 0:
        print 'WARNING! The following fluxes are appear both in exp_flux_bounds and other_flux_bounds: {}. Their bounds in other_flux_bound will be overwritten by those in exp_flux_bounds'.format(flux_intersect)
    all_other_flux_bounds = media_info['other_flux_bounds']
    for rxn in exp_flux_bounds.keys():
        all_other_flux_bounds[rxn] = exp_flux_bounds[rxn]

    # Define the organism
    model_organism = organism(id = organism_info['id'], name = organism_info['name'])

    model = read_sbml_model(file_name = model_info['sbml_filename'], model_id = model_info['id'], model_organism = model_organism, model_type = 'metabolic',import_params = False)

    model.biomass_reaction = model.reactions_by_id[model_info['biomassrxn_id']]

    # Assign the objective function coefficients
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.biomass_reaction.objective_coefficient = 1

    # Growth medium
    set_specific_bounds(model = model, file_name = media_info['media_filename'], flux_bounds = all_other_flux_bounds)

    # Perform FBA for the wild-type
    model.fba(assign_wildType_max_biomass = True, stdout_msgs = stdout_msgs)

def Ecoli_minimalAerobic(): 
    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/'

    check_fluxData(flux_data_filename = '/usr2/postdoc/alizom/work/tools/strain_design/flux_data/Escherichia_coli/iJO1366_minimal_aerobic_glc_PMID_23036703.py', organism_info = {'id':'Ecoli', 'name':'Escherichia_coli'}, model_info = {'id':'iJO1366_minimal_aerobic_glc_PMID_23036703', 'sbml_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, media_info = {'media_filename':model_path + 'iJO1366_minimal_aerobic_glucose.py', 'other_flux_bounds': {'EX_glc(e)':[-10,1000], 'EX_o2(e)':[-20,1000]}}, stdout_msgs = True, warnings = True)

    #model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iAF1260/'

    #check_fluxData(flux_data_filename = '/usr2/postdoc/alizom/work/tools/strain_design/flux_data/Escherichia_coli/iJO1366_minimal_aerobic_glc_PMID_23036703.py', organism_info = {'id':'Ecoli', 'name':'Escherichia_coli'}, model_info = {'id':'iAF1260_minimal_aerobic_glc_PMID_23036703', 'sbml_filename':model_path + 'iAF1260_updated.xml', 'biomassrxn_id':'Ec_biomass_iAF1260_core_59p81M'}, media_info = {'media_filename':model_path + 'iAF1260_minimal_glucose_aerobic.py', 'other_flux_bounds': {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}}, stdout_msgs = True, warnings = True)

#-------------------------
if __name__ == '__main__':
    Ecoli_minimalAerobic()
 
