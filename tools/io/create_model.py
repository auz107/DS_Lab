from __future__ import division
import sys, time
sys.path.append('../../')
from tools.fba.fbaTools import fbaTools
from tools.fba.fba import fba
from tools.userError import userError
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.io.read_sbml_model import read_sbml_model
from tools.io.read_pydict_model import read_pydict_model
from tools.fba.set_specific_bounds import set_specific_bounds
from coopr import pyomo

def create_model(model_organism = {}, model_info = {'id':'', 'file_format':'sbml', 'model_filename':None, 'biomassrxn_id':''}, growthMedium_flux_bounds = {'flux_bounds_filename':None, 'flux_bounds_dict': {}}, import_params = False, validate = True, perform_fba = True, stdout_msgs = True, warnings = True):
    """
    Creates a metabolic model and assigns the flux bounds for a desired growth media and performs FBA.

    INPUTS:
    ------
              organism_info: A dictionary containing the information about he organism (id and name). 'organism_object' can be an instance
                             of class organism. If 'organism_object' is provided, there is no need to provide 'anem' and 'id'
                 model_info: A dictionary containing the metabolic model information with the following keys:
                                          'id': model id
                                 'file_format': Type of the model file (allowed choices are sbml, pydict, pickle, json)
                              'model_filename': Name of and path of the SBML file containing the model
                               'biomassrxn_id': The if of biomass reaction
    growthMedium_flux_bounds: Information about growth media
                             flux_bounds_filename: Name of and path of the file containing the flux bounds
                                                  file for exchange reactions corresponding to compounds in 
                                                  in the growth media
                                 flux_bounds_dict: A dictionary containing the flux bounds for other reactions (such as carbon
                                                  source uptake, oxygen uptake, etc), that have to be set for reactions not 
                                                  in flux_data_filenam or media_filename
                perform_fba: If True, fba is performed using the specified flux bounds
              import_params: If Ture, model parameters (flux bounds, etc) are imported  

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 04-18-2016
    """
    # Perform checks on the input arguments
    if not isinstance(model_organism,organism) and not isinstance(model_organism,dict):
        raise TypeError('model_organism must be an instance of class organism or a dictionary')

    if not isinstance(model_info,dict):
        raise TypeError('model_info must be a dictionary')

    if not isinstance(growthMedium_flux_bounds,dict):
        raise TypeError('growthMedium_flux_bounds must be a dictionary')

    if not isinstance(warnings,bool):
        raise TypeError('warnings must be either True or False')

    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be either True or False')

    if model_info['file_format'].lower() == 'sbml':
        model = read_sbml_model(file_name = model_info['model_filename'], model_id = model_info['id'], model_organism = model_organism, model_type = 'metabolic',import_params = False, validate = validate)
    elif model_info['file_format'].lower() == 'pydict':
        model = read_pydict_model(file_name = model_info['model_filename'], import_params = import_params, validate = validate, warnings = warnings, stdout_msgs = stdout_msgs)
    elif model_info['file_format'].lower() == 'pickle':
        pass
    elif model_info['file_format'].lower() == 'json':
        pass
    else:
        raise userError("Invalid model_info['file_format']: {}. Allowed options are sbml, pydict, pickle and json".format(model_info['file_format'])) 

    if model_info['biomassrxn_id'] != '':
        model.biomass_reaction = model.reactions_by_id[model_info['biomassrxn_id']]

    # Perform FBA for the wild-type
    if perform_fba:
        # Assign the objective function coefficients
        for rxn in model.reactions:
            rxn.objective_coefficient = 0
        model.biomass_reaction.objective_coefficient = 1

        # Growth medium
        set_specific_bounds(model = model, file_name = growthMedium_flux_bounds['flux_bounds_filename'], flux_bounds = growthMedium_flux_bounds['flux_bounds_dict'])

        model.fba(assign_wildType_max_biomass = True, stdout_msgs = stdout_msgs)

    return model

