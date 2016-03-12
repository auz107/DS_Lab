from __future__ import division
import sys, os
sys.path.append('../../')
from tools.userError import userError
from tools.core.compound import compound 
from tools.core.reaction import reaction
from tools.core.model import model
from imp import load_source


def set_specific_bounds(model, flux_bounds = {}, file_name = None, simulation_condition = None, reset_flux_bounds = True, **additional_args): 
    """
    A function to read a gams growth medium information stored 
    in python-formatted and then to convert it into the format 
    that we need

    INPUTS:
    ------
                   model: A model object
             flux_bounds: The rest of the inputs should be in the form of a dicntionary 
                          with reaction ids as keys and a list {LB,UB} containing the lower 
                          and upper bounds on reactions
    simulation_condition: Name of the simulation condition (string)
       reset_flux_bounds: Resets the flux bonds to default if True
         additional_args: The user can use desirable keyword arguments to enter flux bounds, e.g.,
                          limiting_nutrients, excess_nutrients, regulation, etc. The should all 
                          be in the form of a dictionary similar to flux_bounds though

    NOTE: If a file name provided but some reaction fluxes are specified through the input 
          arguments they will overwrite those given in the file

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 03-11-2016
    """
    # Reset all reaction bounds
    if reset_flux_bounds:
        model.reset_flux_bounds()

    model.simulation_condition = simulation_condition

    #--- First set hte flux bounds from the file ---
    # A list containing all flux bounds
    flux_bounds_fromFile = [] 

    # Import the data in the module stored in file_name
    if type(file_name) == str:
        if not os.path.isfile(file_name):
            raise IOError("No such file was found :'" + file_name + "'")
        else:
            # First delete the model dataFile if it already exists. If it is not deleted
            # the new module data is merged with the previous ones
            try:
                del sys.modules['dataFile']
            except:
                pass
            load_source('dataFile',file_name)
            import dataFile 

            for s in [s for s in dir(dataFile) if '__' not in s and "'_" not in s]:
                exec 'flux_bounds_fromFile += dataFile.' + s + '.items()'

    flux_bounds_fromFile = dict(flux_bounds_fromFile)
    for rxn_id in flux_bounds_fromFile.keys(): 
        rxn = model.reactions_by_id[rxn_id]
        if rxn == None:
            raise ValueError('Reaction ' + rxn_id + ' was not found in the model.')
        else:
            if rxn.flux_bounds == []:
                rxn.assign_flux_bounds()
            if flux_bounds_fromFile[rxn_id][0] != None:
                rxn.flux_bounds[0] = flux_bounds_fromFile[rxn_id][0]
            if flux_bounds_fromFile[rxn_id][1] != None:
                rxn.flux_bounds[1] = flux_bounds_fromFile[rxn_id][1]
 
    #--- Next set the flux bounds using flux_bounds and other optimal input arguments ---
    flux_bounds_fromInput = flux_bounds.items()
    argnames = additional_args.keys()
    argvals = additional_args.values()
    for argname in argnames:
      exec "flux_bounds_fromInput += additional_args['" + argname + "'].items()"

    flux_bounds_fromInput = dict(flux_bounds_fromInput)
    for rxn_id in flux_bounds_fromInput.keys(): 
        rxn = model.reactions_by_id[rxn_id]
        if rxn == None:
            raise ValueError('Reaction ' + rxn_id + ' was not found in the model.')
        else:
            if rxn.flux_bounds == []:
                rxn.assign_flux_bounds()
            if flux_bounds_fromInput[rxn_id][0] != None:
                rxn.flux_bounds[0] = flux_bounds_fromInput[rxn_id][0]
            if flux_bounds_fromInput[rxn_id][1] != None:
                rxn.flux_bounds[1] = flux_bounds_fromInput[rxn_id][1]
 
