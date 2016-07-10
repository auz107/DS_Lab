from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../../')
from copy import deepcopy
from tools.userError import userError
from tools.globalVariables import *
from tools.io.read_sbml_model import read_sbml_model
from tools.io.read_pydict_model import read_pydict_model
from tools.gap_filling.create_superModel import create_superModel_from_ModelSEED
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.io.create_model import create_model
from tools.utilities.get_ModelSEED_ids import *
from imp import load_source
from multiprocessing import Process, Manager, Value, Array
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

def test_create_superModel(warnings = True, stdout_msgs = True):
    """
    Creates the updated model after adding external reactions for the production of nsAAs

    Ali R. Zomorrodi - Segre's Lab @ BU
    Last updated: 03-17-2016
    """
    # Model path
    model_path = home_dir + 'work/models/Escherichia_coli/iJO1366/'

    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    flux_bounds_dict = {'EX_glc(e)':[-100,1000], 'EX_o2(e)':[-200,1000]}
    flux_bounds_filename = model_path + 'iJO1366_minimal_glucose_aerobic.py'

    # Orignal iJo1266 model
    model = create_model(model_organism = model_organism, model_info = {'id':'iJO1366', 'file_format':'sbml', 'model_filename':model_path + 'iJO1366_updated.xml', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, validate = True, stdout_msgs = True, warnings = True) 

    model.organism.ModelSEED_type = 'bacteria_GramNegative'

    """
    #rxn_id = 'GAPD'
    #rxn_id = 'ZN2abcpp'
    rxn_id = 'VALtex'
    get_cpds_ModelSEED_id(cpds_list = model.reactions_by_id[rxn_id].compounds)
    for cpd in model.reactions_by_id[rxn_id].compounds:
        print cpd.id,': ',cpd.ModelSEED_id,'   stoic = ',model.reactions_by_id[rxn_id].stoichiometry[cpd]
    rxn_ModelSEED_id, ModelSEED_id_found_by = match_rxn_eqn(rxn = model.reactions_by_id[rxn_id])    
    print 'rxn_ModelSEED_id = {}  , ModelSEED_id_found_by = {}\n'.format(rxn_ModelSEED_id, ModelSEED_id_found_by)
    """

    print '\n----- Getting ModelSEED ids ----'
    get_cpds_ModelSEED_id(cpds_list = model.compounds)
    get_rxns_ModelSEED_id(rxns_list = model.reactions)
    sys.stdout.flush()

    """
    #--- Text exporiting to pydict and importing back ----
    # Save the model into a python dictionary
    model.export(output_format = 'pydict', output_filename = 'iJO1366.py')

    # Re-import the pydict model
    model = create_model(model_info = {'id':'iJO1366', 'file_format':'pydict', 'model_filename':'iJO1366.py', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': flux_bounds_dict}, validate = True, stdout_msgs = True, warnings = True) 
    """

    print '\n----- Creating super_model ----'
    sys.stdout.flush()

    super_model = create_superModel_from_ModelSEED(original_model = model, standard_to_model_compartID_map = {'c':'c','e':'e','p':'p'}, validate = True)
    print 'super_model statistics: # of compounds (total/original/external) = {}/{}/{} ,  # of reactions (total/original/external) = {}/{}/{}'.format(len(super_model.compounds), len([c for c in super_model.compounds if not c.external]), len([c for c in super_model.compounds if c.external]), len(super_model.reactions), len([r for r in super_model.reactions if not r.external]), len([r for r in super_model.reactions if r.external]))    

    print '\n----- fba with super_model ----'
    for rxn in [r for r in super_model.reactions if r.external]:
        rxn.flux_bounds = [0,0]
    for rxn in super_model.reactions:
        rxn.objective_coefficient = 0
    super_model.reactions_by_id['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 1
    super_model.fba()

    print '\n----- Exporting super_model to pydict ----'
    # Export to a pydict model
    super_model.export(output_format = 'pydict', output_filename = 'super_model_iJO1366.py')

    print '\n----- Re-importing super_model from pydict ----'
    # Re-import the pydict model
    super_modelmodel = create_model(model_info = {'id':'super_model_iJO1366', 'file_format':'pydict', 'model_filename':'super_model_iJO1366.py', 'biomassrxn_id':'Ec_biomass_iJO1366_core_53p95M'}, growthMedium_flux_bounds = {'flux_bounds_filename':flux_bounds_filename, 'flux_bounds_dict': dict(flux_bounds_dict.items() + [(r.id,[0,0]) for r in super_model.reactions if r.external])}, validate = True, stdout_msgs = True, warnings = True) 


#------------------------------------
if __name__ == '__main__':
    test_create_superModel(warnings = True, stdout_msgs = True)



