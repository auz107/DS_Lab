import sys
sys.path.append('/data/alizom/')
from tools.io.read_sbml_model import read_sbml_model
from tools.io.read_gams_model import read_gams_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric


if __name__ == "__main__":

    # Model path
    model_path = '/data/alizom/models/Escherichia_coli/iAF1260/'

    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    #-- Read the sbml file --
    # iAF1260_udpated.xml: Contains the following updates: (1) names of essential reactions were udpated
    # These transport reactions for amino acids were made reversible:
    # ALAt2pp, ALAt4pp, ASPt2_2pp, ASPt2pp, GLNabcpp, GLYt2pp, GLYt4pp, METabcpp
    iAF1260 = read_sbml_model(file_name = model_path + 'iAF1260_updated.xml', model_id = 'iAF1260',model_organism = model_organism,model_type = 'metabolic')
    #iAF1260 = read_sbml_model(file_name = model_path + 'iAF1260_fromGAMS.xml', model_id = 'iAF1260',model_organism = model_organism,model_type = 'metabolic',import_params = True,use_cobra_parser = True)
    iAF1260g = read_gams_model(file_name = 'iAF1260ModelData.py',model_name = 'iAF1260',model_organism = model_organism,model_type = 'metabolic') 

    """
    print '\nComparing models ...'
    for rxn in iAF1260.reactions:
        gams_rxn = iAF1260g.get_reactions({rxn.id:'id'})
        if gams_rxn != None: 
            if rxn.type != gams_rxn.type:
                print 'rxn: {}\tSBML:{}\tGAMS: {}'.format(rxn.id,rxn.type,gams_rxn.type)
            elif rxn.flux_bounds != gams_rxn.flux_bounds:
                print 'rxn: {}\tSBML:{}\tGAMS: {}'.format(rxn.id,rxn.flux_bounds,gams_rxn.flux_bounds)
        else:
            print 'rxn {} was not found in iAF1260g.'.format(rxn.id)
    print '\nFinished comparing models ...'
    """

    # Set the growth medium
    set_specific_bounds(model = iAF1260,file_name = model_path + 'iAF1260_minimal_glucose_aerobic.py',simulation_condition = 'minimal_glucose_aerobic')
    #set_specific_bounds(model = iAF1260,file_name = model_path + 'iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign and objective function coefficients
    for rxn in iAF1260.reactions:
        #if rxn.id == 'Ec_biomass_iAF1260_WT_59p81M':
        if rxn.id == 'Ec_biomass_iAF1260_core_59p81M':
            rxn.objective_coefficient = 1
        else:
            rxn.objective_coefficient = 0

    # FBA  
    iAF1260.fba(optimization_solver = optimization_solver) 

    #--- Check a mutant ---
    set_specific_bounds(model = iAF1260,file_name = model_path + 'iAF1260_minimal_glucose_aerobic.py',flux_bounds = {'HISTD':[0,0],'EX_his-L(e)':[-0.-5,1000]})
    iAF1260.fba(optimization_solver = optimization_solver) 
