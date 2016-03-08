import sys
sys.path.append('../../../')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric


if __name__ == "__main__":

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Escherichia_coli/iJO1366/' 

    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655')

    #-- Read the sbml file --
    # iJO1366_udpated.xml: Contains the following updates: (1) names of exchange reactions were udpated
    # These transport reactions for amino acids were made reversible:
    # ALAt2pp, ALAt4pp, ASPt2_2pp, ASPt2pp, GLNabcpp, GLYt2pp, GLYt4pp, METabcpp. Additionally MOX was made
    # irreversible according to Jonathan Monk's email (07-14-2015)
    iJO1366 = read_sbml_model(file_name = model_path + 'iJO1366_updated.xml', model_id = 'iJO1366',model_organism = model_organism,model_type = 'metabolic',import_params = False)
    #iJO1366 = read_sbml_model(file_name = model_path + 'iJO1366_updated.xml', model_id = 'iJO1366',model_organism = model_organism,model_type = 'metabolic',import_params = True)

    # Assign and objective function coefficients
    for rxn in iJO1366.reactions:
        rxn.objective_coefficient = 0
    #iJO1366.biomass_reaction = iJO1366.get_reactions({'Ec_biomass_iJO1366_WT_53p95M':'id'})
    iJO1366.biomass_reaction = iJO1366.get_reactions({'Ec_biomass_iJO1366_core_53p95M':'id'})
    iJO1366.biomass_reaction.objective_coefficient = 1

    """
    # Rxns that must be off according to iJO1366 paper to get a reasonable flux distirbution
    zero_rxns = iJO1366.get_reactions({'CAT':'id','DHPTDNR':'id','DHPTDNRN':'id','FHL':'id','SPODM':'id','SPODMpp':'id','SUCASPtpp':'id','SUCFUMtpp':'id','SUCMALtpp':'id','SUCTARTtpp':'id'})
    for rxn_id in zero_rxns.keys():
        #zero_rxns[rxn_id].assign_flux_bounds()
        print rxn_id,'\t',zero_rxns[rxn_id].flux_bounds
    MOX = iJO1366.get_reactions({'MOX':'id'})
    #print 'MOX type = ',MOX.type
    #MOX.type = 'irreversible'
    #MOX.assign_flux_bounds()
    print 'MOX type = ',MOX.type,' bounds = ',MOX.flux_bounds
    """

    # FBA  
    print '--- Aerobic minimal conditions ----'
    set_specific_bounds(model = iJO1366,file_name = model_path + 'iJO1366_minimal_glucose_aerobic.py',limiting_nutrients = {'EX_glc(e)':[-10,1000],'EX_o2(e)':[-20,1000]},simulation_condition = 'minimal_glucose_aerobic')
    iJO1366.fba() 

    print '--- Anaerobic minimal conditions ----'
    set_specific_bounds(model = iJO1366,file_name = model_path + 'iJO1366_minimal_glucose_anaerobic.py',limiting_nutrients = {'EX_glc(e)':[-10,1000]},simulation_condition = 'minimal_glucose_anaerobic')
    iJO1366.fba() 
