import sys
sys.path.append('/usr2/postdoc/alizom/work/')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
import cobra

def update_model(model_path):
    # Define the organism
    model_organism = organism(id = 'Pputida', name = 'Pseudomonas putida',domain = 'Bacteria', genus = 'Pseudomonas', species = 'putida', strain = 'KT2440')

    # Read the sbml file
    iJP962 = read_sbml_model(file_name = model_path + 'iJP962.xml', model_id = 'iJP962',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    #--- Fix the coefficients of the biomass compounds ---
    # Divide stoichiomteric coefficients of all biomass metabolites by 10000
    biomass_rxn = iJP962.get_reactions({'IR10370':'id'})
    for cmp in biomass_rxn.compounds:
        biomass_rxn.stoichiometry[cmp] = biomass_rxn.stoichiometry[cmp]/1000000

    GAM_rxn = iJP962.get_reactions({'IR10377':'id'})
    for cmp in GAM_rxn.compounds:
        GAM_rxn.stoichiometry[cmp] = GAM_rxn.stoichiometry[cmp]/10000

    # Remove C9324, EC9324 and EX_EC9324 from the model. These removals change he biomass somehow
    #iJP962.remove_reactions([iJP962.get_reactions({'EX_EC9324':'id'})])
    #iJP962.remove_compounds(iJP962.get_compounds({'C9324':'id','EC9324':'id'}).values())

    # export to cobra model
    iJP962_cobra = iJP962.export_to_cobra()

    print '---- Add ATPM and a transport and exchange reaction for Galacturonic Acid to the model ---\n'
    # Load the original model
    #iJP962_cobra = cobra.io.read_sbml_model(model_path + 'iJP962.xml')

    #- First load the SBML file and add a transport and exchange reaction for Galacturonic Acid -
    # Get the metabolite object for Galacturonic Acid
    C0274 = iJP962_cobra.metabolites.get_by_id('C0274')

    #- Create an extracellular compound for Galacturonic Acid - 
    EC0274 = cobra.Metabolite('EC0274', formula='',name='D-Galacturonic Acid', compartment='Extra_organism')

    # Create a transport and an exchange reaction
    EX_EC0274 = cobra.Reaction('EX_EC0274')
    EX_EC0274.name = 'D-Galacturonic Acid exchange'
    EX_EC0274.subsystem = 'exchange'
    EX_EC0274.add_metabolites({EC0274:-1})

    C0274t = cobra.Reaction('C0274t')
    C0274t.name = 'Galacturonic Acid transport'
    C0274t.subsystem = 'transport'
    C0274t.lower_bound = -1000
    C0274t.upper_bound = 1000
    C0274t.add_metabolites({EC0274:-1,C0274:1})

    iJP962_cobra.add_reaction(EX_EC0274)
    iJP962_cobra.add_reaction(C0274t)

    #-- Create an extracellular compound for Phenylethylamine -- 
    # Get the metabolite object for Phenylethylamine 
    C3090 = iJP962_cobra.metabolites.get_by_id('C3090')

    EC3090 = cobra.Metabolite('EC3090', formula='',name='Phenylethylamine', compartment='Extra_organism')

    # Create a transport and an exchange reaction
    EX_EC3090 = cobra.Reaction('EX_EC3090')
    EX_EC3090.name = 'Phenylethylamine exchange'
    EX_EC3090.subsystem = 'exchange'
    EX_EC3090.add_metabolites({EC3090:-1})

    C3090t = cobra.Reaction('C3090t')
    C3090t.name = 'Phenylethylamine transport'
    C3090t.subsystem = 'transport'
    C3090t.lower_bound = -1000
    C3090t.upper_bound = 1000
    C3090t.add_metabolites({EC3090:-1,C3090:1})

    iJP962_cobra.add_reaction(EX_EC3090)
    iJP962_cobra.add_reaction(C3090t)

    #- Add ATP maintenance reaction to the model ATPM	[c] : atp + h2o --> adp + h + pi - 
    # atp: C0002   adp: C0008   h2o: C0001  h: C0065    pi: C0009
    C0002 = iJP962_cobra.metabolites.get_by_id('C0002')
    C0008 = iJP962_cobra.metabolites.get_by_id('C0008')
    C0001 = iJP962_cobra.metabolites.get_by_id('C0001')
    C0065 = iJP962_cobra.metabolites.get_by_id('C0065')
    C0009 = iJP962_cobra.metabolites.get_by_id('C0009')
    ATPM = cobra.Reaction('ATPM')
    ATPM.name = 'NGAM ATP'
    ATPM.subsystem = 'cytosol'
    ATPM.add_metabolites({C0002:-1,C0001:-1,C0008:1,C0065:1,C0009:1})
    iJP962_cobra.add_reaction(ATPM)

    # Export the updated model to a new SBML file
    cobra.io.write_sbml_model(iJP962_cobra,model_path + 'iJP962_updated.xml',use_fbc_package=False)

if __name__ == "__main__":

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Pseudomonas_putida/iJP962/'

    # Update the model
    #update_model(model_path = model_path)

    # Define the organism
    model_organism = organism(id = 'Pputida', name = 'Pseudomonas putida',domain = 'Bacteria', genus = 'Pseudomonas', species = 'putida', strain = 'KT2440')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Read the sbml file
    iJP962 = read_sbml_model(file_name = model_path + 'iJP962_updated.xml', model_id = 'iJP962',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    iJP962.validate()

    # Assign and objective function coefficients
    for rxn in iJP962.reactions:
        rxn.objective_coefficient = 0
    #iJP962.get_reactions({'EX_EC9324':'id'}).objective_coefficient = 0
    iJP962.get_reactions({'IR10370':'id'}).objective_coefficient = 1

    #--- List of carbon sources it should grow on according to experimental data ---
    # A dictionary where the keys are the names of the carbon sources and values are
    # the anem of cooresponding exchange reactions in the model
    carbon_sources = {
    '4_hydroxy_benzoic_acid':'EX_EC0133',
    'D_malic_acid':None,
    'glucose':'EX_EC0027',
    'L_arginine':'EX_EC0050',
    'L_asparagine':'EX_EC0129',
    'L_serine':'EX_EC0053',
    'galacturonic acid':'EX_EC0274',
    'glycyl_L_glutamic_acid':None,
    'itaconic acid':None,
    'phenylalanine':'EX_EC0064',
    'phenylethylamine':'EX_EC3090',
    'putrescine':'EX_EC0116',
    'tween_80':None
    }


    print '\n--- Growth on glucose ---'
    # Set the growth medium
    set_specific_bounds(model = iJP962,flux_bounds = {carbon_sources['glucose']:[-10,1000],'EX_EC0007':[-20,1000]},file_name = model_path + 'iJP962_minimal.py')

    # FBA  
    iJP962.fba(optimization_solver = optimization_solver) 

    #--- Growth on other carbon sources ---
    can_grow = []      # List of carbon sources the model can grow on
    cannot_grow = []   # List of carbon sources the model cannot grow on
    for carbon_src in sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose']):
        print '\n--- Growth on {}  ---'.format(carbon_src)
        set_specific_bounds(model = iJP962,flux_bounds = {carbon_sources[carbon_src]:[-10,1000],'EX_EC0007':[-20,1000]},file_name = model_path + 'iJP962_minimal.py')
        iJP962.fba(optimization_solver = optimization_solver) 
        if iJP962.fba_model.solution['exit_flag'] == 'globallyOptimal':
            can_grow.append(carbon_src)
        else:
            cannot_grow.append(carbon_src)

    print '\nTotal # of carbon sources that can grow on based on experimental data = {}'.format(len(carbon_sources.keys())-1)
    print '\nTotal # of carbon sources that are in the model can grow on based on experimental data = {}'.format(len(sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose'])))
    print '\nTotal # of carbon sourced the model can grow on = {}: {}'.format(len(can_grow),can_grow)
    print '\nTotal # of carbon sourced the model cannot grow on = {}: {}'.format(len(cannot_grow),cannot_grow)

    #----------------------------------------------------------
    print '\n--------- Simulations with KBase updated model ------\n'
    # Read the sbml file
    iJP962_kbase = read_sbml_model(file_name = model_path + 'Pseudomonas_putida_iJP962.5.xml', model_id = 'iJP962',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    iJP962_kbase.validate()

    # Assign and objective function coefficients
    for rxn in iJP962_kbase.reactions:
        rxn.objective_coefficient = 0
    iJP962_kbase.get_reactions({'biomass0':'id'}).objective_coefficient = 1

    #--- List of carbon sources it should grow on according to experimental data ---
    carbon_sources_kbase = {
    '4_hydroxy_benzoic_acid':'EX_cpd00136_e0', 
    'D_malic_acid':'EX_cpd00386_e0',    
    'glucose':'EX_cpd00027_e0',   
    'L_arginine':'EX_cpd00051_e0',
    'L_asparagine':'EX_cpd00132_e0',  
    'L_serine':'EX_cpd00054_e0', 
    'galacturonic acid':'EX_cpd00280_e0',  
    'glycyl_L_glutamic_acid':None,        # cpd11592
    'itaconic acid':None,                 # cpd00380
    'phenylalanine':'EX_cpd00066_e0',   
    'phenylethylamine':'EX_cpd03161_e0', 
    'putrescine':'EX_cpd00118_e0',  
    'tween_80':None                       # cpd13392
    } 


    print '\n--- Growth with no carbon source ---'
    # Set the growth medium
    set_specific_bounds(model = iJP962_kbase,flux_bounds = {'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium_kbase.py')

    # FBA  
    iJP962_kbase.fba(optimization_solver = optimization_solver) 

    print '\n--- Growth on glucose ---'
    # Set the growth medium
    set_specific_bounds(model = iJP962_kbase,flux_bounds = {carbon_sources_kbase['glucose']:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium_kbase.py')

    # FBA  
    iJP962_kbase.fba(optimization_solver = optimization_solver) 

    #--- Growth on other carbon sources ---
    can_grow = []      # List of carbon sources the model can grow on
    cannot_grow = []   # List of carbon sources the model cannot grow on
    for carbon_src in sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose']):
        print '\n--- Growth on {}  ---'.format(carbon_src)
        set_specific_bounds(model = iJP962_kbase,flux_bounds = {carbon_sources_kbase[carbon_src]:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium_kbase.py')
        iJP962_kbase.fba(optimization_solver = optimization_solver) 
        iJP962.fba(optimization_solver = optimization_solver) 
        if iJP962_kbase.fba_model.solution['exit_flag'] == 'globallyOptimal':
            can_grow.append(carbon_src)
        else:
            cannot_grow.append(carbon_src)

    print '\nTotal # of carbon sources that can grow on based on experimental data = {}'.format(len(carbon_sources.keys())-1)
    print '\nTotal # of carbon sources that are in the model can grow on based on experimental data = {}'.format(len(sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose'])))
    print '\nTotal # of carbon sourced the model can grow on = {}: {}'.format(len(can_grow),can_grow)
    print '\nTotal # of carbon sourced the model cannot grow on = {}: {}'.format(len(cannot_grow),cannot_grow)




