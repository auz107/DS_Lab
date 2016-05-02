import sys
sys.path.append('../../../')
from tools.userError import userError
from tools.globalVariables import *
from tools.io.read_sbml_model import read_sbml_model
from tools.io.read_pydict_model import read_pydict_model
from tools.io.create_superModel import create_superModel
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.create_model import create_model
from tools.ancillary.get_ModelSEED_ids import *
from imp import load_source
from multiprocessing import Process, Manager, Value, Array
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

def add_ATPM_rxn(model_name):
    """
    Add ATPM reaction to the cobra model
    """
    # Load the original model
    cobra_model = cobra.io.read_sbml_model(model_name)

    #- Add ATP maintenance reaction to the model ATPM   [c] : atp + h2o --> adp + h + pi - 
    # atp: cpd00002   h2o: cpd00001   adp: cpd00008   h: cpd00067    pi: cpd00012
    cpd00002 = cobra_model.metabolites.get_by_id('cpd00002_c0')
    cpd00001 = cobra_model.metabolites.get_by_id('cpd00001_c0')
    cpd00008 = cobra_model.metabolites.get_by_id('cpd00008_c0')
    cpd00067 = cobra_model.metabolites.get_by_id('cpd00067_c0')
    cpd00012 = cobra_model.metabolites.get_by_id('cpd00012_c0')
    ATPM = cobra.Reaction('ATPM_c0')
    ATPM.name = 'NGAM ATP'
    ATPM.subsystem = 'cytosol'
    ATPM.add_metabolites({cpd00002:-1,cpd00001:-1,cpd00008:1,cpd00067:1,cpd00012:1})
    cobra_model.add_reaction(ATPM)

    # Export the updated model to a new SBML file
    cobra.io.write_sbml_model(cobra_model,model_path + 'Pf_fromGenome_PpBiomass_gapfilled_ATPM_v6.xml',use_fbc_package=False)

def make_manual_fixes(kbase_model):
    kbase_model.get_reactions({'rxn01507_c0':'id'}).flux_bounds = [0,0]

    #--- Make some manual fixes ------
    # cpd10515 (Fe2+)
    kbase_model.get_reactions({'rxn02056_c0':'id'}).flux_bounds = [-1000,0]  # (2) H+ + Siroheme <=> Sirohydrochlorin + Fe2+
    kbase_model.get_reactions({'rxn00224_c0':'id'}).flux_bounds = [0,1000]   # Protoporphyrin + Fe2+ <=> Heme + (2) H+ (Make irreversible based on iAF1260 (FCLT))

    # cpd10516 (Fe3)
    kbase_model.get_reactions({'rxn05557_c0':'id'}).flux_bounds = [0,0]   # H2O + ATP + (2) Citrate[e] + fe3[e] => ADP + Phosphate + H+ + (2) Citrate + fe3
    kbase_model.get_reactions({'rxn00056_c0':'id'}).flux_bounds = [0,1000]  # O2 + (4) H+ + (4) Fe2+ <=> (2) H2O + (4) fe3 (Make irreversible based on iAF1260)

    # These reactions are shown as irreversible in ModelSeed but were reversible in the model
    #kbase_model.get_reactions({'rxn00711_c0':'id'}).flux_bounds = [-1000,0]
    #kbase_model.get_reactions({'rxn01507_c0':'id'}) = [0,1000]
    kbase_model.get_reactions({'rxn00711_c0':'id'}).flux_bounds = [-1000,0]

def fva(model,file_name,stdout_msgs = True, warnings = True):
    """
    Performs flux variability analysis
    INPUTS:
    -------
        model: Metabolic model
    file_name: file name to store the results
    """
    with open(file_name,'w') as f:
        f.write('fva_flux_bounds = {\n')

    # Original objective funciton coefficients of the reactions
    orig_obj_coeff = {}
    for rxn in model.reactions:
        orig_obj_coeff[rxn] = rxn.objective_coefficient

    # First find out which reacitons can reach their wide lower and upper bound
    for rxn in [r for r in model.reactions if r.type.lower() != 'exchange']:
        for r in model.reactions:
            r.objective_coefficient = 0
        rxn.objective_coefficient = 1

        model.fba(build_new_optModel = False, stdout_msgs = stdout_msgs,warnings = warnings)
        if model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            UB = model.fba_model.solution['objective_value']
        else:
            UB = None

        rxn.objective_coefficient = -1
        model.fba(build_new_optModel = False, stdout_msgs = stdout_msgs,warnings = warnings)
        if model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            LB = -model.fba_model.solution['objective_value']
        else:
            LB = None

        with open(file_name,'a') as f:
            f.write("'" + rxn.id + "':" + str([LB,UB]) + ',\n')

    with open(file_name,'a') as f:
        f.write('}')

    # Assign and objective function coefficients
    for rxn in model.reactions:
        rxn.objective_coefficient = orig_obj_coeff[rxn]

if __name__ == "__main__":
 
    # Path to the model
    model_path = '/usr2/postdoc/alizom/work/models/Pseudomonas_fluorescens/kbase/'

    # Add ATPM reaction to the model
    #add_ATPM_rxn(model_id = model_path + 'Pf_fromGenome_PpBiomass_gapfilled_v5.xml')

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    model_organism = organism(id = 'Pfluorescens', name = 'Pseudomonas fluorescens',domain = 'Bacteria', genus = 'Pseudomonas', species = 'fluorescens', strain = 'SBW25')

    #kbase_model = read_sbml_model(file_name = model_path + 'Pf_fromGenome_aeroGluMinimal_gapfilled_v3.xml', model_id = 'Pf_kbase_model',model_organism = model_organism,model_type = 'metabolic',import_params = False)
    #kbase_model = read_sbml_model(file_name = model_path + 'Pf_fromGenome_aeroLSerineMinimal_gapfilled_v4.xml', model_id = 'Pf_kbase_model',model_organism = model_organism,model_type = 'metabolic',import_params = False)
    kbase_model = read_sbml_model(file_name = model_path + 'Pf_fromGenome_PpBiomass_gapfilled_v5.xml', model_id = 'Pf_kbase_model',model_organism = model_organism,model_type = 'metabolic',import_params = False)
    #kbase_model = read_sbml_model(file_name = model_path + 'Pf_fromGenome_PpBiomass_gapfilled_ATPM_v6.xml', model_id = 'Pf_kbase_model',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Biomass reaction
    kbase_model.biomass_reaction = kbase_model.get_reactions({'biomass0':'id'})    

    # Set the objective function coefficients
    for rxn in kbase_model.reactions:
        rxn.objective_coefficient = 0
    kbase_model.biomass_reaction.objective_coefficient = 1

    carbon_sources = {
    'glucose':'EX_cpd00027_e0',
    'L-Asparagine':'EX_cpd00132_e0',
    'L-serine':'EX_cpd00054_e0',
    'Galacturonic Acid':'EX_cpd00280_e0'
    }

    print '\n--- Growth on glucose ---'
    # Set the growth medium
    set_specific_bounds(model = kbase_model,flux_bounds = {carbon_sources['glucose']:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium.py')

    # FBA  
    kbase_model.fba(optimization_solver = optimization_solver)

    #--- Growth on other carbon sources ---
    for carbon_src in sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose']):
        print '\n--- Growth on {}  ---'.format(carbon_src)
        set_specific_bounds(model = kbase_model,flux_bounds = {carbon_sources[carbon_src]:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium.py')
        kbase_model.fba(optimization_solver = optimization_solver)

    #----------------------------------------------------------
    print '\n--------- Simulations with KBase gapfilled model ------\n'
    # Read the sbml file
    Pf_kbase = read_sbml_model(file_name = model_path + 'Pseudomonas_fluorescens.11.xml', model_id = 'Pf_kbase',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    Pf_kbase.validate()

    # Assign and objective function coefficients
    for rxn in Pf_kbase.reactions:
        rxn.objective_coefficient = 0
    Pf_kbase.get_reactions({'biomass0':'id'}).objective_coefficient = 1

    #--- List of carbon sources it should grow on according to experimental data ---
    carbon_sources = {}
    carbon_sources = {
    #'4_hydroxy_benzoic_acid':'EX_cpd00136_e0',
    #'D_malic_acid':'EX_cpd00386_e0',
    'D-mannitol':'EX_cpd00314_e0',
    'glucose':'EX_cpd00027_e0',
    'L_arginine':'EX_cpd00051_e0',
    #'L_asparagine':'EX_cpd00132_e0',
    'L_serine':'EX_cpd00054_e0',
    'N-acetyl-D-glucosamine':'EX_cpd00122_e0',
    #'phenylethylamine':'EX_cpd03161_e0',
    'putrescine':'EX_cpd00118_e0',
    'tween_80':None                       # cpd13392
    }
    print '\n--- Growth with no carbon source ---'
    # Set the growth medium
    set_specific_bounds(model = Pf_kbase,flux_bounds = {'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium.py')

    # FBA  
    Pf_kbase.fba(optimization_solver = optimization_solver)

    print '\n--- Growth on glucose ---'
    # Set the growth medium
    set_specific_bounds(model = Pf_kbase,flux_bounds = {carbon_sources['glucose']:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium.py')

    # FBA  
    Pf_kbase.fba(optimization_solver = optimization_solver)

    #--- Growth on other carbon sources ---
    for carbon_src in sorted([c for c in carbon_sources.keys() if carbon_sources[c] != None and c != 'glucose']):
        print '\n--- Growth on {}  ---'.format(carbon_src)
        set_specific_bounds(model = Pf_kbase,flux_bounds = {carbon_sources[carbon_src]:[-10,1000],'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium.py')
        Pf_kbase.fba(optimization_solver = optimization_solver)


def break_pf_cycles():
    """
    Checks the code for breaking cycles
    """
    from tools.gap_filling.break_cycles import break_cycles

    # Path to the model
    model_path = '/usr2/postdoc/alizom/work/models/Pseudomonas_fluorescens/kbase/'

    model_organism = organism(id = 'Pfluorescens', name = 'Pseudomonas fluorescens',domain = 'Bacteria', genus = 'Pseudomonas', species = 'fluorescens', strain = 'SBW25')

    # Read the sbml file
    Pf_kbase = read_sbml_model(file_name = model_path + 'Pseudomonas_fluorescens.11.xml', model_id = 'Pf_kbase',model_organism = model_organism,model_type = 'metabolic',import_params = False)

    # Check the model after the modifications 
    Pf_kbase.validate()

    # Assign and objective function coefficients
    for rxn in Pf_kbase.reactions:
        rxn.objective_coefficient = 0
    Pf_kbase.get_reactions({'biomass0':'id'}).objective_coefficient = 1

    print '\n--- Growth with no carbon source ---'
    # Set the growth medium
    set_specific_bounds(model = Pf_kbase,flux_bounds = {'EX_cpd00007_e0':[-20,1000]},file_name = model_path + 'minimal_medium.py')

    # FBA  
    Pf_kbase.fba()

    #---- Check cycle breaking code -----
    bc_inst = break_cycles(model = Pf_kbase)

    # Break cycles test
    test_primal_dual = False
    if test_primal_dual:
        print '\n--- Primal ---'
        bc_inst.run(optModel_name = 'primal')

        print '\n--- Dual ---'
        bc_inst.run(optModel_name = 'dual')

        print '\n--- Bilevel ---'
        solution = bc_inst.run(optModel_name = 'bilevel', max_modification_num = 0, max_biomass_thr = 1e-6)


    # Perform flux variablity analysis
    print '\n Performing FVA ...',
    file_name = 'fva_noCarbonSrc_results.py'
    #fva(model = Pf_kbase, file_name = 'fva_noCarbonSrc_results.py',stdout_msgs = False, warnings = False)
    print 'Fone with FVA!\n'

    # Load the FVA results
    from imp import load_source
    import os
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

    fva_flux_bounds = dataFile.fva_flux_bounds
    cycle_counter = 0
    for rxn_id in fva_flux_bounds.keys():
        rxn = Pf_kbase.get_reactions({rxn_id:'id'})
        LB = fva_flux_bounds[rxn_id][0]
        UB = fva_flux_bounds[rxn_id][1]
        if LB == -1000:
            rxn.ignoreLB = False
        else:
            rxn.ignoreLB = True
    
        if UB == 1000:
            rxn.ignoreUB = False
        else:
            rxn.ignoreUB = True
            
        if abs(LB) < 1e-6 and abs(UB) <= 1e-6:
            rxn.blocked = True

        if LB == -1000 or UB == 1000:
            cycle_counter += 1

    print '\nTotal # of reactions participating in cycles = {}\n'.format(cycle_counter)

    solutions = bc_inst.run(optModel_name = 'bilevel', max_modification_num = 1, max_solution_num = 300, max_biomass_thr = 1e-6)
    counter = 0
    print '\n----- Solutions -----\n'
    for soln in solutions:
        counter += 1
        print '\n--- Solution # {}, # of modificaitons = {} , total penalites = {}'.format(counter,soln['modifications_num'],soln['objective_value'])
        print '\nzero_yLBopt_rxns:'
        if len(soln['zero_yLBopt_rxns']) == 0:
            print 'None'
        else:
            for rxn in soln['zero_yLBopt_rxns']:
                print '{} ({}): {}\t{}'.format(rxn.id,rxn.name,rxn.get_equation(ref_type = 'name'),rxn.flux_bounds)

        print '\nzero_yUBopt_rxns:'
        if len(soln['zero_yUBopt_rxns']) == 0:
            print 'None'
        else:
            for rxn in soln['zero_yUBopt_rxns']:
                print '{} ({}): {}\t{}'.format(rxn.id,rxn.name,rxn.get_equation(ref_type = 'name'),rxn.flux_bounds)


