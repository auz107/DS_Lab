from __future__ import division
import sys, time
sys.path.append('../')
from tools.io.read_sbml_model import read_sbml_model
from tools.core.organism import organism
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.compartment import compartment
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.GAMETES.game import game
from tools.userError import userError
from coopr.pyomo import *
from coopr.opt import *
from copy import deepcopy
from multiprocessing import Process, Manager, Value, Array
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

def create_model(stdout_msgs = True, o2_uptake_rate = 2):
    """
    Creates the metabolic model
    """
    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/'

    # Define the organism (see function findHisSatConc for the calculations of gDW_per_cell for yeast)
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '',gDW_per_cell = 4e-11)

    # Load the original model
    iAZ900 = read_sbml_model(file_name = model_path + 'iAZ900_noCycles_03_25_2015.xml', model_id = 'iAZ900', model_name = 'iAZ900 S. cerevisiae model',model_organism = model_organism, model_type = 'metabolic',import_params = False, stdout_msgs = stdout_msgs)

    # Set the objective function
    iAZ900.biomass_reaction = iAZ900.get_reactions({'biomass_core':'id'})
    for rxn in iAZ900.reactions:
        rxn.objective_coefficient = 0
    #biomass_rxn = iAZ900.get_reactions({'biomass_wildType':'id'})
    biomass_rxn = iAZ900.get_reactions({'biomass_core':'id'})
    biomass_rxn.objective_coefficient = 1


    if stdout_msgs:
        print '\n--- FBA for original iAZ900 ---'
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
    iAZ900.fba(stdout_msgs = stdout_msgs)

    return iAZ900

def test_model(stdout_msgs = True, o2_uptake_rate = 2):
    """
    Creates the metabolic model
    """
    model_rxnmap = {}
    model_id = 'iSce926'
    if model_id == 'iAZ900':
        model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/'
        model_filename = 'test.xml'
        model_media_filename = 'iAZ900_minimal.py'
        model_rxnmap['biomass'] = 'biomass_core'
        model_rxnmap['EX_sucr_e_'] = 'EX_sucr_e_'
        model_rxnmap['SUCRe'] = 'SUCRe'
        model_rxnmap['EX_glc_e_'] = 'EX_glc_e_'
        model_rxnmap['GLCt1'] = 'GLCt1'
        model_rxnmap['EX_fru_e_'] = 'EX_fru_e_'
        model_rxnmap['FRUt2'] = 'FRUt2'
        model_rxnmap['EX_o2_e_'] = 'EX_o2_e_'
        model_rxnmap['ATPM'] = 'ATPM'

    elif model_id == 'iSce926':
        model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iSce926/'
        #model_filename = 'iSce926_Jan11_2016.xml'
        model_filename = 'test.xml'
        model_biomass_rxnname = 'r_2133'
        model_media_filename = 'iSce926_minimal.py'
        model_rxnmap['biomass'] = 'r_2133'
        model_rxnmap['EX_sucr_e_'] = 'r_2058'
        model_rxnmap['SUCRe'] = 'r_1024'
        model_rxnmap['EX_glc_e_'] = 'r_1714'
        model_rxnmap['GLCt1'] = 'r_1166'
        model_rxnmap['EX_fru_e_'] = 'r_1709'
        model_rxnmap['FRUt2'] = 'r_1134'
        model_rxnmap['EX_o2_e_'] = 'r_1992'
        model_rxnmap['ATPM'] = 'r_4053'

    print 'Setting bounds ...'
    # Define the organism (see function findHisSatConc for the calculations of gDW_per_cell for yeast)
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '',gDW_per_cell = 4e-11)

    # Load the original model
    model = read_sbml_model(file_name = model_path + model_filename, model_id = model_id, model_name = model_id + ' S. cerevisiae model',model_organism = model_organism, model_type = 'metabolic',import_params = False, stdout_msgs = stdout_msgs)

    # Set the objective function
    model.biomass_reaction = model.reactions_by_id[model_rxnmap['biomass']]
    for rxn in model.reactions:
        rxn.objective_coefficient = 0
    model.biomass_reaction.objective_coefficient = 1

    print '\n--- FBA for original iAZ900 ---'
    #SUCRe_flux = 1.12213158334
    SUCRe_flux = 1.58002065773 
    #set_specific_bounds(model = model,file_name = model_path + model_media_filename,flux_bounds = {model_rxnmap['EX_sucr_e_']:[-10,1000],model_rxnmap['EX_o2_e_']:[-o2_uptake_rate,1000]})
    #set_specific_bounds(model = model,file_name = model_path + model_media_filename,flux_bounds = {model_rxnmap['EX_sucr_e_']:[-10,1000],model_rxnmap['EX_o2_e_']:[-o2_uptake_rate,1000],model_rxnmap['EX_glc_e_']:[-0.7*SUCRe_flux,1000],model_rxnmap['EX_fru_e_']:[-0.7*SUCRe_flux,1000],model_rxnmap['SUCRe']:[0,1000]})
    set_specific_bounds(model = model,file_name = model_path + model_media_filename,flux_bounds = {model_rxnmap['EX_sucr_e_']:[-10,1000],model_rxnmap['EX_o2_e_']:[-o2_uptake_rate,1000],model_rxnmap['EX_glc_e_']:[-0.7*SUCRe_flux,1000],model_rxnmap['EX_fru_e_']:[-0.7*SUCRe_flux,1000],model_rxnmap['SUCRe']:[SUCRe_flux,1000]})
    #set_specific_bounds(model = model,file_name = model_path + model_media_filename,flux_bounds = {model_rxnmap['EX_sucr_e_']:[-10,1000],model_rxnmap['EX_o2_e_']:[-o2_uptake_rate,1000],model_rxnmap['EX_glc_e_']:[-0.7*SUCRe_flux,1000],model_rxnmap['EX_fru_e_']:[-0.7*SUCRe_flux,1000],model_rxnmap['SUCRe']:[SUCRe_flux,1000], model_rxnmap['ATPM']:[1 + SUCRe_flux,1000]})

    model.fba(stdout_msgs = stdout_msgs)

    print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {} ,  FRUt2 = {}\n'.format(model.fba_model.solution['opt_rxnFluxes'][model_rxnmap['EX_sucr_e_']], model.fba_model.solution['opt_rxnFluxes'][model_rxnmap['SUCRe']],model.fba_model.solution['opt_rxnFluxes'][model_rxnmap['EX_glc_e_']],model.fba_model.solution['opt_rxnFluxes'][model_rxnmap['GLCt1']], model.fba_model.solution['opt_rxnFluxes'][model_rxnmap['EX_fru_e_']],model.fba_model.solution['opt_rxnFluxes'][model_rxnmap['FRUt2']])


def create_maxSUCRe_fbaModel(model,stdout_msgs = True):
    """
    Creates an FBA model for maximizing SUCRe reaction flux
    """
    # Find reactions participating in the objective function
    obj_rxns = [r for r in model.reactions if r.objective_coefficient != 0]

    # Save the objective coefficients for these reactions in a dictionary
    obj_coeff = dict([(r,r.objective_coefficient) for r in obj_rxns])

    # Set the objective coefficient for all reactions to zero
    for r in model.reactions:
        r.objective_coefficient = 0

    model.reactions_by_id['SUCRe'].objective_coefficient = 1

    atpm_rxn_bounds = model.reactions_by_id['ATPM'].flux_bounds
    model.reactions_by_id['ATPM'].flux_bounds = [0,1000]

    model.max_SUCRe_fbaModel = fba(model = model, build_new_optModel = True, store_opt_fluxes = False, stdout_msgs = False)
    model.max_SUCRe_fbaModel.run()

    # Set the objective coefficients back to what they were before
    model.reactions_by_id['SUCRe'].objective_coefficient = 0
    for r in obj_rxns:
        r.objective_coefficient = obj_coeff[r]
    model.reactions_by_id['ATPM'].flux_bounds = atpm_rxn_bounds

    return model

def gameCreator(input_data,output_data):

    #--- Creating a shared memory using the manager ---
    model_path = input_data['model_path']
    iAZ900 = input_data['iAZ900']
    capture_eff = input_data['capture_eff']
    atp_coeff = input_data['atp_coeff']
    his_uptake_rate = input_data['his_uptake_rate']
    glc_uptake_rate = input_data['glc_uptake_rate']
    fru_uptake_rate = input_data['fru_uptake_rate']
    o2_uptake_rate = input_data['o2_uptake_rate']
    simulate_his = input_data['simulate_his']
    fixed_o2 = input_data['fixed_o2']
    death_rate = input_data['death_rate']
    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings']
    results_filename = input_data['results_filename']

    if stdout_msgs:
        if stdout_msgs:
            print '\n--- FBA for iAZ900 before adding atp/adp to SUCRe---'

        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs)

    #--- Incorporate atp into SUCRe ---
    # According to Jeff Gore's paper (suppl. figure 5) the cooperators are 2.5% less fit than cheasters
    # Note that the plot in this figure is for ln(X/X0)_C/ln(X/X0)_D but ln(X/X0) = mu*t or ln(X/X0)
    # is proportional to mu.
    # In the absence of a metabolic cost for sucrose hydrolysis (SUCRe) the growth rate for WT is 0.5031.
    # We next add an ATP cost to SUCRe reaction (for WT) and adjust its coefficinet such that the growth
    # becomes 97.5% of 0.5023 = 0.49056. The coefficient of atp was determined to be 0.115
    if atp_coeff > 0:
        SUCRe = iAZ900.get_reactions({'SUCRe':'id'})
        atp = iAZ900.get_compounds({'atp_c':'id'})
        adp = iAZ900.get_compounds({'adp_c':'id'})
        pi = iAZ900.get_compounds({'pi_c':'id'})
        SUCRe.add_compounds({atp:-atp_coeff,adp:atp_coeff,pi:atp_coeff})

        if stdout_msgs:
            print '\n--- FBA for iAZ900 after adding atp/adp to SUCRe---'
            print 'SUCRe:\t',SUCRe.get_equation('id')

            if iAZ900.fixed_o2:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
            iAZ900.fba(stdout_msgs = stdout_msgs)
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])
            if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
                payoff_C_capeff_one = iAZ900.fba_model.solution['objective_value'] 

    #--- Find the metabolic cost of sucrose production and secretion ---
    # invertase_cost = -[(max biomass before the addition of ATP to SUCRe) - (max biomass after the additon of ATP to SUCRe)]
    if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
        # Max biomass before the addition of ATP
        invertase_cost = -(iAZ900.noATPM_max_biomass - iAZ900.fba_model.solution['objective_value'])
    elif iAZ900.noATPM_max_biomass > 0:
        invertase_cost = -iAZ900.noATPM_max_biomass 
    elif iAZ900.noATPM_max_biomass == 0:
        invertase_cost = - atp_coeff
    # Make sure intertase cost <= 0
    if invertase_cost > 0:
        raise userError('Positive invertase cost: ' + str(invertase_cost)) 
    if stdout_msgs:
        print '\nInvertasee production cost = {}\n'.format(invertase_cost)

    #--- Incorporate capture efficiency ---
    # In order to incorporate the capture efficiency into the model, change the SUCRe equation from 
    # SUCRe:   1.0 h2o_e + 1.0 sucr_e --> 1.0 fru_e + 1.0 glc_D_e
    # to (for a capture efficiency of 0.01):
    #          1.0 h2o_e + 1.0 sucr_e --> 0.01 fru_e + 0.01 glc_D_e + 0.99 fru_secreted + 0.99 glc_secreted
    # and then add two reactions exclusively for the secreion of glucose and fructose
    # EX_glc_secreted: glc_D_secreted --> and EX_fru_secreted: fru_secreted -->
    SUCRe = iAZ900.get_reactions({'SUCRe':'id'})
    glc_D_e = iAZ900.get_compounds({'glc_D_e':'id'})
    fru_e = iAZ900.get_compounds({'fru_e':'id'})
    if capture_eff != 0:
        #SUCRe.stoichiometry[glc_D_e] = capture_eff
        #SUCRe.stoichiometry[fru_e] = capture_eff
        SUCRe.set_stoichiometry({glc_D_e:capture_eff, fru_e:capture_eff}, replace = False, model = iAZ900)
    else:
        #SUCRe.del_compounds([glc_D_e,fru_e])
        SUCRe.del_compounds([glc_D_e,fru_e], model = iAZ900)
    glc_D_secreted = compound(id = 'glc_D_secreted',name = 'D-Glucose',compartment = glc_D_e.compartment)
    fru_secreted = compound(id = 'fru_secreted',name = 'D-Fructose',compartment = fru_e.compartment)
    if capture_eff < 1:
        #SUCRe.add_compounds({glc_D_secreted:1 - capture_eff,fru_secreted:1 - capture_eff})
        SUCRe.add_compounds({glc_D_secreted:1 - capture_eff,fru_secreted:1 - capture_eff}, model = iAZ900)
    EX_glc_secreted = reaction(id = 'EX_glc_secreted', name = 'D-Glucose exclusive export',stoichiometry = {glc_D_secreted:-1}, reversibility = 'reversible', is_exchange = True)
    EX_glc_secreted.objective_coefficient = 0
    EX_fru_secreted = reaction(id = 'EX_fru_secreted', name = 'D-Fructose exclusive export',stoichiometry = {fru_secreted:-1}, reversibility = 'reversible', is_exchange = True)
    EX_fru_secreted.objective_coefficient = 0
    iAZ900.add_reactions([EX_glc_secreted,EX_fru_secreted])
    iAZ900.set_cpds_genes_rxns(do_cpds = True, do_gens = False)
    iAZ900.validate(stdout_msgs = stdout_msgs)
    if stdout_msgs:
        print '\n--- FBA for iAZ900 after incorporating the capture efficiency --'
        print 'SUCRe:\t',SUCRe.get_equation('id')
        print 'SUCRe:\t',EX_glc_secreted.get_equation('id')
        print 'SUCRe:\t',EX_fru_secreted.get_equation('id')

        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000],'GLCt1':[3,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  , EX_glc_secreted = {} ,  EX_fru_e_ = {}  ,  EX_fru_secreted = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted'])

    #--- Cheater and cooperator strains ---
    # Cooperator has a defective HIS3 (YOR202W) gene (for histidone production) corresponding to
    # reaction IGPDH in the model
    # Cheater is missing the gene SUC2 (YIL162W) corresponding ot reaction SUCRe
    Cooperator = deepcopy(iAZ900)
    Cooperator.id = 'Cooperator'

    Cheater = deepcopy(iAZ900)
    Cheater.id = 'Cheater'

    #--- Create the players' strategies ---
    players_names = ['player_1','player_2']
    players_strategies = {}
    players_strategies['player_1'] = ['Cooperate','Defect']
    players_strategies['player_2'] = ['Cooperate','Defect']

    payoff_matrix = {}    

    # --- Find out whether the cooperator strains can cooperate and how much ---
    # A parameter showing whether cooperation (performing SUCRe reaction) is technically possible
    # In some cases even if we relax ATPM constriant SUCRe would not be able to carry any flux as
    # its ATP requirements cannot be satisfied. Under these conditions cooperation is not possible
    can_cooperate = True

    if stdout_msgs:
        print '-- Check whether cooperation is possible --'

        print '\nFBA for Cooperator alone '

    if simulate_his:
        if Cooperator.fixed_o2:
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
    else:
        if Cooperator.fixed_o2:
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})

    Cooperator.fba(stdout_msgs = stdout_msgs)
    if Cooperator.fba_model.solution['exit_flag'] == 'globallyOptimal':
        # Cooperator's payoff
        payoff_C_alone = Cooperator.fba_model.solution['objective_value'] 
        glc_secretion_flux = Cooperator.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted']
        fru_secretion_flux = Cooperator.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted']
        coopr_SUCRe_flux = Cooperator.fba_model.solution['opt_rxnFluxes']['SUCRe']
        if stdout_msgs:
            print '\nCooperator alone:  EX_sucr_e_ = {}, SUCRe = {}  ,  glc_secretion_flux = {} , fru_secretion_flux = {}\n'.format(Cooperator.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], Cooperator.fba_model.solution['opt_rxnFluxes']['SUCRe'],glc_secretion_flux,fru_secretion_flux)
    else:
        # Create fba model for maximizing maxSUCRe
        create_maxSUCRe_fbaModel(model = Cooperator)

        # Update the flux bounds for max_SUCRe_fbaModel
        for j in [rxn for rxn in Cooperator.reactions if rxn.type.lower() == 'exchange' or rxn.id == 'IGPDH']: 
            Cooperator.max_SUCRe_fbaModel.optModel.v[j.id].setlb(j.flux_bounds[0])
            Cooperator.max_SUCRe_fbaModel.optModel.v[j.id].setub(j.flux_bounds[1])
        SUCRe_fba_solution = Cooperator.max_SUCRe_fbaModel.run()
        if SUCRe_fba_solution['exit_flag'] == 'globallyOptimal':
            # Note that both death_rate and invertase_cost are negative
            payoff_C_alone = death_rate + invertase_cost 
            glc_secretion_flux = SUCRe_fba_solution['opt_rxnFluxes']['EX_glc_secreted']
            fru_secretion_flux = SUCRe_fba_solution['opt_rxnFluxes']['EX_fru_secreted']
            coopr_SUCRe_flux = SUCRe_fba_solution['opt_rxnFluxes']['SUCRe'] 
            if stdout_msgs:
                print '\nmax_SUCRe_fbaModel was solved sucessfully with  EX_sucr_e_ = {}, SUCRe = {}  ,  glc_secretion_flux = {} , fru_secretion_flux = {}\n'.format(Cooperator.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], Cooperator.fba_model.solution['opt_rxnFluxes']['SUCRe'],glc_secretion_flux,fru_secretion_flux)
        else:
            # Note that both death_rate and invertase_cost are negative
            payoff_C_alone = death_rate            
            glc_secretion_flux = 0
            fru_secretion_flux = 0
            can_cooperate = False # Cooperation is not possible in this case
            if stdout_msgs:
                print '\nmax_SUCRe_fbaModel was not solved sucessfully. Cooperation is not possible.'

    # If cooperation was not possible, they are both defecting (defecting is the only possible strategy of both players)
    if not can_cooperate:
        players_strategies = {}
        players_strategies['player_1'] = ['Defect']
        players_strategies['player_2'] = ['Defect']
        # payoffs when both defect
        payoff_D = death_rate
        payoff_matrix = {}   
        payoff_matrix[(('player_1','Defect'),('player_2','Defect'))] = {'player_1':payoff_D,'player_2':payoff_D}

    else: # If Cooperation is possible
    
        #-- Cooperator vs. Cooperator --
        if stdout_msgs:
            print '-- Cooperator vs. Cooperator --'

        if simulate_his:
            if Cooperator.fixed_o2:
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000], 'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
            else:
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            if Cooperator.fixed_o2:
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000], 'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
            else:
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
    
        Cooperator.fba(stdout_msgs = stdout_msgs)
        if Cooperator.fba_model.solution['exit_flag'] == 'globallyOptimal':
            # Cooperator's payoff
            payoff_C = Cooperator.fba_model.solution['objective_value'] 

            if stdout_msgs:
                print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  glc_secretion_flux = {} , fru_secretion_flux = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(Cooperator.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], Cooperator.fba_model.solution['opt_rxnFluxes']['SUCRe'],glc_secretion_flux,fru_secretion_flux,Cooperator.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],Cooperator.fba_model.solution['opt_rxnFluxes']['GLCt1'], Cooperator.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])
    
        else:
            if warnings:
                print 'WARNING!FBA for CC was not solved to optimality. Zero was assigned as payoffs'
            payoff_C_alone = death_rate + invertase_cost 

        payoff_matrix[(('player_1','Cooperate'),('player_2','Cooperate'))] = {'player_1':payoff_C,'player_2':payoff_C}
    
        #-- Cooperator vs. Cheater --
        if stdout_msgs:
            print '-- Cooperator vs. Cheater --'
            print '\nglc_secretion_flux = {} , fru_secretion_flux = {}\n'.format(glc_secretion_flux,fru_secretion_flux)

        # Cooperator
        payoff_C = payoff_C_alone

        # Cheater 
        if stdout_msgs:
            print 'FBA for Cheater:'
    
        if Cheater.fixed_o2:
            set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[- glc_uptake_rate - glc_secretion_flux,1000],'EX_fru_e_':[- fru_uptake_rate - fru_secretion_flux,1000],'EX_o2_e_':[-o2_uptake_rate,1000],'SUCRe':[0,0]})
        else:
            set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_uptake_rate - glc_secretion_flux,1000],'EX_fru_e_':[- fru_uptake_rate - fru_secretion_flux,1000],'EX_o2_e_':[(- glc_uptake_rate - fru_uptake_rate - glc_secretion_flux - fru_secretion_flux)/5,1000],'SUCRe':[0,0]})
    
        Cheater.fba(stdout_msgs = stdout_msgs)
        if Cheater.fba_model.solution['exit_flag'] == 'globallyOptimal':
            payoff_D = Cheater.fba_model.solution['objective_value']
            if stdout_msgs:
                print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {}  ,  EX_fru_e_ = {}\n'.format(Cheater.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], Cheater.fba_model.solution['opt_rxnFluxes']['SUCRe'],Cheater.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],Cheater.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])
        else:
            if stdout_msgs:
                print 'WARNING! FBA for Cheater in CD was not solved to optimality. Death rate was assigned as payoff'
            payoff_D = death_rate
    
        payoff_matrix[(('player_1','Cooperate'),('player_2','Defect'))] = {'player_1':payoff_C,'player_2':payoff_D}
        payoff_matrix[(('player_1','Defect'),('player_2','Cooperate'))] = {'player_1':payoff_D,'player_2':payoff_C}
    
        #-- Cheater vs.Cheater --
        if stdout_msgs:
            print '-- Cheater vs. Cheater --'
        # Cheater vs. defector not only cannot they grow but also cannot satisfy their ATPM maintenance. 
        # So their payoffs is equal to their death rate.
        payoff_D = death_rate
        payoff_matrix[(('player_1','Defect'),('player_2','Defect'))] = {'player_1':payoff_D,'player_2':payoff_D}

    # Find Nash equilibria
    suc_game = game(game_name = 'yeast_sucrose' + str(capture_eff), players_names = players_names, players_strategies = players_strategies, payoff_matrix = payoff_matrix)
    suc_game.find_NashEq(stdout_msgs = stdout_msgs)

    output_data['suc_game'] = suc_game

    if results_filename != '':
        with open(results_filename,'a') as f:
            f.write('results[' + str((('capture_eff',capture_eff),('atp_coeff',atp_coeff),('his',his_uptake_rate),('glc',glc_uptake_rate),('fru',fru_uptake_rate))) + '] = ' + str(suc_game.pureNash_equilibria) + '\n')

def master_func(start_pos = None, end_pos = None, capture_efficiency = 0.99, SUCRe_atp = 0.115, histidine_uptake_rate = 10, glucose_uptake_rate = 0, fructose_uptake_rate = 0, o2_uptake_rate = 2, simulate_his = True, fixed_o2 = False, task = 'analyze_games', warnings = True, stdout_msgs = True, results_filename = ''):

    """
    Creates the required inputs for gameCreator

    INPUTS:
    -------
           start_pos: Start position of the array containing all possble cases to consider (see all_cases variable)
             end_pos: End position of the array containing all possble cases to consider (see all_cases variable)
        simulate_his: A parameter showing whether we are intereated in simulating the cooperation cost with histidine (True) 
                      or with atp (False)
            fixed_o2: A parameter showing whether the oxygen uptake rate must be fixed at -2 for all simulations (True) or not (False)    
                task: A string showing what the function needs to perform. Eligible cases are analyze_games (finding the NE of games)
                      and analyze_cost (assessing the cooperation cost)
    results_filename: Name of the output results file name 
    """
    # Task
    if task not in ['analyze_games','analyze_cost']:
        raise ValueError("Invalud value for 'taks'! Allowed values are 'analyze_games' and 'analyze_cost'")

    # If a float is provided convert it to a list
    if isinstance(capture_efficiency,int) or isinstance(capture_efficiency,float):
        capture_efficiency = [capture_efficiency]
    elif not isinstance(capture_efficiency,list):
        raise TypeError('capture_efficiency must be an integer, a float or a list of integers or floats.')

    if isinstance(SUCRe_atp,int) or isinstance(SUCRe_atp,float):
        SUCRe_atp = [SUCRe_atp]
    elif not isinstance(SUCRe_atp,list):
        raise TypeError('SUCRe_atp must be an integer, a float or a list of integers or floats. SUCRe_atp = {}'.format(SUCRe_atp))

    if isinstance(histidine_uptake_rate,int) or isinstance(histidine_uptake_rate,float):
        histidine_uptake_rate = [histidine_uptake_rate]
    elif not isinstance(histidine_uptake_rate,list):
        raise TypeError('histidine_uptake_rate must be an integer, a float or a list of integers or floats.')

    if isinstance(glucose_uptake_rate,int) or isinstance(capture_efficiency,float):
        glucose_uptake_rate = [glucose_uptake_rate]
    elif not isinstance(glucose_uptake_rate,list):
        raise TypeError('glucose_uptake_rate must be an integer, a float or a list of integers or floats.')

    if isinstance(fructose_uptake_rate,int) or isinstance(fructose_uptake_rate,float):
        fructose_uptake_rate = [fructose_uptake_rate]
    elif not isinstance(fructose_uptake_rate,list):
        raise TypeError('fructose_uptake_rate must be an integer, a float or a list of integers or floats.')

    # Check that all elements of SUCRe_atp, histidine_uptake_rate, glucose_uptake and furctose_uptake_rate  are non-negative
    if len([a for a in histidine_uptake_rate if a < 0]) > 0:
        raise ValueError('All elements of histidine_uptake_rate must be non-negative')
    if len([a for a in glucose_uptake_rate if a < 0]) > 0:
        raise ValueError('All elements of glucose_uptake_rate must be non-negative')
    if len([a for a in fructose_uptake_rate if a < 0]) > 0:
        raise ValueError('All elements of fructose_uptake_rate must be non-negative')
    if len([a for a in SUCRe_atp if a < 0]) > 0:
        raise ValueError('All elements of SUCRe_atp must be non-negative')

    # O2 uptake
    if o2_uptake_rate < 0:
        raise ValueError('o2_uptake_rate must be non-negative')

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/' 

    # Create the metabolic model object
    iAZ900 = create_model(stdout_msgs = stdout_msgs, o2_uptake_rate = o2_uptake_rate)
    iAZ900.fixed_o2 = fixed_o2
    if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
        # Max biomass before the addition of ATP
        iAZ900.noATPM_max_biomass = iAZ900.fba_model.solution['objective_value']
    else:
        iAZ900.noATPM_max_biomass = 0 


    #---- Calculate the death rate in the absence of glucose and fructose (no invertase) ---
    # mu_death = Yxs*(v_uptake - ms)
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[0,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
    glc_D_e = iAZ900.get_compounds({'glc_D_e':'id'})
    glc_D_e.ms_calc(stdout_msgs = False, warnings = False, model = iAZ900)
    glc_D_e.biomass_yield_calc(stdout_msgs = False, warnings = False, model = iAZ900)
    mu_death_glc_D = - glc_D_e.biomass_yield*glc_D_e.ms 
    fru_e = iAZ900.get_compounds({'fru_e':'id'})
    fru_e.ms_calc(stdout_msgs = False, warnings = False, model = iAZ900)
    fru_e.biomass_yield_calc(stdout_msgs = False, warnings = False, model = iAZ900)
    mu_death_fru = - fru_e.biomass_yield*fru_e.ms 
    death_rate = min(mu_death_glc_D,mu_death_fru)
    if stdout_msgs:
        print 'death_rate = min(glc_death = {:4f},fru_death = {:4f}) = {:4f}'.format(mu_death_glc_D,mu_death_fru,death_rate)

    if results_filename != '':    
        with open(results_filename,'w') as f:
            if task.lower() == 'analyze_games':
                f.write('results = {}\n')
            elif task.lower() == 'analyze_cost':
                f.write('cooperation_cost = {}\n')

    print '\ncapture_efficiency = ',capture_efficiency
    print '\no2_uptake_rate = ',o2_uptake_rate
    print '\nSUCRe_atp = ',SUCRe_atp
    print '\nhistidine_uptake_rate = ',histidine_uptake_rate
    print '\nglucose_uptake_rate = ',glucose_uptake_rate
    print '\nfructose_uptake_rate = ',fructose_uptake_rate
    all_cases = [(capture_eff,atp_coeff,his_uptake_rate,glc_uptake_rate,fru_uptake_rate) for capture_eff in capture_efficiency for atp_coeff in SUCRe_atp for his_uptake_rate in histidine_uptake_rate for glc_uptake_rate in glucose_uptake_rate for fru_uptake_rate in fructose_uptake_rate]
    print '\nThe total # of cases to consider = {}\n'.format(len(all_cases))
    print 'Simulating slice {}\n'.format((start_pos,end_pos))

    if start_pos != None and end_pos != None:
        cases_to_consider = all_cases[start_pos - 1:end_pos]
        counter = start_pos - 1
    else:
        cases_to_consider = all_cases
        counter = 0 

    # Initializations
    sucrose_games = []

    # Check the required time for one run of the loop
    if len(cases_to_consider) == 1:
        # Total processing and wall time required to create the pyomo model, solve it and store the results
        start_pt = time.clock()
        start_wt = time.time()

    for (capture_eff,atp_coeff,his_uptake_rate,glc_uptake_rate,fru_uptake_rate) in cases_to_consider:

        counter += 1
        print '{}. (capture_eff = {}, atp_coeff = {}, his_uptake_rate = {}, glc_uptake_rate = {}, fru_uptake_rate = {})'.format(counter,capture_eff,atp_coeff,his_uptake_rate,glc_uptake_rate,fru_uptake_rate)

        #--- Creating a shared memory using the manager ---
        input_data = {}
        input_data['model_path'] = model_path
        input_data['stdout_msgs'] = stdout_msgs
        input_data['warnings'] = warnings
        input_data['iAZ900'] = iAZ900
        input_data['capture_eff'] = capture_eff
        input_data['atp_coeff'] = atp_coeff
        input_data['his_uptake_rate'] = his_uptake_rate
        input_data['glc_uptake_rate'] = glc_uptake_rate
        input_data['fru_uptake_rate'] = fru_uptake_rate
        input_data['o2_uptake_rate'] = o2_uptake_rate
        input_data['simulate_his'] = simulate_his
        input_data['fixed_o2'] = fixed_o2
        input_data['death_rate'] = death_rate
        input_data['results_filename'] = results_filename
 
        output_data = Manager().dict()
        output_data['suc_game'] = None

        if task.lower() == 'analyze_games':
            p = Process(target = gameCreator, args = (input_data,output_data))
        elif task.lower() == 'analyze_cost':
            p = Process(target = computeCost, args = (input_data,))
        p.start()
        p.join()
        if p.exitcode > 0:
            raise RuntimeError('Error in python subprocess. Please check gameCreator\n')
        else:
            sucrose_games.append(output_data['suc_game'])

    return sucrose_games

    if len(cases_to_consider) == 1:
        # Time required to perform FBA
        elapsed_pt = (time.clock() - start_pt)
        elapsed_wt = (time.time() - start_wt)
        print '\nTotal required time for one run of the loop: Processor time = {} sec/{} min  ,  Walltime = {} sec/{} min\n'.format(elapsed_pt,elapsed_pt/60,elapsed_wt,elapsed_wt/60)

def computeCost(input_data):
    """
    Computes the cost of cooperation for a given histidin uptake and a given capture efficiency.  
    The cost is computed as the reduction in biomass flux compared to a wild-type strain  
    (no histidin mutation) when carrying one unit of flux for SUCRe
    """
    model_path = input_data['model_path']
    iAZ900 = input_data['iAZ900']
    capture_eff = input_data['capture_eff']
    atp_coeff = input_data['atp_coeff']
    his_uptake_rate = input_data['his_uptake_rate']
    glc_uptake_rate = input_data['glc_uptake_rate']
    fru_uptake_rate = input_data['fru_uptake_rate']
    o2_uptake_rate = input_data['o2_uptake_rate']
    fixed_o2 = input_data['fixed_o2']
    simulate_his = input_data['simulate_his']
    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings']
    results_filename = input_data['results_filename']

    #--- Incorporate capture efficiency into the model ---
    # In order to incorporate the capture efficiency into the model, change the SUCRe equation from 
    # SUCRe:   1.0 h2o_e + 1.0 sucr_e --> 1.0 fru_e + 1.0 glc_D_e
    # to (for a capture efficiency of 0.01):
    #          1.0 h2o_e + 1.0 sucr_e --> 0.01 fru_e + 0.01 glc_D_e + 0.99 fru_secreted + 0.99 glc_secreted
    # and then add two reactions exclusively for the secreion of glucose and fructose
    # EX_glc_secreted: glc_D_secreted --> and EX_fru_secreted: fru_secreted -->
    SUCRe = iAZ900.get_reactions({'SUCRe':'id'})
    glc_D_e = iAZ900.get_compounds({'glc_D_e':'id'})
    fru_e = iAZ900.get_compounds({'fru_e':'id'})
    if capture_eff != 0:
        SUCRe.stoichiometry[glc_D_e] = capture_eff
        SUCRe.stoichiometry[fru_e] = capture_eff
    else:
        SUCRe.del_compounds([glc_D_e,fru_e])
    glc_D_secreted = compound(id = 'glc_D_secreted',name = 'D-Glucose',compartment = glc_D_e.compartment)
    fru_secreted = compound(id = 'fru_secreted',name = 'D-Fructose',compartment = fru_e.compartment)
    if capture_eff < 1:
        SUCRe.add_compounds({glc_D_secreted:1 - capture_eff,fru_secreted:1 - capture_eff})
    EX_glc_secreted = reaction(id = 'EX_glc_secreted', name = 'D-Glucose exclusive export',stoichiometry = {glc_D_secreted:-1}, type = 'exchange')
    EX_glc_secreted.objective_coefficient = 0
    EX_fru_secreted = reaction(id = 'EX_fru_secreted', name = 'D-Fructose exclusive export',stoichiometry = {fru_secreted:-1}, type = 'exchange')
    EX_fru_secreted.objective_coefficient = 0
    iAZ900.add_reactions([EX_glc_secreted,EX_fru_secreted])
    iAZ900.validate()
    if stdout_msgs:
        print 'SUCRe:\t',SUCRe.get_equation('id')

    # If analyzing histidine cost, incorporate the current atp into SUCRe and find out biomass flux for 
    # an excess amount of histidine
    if simulate_his:
        #--- Incorporate atp into SUCRe ---
        if atp_coeff > 0:
            SUCRe = iAZ900.get_reactions({'SUCRe':'id'})
            atp = iAZ900.get_compounds({'atp_c':'id'})
            adp = iAZ900.get_compounds({'adp_c':'id'})
            pi = iAZ900.get_compounds({'pi_c':'id'})
            SUCRe.add_compounds({atp:-atp_coeff,adp:atp_coeff,pi:atp_coeff})

        # compute the biomass flux for the reference strain by allowing for an excess amount of histidine uptake 
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-100,1000],'IGPDH':[0,0]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-100,1000],'IGPDH':[0,0]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        fba_reference_soln = iAZ900.fba_model.solution

    # For atp cost, do not incoporate IGPDH mutation and compute the biomass flux for the reference strain by 
    # without incorporaing atp into SUCRe
    else: 
        # compute the biomass flux for the reference strain by allowing for an excess amount of histidine uptake 
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        fba_reference_soln = iAZ900.fba_model.solution

    #--- Now compute the biomass flux for the case cooperation cost (histidine uptake rate or atp in SUCRe is incoporated
    if simulate_his:
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        fba_coopr_cost = iAZ900.fba_model.solution

    # If analyzing atp cost
    else: 
        # Incorporate atp in to SUCRe
        if atp_coeff > 0:
            SUCRe = iAZ900.get_reactions({'SUCRe':'id'})
            atp = iAZ900.get_compounds({'atp_c':'id'})
            adp = iAZ900.get_compounds({'adp_c':'id'})
            pi = iAZ900.get_compounds({'pi_c':'id'})
            SUCRe.add_compounds({atp:-atp_coeff,adp:atp_coeff,pi:atp_coeff})

        # Perform the simulations with no histidine mutations
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        fba_coopr_cost = iAZ900.fba_model.solution

    #--- Compute cost ---   
    if fba_reference_soln['exit_flag'] == 'globallyOptimal' and fba_coopr_cost['exit_flag'] == 'globallyOptimal':
        cooperation_cost = fba_reference_soln['objective_value'] - fba_coopr_cost['objective_value']
        if abs(cooperation_cost) <= 1e-8:
            cooperation_cost = 0
    elif fba_reference_soln['exit_flag'] == 'globallyOptimal' and fba_coopr_cost['exit_flag'] != 'globallyOptimal':
        cooperation_cost = 10  # assign a high value
    elif fba_reference_soln['exit_flag'] != 'globallyOptimal' and fba_coopr_cost['exit_flag'] != 'globallyOptimal':
        cooperation_cost = 15   # assign a high but different value (than 10)
    else: # Infeasible WT but feasible iAZ900
        if simulate_his:
            raise userError('FBA is infeasible FBA for wild-type while feasible for cooperator for capture efficiency = ' + str(capture_eff) + ' histidine uptake rate = ' + str(his_uptake_rate))
        else:
            raise userError('FBA is infeasible FBA for wild-type while feasible for cooperator for capture efficiency = ' + str(capture_eff) + ' atp_coeff = ' + str(atp_coeff))

    #--- Store results ---
    if results_filename != '':
        with open(results_filename,'a') as f:
            if simulate_his:
                 f.write("cooperation_cost[(('capture_eff'," + str(capture_eff) + "),('his'," + str(his_uptake_rate) + '))] = ' + str(cooperation_cost) + '\n')
            else: 
                f.write("cooperation_cost[(('capture_eff'," + str(capture_eff) + "),('atp_coeff'," + str(atp_coeff) + '))] = ' + str(cooperation_cost) + '\n')


def plot_gameResults(results_filename, x, y, x_key, y_key, title = '', xaxis_label = '', yaxis_label = '', set_minor_xticks = True, invert_xaxis = True, output_filename = ''):
    """
    Plots a heatmap of the anaylsis results.
 
    INPUTS:
    -------
    results_filenames: A list of strings containing the names of the files containing the results
                     x: A list containing the elements of horizontal axis 
                     x: A list containing the elements of vertical axis 
         x_key & y_key: Name of the x and y variables in results (Note that results (stored in file_name) is with keys as follows:
                        (('capture_eff', 0.2), ('atp_coeff', 0.115), ('his', 0), ('glc', 0), ('fru', 0))
                        For example, if x if histidine uptake rate and y is capture efficiency, then x_key_name is his
                        and y_key_name is capture_eff. 
           xaxis_label: x-axis label
           yaxis_label: y-axis label
      output_fiilename: Name of the output file containing the plot
    """
    import os
    from imp import load_source
    import numpy as np
    from matplotlib import colors
    from tools.ancillary.plot_heatmap import plot_heatmap 
    from tools.ancillary.plot import plot, axis, color_bar 

    if not isinstance(set_minor_xticks,bool):
        raise TypeError('set_minor_xticks must be boolean')

    #--- Load the data ---
    results = {}
    # Import the data in the module stored in file_name
    if type(results_filename) == str:
        if not os.path.isfile(results_filename):
            raise IOError("No such file was found :'" + results_filename + "'")
        else:
            # First delete the model dataFile if it already exists. If it is not deleted
            # the new module data is merged with the previous ones
            try:
                del sys.modules['dataFile']
            except:
                pass
            load_source('dataFile',results_filename)
            import dataFile

    results = dataFile.results

    #--- Types of Nash equilibria ---
    # Find the Nash equilibria as follows:
    # Possible NE are as follows
    # 0)  [[('player_1', 'D'), ('player_2', 'D')]] --> Prisoner's Dilemma 
    # 1)  [[('player_1', 'C'), ('player_2', 'C')]] --> Mutually beneficial
    # 3)  [[('player_1', 'D'), ('player_2', 'D')], [('player_1', 'C'), ('player_2', 'C')]] --> Harmony
    # 3)  [[('player_1', 'D'), ('player_2', 'C')], [('player_1', 'C'), ('player_2', 'D')]] --> Snowdrift
    res = results
    for k in res.keys():
        res[k] = str(res[k])
    print '\n--- Nash equilibria ---'
    Nash_eq = list(set(res.values()))
    for ne in Nash_eq: 
        print ne
    print '\n'

    # 0 = PD --> red  1 = MB --> green  2 = SD --> cyan  3 = HR --> blue  
    # Source: http://stackoverflow.com/questions/30893483/make-a-heatmap-with-a-specified-discrete-color-mapping-with-matplotlib-in-python
    pd = "[[('player_1', 'Defect'), ('player_2', 'Defect')]]" 
    mb = "[[('player_1', 'Cooperate'), ('player_2', 'Cooperate')]]"
    sd = "[[('player_1', 'Defect'), ('player_2', 'Cooperate')], [('player_1', 'Cooperate'), ('player_2', 'Defect')]]"

    # Number of known Nash equilibria (PD, MB and SD)
    if pd not in Nash_eq and mb not in Nash_eq and sd not in Nash_eq:
        known_NE_num = 0
    elif (pd in Nash_eq and mb not in Nash_eq and sd not in Nash_eq) or (pd not in Nash_eq and mb in Nash_eq and sd not in Nash_eq) or (pd not in Nash_eq and mb not in Nash_eq and sd in Nash_eq): 
        known_NE_num = 1
    elif (pd in Nash_eq and mb in Nash_eq and sd not in Nash_eq) or (pd in Nash_eq and mb not in Nash_eq and sd in Nash_eq) or (pd not in Nash_eq and mb in Nash_eq and sd in Nash_eq): 
        known_NE_num = 2
    elif pd in Nash_eq and mb in Nash_eq and sd in Nash_eq: 
        known_NE_num = 3

    if known_NE_num == 0:
        colormap = colors.ListedColormap(['blue'])
        colorbar_ticklabels = ['Mixed']
        data_value_mixed = 0

    elif known_NE_num == 1:
        if pd in Nash_eq: 
             if len(Nash_eq) == 1:
                colormap = colors.ListedColormap(['red'])
                colorbar_ticklabels = ["Prisoner's \nDilemma"]
                data_value_pd = 0
             else:
                colormap = colors.ListedColormap(['red','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma","Mixed"]
                data_value_pd = 0
                data_value_mixed = 1
        elif mb in Nash_eq: 
             if len(Nash_eq) == 1:
                colormap = colors.ListedColormap(['green'])
                colorbar_ticklabels = ["Mutually \nbeneficial"]
                data_value_mb = 0
             else:
                colormap = colors.ListedColormap(['green','blue'])
                colorbar_ticklabels = ["Mutually \nbeneficial","Mixed"]
                data_value_mb = 0
                data_value_mb = 1
        elif sd in Nash_eq: 
             if len(Nash_eq) == 1:
                colormap = colors.ListedColormap(['cyan'])
                colorbar_ticklabels = ["Snowdrift"]
                data_value_sd = 0
             else:
                colormap = colors.ListedColormap(['cyan','blue'])
                colorbar_ticklabels = ["Snowdrift","Mixed"]
                data_value_sd = 0
                data_value_mixed = 1
        else:
            raise userError('Unknown case for known_NE_num = 1')

    elif known_NE_num == 2:
        if pd in Nash_eq and mb in Nash_eq:
             if len(Nash_eq) == 2:
                colormap = colors.ListedColormap(['red', 'green'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial']
                data_value_pd = 0
                data_value_mb = 1
             else:
                colormap = colors.ListedColormap(['red', 'green','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial','Mixed']
                data_value_pd = 0
                data_value_mb = 1
                data_value_mixed = 2
        elif pd in Nash_eq and sd in Nash_eq: 
             if len(Nash_eq) == 2:
                colormap = colors.ListedColormap(['red', 'cyan'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Snowdrift']
                data_value_pd = 0
                data_value_sd = 1
             else:
                colormap = colors.ListedColormap(['red', 'cyan','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Snowdrift','Mixed']
                data_value_pd = 0
                data_value_sd = 1
                data_value_mixed = 2
        elif mb in Nash_eq and sd in Nash_eq: 
             if len(Nash_eq) == 2:
                colormap = colors.ListedColormap(['green', 'cyan'])
                colorbar_ticklabels = ['Mutually \nBeneficial','Snowdrift']
                data_value_mb = 0
                data_value_sd = 1
             else:
                colormap = colors.ListedColormap(['green', 'cyan','blue'])
                colorbar_ticklabels = ['Mutually \nBeneficial','Snowdrift','Mixed']
                data_value_mb = 0
                data_value_sd = 1
                data_value_mixed = 2
        else:
            raise userError('Unknown case for known_NE_num = 2')

    elif known_NE_num == 3:
        if pd in Nash_eq and mb in Nash_eq and sd in Nash_eq:
            if len(Nash_eq) == 3:
                colormap = colors.ListedColormap(['red', 'green','cyan'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial','Snowdrift']
                data_value_pd = 0
                data_value_mb = 1
                data_value_sd = 2
            else:
                colormap = colors.ListedColormap(['red', 'green','cyan','blue'])
                colorbar_ticklabels = ["Prisoner's \nDilemma",'Mutually \nBeneficial','Snowdrift','Mixed']
                data_value_pd = 0
                data_value_mb = 1
                data_value_sd = 2
                data_value_mixed = 3
        else:
            raise userError('Unknown case for known_NE_num = 3')


    data = np.zeros((len(y),len(x)))

    for i, yy in enumerate(y):
        for j, xx in enumerate(x):
            # Find the dictionary key in results 
            res_key = [k for k in results.keys() if dict(k)[x_key] == xx and dict(k)[y_key] == yy]
            if len(res_key) > 1:
                raise userError('More than one element found in results for (({},{}),({},{})) : {}'.format(x_key,xx,y_key,yy),res_key)
            elif len(res_key) == 0:
                raise userError('No key corresponding to (({},{}),({},{})) found in results'.format(x_key,xx,y_key,yy,res_key))
            else:   
                res_key = res_key[0]
                if results[res_key] == pd: 
                    data[i,j] = data_value_pd
                elif results[res_key] == mb: 
                    data[i,j] = data_value_mb
                elif results[res_key] == sd: 
                    data[i,j] = data_value_sd
                else:
                    data[i,j] = data_value_mixed

    #-- Create a color map --

    #plot_heatmap(x = np.array(x),y = 100*np.array(y),data = data,title = title, xaxis_label = xaxis_label, yaxis_label = yaxis_label, plot_func = 'pcolor', colormap = colormap, colorbar_ticklabels = colorbar_ticklabels, set_minor_xticks = set_minor_xticks, set_minor_yticks = True, x_majorticks_spacing = None, x_minorticks_spacing = None,y_majorticks_spacing = None, y_minorticks_spacing = 10,invert_xaxis = invert_xaxis, invert_yaxis = False, grid = False, figsize = (10,6), dpi = None, output_filename = output_filename)

    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, set_minorticks = set_minor_xticks, limits = (min(x),max(x)), invert = invert_xaxis), yaxis = axis(label = yaxis_label, set_minorticks = True, minorticks_spacing = 10, limits = (100*np.array(y).min(),100*np.array(y).max()), invert = False), plot_gridlines = False, fig_format = {'figsize':(7.5,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = np.array(x), y = 100*np.array(y), data = data, plot_func = 'pcolor', clrbar = color_bar(colormap = colormap, colorlimits = (data.min() - 0.5 ,data.max() + 0.5), label = '',  ticklabels =  colorbar_ticklabels))

    print 'The figure was saved into {}\n'.format(output_filename)

def plot_costResults(results_filename, x, y, x_key, y_key, title = '', xaxis_label = '', yaxis_label = '', x_minorticks_spacing = None, set_minor_xticks = True, cost10_replacement_value = 10, cost15_replacement_value = 15, colorbar_label = '', output_filename = ''):
    """
    Plots a heatmap of the anaylsis results.
 
    INPUTS:
    -------
      results_filenames: A list of strings containing the names of the files containing the results
                       x: A list containing the elements of horizontal axis 
                       x: A list containing the elements of vertical axis 
           x_key & y_key: Name of the x and y variables in results (Note that results (stored in file_name) is with keys as follows:
                          (('capture_eff', 0.2), ('atp_coeff', 0.115), ('his', 0), ('glc', 0), ('fru', 0))
                          For example, if x if histidine uptake rate and y is capture efficiency, then x_key_name is his
                          and y_key_name is capture_eff. 
             xaxis_label: x-axis label
             yaxis_label: y-axis label
        set_minor_xticks: Whethter to set minor xtick_labels
    x_minorticks_spacing: Spacing between minor xtick_labels
        output_fiilename: Name of the output file containing the plot
    """
    import os
    from imp import load_source
    import numpy as np
    from matplotlib import colors
    from tools.ancillary.plot_heatmap import plot_heatmap 
    from tools.ancillary.plot import plot, axis, color_bar 

    #--- Load the data ---
    # Import the data in the module stored in results_filename
    if type(results_filename) == str:
        if not os.path.isfile(results_filename):
            raise IOError("No such file was found :'" + results_filename + "'")
        else:
            # First delete the model dataFile if it already exists. If it is not deleted
            # the new module data is merged with the previous ones
            try:
                del sys.modules['dataFile']
            except:
                pass
            load_source('dataFile',results_filename)
            import dataFile

    results = dataFile.cooperation_cost

    # Fiind the maximum value of the cost other than 10 and 15
    max_cost_no10_15 = max([v for v in results.values() if v != 10 and v != 15])
    print '\nMaximum value of the cost other than 10 and 15 = ',max_cost_no10_15

    data = np.zeros((len(y),len(x)))

    for i, yy in enumerate(y):
        for j, xx in enumerate(x):
            # Find the dictionary key in results 
            res_key = [k for k in results.keys() if dict(k)[x_key] == xx and dict(k)[y_key] == yy]
            if len(res_key) > 1:
                raise userError('More than one element found in results for (({},{}),({},{})) : {}'.format(x_key,xx,y_key,yy),res_key)
            elif len(res_key) == 0:
                raise userError("No key corresponding to (('{}',{}),('{}',{})) found in results".format(x_key,xx,y_key,yy,res_key))
            else:
                res_key = res_key[0]
                if results[res_key] == 10:
                    data[i,j] = cost10_replacement_value 
                elif results[res_key] == 15:
                    data[i,j] = cost15_replacement_value 
                else:
                    data[i,j] = results[res_key] 

    print '\nmin results/data = {}/{}   ,  max results/data = {}/{}\n'.format(min(results.values()),data.min(),max([v for v in results.values() if v != 10 and v != 15]),data.max())

    #plot_heatmap(x = np.array(x),y = 100*np.array(y),data = data, title = title, xaxis_label = xaxis_label, yaxis_label = yaxis_label, plot_func = 'pcolor', colormap = None, colorbar_ticklabels = None, set_minor_xticks = set_minor_xticks, set_minor_yticks = True, x_majorticks_spacing = None, x_minorticks_spacing = x_minorticks_spacing,y_majorticks_spacing = None, y_minorticks_spacing = 10, invert_xaxis = False, invert_yaxis = False, grid = False, figsize = (9,6), dpi = None, colorbar_label = colorbar_label, colorbar_label_format = {'distance_from_ticklabels':25}, output_filename = output_filename)

    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, set_minorticks = set_minor_xticks, minorticks_spacing = x_minorticks_spacing, limits = (min(x),max(x)), invert = False), yaxis = axis(label = yaxis_label, set_minorticks = True, minorticks_spacing = 10, limits = (100*np.array(y).min(),100*np.array(y).max()), invert = False), plot_gridlines = False, fig_format = {'figsize':(7,4.5), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = np.array(x), y = 100*np.array(y), data = data, plot_func = 'pcolor', clrbar = color_bar(label = colorbar_label, label_format = {'distance_from_ticklabels':25}))

    print 'The figure was saved into {}\n'.format(output_filename)

def create_intervals(total_comb_num,interval_size):
    """
      total_comb_num: Length of the variable 'combinations'
    interval_size: Desired size of the iteration intervals

    Results are printed in the output
    """
    slices = []

    # Start positions
    start_positions = range(1,total_comb_num + 1,interval_size)
    for k in start_positions:
            if k + (interval_size - 1) <= total_comb_num:
                slices.append((k,k + (interval_size - 1)))
            else:
                slices.append((k,total_comb_num))

    print '\n'
    print '**Total # of intervals = ',len(slices),'\n'

    for slice in slices:
       print slice

    return slices

def create_job_files(total_comb_num,interval_size, job_filename_base, joboutput_filename_base, results_filename_base, task, simulate_his, fixed_o2, capture_efficiency_str, histidine_uptake_rate_str, SUCRe_atp_str, o2_uptake_rate, max_walltime):
    """
    Creates the job files given:
             total_comb_num: Total number of combinations 
              interval_size: The desired of iteration intervals
          outfile_base_name: The base name of the output files
                       task: The task that should be performed ('analyze_games' or 'analyze_cost')
          job_filename_base: Base name of the output job file
    joboutput_filename_base: The name of the file storing the dynamic output of the job. It is essentially the same as job_filename_base
                             with the file path deleted (in order to create .out files in the same directory as the job file)
     capture_efficiency_str: String containing the capture efficiency
    histidine_uptake_rate_str: String containing the hsitidine uptake rate (for simulate_his = True)
              SUCRe_atp_str: String containing SUCRe_atp (for simuulate_his = False) 
      results_filename_base: Base name of the results file for the main python script
    """
    # Task 
    if task not in ['analyze_cost','analyze_games']:
        raise ValueError('Invalid task value')

    # A parameter showing whether we are intereated in simulating the cooperation cost with histidine (True) or with atp (False)
    if not isinstance(simulate_his,bool):
        raise TypeError('Invalid simulate_his type')

    # A parameter showing whether the oxygen uptake rate is fixed (True) or not (False)
    if not isinstance(fixed_o2,bool):
        raise TypeError('Invalid fixed_o2 type') 

    slices = create_intervals(total_comb_num = total_comb_num,interval_size = interval_size)

    for slice in slices:
        job_filename = job_filename_base + '_' + str(slice[0]) + '_' + str(slice[1]) + '.sh'
        joboutput_filename = joboutput_filename_base + '_' + str(slice[0]) + '_' + str(slice[1]) + '.out'
        print 'creating file ',job_filename,'...'
        with open(job_filename,'w') as outfile:
            outfile.write('#!/bin/bash -l\n')
            outfile.write('#$-l h_rt=' + str(int(max_walltime))+ ':00:00\n\n')
            outfile.write('cd /usr2/postdoc/alizom/work/yeast_sucrose\n\n')

            # This is to merge the .e[jobid] and .o[jobid] files
            outfile.write('#$ -j y\n')

            # Set the output file
            outfile.write('\n#$ -o ' + joboutput_filename + '\n')

            outfile.write('\nsource activate /projectnb/bioinfor/SEGRE/alizom\n\n')

            # python -c "import time;print '\n**Job started at ',time.strftime('%c'),'\n'" > job_dynamic_yeast_stsrt_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\"\n")

            start_pos = slice[0]
            end_pos = slice[1]
            results_filename = results_filename_base + "_" + str(start_pos) + "_" + str(end_pos) + '.py'

            if simulate_his:
                outfile.write('\npython -c "from __future__ import division;from simulate_games import master_func;master_func(start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", capture_efficiency = " + capture_efficiency_str + ", histidine_uptake_rate = " + histidine_uptake_rate_str + ", SUCRe_atp = 0.115, simulate_his = " + str(simulate_his) + ", fixed_o2 = " + str(fixed_o2) + ", o2_uptake_rate = " + str(o2_uptake_rate) + ", task = '" + task + "', stdout_msgs = False, warnings = False, results_filename = '" + results_filename + "')\"\n")

            else:
                outfile.write('\npython -c "from __future__ import division;from simulate_games import master_func;master_func(start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ", capture_efficiency = " + capture_efficiency_str  + ", SUCRe_atp = " +  SUCRe_atp_str + ", simulate_his = " + str(simulate_his) + ",fixed_o2 = " + str(fixed_o2) + ", o2_uptake_rate = " + str(o2_uptake_rate) + ", task = '" + task + "', stdout_msgs = False, warnings = False, results_filename = '" + results_filename + "')\"\n")

            # python -c "import time;print '\n**Job ended at ',time.strftime('%c'),'\n'" >> job_dynamic_yeast_start_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n")

def run_create_job_files(create_cost_his_jobs = False, create_cost_atp_jobs = False, create_games_his_jobs = False, create_games_atp_jobs = False, create_all = False):
    """
    It just runs create_job_files to create pbs files for different cases
    """
    if create_all:
        create_cost_his_jobs = True
        create_cost_atp_jobs = True
        create_games_his_jobs = True
        create_games_atp_jobs = True

    o2_uptake_rate = 1000
    capture_efficiency = [i/50 for i in range(51)]
    capture_efficiency_str = '[i/50 for i in range(51)]'
    histidine_uptake_rate = [i/200 for i in range(21)] 
    histidine_uptake_rate_str = '[i/200 for i in range(21)]'
    SUCRe_atp = [i/10 for i in range(0,301,3)] 
    SUCRe_atp_str = '[i/10 for i in range(0,301,3)]' 


    # Total combinations for simulating histidine
    total_comb_his = len(capture_efficiency)*len(histidine_uptake_rate) 
    total_comb_atp = len(capture_efficiency)*len(SUCRe_atp) 
    his_interval_size = 600
    atp_interval_size = 500
    print '\ntotal_comb_his = {} , his_interval_size = {}     total_comb_atp = {} , atp_interval_size = {}\n'.format(total_comb_his,his_interval_size,total_comb_atp,atp_interval_size)

    #--- Create the jobs for cost analysis ---
    if create_cost_his_jobs:
        # Changes in the file names apply only to fixed O2 case
        if o2_uptake_rate == 1000:
            job_filename_base = 'jobs/job_cost_his_indepO2'
            joboutput_filename_base = 'job_cost_his_indepO2'
            results_filename_base = 'results/cost_results_his_indepO2'
        else:
            job_filename_base = 'jobs/job_cost_his_fixedo2'
            joboutput_filename_base = 'job_cost_his_fixedo2'
            results_filename_base = 'results/cost_results_his_fixedo2'

        # his fixed o2
        print '\n--- Creating job file for cost, his, fixedo2 ---'
        create_job_files(total_comb_num = total_comb_his, interval_size = his_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, task = 'analyze_cost',fixed_o2 = True, simulate_his = True, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = o2_uptake_rate, max_walltime = 48)

        # his variable o2
        print '\n--- Creating job file for cost, his, variable O2 ---'
        create_job_files(total_comb_num = total_comb_his, interval_size = his_interval_size, job_filename_base = 'jobs/job_cost_his_variableO2', joboutput_filename_base = 'job_cost_his_variableO2', results_filename_base = 'results/cost_results_his_variableO2', task = 'analyze_cost',fixed_o2 = False, simulate_his = True, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = 2, max_walltime = 48)       
    
    # atp fixed o2
    if create_cost_atp_jobs:
        # Changes in the file names apply only to fixed O2 case
        if o2_uptake_rate == 1000:
            job_filename_base = 'jobs/job_cost_atp_indepO2' 
            joboutput_filename_base = 'job_cost_atp_indepO2' 
            results_filename_base = 'results/cost_results_atp_indepO2'
        else:
            job_filename_base = 'jobs/job_cost_atp_fixedo2' 
            joboutput_filename_base = 'job_cost_atp_fixedo2' 
            results_filename_base = 'results/cost_results_atp_fixedo2'

        print '\n--- Creating job file for cost, atp, fixedo2 ---'
        create_job_files(total_comb_num = total_comb_atp, interval_size = atp_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, task = 'analyze_cost',fixed_o2 = True, simulate_his = False, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = o2_uptake_rate, max_walltime = 48)

        # his variable o2
        print '\n--- Creating job file for cost, atp, variable O2 ---'
        create_job_files(total_comb_num = total_comb_atp, interval_size = atp_interval_size, job_filename_base = 'jobs/job_cost_atp_variableO2', joboutput_filename_base = 'job_cost_atp_variableO2', results_filename_base = 'results/cost_results_atp_variableO2', task = 'analyze_cost',fixed_o2 = False, simulate_his = False, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = 2, max_walltime = 48)

    #--- Create the jobs for game analysis ---
    if create_games_his_jobs:
        if o2_uptake_rate == 1000:
            job_filename_base = 'jobs/job_ss_games_his_indepO2' 
            joboutput_filename_base = 'job_ss_games_his_indepO2' 
            results_filename_base = 'results/ss_game_results_his_indepO2'
        else:
            job_filename_base = 'jobs/job_ss_games_his_fixedo2' 
            joboutput_filename_base = 'job_ss_games_his_fixedo2' 
            results_filename_base = 'results/ss_game_results_his_fixedo2'

        # his fixed o2
        print '\n--- Creating job file for games, his, fixedo2 ---'
        create_job_files(total_comb_num = total_comb_his, interval_size = his_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, task = 'analyze_games',fixed_o2 = True, simulate_his = True, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = o2_uptake_rate, max_walltime = 48)

        # his variable o2
        print '\n--- Creating job file for games, his, variable O2 ---'
        create_job_files(total_comb_num = total_comb_his, interval_size = his_interval_size, job_filename_base = 'jobs/job_ss_games_his_variableO2', joboutput_filename_base = 'job_ss_games_his_variableO2', results_filename_base = 'results/ss_game_results_his_variableO2', task = 'analyze_games',fixed_o2 = False, simulate_his = True, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = 2, max_walltime = 48)       
    
    if create_games_atp_jobs:
        if o2_uptake_rate == 1000:
            job_filename_base = 'jobs/job_ss_games_atp_indepO2'
            joboutput_filename_base = 'job_ss_games_atp_indepO2'
            results_filename_base = 'results/ss_game_results_atp_indepO2'
        else:
            job_filename_base = 'jobs/job_ss_games_atp_fixedo2'
            joboutput_filename_base = 'job_ss_games_atp_fixedo2'
            results_filename_base = 'results/ss_game_results_atp_fixedo2'

        # atp fixed o2
        print '\n--- Creating job file for games, atp, fixedo2 ---'
        create_job_files(total_comb_num = total_comb_atp, interval_size = atp_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, task = 'analyze_games',fixed_o2 = True, simulate_his = False, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = o2_uptake_rate, max_walltime = 48)

        # atp variable o2
        print '\n--- Creating job file for games, atp, variable O2 ---'
        create_job_files(total_comb_num = total_comb_atp, interval_size = atp_interval_size, job_filename_base = 'jobs/job_ss_games_atp_variableO2', joboutput_filename_base = 'job_ss_games_atp_variableO2', results_filename_base = 'results/ss_game_results_atp_variableO2', task = 'analyze_games',fixed_o2 = False, simulate_his = False, capture_efficiency_str = capture_efficiency_str, histidine_uptake_rate_str = histidine_uptake_rate_str, SUCRe_atp_str = SUCRe_atp_str, o2_uptake_rate = 2, max_walltime = 48)

def run_plot_results():
    """
    Plots the results
    """
    capture_efficiency = [i/50 for i in range(51)]
    histidine_uptake_rate = [i/200 for i in range(21)] 
    SUCRe_atp = [i/10 for i in range(0,301,3)] 

    cost_his_fixedo2 = True
    cost_his_indepO2 = True
    cost_his_variableO2 = True

    game_his_fixedo2 = True
    game_his_indepO2 = True
    game_his_variableO2 = True

    cost_atp_fixedo2 = True
    cost_atp_indepO2 = True
    cost_atp_variableO2 = True

    game_atp_fixedo2 = True
    game_atp_indepO2 = True
    game_atp_variableO2 = True

    # ----- Plot the cost results -----
    # Use \rm to remove italic from matplotlib LaTex text.
    # Source: http://stackoverflow.com/questions/19671659/remove-italics-in-latex-subscript-in-matplotlib
    # You can also use \mathregular{} to make the font regular
    #-- his --
    if cost_his_fixedo2: 
        plot_costResults(results_filename = 'results/cost_results_his_fixedo2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Fixed limiting O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', x_minorticks_spacing = None, cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'Histidine production cost', output_filename = 'results/cost_results_his_fixedo2.pdf')

    if cost_his_indepO2: 
        plot_costResults(results_filename = 'results/cost_results_his_indepO2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Excess O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', x_minorticks_spacing = None, cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'Histidine production cost', output_filename = 'results/cost_results_his_indepO2.pdf')

    if cost_his_variableO2: 
        plot_costResults(results_filename = 'results/cost_results_his_variableO2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Variable O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'Histidine production cost', output_filename = 'results/cost_results_his_variableO2.pdf')

    #-- ATP --
    if cost_atp_fixedo2: 
        plot_costResults(results_filename = 'results/cost_results_atp_fixedo2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Fixed limiting O$_{\rm 2}$', xaxis_label = 'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'SUCRe ATP Cost', output_filename = 'results/cost_results_atp_fixedo2.pdf')

    if cost_atp_indepO2: 
        plot_costResults(results_filename = 'results/cost_results_atp_indepO2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Excess O$_{\rm 2}$', xaxis_label = 'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, cost10_replacement_value = 1.6, cost15_replacement_value = 1.65, colorbar_label = 'SUCRe ATP Cost', output_filename = 'results/cost_results_atp_indepO2.pdf')

    if cost_atp_variableO2: 
        plot_costResults(results_filename = 'results/cost_results_atp_variableO2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Variable O$_{\rm 2}$', xaxis_label = 'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'SUCRe ATP Cost', output_filename = 'results/cost_results_atp_variableO2.pdf')

    # ------ Plot the game results ------
    # -- his --
    if game_his_fixedo2:
        plot_gameResults(results_filename = 'results/ss_game_results_his_fixedo2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Fixed limitting O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = True, output_filename = 'results/ss_game_results_his_fixedo2.pdf')
        #plot_gameResults(results_filename = 'results/ss_game_results_his_fixedo2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'', xaxis_label = r'Histidine uptake rate ($\rm{\frac{mmol}{gDW.h}}$)', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = True, output_filename = 'results/ss_game_results_his_fixedo2.pdf')

    if game_his_indepO2:
        plot_gameResults(results_filename = 'results/ss_game_results_his_indepO2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Excess O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = True, output_filename = 'results/ss_game_results_his_indepO2.pdf')

    if game_his_variableO2:
        plot_gameResults(results_filename = 'results/ss_game_results_his_variableO2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Variable O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = True, output_filename = 'results/ss_game_results_his_variableO2.pdf')

    #-- ATP --
    if game_atp_fixedo2:
        plot_gameResults(results_filename = 'results/ss_game_results_atp_fixedo2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Fixed limitting O$_{\rm 2}$', xaxis_label = r'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = False, output_filename = 'results/ss_game_results_atp_fixedo2.pdf')
        #plot_gameResults(results_filename = 'results/ss_game_results_atp_fixedo2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'', xaxis_label = r'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = False, output_filename = 'results/ss_game_results_atp_fixedo2.pdf')

    if game_atp_indepO2:
        plot_gameResults(results_filename = 'results/ss_game_results_atp_indepO2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Excess O$_{\rm 2}$', xaxis_label = r'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = False, output_filename = 'results/ss_game_results_atp_indepO2.pdf')

    if game_atp_variableO2:
        plot_gameResults(results_filename = 'results/ss_game_results_atp_variableO2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Variable O$_{\rm 2}$', xaxis_label = r'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = False, output_filename = 'results/ss_game_results_atp_variableO2.pdf')

#-----------------------------------
if __name__ == '__main__':
    print '\n------- Test -------- '
    # Create the payoff matrix of the game
    #sucrose_game = master_func(capture_efficiency = 0.01, histidine_uptake_rate = 1, stdout_msgs = True)
    #sucrose_game = master_func(capture_efficiency = 0.9, histidine_uptake_rate = [0.003,0.004,0.006,0.008,0.01], stdout_msgs = True, results_file = 'a.py')
    sucrose_game = master_func(capture_efficiency = 0.2, histidine_uptake_rate = 0, glucose_uptake_rate = 0, fructose_uptake_rate = 0, o2_uptake_rate = 2, simulate_his = True, fixed_o2 = True, task = 'analyze_games', warnings = True, stdout_msgs = True, results_filename = 't2.py')
    sucrose_game = sucrose_game[0] 
    print '-- The payoff matrix is as follows --'
    for s in sorted(sucrose_game.payoff_matrix.keys()):
        print '{}: {}'.format(s,sucrose_game.payoff_matrix[s])    

    print '\n-- Nash Equilibria of the game --'
    print 'exit flag = ',sucrose_game.pureNashEq_exitflag
    for eq in sucrose_game.pureNash_equilibria:
        print eq
   
