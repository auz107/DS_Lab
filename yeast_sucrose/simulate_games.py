from __future__ import division
import sys, time
sys.path.append('../')
from tools.io.read_sbml_model import read_sbml_model
from tools.io.create_model import create_model
from tools.core.organism import organism
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.compartment import compartment
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.GAMETES.game import game
from tools.GAMETES.replicator_dynamics import replicator_dynamics
from tools.userError import userError
from coopr.pyomo import *
from coopr.opt import *
from copy import deepcopy
from multiprocessing import Process, Manager, Value, Array
import numpy as np
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

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
    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/'
    # Define the organism (see function findHisSatConc for the calculations of gDW_per_cell for yeast)
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '',gDW_per_cell = 4e-11)

    # Create the metabolic model object
    model = create_model(model_organism = model_organism, model_info = {'id':'iAZ900', 'file_format':'sbml', 'model_filename':model_path + 'iAZ900_noCycles_03_25_2015.xml', 'biomassrxn_id':'biomass_core'}, growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iAZ900_minimal.py', 'flux_bounds_dict': {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]}}, perform_fba = True)

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
    sucrose_uptake_rate = input_data['sucrose_uptake_rate']
    his_uptake_rate = input_data['his_uptake_rate']
    glc_uptake_rate = input_data['glc_uptake_rate']
    fru_uptake_rate = input_data['fru_uptake_rate']
    o2_uptake_rate = input_data['o2_uptake_rate']
    simulate_his = input_data['simulate_his']
    fixed_o2 = input_data['fixed_o2']
    death_rate = input_data['death_rate']
    relax_ATPM = input_data['relax_ATPM']
    simulate_rep_dynamics = input_data['simulate_rep_dynamics']
    t0 = input_data['t0']
    tf = input_data['tf']
    dt = input_data['dt']
    save_details = input_data['save_details']
    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings']
    results_filename = input_data['results_filename']

    if stdout_msgs:
        if stdout_msgs:
            print '\n--- FBA for iAZ900 before adding atp/adp to SUCRe---'

        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)

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
        SUCRe.set_stoichiometry(stoichiometry = dict(SUCRe.stoichiometry.items() + [(atp,-atp_coeff),(adp,atp_coeff),(pi,atp_coeff)]), replace = True, model = iAZ900)
        iAZ900.set_cpds_genes_rxns(do_cpds = True, do_gens = False)

        if stdout_msgs:
            print '\n--- FBA for iAZ900 after adding atp/adp to SUCRe---'
            print 'SUCRe:\t',SUCRe.get_equation('id')

            if iAZ900.fixed_o2:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
            iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = True)

            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])

    #--- Find the metabolic cost of sucrose production and secretion ---
    # invertase_cost = -[(max biomass before the addition of ATP to SUCRe) - (max biomass after the additon of ATP to SUCRe)]
    # Compute the invertase production cost for the case where sucrose is the only carbon source  
    # (do not allow for the uptake of extra glucose and fructose from the medium)
    if iAZ900.fixed_o2:
        set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
    else:
        set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000], 'EX_o2_e_':[-10/5,1000]})
    iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = True)

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
    if atp_coeff > 0 and invertase_cost == 0:
        raise userError('Positive invertase cost is zero') 
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
        SUCRe.set_stoichiometry({glc_D_e:capture_eff, fru_e:capture_eff}, replace = False, model = iAZ900)
    else:
        SUCRe.del_compounds([glc_D_e,fru_e], model = iAZ900)
    glc_D_secreted = compound(id = 'glc_D_secreted',name = 'D-Glucose',compartment = glc_D_e.compartment)
    fru_secreted = compound(id = 'fru_secreted',name = 'D-Fructose',compartment = fru_e.compartment)
    if capture_eff < 1:
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
        if capture_eff < 1:
            print 'EX_glc_secreted:\t',EX_glc_secreted.get_equation('id')
            print 'EX_fru_secreted:\t',EX_fru_secreted.get_equation('id')

        if simulate_his:
            if iAZ900.fixed_o2:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            if iAZ900.fixed_o2:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})

        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = True)
        print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  , EX_glc_secreted = {} ,  EX_fru_e_ = {}  ,  EX_fru_secreted = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted'])

    #--- Cheater and cooperator strains ---
    # Cooperator has a defective HIS3 (YOR202W) gene (for histidone production) corresponding to
    # reaction IGPDH in the model
    # Cheater is missing the gene SUC2 (YIL162W) corresponding ot reaction SUCRe

    #--- Create the players' strategies ---
    players_names = ['player_1','player_2']
    players_strategies = {}
    players_strategies['player_1'] = ['Cooperator','Defector']
    players_strategies['player_2'] = ['Cooperator','Defector']

    payoff_matrix = {}    

    # --- Find out glc and fru secretion rates ---
    if stdout_msgs:
        print '\nFBA for Cooperator alone '

    # NOTE: In the following we se the LB for SUCRe flux to sucrose_uptake_rate. This is because
    #       if glc_uptake_rate or fru_uptake_rate are not zero, then FBA picks to not take up
    #       any sucrose at all
    if glc_uptake_rate > 0 or fru_uptake_rate > 0:
        SUCRe_lb = sucrose_uptake_rate
    else:
        SUCRe_lb = 0 

    if simulate_his:
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[SUCRe_lb, 1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[sucrose-uptake_rate, 1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
    else:
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[SUCRe_lb, 1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[SUCRe_lb, 1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})

    iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = True)
    if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
        # Cooperator's payoff
        glc_secretion_flux = iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted']
        fru_secretion_flux = iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted']
        coopr_SUCRe_flux = iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe']
        if stdout_msgs:
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  glc_secretion_flux = {} , fru_secretion_flux = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],glc_secretion_flux,fru_secretion_flux,iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])
    else:
        coopr_SUCRe_flux = sucrose_uptake_rate 
        glc_secretion_flux = (1 - capture_eff)*sucrose_uptake_rate 
        fru_secretion_flux = (1 - capture_eff)*sucrose_uptake_rate 

    #-- Cooperator vs. Cooperator --
    if stdout_msgs:
        print '-- Cooperator vs. Cooperator --'

    if simulate_his:
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000], 'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
    else:
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000], 'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_secretion_flux-glc_uptake_rate,1000],'EX_fru_e_':[-fru_secretion_flux-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})

    iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
    if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
        # Cooperator's payoff
        if simulate_his and his_uptake_rate == 0:
            payoff_C = death_rate + invertase_cost
        else:
            payoff_C = iAZ900.fba_model.solution['objective_value'] 

        if stdout_msgs:
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  glc_secretion_flux = {} , fru_secretion_flux = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],glc_secretion_flux,fru_secretion_flux,iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])

    else:
        if warnings:
            print 'WARNING!FBA for CC was not solved to optimality. Zero was assigned as payoffs'
        payoff_C = death_rate + invertase_cost 

    payoff_matrix[(('player_1','Cooperator'),('player_2','Cooperator'))] = {'player_1':payoff_C,'player_2':payoff_C}

    #-- Cooperator vs. Cheater --
    payoff_C, payoff_D = None, None
    if stdout_msgs:
        print '-- Cooperator vs. Cheater --'
        print '\nglc_secretion_flux = {} , fru_secretion_flux = {}\n'.format(glc_secretion_flux,fru_secretion_flux)
        print 'FBA for Cooperator:'

    # Cooperator
    if simulate_his and his_uptake_rate == 0:
        payoff_C = death_rate + invertase_cost
    else:
        if simulate_his:
            if iAZ900.fixed_o2:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000], 'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            if iAZ900.fixed_o2:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000], 'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'SUCRe':[coopr_SUCRe_flux,1000], 'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
    
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
        if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
            payoff_C = iAZ900.fba_model.solution['objective_value'] 
    
        else:
            payoff_C = death_rate + invertase_cost 
            if warnings:
                print 'WARNING!FBA for CC was not solved to optimality. Zero was assigned as payoffs'
    

    # Cheater 
    if stdout_msgs:
        print 'FBA for Cheater:'

    if iAZ900.fixed_o2:
        set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[- glc_uptake_rate - glc_secretion_flux,1000],'EX_fru_e_':[- fru_uptake_rate - fru_secretion_flux,1000],'EX_o2_e_':[-o2_uptake_rate,1000],'SUCRe':[0,0]})
    else:
        set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_uptake_rate - glc_secretion_flux,1000],'EX_fru_e_':[- fru_uptake_rate - fru_secretion_flux,1000],'EX_o2_e_':[(- glc_uptake_rate - fru_uptake_rate - glc_secretion_flux - fru_secretion_flux)/5,1000],'SUCRe':[0,0]})

    iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
    if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
        payoff_D = iAZ900.fba_model.solution['objective_value']
        if stdout_msgs:
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,   EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])
    else:
        if stdout_msgs:
            print 'WARNING! FBA for Cheater in CD was not solved to optimality. Death rate ({}) was assigned as payoff'.format(death_rate)
        payoff_D = death_rate

    payoff_matrix[(('player_1','Cooperator'),('player_2','Defector'))] = {'player_1':payoff_C,'player_2':payoff_D}
    payoff_matrix[(('player_1','Defector'),('player_2','Cooperator'))] = {'player_1':payoff_D,'player_2':payoff_C}

    #-- Cheater vs.Cheater --
    payoff_D = None
    if stdout_msgs:
        print '-- Cheater vs. Cheater --'
    # Defector vs. Defector not only cannot they grow but also cannot satisfy their ATPM maintenance. 
    # So their payoffs is equal to their death rate.
    if glc_uptake_rate == 0 and fru_uptake_rate == 0:
        payoff_D = death_rate
    else:
        if stdout_msgs:
            print 'FBA for Cheater:'
    
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[- glc_uptake_rate,1000],'EX_fru_e_':[- fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000],'SUCRe':[0,0]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[- fru_uptake_rate,1000],'EX_o2_e_':[(- glc_uptake_rate - fru_uptake_rate)/5,1000],'SUCRe':[0,0]})
    
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
        if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
            payoff_D = iAZ900.fba_model.solution['objective_value']
            if stdout_msgs:
                print '\nEX_sucr_e_ = {}, SUCRe = {}  ,   EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])
        else:
            if stdout_msgs:
                print 'WARNING! FBA for Cheater in CD was not solved to optimality. Death rate ({}) was assigned as payoff'.format(death_rate)
            payoff_D = death_rate

    payoff_matrix[(('player_1','Defector'),('player_2','Defector'))] = {'player_1':payoff_D,'player_2':payoff_D}

    # Find Nash equilibria
    suc_game = game(game_name = 'yeast_sucrose' + str(capture_eff), players_names = players_names, players_strategies = players_strategies, payoff_matrix = payoff_matrix)
    suc_game.find_NashEq(stdout_msgs = stdout_msgs)

    # Find the payoff matrix of the symmetric game (to save in file)
    suc_game.create_symmetric_payoff_matrix()

    output_data['suc_game'] = suc_game

    game_info = {'NashEq':list(set([tuple(sorted(neq)) for neq in [dict(e).values() for e in suc_game.pureNash_equilibria]])), 'payoff matrix':suc_game.symmetric_payoff_matrix}

    results_key = (('capture_eff',capture_eff),('atp_coeff',atp_coeff),('suc', sucrose_uptake_rate), ('his',his_uptake_rate),('glc',glc_uptake_rate),('fru',fru_uptake_rate), ('o2', o2_uptake_rate))

    #--- Simulate replicator's dynamics ---
    if simulate_rep_dynamics:
        input_data_repdyn = {}
        input_data_repdyn['results_key'] = results_key 
        input_data_repdyn['game_info'] = game_info
        input_data_repdyn['t0'] = t0
        input_data_repdyn['tf'] = tf
        input_data_repdyn['dt'] = dt
        input_data_repdyn['save_details'] = save_details
        input_data_repdyn['results_filename'] = results_filename
        input_data_repdyn['stdout_msgs'] = stdout_msgs

        output_data_repdyn = Manager().dict()
        output_data_repdyn['replicator_dynamics_info'] = None

        p = Process(target = do_replicator_dynamics, args = (input_data_repdyn, output_data_repdyn))
        p.start()
        p.join()
        if p.exitcode > 0:
            raise userError('Error in python subprocess. Please check do_replicator_dynamics\n')
        else:
            replicator_dynamics_info = output_data_repdyn['replicator_dynamics_info']

    else:
        if stdout_msgs:
            print '\nNash equilibria:'
            for neq in game_info['NashEq']:
                print neq

        if results_filename != '':
            if not save_details:
                del game_info['payoff matrix']

            replicator_dynamics_info = {}
            results = {'game_info':game_info, 'replicator_dynamics_info':replicator_dynamics_info}

            with open(results_filename,'a') as f:
                f.write('results[' + str(results_key) + '] = ' + str(results) + '\n')


def do_replicator_dynamics(input_data_repdyn, output_data_repdyn):
    """
    Performs replicator's dynamic simulations
    """
    results_key = input_data_repdyn['results_key'] 
    game_info = input_data_repdyn['game_info']
    t0 = input_data_repdyn['t0']
    tf = input_data_repdyn['tf']
    dt = input_data_repdyn['dt']
    save_details = input_data_repdyn['save_details']
    results_filename = input_data_repdyn['results_filename']
    stdout_msgs = input_data_repdyn['stdout_msgs']

    replicator_dynamics_info = {}

    # Check if any of the Cheater can invade the Cooperator 
    x_init = {'Cooperator':0.99,'Defector':0.01} 
    strains_fracs = replicator_dynamics(payoff_matrix = dict([(k,v) for (k, v) in game_info['payoff matrix'].iteritems()]), x_init = x_init, time_range = [t0,dt,tf])
    # Strain fractions in the last time point
    replicator_dynamics_info[tuple(x_init.items())] = dict([(strain_name, strains_fracs[strain_name][tf]) for strain_name in x_init.keys()])

    # Check if wild-type can invade the population of mutants 
    x_init = {'Cooperator':0.01,'Defector':0.99} 
    strains_fracs = replicator_dynamics(payoff_matrix = dict([(k,v) for (k, v) in game_info['payoff matrix'].iteritems()]), x_init = x_init, time_range = [t0,dt,tf])
    # Strain fractions in the last time point
    replicator_dynamics_info[tuple(x_init.items())] = dict([(strain_name, strains_fracs[strain_name][tf]) for strain_name in x_init.keys()])

    output_data_repdyn['replicator_dynamics_info'] = replicator_dynamics_info

    if stdout_msgs:
        print "\n-- Strain fractions according to Replicator's dynamics ---\n"
        for x0 in replicator_dynamics_info.keys():
            print '\nx_init = {} --> x_tf = {}\n'.format(x0,replicator_dynamics_info[x0])

        print '\nNash equilibria:'
        for neq in game_info['NashEq']:
            print neq

    if results_filename != '':
        if not save_details:
            del game_info['payoff matrix']
        results = {'game_info':game_info, 'replicator_dynamics_info':replicator_dynamics_info}
        with open(results_filename,'a') as f:
            f.write('results[' + str(results_key) + '] = ' + str(results) + '\n')


def master_func(start_pos = None, end_pos = None, capture_efficiency = 0.99, SUCRe_atp = 0.115, sucrose_uptake_rate = 10, histidine_uptake_rate = 10, glucose_uptake_rate = 0, fructose_uptake_rate = 0, o2_uptake_rate = 2, simulate_his = True, fixed_o2 = True, relax_ATPM = True, task = 'analyze_games', simulate_rep_dynamics = True, t0 = 0, dt = 1, tf = 5000, save_details = False, results_filename_base = '', results_filename = '', warnings = True, stdout_msgs = True,):

    """
    Creates the required inputs for gameCreator

    INPUTS:
    -------
    start_pos & end_pos: 
    Start and end positions of the array containing all possble cases to consider (see all_cases variable)
    NOTE: The user can enter the start and end positions numbers assuming that they start at one.
          The code will take care of array indexing convention of python (indexing starts at zero) 


    simulate_his: 
    A parameter showing whether we are intereated in simulating the cooperation cost with histidine (True) 
    or with atp (False)

    fixed_o2: 
    A parameter showing whether the oxygen uptake rate must be fixed at -2 for all simulations (True) or not (False)    

    relax_ATPM:
    If True, the constraint on ATPM is relaxed when the fba problem for a cooperator alone
    is infeasible and the fba problem is resolved to check whether the sucrose hydrolysis
    reaction can carry any flux an dproduce glc and fru. . 

    task: 
    A string showing what the function needs to perform. Eligible cases are analyze_games (finding the NE of games)
    and analyze_cost (assessing the cooperation cost)

    simulate_rep_dynamics:
    If True, simulates the replicator's dynamics
   
    t0, tf, dt:
    Start and final simulation times and delta for discretizing the ODE 

    results_filename_base: 
    The base for the file name storing the results
    Example: results_filename_base = 'results/emc_results'. 
    The code will add the start and end positions to the file name.
    Example: 'results/emc_results_1_500.txt'

    results_filename:
    The same as results_filename_base except that if this parameter is provided
    the start and end positions are not added to the file name. Provide only 
    results_filename_base or results_filename 
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

    if not isinstance(sucrose_uptake_rate,int) and not isinstance(sucrose_uptake_rate, float):
        raise TypeError('sucrose_uptake_rate must be a float or an integer')

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

    if not isinstance(relax_ATPM,bool):
        raise TypeError('relax_ATPM must be either True to False')

    if not isinstance(simulate_rep_dynamics,bool):
        raise TypeError('simulate_rep_dynamics must be either True to False')

    if not isinstance(save_details,bool):
        raise TypeError('save_details must be either True to False')

    if not isinstance(t0,int) and not isinstance(t0, float):
        raise TypeError('t0 must be either an integer or a float')
    if not isinstance(tf,int) and not isinstance(tf, float):
        raise TypeError('tf must be either an integer or a float')
    if not isinstance(dt,int) and not isinstance(dt, float):
        raise TypeError('dt must be either an integer or a float')

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
    # Define the organism (see function findHisSatConc for the calculations of gDW_per_cell for yeast)
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '',gDW_per_cell = 4e-11)

    # Create the metabolic model object
    iAZ900 = create_model(model_organism = model_organism, model_info = {'id':'iAZ900', 'file_format':'sbml', 'model_filename':model_path + 'iAZ900_noCycles_03_25_2015.xml', 'biomassrxn_id':'biomass_core'}, growthMedium_flux_bounds = {'flux_bounds_filename':model_path + 'iAZ900_minimal.py', 'flux_bounds_dict': {'EX_sucr_e_':[-sucrose_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]}}, perform_fba = True)

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
    glc_D_e.ms_calc(model = iAZ900, stdout_msgs = False, warnings = False)
    glc_D_e.biomass_yield_calc(model = iAZ900, stdout_msgs = False, warnings = False)
    if stdout_msgs:
        print 'glc ms = {} , glc biomass yield = {}'.format(glc_D_e.ms, glc_D_e.biomass_yield)
    mu_death_glc_D = - glc_D_e.biomass_yield*glc_D_e.ms 
    if mu_death_glc_D == 0:
        print '**WARNING! mu_death_glc_D is zero!'
    fru_e = iAZ900.get_compounds({'fru_e':'id'})
    fru_e.ms_calc(model = iAZ900, stdout_msgs = False, warnings = False)
    fru_e.biomass_yield_calc(model = iAZ900, stdout_msgs = False, warnings = False)
    if stdout_msgs:
        print 'fru ms = {} , fru biomass yield = {}'.format(fru_e.ms, fru_e.biomass_yield)
    mu_death_fru = - fru_e.biomass_yield*fru_e.ms 
    if mu_death_fru == 0:
        print '**WARNING! mu_death_fru is zero!'

    # -0.001 is included in the following to make sure of having a negative death rate
    # in case the glc and fru death rates turn out to be zero for any reasosns
    death_rate = min(mu_death_glc_D,mu_death_fru,-0.001)
    if stdout_msgs:
        print 'death_rate = min(glc_death = {:4f},fru_death = {:4f}) = {:4f}'.format(mu_death_glc_D,mu_death_fru,death_rate)

    # Name of the output file storing the results
    if results_filename_base != '' and results_filename == '':
        results_filename = results_filename_base + '_' + str(leaky_traits_num) + 'LT_' + str(start_pos) + '_' + str(end_pos) + '.py'

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
        input_data['sucrose_uptake_rate'] = sucrose_uptake_rate
        input_data['his_uptake_rate'] = his_uptake_rate
        input_data['glc_uptake_rate'] = glc_uptake_rate
        input_data['fru_uptake_rate'] = fru_uptake_rate
        input_data['o2_uptake_rate'] = o2_uptake_rate
        input_data['relax_ATPM'] = relax_ATPM
        input_data['simulate_his'] = simulate_his
        input_data['fixed_o2'] = fixed_o2
        input_data['death_rate'] = death_rate
        input_data['simulate_rep_dynamics'] = simulate_rep_dynamics
        input_data['t0'] = 50
        input_data['tf'] = tf
        input_data['dt'] = dt
        input_data['save_details'] = save_details
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

    print '\nJob finished successfully!\n' 

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
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
        fba_reference_soln = iAZ900.fba_model.solution

    # For atp cost, do not incoporate IGPDH mutation and compute the biomass flux for the reference strain by 
    # without incorporaing atp into SUCRe
    else: 
        # compute the biomass flux for the reference strain by allowing for an excess amount of histidine uptake 
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
        fba_reference_soln = iAZ900.fba_model.solution

    #--- Now compute the biomass flux for the case cooperation cost (histidine uptake rate or atp in SUCRe is incoporated
    if simulate_his:
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[-o2_uptake_rate,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate,1000],'EX_fru_e_':[-fru_uptake_rate,1000],'EX_o2_e_':[(-10 - glc_uptake_rate - fru_uptake_rate)/5,1000], 'EX_his_L_e_':[-his_uptake_rate,1000],'IGPDH':[0,0]})
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
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
        iAZ900.fba(stdout_msgs = stdout_msgs, build_new_optModel = False)
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


def plot_gameResults(results_filename, x, y, x_key, y_key, title = '', xaxis_label = '', yaxis_label = '', custom_xaxis_ticklabels = None, custom_yaxis_ticklabels = None, set_minor_xticks = True, invert_xaxis = True, mixed_NashEq_label = 'Mixed', integrate_unknown_NashEq = True, output_filename_base = '', output_filetype = 'pdf'):
    """
    Plots a heatmap of the anaylsis results.
 
    INPUTS:
    -------
    results_filename: 
    A string containing the name of the file containing the results
    (use integrate_results_files function to integrate the results
    in different files

    x & y: 
    A list containing the elements of horizontal and vertical axes 

    x_key & y_key: 
    Name of the x and y variables in results (Note that results (stored in file_name) is with keys as follows:
    (('capture_eff', 0.2), ('atp_coeff', 0.115), ('his', 0), ('glc', 0), ('fru', 0))
    For example, if x if histidine uptake rate and y is capture efficiency, then x_key_name is his
    and y_key_name is capture_eff. 

    output_filename_base: 
    Base name of the output file containing the plot. This funciton then creates 
    two files output_filename_base + '_games.pdf' and output_filename_base + '_freq.pdf'
    The first stores the plot for the games and the second for species frequencies. 

    output_filetype:
    File type of the output file (e.g., pdf, png, etc). Do not include dot "." here
    """
    from tools.utilities.load_data_fromFile import load_data_from_python_file
    from matplotlib import colors
    from tools.utilities.plot import plot, axis, color_bar

    if not isinstance(set_minor_xticks,bool):
        raise TypeError('set_minor_xticks must be boolean')

    #--- Load the data ---
    results = load_data_from_python_file(file_name = results_filename, var_names = ['results'])['results']

    #--- Types of Nash equilibria ---
    # Find the Nash equilibria as follows:
    # Possible NE are as follows
    # 0)  ('Defector','Defector') --> Prisoner's Dilemma 
    # 1)  ('Cooperator','Cooperator') --> Mutually beneficial
    # 3)  (('Cooperator','Cooperator'),('Defector','Defector')) --> Harmony
    # 3)  ('Cooperator','Defector') --> Snowdrift
    print '\n--- Nash equilibria ---'
    Nash_eq_allResults = [tuple(results[res_key]['game_info']['NashEq']) for res_key in results.keys()]
    #Nash_eq = list(set([tuple(results[res_key]['game_info']['NashEq']) for res_key in results.keys()]))
    Nash_eq = list(set(Nash_eq_allResults))

    # 0 = PD --> red  1 = MB --> green  2 = SD --> cyan  3 = HR --> blue  
    # Source: http://stackoverflow.com/questions/30893483/make-a-heatmap-with-a-specified-discrete-color-mapping-with-matplotlib-in-python
    pd = (('Defector','Defector'),)                                # Red
    mb = (('Cooperator','Cooperator'),)                            # Green
    sd = (('Cooperator','Defector'),)                              # Cyan
    hr = (('Cooperator','Cooperator'),('Defector','Defector'))  # Blue

    known_Nash_eq = [pd, mb, sd, hr]
    unknown_Nash_eq = [n for n in Nash_eq if n not in known_Nash_eq]

    if integrate_unknown_NashEq:
        neq_color_map = {pd:'Red', mb:'Green', sd: 'Cyan', hr:'Blue', 'mixed':np.array([49,79,79])/255}
        neq_name_map = {pd:"Prisoner's\nDilemma", mb:'Mutually\nBeneficial', sd: 'Snowdrift', 'mixed':mixed_NashEq_label}
    else:
        neq_color_map = {pd:'Red', mb:'Green', sd: 'Cyan', hr:'Blue'}
        neq_name_map = {pd:"Prisoner's\nDilemma", mb:'Mutually\nBeneficial', sd: 'Snowdrift'}

        # Color and name for unknonw Nash equilibria    
        counter = 0
        gray_rgbs = [np.array([49,79,79])/255, np.array([105,105,105])/255, np.array([112, 138, 144])/255, np.array([190, 190, 190])/255, np.array([211, 211, 211])/255, np.array([238, 233, 233])/255, np.array([205, 201, 201])/255, np.array([220,220,220])/255, np.array([255,222,173])/255, np.array([238,203,173])/255, np.array([193,205,193])/255, np.array([100,149,237])/255]
        for neq in unknown_Nash_eq:
            neq_color_map[neq] = gray_rgbs[counter]
            neq_name_map[neq] = mixed_NashEq_label + str(counter + 1)
            counter += 1
            print neq_name_map[neq],' --> ',neq

    print '\nNash equilibria of the systems are:'
    print 'Known:'
    for neq in [n for n in Nash_eq if n in known_Nash_eq]:
        print neq,  ' --> ',Nash_eq_allResults.count(neq)
    print '\nUnknown:'
    for neq in unknown_Nash_eq:
        print neq, ' --> ',Nash_eq_allResults.count(neq)

    data_value = {}

    # Number of known Nash equilibria (m1m1,m1m2, m1m2, wtwt, m1wt,m2wt)
    known_neq = [kneq for kneq in known_Nash_eq if kneq in Nash_eq]
    known_NE_num = len(known_neq) 

    if integrate_unknown_NashEq:

        if known_NE_num == 0:
            colormap = colors.ListedColormap([neq_color_map['mixed']])
            colorbar_ticklabels = [neq_name_map['mixed']]
            data_value['mixed'] = 0
        else:
            # Find which known Nash equilibria appear in Nash_eq
            if len(Nash_eq) == known_NE_num:
                colormap = colors.ListedColormap([neq_color_map[kneq] for kneq in known_neq])
                colorbar_ticklabels = [neq_name_map[kneq] for kneq in known_neq]
                for kneq in known_neq:
                    data_value[neq_name_map[kneq]] = known_neq.index(kneq)

            elif len(Nash_eq) > known_NE_num:
                colormap = colors.ListedColormap([neq_color_map[kneq] for kneq in known_neq] + [neq_color_map['mixed']])
                colorbar_ticklabels = [neq_name_map[kneq] for kneq in known_neq] + [neq_name_map['mixed']]
                for kneq in known_neq:
                    data_value[neq_name_map[kneq]] = known_neq.index(kneq)
                data_value['mixed'] = len(known_neq)
            else:
                 raise userError('len(Nash_eq) = {} is less than known_NE_num = {}'.format(len(Nash_eq), known_NE_num))
    else:
        colormap = colors.ListedColormap([neq_color_map[neq] for neq in Nash_eq])
        colorbar_ticklabels = [neq_name_map[neq] for neq in Nash_eq]
        for neq in Nash_eq:
           data_value[neq_name_map[neq]] = Nash_eq.index(neq)

    print '\ndata_value = ', data_value
 
    dataG = np.zeros((len(y),len(x)))

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

            if integrate_unknown_NashEq:
                # Find out whether results_neq[res_key] is a known Nash equilibirum
                neqs = [kneq for kneq in known_neq if kneq == tuple(results[res_key]['game_info']['NashEq'])]
                if len(neqs) == 1:
                    dataG[i,j] = data_value[neq_name_map[neqs[0]]]
                elif len(neqs) == 0:
                    dataG[i,j] = data_value['mixed']
                else:
                    raise userError('len(neqs) = {} is greater than 1'.format(len(neqs)))
            else:
                neqs = [kneq for kneq in Nash_eq if kneq == tuple(results[res_key]['game_info']['NashEq'])]
                if len(neqs) == 1:
                    dataG[i,j] = data_value[neq_name_map[neqs[0]]]
                else:
                    raise userError('len(neqs) = {} is not 1'.format(len(neqs)))

    print '\nmin(data) = {}  , max(data) = {}'.format(dataG.min(), dataG.max())

    if integrate_unknown_NashEq:
       output_filename = output_filename_base + '_NEs.' + output_filetype
    else:
        output_filename = output_filename_base + '_NEs_notMixed.' + output_filetype

    curr_plt = plot(title = title, xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, limits = (min(x),max(x)), plot_gridlines = False, invert = invert_xaxis), yaxis = axis(label = yaxis_label,  set_minorticks = True, minorticks_spacing = 10, limits = (100*np.array(y).min(),100*np.array(y).max()), invert = False), fig_format = {'figsize':(7,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = np.array(x), y = 100*np.array(y), data = dataG, plot_func = 'pcolor', clrbar = color_bar(colormap = colormap, colorlimits = (dataG.min() ,dataG.max()), label = '',  custom_ticklabels =  colorbar_ticklabels, ticklabels_format = {'ticks_position':'middle'}))

    print 'The figure was saved into {}\n'.format(output_filename)

    #--------------- Evolutionary dynamics ---------------------------
    # Species frequences when wild-type invades (dataX_C) or mutants invade (dataX_D)
    dataX_C = np.zeros((len(y),len(x)))
    dataX_D = np.zeros((len(y),len(x)))

    # Here, x is the amino acid (or cost) and y is capture efficiency
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

            # Cooperator invades
            if (('Defector', 0.99), ('Cooperator', 0.01)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_C[i,j] = results[res_key]['replicator_dynamics_info'][(('Defector', 0.99), ('Cooperator', 0.01))]['Cooperator']
            elif (('Cooperator', 0.01), ('Defector', 0.99)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_C[i,j] = results[res_key]['replicator_dynamics_info'][(('Cooperator', 0.01), ('Defector', 0.99))]['Cooperator']
            elif (('Defector', 1),) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_C[i,j] = 0
            else:
                raise userError("Neither {} nor {} is in results[res_key]['replicator_dynamics_info'].keys()".format((('Defector', 0.01), ('Cooperator', 0.99)), (('Cooperator', 0.99), ('Defector', 0.01))))

            # Defector invades
            if (('Defector', 0.01), ('Cooperator', 0.99)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_D[i,j] = results[res_key]['replicator_dynamics_info'][(('Defector', 0.01), ('Cooperator', 0.99))]['Cooperator']
            elif (('Cooperator', 0.99), ('Defector', 0.01)) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_D[i,j] = results[res_key]['replicator_dynamics_info'][(('Cooperator', 0.99), ('Defector', 0.01))]['Cooperator']
            elif (('Defector', 1),) in results[res_key]['replicator_dynamics_info'].keys():
                dataX_C[i,j] = 0
            else:
                raise userError("Neither {} nor {} is in results[res_key]['replicator_dynamics_info'].keys()".format((('Defector', 0.01), ('Cooperator', 0.99)), (('Cooperator', 0.99), ('Defector', 0.01))))

    print 'min(dataX_C) = {}  , max(dataX_C) = {}'.format(dataX_C.min(), dataX_C.max())
    print 'min(dataX_D) = {}  , max(dataX_D) = {}'.format(dataX_D.min(), dataX_D.max())

    custom_colorbar_ticks = np.arange(0,1+0.2,0.2)

    output_filename = output_filename_base + '_freq_Cinvades.' + output_filetype
    curr_plt = plot(title = title + ' Cooperator invades', xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, custom_ticklabels = custom_xaxis_ticklabels, plot_gridlines = False, invert = invert_xaxis), yaxis = axis(label = yaxis_label, custom_ticklabels = custom_yaxis_ticklabels, invert = True), plot_gridlines = False, fig_format = {'figsize':(7,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = np.array(x), y = 100*np.array(y), data = dataX_C, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_colorbar_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
    print 'The figure was saved into {}\n'.format(output_filename)

    output_filename = output_filename_base + '_freq_Dinvades.' + output_filetype
    curr_plt = plot(title = title + ' Defector invades', xaxis = axis(label = xaxis_label, ticklabels_format = {'fontsize':20, 'position':'bottom','rotation':0}, custom_ticklabels = custom_xaxis_ticklabels, plot_gridlines = False, invert = invert_xaxis), yaxis = axis(label = yaxis_label, custom_ticklabels = custom_yaxis_ticklabels, invert = True), plot_gridlines = False,fig_format = {'figsize':(7,4), 'mathtext_fontname':'Arial'}, output_filename = output_filename)
    curr_plt.heatmap(x = np.array(x), y = 100*np.array(y), data = dataX_D, plot_func = 'matshow', clrbar = color_bar(colorlimits = (0,1), custom_ticks = custom_colorbar_ticks, label = 'Frequency', label_format = {'distance_from_ticklabels':30}), interpolate = True)
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

def create_jobfiles(create_cost_his_jobs = False, create_cost_atp_jobs = False, create_games_his_jobs = False, create_games_atp_jobs = False, create_all = False, glucose_uptake_rate = None):
    """
    Creates the job files given:
    task: The task that should be performed ('analyze_games' or 'analyze_cost')
    """
    from tools.utilities.create_job_files import create_intervals, create_job_files

    if create_all:
        create_cost_his_jobs = True
        create_cost_atp_jobs = True
        create_games_his_jobs = True
        create_games_atp_jobs = True

    sucrose_uptake_rate = 10
    o2_uptake_rate = sucrose_uptake_rate/5
    capture_efficiency = [i/50 for i in range(51)]
    capture_efficiency_str = '[i/50 for i in range(51)]'
    histidine_uptake_rate = [i/200 for i in range(21)] 
    histidine_uptake_rate_str = '[i/200 for i in range(21)]'
    #SUCRe_atp = [i/10 for i in range(0,301,3)] 
    #SUCRe_atp_str = '[i/10 for i in range(0,301,3)]' 
    #SUCRe_atp = [i/100 for i in range(0,3001,25)] 
    #SUCRe_atp_str = '[i/100 for i in range(0,3001,25)]' 
    SUCRe_atp = [i/10 for i in range(0,301,5)] 
    SUCRe_atp_str = '[i/10 for i in range(0,301,5)]' 

    # Total combinations for simulating histidine
    total_comb_his = len(capture_efficiency)*len(histidine_uptake_rate) 
    total_comb_atp = len(capture_efficiency)*len(SUCRe_atp) 
    his_interval_size = 600
    atp_interval_size = 500
    print '\ntotal_comb_his = {} , his_interval_size = {}     total_comb_atp = {} , atp_interval_size = {}\n'.format(total_comb_his,his_interval_size,total_comb_atp,atp_interval_size)

    #--- Create the jobs for cost analysis ---
    if create_cost_his_jobs:
        job_filename_base = 'jobs/job_cost_his'
        joboutput_filename_base = 'job_cost_his'
        results_filename_base = 'results/cost_results_his'

        commands = ['python -c "from __future__ import division;from simulate_games import master_func;master_func(start_pos = None , end_pos = None, capture_efficiency = ' + capture_efficiency_str + ", histidine_uptake_rate = " + histidine_uptake_rate_str + ", SUCRe_atp = 0.115, simulate_his = True, fixed_o2 = True, o2_uptake_rate = 2, task = 'analyze_cost', stdout_msgs = False, warnings = True, results_filename = None)" + '"']

        create_job_files(total_cases_num = total_comb_his, interval_size = his_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos = None', end_pos_input_format = 'end_pos = None', main_code_dir = '/usr2/postdoc/alizom/work/yeast_sucrose', commands = commands, max_walltime = max_walltime)

    if create_cost_atp_jobs:
        job_filename_base = 'jobs/job_cost_atp' 
        joboutput_filename_base = 'job_cost_atp' 
        results_filename_base = 'results/cost_results_atp'

        commands = ['python -c "from __future__ import division;from simulate_games import master_func;master_func(start_pos = None , end_pos = None, capture_efficiency = ' + capture_efficiency_str + ', histidine_uptake_rate = 0, SUCRe_atp = ' + SUCRe_atp_str + ", simulate_his = False, fixed_o2 = True, o2_uptake_rate = 2, task = 'analyze_cost', stdout_msgs = False, warnings = True, results_filename = None)" + '"']

        create_job_files(total_cases_num = total_comb_atp, interval_size = atp_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos = None', end_pos_input_format = 'end_pos = None', main_code_dir = '/usr2/postdoc/alizom/work/yeast_sucrose', commands = commands, max_walltime = max_walltime)

    #--- Create the jobs for game analysis ---
    if create_games_his_jobs or create_games_atp_jobs:
        max_walltime = 48
        # Additions to the file names
        if glucose_uptake_rate == 0:
            if sucrose_uptake_rate == 10: 
                filename_additions = '_noGlc'
            else:
                filename_additions = '_noGlc_lowSuc'
        elif glucose_uptake_rate > 0: 
            if sucrose_uptake_rate == 10: 
                filename_additions = '_glc'
            else:
                filename_additions = '_glc_lowSuc'
        else:
            raise userError('Provide a non-zero value for glucose_uptake_rate')

    if create_games_his_jobs:
        max_walltime = 48
        job_filename_base = 'jobs/job_ss_games_his' + filename_additions 
        joboutput_filename_base = 'job_ss_games_his' + filename_additions 
        results_filename_base = 'results/ss_game_results_his' + filename_additions

        commands = ['python -c "from __future__ import division;from simulate_games import master_func;master_func(start_pos = None , end_pos = None, capture_efficiency = ' + capture_efficiency_str + ", sucrose_uptake_rate = " + str(sucrose_uptake_rate) + ", o2_uptake_rate = " + str(o2_uptake_rate) + ", histidine_uptake_rate = " + histidine_uptake_rate_str + ', SUCRe_atp = 0.115, glucose_uptake_rate = ' + str(glucose_uptake_rate) + ", simulate_his = True, fixed_o2 = True, task = 'analyze_games', simulate_rep_dynamics = True, dt = 1, tf = 5000, save_details = False, stdout_msgs = False, warnings = True, results_filename = None)" + '"']

        create_job_files(total_cases_num = total_comb_his, interval_size = his_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos = None', end_pos_input_format = 'end_pos = None', main_code_dir = '/usr2/postdoc/alizom/work/yeast_sucrose', commands = commands, max_walltime = max_walltime)

    if create_games_atp_jobs:
        max_walltime = 48
        job_filename_base = 'jobs/job_ss_games_atp' + filename_additions
        joboutput_filename_base = 'job_ss_games_atp' + filename_additions
        results_filename_base = 'results/ss_game_results_atp' + filename_additions

        commands = ['python -c "from __future__ import division;from simulate_games import master_func;master_func(start_pos = None , end_pos = None, capture_efficiency = ' + capture_efficiency_str + ", sucrose_uptake_rate = " + str(sucrose_uptake_rate) + ", o2_uptake_rate = " + str(o2_uptake_rate) +  ', histidine_uptake_rate = 0, SUCRe_atp = ' + SUCRe_atp_str + ", glucose_uptake_rate = " + str(glucose_uptake_rate) + ", simulate_his = False, fixed_o2 = True, task = 'analyze_games', simulate_rep_dynamics = True, dt = 1, tf = 5000, save_details = False, stdout_msgs = False, warnings = True, results_filename = None)" + '"']

        create_job_files(total_cases_num = total_comb_atp, interval_size = atp_interval_size, job_filename_base = job_filename_base, joboutput_filename_base = joboutput_filename_base, results_filename_base = results_filename_base, results_filename_input_format = 'results_filename = None', start_pos_input_format = 'start_pos = None', end_pos_input_format = 'end_pos = None', main_code_dir = '/usr2/postdoc/alizom/work/yeast_sucrose', commands = commands, max_walltime = max_walltime)

def run_create_job_files():
    """
    It just runs create_job_files to create pbs files for different cases
    """
    create_games_his_jobs = True
    create_games_atp_jobs = True

    # Create jobs to analyze games for simulations iwth histidne uptake rate 
    if create_games_his_jobs:
        create_jobfiles(create_games_his_jobs = True, glucose_uptake_rate = 0)
        create_jobfiles(create_games_his_jobs = True, glucose_uptake_rate = 5)

    # Create jobs to analyze games for simulations iwth histidne SUCRe atp cost 
    if create_games_atp_jobs:
        create_jobfiles(create_games_atp_jobs = True, glucose_uptake_rate = 0)
        create_jobfiles(create_games_atp_jobs = True, glucose_uptake_rate = 5)

def run_plot_results():
    """
    Plots the results
    """
    capture_efficiency = [i/50 for i in range(51)]
    histidine_uptake_rate = [i/200 for i in range(21)] 
    #SUCRe_atp = [i/10 for i in range(0,301,3)] 
    SUCRe_atp = [i/10 for i in range(0,301,5)] 

    cost_his = False
    cost_atp = False

    game_his_noGlc = True
    game_atp_noGlc = True

    game_his_glc = True
    game_atp_glc = True

    game_his_noGlc_lowSuc = False
    game_atp_noGlc_lowSuc = False

    game_his_glc_lowSuc = False
    game_atp_glc_lowSuc = False

    # ----- Plot the cost results -----
    # Use \rm to remove italic from matplotlib LaTex text.
    # Source: http://stackoverflow.com/questions/19671659/remove-italics-in-latex-subscript-in-matplotlib
    # You can also use \mathregular{} to make the font regular
    #-- his --
    if cost_his: 
        plot_costResults(results_filename = 'results/cost_results_his_fixedo2.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = r'Fixed limiting O$_{\rm 2}$', xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', x_minorticks_spacing = None, cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'Histidine production cost', output_filename = 'results/figures/cost_results_his_fixedo2.pdf')

    #-- ATP --
    if cost_atp: 
        plot_costResults(results_filename = 'results/cost_results_atp_fixedo2.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = r'Fixed limiting O$_{\rm 2}$', xaxis_label = 'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, cost10_replacement_value = 0.55, cost15_replacement_value = 0.6, colorbar_label = 'SUCRe ATP Cost', output_filename = 'results/figures/cost_results_atp_fixedo2.pdf')

    # ------ Plot the game results ------
    # custom_xaxis_ticklabels and custom_yaxis_ticklabels are used only to plot fractions where use matshow

    # -- his --
    if game_his_noGlc:
        print '\n--------------- Histin games, no glc ------------\n'
        custom_xaxis_ticklabels = [str(xx) if xx/0.02 == int(xx/0.02) else '' for xx in histidine_uptake_rate]
        custom_yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) else '' for yy in capture_efficiency]
        title = 'no glc' 
        title = '' 
        plot_gameResults(results_filename = 'results/ss_game_results_his_noGlc.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = title, xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = True, integrate_unknown_NashEq = True, custom_xaxis_ticklabels = custom_xaxis_ticklabels, custom_yaxis_ticklabels = custom_yaxis_ticklabels, output_filename_base = 'results/figures/ss_game_results_his_noGlc', output_filetype = 'pdf')

    if game_his_glc:
        print '\n--------------- Histin games, glc ------------\n'
        title = 'glc' 
        title = '' 
        custom_xaxis_ticklabels = [str(xx) if xx/0.02 == int(xx/0.02) else '' for xx in histidine_uptake_rate]
        custom_yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) else '' for yy in capture_efficiency]
        plot_gameResults(results_filename = 'results/ss_game_results_his_glc.py', x = histidine_uptake_rate, y = capture_efficiency, x_key = 'his', y_key = 'capture_eff', title = title, xaxis_label = r'Histidine uptake rate $\left(\frac{mmol}{gDW.h}\right)$', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = True, integrate_unknown_NashEq = True, custom_xaxis_ticklabels = custom_xaxis_ticklabels, custom_yaxis_ticklabels = custom_yaxis_ticklabels, output_filename_base = 'results/figures/ss_game_results_his_glc')

    #-- ATP --
    if game_atp_noGlc:
        print '\n--------------- ATP games, no glc ------------\n'
        title = 'no glc'
        title = '' 
        custom_xaxis_ticklabels = [str(int(xx)) if xx/6 == int(xx/6) else '' for xx in SUCRe_atp]
        custom_yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) else '' for yy in capture_efficiency]
        plot_gameResults(results_filename = 'results/ss_game_results_atp_noGlc.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = title, xaxis_label = r'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = False, integrate_unknown_NashEq = True, custom_xaxis_ticklabels = custom_xaxis_ticklabels, custom_yaxis_ticklabels = custom_yaxis_ticklabels, output_filename_base = 'results/figures/ss_game_results_atp_noGlc')

    if game_atp_glc:
        print '\n--------------- ATP games, glc ------------\n'
        title = 'glc'
        title = '' 
        custom_xaxis_ticklabels = [str(int(xx)) if xx/3 == int(xx/3) else '' for xx in SUCRe_atp]
        custom_yaxis_ticklabels = [str(int(100*yy)) if 100*yy/20 == int(100*yy/20) else '' for yy in capture_efficiency]
        plot_gameResults(results_filename = 'results/ss_game_results_atp_glc.py', x = SUCRe_atp, y = capture_efficiency, x_key = 'atp_coeff', y_key = 'capture_eff', title = title, xaxis_label = r'mol of ATP', yaxis_label = 'Capture efficiency (%)', set_minor_xticks = False, invert_xaxis = False, integrate_unknown_NashEq = True, custom_xaxis_ticklabels = custom_xaxis_ticklabels, custom_yaxis_ticklabels = custom_yaxis_ticklabels, output_filename_base = 'results/figures/ss_game_results_atp_glc')


def run_integrate_resultsfiles():

    from tools.utilities.integrate_results_files import integrate_results_files, del_results_files

    games_his_noGlc = True
    games_his_glc = True
    games_atp_noGlc = True
    games_atp_glc = True

    games_his_noGlc_lowSuc = False
    games_his_glc_lowSuc = False
    games_atp_noGlc_lowSuc = False
    games_atp_glc_lowSuc = False
    
    #---- Games, his, no glc -----
    if games_his_noGlc: 
        results_filenames = ['results/ss_game_results_his_noGlc_1_600.py','results/ss_game_results_his_noGlc_601_1071.py']

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_his_noGlc.py')

        #del_results_files(results_filenames)

    if games_his_glc: 
        results_filenames = ['results/ss_game_results_his_glc_1_600.py','results/ss_game_results_his_glc_601_1071.py']

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_his_glc.py')

        #del_results_files(results_filenames)

    if games_his_noGlc_lowSuc: 
        results_filenames = ['results/ss_game_results_his_noGlc_lowSuc_1_600.py','results/ss_game_results_his_noGlc_lowSuc_601_1071.py']

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_his_noGlc_lowSuc.py')

        #del_results_files(results_filenames)

    if games_his_glc_lowSuc: 
        results_filenames = ['results/ss_game_results_his_glc_lowSuc_1_600.py','results/ss_game_results_his_glc_lowSuc_601_1071.py']

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_his_glc_lowSuc.py')

        #del_results_files(results_filenames)


    #---- Games, atp, no glc -----
    if games_atp_noGlc: 
        results_filenames = ['results/ss_game_results_atp_noGlc_1_500.py', 'results/ss_game_results_atp_noGlc_501_1000.py', 'results/ss_game_results_atp_noGlc_1001_1500.py', 'results/ss_game_results_atp_noGlc_1501_2000.py', 'results/ss_game_results_atp_noGlc_2001_2500.py', 'results/ss_game_results_atp_noGlc_2501_3000.py', 'results/ss_game_results_atp_noGlc_3001_3111.py' ]

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_atp_noGlc.py')

        #del_results_files(results_filenames)

    if games_atp_glc: 
        results_filenames = ['results/ss_game_results_atp_glc_1_500.py', 'results/ss_game_results_atp_glc_501_1000.py', 'results/ss_game_results_atp_glc_1001_1500.py', 'results/ss_game_results_atp_glc_1501_2000.py', 'results/ss_game_results_atp_glc_2001_2500.py', 'results/ss_game_results_atp_glc_2501_3000.py', 'results/ss_game_results_atp_glc_3001_3111.py' ]

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_atp_glc.py')

        #del_results_files(results_filenames)

    if games_atp_noGlc_lowSuc: 
        results_filenames = ['results/ss_game_results_atp_noGlc_lowSuc_1_500.py', 'results/ss_game_results_atp_noGlc_lowSuc_501_1000.py', 'results/ss_game_results_atp_noGlc_lowSuc_1001_1500.py', 'results/ss_game_results_atp_noGlc_lowSuc_1501_2000.py', 'results/ss_game_results_atp_noGlc_lowSuc_2001_2500.py', 'results/ss_game_results_atp_noGlc_lowSuc_2501_3000.py', 'results/ss_game_results_atp_noGlc_lowSuc_3001_3111.py' ]

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_atp_noGlc_lowSuc.py')

        #del_results_files(results_filenames)

    if games_atp_glc_lowSuc: 
        results_filenames = ['results/ss_game_results_atp_glc_lowSuc_1_500.py', 'results/ss_game_results_atp_glc_lowSuc_501_1000.py', 'results/ss_game_results_atp_glc_lowSuc_1001_1500.py', 'results/ss_game_results_atp_glc_lowSuc_1501_2000.py', 'results/ss_game_results_atp_glc_lowSuc_2001_2500.py', 'results/ss_game_results_atp_glc_lowSuc_2501_3000.py', 'results/ss_game_results_atp_glc_lowSuc_3001_3111.py' ]

        integrate_results_files(results_filenames = results_filenames, results_varname = 'results', output_filename = 'results/ss_game_results_atp_glc_lowSuc.py')

        #del_results_files(results_filenames)


#-----------------------------------
if __name__ == '__main__':
    print '\n------- Test -------- '
    # Create the payoff matrix of the game
    #sucrose_game = master_func(capture_efficiency = 0.3, histidine_uptake_rate = 0.02, sucrose_uptake_rate = 10, o2_uptake_rate = 10/5, glucose_uptake_rate = 0, simulate_his = True, fixed_o2 = True, task = 'analyze_games', save_details = True,  warnings = True, stdout_msgs = True, results_filename = 't.py')
    sucrose_game = master_func(capture_efficiency = 0.0, histidine_uptake_rate = 0.0, SUCRe_atp = 0, sucrose_uptake_rate = 10, o2_uptake_rate = 2, simulate_his = False, fixed_o2 = True, task = 'analyze_games', save_details = True,  warnings = True, stdout_msgs = True, results_filename = 't.py')
   
