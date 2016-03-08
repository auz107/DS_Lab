from __future__ import division
import sys, time
sys.path.append('../')
import numpy as np
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.compartment import compartment
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from DMMM_yeast import DMMM_yeast
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.GAMETES.game import game
from tools.userError import userError
from multiprocessing import Process, Manager
from copy import deepcopy
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

"""
This module contains the following functions:

      assign_uptakeKinetics: Assigns uptake kinetics for shared compounds to the desired reactions in the metabolic model
             findHisSatConc: Finds the saturaated concetration of histidine according to uptake kinetics
             findSucSatConc: Finds the saturaated concetration of sucrose according to uptake kinetics
    find_min_sucrose_uptake: Finds the minimum required sucrose uptake for model to be feasible
plot_min_suc_uptake_results: Plots the results of plot_min_suc_uptake_results
               create_model: Create a metabolic model
                performDMMM: Performs Dynamics multi-species metabolic modeling (DMMM)
          master_func: Runs performDMM for various conditions
      run_master_func: Runs master_func for various input parameters
  plot_dynamic_game_results: Plots the results of run_master_func() 
           create_intervals: Create intervals for creating pbs jobs
           create_job_files: Creates pbs job files
       run_create_job_files: Runs create_job_files
 coopr_density_effect_on_mu: Investigates the impact of cooperators' density on their growth rate
               plot_density: Plots the results of coopr_density_effect_on_mu
       coopr_cheater_invade: Investigates the invation of a population of cooperators by cheaters and vice versa
              test_dynamics: Tests master_func
"""

def assign_uptakeKinetics(model):
    """
    Assigns uptake kinetics for glucose, fructose, sucrose and histidine

    Sucrose:
    -------
    1. (http://www.sciencedirect.com/science/article/pii/0003986182902557):
       Km = 6 mM ,  Vmax = 2 nano-mol/(mg DW.min) --> see page 655, left column
       It should be noted that these constants do not have the significance of K, and V as 
       there is probably more than one system involved in the transport of sucrose. In fact,
       there are two different kinetics components with very different apparent affinities, 
       at higher concentration of substrate. The point of inflection showing the separation 
       between both components appears around 13 mM sucrose. Thus, it seems that a simple 
       Michaelis-Menten kinetics does not apply for this transport system.  
    2. (http://www.sciencedirect.com/science/article/pii/0922338X94901546): 
       The sucrose transport system does not obey the basic Michaelis-Menten kinetics and consists of at least 
       two transport mechanisms for concentrations below and above 12-13 mM. For this transport system, 
       Therefore, it is very difficult to obtain the usual saturation kinetics parameters of
    3. (http://www.sciencedirect.com/science/article/pii/0922338X9685029X):
       Km = 6.41 mM  , Vmax = 18.36 mmol/(mg DW.min) (see Table 3)
    4. (PMID: 16232731):
       The kinetic analysis of active sucrose-H+ uptake by Saccharomyces cerevisiae revealed the presence 
       of two transport systems with high and low af6nlty for sucrose. The MALZT permease has a low aliinity 
       (Km= 120 +/- 20 mM) for sucrose, while the a-glucoside transporter encoded by the AGTl gene is a 
       high aflinlty sucrose-H+ symporter (Km = 7.95 +/- 0.8 mM) that increases the specificgrowth rate of 
       cellsgrowing on sucrose.       

    Assuming that Vmax is the same for both low- and high- affinity systems:
    Vmax = 2 nmol/(mg DW.min)*(10^-6 mmol/1 nmol)*(1000 mg/1 g)*(60 min/1 h) = 0.12 mmol/(gDW.h) 
    Vmax = 18.36 mmol/(mg DW.min)*(1000 mg/1 g)*(60 min/1 h) = 1101600.0 mmol/(gDW.h) --> Unlikely

    According to the papers above, the uptake kinetics of sucrose is as follows:
    if C <= (12 + 13)/2 = 12.5 mM,  Km = 120 mM
    if C > 12.5 mM,                 Km = 7.95 mM                

    Glucose:
    --------
   1.(PMID: 10336421):
    Table II: Protein content: 0.42 g/g DW
    Table III (Exponential): Km = 18 mM  , Vmax = 16 - 19 (nmol/mg Protein.min)
    Vmax = (16 + 19)/2 = 17.5(nmol/mg Protein.min)*(10^-6 mmol/1 nol)*(0.42 mg Protein/1 mg DW)*
                         (1000 mg DW/1 gDW)*(60 min/1 h) = 0.441 mmol/(gDW.h)

    2. (http://www.pnas.org/content/80/6/1730.full.pdf):
       High affinity: Km = 1.5 +/- 0.25 mM  and Low affinity: 20 +/- 8
    
    According to these papers and assuming Vmax is the same for both low and high affinity enzymes:
    Vmax = 0.441 mmol/gDW.h and Km = 1.5 mM and Km = 20 mM

    Fructose:
    ---------
    1. (http://www.pnas.org/content/80/6/1730.full.pdf):
       High affinity: Km = 6 +/- 2 mM  and Low affinity: 40 +/- 15
    Take Vmax as that for glucose Vmax = 0.441 mmol/(gDW.h)


    Histidine:
    ----------
    According to PMID: 9473505 the kinetics of histidine uptake in yeast is as follows: 
    V = Vmax_HIP1*S/(Km_HIP1 + S) + Vmax_TAT1/(Km_TAT1 + S)
    where,
    Vmax_HIP1 = 44 pmol/min/(10^6 cells)  ,  Km_HIP1 = 17 muM
    Vmax_TAT1 = 35 pmol/min/(10^6 cells)  ,  Km_TAT1 = 370 muM

    Cell dry weight is about 4 times the buoyant mass (PMID: 20383132).
    the buoyant mass of yeast ranges from 3-10pg (pico = 10^-12)
    (Source: http://www.weizmann.ac.il/plants/Milo/images/YeastSize-Feb2010.pdf)
    Therefore, the dry cell weight ~ 40*10^-12 = 4*10^-11 g/cell
    1 pico = 10^-9 mili 
    Convert pmol/min/(10^6 cells) to mmol/gDW.h:
    pmol/(min*10^6 cells) * (1 mmol/10^9 pmol)*(60 min/1 h)*(1 cell/(4*10^-11 g)) 
    pmol/(min*10^6 cells) * (15*10^-4) = 1.5e-3 --> Conversion factor
    Therefore,
    Vmax_HIP1 = 44*1.5e-3 = 0.066 mmol/(gDW.h)  ,  Km_HIP1 = 0.017 mM
    Vmax_TAT1 = 35*1.5-3 = 0.0525 mmol/(gDW.h)  ,  Km_TAT1 = 0.370 mM

    From Jeff Gore's paper the saturation concentration is 20 mu-g/ml. We want to convert it
    to mmol/L or mM. The molecular weight of histidine is 155.15 g/mol
    20 mu-g/ml *(10^-6 g/1 mu-g)*(1000 ml/1 l)*(1 mol/155.15 g)*(1000 mmol/1 mol) = 0.129 mmol/l 
    """
    # Exchange reactions for glucose, fructose, sucrose and histidine
    glc_exchrxn = model.get_reactions({'EX_glc_e_':'id'})
    fru_exchrxn = model.get_reactions({'EX_fru_e_':'id'})
    suc_exchrxn = model.get_reactions({'EX_sucr_e_':'id'})
    his_exchrxn = model.get_reactions({'EX_his_L_e_':'id'})
    
    # Glucose
    glc_e_id = glc_exchrxn.reactants[0].id
    glc_exchrxn.kinetics = "0.441*C['" + glc_e_id + "']/(1.5 + C['" + glc_e_id + "']) +" + "0.441*C['" + glc_e_id + "']/(20 + C['" + glc_e_id + "'])"

    # Fructose
    fru_e_id = fru_exchrxn.reactants[0].id
    fru_exchrxn.kinetics = "0.441*C['" + fru_e_id + "']/(6 + C['" + fru_e_id + "']) +" + "0.441*C['" + fru_e_id + "']/(40 + C['" + fru_e_id + "'])"

    # Sucrose
    affinity_func = lambda C: 1 if C <= 12.5 else 0
    suc_e_id = suc_exchrxn.reactants[0].id
    suc_exchrxn.affinity_func = affinity_func
    #suc_exchrxn.kinetics = "self.affinity_func(C['" + suc_e_id + "'])*0.12*C['" + suc_e_id + "']/(120 + C['" + suc_e_id + "']) + " + "(1 - self.affinity_func(C['" + suc_e_id + "']))*0.12*C['" + suc_e_id + "']/(7.95 + C['" + suc_e_id + "'])"
    suc_exchrxn.kinetics = "0.12*C['" + suc_e_id + "']/(120 + C['" + suc_e_id + "']) +  0.12*C['" + suc_e_id + "']/(7.95 + C['" + suc_e_id + "'])"

    # Histidine
    his_e_id = his_exchrxn.reactants[0].id
    his_exchrxn.kinetics = "0.066*C['" + his_e_id + "']/(0.017 + C['" + his_e_id + "']) +" + "0.052*C['" + his_e_id + "']/(0.370 + C['" + his_e_id + "'])"

def findHisSatConc(Smax = 100, S_step = 0.1, tol = 0.95, output_filename = '',stdout_msgs = True):
    """
    Finds the saturating histidine concentration of a strain lacking histidin synthesizing genes

    See assign_uptakeKinetics for details of uptake kinetics
    """    
    from tools.ancillary.plot import plot, axis

    Vmax = 0.066 + 0.0525  # According to uptake kinetics
    V = 0 
    S = 0
    S_vec, V_vec = [0], [0]  # Vectors containing S and V values
    while V < tol*Vmax and S < Smax:
         S += S_step 
         S_vec.append(S)
         V = 0.066*S/(0.017 + S) + 0.0525*S/(0.370 + S) 
         V_vec.append(V)

    S_Jeff = 0.129
    V_Jeff = 0.066*S_Jeff/(0.017 + S_Jeff) + 0.0525*S_Jeff/(0.370 + S_Jeff)

    if stdout_msgs:
        print 'Saturated histidine concentration according to uptake kinetics = {} --> Uptake rate = {}\n'.format(S,V)
        print "Reported saturated histidiine concentration in Gore's paper = {} --> Uptake rate = {}".format(S_Jeff, V_Jeff)

    # Plot the results
    linplot = plot(xaxis = axis(label = 'Substrate conc. (mM)', limits = (0,max(S_vec)), majorticks_spacing = 0.2, spines_format = {'bottom':{'linewidth':2},'top':{'linewidth':2}}), yaxis = axis(label = 'Uptake rate \n(mmol/gDW.h)', spines_format = {'left':{'linewidth':2}, 'right':{'linewidth':2}}), plot_gridlines = False, figsize = None, output_filename = output_filename)
 
    linplot.plot2D(x = S_vec, y = V_vec, sort_data = False, save_current = False)

    # Create the lines for Jeff Gore data
    linplot.plot2D(x = [0.12,0.12], y = [0,0.08], line_format = {'width':1,'color':'black'}, sort_data = False, save_current = False)
    linplot.add_text(text = 'Reported saturation \nconc. (' + str(S_Jeff) + ' (mM)', text_format = {'x_pos':0.15, 'y_pos':0.02, 'fontweight':'regular'})    
    linplot.customize_and_save()

    return S

def findSucSatConc(uptake_kinetics_type = 'conditional', Smax = 100, S_step = 0.1, tol = 0.95, output_filename = '', stdout_msgs = True):
    """
    INPUTS:
    ------
    uptake_kinetics_type: Can be conditional or combined, where
    conditional: V = 0.12*S/(120 + S), if S <= 12.5  &  V = 0.12*S/(7.95 + S),  otherwise 
      combined:  V = 0.12*S/(120 + S) +  0.12*S/(7.95 + S)    
    """
    from tools.ancillary.plot import plot, axis

    #--- FInd the maximum sucrose uptake rate according to kinetic equations ---
    if uptake_kinetics_type.lower() == 'conditional':
        Vmax = 0.12 # According to uptake kinetics
    elif uptake_kinetics_type.lower() == 'combined':
        Vmax = 0.12 + 0.12 # According to uptake kinetics
    else:
        raise ValueError('Invalid uptake_kinetics_type value!')

    Vold = 0
    V = 0 
    S = 0
    S_vec, V_vec = [0], [0]  # Vectors containing S and V values
    
    while V < tol*Vmax and S < Smax:
        Vold = V
        S += S_step 
        S_vec.append(S)
        if uptake_kinetics_type.lower() == 'conditional':
            if S <= 12.5:
                 V = 0.12*S/(120 + S) 
            else:
                 V = 0.12*S/(7.95 + S) 
        else:
             V = 0.12*S/(120 + S) + 0.12*S/(7.95 + S) 

        V_vec.append(V)

    S_sat = S # Saturated concentration of glucose 
    V_sat = V # Uptake rate at the saturated sucrose concentration

    S_Jeff = 1.46 
    if uptake_kinetics_type.lower() == 'conditional':
        V_Jeff = 0.12*S_Jeff/(120 + S_Jeff) 
    else:
        V_Jeff = 0.12*S_Jeff/(120 + S_Jeff) + 0.12*S_Jeff/(7.95 + S_Jeff) 

    if stdout_msgs:
        print 'Saturated sucrose concentration according to uptake kinetics = {} --> Uptake rate = {} ({:0.2f}% of Vmax = {})\n'.format(S_sat,V_sat, 100*V_sat/Vmax,Vmax)
        print "Reported used sucrose concentration in Gore's paper = {} --> Uptake rate = {} ({:0.2f}% of Vmax = {})".format(S_Jeff, V_Jeff, 100*V_Jeff/Vmax,Vmax)

    # Plot the results
    linplot = plot(xaxis = axis(label = 'Sucrose conc. (mM)', limits = (0,max(S_vec)), majorticks_spacing = None, spines_format = {'bottom':{'linewidth':2},'top':{'linewidth':2}}), yaxis = axis(label = 'Uptake rate \n(mmol/gDW.h)', spines_format = {'left':{'linewidth':2}, 'right':{'linewidth':2}}), plot_gridlines = False, figsize = None, output_filename = output_filename)
 
    linplot.plot2D(x = S_vec, y = V_vec, sort_data = False, save_current = False)

    # Create the lines for Jeff Gore data
    linplot.plot2D(x = [1.46,1.46], y = [0,0.06], line_format = {'width':1,'color':'black'}, sort_data = False, save_current = False)
    linplot.add_text(text = "conc. in Gore's \nexperiments (" + str(S_Jeff) + ' mM)', text_format = {'x_pos':2.5, 'y_pos':0.01, 'fontweight':'regular'})    
    linplot.customize_and_save()

def find_min_sucrose_uptake(uptake_kinetics_type = 'conditional', output_filename = '', stdout_msgs = True):
    """
    Finds the minimum requried sucrose uptake for the cell to growth after incorporating capture efficiency
    and uptake kinetics 
    """
    from tools.ancillary.plot import plot, axis

    #--- Metabolic models ---
    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/'

    # Create the metabolic model object
    iAZ900 = create_model()

    capture_eff_vec = []
    ms_vec = []
    
    #--- Incorporate atp into SUCRe ---
    # According to Jeff Gore's paper (suppl. figure 5) the cooperators are 2.5% less fit than cheasters
    # Note that the plot in this figure is for ln(X/X0)_C/ln(X/X0)_D but ln(X/X0) = mu*t or ln(X/X0)
    # is proportional to mu.
    # In the absence of a metabolic cost for sucrose hydrolysis (SUCRe) the growth rate for WT is 0.5031.
    # We next add an ATP cost to SUCRe reaction (for WT) and adjust its coefficinet such that the growth
    # becomes 97.5% of 0.5023 = 0.49056. The coefficient of atp was determined to be 0.115
    atp_coeff = 0.115
    SUCRe = iAZ900.get_reactions({'SUCRe':'id'})
    atp = iAZ900.get_compounds({'atp_c':'id'})
    adp = iAZ900.get_compounds({'adp_c':'id'})
    pi = iAZ900.get_compounds({'pi_c':'id'})
    SUCRe.add_compounds({atp:-atp_coeff,adp:atp_coeff,pi:atp_coeff})
    if stdout_msgs: 
        print '\n--- FBA for iAZ900 after adding atp/adp to SUCRe---'  
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-2,1000]})
    iAZ900.fba(stdout_msgs = stdout_msgs)

    #--- Compute ms for sucrose for a capture efficiency of one ---
    # Compute ms and biomass yield for sucrose, glucose, fructose and histidine
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_o2_e_':[-2,1000]})
    iAZ900.compounds_by_id['sucr_e'].ms_calc()
    ms_capture_eff_one = iAZ900.compounds_by_id['sucr_e'].ms

    #--- Incorporate capture efficiency ---
    capture_eff = 0.1
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
        print '\n--- FBA for modified iAZ900 (after adding glc_secreted and fru_secreted to SUCRe ---'  
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-2,1000]})
    iAZ900.fba(stdout_msgs = stdout_msgs)


    capture_eff = 0.00
    print '(capture_eff, ms) = ',
    while capture_eff < 0.99:
        capture_eff += 0.01
        capture_eff_vec.append(100*capture_eff)

        SUCRe.stoichiometry[glc_D_e] = capture_eff
        SUCRe.stoichiometry[fru_e] = capture_eff
        SUCRe.stoichiometry[glc_D_secreted] = 1 - capture_eff
        SUCRe.stoichiometry[fru_secreted] = 1 - capture_eff

        # Compute ms and biomass yield for sucrose, glucose, fructose and histidine
        set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_o2_e_':[-2,1000]})
        iAZ900.compounds_by_id['sucr_e'].ms_calc()
        ms_vec.append(iAZ900.compounds_by_id['sucr_e'].ms)
        print '({} ,{}), '.format(capture_eff,iAZ900.compounds_by_id['sucr_e'].ms),

    capture_eff_vec.append(100)
    ms_vec.append(ms_capture_eff_one)
    print '(1 ,{})\n'.format(ms_capture_eff_one)

    return capture_eff_vec, ms_vec

def plot_min_suc_uptake_results(capture_eff_vec, ms_vec, uptake_kinetics_type = 'conditional', output_filename = ''):
    # Plot the results
    from tools.ancillary.plot import plot, axis

    linplot = plot(xaxis = axis(label = 'Capture efficiency (%)', limits = None, majorticks_spacing = None, spines_format = {'bottom':{'linewidth':2},'top':{'linewidth':2}}), yaxis = axis(label = 'min sucrose uptake rate \n(mmol/gDW.h)', spines_format = {'left':{'linewidth':2}, 'right':{'linewidth':2}}), plot_gridlines = False, figsize = None, output_filename = output_filename)
 
    linplot.plot2D(x = capture_eff_vec, y = ms_vec, sort_data = False, save_current = False)

    # Create the lines for Jeff Gore data
    S_Jeff = 1.46 
    if uptake_kinetics_type.lower() == 'conditional':
        V_Jeff = 0.12*S_Jeff/(120 + S_Jeff) 
    elif uptake_kinetics_type.lower() == 'combined':
        V_Jeff = 0.12*S_Jeff/(120 + S_Jeff) + 0.12*S_Jeff/(7.95 + S_Jeff) 
    else:
        raise ValueError('Invalid uptake_kinetics_type value!')

    if len([(cef,ms) for cef in capture_eff_vec for ms in ms_vec if ms_vec.index(ms) == capture_eff_vec.index(cef) and ms <= V_Jeff]) > 0:
        (min_cef_Jeff,ms_Jeff) = min([(cef,ms) for cef in capture_eff_vec for ms in ms_vec if ms_vec.index(ms) == capture_eff_vec.index(cef) and ms <= V_Jeff],key=lambda x:x[0])
        print '\V_Jeff = {} --> (min capture eff,ms) = {}\n'.format(V_Jeff, (min_cef_Jeff,ms_Jeff))

        linplot.plot2D(x = [0,min_Cef_Jeff], y = [V_Jeff,V_Jeff], line_format = {'width':1,'color':'black'}, sort_data = False, save_current = False)
        #linplot.add_text(text = "conc. in Gore's \nexperiments (" + str(S_Jeff) + ' mM)', text_format = {'x_pos':15, 'y_pos':V_Jeff + 0.5})    
    else:
        print 'V_Jeff = {} --> No capture efficiency is less than V_Jeff (min ms = {}'.format(V_Jeff, min(ms_vec))

    # Create the lines for S = 100 mM 
    S_high = 100
    if uptake_kinetics_type.lower() == 'conditional':
        V_high = 0.12*S_high/(120 + S_high) 
    elif uptake_kinetics_type.lower() == 'combined':
        V_high = 0.12*S_high/(120 + S_high) + 0.12*S_high/(7.95 + S_high) 
    else:
        raise ValueError('Invalid uptake_kinetics_type value!')

    if len([(cef,ms) for cef in capture_eff_vec for ms in ms_vec if ms_vec.index(ms) == capture_eff_vec.index(cef) and ms <= V_high]) > 0:
        (min_cef_high, ms_high) = min([(cef,ms) for cef in capture_eff_vec for ms in ms_vec if ms_vec.index(ms) == capture_eff_vec.index(cef) and ms <= V_high],key=lambda x:x[0])
        print 'V_high = {} --> (min capture eff,ms) = {}\n'.format(V_high,(min_cef_high, ms_high))

        linplot.plot2D(x = [0,min_cef_high], y = [V_high,V_high], line_format = {'width':1,'color':'black'}, sort_data = False, save_current = False)
        #linplot.add_text(text = 'Uptake rate at \nsucrose conc = 100 mM', text_format = {'x_pos':10, 'y_pos':V_Jeff + 0.5, 'fontweight':'regular'})    
    else:
        print 'V_high = {} --> No capture efficiency is less than V_high (min ms = {}'.format(V_high, min(ms_vec))

    # Shade the area showing the allowed capture efficiency
    linplot.shade(x_range = (min_cef_high,max(capture_eff_vec)), y_range = (),  shade_format = {'color':'green', 'transparency':0.5})

    linplot.customize_and_save()


def performDMMM(input_data):
    """
    Performs DMMM
    """
    model_path = input_data['model_path'] 
    t0 = input_data['t0'] 
    tf = input_data['tf'] 
    delta_t = input_data['delta_t']
    iAZ900 = deepcopy(input_data['iAZ900'])
    capture_eff = deepcopy(input_data['capture_eff'])
    coopr_frac_init = deepcopy(input_data['coopr_frac_init'])
    atp_coeff = deepcopy(input_data['atp_coeff']) 
    init_concs = deepcopy(input_data['init_concs'])
    reactor_type = input_data['reactor_type'] 
    serial_dilution_params = input_data['serial_dilution_params'] 
    simulate_his = input_data['simulate_his']
    his_uptake_byCheater = input_data['his_uptake_byCheater']
    stdout_msgs = input_data['stdout_msgs']
    warnings = input_data['warnings'] 
    results_filename = input_data['results_filename'] 
    save_details = input_data['save_details'] 

    # Compute the initial uptake rates (for demonstration purposes only)
    if stdout_msgs:
        sucr_uptake_rate_init = 0.12*init_concs['sucr_mM']/(120 + init_concs['sucr_mM']) + 0.12*init_concs['sucr_mM']/(7.95 + init_concs['sucr_mM'])
        glc_uptake_rate_init = 0.441*init_concs['glc_mM']/(1.5 + init_concs['glc_mM']) + 0.441*init_concs['glc_mM']/(20 + init_concs['glc_mM'])
        fru_uptake_rate_init = 0.441*init_concs['fru_mM']/(6 + init_concs['fru_mM']) + 0.441*init_concs['fru_mM']/(40 + init_concs['fru_mM'])
        if simulate_his:
            his_uptake_rate_init = 0.066*init_concs['his_mM']/(0.017 + init_concs['his_mM']) + 0.0525*init_concs['his_mM']/(0.370 + init_concs['his_mM'])

        print 'Sucrose: conc = ',init_concs['sucr_mM'],'\tuptake rate = ',sucr_uptake_rate_init
        print 'Glucose: conc = ',init_concs['glc_mM'],'\tuptake rate = ',glc_uptake_rate_init
        print 'Fructose: conc = ',init_concs['fru_mM'],'\tuptake rate = ',fru_uptake_rate_init
        if simulate_his:
            print 'histidine: conc = ',init_concs['his_mM'],'\tuptake rate = ',his_uptake_rate_init

        print '\n--- FBA for iAZ900 before adding atp/adp to SUCRe and incorporating the capture efficiency ---'
        print 'SUCRe:\t',iAZ900.reactions_by_id['SUCRe'].get_equation('id')
        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[-2,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[(-sucr_uptake_rate_init - glc_uptake_rate_init - fru_uptake_rate_init)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
            payoff_C_capeff_one = iAZ900.fba_model.solution['objective_value']
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])


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
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[-2,1000]})
            else:
                set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[(-sucr_uptake_rate_init - glc_uptake_rate_init - fru_uptake_rate_init)/5,1000]})
            iAZ900.fba(stdout_msgs = stdout_msgs)
            if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
                payoff_C_capeff_one = iAZ900.fba_model.solution['objective_value']
                print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  ,  EX_fru_e_ = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'])

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
        print '\n--- FBA for iAZ900 after incorporating the capture efficiency --'
        print 'SUCRe:\t',SUCRe.get_equation('id')
        print 'SUCRe:\t',EX_glc_secreted.get_equation('id')
        print 'SUCRe:\t',EX_fru_secreted.get_equation('id')

        if iAZ900.fixed_o2:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[-2,1000]})
        else:
            set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[(-sucr_uptake_rate_init - glc_uptake_rate_init - fru_uptake_rate_init)/5,1000]})
        iAZ900.fba(stdout_msgs = stdout_msgs)
        if iAZ900.fba_model.solution['exit_flag'] == 'globallyOptimal':
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  , EX_glc_secreted = {} ,  EX_fru_e_ = {}  ,  EX_fru_secreted = {}\n'.format(iAZ900.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['SUCRe'],iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],iAZ900.fba_model.solution['opt_rxnFluxes']['GLCt1'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'], iAZ900.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted'])

    #--- Compute ms and biomass yield for sucrose, glucose, fructose and histidine ---
    carbon_srcs = iAZ900.get_compounds({'sucr_e':'id','glc_D_e':'id','fru_e':'id','his_L_e':'id'})
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_o2_e_':[-2,1000]})
    carbon_srcs['sucr_e'].ms_calc(warnings = False)
    carbon_srcs['sucr_e'].biomass_yield_calc(warnings = False)
    if stdout_msgs:
        print 'sucr: ms = {} ,  biomass_yield = {}'.format(carbon_srcs['sucr_e'].ms,carbon_srcs['sucr_e'].biomass_yield)

    carbon_srcs['glc_D_e'].ms_calc(warnings = False)
    carbon_srcs['glc_D_e'].biomass_yield_calc(warnings = False)
    if stdout_msgs:
        print 'glc_D: ms = {} ,  biomass_yield = {}'.format(carbon_srcs['glc_D_e'].ms,carbon_srcs['glc_D_e'].biomass_yield)

    carbon_srcs['fru_e'].ms_calc(warnings = False)
    carbon_srcs['fru_e'].biomass_yield_calc(warnings = False)
    if stdout_msgs:
        print 'fru_e: ms = {} ,  biomass_yield = {}'.format(carbon_srcs['fru_e'].ms,carbon_srcs['fru_e'].biomass_yield)

    carbon_srcs['his_L_e'].ms_calc(warnings = False)
    carbon_srcs['his_L_e'].biomass_yield_calc(warnings = False)
    if stdout_msgs:
        print 'his_L_e: ms = {} ,  biomass_yield = {}'.format(carbon_srcs['his_L_e'].ms,carbon_srcs['his_L_e'].biomass_yield)

    #--- Cheater and cooperator strains ---
    # Cooperator has a defective HIS3 (YOR202W) gene (for histidone production) corresponding to
    # reaction IGPDH in the model
    # Cheater is missing the gene SUC2 (YIL162W) corresponding ot reaction SUCRe
    Cooperator = deepcopy(iAZ900)
    Cooperator.id = 'Cooperator'
    Cooperator.reset_flux_bounds()
    Cooperator.organism.cells_per_ml = {0:init_concs['coopr_cells_per_ml']}

    Cheater = deepcopy(iAZ900)
    Cheater.id = 'Cheater'
    Cheater.reset_flux_bounds()
    Cheater.organism.cells_per_ml = {0:init_concs['cheater_cells_per_ml']}

    #--- Test cheater and cooperator strains ---
    if stdout_msgs:
        print '\n--- FBA for Cooperator ---' 
        print 'Sucrose: conc = ',init_concs['sucr_mM'],'\tuptake rate = ',sucr_uptake_rate_init
        print 'Glucose: conc = ',init_concs['glc_mM'],'\tuptake rate = ',glc_uptake_rate_init
        print 'Fructose: conc = ',init_concs['fru_mM'],'\tuptake rate = ',fru_uptake_rate_init
        if simulate_his:
            print 'histidine: conc = ',init_concs['his_mM'],'\tuptake rate = ',his_uptake_rate_init

        if Cooperator.fixed_o2:
            if simulate_his:
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_his_L_e_':[-his_uptake_rate_init,1000],'EX_o2_e_':[-2,1000],'IGPDH':[0,0]})
            else: 
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[-2,1000]})
        else:
            if simulate_his:
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_his_L_e_':[-his_uptake_rate_init,1000],'EX_o2_e_':[(-sucr_uptake_rate_init - glc_uptake_rate_init - fru_uptake_rate_init)/5,1000], 'IGPDH':[0,0]})
            else: 
                set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-sucr_uptake_rate_init,1000],'EX_glc_e_':[-glc_uptake_rate_init,1000],'EX_fru_e_':[-fru_uptake_rate_init,1000],'EX_o2_e_':[(-sucr_uptake_rate_init - glc_uptake_rate_init - fru_uptake_rate_init)/5,1000]})

        Cooperator.fba(store_opt_fluxes = False, stdout_msgs = stdout_msgs)
        if Cooperator.fba_model.solution['exit_flag'] == 'globallyOptimal':
            print '\nEX_sucr_e_ = {}, SUCRe = {}  ,  EX_glc_e = {} , GLCt1 = {}  , EX_glc_secreted = {} ,  EX_fru_e_ = {}  ,  EX_fru_secreted = {}\n'.format(Cooperator.fba_model.solution['opt_rxnFluxes']['EX_sucr_e_'], Cooperator.fba_model.solution['opt_rxnFluxes']['SUCRe'],Cooperator.fba_model.solution['opt_rxnFluxes']['EX_glc_e_'],Cooperator.fba_model.solution['opt_rxnFluxes']['GLCt1'], Cooperator.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted'], Cooperator.fba_model.solution['opt_rxnFluxes']['EX_fru_e_'], Cooperator.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted'])
            glc_secreted = Cooperator.fba_model.solution['opt_rxnFluxes']['EX_glc_secreted']
            fru_secreted = Cooperator.fba_model.solution['opt_rxnFluxes']['EX_fru_secreted'] 
        else:
            glc_secreted = 0 
            fru_secreted = 0 

        # Test Cheater 
        print '\n--- FBA for Cheater ---' 
        print '\nglc_secreted = {} , fru_secreted = {}\n'.format(glc_secreted,fru_secreted)
        if Cheater.fixed_o2:
            if simulate_his and his_uptake_byCheater:
                set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_secreted,1000],'EX_fru_e_':[-fru_secreted,1000],'EX_his_L_e_':[-his_uptake_rate_init,1000],'EX_o2_e_':[-2,1000], 'SUCRe':[0,0]})
            else:
                set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_secreted,1000],'EX_fru_e_':[-fru_secreted,1000],'EX_o2_e_':[-2,1000], 'SUCRe':[0,0]})
        else:
            if simulate_his and his_uptake_byCheater:
                set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_secreted,1000],'EX_fru_e_':[-fru_secreted,1000],'EX_his_L_e_':[-his_uptake_rate_init,1000],'EX_o2_e_':[(- glc_secreted - fru_secreted)/5,1000],'SUCRe':[0,0]})
            else:
                set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_glc_e_':[-glc_secreted,1000],'EX_fru_e_':[-fru_secreted,1000],'EX_o2_e_':[(-- glc_secreted - fru_secreted)/5,1000],'SUCRe':[0,0]})

        Cheater.fba(store_opt_fluxes = False, stdout_msgs = stdout_msgs)

    #--- Create the list of shared compounds ---
    # Shared compounds include sucrose, glucose and fructorse and histidine
    #-- Sucrose --
    reactant_reactions = [Cooperator.get_reactions({'EX_sucr_e_':'id'})] 
    sucr_shared = compound(id = 'sucr_e', name = 'Sucrose', ModelSEED_id = ['cpd00076'], KEGG_id = 'C00089', reactant_reactions = reactant_reactions,reactions = reactant_reactions,concentration = {0:init_concs['sucr_mM']})

    #-- Glucose --
    reactant_reactions = [Cooperator.get_reactions({'EX_glc_e_':'id'}),Cheater.get_reactions({'EX_glc_e_':'id'})]
    product_reactions = [Cooperator.get_reactions({'EX_glc_secreted':'id'})]
    glc_shared = compound(id = 'glc_D_e', name = 'D-Glucose', ModelSEED_id = 'cpd00027', KEGG_id = ['C00031','C00267'], reactant_reactions = reactant_reactions,product_reactions = product_reactions, reactions = reactant_reactions + product_reactions, concentration = {0:init_concs['glc_mM']})

    #-- Fructose --
    reactant_reactions = [Cooperator.get_reactions({'EX_fru_e_':'id'}),Cheater.get_reactions({'EX_fru_e_':'id'})]
    product_reactions = [Cooperator.get_reactions({'EX_fru_secreted':'id'})]
    fru_shared = compound(id = 'fru_e', name = 'D-Fructose', ModelSEED_id = 'cpd00082', KEGG_id = ['C10906','C01496','C00095','C02336'], reactant_reactions = reactant_reactions,product_reactions = product_reactions, reactions = reactant_reactions + product_reactions, concentration = {0:init_concs['fru_mM']})

    #-- Histidine --
    if simulate_his:
        if his_uptake_byCheater:
            #reactant_reactions = [Cooperator.get_reactions({'EX_his_L_e_'}),Cheater.get_reactions({'EX_his_L_e_'})]
            reactant_reactions = [Cooperator.reactions_by_id['EX_his_L_e_'],Cheater.reactions_by_id['EX_his_L_e_']]
        else:
            #reactant_reactions = [Cooperator.get_reactions({'EX_his_L_e_'})]
            reactant_reactions = [Cooperator.reactions_by_id['EX_his_L_e_']]

        his_shared = compound(id = 'his_e', name = 'Histidine', ModelSEED_id = 'cpd00119', KEGG_id = 'C00135', reactant_reactions = reactant_reactions,product_reactions = [], reactions = reactant_reactions, concentration = {0:init_concs['his_mM']})
        shared_cmps = [sucr_shared,glc_shared,fru_shared,his_shared]

    #--- Reset the bounds for Cooperator and Cheater ---
    if Cooperator.fixed_o2:
        if simulate_his:
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_o2_e_':[-2,1000],'IGPDH':[0,0]})
        else: 
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_o2_e_':[-2,1000]})
    else:
        # In this case, the bounds on oxygen uptake are updated at each time inside DMMM_yeast.py
        if simulate_his:
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'IGPDH':[0,0]})
        else: 
            set_specific_bounds(model = Cooperator,file_name = model_path + 'iAZ900_minimal.py')

    if Cheater.fixed_o2:
        set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_o2_e_':[-2,1000], 'SUCRe':[0,0]})
    else:
        # In this case, the bounds on oxygen uptake are updated at each time inside DMMM_yeast.py
        set_specific_bounds(model = Cheater,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'SUCRe':[0,0]})

    # Lag phase time was set according to PMID: 23637571 (see the table in supplementary material)
    DMMM_CooprDefec = DMMM_yeast(community_members = [Cooperator,Cheater],shared_compounds = shared_cmps,time_range = [t0,delta_t,tf],reactor_type = reactor_type, carrying_capacity = {'cells_per_ml':1e20,'compounds_mM':1e15}, serial_dilution_params = serial_dilution_params, lag_phase_time = 3, warnings = warnings, stdout_msgs = stdout_msgs)
    DMMM_CooprDefec.run()

    # The fraction of Cooperator at the end of the experiment
    cooperator_frac = Cooperator.organism.cells_per_ml[tf]/(Cooperator.organism.cells_per_ml[tf] + Cheater.organism.cells_per_ml[tf])

    # A dictionary with the following keys: cooperator_conc = cell conc of cooperators, defectors_conc = cell conc of defectors
    # shared_cpds_conc = A dictionary containing the conc of shared compounds. The keys are shared compounds ids and values are
    # shared compounds conentrations (a dictionary with keys being time point and values concetrations). 
    # cooperator_frac: Fraciton of cooperators at the end of simulaiton
    # shared_cpds_rxns = Flux of exchange reactions for shared compounds. This is a dictionary where keys are a tuple with the 
    # first and second element being the id of the respective community member and that of the exchange reaction, respectively.
    # Write the final fraction of Cooperator into a file
    results = {}
    results['cooperator_conc'] = Cooperator.organism.cells_per_ml    
    results['cheater_conc'] = Cheater.organism.cells_per_ml    
    results['cooperator_frac'] = cooperator_frac
    results['shared_cpds_conc'] = {}
    for shared_cmp in shared_cmps:
        results['shared_cpds_conc'][shared_cmp.id] = shared_cmp.concentration 

    if save_details:
        results['cooperator_mu'] = Cooperator.organism.mu
        results['cheater_mu'] = Cheater.organism.mu

        results['shared_cpds_rxns'] = {}
        for rxn in list(set([r for c in shared_cmps for r in c.reactions])):
            results['shared_cpds_rxns'][(rxn.model.id,rxn.id)] = rxn.flux
 
    if results_filename != '':
        with open(results_filename,'a') as f:
            if simulate_his:
                f.write('results[' + str((('coopr_frac',coopr_frac_init), ('capture_eff',capture_eff), ('atp_coeff',atp_coeff),('his',init_concs['his_mM']),('glc',init_concs['glc_mM']),('fru',init_concs['fru_mM']))) + '] = ' + str(results) + '\n')
            else:
                f.write('results[' + str((('coopr_frac',coopr_frac_init), ('capture_eff',capture_eff), ('atp_coeff',atp_coeff),('glc',init_concs['glc_mM']),('fru',init_concs['fru_mM']))) + '] = ' + str(results) + '\n')

def create_model(stdout_msgs = True):
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

    # Assign uptake kinetics for sucrose, glucose, fructose and histidin
    assign_uptakeKinetics(iAZ900)
 
    if stdout_msgs: 
        print '\n--- FBA for original iAZ900 ---'  
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-2,1000]})
    iAZ900.fba(stdout_msgs = stdout_msgs)

    return iAZ900

def master_func(t0,delta_t,tf, start_pos = None, end_pos = None, capture_efficiency = 0.99, SUCRe_atp = 0.115, serial_dilution_params = {'dilution_factor':1000,'time_between_dilutions':24}, reactor_type = 'batch', initial_concs = {'sucr_mM':10, 'glc_mM':0,'fru_mM':0,'his_mM':10,'total_cells_per_ml':1e5,'cooperator_frac_init':0.5}, simulate_his = True, his_uptake_byCheater = True, fixed_o2 = True, warnings = True, stdout_msgs = True, results_filename = '', save_details = True):
    """
    Performs DMMM

    INPUTS:
    ------
                      t0: Initial time
                 delta_t: Time interavl
                      tf: Final time point
               start_pos: Start position of the array containing all possble cases to consider (see all_cases variable)
                 end_pos: End position of the array containing all possble cases to consider (see all_cases variable)
               SUCRe_atp: The atp requirement for SUCRe reaciton (must be positive)
           initial_concs: A dictionary containing the initial concentrations of cells and shared compounds 
            simulate_his: A parameter showing whether we are intereated in simulating the cooperation cost with histidine (True) 
                          or with atp (False)
    his_uptake_byCheater: Indicates whether to allow histidine uptake by cheater (True) or not (False)
                fixed_o2: A parameter showing whether the oxygen uptake rate is fixed (True) or not (False)
       results_filesname: Name of the file containing the results (with NO ".py" extension
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError('stdout_msgs must be True or False')
    if not isinstance(warnings,bool):
        raise TypeError('warnings must be True or False')
    if not isinstance(save_details,bool):
        raise TypeError('save_details must be True or False')

    if not isinstance(capture_efficiency,float) and not isinstance(capture_efficiency,int) and not isinstance(capture_efficiency,list):
        raise TypeError('capture_efficiency must be either a float or list')
    if not isinstance(SUCRe_atp,float) and not isinstance(SUCRe_atp,int) and not isinstance(SUCRe_atp,list):
        raise TypeError('SUCRe_atp must be either a float or list')
    if not isinstance(initial_concs['his_mM'],float) and not isinstance(initial_concs['his_mM'],int) and not isinstance(initial_concs['his_mM'],list):
        raise TypeError("initial_concs['sucr_mM'] must be either a float or list")
    if not isinstance(initial_concs['sucr_mM'],float) and not isinstance(initial_concs['sucr_mM'],int) and not isinstance(initial_concs['sucr_mM'],list):
        raise TypeError("initial_concs['sucr_mM'] must be either a float or list")
    if not isinstance(initial_concs['glc_mM'],float) and not isinstance(initial_concs['glc_mM'],int) and not isinstance(initial_concs['glc_mM'],list):
        raise TypeError("initial_concs['glc_mM'] must be either a float or list")
    if not isinstance(initial_concs['fru_mM'],float) and not isinstance(initial_concs['fru_mM'],int) and not isinstance(initial_concs['fru_mM'],list):
        raise TypeError("initial_concs['fru_mM'] must be either a float or list")
    if not isinstance(initial_concs['cooperator_frac_init'],float) and not isinstance(initial_concs['cooperator_frac_init'],int) and not isinstance(initial_concs['cooperator_frac_init'],list):
        raise TypeError("initial_concs['cooperator_frac_init'] must be either a float or list")
   
    # If a float is provided for capture efficiency, convert it to a list
    if isinstance(capture_efficiency,float) or isinstance(capture_efficiency,int):
        capture_efficiency = [capture_efficiency]
    if isinstance(SUCRe_atp,float) or isinstance(SUCRe_atp,int):
        SUCRe_atp = [SUCRe_atp]
    if isinstance(initial_concs['glc_mM'],float) or isinstance(initial_concs['glc_mM'],int):
        initial_concs['glc_mM'] = [initial_concs['glc_mM']]
    if isinstance(initial_concs['fru_mM'],float) or isinstance(initial_concs['fru_mM'],int):
        initial_concs['fru_mM'] = [initial_concs['fru_mM']]
    if isinstance(initial_concs['his_mM'],float) or isinstance(initial_concs['his_mM'],int):
        initial_concs['his_mM'] = [initial_concs['his_mM']]
    if isinstance(initial_concs['cooperator_frac_init'],float) or isinstance(initial_concs['cooperator_frac_init'],int):
        initial_concs['cooperator_frac_init'] = [initial_concs['cooperator_frac_init']]

    # Check that all elements of SUCRe_atp are non-negative
    if len([a for a in SUCRe_atp if a < 0]) > 0:
        raise ValueError('All elements of SUCRe_atp must be non-negative')

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]

    # Create a compartment for the shared compounds
    shared_camprt = compartment(id = 'shared', name = 'Shared compounds pool')

    #-- Convert the concentration of sucrose from % solution to mM --
    # See: http://abacus.bates.edu/~ganderso/biology/resources/dilutions.html
    # To convert X% solution of a compound with molecular weight MW:
    # molarity = 10*X g/(100 ml) * (1 mol/ MW g)
    # Example: Sucrose concentration = 5% = 10 *5 g/(100 ml) =  0.5 g/ml  
    # 0.5 g/ml * 1 mol/342.2965 g = 0.00146 mol/ml = 1.46 mmol/L = 1.46 mM
    if initial_concs['sucr_percent'] != None and initial_concs['sucr_mM'] == None:
        initial_concs['sucr_mM'] = (10*initial_concs['sucr_percent']/100)*(1/342.2965)*1000
    elif initial_concs['sucr_percent'] != None and initial_concs['sucr_mM'] != None:
        raise ValueError('Values are provided for both "sucr_percent" and "sucr_mM". Set the one that should not be used to "None" and rerun the script')

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/'

    # Create the metabolic model object
    iAZ900 = create_model(stdout_msgs = False)

    iAZ900.fixed_o2 = fixed_o2

    # All possible cases to consider
    all_cases = [(coopr_frac_init,capture_eff,atp_coeff,his_init_conc,glc_init_conc,fru_init_conc) for coopr_frac_init in initial_concs['cooperator_frac_init'] for atp_coeff in SUCRe_atp for capture_eff in capture_efficiency for his_init_conc in initial_concs['his_mM'] for glc_init_conc in initial_concs['glc_mM'] for fru_init_conc in initial_concs['fru_mM']]
    print '\nThe total # of cases to consider = {}'.format(len(all_cases))
    print 'Simulating slice {}\n'.format((start_pos,end_pos))

    if results_filename != '':
        with open(results_filename,'w') as f:
            f.write('results = {}\n')

    if start_pos != None and end_pos != None:
        cases_to_consider = all_cases[start_pos - 1:end_pos]
        counter = start_pos - 1
    else:
        cases_to_consider = all_cases
        counter = 0 


    # Check the required time for one run of the loop
    if len(cases_to_consider) == 1:
        # Total processing and wall time required to create the pyomo model, solve it and store the results
        start_pt = time.clock()
        start_wt = time.time()

    for (coopr_frac_init,capture_eff,atp_coeff,his_init_conc,glc_init_conc,fru_init_conc) in cases_to_consider: 
                       
        counter += 1
        print '{}. coopr_frac = {}, capture_eff = {}, atp_coeff = {}, (his_conc, glc_conc, fru_conc) = ({}, {}, {})'.format(counter, coopr_frac_init, capture_eff, atp_coeff, his_init_conc, glc_init_conc, fru_init_conc)
        
        init_concs = {'sucr_mM':initial_concs['sucr_mM'], 'his_mM':his_init_conc, 'glc_mM':glc_init_conc, 'fru_mM':fru_init_conc, 'coopr_cells_per_ml':initial_concs['total_cells_per_ml']*coopr_frac_init, 'cheater_cells_per_ml':initial_concs['total_cells_per_ml']*(1 - coopr_frac_init)}
    
        #--- Creating a shared memory using the manager ---
        input_data = {}
        input_data['model_path'] = model_path
        input_data['t0'] = t0
        input_data['tf'] = tf
        input_data['delta_t'] = delta_t
        input_data['iAZ900'] = iAZ900
        input_data['coopr_frac_init'] = coopr_frac_init
        input_data['capture_eff'] = capture_eff
        input_data['atp_coeff'] = atp_coeff
        input_data['init_concs'] = init_concs
        input_data['reactor_type'] = reactor_type 
        input_data['serial_dilution_params'] = serial_dilution_params 
        input_data['simulate_his'] = simulate_his
        input_data['his_uptake_byCheater'] = his_uptake_byCheater
        input_data['stdout_msgs'] = stdout_msgs
        input_data['warnings'] = warnings
        input_data['results_filename'] = results_filename 
        input_data['save_details'] = save_details 
                
        p = Process(target = performDMMM, args = (input_data,))
        p.start()
        p.join()
        if p.exitcode > 0:
            raise RuntimeError('Error in python subprocess. Please check performDMMM\n')


    if len(cases_to_consider) == 1:
        # Time required to perform FBA
        elapsed_pt = (time.clock() - start_pt)
        elapsed_wt = (time.time() - start_wt)
        print '\nTotal required time for one run of the loop: Processor time = {} sec/{} min  ,  Walltime = {} sec/{} min\n'.format(elapsed_pt,elapsed_pt/60,elapsed_wt,elapsed_wt/60)

def run_master_func(start_pos = None, end_pos = None, simulate_his = True, fixed_o2 = False, results_filename = ''):
    """
    This funciton runs function master_func for various capture efficiencies and histidine uptake rates to identify
    their impact on the type of game

    INPUTS:
    -------
    simulate_his: A parameter showing whether we are intereated in simulating the cooperation cost with histidine (True) 
                  or with atp (False)
        fixed_02: A parameter showing whether the oxygen uptake rate is fixed (True) or not (False)
    """
    # A [arameter showing whether the details of analysis (such time course concentration of
    # shared metabolites must be saved
    save_details = False

    #--- Perform dynamic simulations --- 
    capture_efficiency = [i/50 for i in range(51)] 

    if simulate_his:
        histidine_init_conc = [i/10 for i in range(14)] # Saturating concentration is 1.2. We go up to 1.3

        # Cooperation cost
        SUCRe_atp = 0.115
    else:
        histidine_init_conc = 0 

        # Cooperation cost
        SUCRe_atp = [i/10 for i in range(101)]

    # Percent concentration of cursoe
    sucr_percent = None
    sucr_mM = 100

    # Total initial cell conc per ml = 150000 (tests were prformed in 5 ml batch culture)
    total_cell_conc_init = 150000

    # Initial fraction of cooperator cells
    coopr_frac = 0.5

    # Number of cycles of serial dilution
    cycles_num = 10

    # Initial concentration of cooperator and cheaters
    coopr_conc_init = coopr_frac*total_cell_conc_init
    cheater_conc_init = (1 - coopr_frac)*total_cell_conc_init

    master_func(t0 = 0, delta_t = 0.25, tf = cycles_num*24, start_pos = start_pos, end_pos = end_pos, capture_efficiency = capture_efficiency, SUCRe_atp = SUCRe_atp, reactor_type = 'serial_dilution', serial_dilution_params = {'dilution_factor':None,'time_between_dilutions':24,'total_cell_conc_beginCycle':total_cell_conc_init}, initial_concs = {'sucr_percent':sucr_percent,'sucr_mM':sucr_mM,'his_mM':histidine_init_conc, 'glc_mM':0,'fru_mM':0,'total_cells_per_ml':total_cell_conc_init, 'cooperator_frac_init':coopr_frac}, simulate_his = simulate_his, fixed_o2 = fixed_o2, warnings = True, stdout_msgs = False, results_filename = results_filename, save_details = save_details)

def plot_dynamic_game_results(results_file_names):
    """
    Plots the results of run_master_func()

    results_file_names: A list of strings containing the names of the results files
    """
    import os
    from imp import load_source
    import numpy as np

    all_results = {}

    for file_name in results_file_names: 
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
    
            results = dataFile.results
            all_results_keys = all_results.keys()
            for k in results.keys():
                if k in all_results_keys:
                    raise userError(str(k) + ' in ' + file_name + ' already in all_results_keys\n')
                else:
                    all_results[k] = results[k]

    if len(all_results.keys()) != 294:
        print '\n**WARNING! Number of integrated results is not equal to 294\n'

    # coopr_frac_init
    coopr_frac_init = list(set([dict(k)['coopr_frac'] for k in all_results.keys()]))
    if len(coopr_frac_init) > 1:
        raise userError('More than one coopr_frac_init found: {} '.format(coopr_frac_init))
    elif len(coopr_frac_init) == 0:
        raise ValueError('No coopr_frac_init found!')
    elif len(coopr_frac_init) == 1:
        coopr_frac_init = coopr_frac_init[0]
    
    # atp_coeff
    atp_coeff = list(set([dict(k)['atp_coeff'] for k in all_results.keys()]))
    if len(atp_coeff) > 1:
        raise userError('More than one atp_coeff found: {} '.format(atp_coeff))
    elif len(atp_coeff) == 0:
        raise ValueError('No atp_coeff found!')
    elif len(atp_coeff) == 1:
        atp_coeff = atp_coeff[0]

    # glc_conc
    glc_conc = list(set([dict(k)['glc'] for k in all_results.keys()]))
    if len(glc_conc) > 1:
        raise userError('More than one glc_conc found: {} '.format(glc_conc))
    elif len(glc_conc) == 0:
        raise ValueError('No glc found!')
    elif len(glc_conc) == 1:
        glc_conc = glc_conc[0]

    # fru_conc
    fru_conc = list(set([dict(k)['fru'] for k in all_results.keys()]))
    if len(fru_conc) > 1:
        raise userError('More than one fru_conc found: {} '.format(fru_conc))
    elif len(fru_conc) == 0:
        raise ValueError('No fru found!')
    elif len(fru_conc) == 1:
        fru_conc = fru_conc[0]


    # Array of capture efficiencies and histidine uptake rate
    capture_efficiencies = sorted(set([dict(k)['capture_eff'] for k in all_results.keys()]))
    his_uptake_rates = sorted(set([dict(k)['his'] for k in all_results.keys()]))

    data = np.zeros((len(his_uptake_rates),len(capture_efficiencies)))

    for i, capture_eff in enumerate(capture_efficiencies):
        for j, his_uptake in enumerate(his_uptake_rates):
            data[i,j] = all_results[(('coopr_frac',coopr_frac_init), ('capture_eff',capture_eff), ('atp_coeff',atp_coeff),('his',his_uptake),('glc',glc_conc),('fru',fru_conc))]['cooperator_frac']

    plot_heatmap(x = np.array(his_uptake_rates), y = 100*np.array(capture_efficiencies), data = data, title = '', xaxis_label = 'Histidine uptake rate (mmol/gDW.h)', yaxis_label = 'Capture efficiency (%)',plot_func = 'pcolor', color_map = None, colormap_labels = None,set_minor_xticks = True, x_minorticks_spacing = 0.1, x_majorticks_spacing = None, set_minor_yticks = True, y_minorticks_spacing = 10, y_majorticks_spacing = None,invert_xaxis = True, invert_yaxis = False, interpolate = True, grid = False, figsize = [25,15], dpi = None, output_filename = 'results/dynamic_games_his.pdf')


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

def create_job_files(total_comb_num,interval_size, job_filename_base, joboutput_filename_base, results_filename_base, simulate_his, fixed_o2, max_walltime):
    """
    Creates the job files given:
             total_comb_num: Total number of combinations 
              interval_size: The desired of iteration intervals
          outfile_base_name: The base name of the output files
          job_filename_base: Base name of the output job file
    joboutput_filename_base: The name of the file storing the dynamic output of the job. It is essentially the same as job_filename_base
                             with the file path deleted (in order to create .out files in the same directory as the job file)
      results_filename_base: Base name of the results file for the main python script
    """
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
        results_filename = results_filename_base + "_" + str(slice[0]) + "_" + str(slice[1]) + '.py'
        print 'creaitng file ',file_name,'...'
        with open(job_filename,'w') as outfile:
            outfile.write('#!/bin/bash -l\n')
            outfile.write('#$-l h_rt=' + str(int(max_walltime))+ ':00:00\n\n')
            outfile.write('cd /usr2/postdoc/alizom/work/yeast_sucrose\n\n')

            # This is to merge the .e[jobid] and .o[jobid] files
            outfile.write('#$ -j y\n')  

            # Set the output file
            outfile.write('\n#$ -o ' + joboutput_filename + '\n') 

            outfile.write('\nsource activate /projectnb/bioinfor/SEGRE/alizom\n\n')

            # python -c "import time;print '\n**Job started at ',time.strftime('%c'),'\n'" 
            outfile.write("\npython -c \"import time;print '\\n**Job started at ',time.strftime('%c'),'\\n'\"\n") 

            # python -c "from dynamicAnalysis import run_master_func;run_master_func(start_pos = start_pos, end_pos = end_pos, siimulate_his = simulate_his, fixed_o2 = fixed_o2, results_filename = results_filename)" 
            outfile.write('\npython -c "from dynamicAnalysis import run_master_func;run_master_func(start_pos = ' + str(slice[0]) + ', end_pos = ' + str(slice[1]) + ', simulate_his = ' + str(simulate_his) + ', fixed_o2 = ' + str(fixed_o2) + ", results_filename = '" + results_filename + "')\"\n") 

            # python -c "import time;print '\n**Job ended at ',time.strftime('%c'),'\n'" >> job_dynamic_yeast_start_pos_end_pos.out 2>&1
            outfile.write("\npython -c \"import time;print '\\n**Job ended at ',time.strftime('%c'),'\\n'\"\n")

def run_create_job_files():
    """
    runs funciton create_job_files to create job files for different cases
    """
    capture_efficiency = [i/20 for i in range(21)] 
    histidine_init_conc = [i/10 for i in range(14)] # Saturating concentration is 1.2. We go up to 1.3
    SUCRe_atp = [i/2 for i in range(51)] 

    # Total combinations for simulating histidine
    total_comb_his = len(capture_efficiency)*len(histidine_uptake_rate)
    total_comb_atp = len(capture_efficiency)*len(SUCRe_atp)
    interval_size_his = 100
    interval_size_atp = 100
    print '\ntotal_comb_his = {} , interval_size_his = {}     total_comb_atp = {} , interval_size_atp = {}\n'.format(total_comb_his,interval_size_his,total_comb_atp,interval_size_atp)
    
    # his, FixedO2
    print '\n--- Creating job files for his, fixedo2 ---'
    create_job_files(total_comb_num = total_comb_his, interval_size = interval_size_his, job_filename_base = 'jobs/job_dyn_his_fixedo2', joboutput_filename_base = 'job_dyn_his_fixedo2', results_filename_base = 'results/game_results_dyn_his_fixedo2', simulate_his = True, fixed_o2 = True, max_walltime = 72)
    
    # his, variable o2
    print '\n--- Creating job files for his, variable O2 ---'
    create_job_files(total_comb_num = total_comb_his, interval_size = interval_size_his, job_filename_base = 'jobs/job_dyn_his_variableO2', joboutput_filename_base = 'job_dyn_his_variableO2', results_filename_base = 'results/game_results_dyn_his_variableO2', simulate_his = True, fixed_o2 = False, max_walltime = 72)
    
    # atp, FixedO2
    print '\n--- Creating job files for atp, fixedo2 ---'
    create_job_files(total_comb_num = total_comb_atp, interval_size = interval_size_atp, job_filename_base = 'jobs/job_dyn_atp_fixedo2', joboutput_filename_base = 'job_dyn_atp_fixedo2', results_filename_base = 'results/game_results_dyn_atp_fixedo2', simulate_his = False, fixed_o2 = True, max_walltime = 72)
    
    # atp, variable O2
    print '\n--- Creating job files for atp, variableO2 ---'
    create_job_files(total_comb_num = total_comb_atp, interval_size = interval_size_atp, job_filename_base = 'jobs/job_dyn_atp_variableO2', joboutput_filename_base = 'job_dyn_atp_variableO2', results_filename_base = 'results/game_results_dyn_atp_variableO2', simulate_his = False, fixed_o2 = False, max_walltime = 72)
    

def coopr_density_effect_on_mu():
    """
    This function assesses the observation that if sugar loss by diffusion is significant 
    then we might expect high-density pure cultures of cooperators to grow faster than low-density cultures, 
    because cells at high density benefit from their hydrolysis products and those of their 
    abundant neighbours.
    """
    #--- Perform dynamic simulations --- 
    capture_efficiency = [0.25]
    #histidine_init_conc = [1.5]
    histidine_init_conc = [5]

    # Cooperation cost
    SUCRe_atp = 0.115

    # Percent concentration of sucrose
    sucr_percent = None
    sucr_mM = 100

    # Number of cycles of serial dilution
    cycles_num = 1

    # Initial fraction of cooperator cells
    coopr_frac = 1

    print "\nsimulation parameters:  capture_efficiency = {}%, initial sucrose conc. = {} mM, initial histidine conc. = {} mM  ,  Cooperators' fraction = {}\n".format(100*capture_efficiency[0], sucr_mM, histidine_init_conc[0], coopr_frac)

    print '\n---- Simulating low cell density ----\n'
    # Total initial cell conc per ml = 150000 (tests were prformed in 5 ml batch culture). 
    for total_cell_conc_init in [0.5*150000, 150000, 2*150000]:
        if total_cell_conc_init == 0.5*150000:
            results_filename = 'results/dynamic_coopr_density_effect_on_mu_low.py'
        elif total_cell_conc_init == 150000:
            results_filename = 'results/dynamic_coopr_density_effect_on_mu_base.py'
        elif total_cell_conc_init == 2*150000:
            results_filename = 'results/dynamic_coopr_density_effect_on_mu_high.py'

        master_func(t0 = 0, delta_t = 0.25, tf = cycles_num*24, capture_efficiency = capture_efficiency, SUCRe_atp = SUCRe_atp, reactor_type = 'serial_dilution', serial_dilution_params = {'dilution_factor':None,'time_between_dilutions':24,'total_cell_conc_beginCycle':total_cell_conc_init}, initial_concs = {'sucr_percent':sucr_percent,'sucr_mM':sucr_mM, 'his_mM':histidine_init_conc, 'glc_mM':0,'fru_mM':0,'total_cells_per_ml':total_cell_conc_init, 'cooperator_frac_init':coopr_frac}, warnings = True, stdout_msgs = True, results_filename = results_filename)


def plot_density(results_filenames, output_filename):
    """
    Makes a plot of the cell densities for two low and high starting cell densities

    INPUTS:
        output_file_name: Name of the output figure
    """
    import numpy as np
    from imp import load_source
    from tools.ancillary.plot import plot, axis
 
    linplot = plot(xaxis = axis(label = 'Time ($h$}', spines_format = {'bottom':{'linewidth':2},'top':{'linewidth':2}}), yaxis = axis(label = '$\mu \, (h^{-1})$', spines_format = {'left':{'linewidth':2}, 'right':{'linewidth':2}}), plot_gridlines = True, show_legend = True, legend_format = {'location':'center left', 'bbox_to_anchor':(1.05,0.5)}, fig_format = {'mathtext_fontname':'Arial'}, output_filename = output_filename)
    
    counter = 0

    for results_filename in results_filenames:
        counter += 1

        # Load the data from the files
        load_source('dataFile',results_filename)
        import dataFile
        results = dataFile.results
    
        if counter == 1:
            # Time points
            time_points = sorted(results[results.keys()[0]]['cooperator_mu'].keys())
            tf = max(time_points)
 
            linplot.xaxis.custom_ticks = range(0,int(tf)+1,24)
            linplot.xaxis.limits = (0,tf)
            create_new_figure = True
        else:
            create_new_figure = False

        current_label = str('{:.1E}'.format(int(results[results.keys()[0]]['cooperator_conc'][0])))
        linplot.plot2D(x = time_points, y = [results[results.keys()[0]]['cooperator_mu'][t] for t in time_points], sort_data = False, label = current_label, line_format = {'width':3}, create_new_figure = create_new_figure, save_current = False)
    
    linplot.customize_and_save()

def coopr_cheater_invade():
    """
    This function checks the obervation that at saturated histidine concentration a cheater
    can invade a cooperator and a cooperator can invade a cheater as well.  
    Under the hsitdine saturation condition, the equilibrium fraction of cooperators depend on the 
    iniital fraction 
    """
    #--- Perform dynamic simulations --- 
    capture_efficiency = [0.25]
    histidine_init_conc = [0.01,1.2]

    # Cooperation cost
    SUCRe_atp = 0.115

    # Percent concentration of sucrose
    sucr_percent = None
    sucr_mM = 100

    # Number of cycles of serial dilution
    cycles_num = 10

    # Total initial cell conc per ml = 150000 (tests were prformed in 5 ml batch culture)
    total_cell_conc_init = 150000

    # Initial fraction of cooperator cells
    coopr_frac = [0.99,0.01]

    print '\n---- Simulating low cell density ----\n'
    master_func(t0 = 0, delta_t = 0.25, tf = cycles_num*24, capture_efficiency = capture_efficiency, SUCRe_atp = SUCRe_atp, reactor_type = 'serial_dilution', serial_dilution_params = {'dilution_factor':None,'time_between_dilutions':24,'total_cell_conc_beginCycle':total_cell_conc_init}, initial_concs = {'sucr_percent':sucr_percent,'sucr_mM':sucr_mM, 'his_mM':histidine_init_conc, 'glc_mM':0,'fru_mM':0,'total_cells_per_ml':total_cell_conc_init, 'cooperator_frac_init':coopr_frac}, warnings = True, stdout_msgs = True, results_filename = 'results/dynamic_coopr_density_effect_on_mu_low')

    # Total initial cell conc per ml = 150000 (tests were prformed in 5 ml batch culture)
    total_cell_conc_init = 10*150000
    
    # Initial concentration of cooperator and cheaters
    coopr_conc_init = coopr_frac*total_cell_conc_init
    cheater_conc_init = (1 - coopr_frac)*total_cell_conc_init

    print '\n---- Simulating high cell density ----\n'

def test_dynamics():
    #--- Perform dynamic simulations --- 
    capture_efficiency = [0.25]
    histidine_init_conc = [1.5]

    # Cooperation cost
    SUCRe_atp = 0.115

    # Percent concentration of cursoe
    sucr_percent = None
    sucr_mM = 100

    # Number of cycles of serial dilution
    cycles_num = 1

    # Total initial cell conc per ml = 150000 (tests were prformed in 5 ml batch culture)
    total_cell_conc_init = 150000

    # Initial fraction of cooperator cells
    coopr_frac = 1 

    master_func(t0 = 0, delta_t = 0.5, tf = cycles_num*24, start_pos = 1, end_pos = 1, capture_efficiency = capture_efficiency, SUCRe_atp = SUCRe_atp, reactor_type = 'serial_dilution', serial_dilution_params = {'dilution_factor':None,'time_between_dilutions':24,'total_cell_conc_beginCycle':150000}, initial_concs = {'sucr_percent':sucr_percent,'sucr_mM':sucr_mM, 'his_mM':histidine_init_conc, 'glc_mM':0,'fru_mM':0,'total_cells_per_ml':total_cell_conc_init, 'cooperator_frac_init':coopr_frac}, warnings = True, stdout_msgs = True, results_filename = 'results/dynamicAnalysis_results_test')


if __name__ == '__main__':
    # Test whether master_func works OK
    test_dynamics()

    # Find the minimum capture efficiency for cooperator to grow on sucrose
    #find_min_sucrose_uptake()

    # Assess the impact of cell density on the growth rates 
    #coopr_density_effect_on_mu()
    #plot_results(results_filenames = 'results/dynamic_coopr_density_effect_on_mu_low.py',input_file_name_high = 'results/dynamic_coopr_density_effect_on_mu_high.py' ,output_file_name = 'results/coopr_density_effect_on_mu.pdf')

