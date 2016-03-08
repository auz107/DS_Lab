from __future__ import division
import sys,os, time
import numpy as np
sys.path.append('../')
sys.path.append('results/')
from copy import deepcopy
import cPickle as pk
import shelve
import itertools
from tools.userError import *
from tools.io.read_gams_model import read_gams_model
from tools.io.read_sbml_model import read_sbml_model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.organism import organism
from tools.core.model import model
from tools.fba.fba import fba
from tools.fba.DMMM import DMMM
from tools.fba.set_specific_bounds import set_specific_bounds
from read_exp_data import read_exp_data
import re
# Increse the recursion limit, otherwise deepcopy will complain
sys.setrecursionlimit(10000)

def  calibrate(t0,delta_t,tf):
    """
    This function compares the growth for wild-type based on experimental data
    with that of simulations
    """
    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]    
    print '\ntime_points = ',time_points

    #--- E. coli iAF1260 model ---
    print '\n--- Wild-type E.coli (iAF1260 model) ----'
    # Define the organism
    model_organism = organism(id = 'Ecoli', name = 'Escherichia coli',domain = 'Bacteria', genus = 'Escherichia', species = 'coli', strain = 'MG1655',gDW_per_cell = 2.8e-13)

    WT = read_gams_model(file_name = '../models/Escherichia_coli/iAF1260/iAF1260ModelData.py',model_name = 'iAF1260',model_organism = model_organism,model_type = 'metabolic')

    WT.biomass_reaction = WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'})
    WT.all_biomass_reactions = {'core':WT.get_reactions({'Ec_biomass_iAF1260_core_59p81M':'id'}),'WT':WT.get_reactions({'Ec_biomass_iAF1260_WT_59p81M':'id'})}
    WT.organism.cells_per_mL = {0:7.5e6}

    # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
    # Example: EX_glc(E): glc-D[e] <==>    Vmax*C['glc-D[e]']/(Km + C['glc-D[e]']) 
    # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
    for reaction in [r for r in WT.reactions if r.type.lower() == 'exchange']:
        # The id of compound participating in the exchange reaction
        metab_id = [m.id for m in reaction.compounds][0]
        reaction.kinetics = "10*C['" + metab_id + "']/(10 + C['" + metab_id + "'])"

    # Glucose uptake kinetics 
    exch_rxns = WT.get_reactions({'EX_glc(e)':'id','EX_lys-L(e)':'id','EX_ile-L(e)':'id'})
    exch_rxns['EX_glc(e)'].kinetics = "10*C['glc-D[e]']/(0.15 + C['glc-D[e]'])"
    exch_rxns['EX_lys-L(e)'].kinetics = "0.1964*C['lys-L[e]']/(5e-4 + C['lys-L[e]']) + 0.3055*C['lys-L[e]']/(1e-2 + C['lys-L[e]'])"
    exch_rxns['EX_ile-L(e)'].kinetics = "0.0346*C['ile-L[e]']/(1.22e-3 + C['ile-L[e]'])"
   
    # Growth medium
    set_specific_bounds(WT,file_name = '../models/Escherichia_coli/iAF1260/iAF1260_minimal_glucose_anaerobic.py',simulation_condition = 'minimal_glucose_anaerobic')

    # Assign the objective function coefficients
    for rxn in WT.reactions:
        rxn.objective_coefficient = 0
    WT.biomass_reaction.objective_coefficient = 1

    # Perform FBA for the wild-type
    print '\n------- Perform FBA for the wild-type ------'
    WT.fba(assign_wildType_max_biomass = True)

    # Compute ms and biomass yield for glucose
    glc_D = WT.get_compounds({'glc-D[e]':'id'})
    glc_D.ms_calc() 
    glc_D.biomass_yield_calc() 

    #-------- pfs mutant (positive control in Pam Silver's experiments ------
    pfs = deepcopy(WT)
    pfs.id = WT.id + '_pfs'
    pfs.organism.id = 'pfs_Ecoli'
    pfs.organism.cells_per_mL = {0:7.5e6}
    pfs.organism.gDW_per_mL = {0:7.5e6*pfs.organism.gDW_per_cell}
    # Set the LB and UB corresponding to reacitons for pfs deletion to zero
    for rxn_id in ['5DOAN','AHCYSNS','MTAN']:
        rxn = pfs.get_reactions({rxn_id:'id'})
        rxn.flux_bounds = [0,0]
    print '\n----- Performing fba for the positive control (pfs) -----'
    pfs.fba(create_model = False, store_opt_fluxes = False)

    #---- Perform dynamic simulations for WT because pfs is not viable in silico ----
    # Create the list of shared compounds
    shared_cmps = []

    # Glucose as a shared compound (concentration in mM)
    glucose = compound(id = 'glc-D[e]', name = 'D-Glucose', Kegg_id = 'C00031', reactant_reactions = [r for r in WT.reactions if r.id == 'EX_glc(e)'],reactions = [r for r in WT.reactions if r.id == 'EX_glc(e)'],concentration = {0:111.01})

    shared_cmps.append(glucose)

    DMMM_WT = DMMM(community_members = [WT],shared_compounds = shared_cmps, time_range = [t0,delta_t,tf],store_dynamic_fluxes = False, screen_output = 'on')
    DMMM_WT.run()

    # Load the experimental growth data
    [day1Ave,day2Ave,day3Ave,day4Ave] = read_exp_data()

    # Write the results into a MATLAB file to plot
    with open('results/calibrate_results_toPlot_withMatlab.m','w') as f:
        # Write simulation results
        f.write('model = [\n')
        for t in time_points:
            f.write(str(t) + '\t' + str(WT.organism.cells_per_mL[t]) + '\n')
        f.write('];\n')

        # Write experimental data
        f.write('\nexperiment = [\n')
        f.write('0\t7.5e6\n')
        f.write('24\t' + str(day1Ave[('Positive','Positive')]/2)+'\n')
        f.write('48\t' + str(day2Ave[('Positive','Positive')]/2)+'\n')
        f.write('72\t' + str(day3Ave[('Positive','Positive')]/2)+'\n')
        f.write('96\t' + str(day4Ave[('Positive','Positive')]/2)+'\n')
        f.write('];\n')

#---------------------
if __name__ == '__main__':
    calibrate(t0 = 0,delta_t = 0.5,tf = 96)

