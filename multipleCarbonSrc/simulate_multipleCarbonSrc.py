from __future__ import division
import sys, os
sys.path.append('../')
from tools.io.read_sbml_model import read_sbml_model
from tools.core.model import model
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.compartment import compartment
from tools.core.organism import organism
from tools.userError import userError
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from DMMM_multipleCarbonSrc import DMMM_multipleCarbonSrc
from tools.ancillary.get_ModelSeed_ids import get_ModelSeed_ids
from models.model_seed_database.ModelSeed_compounds import ModelSeed_compounds
from tools.ancillary.remove_non_alphanumeric import remove_non_alphanumeric
import re, math, itertools, cobra


def addATPMrxn(model):
    """
    Adds ATPM reaction to a metabolic model
    """
    #- ATPM   [c] : atp + h2o --> adp + h + pi - 
    # atp: cpd00002   h2o: cpd00001   adp: cpd00008   h: cpd00067    pi: cpd00012
    cpd00002 = [cmp for cmp in model.compounds if cmp.id.lower() in ['cpd00002_c0','atp_c','atp[c]']]
    if len(cpd00002) == 1:
       cpd00002 = cpd00002[0]
    elif len(cpd00002) > 1:
        raise userError('More than one compounds found for atp')
    else:
        raise userError('No compound for atp')
    cpd00001 = [cmp for cmp in model.compounds if cmp.id.lower() in ['cpd00001_c0','h2o_c','h2o[c]']]
    if len(cpd00001) == 1:
       cpd00001 = cpd00001[0]
    elif len(cpd00001) > 1:
        raise userError('More than one compounds found for h2o')
    else:
        raise userError('No compound for h2o')
    cpd00008 = [cmp for cmp in model.compounds if cmp.id.lower() in ['cpd00008_c0','adp_c','adp[c]']]
    if len(cpd00008) == 1:
       cpd00008 = cpd00008[0]
    elif len(cpd00008) > 1:
        raise userError('More than one compounds found for adp')
    else:
        raise userError('No compound for adp')
    cpd00067 = [cmp for cmp in model.compounds if cmp.id.lower() in ['cpd00067_c0','h_c','h[c]']]
    if len(cpd00067) == 1:
       cpd00067 = cpd00067[0]
    elif len(cpd00067) > 1:
        raise userError('More than one compounds found for h')
    else:
        raise userError('No compound for h')
    cpd00012 = [cmp for cmp in model.compounds if cmp.id.lower() in ['cpd00012_c0','pi_c','pi[c]']]
    if len(cpd00012) == 1:
        cpd00012 = cpd00012[0]
    elif len(cpd00012) > 1:
        raise userError('More than one compounds found for pi')
    else:
        raise userError('No compound for pi')
    
    # atp: cpd00002   h2o: cpd00001   adp: cpd00008   h: cpd00067    pi: cpd00012
    ATPM = reaction(id = 'ATPM', name = 'Maintenance ATP', stoichiometry = {cpd00002:-1,cpd00001:-1,cpd00008:1,cpd00067:1,cpd00012:1},type = 'irreversible')
    
    model.add_reactions(new_reactions = [ATPM])

def simulate_community(t0,delta_t,tf,members_info,carbon_sources,media_cmp_concs,serial_dilution_params,exchrxns_objcoeff,results_filename_base):
    """
    Simulates the growth of multiple species

    INPUTS:
    ------
              t0: Initial time
         delta_t: Time interavl
              tf: Final time point
    members_info: A dictionary where the keys are the names of community members and values are another
                  dictionary with the following keys: 
                          'model_filename': The model filename (including the path) for that community member 
                                            (some exchange reactions in the minimal medium do not exist in the
                                            models for some species)
                                'model_id': The model id
                  'growth_medium_filename': The nmae of the file (including the path) containing the 
                                            flux bounds for reactionsn in the medium other than the carbon 
                                            source and oxygen (some exchange reactions in the minimal medium
                                            do not exist in the models for some species)
                             biomassrxn_id: Id of the biomass reaction in the model
                         oxygen_exchrxn_id: Id of exchange reaction for exygen in the model
                            'gDW_per_cell': Gram of dry cell weight per cell
                       'cells_per_mL_init': Initial cell concentration (cells/mL)
    media_cmp_concs: A dictionary where keys are the ModelSeed ids of compounds in the growth medium
                     and values are their concentrations in mM, M or g/l 
    carbon_sources: A dictionary with keys being the names of the limiting carbon sources in the medium and 
                    values their ModelSeed ids
    exchrxns_weight: Objective coefficient of exchange reactions in the model (in order to force metabolites
                     to be produced)
    serial_dilution_params: Same as serial_dilution_params in DMMM
    results_filename_base: Base file name to store the results
    """
    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Generate all time points
    time_points = [k/10 for k in range(t0,int(tf*10 + delta_t*10),int(delta_t*10))]

    # Create a compartment for the shared compounds
    shared_camprt = compartment(id = 'shared', name = 'Shared compounds pool')

    # A list of model objects for community members
    models = []

    for org_name in members_info.keys():

        print '\n----- ' + org_name + ' ------'

        # Define the organism
        model_organism = organism(id = org_name, gDW_per_cell = members_info[org_name]['gDW_per_cell'])

        # Read the sbml file
        model = read_sbml_model(file_name = members_info[org_name]['model_filename'], model_id = members_info[org_name]['model_id'],model_organism = model_organism,model_type = 'cmpolic',import_params = False)

        # Check whether the model has ATPM reaction
        if model.get_reactions({'ATPM':'id'}) == None:
            print 'Adding ATPM reaction to the model ...'
            addATPMrxn(model)

        # Get ModelSeed ids for reactions and compounds
        get_ModelSeed_ids(model = model,screen_output = 'off') 

        # Check the model after the modifications 
        model.validate()

        # Assign the coefficient of the biomass reaction
        model.biomass_reaction = model.get_reactions({members_info[org_name]['biomassrxn_id']:'id'})
        for rxn in model.reactions:
            rxn.objective_coefficient = 0
        model.biomass_reaction.objective_coefficient = 1

        # Find out the growth without any carbon source
        print 'FBA with no carbon source ...'
        set_specific_bounds(model = model,flux_bounds = {members_info[org_name]['oxygen_exchrxn_id']:[-2*min([10*media_cmp_concs[carbonSrc]/(10 + media_cmp_concs[carbonSrc])for carbonSrc in carbon_sources.values()]),1000],'ATPM':[1.5,1.5]},file_name = members_info[org_name]['growth_medium_filename'],simulation_condition = 'minimal_aerobic_' + '_'.join(carbon_sources.keys()))
        model.fba() 
        if model.fba_model.solution['exit_flag'] == 'globallyOptimal':
            model.noCarbonSrc_maxBiomass = model.fba_model.solution['objective_value']
        else:
            model.noCarbonSrc_maxBiomass = 0 

        # Set the growth medium. Note that the bound on oxygen and carbon source will change duriing the 
        # dynamic simulaitons 
        print 'FBA with supplied carbon sources ...'
        flux_bounds = dict([(model.get_reactions({'EX_' + carbonSrc + '_e0':'ModelSeed_id'}).id,[-10*media_cmp_concs[carbonSrc]/(10 + media_cmp_concs[carbonSrc]),1000]) for carbonSrc in carbon_sources.values()] + [(members_info[org_name]['oxygen_exchrxn_id'],[-2*min([10*media_cmp_concs[carbonSrc]/(10 + media_cmp_concs[carbonSrc]) for carbonSrc in carbon_sources.values()]),1000]),('ATPM',[1.5,1.5])])
        set_specific_bounds(model = model,flux_bounds = flux_bounds,file_name = members_info[org_name]['growth_medium_filename'],simulation_condition = 'minimal_aerobic_')
        for r in model.reactions:
            r.flux_bounds_init = r.flux_bounds
        model.fba() 

        # Assign a general Michaelis-Menten type uptake kinetics to all exchange reactions
        # Example: EX_glc(E): glc-D[e] <==>    Vmax*C['glc-D[e]']/(Km + C['glc-D[e]']) 
        # Use a Vmax value of 10 mmole/gDW.h and a Km value of 10 micro-M 
        for exch_rxn in [r for r in model.reactions if r.type.lower() == 'exchange']:
            # The id of compound participating in the exchange reaction
            cmp_id = exch_rxn.reactants[0].id
            exch_rxn.kinetics = "10*C['" + cmp_id + "']/(10 + C['" + cmp_id + "'])"

        # Compute ms and biomass yield for limiting carbon sources in the extracelluular environment
        for carbonSrc in carbon_sources.values():
            #carbonSrc_cmp = model.get_compounds({carbonSrc + '_e0':'ModelSeed_id_compart'}) 
            carbonSrc_cmp = model.get_reactions({'EX_' + carbonSrc + '_e0':'ModelSeed_id'}).reactants[0]
            carbonSrc_cmp.ms_calc()
            carbonSrc_cmp.biomass_yield_calc()

        # Set the initial cell concentrations
        model.organism.cells_per_mL = {0:members_info[org_name]['cells_per_mL_init']}
        model.organism.gDW_per_mL = {0:members_info[org_name]['cells_per_mL_init']*model.organism.gDW_per_cell}

        # Exchange reaction for the limiting carbon source and oxygen
        model.limiting_carbonSrc_exchrxns = []
        for carbonSrc in carbon_sources.values(): 
            carbonSrc_exchrxn = [r for r in model.reactions if r.type.lower() == 'exchange' and r.reactants[0].ModelSeed_id == carbonSrc]
            if len(carbonSrc_exchrxn) != 1:
                raise userError(str(len(carbonSrc_exchrxn)) + ' exchange reactions for the limiting carbon source in model ' + model.id)
            carbonSrc_exchrxn = carbonSrc_exchrxn[0]
            carbonSrc_exchrxn.flux_bounds = [-1000,1000]
            model.limiting_carbonSrc_exchrxns.append(carbonSrc_exchrxn)

        oxygen_exchrxn = model.get_reactions({members_info[org_name]['oxygen_exchrxn_id']:'id'}) 
        if isinstance(oxygen_exchrxn,list) and len(oxygen_exchrxn) != 1:
            raise userError(str(len(oxygen_exchrxn)) + ' exchange reactions for oxygen in model ' + model.id)
        oxygen_exchrxn.flux_bounds = [-1000,1000]
        model.oxygen_exchrxn = oxygen_exchrxn

        models.append(model)

    print '\n----- Find the list of shared compounds in the medium ------'
    shared_cmps = []
     
    #-- Add the limiting resource given in media_cmp_concs --
    for cmp in media_cmp_concs.keys():
        reactant_reactions = [exchrxn for exchrxn in [r for model in models for r in model.reactions if r.type.lower() == 'exchange' and (len(r.reactants) > 0 and cmp in [re.sub('_e0','',r.reactants[0].id),r.reactants[0].ModelSeed_id])]]
        if len(reactant_reactions) == 0:
            raise userError('No reactant reactions for limiting carbon source ' + cmp)
        cmp_shared = compound(id = cmp, compartment = shared_camprt, name = ModelSeed_compounds[cmp]['name'], ModelSeed_id = media_cmp_concs[cmp], Kegg_id = ModelSeed_compounds[cmp]['Kegg_id'],formula = ModelSeed_compounds[cmp]['formula'], reactant_reactions = reactant_reactions, concentration = {0:media_cmp_concs[cmp]})
        shared_cmps.append(cmp_shared)

    #-- All other common compounds --
    # We should also allow the exchange of any metabolites that might produced by some community membersi
    # However, we need to exchlude from this list, exchange reactions for metabolites present in excess
    # in the growth medium and those for biomass
    # List of exchange rxns for all other compounds in the minimal medium (e.g., those available in excess)
    media_excess_exchrxns = [r for model in models for r in model.reactions if r.type.lower() == 'exchange' and r.flux_bounds[0] < 0 if r.reactants[0].id not in media_cmp_concs.keys()]

    # Exchange reaction for biomass
    biomass_exchrxns = [r for model in models for r in model.reactions if r.type.lower() == 'exchange' and (r.reactants[0].ModelSeed_id == 'cpd11416' or 'biomass' in r.reactants[0].id.lower() or (r.reactants[0].name != None and 'biomass' in r.reactants[0].name.lower()) or 'cpd11416' in r.reactants[0].id.lower()  or len([syn for syn in r.reactants[0].synonyms if 'biomass' in syn.lower()]) > 0 or (r.reactants[0].name != None and 'cpd11416' in  r.reactants[0].name))] 

    # Reactions that must be excluded from the list of shared metabolites
    excluded_exchrxns = media_excess_exchrxns + biomass_exchrxns

    # Add a filed named "matched" to all exchange reactions that may serve for the uptake/export of any 
    # potential shared compound. This field geta a value of True or False indicating whether that reactions 
    # has been matached with any other reaction in the model or not
    for exchrxn in [r for model in models for r in model.reactions if r.type.lower() == 'exchange' and r not in excluded_exchrxns]: 
        exchrxn.matched = False

    # List of exchange reactions for any other compound that is not in the growth medium initially but 
    # can be potentially shared compounds through inter-species interactions
    for model in models:
        model.alreadyConsidered = False 
    for model in models:
        # Remove a model being considered in this iteration from the list to avoid repetitions
        model.alreadyConsidered = True 
        for exchrxn in [r for r in model.reactions if r.type.lower() == 'exchange' and r not in excluded_exchrxns and not r.matched]: 
            # Check whether this reaction has any match in any other model
            matches = [r for othermodel in models if not othermodel.alreadyConsidered for r in othermodel.reactions if r.type.lower() == 'exchange' and r not in excluded_exchrxns and (r.id == exchrxn.id or (r.ModelSeed_id != None and r.ModelSeed_id == exchrxn.ModelSeed_id) or (r.Kegg_id != None and r.Kegg_id == exchrxn.Kegg_id))]
            # Create a shared compound if at least one match is found
            if len(matches) >= 1:

                exchrxn.matched = True
                exchrxn.objective_coefficient = exchrxns_objcoeff 
                for r in matches:
                    r.matched = True
                    r.objective_coefficient = exchrxns_objcoeff 

                cmp_ModelSeed_id = exchrxn.reactants[0].ModelSeed_id
                shared_cmp = compound(id = cmp_ModelSeed_id, compartment = shared_camprt, name = ModelSeed_compounds[cmp_ModelSeed_id]['name'], Kegg_id = ModelSeed_compounds[cmp_ModelSeed_id]['Kegg_id'],formula = ModelSeed_compounds[cmp_ModelSeed_id]['formula'], reactant_reactions = [exchrxn] + matches, product_reactions = [exchrxn] + matches, concentration = {0:0.00})
                shared_cmps.append(shared_cmp)

    print 'Total # of exchange rxns in all models = ',len([r for model in models for r in model.reactions if r.type.lower() == 'exchange'])
    print 'Total # of shared compounds = ',len(shared_cmps)

    print '\n----- Performing DMMM ---------'
    DMMM_m1m2 = DMMM_multipleCarbonSrc(community_members = models,shared_compounds = shared_cmps, time_range = [t0,delta_t,tf], reactor_type = 'serial_dilution', serial_dilution_params = serial_dilution_params, carrying_capacity = {'cells_per_mL':1e14,'compounds_mM':1e15},store_dynamic_fluxes = False, screen_output = 'on')
    DMMM_m1m2.run()

    print '\n----- Storing the results for ',tuple([model.organism.id for model in models]),'---------'
    storeDMMMresults(models = models,shared_cmps = [c for c in shared_cmps if c.id not in media_cmp_concs.keys()],carbon_sources = carbon_sources, results_filename_base = results_filename_base)

def storeDMMMresults(models,shared_cmps,carbon_sources,results_filename_base): 
    """
    This function plots the DMMM results into a matlab file
    INPUTS:
    -------
         models: A list of models for community members
    shared_cmps: The same shared_cmps in the input of DMMM
    time_points: The set of time points computed in performDMMM
    carbon_sources
    results_filename_base: The same as in performDMMM
    """
    # Initialize the file
    with open(results_filename_base + '_' + '_'.join(carbon_sources.keys()) + '.py','w') as f:
        f.write('results = []\n')
    with open(results_filename_base + '_' + '_'.join(carbon_sources.keys()) + '.m','w') as f:
        f.write('results = struct([]);\n')

    time_points = sorted(shared_cmps[0].concentration.keys())

    print '\nCell concentrations (',[model.id for model in models],'):'
    #for t in time_points:
    #    print t,'\t',[(model.organism.mu[t],model.organism.cells_per_mL[t]) for model in models]
    print '\n'

    # Threshould for metabolite concentrations
    conc_thr = 1e-6
    shared_cmps_concs_stored = {}
    shared_cmps_nonzero_conc = [c for c in shared_cmps if max(c.concentration.values()) > conc_thr]
    for shared_cmp in shared_cmps_nonzero_conc: 
        shared_cmps_concs_stored[shared_cmp.id] = shared_cmp.concentration
        print '\n{}\t{}\t{:6f}'.format(shared_cmp.id,ModelSeed_compounds[shared_cmp.id]['name'][0],max(shared_cmp.concentration.values()))
        #for t in sorted(shared_cmp.concentration.keys()):
        #    print '{}\t{}'.format(t,shared_cmp.concentration[t]) 

    print '\nExchange rxns for shared compounds with non-zero concentration = '
    # Threshould for reaction flux
    rxn_thr = 1e-6
    exch_rxns_stored = {}
    for rxn in sorted(list(set([r for c in shared_cmps_nonzero_conc for r in c.reactions if max([abs(v) for v in r.flux.values()]) > rxn_thr])),key=lambda x:x.model.organism.id):
        exch_rxns_stored[(rxn.model.organism.id,rxn.id)] = rxn.flux
        print '\n{}\t{}\t{}\t({:6f},{:6f})'.format(rxn.model.organism.id,rxn.id,rxn.name,min([v for v in rxn.flux.values()]),max([v for v in rxn.flux.values()]))
        #for t in sorted(r.flux.keys()):
        #    print '{}\t{}'.format(t,r.flux[t])

    # The following outputs are saved: shared_cmps and organism object 
    # for each community member 
    with open(results_filename_base  + '_' + '_'.join(carbon_sources.keys())  + '.py','a') as f:
        f.write('results.append(' + repr({'name':tuple([model.organism.id for model in models]),'cell_concs':dict([(model.organism.id,model.organism.cells_per_mL) for model in models]),'shared_cmp_concs':shared_cmps_concs_stored,'exch_rxns':exch_rxns_stored}) + ')\n')

    # Write the results into a MATLAB file to plot later
    with open(results_filename_base  + '_' + '_'.join(carbon_sources.keys())  + '.m','a') as f:
        f.write('\ncurrent_index = length(results) + 1;\n')
        f.write("results(current_index).organism_names = {'" + "','".join([model.organism.id for model in models]) + "'};\n")
        f.write('results(current_index).time = ' + str(time_points) + ';\n')

        # Cell concentrations
        model_cell_concs = []
        for t in time_points:
            for model in models:
                model_cell_concs.append(model.organism.cells_per_mL[t])
        f.write('results(current_index).cell_concs.' + model.organism.id + ' = ' + str(model_cell_concs) + ';\n')

        # Shared cmp concnntrations
        for shared_cmp_name in shared_cmps_concs_stored.keys():
            shared_cmp_concs = []
            for t in time_points:
                shared_cmp_concs.append(shared_cmps_concs_stored[shared_cmp_name][t])
            f.write('results(current_index).shared_cmp_concs.' + remove_non_alphanumeric(shared_cmp_name) + ' = ' + str(shared_cmp_concs) + ';\n')

        # Exchange reactions for shared cmps 
        for (org_id,exch_rxn_id) in exch_rxns_stored.keys():
            exch_rxn_fluxes = []
            for t in time_points:
                exch_rxn_fluxes.append(exch_rxns_stored[(org_id,exch_rxn_id)][t])
            f.write('results(current_index).exch_rxn_fluxes.' + re.sub('_LParen_e_RParen','',remove_non_alphanumeric(exch_rxn_id)) + '_' + remove_non_alphanumeric(org_id) + ' = ' + str(exch_rxn_fluxes) + ';\n')

def plotResults():
    """
    Plots the DMMM results

    INPUTS:
    -------
    tf: Final time point
    """
    import numpy as np
    from imp import load_source
    import matplotlib
    import matplotlib.pyplot as plt

    res_files = ['Pseudomonas_chlororaphis_Enterobacter_aerogenes_galacturonic_acid_L_serine.py', 'Pseudomonas_chlororaphis_Pseudomonas_fluorescens_galacturonic_acid_L_serine.py', 'Pseudomonas_chlororaphis_Serratia_marcescens_galacturonic_acid_L_serine.py', 'Pseudomonas_fluorescens_Enterobacter_aerogenes_galacturonic_acid_L_serine.py', 'Pseudomonas_putida_Enterobacter_aerogenes_galacturonic_acid_L_serine.py', 'Pseudomonas_putida_Pseudomonas_chlororaphis_galacturonic_acid_L_serine.py', 'Pseudomonas_putida_Pseudomonas_fluorescens_galacturonic_acid_L_serine.py', 'Pseudomonas_putida_Serratia_marcescens_galacturonic_acid_L_serine.py', 'Pseudomonas_veronii_Enterobacter_aerogenes_galacturonic_acid_L_serine.py', 'Pseudomonas_veronii_Pseudomonas_chlororaphis_galacturonic_acid_L_serine.py', 'Pseudomonas_veronii_Pseudomonas_fluorescens_galacturonic_acid_L_serine.py', 'Pseudomonas_veronii_Pseudomonas_putida_galacturonic_acid_L_serine.py', 'Pseudomonas_veronii_Serratia_marcescens_galacturonic_acid_L_serine.py', 'Serratia_marcescens_Enterobacter_aerogenes_galacturonic_acid_L_serine.py', 'Serratia_marcescens_Pseudomonas_fluorescens_galacturonic_acid_L_serine.py']

    #for res_file in res_files:
    for res_file in  ['Pseudomonas_putida_Pseudomonas_fluorescens_galacturonic_acid_L_serine.py']:

        #fig, ax = plt.subplots(figsize = (16,10))
        fig, ax = plt.subplots(figsize = (25,15))
        #ax.set_aspect('auto')

        # Update the matplotlib configuration parameters:
        matplotlib.rcParams.update({'font.size': 30,'font.weight':'bold'})

        results = []

        # Load the data from the files
        load_source('dataFile','results/' + res_file) 
        import dataFile
        results = dataFile.results

        # Names of the organisms
        org_names = results[0]['cell_concs'].keys()

        # Time points
        time_points = [t for t in sorted(results[0]['cell_concs'][org_names[0]].keys())]
        tf = max(time_points)

        # Total cell concentration at each time point
        total_cell_conc = dict([(t,sum([results[0]['cell_concs'][org][t] for org in org_names])) for t in time_points])

        # Fraction of each microorganism over time
        org_fracs = dict([(org,dict([(t,results[0]['cell_concs'][org][t]/total_cell_conc[t]) for t in time_points])) for org in org_names]) 
        # Print the fractions in the final time piont to the output
        print '{}: {:0.2f}  ,  {}: {:0.2f}'.format(org_names[0],org_fracs[org_names[0]][tf],org_names[1],org_fracs[org_names[1]][tf])

        # Labels for organisms
        org_labels = {'Pseudomonas_chlororaphis':'Pc','Pseudomonas_fluorescens':'Pf','Pseudomonas_putida':'Pp','Pseudomonas_veronii':'Pv','Enterobacter_aerogenes':'Ea','Serratia marcescens':'Sm'}

        #--- plot the results ---
        for org in org_names:
            print org
            log10_cell_conc = dict([(t,np.log10(results[0]['cell_concs'][org][t])) for t in time_points])
            ax.plot(time_points, [results[0]['cell_concs'][org][t] for t in time_points], label=org_labels[org], linewidth = 2)
            #ax.plot(time_points, [log10_cell_conc[t] for t in time_points], label=org_labels[org], linewidth = 2)

        ax.set_xlabel('Time (h)',{'weight':'bold','size':30})
        ax.set_ylabel('Cell conc (cells/ml)',{'weight':'bold','size':30})
     
        ax.legend(loc=0)

        ax.set_xlim([0,tf])
        ax.set_xticks(range(0,240+1,24))

        for ticklabel in ax.get_xmajorticklabels():
            ticklabel.set_fontsize(30)
            ticklabel.set_fontweight('bold')
        for ticklabel in ax.get_ymajorticklabels():
            ticklabel.set_fontsize(30)
            ticklabel.set_fontweight('bold')

        ax.grid(True, color = 'k', linewidth = 2)

        fig.savefig('results/plot_'+ org_labels[org_names[0]] + '_' + org_labels[org_names[1]] + '.pdf')

        plt.tight_layout()
        plt.show()

#------------------------

if __name__ == '__main__':

    # Path to the directory where the models exist
    models_path = '/data/alizom/models/'

    all_organisms = {}
    all_organisms['Enterobacter_aerogenes'] = {'model_filename':models_path + 'Enterobacter_aerogenes/kbase/Enterobacter_aerogenes.20.xml','model_id':'Enterobacter_aerogenes.20','growth_medium_filename':models_path + 'Enterobacter_aerogenes/kbase/minimal_medium.py', 'biomassrxn_id':'biomass0','oxygen_exchrxn_id':'EX_cpd00007_e0','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}
    all_organisms['Pseudomonas_chlororaphis'] = {'model_filename':models_path + 'Pseudomonas_chlororaphis/kbase/Pseudomonas_chlororaphis.13.xml','model_id':'Pseudomonas_chlororaphis.13','growth_medium_filename':models_path + 'Pseudomonas_chlororaphis/kbase/minimal_medium.py', 'biomassrxn_id':'biomass0','oxygen_exchrxn_id':'EX_cpd00007_e0','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}
    all_organisms['Pseudomonas_fluorescens'] = {'model_filename':models_path + 'Pseudomonas_fluorescens/kbase/Pseudomonas_fluorescens.11.xml','model_id':'Pseudomonas_fluorescens.11','growth_medium_filename':models_path + 'Pseudomonas_fluorescens/kbase/minimal_medium.py', 'biomassrxn_id':'biomass0','oxygen_exchrxn_id':'EX_cpd00007_e0','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}
    #all_organisms['Pseudomonas_putida'] = {'model_filename':models_path + 'Pseudomonas_putida/iJP962/Pseudomonas_putida_iJP962.5.xml','model_id':'Pseudomonas_putida_iJP962.5','growth_medium_filename':models_path + 'Pseudomonas_putida/iJP962/minimal_medium_kbase.py', 'biomassrxn_id':'biomass0','oxygen_exchrxn_id':'EX_cpd00007_e0','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}
    all_organisms['Pseudomonas_putida'] = {'model_filename':models_path + 'Pseudomonas_putida/iJP962/iJP962_updated.xml','model_id':'Pseudomonas_putida_iJP962_updated','growth_medium_filename':models_path + 'Pseudomonas_putida/iJP962/iJP962_minimal.py', 'biomassrxn_id':'IR10370','oxygen_exchrxn_id':'EX_EC0007','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}
    all_organisms['Pseudomonas_veronii'] = {'model_filename':models_path + 'Pseudomonas_veronii/kbase/Pseudomonas_veronii.13.xml','model_id':'Pseudomonas_veronii.13','growth_medium_filename':models_path + 'Pseudomonas_veronii/kbase/minimal_medium.py', 'biomassrxn_id':'biomass0','oxygen_exchrxn_id':'EX_cpd00007_e0','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}
    all_organisms['Serratia_marcescens'] = {'model_filename':models_path + 'Serratia_marcescens/kbase/Serratia_marcescens.12.xml','model_id':'Serratia_marcescens.12','growth_medium_filename':models_path + 'Serratia_marcescens/kbase/minimal_medium.py', 'biomassrxn_id':'biomass0','oxygen_exchrxn_id':'EX_cpd00007_e0','gDW_per_cell':2.8e-13,'cells_per_mL_init':1e7}

    # All possible pairs of these organisms
    org_pairs = list(itertools.combinations(all_organisms.keys(),r=2))
    print '\nTotal # of pairs to consider = ',len(org_pairs)
    for org_pair in org_pairs:
    #for org_pair in [('Pseudomonas_putida','Pseudomonas_fluorescens')]:
        print '\n************ ',org_pairs.index(org_pair) + 1,'. ',org_pair,' ************'
        members_info = {}
        for org in org_pair:
            members_info[org] = all_organisms[org]
        simulate_community(t0 = 0,delta_t = 0.5,tf = 5*48,members_info = members_info,carbon_sources = {'L_serine':'cpd00054','galacturonic_acid':'cpd00280'},media_cmp_concs = {'cpd00054':10,'cpd00280':10},serial_dilution_params = {'dilution_factor':1500,'time_between_dilutions':48},exchrxns_objcoeff = 0.001,results_filename_base = 'results/' + '_'.join(org_pair))

