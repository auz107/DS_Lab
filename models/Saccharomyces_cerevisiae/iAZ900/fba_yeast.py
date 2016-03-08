import sys
sys.path.append('../../../')
from tools.io.read_sbml_model import read_sbml_model
from tools.fba.set_specific_bounds import set_specific_bounds
from tools.fba.fba import fba
from tools.core.model import model
from tools.core.organism import organism
import cobra

def write_gams_model():
    with open('metabolites.txt','w') as f:
        f.write('\/\n')
        for cmp in iAZ900_noCycles.compounds:
            f.write("'" + cmp.id + "'\n")
        f.write('\/\n')

    with open('reactions.txt','w') as f:
        f.write('\/\n')
        for rxn in iAZ900_noCycles.reactions:
            f.write("'" + rxn.id + "'\n")
        f.write('\/\n')

    with open('reaction_types.txt','w') as f:
        f.write('\/\n')
        for rxn in iAZ900_noCycles.reactions:
            if rxn.type.lower() == 'irreversible':
                f.write("'" + rxn.id + "'\t0\n")
            elif rxn.type.lower() == 'reversible':
                f.write("'" + rxn.id + "'\t1\n")
            elif rxn.type.lower() == 'exchange':
                f.write("'" + rxn.id + "'\t3\n")
            else:
                raise userError('unknonw reaction name')
        f.write('\/\n')

    with open('stoichiometric_matrix.txt','w') as f:
        for cmp in iAZ900_noCycles.compounds:
            for rxn in cmp.reactant_reactions:
                f.write("S('" + cmp.id + "','" + rxn.id + "') = " + rxn.stoichiometry[cmp] +';\n')
            for rxn in cmp.product_reactions:
                f.write("S('" + cmp.id + "','" + rxn.id + "') = " + rxn.stoichiometry[cmp] +';\n')
            f.write('\n')

def compare_models():
    print '\nrxn id\t\tiAZ900\t\tiAZ900_noCycles'
    for iAZ900_rxn in iAZ900.reactions:
        iAZ900_noCycles_rxn = iAZ900_noCycles.get_reactions({iAZ900_rxn.id:'id'})
        if iAZ900_rxn.type != iAZ900_noCycles_rxn.type:
            print iAZ900_rxn.id,'\t\t',iAZ900_rxn.type,'\t\t',iAZ900_noCycles_rxn.type
          
    print '\n'

    print '\nCompartments'
    for compart in iAZ900.compartments:
        print compart.id,'  ',compart.name

    print '\niAZ900'
    print '# of reactions:',len(iAZ900.reactions)
    print '# of metabolites:',len(iAZ900.compounds)
    print '# of genes:',len(iAZ900.genes)
    print '# of compartments:',len(iAZ900.compartments)

    print '\niAZ900_noCycles'
    print '# of reactions:',len(iAZ900_noCycles.reactions)
    print '# of metabolites:',len(iAZ900_noCycles.compounds)
    print '# of genes:',len(iAZ900_noCycles.genes)
    print '# of compartments:',len(iAZ900.compartments)


if __name__ == "__main__":

    # Optimizaiton solver
    optimization_solver = 'gurobi'

    # Model path
    model_path = '/usr2/postdoc/alizom/work/models/Saccharomyces_cerevisiae/iAZ900/'

    # Define the organism
    model_organism = organism(id = 'Scerevisiae', name = 'Saccharomyces cerevisiae',domain = 'Eukaryotes', genus = 'Saccharomyces', species = 'cerevisiae', strain = '')

    #------- Original iAZ900 model -----
    print '\n----- Original iAZ900 model -----\n'
    # Load the original model
    iAZ900 = read_sbml_model(file_name = model_path + 'iAZ900.xml', model_name = 'iAZ900',model_organism = model_organism, model_type = 'metabolic',import_params = False)

    # Set the growth medium (minmal aerobic glucose
    set_specific_bounds(model = iAZ900,file_name = model_path + 'iAZ900_minimal.py',limiting_nutrients = {'EX_glc_e_':[-10,1000],'EX_o2_e_':[-2,1000]})

    for rxn in iAZ900.reactions:
        if rxn.id == 'biomass_wildType':
        #if rxn.id == 'biomass_core':
            rxn.objective_coefficient = 1
        else:
            rxn.objective_coefficient = 0

    # FBA  
    iAZ900.fba(optimization_solver = optimization_solver) 

    #------- iAZ900 model modified with cycles removed (by Anupam) -----
    print '\n----- iAZ900_noCycles model modified with cycles removed (by Anupam) -----\n'
    # Load the original model
    iAZ900_noCycles = read_sbml_model(file_name = model_path + 'iAZ900_noCycles_03_25_2015.xml', model_name = 'iAZ900_noCycles',model_organism = model_organism, model_type = 'metabolic',import_params = False)

    # Set the growth medium (minimal aerobic glucose)
    set_specific_bounds(model = iAZ900_noCycles,file_name = model_path + 'iAZ900_minimal.py',limiting_nutrients = {'EX_glc_e_':[-10,1000],'EX_o2_e_':[-2,1000]})

    for rxn in iAZ900_noCycles.reactions:
        rxn.objective_coefficient = 0

    biomass_rxn = iAZ900_noCycles.get_reactions({'biomass_wildType':'id'}) 
    #biomass_rxn = iAZ900_noCycles.get_reactions({'biomass_core':'id'}) 
    biomass_rxn.objective_coefficient = 1


    # FBA  
    iAZ900_noCycles.fba(optimization_solver = optimization_solver) 
  
    #-------- Change sucrose hydrolyzing reaciton (add ATP cost) ---
    SUCRe = iAZ900_noCycles.get_reactions({'SUCRe':'id'})
    atp = iAZ900_noCycles.get_compounds({'atp_c':'id'})
    adp = iAZ900_noCycles.get_compounds({'adp_c':'id'})
    pi = iAZ900_noCycles.get_compounds({'pi_c':'id'})
    # According to Jeff Gore's paper (suppl. figure 5) the cooperators are 2.5% less fit than cheasters
    # In the abscne of a metabolic cost for sucrose hydrolysis (SUCRe) the growth rate is 0.5023. We adjust
    # the coefficinet of ATP in SUCRe such the growth becomes 97.5% of 0.5023 = 0.48974
    SUCRe.add_compounds({atp:-0.12,adp:0.12,pi:0.12})

    # Set the growth medium (minimal aerobic glucose)
    set_specific_bounds(model = iAZ900_noCycles,file_name = model_path + 'iAZ900_minimal.py',limiting_nutrients = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-2,1000]})

    iAZ900_noCycles.fba(optimization_solver = optimization_solver, create_model = True) 

    #--- Simulate histidone defect ---
    print '\n--- Simulation histidine synthesis defect ---'
    set_specific_bounds(model = iAZ900_noCycles,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_o2_e_':[-2,1000],'IGPDH':[0,0]})
    iAZ900_noCycles.fba(optimization_solver = optimization_solver, create_model = True) 

    set_specific_bounds(model = iAZ900_noCycles,file_name = model_path + 'iAZ900_minimal.py',flux_bounds = {'EX_sucr_e_':[-10,1000],'EX_his_L_e_':[-1,1000],'EX_o2_e_':[-2,1000],'IGPDH':[0,0]})
    iAZ900_noCycles.fba(optimization_solver = optimization_solver, create_model = True) 

