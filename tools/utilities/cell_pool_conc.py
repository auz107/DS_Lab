from __future__ import division
import sys, time
sys.path.append('../../')
from tools.userError import userError
from tools.core.compound import compound
from tools.core.reaction import reaction
from tools.core.model import model
from remove_non_alphanumeric import remove_non_alphanumeric
from get_ModelSEED_ids import get_cmp_ModelSEED_id
import re

def cell_pool_conc(cmp_org, compart_list = None, stdout_msgs = True, warnings = True):
    """
    This functions returns the cell pool concentration of a a compound or a list
    of compounds 

    INPUTS:
    -------
          cmp_org: A dictionary where the keys are instances of the object compound
                   and values are the cooresponding organism (an instance of object
                   organism) 
      stdout_msgs: Prints a summary of the FBA and a number of other messages in the output if 'on.
                   Eligible values are True and False and the default is True 
         warnings: Can be True or False shwoing whether the warnings should be writtten to the 
                   screen or not. The default is True  


    OUTPUTS:
    --------
    The cell pool concentration of the compound. It returns zero if the compound
    is not in the list 

    Ali R. Zomorrodi - Daniel Segre Lab @ BU
    Last updated: 03-31-2015
    """
    if not isinstance(stdout_msgs,bool):
        raise TypeError("stdout_msgs must be either True or False")

    if not isinstance(warnings,bool):
        raise TypeError("warnings must be either True or False")

    # List of all compartments in the input compounds
    if compart_list == None:
        compart_list = list(set([c.compartment.id for c in cmp_org.keys() if c.compartment != None]))

    # Try to get the ModelSEED ids, if they are not given
    get_cmp_ModelSEED_id([c for c in cmp_org.keys() if c.ModelSEED_id == None],compart_list = compart_list,stdout_msgs = False)

    #---- Pool of compounds in E. coli ----
    # Reference: R. Raunio and H Rosenqvist (1970) "Amino acid pool of Escherchia coli during
    # the different phases of growth", Acta Chemica Scandinavica, 24, 2737-2744"
    # Concentrations are in mmol/gDW
    Ecoli = {}
   
    # L-Alanine 
    Ecoli['cpd00035'] = 6.7e-3

    # L-Valine
    Ecoli['cpd00156'] = 9.2e-3

    # Glycine
    Ecoli['cpd00033'] = 1.7e-3

    # L-isoleucine
    Ecoli['cpd00322'] = 1.6e-3

    # L-threonine
    Ecoli['cpd00161'] = 0.2e-3

    # L-leucine
    Ecoli['cpd00107'] = 0.4e-3

    # L-serine
    Ecoli['cpd00054'] = 0

    # L-proline 
    Ecoli['cpd00129'] = 0.4e-3

    # L-aspartate 
    Ecoli['cpd00041'] = 0

    # L-cystein 
    Ecoli['cpd00084'] = 1.8e-3

    # L-methionine 
    Ecoli['cpd00060'] = 0.3e-3

    # L-glutamate 
    Ecoli['cpd00023'] = 25.5e-3

    # L-phenylalanine 
    Ecoli['cpd00066'] = 6.8e-3

    # L-tyrosine 
    Ecoli['cpd00069'] = 0.1e-3

    # L-ornithine 
    Ecoli['cpd00064'] = 0

    # L-lysine 
    Ecoli['cpd00039'] = 0.2e-3

    # L-tryptophan 
    Ecoli['cpd00065'] = 0

    # L-arginine 
    Ecoli['cpd00051'] = 0

    # L-cystine 
    Ecoli['cpd00381'] = 0

    # compounds not found 
    not_found_cmp_ids = []

    for cmp in cmp_org.keys():
        cmp.cell_pool_concentration = {}
        for org in cmp_org[cmp]:
            if org.id in ['Ecoli', 'E. coli', 'Escherichia coli'] or org.name in ['Ecoli', 'E. coli', 'Escherichia coli']: 
                if cmp.ModelSEED_id in Ecoli.keys():
                    cmp.cell_pool_concentration[org.id] = Ecoli[cmp.ModelSEED_id]
                else:
                    cmp.cell_pool_concentration[org.id] = 0 
                    not_found_cmp_ids.append(cmp.id)
            else:
                cmp.cell_pool_concentration[org.id] = 0 
                if warnings: 
                    print '**WARNING! organism ',org.id,' was not found. Zero assigned too to cell_pool_concentration'

    if warnings: 
        print '**WARNING! These compounds were not found in the list of cell free compounds and zero assigned to their cell_pool_concentration:',list(set(not_found_cmp_ids)),'\n'
