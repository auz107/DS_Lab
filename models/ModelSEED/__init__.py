__all__ = ['KEGG_rxns_pathways', 'ModelSEED_compartments', 'ModelSEED_cpds_BiGG_aliases', 'ModelSEED_cpds_KEGG_aliases', 'ModelSEED_cpds_formula_aliases', 'ModelSEED_cpds_master', 'ModelSEED_cpds_name_aliases', 'ModelSEED_cpds_website', 'ModelSEED_rxns_BiGG_aliases', 'ModelSEED_rxns_EC_number_aliases', 'ModelSEED_rxns_GramNegative', 'ModelSEED_rxns_GramNegative_Biomass', 'ModelSEED_rxns_GramPositive', 'ModelSEED_rxns_GramPositive_Biomass', 'ModelSEED_rxns_Human', 'ModelSEED_rxns_Human_Biomass', 'ModelSEED_rxns_KEGG_aliases', 'ModelSEED_rxns_master', 'ModelSEED_rxns_name_aliases', 'ModelSEED_rxns_website']  

# Add Directories to PYTHONPATH
import sys, os
dir = os.path.dirname(__file__)

for subpackage in __all__:
   sys.path.append(os.path.join(dir,subpackage))

