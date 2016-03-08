% NOTES:
% iIB711 is an older version of the model for Streptomyces coelicolor and was not used
% iAN840m Unclear to what orgganism it belongs
% iJN746 and iJP815 are older versions of the Pseudomonas putida model
% iJR904 is an older version of the E. coli model
% iLC915 is an older version of the Pichia pastoris model
% iMB745 is an older version of the Methanosarcina acetivorans model
% iMM904 and iND750 are older versions of the S. cerevisiae model
% iNJ661 is an older version of the Mycobacterium tuberculosis model
% AraGEM is an older version of the Arabidopsis thaliana model
% iSB619 is an older version of the Staphylococcus aureus model
% PpuMBEL1071 is an older version of the model for Pseudomonas putida 
% iSyn669 is an older version of the Synechocystis sp. PCC6803 model 
% iCR744 is an older version of the Rhodoferax ferrireducens model
% iRM588 and iJS747 are older versions of the Geobacter sulfurreducens model
% mus_musculus is an older version of the Mus musculus model

clear all

% Load all models from the file
load('all_models.mat');

bigg_models= struct();

model_names = {};

%------ Bacteria ---------
% Acinetobacter baumannii
Acinetobacter_baumannii = AbyMBEL891;
Acinetobacter_baumannii.organism = {'Acinetobacter baumannii','Bacteria'};
Acinetobacter_baumannii.description = 'AbyMBEL891';

% Salmonella enterica%
Salmonella_enterica = STM_v1_0; 
Salmonella_enterica.organism = {'Salmonella enterica','Bacteria'};
Salmonella_enterica.description = 'STM_v1_0'; 

% Streptomyces coelicolor (A newer model called iMK1208 is available at PMID: 24623710)
Streptomyces_coelicolor = S_coilicolor_fixed;
Streptomyces_coelicolor.organism = {'Streptomyces coelicolor','Bacteria'};
Streptomyces_coelicolor.description = 'S_coilicolor';

% Thermotoga maritima
Thermotoga_maritima = T_Maritima;
Thermotoga_maritima.organism = {'Thermotoga_maritima','Bacteria'};
Thermotoga_maritima.description = 'Thermotoga_maritima';

% Vibrio vulnificus
Vibrio_vulnificus = VvuMBEL943;
Vibrio_vulnificus.organism = {'Vibrio_vulnificus','Bacteria'};
Vibrio_vulnificus.description = 'VvuMBEL943';

% Escherichia coli
Escherichia_coli_iAF1260 = iAF1260;
Escherichia_coli_iAF1260.organism = {'Escherichia coli','Bacteria'};
Escherichia_coli_iAF1260.description = 'iAF1260';
Escherichia_coli_iJO1366 = iJO1366;
Escherichia_coli_iJO1366.organism = {'Escherichia coli','Bacteria'};
Escherichia_coli_iJO1366.description = 'iJO1366';

% Escherichia coli W
Escherichia_coli_W  = iCA1273;
Escherichia_coli_W.organism  = {'Escherichia coli W','Bacteria'};
Escherichia_coli_W.description  = 'iCA1273'; 

% Dehalococcoides ethenogenes
Dehalococcoides_ethenogenes = iAI549;
Dehalococcoides_ethenogenes.organism = {'Dehalococcoides ethenogenes','Bacteria'};
Dehalococcoides_ethenogenes.description = 'iAI549'; 

% Bacillus subtilis
Bacillus_subtilis = iBsu1103;
Bacillus_subtilis.organism = {'Bacillus subtilis','Bacteria'};
Bacillus_subtilis.description = 'iBsu1103';

% Clostridium beijerinckii
Clostridium_beijerinckii = iCB925;
Clostridium_beijerinckii.organism = {'Clostridium beijerinckii','Bacteria'};
Clostridium_beijerinckii.description = 'iCB925'; 

% Clostridium acetobutylicum 
Clostridium_acetobutylicum  = iCac802;
Clostridium_acetobutylicum.organism = {'Clostridium_acetobutylicum','Bacteria'};
Clostridium_acetobutylicum.description = 'iCac802';

% Helicobacter pylori
Helicobacter_pylori = iIT341;
Helicobacter_pylori.organism = {'Helicobacter pylori','Bacteria'};
Helicobacter_pylori.description = 'iIT341';

% Synechocystis sp. PCC6803
Synechocystis_sp = iJN678;
Synechocystis_sp.organism = {'Synechocystis sp. PCC6803','Bacteria'};
Synechocystis_sp.description = 'iJN678';

% Burkholderia cenocepacia
Burkholderia_cenocepacia = iKF1028;
Burkholderia_cenocepacia.organism = {'Burkholderia cenocepacia','Bacteria'};
Burkholderia_cenocepacia.description = 'iKF1028';

% Pseudomonas aeruginosa
Pseudomonas_aeruginosa = iMO1056;
Pseudomonas_aeruginosa.organism = {'Pseudomonas aeruginosa','Bacteria'};
Pseudomonas_aeruginosa.description = 'iMO1056';

% Mycobacterium tuberculosis
Mycobacterium_tuberculosis = iNJ661m;
Mycobacterium_tuberculosis.organism = {'Mycobacterium tuberculosis','Bacteria'};
Mycobacterium_tuberculosis.description = 'iNJ661m';

% Mycoplasma genitalium 
Mycoplasma_genitalium = iPS189_fixed;
Mycoplasma_genitalium.organism = {'Mycoplasma genitalium','Bacteria'};
Mycoplasma_genitalium.description = 'iPS189';

% Rhodobacter sphaeroides
Rhodobacter_sphaeroides = iRsp1095;
Rhodobacter_sphaeroides.organism = {'Rhodobacter sphaeroides','Bacteria'};
Rhodobacter_sphaeroides.description = 'iRsp1095';

% Clostridium thermocellum
Clostridium_thermocellum = iSR432;
Clostridium_thermocellum.organism = {'Clostridium thermocellum','Bacteria'};
Clostridium_thermocellum.description = 'iSR432';

% Plasmodium falciparum HB3 
Plasmodium_falciparum = iTH366;
Plasmodium_falciparum.organism = {'Plasmodium falciparum HB3','Bacteria','Parasite'};
Plasmodium_falciparum.description = 'iTH366';

% Klebsiella pneumoniae
Klebsiella_pneumoniae = iYL1228;
Klebsiella_pneumoniae.organism = {'Klebsiella pneumoniae','Bacteria'};
Klebsiella_pneumoniae.description = 'iYL1228';

% Rhizobium etli
Rhizobium_etli = iOR363;
Rhizobium_etli.organism = {'Rhizobium etli','Bacteria'};
Rhizobium_etli.description = 'iOR363';

% Shewanella oneidensis 
Shewanella_oneidensis = iSO783;
Shewanella_oneidensis.organims = {'Shewanella oneidensis','Bacteria'};
Shewanella_oneidensis.description = 'iSO783';

% Ketogulonicigenium vulgare 
Ketogulonicigenium_vulgare = iWZ663;
Ketogulonicigenium_vulgare.organism = {'Ketogulonicigenium vulgare','Bacteria'};
Ketogulonicigenium_vulgare.description = 'iWZ663';

%------ Archaea -------
% Methanosarcina barkeri
Methanosarcina_barkeri = iAF692;
Methanosarcina_barkeri.organism = {'Methanosarcina barkeri','Archaea'};
Methanosarcina_barkeri.description = 'iAF692';

% Methanosarcina acetivorans
Methanosarcina_acetivorans = iVS941_fixed;
Methanosarcina_acetivorans.organism = {'Methanosarcina acetivorans','Archaea'};
Methanosarcina_acetivorans.description = 'iVS941';

%------ Eukaryota -------
% Aspergillus oryzae
Aspergillus_oryzae = AORYZAE_COBRA;
Aspergillus_oryzae.organism = {'Aspergillus oryzae','Eukaryota','Fungus'};
Aspergillus_oryzae.description = 'AORYZAE';

% Arabidopsis thaliana
Arabidopsis_thaliana = iRS1597;
Arabidopsis_thaliana.organism = {'Arabidopsis thaliana','Eukaryota','Plant'};
Arabidopsis_thaliana.description = 'iRS1597';

% Pichia pastoris
Pichia_pastoris = PpaMBEL1254;
Pichia_pastoris.organism = {'Pichia pastoris','Eukaryota','Yeast'};
Pichia_pastoris.description = 'PpaMBEL1254';

% Schizosaccharomyces pombe
Schizosaccharomyces_pombe = SpoMBEL1693;
Schizosaccharomyces_pombe.organism = {'Schizosaccharomyces pombe','Eukaryota','Yeast'};
Schizosaccharomyces_pombe.description = 'SpoMBEL1693';

% Leishmania major
Leishmania_major = iAC560;
Leishmania_major.organism = {'Leishmania major','Eukaryota','Unicellular'};
Leishmania_major.description = 'iAC560';

% Aspergillus niger
Aspergillus_niger = iMA871;
Aspergillus_niger.organism = {'Aspergillus niger','Eukaryota','Fungus'};
Aspergillus_niger.description = 'iMA871';

% Mus musculus
Mus_musculus = iMM1415;
Mus_musculus.organism = {'Mus musculus','Eukaryota','Mouse'};
Mus_musculus.description = 'iMM1415';

% Chlamydomonas reinhardtii
Chlamydomonas_reinhardtii = iRC1080;
Chlamydomonas_reinhardtii.organism = {'Chlamydomonas reinhardtii','Eukaryota','Unicellular_Alga'};
Chlamydomonas_reinhardtii.description = 'iRC1080';

% Zea mays
Zea_mays = iRS1563;
Zea_mays.organism = {'Zea mays','Eukaryota','Plant'};
Zea_mays.description = 'iRS1563';

% Pichia stipitis
Pichia_stipitis = iSS884;
Pichia_stipitis.organism = {'Pichia stipitis','Eukaryota','Yeast'};
Pichia_stipitis.description = 'iSS884';

% Cryptosporidium hominis 
Cryptosporidium_hominis = iNV213;
Cryptosporidium_hominis.organism = {'Cryptosporidium hominis','Microorganism'};
Cryptosporidium_hominis.description = 'iNV213';

% Scheffersomyces stipitis (Pichia stipitis)
Scheffersomyces_stipitis = iTL885;
Scheffersomyces_stipitis.organism = {'Scheffersomyces stipitis','Yeast'};
Scheffersomyces_stipitis.description = 'iTL885';
;
all_models_names = {'Acinetobacter_baumannii','Salmonella_enterica','Streptomyces_coelicolor','Thermotoga_maritima','Vibrio_vulnificus','Escherichia_coli_iAF1260','Escherichia_coli_iJO1366','Escherichia_coli_W','Dehalococcoides_ethenogenes','Bacillus_subtilis','Clostridium_beijerinckii','Clostridium_acetobutylicum','Helicobacter_pylori','Synechocystis_sp_PCC6803','Burkholderia_cenocepacia','Pseudomonas_aeruginosa','Mycobacterium_tuberculosis','Rhodobacter_sphaeroides','Clostridium_thermocellum','Plasmodium_falciparum','Klebsiella_pneumoniae','Rhizobium_etli','Shewanella_oneidensis','Ketogulonicigenium_vulgare','Methanosarcina_barkeri','Methanosarcina_acetivorans','Aspergillus_oryzae','Pichia_pastoris','Schizosaccharomyces_pombe','Leishmania_major','Aspergillus_niger','Mus_musculus','Chlamydomonas_reinhardtii','Zea_mays','Pichia_stipitis','Cryptosporidium_hominis','Scheffersomyces_stipitis'};

fprintf('\nThe total number of models = %i\n',length(all_models_names))

save('all_models_refined.mat','Acinetobacter_baumannii','Salmonella_enterica','Streptomyces_coelicolor','Thermotoga_maritima','Vibrio_vulnificus','Escherichia_coli_iAF1260','Escherichia_coli_iJO1366','Escherichia_coli_W','Dehalococcoides_ethenogenes','Bacillus_subtilis','Clostridium_beijerinckii','Clostridium_acetobutylicum','Helicobacter_pylori','Synechocystis_sp_PCC6803','Burkholderia_cenocepacia','Pseudomonas_aeruginosa','Mycobacterium_tuberculosis','Rhodobacter_sphaeroides','Clostridium_thermocellum','Plasmodium_falciparum','Klebsiella_pneumoniae','Rhizobium_etli','Shewanella_oneidensis','Ketogulonicigenium_vulgare','Methanosarcina_barkeri','Methanosarcina_acetivorans','Aspergillus_oryzae','Pichia_pastoris','Schizosaccharomyces_pombe','Leishmania_major','Aspergillus_niger','Mus_musculus','Chlamydomonas_reinhardtii','Zea_mays','Pichia_stipitis','Cryptosporidium_hominis','Scheffersomyces_stipitis');


