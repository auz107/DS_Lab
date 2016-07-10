mutants_rxn_info_PamSilver = {} 
mutants_rxn_info_PamSilver['ppc'] = ['PPC']
mutants_rxn_info_PamSilver['argB'] = ['ACGK']
mutants_rxn_info_PamSilver['argC'] = ['AGPR']
mutants_rxn_info_PamSilver['argE'] = ['ACODA','NACODA']
mutants_rxn_info_PamSilver['argH'] = ['ARGSL']
mutants_rxn_info_PamSilver['proC'] = ['P5CR']
mutants_rxn_info_PamSilver['icd'] = ['ICDHyr']
mutants_rxn_info_PamSilver['nadC'] = ['NNDPR']
mutants_rxn_info_PamSilver['panB'] = ['MOHMT']
mutants_rxn_info_PamSilver['panD'] = ['ASP1DC']
mutants_rxn_info_PamSilver['serC'] = ['OHPBAT','PSERT']
mutants_rxn_info_PamSilver['cysC'] = ['ADSK']
mutants_rxn_info_PamSilver['cysE'] = ['SERAT']
mutants_rxn_info_PamSilver['glnA'] = ['GLNS']
mutants_rxn_info_PamSilver['hisB'] = ['HISTP','IGPDH']
mutants_rxn_info_PamSilver['hisC'] = ['HSTPT']
mutants_rxn_info_PamSilver['hisD'] = ['HISTD']
mutants_rxn_info_PamSilver['hisI'] = ['PRAMPC','PRATPP']
mutants_rxn_info_PamSilver['metA'] = ['HSST']
mutants_rxn_info_PamSilver['guaA'] = ['GMPS2']
mutants_rxn_info_PamSilver['purA'] = ['ADSS','LPADSS']
mutants_rxn_info_PamSilver['purC'] = ['PRASCSi']
mutants_rxn_info_PamSilver['purD'] = ['PRAGSr']
mutants_rxn_info_PamSilver['purF'] = ['GLUPRT']
mutants_rxn_info_PamSilver['purL'] = ['PRFGS']
mutants_rxn_info_PamSilver['purM'] = ['PRAIS']
mutants_rxn_info_PamSilver['pyrB'] = ['ASPCT']
mutants_rxn_info_PamSilver['pyrC'] = ['DHORTS']
mutants_rxn_info_PamSilver['pyrD'] = ['DHORD2','DHORD5']
mutants_rxn_info_PamSilver['pyrE'] = ['ORPT']
mutants_rxn_info_PamSilver['pyrF'] = ['OMPDC']
mutants_rxn_info_PamSilver['lysA'] = ['DAPDC']
mutants_rxn_info_PamSilver['tyrA'] = ['PPND']
mutants_rxn_info_PamSilver['trpA'] = ['TRPS1','TRPS2','TRPS3']
mutants_rxn_info_PamSilver['trpB'] = ['TRPS1','TRPS2','TRPS3']
mutants_rxn_info_PamSilver['trpC'] = ['IGPS','PRAIi']
mutants_rxn_info_PamSilver['trpD'] = ['ANPRT','ANS']
mutants_rxn_info_PamSilver['trpE'] = ['ANS']
mutants_rxn_info_PamSilver['ilvC'] = ['KARA1','KARA2']
mutants_rxn_info_PamSilver['ilvE'] = ['ILETA','VALTA']
mutants_rxn_info_PamSilver['leuB'] = ['IPMD','OMCDC']
mutants_rxn_info_PamSilver['leuD'] = ['IPPMIa','IPPMIb']
mutants_rxn_info_PamSilver['leuC'] = ['IPPMIa','IPPMIb']

#----------------------------------------------
mutants_rxn_info_AAs = {} 
genes_AA_map = {}

# ala_L
genes_AA_map['alaA_alaC_avtA'] = 'ala_L'
mutants_rxn_info_AAs['alaA_alaC_avtA'] = ['ALATA_L','VPAMTr']  

# arg_L
genes_AA_map['argH'] = 'arg_L'
mutants_rxn_info_AAs['argH'] = ['ARGSL']                 

# asn_L
genes_AA_map['asnA_asnB'] = 'asn_L'
mutants_rxn_info_AAs['asnA_asnB'] = ['ASNS1','ASNS2']   

# asp_L
genes_AA_map['aspC'] = 'asp_L'
mutants_rxn_info_AAs['aspC'] = ['ASPTA']               

# cys_L
genes_AA_map['cysC'] = 'cys_L'
mutants_rxn_info_AAs['cysC'] = ['ADSK']               

# gln_L
genes_AA_map['glnA'] = 'gln_L'
mutants_rxn_info_AAs['glnA'] = ['GLNS']              

# glu_L
genes_AA_map['gltBD_gdhA'] = 'glu_L'
mutants_rxn_info_AAs['gltBD_gdhA'] = ['GLUSy','GLUDy']  

# gly
genes_AA_map['glyA_ltaE'] = 'gly'
mutants_rxn_info_AAs['glyA_ltaE'] = ['GHMT2r', 'GLYAT', 'THRAi','THRA2i']  

# his_L
genes_AA_map['hisB'] = 'his_L'
mutants_rxn_info_AAs['hisB'] = ['HISTP','IGPDH'] 
genes_AA_map['hisC'] = 'his_L'
mutants_rxn_info_AAs['hisC'] = ['HSTPT']
genes_AA_map['hisD'] = 'his_L'
mutants_rxn_info_AAs['hisD'] = ['HISTD']
genes_AA_map['hisI'] = 'his_L'
mutants_rxn_info_AAs['hisI'] = ['PRAMPC','PRATPP']

# ile_L
genes_AA_map['ilvE'] = 'ile_L'
mutants_rxn_info_AAs['ilvE'] = ['ILETA','VALTA']       

# leu_L
genes_AA_map['leuB'] = 'leu_L'
mutants_rxn_info_AAs['leuB'] = ['IPMD','OMCDC']       
genes_AA_map['leuD'] = 'leu_L'
mutants_rxn_info_AAs['leuD'] = ['IPPMIa','IPPMIb']   
genes_AA_map['leuC'] = 'leu_L'
mutants_rxn_info_AAs['leuC'] = ['IPPMIa','IPPMIb']  

# lys_L
genes_AA_map['lysA'] = 'lys_L'
mutants_rxn_info_AAs['lysA'] = ['DAPDC']               

# met_L
genes_AA_map['metC_malY'] = 'met_L'
mutants_rxn_info_AAs['metC_malY'] = ['CYSTL']            

# phe_L
genes_AA_map['tyrB_ilvE_aspC'] = 'phe_L'
mutants_rxn_info_AAs['tyrB_ilvE_aspC'] = ['PHETA1']     

# pro_L
genes_AA_map['proC'] = 'pro_L'
mutants_rxn_info_AAs['proC'] = ['P5CR']                

# ser_L
genes_AA_map['serB_glyA'] = 'ser_L'
mutants_rxn_info_AAs['serB_glyA'] = ['PSP_L','GHMT2r']             

# thr_L
genes_AA_map['thrC'] = 'thr_L'
mutants_rxn_info_AAs['thrC'] = ['THRS', '4HTHRS']     

# trp_L
genes_AA_map['trpA'] = 'trp_L'
mutants_rxn_info_AAs['trpA'] = ['TRPS1','TRPS2','TRPS3'] 
genes_AA_map['trpB'] = 'trp_L'
mutants_rxn_info_AAs['trpB'] = ['TRPS1','TRPS2','TRPS3']

# tyr_L
genes_AA_map['tyrA'] = 'tyr_L'
mutants_rxn_info_AAs['tyrA'] = ['PPND']                

# val_L
genes_AA_map['ilvE_avtA'] = 'val_L'
mutants_rxn_info_AAs['ilvE_avtA'] = ['VALTA', 'VPAMTr'] 

