excess_nutrients = {
#-- Compounds in the minimal medium --
'EX_ca2(e)':[-1000,1000],
'EX_cbl1(e)':[-1000,1000],
'EX_cl(e)':[-1000,1000],
'EX_co2(e)':[-1000,1000],
'EX_cobalt2(e)':[-1000,1000],
'EX_cu2(e)':[-1000,1000],
'EX_fe2(e)':[-1000,1000],
'EX_fe3(e)':[-1000,1000],
'EX_h(e)':[-1000,1000],
'EX_h2o(e)':[-1000,1000],
'EX_k(e)':[-1000,1000],
'EX_mg2(e)':[-1000,1000],
'EX_mn2(e)':[-1000,1000],
'EX_mobd(e)':[-1000,1000],
'EX_na1(e)':[-1000,1000],
'EX_nh4(e)':[-1000,1000],
'EX_ni2(e)':[-1000,1000],
'EX_pi(e)':[-1000,1000],
'EX_sel(e)':[-1000,1000],
'EX_slnt(e)':[-1000,1000],
'EX_so4(e)':[-1000,1000],
'EX_tungs(e)':[-1000,1000],
'EX_zn2(e)':[-1000,1000]

#-- Compounds in the rich medium --
# All amino acids
'EX_ala_L(e)':[-20,1000],
'EX_arg_L(e)':[-20,1000],
'EX_asn_L(e)':[-20,1000],
'EX_asp_L(e)':[-20,1000],
'EX_cys_L(e)':[-20,1000],
'EX_gln_L(e)':[-20,1000],
'EX_glu_L(e)':[-20,1000],
'EX_his_L(e)':[-20,1000],
'EX_ile_L(e)':[-20,1000],
'EX_leu_L(e)':[-20,1000],
'EX_lys_L(e)':[-20,1000],
'EX_met_L(e)':[-20,1000],
'EX_phe_L(e)':[-20,1000],
'EX_pro_L(e)':[-20,1000],
'EX_ser_L(e)':[-20,1000],
'EX_thr_L(e)':[-20,1000],
'EX_trp_L(e)':[-20,1000],
'EX_tyr_L(e)':[-20,1000],
'EX_val_L(e)':[-20,1000],
'EX_gly(e)':[-20,1000],

# Four nucleotides bases
'EX_ade(e)':[-20,1000],
'EX_thym(e)':[-20,1000],
'EX_gua(e)':[-20,1000],
'EX_csn(e)':[-20,1000]

}

regulation = {
'CADVtpp':[0,0],
'DMSOR1':[0,0],
'DMSOR2':[0,0],
'FRD2':[0,0],
'FRD3':[0,0],
'HYD1pp':[0,0],
'HYD2pp':[0,0],
'HYD3pp':[0,0],
'NO2t2rpp':[0,0],
'NTRIR2x':[0,0],
'PFL':[0,0],
'RNTR1c':[0,0],
'RNTR2c':[0,0],
'RNTR3c':[0,0],
'RNTR4c':[0,0],
'TMAOR1':[0,0],
'TMAOR2':[0,0],
}

ngam_atp = {
'ATPM':[8.39000000,8.39000000],
}

others = {
'ALAt2pp':[0,0],    # Creates a thermodynamically infeasible loop with ALAt2rpp 
'GLYt2pp':[0,0]    # Creates a thermodynamically infeasible loop with GLYt2rpp
}

