#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%% ===========================================================================
# Description
# =============================================================================

'''
This script serves as a module to contain commonly used functions across other scripts.
Functions in this module:
    - KTC_pos_to_gene
        - Takes as input the path to a gtf file and a list of genomic positions (a tuple of (chr, pos))
        - It returns a list of gene names for those positions
'''

from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from adjustText import adjust_text
from gseapy import Msigdb
import gseapy as gp
import requests


#%% ===========================================================================
# KTC_pos_to_gene
# =============================================================================
# Takes as input the path to a gtf file and a list of genomic positions (tuples of (chr, pos))
# e.g. [(chr1, 1234),(chr2, 2345),(chr3, 3456)]
# Function should handle cases with and without 'chr' prefix
# Returns a list of gene names for those positions

def KTC_pos_to_gene(path_GTF, list_pos):
    import pandas as pd

    df_gtf            = pd.read_csv(path_GTF, sep='\t', comment='#', header=None)
    df_gtf.columns    = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df_gtf['seqname'] = df_gtf['seqname'].str.lower().str.replace('chr', '')
    df_gtf            = df_gtf[df_gtf['feature'] == 'transcript']
    gene_names        = []

    print('Annotating genes...')
    for chromosome, position in tqdm(list_pos, unit=' positions'):
        normalized_chromosome = chromosome.lower().replace('chr', '')
        gene_at_position = df_gtf[
            (df_gtf['seqname'] == normalized_chromosome) &
            (df_gtf['start'] <= int(position)) &
            (df_gtf['end'] >= int(position))
            ]
        if not gene_at_position.empty:
            gene_name = gene_at_position['attribute'].str.extract(r'gene_name "([^"]+)"')[0].values[0]
            gene_names.append(gene_name)
        else:
            gene_names.append(None)

    return(gene_names)

#%%
#%% ===========================================================================
# KTC_GetGeneSet
# =============================================================================

gene_sets = {
    'None' : [],
    'm6a_readers'      : ['HNRNPC', 'YTHDF1', 'YTHDF2', 'YTHDF3', 'YTHDC1', 'YTHDC2', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'HNRNPA2B1'],
    'splicing_factors' : ['SFRS7', 'CUGBP1', 'DAZAP1', 'CUGBP2', 'FMR1', 'A2BP1', 'RBFOX2', 'HNRNPA0', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPC', 'HNRNPC', 'HNRNPD', 'HNRNPD', 'HNRPDL', 'PCBP1', 'PCBP2', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPH3', 'PTBP1', 'HNRNPK', 'HNRNPK', 'HNRNPL', 'HNRPLL', 'HNRNPM', 'FUS', 'HNRNPU', 'TRA2A', 'TRA2B', 'ELAVL2', 'ELAVL4', 'ELAVL1', 'KHSRP', 'MBNL1', 'NOVA1', 'NOVA2', 'PTBP2', 'SFPQ', 'RBM25', 'RBM4', 'KHDRBS1', 'SF3B1', 'SFRS2', 'SF1', 'SFRS1', 'KHDRBS2', 'KHDRBS3', 'SFRS3', 'SFRS9', 'SFRS13A', 'SFRS5', 'SFRS11', 'SFRS6', 'SFRS4', 'TARDBP', 'TIA1', 'TIAL1', 'YBX1', 'ZRANB2', 'ELAVL3', 'RBM5', 'SYNCRIP', 'HNRNPA3', 'QKI', 'RBMX', 'SRRM1', 'ESRP1', 'ESRP2'], # From SpliceAid
    'EpiFactors'       : ['A1CF', 'ACINU', 'ACTB', 'ACTL6A', 'ACTL6B', 'ACTR3B', 'ACTR5', 'ACTR6', 'ACTR8', 'ADNP', 'AEBP2', 'AICDA', 'AIRE', 'ALKBH1', 'ALKBH1', 'ALKBH4', 'ALKBH5', 'ANKRD32', 'ANP32A', 'ANP32B', 'ANP32E', 'APBB1', 'APEX1', 'APOBEC1', 'APOBEC2', 'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H', 'ARID1A', 'ARID1B', 'ARID2', 'ARID4A', 'ARID4B', 'ARNTL', 'ARRB1', 'ASF1A', 'ASF1B', 'ASH1L', 'ASH2L', 'ASXL1', 'ASXL2', 'ASXL3', 'ATAD2', 'ATAD2B', 'ATF2', 'ATF7IP', 'ATM', 'ATN1', 'ATR', 'ATRX', 'ATXN7', 'ATXN7L3', 'AURKA', 'AURKB', 'AURKC', 'BABAM1', 'BAHD1', 'BANP', 'BAP1', 'BARD1', 'BAZ1A', 'BAZ1B', 'BAZ2A', 'BAZ2B', 'BCOR', 'BCORL1', 'BMI1', 'BPTF', 'BRCA1', 'BRCA2', 'BRCC3', 'BRD1', 'BRD2', 'BRD3', 'BRD4', 'BRD7', 'BRD8', 'BRD9', 'BRDT', 'BRE', 'BRMS1', 'BRMS1L', 'BRPF1', 'BRPF3', 'BRWD1', 'BRWD3', 'BUB1', 'C11orf30', 'C14orf169', 'C17orf49', 'CARM1', 'CBLL1', 'CBX1', 'CBX2', 'CBX3', 'CBX4', 'CBX5', 'CBX6', 'CBX7', 'CBX8', 'CCDC101', 'CDC6', 'CDC73', 'CDK1', 'CDK17', 'CDK2', 'CDK3', 'CDK5', 'CDK7', 'CDK9', 'CDY1', 'CDY1B', 'CDY2A', 'CDY2B', 'CDYL', 'CDYL2', 'CECR2', 'CELF1', 'CELF2', 'CELF3', 'CELF4', 'CELF5', 'CELF6', 'CENPC', 'CHAF1A', 'CHAF1B', 'CHD1', 'CHD1L', 'CHD2', 'CHD3', 'CHD4', 'CHD5', 'CHD6', 'CHD7', 'CHD8', 'CHD9', 'CHEK1', 'CHRAC1', 'CHTOP', 'CHUK', 'CIR1', 'CIT', 'CLNS1A', 'CLOCK', 'CRB2', 'CREBBP', 'CSNK2A1', 'CSRP2BP', 'CTBP1', 'CTBP2', 'CTCF', 'CTCFL', 'CTR9', 'CUL1', 'CUL2', 'CUL3', 'CUL4A', 'CUL4B', 'CUL5', 'CXXC1', 'DAPK3', 'DAXX', 'DDB1', 'DDB2', 'DDX17', 'DDX21', 'DDX5', 'DDX50', 'DEK', 'DHX9', 'DMAP1', 'DNAJC1', 'DNAJC2', 'DND1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DNMT3L', 'DNTTIP2', 'DOT1L', 'DPF1', 'DPF2', 'DPF3', 'DPPA3', 'DPY30', 'DR1', 'DTX3L', 'DZIP3', 'E2F6', 'EED', 'EEF1AKMT3', 'EEF1AKMT4', 'EEF1AKNMT', 'EHMT1', 'EHMT2', 'EID1', 'EID2', 'EID2B', 'EIF4A3', 'ELP2', 'ELP3', 'ELP4', 'ELP5', 'ELP6', 'ENY2', 'EP300', 'EP400', 'EPC1', 'EPC2', 'ERBB4', 'ERCC6', 'EXOSC1', 'EXOSC2', 'EXOSC3', 'EXOSC4', 'EXOSC5', 'EXOSC6', 'EXOSC7', 'EXOSC8', 'EXOSC9', 'EYA1', 'EYA2', 'EYA3', 'EYA4', 'EZH1', 'EZH2', 'FAM175A', 'FAM175B', 'FBL', 'FBRS', 'FBRSL1', 'FOXA1', 'FOXO1', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXP4', 'FTO', 'GADD45A', 'GADD45B', 'GADD45G', 'GATAD1', 'GATAD2A', 'GATAD2B', 'GFI1', 'GFI1B', 'GLYR1', 'GSE1', 'GSG2', 'GTF2I', 'GTF3C4', 'HAT1', 'HCFC1', 'HCFC2', 'HDAC1', 'HDAC10', 'HDAC11', 'HDAC2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6', 'HDAC7', 'HDAC8', 'HDAC9', 'HDGF', 'HDGFL2', 'HELLS', 'HIF1AN', 'HINFP', 'HIRA', 'HIRIP3', 'HJURP', 'HLCS', 'HLTF', 'HMG20A', 'HMG20B', 'HMGB1', 'HMGN1', 'HMGN2', 'HMGN3', 'HMGN4', 'HMGN5', 'HNRNPU', 'HNRPL', 'HNRPM', 'HP1BP3', 'HR', 'HSFX3', 'HSPA1A', 'HSPA1A', 'HSPA1B', 'HSPA1B', 'HUWE1', 'IKBKAP', 'IKZF1', 'IKZF3', 'ING1', 'ING2', 'ING3', 'ING4', 'ING5', 'INO80', 'INO80B', 'INO80C', 'INO80D', 'INO80E', 'JADE1', 'JADE2', 'JADE3', 'JAK2', 'JARID2', 'JDP2', 'JMJD1C', 'JMJD6', 'KANSL1', 'KANSL2', 'KANSL3', 'KAT2A', 'KAT2B', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'KAT8', 'KDM1A', 'KDM1B', 'KDM2A', 'KDM2B', 'KDM3A', 'KDM3B', 'KDM4A', 'KDM4B', 'KDM4C', 'KDM4D', 'KDM4E', 'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'KDM6A', 'KDM6B', 'KDM7A', 'KDM8', 'KEAP1', 'KHDRBS1', 'KLF18', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'KMT2E', 'L3MBTL1', 'L3MBTL2', 'L3MBTL3', 'L3MBTL4', 'LAS1L', 'LBR', 'LEO1', 'LRWD1', 'MAGOH', 'MAP3K7', 'MAPKAPK3', 'MASTL', 'MAX', 'MAZ', 'MBD1', 'MBD2', 'MBD3', 'MBD4', 'MBD5', 'MBD6', 'MBIP', 'MBNL1', 'MBNL3', 'MBTD1', 'MCRS1', 'MDC1', 'MEAF6', 'MECP2', 'MEN1', 'METTL11B', 'METTL14', 'METTL16', 'METTL21A', 'METTL3', 'METTL4', 'MGA', 'MGEA5', 'MINA', 'MIS18A', 'MIS18BP1', 'MLLT1', 'MLLT10', 'MLLT6', 'MORF4L1', 'MORF4L2', 'MOV10', 'MPHOSPH8', 'MPND', 'MRGBP', 'MSH6', 'MSL1', 'MSL2', 'MSL3', 'MST1', 'MTA1', 'MTA2', 'MTA3', 'MTF2', 'MUM1', 'MYBBP1A', 'MYO1C', 'MYSM1', 'NAA60', 'NAP1L1', 'NAP1L2', 'NAP1L4', 'NASP', 'NAT10', 'NAT10', 'NBN', 'NCL', 'NCOA1', 'NCOA2', 'NCOA3', 'NCOA6', 'NCOR1', 'NCOR2', 'NEK6', 'NEK9', 'NFRKB', 'NFYB', 'NFYC', 'NIPBL', 'NOC2L', 'NPAS2', 'NPM1', 'NPM2', 'NSD1', 'NSL1', 'NSRP1', 'NSUN2', 'NSUN6', 'NTMT1', 'NUP98', 'OGT', 'OIP5', 'PADI1', 'PADI2', 'PADI3', 'PADI4', 'PAF1', 'PAGR1', 'PAK2', 'PARG', 'PARP1', 'PARP2', 'PARP3', 'PAXIP1', 'PBK', 'PBRM1', 'PCGF1', 'PCGF2', 'PCGF3', 'PCGF5', 'PCGF6', 'PCNA', 'PDP1', 'PELP1', 'PHC1', 'PHC2', 'PHC3', 'PHF1', 'PHF10', 'PHF12', 'PHF13', 'PHF14', 'PHF19', 'PHF2', 'PHF20', 'PHF20L1', 'PHF21A', 'PHF8', 'PHIP', 'PIWIL4', 'PKM', 'PKN1', 'POGZ', 'POLE3', 'PPARGC1A', 'PPM1G', 'PPP2CA', 'PPP4C', 'PPP4R2', 'PQBP1', 'PRDM1', 'PRDM11', 'PRDM12', 'PRDM13', 'PRDM14', 'PRDM16', 'PRDM2', 'PRDM4', 'PRDM5', 'PRDM6', 'PRDM7', 'PRDM8', 'PRDM9', 'PRKAA1', 'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKAG3', 'PRKCA', 'PRKCB', 'PRKCD', 'PRKDC', 'PRMT1', 'PRMT2', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9', 'PRPF31', 'PRR14', 'PSIP1', 'PTBP1', 'PTBP1', 'PUF60', 'RAD51', 'RAD54B', 'RAD54L', 'RAD54L2', 'RAG1', 'RAG2', 'RAI1', 'RARA', 'RB1', 'RBBP4', 'RBBP5', 'RBBP7', 'RBFOX1', 'RBM11', 'RBM15', 'RBM15B', 'RBM17', 'RBM24', 'RBM25', 'RBM4', 'RBM5', 'RBM7', 'RBM8A', 'RBMY1A1', 'RBX1', 'RCC1', 'RCOR1', 'RCOR3', 'REST', 'RFOX1', 'RING1', 'RLIM', 'RMI1', 'RNF168', 'RNF2', 'RNF20', 'RNF40', 'RNF8', 'RNPS1', 'RPS6KA3', 'RPS6KA4', 'RPS6KA5', 'RPUSD3', 'RRP8', 'RSF1', 'RSRC1', 'RUVBL1', 'RUVBL2', 'RYBP', 'SAFB', 'SAP130', 'SAP18', 'SAP25', 'SAP30', 'SAP30L', 'SATB1', 'SATB2', 'SCMH1', 'SCML2', 'SCML4', 'SENP1', 'SENP3', 'SET', 'SETD1A', 'SETD1B', 'SETD2', 'SETD3', 'SETD5', 'SETD6', 'SETD7', 'SETD8', 'SETDB1', 'SETDB2', 'SETMAR', 'SF3B1', 'SF3B3', 'SFMBT1', 'SFMBT2', 'SFPQ', 'SFSWAP', 'SHPRH', 'SIN3A', 'SIN3B', 'SIRT1', 'SIRT2', 'SIRT6', 'SIRT7', 'SKP1', 'SLU7', 'SMARCA1', 'SMARCA2', 'SMARCA4', 'SMARCA5', 'SMARCAD1', 'SMARCAL1', 'SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2', 'SMARCD3', 'SMARCE1', 'SMEK1', 'SMEK2', 'SMYD1', 'SMYD2', 'SMYD3', 'SMYD4', 'SNAI2', 'SP1', 'SP100', 'SP140', 'SPEN', 'SPOP', 'SRCAP', 'SRRM4', 'SRSF1', 'SRSF10', 'SRSF12', 'SRSF3', 'SRSF6', 'SS18L1', 'SS18L2', 'SSRP1', 'STK4', 'SUDS3', 'SUPT16H', 'SUPT3H', 'SUPT6H', 'SUPT7L', 'SUV39H1', 'SUV39H2', 'SUV420H1', 'SUV420H2', 'SUZ12', 'SYNCRIP', 'TADA1', 'TADA2A', 'TADA2B', 'TADA3', 'TAF1', 'TAF10', 'TAF12', 'TAF1L', 'TAF2', 'TAF3', 'TAF4', 'TAF5', 'TAF5L', 'TAF6', 'TAF6L', 'TAF7', 'TAF8', 'TAF9', 'TAF9B', 'TBL1XR1', 'TDG', 'TDRD3', 'TDRD7', 'TDRKH', 'TET1', 'TET2', 'TET3', 'TEX10', 'TFDP1', 'TFPT', 'THRAP3', 'TLE1', 'TLE2', 'TLE4', 'TLK1', 'TLK2', 'TNP1', 'TNP2', 'TONSL', 'TOP2A', 'TOP2B', 'TP53', 'TP53BP1', 'TRA2B', 'TRIM16', 'TRIM24', 'TRIM27', 'TRIM28', 'TRIM33', 'TRRAP', 'TRUB2', 'TSSK6', 'TTK', 'TYW5', 'U2AF2', 'UBE2A', 'UBE2B', 'UBE2D1', 'UBE2D3', 'UBE2E1', 'UBE2H', 'UBE2N', 'UBE2T', 'UBN1', 'UBR2', 'UBR5', 'UBR7', 'UCHL5', 'UHRF1', 'UHRF2', 'UIMC1', 'USP11', 'USP12', 'USP15', 'USP16', 'USP17L2', 'USP21', 'USP22', 'USP3', 'USP36', 'USP44', 'USP46', 'USP49', 'USP7', 'UTY', 'VDR', 'VIRMA', 'VPS72', 'VRK1', 'WAC', 'WDR5', 'WDR77', 'WDR82', 'WHSC1', 'WHSC1L1', 'WSB2', 'WTAP', 'YAF2', 'YEATS2', 'YEATS4', 'YTHDC1', 'YWHAB', 'YWHAE', 'YWHAZ', 'YY1', 'ZBTB16', 'ZBTB33', 'ZBTB7A', 'ZBTB7C', 'ZC3H13', 'ZCWPW1', 'ZFP57', 'ZGPAT', 'ZHX1', 'ZMYM2', 'ZMYM3', 'ZMYND11', 'ZMYND8', 'ZNF217', 'ZNF516', 'ZNF532', 'ZNF541', 'ZNF592', 'ZNF687', 'ZNF711', 'ZNHIT1', 'ZRANB3', 'ZZZ3'],
    'm6a_writers'      : ['METTL3', 'METTL14', 'METTL16', 'KIAA1429','RBM15', 'WTAP'],
    'm6a_erasers'      : ['FTO', 'ALKBH5'],
    # 'm6a_readers'      : ['YTHS', 'EIF3', 'HNRNPC', 'HNRNPA2B1', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3'],
    'm6a_re_wr_er'     : ['METTL3', 'METTL14',  'KIAA1429','RBM15', 'WTAP', 'FTO', 'ALKBH5', 'YTHS', 'EIF3', 'HNRNPC', 'HNRNPA2B1', 'YTHDF1', 'YTHDF2', 'YTHDC1', 'YTHDC2', 'TYSND1', 'SND1', 'PRRC2A', 'LRPPRC', 'FMR1','FMR1NB', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3'],
    'PRC2'             : ['EZH1', 'EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7', 'JARID2', 'PCGF1', 'PCGF2', 'RING1', 'BMI1'],
    'Freya'            : ['ALKBH5', 'CYP51A1', 'DHCR7', 'DHCR24', 'EBP', 'FDFT1', 'FDPS', 'FTO', 'GGPS1', 'HNRNPC', 'HSD17B7', 'IDI1', 'IGF2BP2', 'LDLR', 'LSS', 'METTL3', 'METTL14', 'MSMO1', 'MVD', 'MVK', 'NSDHL', 'PMVK', 'SC5D', 'SQLE', 'YTHDF1', 'YTHDF2'],
    'Laura'            : ['NAMPT', 'NAPRT', 'IDO', 'DHFR', 'NMNAT1', 'NMNAT2', 'NMNAT3', 'QPRT', 'MAT2A', 'MTAP', 'WTAP', 'E2F1', 'NRK', 'TDO', 'NADSYN', 'MYC', 'BRCA1', 'PCNA', 'RAD51', 'CCNA2', 'CCNE1'],
    'Kevin'            : ['PRPF8', 'SRRM1', 'SRRM2', 'ACIN1', 'RNPS1', 'CLK1', 'CLK2', 'CLK3', 'CLK4'],
    'CM-specific'      : ['Top2a', 'Mki67', 'Prc1', 'Mis18bp1', 'Kif23', 'Kif15', 'Neil3', 'Cenpe', 'Knl1', 'Iqgap3', 'Aspm', 'Racgap1', 'Foxm1', 'Lockd', 'Esco2', 'Cdca8', 'Anln', 'Cdca2', 'Diaph3', 'Ckap2l', 'Birc5', 'Kif11', 'Sgo1', 'Cdca3', 'Knstrn', 'Uhrf1', 'Ttk', 'Cdk1', 'Ankle1', 'Ncapg', 'Spag5', 'Tpx2', 'E2f8', 'Kif4', 'Kif20b', 'Shcbp1', 'Ncapg2', 'Nusap1', 'E2f7', 'Cit', 'Rad51ap1', 'Cenpf', 'Trim59', 'Ncapd2', 'Stmn1', 'Gm42047', 'Smc2', 'Incenp', 'Lmnb1', 'Atad2', 'Plk4', 'Cep128', 'Arl4c', 'Hmgb2', 'Sgo2a', 'Ezh2', 'Lair1', 'Cep192', 'Smc4', 'Hirip3', 'Lcp1', 'Ncapd3', 'Topbp1', 'Laptm5', 'Cdkn2c', 'Runx1', 'Nsd2', 'Hjurp', 'Fam111a', 'H2afz', 'Arhgap30', 'C1qa', 'Ptprc', 'Cenpa', 'Maf', 'F630028O10Rik', 'C1qc', 'Tmpo', 'Lyz2', 'Mrc1', 'Slbp', 'Ctsc', 'Lbr', 'Gatm', 'F13a1', 'Arhgap45', 'C1qb', 'Ccdc82', 'Nav2', 'Tm6sf1', 'Smchd1', 'Nucks1', 'G2e3', 'Apobec3', 'Tubb5', 'Myo5a', 'Dab2', 'Smc1a', 'Git2', 'Arl6ip1', 'Mhrt', 'Gm31251', 'Cpeb3', 'Vegfa', 'Atcayos', 'Prune2', 'Sorbs1', 'Sorbs2', 'Pde4d', 'Rnf207', 'Lmo7', 'Coro6', 'Ryr2', 'Ppip5k2', 'Ivns1abp', 'Myh7b', 'Rbm24', 'Smtn', 'Trim63', 'Dmd', 'D830005E20Rik', 'Rbm20', 'D830024N08Rik', 'Nav2', 'Ppargc1a', 'Ccdc141', 'Pde4dip', 'Dmpk', 'Tnnt2', 'Pde7a', 'Cacna1c', 'Slc4a3', 'Clasp1', 'Gja1', 'Mybpc3', 'Nexn', 'Pcdh7', 'Gm36827', 'Lgals4', 'Obscn', 'Cacnb2', 'Camk2d', 'Palld', 'Neat1', 'Ggnbp1', 'Ttn', 'Tacc2', 'Mapt', 'Pfkfb2', 'Nnt', 'Slc27a1', 'Trim7', 'Ank3', 'Carns1', 'Clip1', 'Ralgapa2', 'Vldlr', 'Ppargc1b', 'Alpk3', '5430431A17Rik', 'Ank2', 'Acacb', 'Myom2', 'Art1', 'Lrrc2', 'Cdh2', 'Enah', 'Dtna', 'Pkp2', 'Lrrfip2', 'Myzap', 'Dot1l', 'Kidins220', 'Tbc1d4', 'Magi2', 'Mlxipl', 'Gpcpd1', 'Hk2', 'Corin', 'Ctnna3', 'Kbtbd12', 'Asph', 'Mlip', 'Retreg1', 'Fblim1', 'Pacsin3', 'Mical3', 'Speg', 'Ldb3', 'Tango2', 'Nebl', 'Agl', 'Grip2', 'Lmod2', 'Ppip5k1', 'Phkg1', 'Cenpa', 'Fhl2', 'Car14', 'Lrtm1', 'D830005E20Rik', 'Rnf207', 'Mhrt', 'Gm31251', 'Atcayos', 'Rbm24', 'Cpeb3', 'Cacnb2', 'Pde7a', 'Myh7b', 'D830024N08Rik', 'Cacna1c', 'Trim63', 'Ppargc1a', 'Mlip', 'Rbm20', 'Pfkfb2', 'Ctnna3', 'Ppip5k2', 'Pkp2', 'Lgals4', 'Vegfa', 'Trim55', 'Ccdc141', 'Mlxipl', 'Coro6', 'Nav2', 'Pde4d', 'Ppp1r3a', 'Sgcd', 'Gm47101', 'Mfn1', 'Nabp1', 'Ank3', 'Dtna', 'Alpk3', 'Ggnbp1', 'Nnt', 'Gm36827', 'Sorbs2', 'Slc8a1', 'Sox6', 'Lmo7', 'Palld', 'Plin4', 'Nexn', 'Lmod2', 'Dmd', 'Gja1', 'Corin', 'Akap6', 'Smtn', 'Slc4a3', 'Art1', 'Carns1', 'Sorbs1', 'Grb14', 'Speg', 'Trim7', 'Prune2', 'Kcnip2', 'Cdh2', 'Ryr2', 'Gm47547', 'Phkg1', '3222401L13Rik', 'Enah', 'Vldlr', 'Lrrc10', 'C130080G10Rik', 'Nebl', 'Pcdh7', 'Agl', 'Dmpk', 'Myzap', 'Clip1', 'Tnnt2', 'Pde4dip', 'Rrad', 'Lrrc2', 'Ppip5k1', 'Neat1', 'Pcgf5', 'Nceh1', 'Obscn', 'Ndufaf4', 'Svil', 'Ank2', 'Tacc2', 'Clasp1', 'Ivns1abp', 'Ttn', 'Asph', '2010111I01Rik', 'Malat1', 'Camk2d', 'Mybpc3', 'Tnnc1', 'Fhl2', 'Ptgds', 'Mb', 'Myl2', 'Actc1', 'Fabp3', 'Cox6a2', 'Slc25a4', 'Myl3', 'Ckm', 'Tpm1', 'Myh6', 'Tnni3', 'Cox7a1', 'Hspb7', 'Des', 'Atp5a1', 'Atp5b', 'Atp5g3', 'Ckmt2', 'Oxct1', 'Cox6c', 'Cryab', 'Cox5a', 'Atp5g1', 'Atp5k', 'Pgam2', 'Cox4i1', 'Ndufa4', 'Chchd10', 'Tcap', 'Atp2a2', 'Pln', 'Cox7c', 'Ttn', 'Tuba4a', 'Cox8b', 'Uqcrq', 'Idh2', 'Ndufa5', 'Uqcr11', 'Myh7', 'Atp5h', 'Actn2', 'Atp5f1', 'Ldhb', 'S100a1', 'Csrp3', 'Aldoa', 'Atp5e', 'Atp5o.1', 'Uqcrfs1', 'Mdh2', 'Sod2', 'Ndufb9', 'Ankrd1', 'Slc25a3', 'Ndufc1', 'Pdha1', 'Mdh1', 'Ech1', 'Cox7b', 'Ndufa13', 'Atp5c1', 'Hrc', 'Atp5j2', 'Xirp2', 'Ndufs2', 'Uqcrh', 'Uqcrb', '2010107E04Rik', 'Atp5j', 'Ndufa2', 'Ndufs6', 'Uqcrc1', 'Acadl', 'Etfb', 'Uqcrc2', 'Usmg5', 'Ndufb8', 'Ndufb7', 'Ndufa1', 'Got1', 'Srl', 'Aco2', 'Dsp', 'Cox6b1', 'Cox5b', 'Cyc1', 'Uqcr10', 'Jph2', 'Acaa2', 'Acadvl', 'Ndufs7', 'Eno3', 'Ndufab1', 'Etfa', 'Trp53inp2', 'Cmya5', 'Zfp106', 'Acta1', 'Dynll2', 'Myl4', 'Myl7', 'Sln', 'Dkk3', 'Stard10', 'Mybphl', 'Sbk3', 'Bmp10', 'Myl1', 'Clu', 'Nppa', 'Myl9', 'Kcnj3', 'Pam', 'Tbx5', 'Gpx3', 'Nudt4', 'Smpx', 'Atp2a2', 'Cryab', 'Ankrd1', 'Chchd10', 'Tnni3', 'Kcnk3', 'Cox6a2', 'Actc1', 'Myh6', 'Cox7c', 'Tpm1', 'Ndufa1', 'Cox6c', 'Cox7a1', 'Atp5b', 'Ndufa4', 'Cox8b', 'Mb', 'Slc25a4', 'Uqcr11', 'Atp5g1', 'Mdh1', 'Uqcrh', 'Atp5a1', 'Csrp3', 'Cycs', 'Atp5e', 'Mtus2', 'Tmod1', 'S100a1', 'Atp5j', 'Cox5a', 'Des', 'Chchd2', 'Cox4i1', 'Atp5g3', 'Cox6b1', 'Mylk3', 'Tcap', 'Atp5k', 'Uqcrq', 'Atp5f1', 'Corin', 'Cox7b', 'Mdh2', 'Fndc5', 'Aldoa', 'Ehd4', 'Atp5d', 'Atp5h', 'Doc2g', 'Uqcrb', 'Uqcrfs1', 'Nppb', 'Hspb7', 'Uqcr10', 'Atp5l', 'Srl', 'Ndufa5', 'Myoz2', 'Angpt1', 'Ndufc1', 'Ndufa11', 'Ndufb9', 'Map1lc3a', 'Ndufa13', 'Gapdh', 'Cox5b', 'Ndufb4', '2010107E04Rik', 'Chrm2', 'Ndufb11', 'Vdac1', 'Atp5j2', 'Cox7a2', 'Atp5c1', 'Eif1', 'Ndufb10', 'Aes', 'Usmg5', 'Ctgf', 'Slc8a1'],
    'PRC2_consistent'  : ['IGF2BP2', 'PLEK', 'KIF21A', 'ID1', 'AFDN', 'SPATS2L', 'CTBP2', 'LMNA', 'RETN', 'BIN1', 'PTK2', 'CST7', 'AHNAK', 'ANXA5', 'CKAP4', 'GOLM1', 'MAT1A', 'SHTN1', 'KCTD12', 'PRKAR2B', 'AP3B2', 'ANXA1', 'HLA-B', 'SMARCA1', 'SPART', 'YPEL5', 'GLRX', 'ANO6', 'FNBP1', 'ACSF2', 'FLNB', 'ALOX5AP', 'SGSH', 'CD38', 'TMEM63A', 'LGALS9', 'ARHGAP25', 'GIMAP2', 'CD2', 'CAMK4', 'CD1A', 'LGALS3BP', 'CD28', 'GNA15', 'UBA7', 'CD1C', 'FYB1', 'TRAC', 'CD1E', 'MAGEA4'],
    'm6a_story'        : ['NFKB', 'IKZF3', 'IRF4', 'MTAP', 'MYC', 'FDFT1', 'BCL2', 'PTEN', 'ASB2', 'RARA', 'FTO', 'METTL3', 'ALKBH', 'YTHDF', 'HNRNP', 'HNRNPC', 'PTEN', 'SP1', 'CEBPA', 'DHCR7', 'HMGCS1', 'SOCS1', 'SOCS3', 'NOTCH1', 'MSI2', 'BRD4', ],
    'Cristina'         : ['KDM6B', 'SPI1', 'PU1', 'CD44', 'CEBPA', 'CEBPB', 'CEBPD', 'CEBPE', 'CEBPG', 'CEBPZ', 'LTF', 'TK1', 'DNMT1'],
    'Cristina_extra'   : ['XIST', 'MLL1', 'CHD4', 'PU1', 'CEBPA', 'CEBPB', 'LSD1', 'KDM6A', 'LTF', 'SMARCA4', 'ARID1A', 'CEBPD', 'HDAC2', 'MIR21', 'CD44', 'EP300', 'KDM5B', 'TET1', 'KDM6B', 'TET2', 'lncRNA', 'SMARCB1', 'SUZ12', 'DNMT1', 'KDM4A', 'RAD21', 'NIPBL', 'HDAC3', 'CTCF', 'TK1', 'KDM1A', 'CBP', 'BRG1', 'CEBPE', 'BAF47', 'KMT2A', 'p300', 'EZH2', 'EED', 'CEBPZ', 'CREBBP', 'DNMT3A', 'DNMT3B', 'HDAC1', 'TET3', 'SPI1', 'CEBPG', 'BRD4'],
    'IGF2BP2_targets'  : ['CD6', 'STAT3', 'ATP6V1A', 'HMGA1', 'CCL20', 'MYC', 'GLUT1', 'NRAS', 'FEN1', 'EPHA2', 'SIRT1', 'SOX2', 'TNS1', 'CCND1', 'MEIS2', 'YAP1', 'SNPRD1', 'NOTCH1', 'SLC1A5', 'GPT2', 'LAMB2', 'GATA6', 'PDGFA', 'HK2', 'VEGFA', 'BCL2', 'NANOG', 'OCT4', 'ABCB1', 'DDX21', 'PRMT5'],
    'Tibo'             : ["PSMG1", "PSMG2", "PSMG3", "PSMG4", "POMP", "PSMB8", "PSMB10", "PSMB9", "HSPB1", "PSMA4", "PSMA6", "PSMA5", "PSMD4", "PSME1", "PSME2", "PSMD10", "PSME3", "PSMF1", "RAD23A", "PSMB5", "ADRM1", "PSMC3", "PSMB3", "PSMD8", "PSMB6", "PSMA7", "PSMD13", "PSMA3", "PSMB1", "PSMA2", "PSMB7", "PSMA1", "PSMB4", "PSMC1", "PSMC2", "PSMD3", "PSMD7", "PSMB2", "VCP", "TXNL1", "UBQLN1", "PSMD9", "ZFAND2A", "PSMD14", "PSMD1", "RAD23B", "PSMC5", "PSMC4", "PSMD2", "UBR1", "PSMC6", "PSMD6", "PSMD5", "UBE3A", "USP14", "PSMD11", "PSME4", "PSMD12", "PSMA8", "UCHL5", "PAAF1", "UBE3C", "ECPAS", "PRICKLE1"],
    'COSMIC_TALL'      : ['ABL1', 'AFF3', 'BCL11B', 'CCNC', 'CNOT3', 'DNM2', 'FBXW7', 'IRS4', 'KDM6A', 'LCK', 'LEF1', 'LMO1', 'LMO2', 'LYL1', 'NOTCH1', 'NUP214', 'OLIG2', 'PICALM', 'PTPRC', 'RAP1GDS1', 'RPL10', 'RPL5', 'RUNX1', 'SET', 'STIL', 'TAL2', 'TLX1', 'TLX3', 'TRA', 'TRB'],
    'SPI1_TRRUST_mmus' : ['ACP5', 'AIF1', 'BCL2L11', 'CCL5', 'CD1D1', 'CSF1R', 'CSF1R', 'DCSTAMP', 'ELANE', 'ELF1', 'FCER1G', 'FES', 'FLT3', 'H2-AB1', 'IGH', 'IGHM', 'IGL', 'IL1B', 'IL7R', 'IL9', 'IRF8', 'ITGA2B', 'ITGB3', 'MED1', 'MEF2C', 'MMP13', 'MRC1', 'NCF1', 'NCF4', 'NFKB1', 'OPRM1', 'PIRB', 'PTGDS', 'PTGDS', 'PTPRC', 'PTPRO', 'SCN8A', 'TAL1', 'TAL2', 'TEC', 'TLR4', 'TNF'],
    'SPI1_TRRUST_hsap' : ['ACP5', 'ALOX15', 'BCL6', 'BPI', 'BTK', 'CCL2', 'CCL5', 'CD163', 'CD22', 'CD40', 'CD68', 'CLEC4G', 'CSF1', 'CSF2RA', 'CSF3R', 'CSF3R', 'CTSG', 'CTSK', 'CTSS', 'CXCR1', 'CYBB', 'CYBB', 'DAPK2', 'ELANE', 'ERAP2', 'FCER1A', 'FCER1A', 'FCGR1A', 'FES', 'FLI1', 'GATA1', 'GATA1', 'HCK', 'IFIT3', 'JCHAIN', 'IL12B', 'IL18', 'IL1B', 'IL1B', 'IL5', 'ITGA2B', 'ITGAM', 'ITGAX', 'ITGAX', 'ITGB2', 'ITGB2', 'MACROD1', 'MAPK1', 'MME', 'MNDA', 'MS4A1', 'MSR1', 'NCF2', 'NCF2', 'NSFL1C', 'P2RY10', 'PARG', 'PEBP1', 'PRG2', 'PRTN3', 'PTGIR', 'RNASE2', 'SCARB1', 'SCARB2', 'STAT3', 'TLR4', 'TLR4', 'TNF', 'TNFRSF11A', 'WAS', 'ZNF300'],
    'Polonen_SPI1_ass' : ['CD68', 'LYZ', 'SPI1', 'CSF3R', 'CYBB', 'TYROBP', 'PLEK', 'CD74', 'SIRPA', 'S100A9', 'NCF2', 'HLA-DRA', 'S100A8', 'THBS1', 'VCAN', 'SERPINA1', 'LILRA2', 'MPO', 'SRGN', 'LILRB2', 'CEBPD', 'TREM1', 'ITGAX', 'FCAR', 'FCN1', 'FCGR2A', 'CEBPA', 'MPEG1', 'AC090559.1', 'MNDA', 'LYN', 'RAB31', 'CFD', 'HLA-DPA1', 'HCK', 'IRAK3', 'LILRA5', 'ITGAM', 'AQP9', 'CTSZ', 'LRRC25', 'LTF', 'CLEC7A', 'IFI30', 'CSF1R', 'PLBD1', 'FPR1', 'HLA-DRB1', 'RBM47', 'S100A12', 'G0S2', 'NFAM1', 'TGFBI', 'TLR4', 'MYO1F', 'HLA-DPB1', 'IL6R', 'C5AR1', 'MS4A6A', 'TLR2', 'TNFRSF1B', 'CD300E', 'CSTA', 'LRP1', 'PTAFR', 'CHST15', 'SEMA6B', 'CD14', 'SIRPB1', 'HK3', 'CD86', 'CLEC12A', 'ZBTB7B', 'ZNF516', 'IL1RN', 'THBD', 'CD300LF', 'SLC11A1', 'FCER1G', 'ANPEP', 'PRAM1', 'HMOX1', 'ELANE', 'SLC7A7', 'CSF2RB', 'NOD2', 'NCF4', 'CXCL8', 'CCR1', 'ICAM1', 'C19orf38', 'CXCL16', 'LRRK1', 'LILRB4', 'FGL2', 'THEMIS2', 'AC005840.1', 'PADI2', 'SIGLEC14', 'AHR', 'IL13RA1', 'NLRP3', 'TYMP', 'PYGL', 'CD300LB', 'FAM49A', 'SAMHD1', 'TRIB1', 'LILRB3', 'TLR8', 'OGFRL1', 'DEFA3', 'CTSS', 'ALDH3B1', 'RIPK2', 'BASP1', 'LRRK2', 'METRNL', 'ADM', 'DEFA1', 'FBP1', 'IGKC', 'BTK', 'SERPINB2', 'PDE4B', 'CD33', 'MEFV', 'PPBP', 'DEFA1B', 'TNFSF13', 'LGALS1', 'LILRA1', 'P2RY13', 'ALOX5', 'SIGLEC5', 'ACSL1', 'CIITA', 'MMP9', 'UPP1', 'CYTIP', 'GPR183', 'NEK6', 'LILRA6', 'OSCAR', 'LILRB1', 'HLX', 'PLXNB2', 'CD163', 'PHACTR1', 'NCF1', 'APP', 'EPB41L3', 'RNASE2', 'FGD2', 'AZU1', 'RASGRP4', 'GSN', 'STAB1', 'CYP1B1', 'CSF2RA.1', 'WDFY4', 'CD302', 'CST3', 'LGALS3', 'CDA', 'SERPINB6', 'NCF1C', 'RMRP', 'DYSF', 'NAPSB', 'NID1', 'COTL1', 'SMIM25', 'TMEM176B', 'IL10RA', 'PRTN3', 'SLC15A3', 'SLCO3A1', 'ATP8B4', 'CD93', 'IRF5', 'APOBEC3A', 'IL1B', 'LY86', 'ADAM15', 'FCGR3A', 'PDK4', 'BCL3', 'STX11', 'MS4A7', 'SECTM1', 'LCN2', 'IRF4', 'VDR', 'BPI', 'KLF4', 'CFP', 'GASK1B', 'LFNG', 'PIK3AP1', 'HLA-DMB', 'FOSL2', 'EREG', 'IGHG1', 'ST14', 'SULF2', 'C3AR1', 'HLA-DQB1', 'TMEM176A', 'LPCAT2', 'PTPRJ', 'ENTPD1', 'TENT5A', 'SNTB1', 'NRGN', 'CD180', 'CD36', 'MMP8', 'LRG1', 'HLA-DMA', 'IGFBP7', 'FES', 'PLD4', 'TFEC', 'KLF9', 'JDP2', 'ODF3B', 'SNX9', 'P2RY2', 'GRN', 'ELL2', '1-Mar', 'SIGLEC10', 'MEF2C', 'AC020916.1', 'ACPP', 'SPOPL', 'PREX1', 'HCAR3', 'PTGES', 'IGHA1', 'NACC2', 'SIRPB2', 'LMNA', 'STX3', 'TLR1', 'MMP25', 'KCNE3', 'PARVB', 'PADI4', 'NLRP12', 'VSIR', 'SCPEP1', 'PTGS1', 'QPCT', 'ADAMTSL4', 'TMCC3', 'MTSS1', 'SIGLEC9', 'KLF11', 'KCNQ1', 'JAML', 'MYADM', 'NUAK2', 'S100A11', 'CES1', 'GPR132', 'MYOF', 'PROK2', 'HLA-DQA1', 'CEACAM8', 'ZNF467', 'S100A6', 'PLXDC2', 'TES', 'FFAR2', 'C1orf162', 'FCGR1A', 'TNFSF13B', 'KYNU', 'ADAM28', 'ARHGEF40', 'LRP3', 'FGD4', 'MS4A3', 'SLC24A4', 'NCF1B', 'KIF13A', 'PTGER2', 'OLR1', 'S100A4', 'UNC93B1', 'HAVCR2', 'HP', 'B3GNT5', 'PTGS2', 'CAVIN2', 'MGST1', 'GAPT', 'MARCKS', 'TNFRSF10D', 'HRH2', 'MGAM', 'F5', 'F13A1', 'STARD8', 'HLA-F', 'IRAK2', 'MTMR11', 'PLAUR', 'ZEB2', 'IQSEC2', 'CD101', 'IGLC2', 'ASGR2', 'MAP3K5', 'AC139495.3', 'CTSL', 'MAP3K8', 'GBP2', 'KLF3'],
    'SG'               : ['MEN1', 'KMT2A', 'HOXA9', 'NPM1', 'MEIS1', 'LMO2', 'PTEN', 'CD2', 'LCK', 'EZH2', 'SUZ12', 'EED', 'LMO2'],
    'Master'           : [
        #PRC1/2
        'RING1', 'BMI1', 'PCGF1', 'PCGF2', 'EED', 'SUZ12', 'EZH1', 'EZH2', 'RBBP4', 'RBBP7', 'JARID2',
        #Chromatin stuff
        'CTCF', 'ARID1A', 'PBRM1', 'ARID1B', 'SMARCA2', 'SMARCA4', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2', 'SMARCD3', 'ACTL6A', 'ACTL6B', 'SMARCB1', 'SA1', 'SA2', 'SMC1A', 'SMC3', 'RAD21',
        #NOTCH1 upregulated
        'MYC', 'HEY1', 'NFKB1', 'NFKB2', 'CCND1', 'CDKN1A', 'CDKN1B', 'NOTCH1', #Proteins of interest in lab below
        #Laura
        'NAMPT', 'NAPRT', 'NNMAT1', 'NNMAT2', 'NNMAT3', 'IDO1', 'QPRT', 'DHFR', 'MAT2A' 'HES1', 'BRD4', 'DTX1', 'FOXO3', 'SESN1', 'U2AF1', 'U2AF2',
        #Involved in monocyte interactions with vascular endothelial repair. From Wtikowski, 2020
        'PECAM1', 'CD44', 'ITGA4', 'CX3CR1', 'TNFSF10', 'CSF1R ',
        #m6a proteins
        'METTL3', 'METTL14', 'METTL16', 'WTAP', 'RBM15', 'RBM15B' 'KIAA429', 'ZC3HI3', 'CBLL1', 'VIRMA', 'FTO', 'ALKBH5', 'HNRNPC', 'YTHDF1', 'YTHDF2', 'YTHDF3', 'YTHDC1', 'YTHDC2', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'HNRNPA2B1',
        #Frequently mutated in AML
        'DNMT3A', 'TET2', 'ASXL1','TP53', 'KRAS', 'IDH',
        #Deaminase stuff
        'APOBEC1', 'APOBEC2', 'APOBEC3A','APOBEC3B','APOBEC3C','APOBEC3D','APOBEC3F','APOBEC3G','APOBEC3H','APOBEC4', 'AID',
        #Mentioned alongside CCR9 in Mansour's paper, expressed in MOLT4 but not healthy cells
        'SD06', 'OR10R2', 'SKAJ', 'DNAH17', 'LZTFL1', 'TCF19', 'MAD2L1', 'UBASH3A', 'LCK', 'ITGAE',
        #Gene is in depmap as having the lowest median score among all cancers for T-ALL
        'AHCY', 'CHERP','CTPS1','DCDC2','H2BC11','NPTXR','RHOH','SUCLG1','TAS1R2','TONSL','TPT1','ZAP70'
        #Alterations in these genes across 15 subtypes of T-ALL
        'BCL11B', 'KMT2A', 'MLLT10','HOXA9','NKX2','CUL1','CDKN2A','PHF6','LEF1','PSIP1','NOL4L','KMT2E','KAT6A','DDX39B','MYL1','CD99','E2F1','GEMIN8','FBXW7',
        #Miscellaneous
        'FBXW7', 'SF3B1', 'BRD9', 'NT5C2', 'FOXO3', 'USP7', 'BCL2L11', 'MYB', 'RUNX1', 'IKZF3', 'KAT2B', 'PTPN4', 'CYLB', 'ITGB1', 'CCND3', 'PSM3IP', 'CREBBP', 'PSMB9', 'SRSF2', 'SRSF3', 'SRSF6', 'SRSF7',
        'SRSF11', 'NUP85', 'PSMD4', 'NUP85', 'UPF1', 'KDM6B', 'IDH2','INTS3','SPI1', 'TCF7','DHCR7','HMGCS1','WT1','FLT3','KIT','MDA5','PKR','NCOR2''FMR1','MATR3','EP300','LMO1','IRF4','KDM4C',
        'JQ1', 'MXD1', 'MTAP', 'WTAP', 'CHEK1', 'RBM39', 'CCR9', 'CCR7', 'CD19', 'CDK6', 'BRD1', 'BRD2', 'BRD3', 'BRD4', 'APCDD1', "FOXO1", 'PRPF39','MSI2', 'NFKIA','SERPINB1','CD69','CCL25', 'DCAF15',
        'TNFSF9', 'TNFRSF9', 'CD38', 'TRBC1', 'TRBC2', 'ZBTB7B', 'RUNX3','ZBP1', 'FBXW7', 'KMT2D', 'IGF2','TAL1','PTEN', 'IDH1','PIK3CA','PIK3R1','STAT5B','JAK1','JAK2','DNM2','IL7R','NGN3','CASC3', 'IDO',
        'MEF2C', 'TCF3','PBX1','HLF','ITGAM','AURKB','TLX3','STAG2','SOX11','XBP1','ID1','ID3','XBP1','CD34','LMO2', 'JAK3', 'RIGI','ADAR','ADARB1','ADARB2','IGHV','FOXM1','AURKA', 'BCL2', 'FOXP3',
        'RELA','TERT','ARNT','HDAC1','HDAC2','HDAC6','HDAC10','NPM1','KDM1B'
        ],
	}


# A set of genes is selected based on the string used as input.
# First, the string is tested in the custom set of gene sets in gene_sets dictionary above.
# If that fails, it will use the string to search for a public gene set on Msigdb
# If that fails, it will default to interpreting the input string as a set with a single gene name (the input string)

def KTC_GetGeneSet(name_gene_set, db_version='2024.1.Hs'):
    try:
        # If the input is a list, capitalize all items and use it as the gene set
        if isinstance(name_gene_set, list):
            gene_set = [gene.upper() for gene in name_gene_set]
            print(f'\n -- KTC_GetGeneSet : Using provided list of genes as the gene set: {gene_set}\n')
            return gene_set

        # If the input is a string, proceed with existing logic
        # First, try to find the gene set in the `gene_sets` dictionary
        gene_set = gene_sets[name_gene_set]
        print(f'\n -- KTC_GetGeneSet : Using custom set of genes from gene_sets dictionary: {name_gene_set}\n')
    except KeyError:
        # If not found, try to fetch it from msigdb, prioritizing 'HALLMARK' over others
        if db_version == '2024.1.Hs':
            categories = ['h.all', 'c1.all', 'c2.all', 'c3.all', 'c4.all', 'c5.all', 'c6.all', 'c7.all', 'c8.all', 'msigdb']
        elif db_version == '2024.1.Mm':
            categories = ['mh.all', 'm1.all', 'm2.all', 'm3.all', 'm5.all', 'm8.all', 'msigdb']
        print(f'\n -- KTC_GetGeneSet : {name_gene_set} Searching {db_version} for gene set {name_gene_set} using Msigdb...\n')
        msigdb = Msigdb()
        gmt_cache = {}
        gene_set = None  # Initialize gene_set in case no match is found

        for category in categories:
            try:
                if category not in gmt_cache:
                    gmt_cache[category] = msigdb.get_gmt(category=category, dbver=db_version)
                gmt = gmt_cache[category]

                if gmt and name_gene_set in gmt:
                    print(f'\n -- KTC_GetGeneSet : Using gene set from msigdb: {name_gene_set} from {category}\n')
                    gene_set = gmt[name_gene_set]
                    break  # Exit loop once a match is found
            except Exception as e:
                print(f"Error fetching or searching category {category}: {e}\n")

        if gene_set is None:
            # If no match is found, use string as a gene name for a set with a single gene
            print(f'\n -- KTC_GetGeneSet : Gene set: {name_gene_set} was not a list of genes, was not found in gene_sets, or using Msigdb. Defaulting to interpreting {name_gene_set} as a set with a single gene\n')
            gene_set = [name_gene_set.upper()]
    return gene_set


#%% ===========================================================================
# KTC_orthologue
# =============================================================================
from pybiomart import Server

def KTC_orthologue(gene_list, organism_origin='hsapiens', organism_target='mmusculus'):
    print(f'\n -- KTC_orthologue : Converting genes from {organism_origin} to {organism_target} orthologues...')

    # Connect to Ensembl BioMart
    server = Server(host='http://www.ensembl.org')
    origin_dataset = f'{organism_origin}_gene_ensembl'
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets[origin_dataset]

    # Attributes for query
    origin_attr = 'external_gene_name'
    target_attr = f'{organism_target}_homolog_associated_gene_name'

    # Fetch orthologs
    results = dataset.query(attributes=[origin_attr, target_attr])

    # Inspect column names (to handle variations)
    print("Columns in results:", results.columns)

    # Match the origin column explicitly for the 'origin' organism
    origin_col = None
    target_col = None

    # Check for the column names based on the organisms
    if organism_origin == 'hsapiens' and organism_target == 'mmusculus':
        origin_col = 'Gene name'
        target_col = 'Mouse gene name'
    elif organism_origin == 'mmusculus' and organism_target == 'hsapiens':
        origin_col = 'Gene name'
        target_col = 'Human gene name'
    else:
        raise ValueError(f"Unsupported organisms: {organism_origin} to {organism_target}")

    # Make sure we found the correct columns
    if not origin_col or not target_col:
        raise ValueError(f"Could not identify origin and target columns. Found columns: {results.columns}")

    print(f'Origin column: {origin_col} | Target column: {target_col}')

    # Filter and convert gene names
    converted_genes = results[results[origin_col].isin(gene_list)][target_col].dropna().tolist()
    
    return converted_genes

#%% ===========================================================================
# KTC_Splice_Mapper
# =============================================================================

def KTC_preParse_gtf(gtf_file):
    import pandas as pd
    print(' -- KTC_preParse_gtf : Pre-parsing GTF...')

    # Load GTF with standard column names
    gtf = pd.read_csv(
        gtf_file,
        sep='\t',
        comment='#',
        header=None,
        names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
        dtype={0: str}
    )

    # Extract gene_name, transcript_id, and exon_number from the attribute field
    gtf['gene_name'] = gtf['attribute'].str.extract(r'gene_name "([^"]+)"')
    gtf['transcript_id'] = gtf['attribute'].str.extract(r'transcript_id "([^"]+)"')
    gtf['exon_number'] = gtf['attribute'].str.extract(r'exon_number "(\d+)"')

    print(' -- KTC_preParse_gtf : Done')
    return gtf

def KTC_splice_map(gtf_file, rmats_file, event_id, splicing_event_type):
    import pandas as pd
    import matplotlib.pyplot as plt
    print(' -- KTC_SpliceMap : Reading gtf...')

    # Check if gtf_file is already a DataFrame (i.e., a pre-parsed GTF)
    if isinstance(gtf_file, pd.DataFrame):
        gtf = gtf_file
        print(' -- KTC_SpliceMap : Using provided pre-parsed GTF.')
    else:
        # Read GTF file if it's not a pre-parsed DataFrame
        gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                          names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                          dtype={0: int})
    
    print(' -- KTC_SpliceMap : Reading rMATS...')
    # Read rMATS file
    rmats = pd.read_csv(rmats_file, sep='\t')
    rmats = rmats[rmats['Splicing Event'] == splicing_event_type]

    # print("Event ID:", event_id)
    
    # Use the first unnamed column (index) for filtering
    rmats_index = rmats.iloc[:, 0]  # First unnamed column as index
    rmats_index = rmats['ID']
    
    # Filter the rMATS data based on the event ID from the unnamed index column
    filtered_event = rmats[rmats_index == event_id]

    if filtered_event.empty:
        print(f' -- KTC_SpliceMap : No event found with ID {event_id}')
        # print(f"")
        return
    
    # Get the first matching event
    event = filtered_event.iloc[0]
    gene_id = event['geneSymbol']
    # print(gene_id)
    
    # Check if transcript_id exists in the rMATS data
    if 'transcript_id' in filtered_event.columns:
        transcript_id = event['transcript_id']  # If available, use the transcript_id for filtering
        print(f' -- KTC_SpliceMap : Using transcript_id: {transcript_id}')
    else:
        transcript_id = None
        print(' -- KTC_SpliceMap : No specific transcript_id in rMATS event, using default transcript version.')
    global exons

    print(' -- KTC_SpliceMap : Filtering gtf...')
    # Filter for exons of the specified gene
    if transcript_id:
        exons = gtf[(gtf['feature'] == 'exon') & (gtf['attribute'].str.contains(f'gene_name "{gene_id}"')) & 
                    (gtf['attribute'].str.contains(f'transcript_id "{transcript_id}"'))]
    else:
        exons = gtf[(gtf['feature'] == 'exon') & (gtf['attribute'].str.contains(f'gene_name "{gene_id}"'))]
    
    exons = exons.sort_values(by='start')  # Sort exons by their start position
    
    # Extract transcript version info from attributes
    exons['transcript_version'] = exons['attribute'].str.extract(r'transcript_version "(\d+)"')
    exons['gene_id'] = exons['attribute'].str.extract(r'gene_id "([^"]+)"')
    # Identify the most frequently found transcript_id
    most_common_transcript = exons['attribute'].str.extract(r'transcript_id "([^"]+)"')[0].mode()[0]
    # print(f"Most common transcript_id: {most_common_transcript}")
    
    # Filter exons to keep only the most common transcript
    exons = exons[exons['attribute'].str.contains(f'transcript_id "{most_common_transcript}"')]
    
    # Extract exon numbers, handling both human and mouse GTF formats
    exons['exon_number'] = exons['attribute'].str.extract(r'exon_number "(\d+)"')
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 2), dpi=100)
    y = 0.5  # Vertical placement of the gene track
    prev_end = None
    
    print(' -- KTC_SpliceMap : Plotting elements...')
    for _, row in exons.iterrows():
        # Use the assigned exon_number instead of extracting it from the attribute
        exon_number = row['exon_number']
        
        # Draw intron as a line
        if prev_end is not None:
            ax.plot([prev_end, row['start']], [y, y], color='black', linewidth=0.8)
        
        # Draw exon as a block
        ax.add_patch(plt.Rectangle((row['start'], y - 0.2), row['end'] - row['start'], 0.4, color='black'))
        
        # Label exons
        # Alternate exon number labels above and below
        label_y = y + 0.35 if int(exon_number) % 2 == 1 else y - 0.35
        exon_center = (row['start'] + row['end']) / 2
        ax.text(exon_center, label_y, f'{exon_number}', ha='center', va='center', fontsize=8)
        prev_end = row['end']

    # Highlight splicing event locations
    highlight_cols = ['longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE',
                      'exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
                      '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd', 'riExonStart_0base', 'riExonEnd']
    highlight_cols = [col for col in highlight_cols if col in event and pd.notna(event[col])]


    colors = plt.colormaps['Dark2'](range(len(highlight_cols)))
    legend_elements = []

    for i, col in enumerate(highlight_cols):
        if pd.notna(event[col]):
            ax.annotate('', (event[col], y - 0.4), (event[col], y - 0.75),
                        arrowprops=dict(arrowstyle='-|>', color=colors[i], lw=1.5))  # Access color using colors[i]
            legend_elements.append(plt.Line2D([0], [0], color=colors[i], lw=1.5, label=col))  # Access color using colors[i]
    
    ax.set_xlim(exons['start'].min() - 100, exons['end'].max() + 100)
    ax.set_ylim(-1, 1)
    ax.axis('off')
    fdr = event['FDR']
    inc_level_diff = event['IncLevelDifference']
    
    # Adjust the title to show event information
    plt.title(f'Gene Map: {gene_id} (Splicing Event: {splicing_event_type} {event_id})', fontsize=14)
    plt.text(0.5, 0.2, f'FDR: {fdr:.4e} | IncLevelDifference: {inc_level_diff:.4f}',
         ha='center', va='center', transform=ax.transAxes, fontsize=10, color='darkred')
    plt.legend(handles=legend_elements, title='Splicing Event Locations', fontsize=8, title_fontsize=10, loc='lower right', bbox_to_anchor=(1.22, 0.15))
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    print(' -- KTC_SpliceMap : Done')


## Example usage - parsing new gtf every time:
# KTC_splice_map('/Users/kasperthorhaugechristensen/Downloads/Homo_sapiens.GRCh38.110.chr.gtf', '/Users/kasperthorhaugechristensen/Downloads/KO1_rMATS_compiled.tsv', 80549, 'SE')

## Example usage - pre-parsing a gtf file once and feeding that dataframe to each iteration (saving time)
# preparsed_gtf = KTC_preParse_gtf('/Users/kasperthorhaugechristensen/Downloads/Mus_musculus.GRCm39.110.chr.gtf')
# KTC_splice_map(preparsed_gtf, '/Volumes/cmgg_pnlab/Kasper/Analyses/Joao/2025_TF_analysis/CD2_v_Vav_rMATS_compiled.tsv', 53351, 'SE')

#%% ===========================================================================
# KTC_Volcano
# =============================================================================

def KTC_PlotVolcano(
    df,
    l2fc_col='log2FoldChange',
    padj_col='padj',
    label_col=None,
    padj_thresh=0.05,
    l2fc_thresh=1.0,
    top_n=10,
    highlight_genes=None,
    label_highlight=True,   # 👈 NEW
    title="Volcano Plot",
    figsize=(8, 6),
    dpi=200,
    alpha=0.7,
    point_size=20,
    show_plot=True,
    save_path=None,
    ax=None,
    x_limits=None
):
    df = df.copy()

    # Ensure numeric
    df[l2fc_col] = pd.to_numeric(df[l2fc_col], errors="coerce")
    df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce")

    # Guard against zeros / negatives / NaNs in padj for -log10
    df["-log10(padj)"] = -np.log10(df[padj_col].clip(lower=np.finfo(float).tiny))

    # Default significance logic
    df["significant"] = (df[padj_col] < padj_thresh) & (df[l2fc_col].abs() > l2fc_thresh)

    # Highlight logic
    if highlight_genes is not None and label_col:
        df["highlight"] = df[label_col].isin(highlight_genes)
    else:
        df["highlight"] = df["significant"]

    # Create axis if needed
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        created_fig = True
    else:
        fig = ax.figure

    # Base (non-highlight) points
    base_df = df[~df["highlight"]]
    sns.scatterplot(
        data=base_df,
        x=l2fc_col,
        y="-log10(padj)",
        color="lightgrey",
        alpha=alpha,
        s=point_size,
        edgecolor=None,
        ax=ax
    )

    # Highlighted points
    highlight_df = df[df["highlight"]]
    if not highlight_df.empty:
        sns.scatterplot(
            data=highlight_df,
            x=l2fc_col,
            y="-log10(padj)",
            color="steelblue",
            alpha=1,
            s=point_size * 1.2,
            edgecolor="black",
            linewidth=0.5,
            ax=ax
        )

    if x_limits is not None:
        ax.set_xlim(x_limits)

    # Labels
    texts = []
    if label_col and not df.empty:
    
        if highlight_genes is not None and not highlight_df.empty:
    
            if label_highlight:
                label_df = highlight_df
            else:
                label_df = pd.DataFrame()  # 👈 no labels
    
        elif top_n > 0:
            sig_df = df[df["significant"]].copy()
            if not sig_df.empty:
                label_df = sig_df.loc[sig_df[l2fc_col].abs().nlargest(top_n).index]
            else:
                label_df = pd.DataFrame()
    
        else:
            label_df = pd.DataFrame()
    
        for _, row in label_df.iterrows():
            texts.append(
                ax.text(
                    row[l2fc_col],
                    row["-log10(padj)"],
                    str(row[label_col]),
                    fontsize=8
                )
            )
    
        if texts:
            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            )

    # Threshold lines
    ax.axhline(-np.log10(padj_thresh), color='black', linestyle='--', lw=1)
    ax.axvline(-l2fc_thresh, color='black', linestyle='--', lw=1)
    ax.axvline(l2fc_thresh, color='black', linestyle='--', lw=1)

    # Final formatting
    ax.set_title(title)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10 Adjusted p-value")
    sns.despine(ax=ax)

    # Save/show only if we created the figure
    if created_fig:
        if save_path:
            fig.savefig(save_path, bbox_inches="tight")
        if show_plot:
            plt.show()
        else:
            plt.close(fig)

    return ax


#%% ===========================================================================
# 
# =============================================================================


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gseapy as gp


def KTC_multi_enrichment_plot(
    gene_list,
    gene_sets=[
        'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
        "TRRUST_Transcription_Factors_2019",
        "GO_Biological_Process_2021",
        "KEGG_2019_Human",
        "ChEA_2016",
        "MSigDB_Hallmark_2020",
        'Reactome_2022',
        'ARCHS4_Tissues'
    ],
    organism="human",
    top_n=3,
    show_all_significant=False,
    figsize=(12, 8),
    out_dir=None,
    title=None,
    title_suffix='',
    filename_suffix='',
    palette='Set2',
    significant_only=True,
    filter_term_keywords=None,
    show_legend=True,
    padj_cutoff=0.05,
    sort_by="gene_ratio" 
):
    if not isinstance(gene_list, list) or not gene_list:
        raise ValueError("gene_list must be a non-empty list of gene symbols.")

    # Clean input list
    gene_list = sorted({g.strip().upper() for g in gene_list if isinstance(g, str) and g.strip()})
    N = len(gene_list)

    print(f"\n -- Performing batch enrichment for {len(gene_sets)} gene sets...")

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            organism=organism,
            gene_sets=gene_sets,
            no_plot=True,
            outdir=None
        )
    except Exception as e:
        print(f"  [ERROR] Batch enrichment request failed: {e}")
        return pd.DataFrame()

    if enr.results is None or enr.results.empty:
        print("⚠ No enrichment results found.")
        return pd.DataFrame()

    if 'Gene_set' not in enr.results.columns:
        enr.results['Gene_set'] = enr.results.index.get_level_values(0)

    all_results = []

    for gs in gene_sets:
        df = enr.results[enr.results['Gene_set'] == gs].copy()

        if df.empty:
            print(f"  [INFO] No results for {gs}, skipping.")
            continue

        if significant_only:
            df = df[df["Adjusted P-value"] < padj_cutoff]
            if df.empty:
                print(f"  [INFO] No significant terms (adj p < {padj_cutoff}) for {gs}, skipping.")
                continue

        df = df.sort_values("Adjusted P-value")
        if not show_all_significant:
            df = df.head(top_n)

        all_results.append(df)

    if not all_results:
        print("No enrichment results to display after filtering.")
        return pd.DataFrame()

    combined_df = pd.concat(all_results, ignore_index=True)

    # Optional keyword filtering
    if filter_term_keywords is not None:
        import re
        regex_pattern = "|".join(map(re.escape, filter_term_keywords))
        combined_df = combined_df[combined_df["Term"].str.contains(regex_pattern, case=False, na=False)]
        if combined_df.empty:
            print(f"No terms matched the filter keywords: {filter_term_keywords}")
            return pd.DataFrame()

        if significant_only:
            combined_df = combined_df[combined_df["Adjusted P-value"] < padj_cutoff]
            if combined_df.empty:
                print("No significant terms remain after keyword filtering.")
                return pd.DataFrame()

    # Compute effect size from Overlap like "12/187"
    overlap = combined_df["Overlap"].astype(str).str.split("/", expand=True)
    combined_df["k_overlap"] = pd.to_numeric(overlap[0], errors="coerce")
    combined_df["K_term"] = pd.to_numeric(overlap[1], errors="coerce")
    combined_df["gene_ratio"] = combined_df["k_overlap"] / max(N, 1)

    combined_df["Term_short"] = combined_df["Term"].astype(str).str.slice(0, 60)
    
    combined_df["y_label"] = (
        combined_df["Term_short"]
        + "  [k="
        + combined_df["k_overlap"].fillna(0).astype(int).astype(str)
        + "]"
    )


    if sort_by == "gene_ratio":
        combined_df = combined_df.sort_values(
            ["gene_ratio", "Adjusted P-value"],
            ascending=[False, True]
        )
    elif sort_by in ["padj", "Adjusted P-value"]:
        combined_df = combined_df.sort_values(
            ["Adjusted P-value", "gene_ratio"],
            ascending=[True, False]
        )
    else:
        raise ValueError("sort_by must be 'gene_ratio' or 'padj'")

    # Alphabetic ordering for deterministic colors
    # present_sets = sorted(combined_df["Gene_set"].unique())
    
    # Build a deterministic palette list with one color per gene set
    # palette_list = sns.color_palette(palette, n_colors=len(present_sets))
    # print(palette_list)


    # Dynamic figure height
    n_terms = combined_df.shape[0]
    buffer = 1
    height_per_row = 0.35
    max_height = 24
    dynamic_height = min(buffer + (n_terms * height_per_row), max_height)

    # Plot:
    #   x-axis = Adjusted P-value (log scale)
    #   color  = Gene_set
    #   size   = gene_ratio
    #   y      = "Gene_set [k=..] Term"
    plt.figure(figsize=(figsize[0], dynamic_height), dpi=200)

    ax = sns.scatterplot(
        data=combined_df,
        x="Adjusted P-value",
        y="y_label",
        hue="Gene_set",
        size="gene_ratio",
        sizes=(40, 350),
        alpha=0.9,
        edgecolor="none",
        hue_order=gene_sets,      # <--- key for reproducibility
        # palette=palette_list,        # <--- list aligned to hue_order
    )

    # Log scale for p-values so you can see range properly
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.axvline(padj_cutoff, linestyle="--", linewidth=1.5)
    ax.set_xlabel("Adjusted p-value (log scale)")
    ax.set_ylabel("")

    if title is None:
        title = f"Enrichment overview ({N} genes)"
    ax.set_title(title + title_suffix)

    if show_legend:
        # Put legend outside if it gets crowded
        ax.legend(loc="lower left", bbox_to_anchor=(1.02, 0.0), borderaxespad=0.0, frameon=True)
    else:
        if ax.legend_ is not None:
            ax.legend_.set_visible(False)
            
    handles, labels = ax.get_legend_handles_labels()

    new_labels = []
    for lab in labels:
        try:
            new_labels.append(f"{float(lab):.2f}")
        except ValueError:
            new_labels.append(lab)
    
    ax.legend(handles, new_labels,
              loc="lower left",
              bbox_to_anchor=(1.02, 0.0),
              borderaxespad=0.0,
              frameon=True)


    plt.tight_layout()

    if out_dir:
        filename = f"multi_enrichment_dotplot_padj{filename_suffix}.svg"
        path = os.path.join(out_dir, filename)
        plt.savefig(path)
        print(f"Saved plot to {path}")

    plt.show()
    return combined_df


#%% ===========================================================================
# KTC_plot_transcript_splicing
# =============================================================================

def KTC_get_canonical_transcript(gene_symbol, species='human'):
    """
    Query Ensembl REST API to get the canonical transcript for a gene.
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/{species}/{gene_symbol}?expand=1"
    headers = {"Content-Type": "application/json"}

    r = requests.get(server + ext, headers=headers)
    if not r.ok:
        print(f"❌ Error fetching data for {gene_symbol}")
        return None

    data = r.json()
    if "Transcript" in data:
        for transcript in data["Transcript"]:
            if transcript.get("is_canonical"):
                return transcript["id"]
    return None


def KTC_plot_transcript_splicing(gtf_file, rmats_file, event_id, splicing_event_type):
    print(' -- KTC_plot_transcript_splicing : Reading GTF...')

    # Load or use pre-parsed GTF DataFrame
    if isinstance(gtf_file, pd.DataFrame):
        gtf = gtf_file.copy()
        print(' -- KTC_plot_transcript_splicing : Using pre-parsed GTF')
    else:
        gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                          names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
        gtf['gene_name'] = gtf['attribute'].str.extract(r'gene_name "([^"]+)"')
        gtf['transcript_id'] = gtf['attribute'].str.extract(r'transcript_id "([^"]+)"')
        gtf['exon_number'] = gtf['attribute'].str.extract(r'exon_number "(\d+)"')

    print(' -- KTC_plot_transcript_splicing : Reading rMATS...')
    rmats = pd.read_csv(rmats_file, sep='\t')
    rmats = rmats[rmats['Splicing Event'] == splicing_event_type]
    event = rmats[rmats.iloc[:, 0] == event_id]
    if event.empty:
        print(f' -- KTC_plot_transcript_splicing : Event ID {event_id} not found')
        return
    event = event.iloc[0]

    gene_id = event['geneSymbol']
    strand = event['strand']
    chrom = event['chr'].replace('chr', '')

    skipped_start = int(event['exonStart_0base']) + 1
    skipped_end = int(event['exonEnd'])

    print(f' -- Plotting transcripts for gene {gene_id} with skipped region {skipped_start}-{skipped_end}')

    # Filter exons for this gene
    gene_exons = gtf[
        (gtf['gene_name'] == gene_id) &
        (gtf['seqname'].str.replace('chr', '') == chrom) &
        (gtf['strand'] == strand) &
        (gtf['feature'] == 'exon')
    ].copy()

    gene_exons.sort_values(['transcript_id', 'start'], inplace=True)
    transcripts = gene_exons.groupby('transcript_id')

    # Identify canonical transcript
    canonical_tx = KTC_get_canonical_transcript(gene_id)
    if canonical_tx:
        print(f" -- Canonical transcript: {canonical_tx}")
    else:
        print(" -- Canonical transcript: ❌ Not found")

    fig, ax = plt.subplots(figsize=(12, max(4, len(transcripts) * 0.8)), dpi=200)

    def causes_frameshift(exons, skipped, tx_id):
        full_length = (exons['end'] - exons['start'] + 1).sum()
        skipped_exon_length = skipped[1] - skipped[0] + 1
        skipped_length = full_length - skipped_exon_length
    
        print()
        print(f" {tx_id}")
        print(f" -- Length of full transcript: {full_length} | Divisible by 3: {full_length % 3 == 0}")
        print(f" -- Skipped exon length: {skipped_exon_length} | Divisible by 3: {skipped_exon_length % 3 == 0}")
        print(f" -- Length without skipped exon: {skipped_length} | Divisible by 3: {skipped_length % 3 == 0}")
    
        return (skipped_length % 3) != 0

    for idx, (tx_id, exons) in enumerate(transcripts):
        y = -idx
        exons = exons.sort_values('start')
        prev_end = None
        is_canonical = (tx_id == canonical_tx)

        frame_shift = causes_frameshift(exons, (skipped_start, skipped_end), tx_id)

        # Background highlight if canonical
        if is_canonical:
            ax.add_patch(plt.Rectangle(
                (gene_exons['start'].min() - 800, y - 0.4),
                gene_exons['end'].max() - gene_exons['start'].min() + 1000,
                0.8, color='yellow', alpha=0.2
            ))

        for _, exon in exons.iterrows():
            start, end = exon['start'], exon['end']
            color = 'black'
            if start == skipped_start and end == skipped_end:
                color = 'red' if frame_shift else 'green'

            ax.add_patch(plt.Rectangle((start, y - 0.2), end - start + 1, 0.4, color=color))

            if prev_end is not None:
                ax.plot([prev_end, start], [y, y], color='black', linewidth=1.0)

            prev_end = end

        label_style = {'fontsize': 9, 'va': 'center', 'ha': 'right'}
        # if is_canonical:
        #     label_style['fontweight'] = 'bold'
            # label_style['color'] = 'darkblue'

        ax.text(gene_exons['start'].min() - 300, y, tx_id, **label_style)

    ax.set_ylim(y - 1, 1)
    ax.set_xlim(gene_exons['start'].min() - 500, gene_exons['end'].max() + 500)
    ax.axis('off')

    plt.title(f'Gene: {gene_id} | Skipped exon: {skipped_start}-{skipped_end}')
    plt.tight_layout()
    plt.show()

#%% =============================================================================
# 
# =============================================================================
def KTC_make_list(not_list):
    # Split by lines, strip whitespace, and remove empty entries
    lines = [line.strip() for line in not_list.strip().splitlines() if line.strip()]
    # Return as list of strings
    return lines


#%% =============================================================================
# 
# =============================================================================

# =========================
# KTC_ncbi_gene_scraper
# =========================
from typing import List, Union, Optional
import time
import pandas as pd

try:
    from Bio import Entrez
except Exception as e:
    raise ImportError("Biopython is required. Install with: pip install biopython") from e


def KTC_ncbi_gene_scraper(
    genes: Union[str, List[str]],
    email: str = "kasperthorhauge.christensen@ugent.be",
    organism: str = "Homo sapiens",
    print_summaries: bool = True,
    sleep_seconds: float = 0.4,
) -> pd.DataFrame:
    """
    Query NCBI Gene (Entrez) for summaries and metadata.

    Parameters
    ----------
    genes : str | list[str]
        Either a multiline string with one gene per line, or a list of gene symbols.
    email : str
        Email required by NCBI usage policy.
    organism : str
        Organism filter for the query (default Homo sapiens).
    print_summaries : bool
        If True, prints 'SYMBOL\tSummary' lines while fetching.
    sleep_seconds : float
        Pause between requests to respect NCBI rate limits.

    Returns
    -------
    pd.DataFrame
        Columns: ['query', 'gene_id', 'symbol', 'full_name', 'description',
                  'summary', 'chromosome', 'map_location', 'other_aliases', 'status', 'organism']
    """

    # Configure Entrez
    Entrez.email = email
    Entrez.tool = "KTC_functions_ncbi_gene_scraper"

    def _normalize_gene_list(x: Union[str, List[str]]) -> List[str]:
        if isinstance(x, str):
            raw = [g.strip() for g in x.splitlines()]
        else:
            raw = [str(g).strip() for g in x]
        # keep order, drop empties, drop 'nan' literals, uppercase, deduplicate preserving order
        seen = set()
        clean = []
        for g in raw:
            if not g or g.lower() == "nan":
                continue
            g_up = g.upper()
            if g_up not in seen:
                seen.add(g_up)
                clean.append(g_up)
        return clean

    def _safe_get(d: dict, key: str, default: Optional[str] = "") -> str:
        try:
            return d.get(key, default) or default
        except Exception:
            return default

    gene_list = _normalize_gene_list(genes)
    records = []

    for symbol in gene_list:
        # Search for the gene ID
        term = f"{symbol}[Gene] AND {organism}[orgn]"
        try:
            with Entrez.esearch(db="gene", term=term, retmax=1) as handle:
                search = Entrez.read(handle)
        except Exception as e:
            # Record a failed lookup line and continue
            records.append({
                "query": symbol,
                "gene_id": "",
                "symbol": "",
                "full_name": "",
                "description": "",
                "summary": f"ERROR during esearch: {e}",
                "chromosome": "",
                "map_location": "",
                "other_aliases": "",
                "status": "",
                "organism": organism,
            })
            time.sleep(sleep_seconds)
            continue

        idlist = search.get("IdList", [])
        if not idlist:
            msg = f"No gene found for {symbol} in {organism}"
            if print_summaries:
                print(f"{symbol}\t{msg}")
                print()
            records.append({
                "query": symbol,
                "gene_id": "",
                "symbol": "",
                "full_name": "",
                "description": "",
                "summary": msg,
                "chromosome": "",
                "map_location": "",
                "other_aliases": "",
                "status": "",
                "organism": organism,
            })
            time.sleep(sleep_seconds)
            continue

        gene_id = idlist[0]

        # Fetch summary for that ID
        try:
            with Entrez.esummary(db="gene", id=gene_id) as handle:
                esum = Entrez.read(handle)
            doc = esum["DocumentSummarySet"]["DocumentSummary"][0]
        except Exception as e:
            records.append({
                "query": symbol,
                "gene_id": gene_id,
                "symbol": "",
                "full_name": "",
                "description": "",
                "summary": f"ERROR during esummary: {e}",
                "chromosome": "",
                "map_location": "",
                "other_aliases": "",
                "status": "",
                "organism": organism,
            })
            time.sleep(sleep_seconds)
            continue

        # Extract fields safely
        official_symbol = _safe_get(doc, "NomenclatureSymbol") or _safe_get(doc, "Name")
        full_name       = _safe_get(doc, "NomenclatureName")
        description     = _safe_get(doc, "Description")
        summary_text    = _safe_get(doc, "Summary")
        # Trim any trailing bracketed references for your preferred printout
        printable_sum   = summary_text.split("[")[0].strip() if summary_text else ""

        if print_summaries:
            print(f"{symbol}\t{printable_sum if printable_sum else 'No summary available'}")
            print()

        records.append({
            "query": symbol,
            "gene_id": gene_id,
            "symbol": official_symbol,
            "full_name": full_name,
            "description": description,
            "summary": summary_text,
            "chromosome": _safe_get(doc, "Chromosome"),
            "map_location": _safe_get(doc, "MapLocation"),
            "other_aliases": _safe_get(doc, "OtherAliases"),
            "status": _safe_get(doc, "Status"),
            "organism": organism,
        })

        time.sleep(sleep_seconds)

    df = pd.DataFrame.from_records(records)
    # If a query matched, prefer the NCBI official symbol; otherwise keep the query
    df["symbol"] = df["symbol"].mask(df["symbol"].eq(""), df["query"])
    return df

#%% =============================================================================
# Extracting gene names from a URL
# =============================================================================

def KTC_GenesFromURL(url, write_file=False, out_path=None):
    """
    Fetch Harmonizome gene set from URL and extract gene symbols.
    """

    response = requests.get(url)
    response.raise_for_status()  # fail loudly if something is wrong

    data = response.json()

    genes = [
        entry["gene"]["symbol"]
        for entry in data.get("associations", [])
        if "gene" in entry and "symbol" in entry["gene"]
    ]

    # remove duplicates while preserving order
    genes = list(dict.fromkeys(genes))

    if write_file and out_path is not None:
        with open(out_path, "w") as f:
            for g in genes:
                f.write(g + "\n")

    return genes

# print(KTC_GenesFromURL('https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/NAMPT/Pathway+Commons+Protein-Protein+Interactions'))


#%% =============================================================================
# GSEA
# =============================================================================
from gseapy.plot import gseaplot


def KTC_GSEA(
    deseq_path,
    genes,                  # your list of HGNC symbols
    set_name="My_Gene_Set",
    rank_by="stat",         # recommended: "stat"
    permutations=2000,
    min_size=5,
    max_size=500,
    out_png="gsea_plot.png",
    seed=0,
    ds = None
):
    """
    Runs preranked GSEA on a single custom gene list and produces
    a classic enrichment plot.
    """

    # -----------------------------
    # 1. Load DESeq2 results
    # -----------------------------
    df = pd.read_csv(deseq_path, sep=None, engine="python")

    required_cols = {"hgnc_symbol", rank_by}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"DESeq file must contain: {required_cols}")

    df["hgnc_symbol"] = df["hgnc_symbol"].astype(str).str.strip()

    # Remove bad rows
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["hgnc_symbol", rank_by])

    # -----------------------------
    # 2. Build ranked list
    # -----------------------------
    rnk = df[["hgnc_symbol", rank_by]].copy()

    # If duplicated gene symbols exist, keep strongest absolute score
    rnk["abs_score"] = rnk[rank_by].abs()
    rnk = rnk.sort_values("abs_score", ascending=False)
    rnk = rnk.drop_duplicates("hgnc_symbol")

    rnk = rnk.sort_values(rank_by, ascending=False)
    rnk = rnk[["hgnc_symbol", rank_by]]
    rnk.columns = ["gene", "score"]

    # -----------------------------
    # 3. Run preranked GSEA
    # -----------------------------
    gene_sets = {set_name: list(set(genes))}

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutations,
        seed=seed,
        outdir=None,  # don't write files
        verbose=False
    )

    # -----------------------------
    # 4. Extract results
    # -----------------------------
    res = pre.results[set_name]

    # -----------------------------
    # 5. Plot classic GSEA figure
    # -----------------------------
    plot_obj = gseaplot(
        rank_metric=rnk["score"].values,
        term=set_name,
        hits=res["hits"],
        RES=res["RES"],
        nes=res["nes"],
        pval=res["pval"],
        fdr=res["fdr"],
        figsize=(6,5)
    )

    # Handle different return types across gseapy versions
    if isinstance(plot_obj, tuple):
        fig = plot_obj[0]
    elif isinstance(plot_obj, list):
        fig = plot_obj[0].figure
    else:
        fig = plot_obj.figure

    fig.tight_layout()
    # fig.savefig(out_png, dpi=300)
    if ds != None:
        plt.title(f'{ds}')
    plt.show()

    return pre, res

