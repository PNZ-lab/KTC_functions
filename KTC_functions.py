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
    'm6a readers'      : ['HNRNPC', 'YTHDF1', 'YTHDF2', 'YTHDF3', 'YTHDC1', 'YTHDC2', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'HNRNPA2B1'],
    'splicing_factors' : ['SFRS7', 'CUGBP1', 'DAZAP1', 'CUGBP2', 'FMR1', 'A2BP1', 'RBFOX2', 'HNRNPA0', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPC', 'HNRNPC', 'HNRNPD', 'HNRNPD', 'HNRPDL', 'PCBP1', 'PCBP2', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPH3', 'PTBP1', 'HNRNPK', 'HNRNPK', 'HNRNPL', 'HNRPLL', 'HNRNPM', 'FUS', 'HNRNPU', 'TRA2A', 'TRA2B', 'ELAVL2', 'ELAVL4', 'ELAVL1', 'KHSRP', 'MBNL1', 'NOVA1', 'NOVA2', 'PTBP2', 'SFPQ', 'RBM25', 'RBM4', 'KHDRBS1', 'SF3B1', 'SFRS2', 'SF1', 'SFRS1', 'KHDRBS2', 'KHDRBS3', 'SFRS3', 'SFRS9', 'SFRS13A', 'SFRS5', 'SFRS11', 'SFRS6', 'SFRS4', 'TARDBP', 'TIA1', 'TIAL1', 'YBX1', 'ZRANB2', 'ELAVL3', 'RBM5', 'SYNCRIP', 'HNRNPA3', 'QKI', 'RBMX', 'SRRM1', 'ESRP1', 'ESRP2'], # From SpliceAid
    'EpiFactors'       : ['A1CF', 'ACINU', 'ACTB', 'ACTL6A', 'ACTL6B', 'ACTR3B', 'ACTR5', 'ACTR6', 'ACTR8', 'ADNP', 'AEBP2', 'AICDA', 'AIRE', 'ALKBH1', 'ALKBH1', 'ALKBH4', 'ALKBH5', 'ANKRD32', 'ANP32A', 'ANP32B', 'ANP32E', 'APBB1', 'APEX1', 'APOBEC1', 'APOBEC2', 'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H', 'ARID1A', 'ARID1B', 'ARID2', 'ARID4A', 'ARID4B', 'ARNTL', 'ARRB1', 'ASF1A', 'ASF1B', 'ASH1L', 'ASH2L', 'ASXL1', 'ASXL2', 'ASXL3', 'ATAD2', 'ATAD2B', 'ATF2', 'ATF7IP', 'ATM', 'ATN1', 'ATR', 'ATRX', 'ATXN7', 'ATXN7L3', 'AURKA', 'AURKB', 'AURKC', 'BABAM1', 'BAHD1', 'BANP', 'BAP1', 'BARD1', 'BAZ1A', 'BAZ1B', 'BAZ2A', 'BAZ2B', 'BCOR', 'BCORL1', 'BMI1', 'BPTF', 'BRCA1', 'BRCA2', 'BRCC3', 'BRD1', 'BRD2', 'BRD3', 'BRD4', 'BRD7', 'BRD8', 'BRD9', 'BRDT', 'BRE', 'BRMS1', 'BRMS1L', 'BRPF1', 'BRPF3', 'BRWD1', 'BRWD3', 'BUB1', 'C11orf30', 'C14orf169', 'C17orf49', 'CARM1', 'CBLL1', 'CBX1', 'CBX2', 'CBX3', 'CBX4', 'CBX5', 'CBX6', 'CBX7', 'CBX8', 'CCDC101', 'CDC6', 'CDC73', 'CDK1', 'CDK17', 'CDK2', 'CDK3', 'CDK5', 'CDK7', 'CDK9', 'CDY1', 'CDY1B', 'CDY2A', 'CDY2B', 'CDYL', 'CDYL2', 'CECR2', 'CELF1', 'CELF2', 'CELF3', 'CELF4', 'CELF5', 'CELF6', 'CENPC', 'CHAF1A', 'CHAF1B', 'CHD1', 'CHD1L', 'CHD2', 'CHD3', 'CHD4', 'CHD5', 'CHD6', 'CHD7', 'CHD8', 'CHD9', 'CHEK1', 'CHRAC1', 'CHTOP', 'CHUK', 'CIR1', 'CIT', 'CLNS1A', 'CLOCK', 'CRB2', 'CREBBP', 'CSNK2A1', 'CSRP2BP', 'CTBP1', 'CTBP2', 'CTCF', 'CTCFL', 'CTR9', 'CUL1', 'CUL2', 'CUL3', 'CUL4A', 'CUL4B', 'CUL5', 'CXXC1', 'DAPK3', 'DAXX', 'DDB1', 'DDB2', 'DDX17', 'DDX21', 'DDX5', 'DDX50', 'DEK', 'DHX9', 'DMAP1', 'DNAJC1', 'DNAJC2', 'DND1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DNMT3L', 'DNTTIP2', 'DOT1L', 'DPF1', 'DPF2', 'DPF3', 'DPPA3', 'DPY30', 'DR1', 'DTX3L', 'DZIP3', 'E2F6', 'EED', 'EEF1AKMT3', 'EEF1AKMT4', 'EEF1AKNMT', 'EHMT1', 'EHMT2', 'EID1', 'EID2', 'EID2B', 'EIF4A3', 'ELP2', 'ELP3', 'ELP4', 'ELP5', 'ELP6', 'ENY2', 'EP300', 'EP400', 'EPC1', 'EPC2', 'ERBB4', 'ERCC6', 'EXOSC1', 'EXOSC2', 'EXOSC3', 'EXOSC4', 'EXOSC5', 'EXOSC6', 'EXOSC7', 'EXOSC8', 'EXOSC9', 'EYA1', 'EYA2', 'EYA3', 'EYA4', 'EZH1', 'EZH2', 'FAM175A', 'FAM175B', 'FBL', 'FBRS', 'FBRSL1', 'FOXA1', 'FOXO1', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXP4', 'FTO', 'GADD45A', 'GADD45B', 'GADD45G', 'GATAD1', 'GATAD2A', 'GATAD2B', 'GFI1', 'GFI1B', 'GLYR1', 'GSE1', 'GSG2', 'GTF2I', 'GTF3C4', 'HAT1', 'HCFC1', 'HCFC2', 'HDAC1', 'HDAC10', 'HDAC11', 'HDAC2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6', 'HDAC7', 'HDAC8', 'HDAC9', 'HDGF', 'HDGFL2', 'HELLS', 'HIF1AN', 'HINFP', 'HIRA', 'HIRIP3', 'HJURP', 'HLCS', 'HLTF', 'HMG20A', 'HMG20B', 'HMGB1', 'HMGN1', 'HMGN2', 'HMGN3', 'HMGN4', 'HMGN5', 'HNRNPU', 'HNRPL', 'HNRPM', 'HP1BP3', 'HR', 'HSFX3', 'HSPA1A', 'HSPA1A', 'HSPA1B', 'HSPA1B', 'HUWE1', 'IKBKAP', 'IKZF1', 'IKZF3', 'ING1', 'ING2', 'ING3', 'ING4', 'ING5', 'INO80', 'INO80B', 'INO80C', 'INO80D', 'INO80E', 'JADE1', 'JADE2', 'JADE3', 'JAK2', 'JARID2', 'JDP2', 'JMJD1C', 'JMJD6', 'KANSL1', 'KANSL2', 'KANSL3', 'KAT2A', 'KAT2B', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'KAT8', 'KDM1A', 'KDM1B', 'KDM2A', 'KDM2B', 'KDM3A', 'KDM3B', 'KDM4A', 'KDM4B', 'KDM4C', 'KDM4D', 'KDM4E', 'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'KDM6A', 'KDM6B', 'KDM7A', 'KDM8', 'KEAP1', 'KHDRBS1', 'KLF18', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'KMT2E', 'L3MBTL1', 'L3MBTL2', 'L3MBTL3', 'L3MBTL4', 'LAS1L', 'LBR', 'LEO1', 'LRWD1', 'MAGOH', 'MAP3K7', 'MAPKAPK3', 'MASTL', 'MAX', 'MAZ', 'MBD1', 'MBD2', 'MBD3', 'MBD4', 'MBD5', 'MBD6', 'MBIP', 'MBNL1', 'MBNL3', 'MBTD1', 'MCRS1', 'MDC1', 'MEAF6', 'MECP2', 'MEN1', 'METTL11B', 'METTL14', 'METTL16', 'METTL21A', 'METTL3', 'METTL4', 'MGA', 'MGEA5', 'MINA', 'MIS18A', 'MIS18BP1', 'MLLT1', 'MLLT10', 'MLLT6', 'MORF4L1', 'MORF4L2', 'MOV10', 'MPHOSPH8', 'MPND', 'MRGBP', 'MSH6', 'MSL1', 'MSL2', 'MSL3', 'MST1', 'MTA1', 'MTA2', 'MTA3', 'MTF2', 'MUM1', 'MYBBP1A', 'MYO1C', 'MYSM1', 'NAA60', 'NAP1L1', 'NAP1L2', 'NAP1L4', 'NASP', 'NAT10', 'NAT10', 'NBN', 'NCL', 'NCOA1', 'NCOA2', 'NCOA3', 'NCOA6', 'NCOR1', 'NCOR2', 'NEK6', 'NEK9', 'NFRKB', 'NFYB', 'NFYC', 'NIPBL', 'NOC2L', 'NPAS2', 'NPM1', 'NPM2', 'NSD1', 'NSL1', 'NSRP1', 'NSUN2', 'NSUN6', 'NTMT1', 'NUP98', 'OGT', 'OIP5', 'PADI1', 'PADI2', 'PADI3', 'PADI4', 'PAF1', 'PAGR1', 'PAK2', 'PARG', 'PARP1', 'PARP2', 'PARP3', 'PAXIP1', 'PBK', 'PBRM1', 'PCGF1', 'PCGF2', 'PCGF3', 'PCGF5', 'PCGF6', 'PCNA', 'PDP1', 'PELP1', 'PHC1', 'PHC2', 'PHC3', 'PHF1', 'PHF10', 'PHF12', 'PHF13', 'PHF14', 'PHF19', 'PHF2', 'PHF20', 'PHF20L1', 'PHF21A', 'PHF8', 'PHIP', 'PIWIL4', 'PKM', 'PKN1', 'POGZ', 'POLE3', 'PPARGC1A', 'PPM1G', 'PPP2CA', 'PPP4C', 'PPP4R2', 'PQBP1', 'PRDM1', 'PRDM11', 'PRDM12', 'PRDM13', 'PRDM14', 'PRDM16', 'PRDM2', 'PRDM4', 'PRDM5', 'PRDM6', 'PRDM7', 'PRDM8', 'PRDM9', 'PRKAA1', 'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKAG3', 'PRKCA', 'PRKCB', 'PRKCD', 'PRKDC', 'PRMT1', 'PRMT2', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9', 'PRPF31', 'PRR14', 'PSIP1', 'PTBP1', 'PTBP1', 'PUF60', 'RAD51', 'RAD54B', 'RAD54L', 'RAD54L2', 'RAG1', 'RAG2', 'RAI1', 'RARA', 'RB1', 'RBBP4', 'RBBP5', 'RBBP7', 'RBFOX1', 'RBM11', 'RBM15', 'RBM15B', 'RBM17', 'RBM24', 'RBM25', 'RBM4', 'RBM5', 'RBM7', 'RBM8A', 'RBMY1A1', 'RBX1', 'RCC1', 'RCOR1', 'RCOR3', 'REST', 'RFOX1', 'RING1', 'RLIM', 'RMI1', 'RNF168', 'RNF2', 'RNF20', 'RNF40', 'RNF8', 'RNPS1', 'RPS6KA3', 'RPS6KA4', 'RPS6KA5', 'RPUSD3', 'RRP8', 'RSF1', 'RSRC1', 'RUVBL1', 'RUVBL2', 'RYBP', 'SAFB', 'SAP130', 'SAP18', 'SAP25', 'SAP30', 'SAP30L', 'SATB1', 'SATB2', 'SCMH1', 'SCML2', 'SCML4', 'SENP1', 'SENP3', 'SET', 'SETD1A', 'SETD1B', 'SETD2', 'SETD3', 'SETD5', 'SETD6', 'SETD7', 'SETD8', 'SETDB1', 'SETDB2', 'SETMAR', 'SF3B1', 'SF3B3', 'SFMBT1', 'SFMBT2', 'SFPQ', 'SFSWAP', 'SHPRH', 'SIN3A', 'SIN3B', 'SIRT1', 'SIRT2', 'SIRT6', 'SIRT7', 'SKP1', 'SLU7', 'SMARCA1', 'SMARCA2', 'SMARCA4', 'SMARCA5', 'SMARCAD1', 'SMARCAL1', 'SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2', 'SMARCD3', 'SMARCE1', 'SMEK1', 'SMEK2', 'SMYD1', 'SMYD2', 'SMYD3', 'SMYD4', 'SNAI2', 'SP1', 'SP100', 'SP140', 'SPEN', 'SPOP', 'SRCAP', 'SRRM4', 'SRSF1', 'SRSF10', 'SRSF12', 'SRSF3', 'SRSF6', 'SS18L1', 'SS18L2', 'SSRP1', 'STK4', 'SUDS3', 'SUPT16H', 'SUPT3H', 'SUPT6H', 'SUPT7L', 'SUV39H1', 'SUV39H2', 'SUV420H1', 'SUV420H2', 'SUZ12', 'SYNCRIP', 'TADA1', 'TADA2A', 'TADA2B', 'TADA3', 'TAF1', 'TAF10', 'TAF12', 'TAF1L', 'TAF2', 'TAF3', 'TAF4', 'TAF5', 'TAF5L', 'TAF6', 'TAF6L', 'TAF7', 'TAF8', 'TAF9', 'TAF9B', 'TBL1XR1', 'TDG', 'TDRD3', 'TDRD7', 'TDRKH', 'TET1', 'TET2', 'TET3', 'TEX10', 'TFDP1', 'TFPT', 'THRAP3', 'TLE1', 'TLE2', 'TLE4', 'TLK1', 'TLK2', 'TNP1', 'TNP2', 'TONSL', 'TOP2A', 'TOP2B', 'TP53', 'TP53BP1', 'TRA2B', 'TRIM16', 'TRIM24', 'TRIM27', 'TRIM28', 'TRIM33', 'TRRAP', 'TRUB2', 'TSSK6', 'TTK', 'TYW5', 'U2AF2', 'UBE2A', 'UBE2B', 'UBE2D1', 'UBE2D3', 'UBE2E1', 'UBE2H', 'UBE2N', 'UBE2T', 'UBN1', 'UBR2', 'UBR5', 'UBR7', 'UCHL5', 'UHRF1', 'UHRF2', 'UIMC1', 'USP11', 'USP12', 'USP15', 'USP16', 'USP17L2', 'USP21', 'USP22', 'USP3', 'USP36', 'USP44', 'USP46', 'USP49', 'USP7', 'UTY', 'VDR', 'VIRMA', 'VPS72', 'VRK1', 'WAC', 'WDR5', 'WDR77', 'WDR82', 'WHSC1', 'WHSC1L1', 'WSB2', 'WTAP', 'YAF2', 'YEATS2', 'YEATS4', 'YTHDC1', 'YWHAB', 'YWHAE', 'YWHAZ', 'YY1', 'ZBTB16', 'ZBTB33', 'ZBTB7A', 'ZBTB7C', 'ZC3H13', 'ZCWPW1', 'ZFP57', 'ZGPAT', 'ZHX1', 'ZMYM2', 'ZMYM3', 'ZMYND11', 'ZMYND8', 'ZNF217', 'ZNF516', 'ZNF532', 'ZNF541', 'ZNF592', 'ZNF687', 'ZNF711', 'ZNHIT1', 'ZRANB3', 'ZZZ3'],
    'm6a_writers'      : ['METTL3', 'METTL14', 'METTL16', 'KIAA1429','RBM15', 'WTAP'],
    'm6a_erasers'      : ['FTO', 'ALKBH5'],
    'm6a_readers'      : ['YTHS', 'EIF3', 'HNRNPC', 'HNRNPA2B1', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3'],
    'm6a_re_wr_er'     : ['METTL3', 'METTL14',  'KIAA1429','RBM15', 'WTAP', 'FTO', 'ALKBH5', 'YTHS', 'EIF3', 'HNRNPC', 'HNRNPA2B1', 'YTHDF1', 'YTHDF2', 'YTHDC1', 'YTHDC2', 'TYSND1', 'SND1', 'PRRC2A', 'LRPPRC', 'FMR1','FMR1NB', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3'],
    'PRC2'             : ['EZH1', 'EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7', 'JARID2', 'PCGF1', 'PCGF2', 'RING1', 'BMI1'],
    'Freya'            : ['ALKBH5', 'CYP51A1', 'DHCR7', 'DHCR24', 'EBP', 'FDFT1', 'FDPS', 'FTO', 'GGPS1', 'HNRNPC', 'HSD17B7', 'IDI1', 'IGF2BP2', 'LDLR', 'LSS', 'METTL3', 'METTL14', 'MSMO1', 'MVD', 'MVK', 'NSDHL', 'PMVK', 'SC5D', 'SQLE', 'YTHDF1', 'YTHDF2'],
    'Laura'            : ['NAMPT', 'NAPRT', 'IDO', 'DHFR', 'NMNAT1', 'NMNAT2', 'NMNAT3', 'QPRT', 'MAT2A', 'MTAP', 'WTAP', 'E2F1', 'NRK', 'TDO', 'NADSYN'],
    'Kevin'            : ['PRPF8', 'SRRM1', 'SRRM2', 'ACIN1', 'RNPS1', 'CLK1', 'CLK2', 'CLK3', 'CLK4'],
    'CM-specific'      : ['Top2a', 'Mki67', 'Prc1', 'Mis18bp1', 'Kif23', 'Kif15', 'Neil3', 'Cenpe', 'Knl1', 'Iqgap3', 'Aspm', 'Racgap1', 'Foxm1', 'Lockd', 'Esco2', 'Cdca8', 'Anln', 'Cdca2', 'Diaph3', 'Ckap2l', 'Birc5', 'Kif11', 'Sgo1', 'Cdca3', 'Knstrn', 'Uhrf1', 'Ttk', 'Cdk1', 'Ankle1', 'Ncapg', 'Spag5', 'Tpx2', 'E2f8', 'Kif4', 'Kif20b', 'Shcbp1', 'Ncapg2', 'Nusap1', 'E2f7', 'Cit', 'Rad51ap1', 'Cenpf', 'Trim59', 'Ncapd2', 'Stmn1', 'Gm42047', 'Smc2', 'Incenp', 'Lmnb1', 'Atad2', 'Plk4', 'Cep128', 'Arl4c', 'Hmgb2', 'Sgo2a', 'Ezh2', 'Lair1', 'Cep192', 'Smc4', 'Hirip3', 'Lcp1', 'Ncapd3', 'Topbp1', 'Laptm5', 'Cdkn2c', 'Runx1', 'Nsd2', 'Hjurp', 'Fam111a', 'H2afz', 'Arhgap30', 'C1qa', 'Ptprc', 'Cenpa', 'Maf', 'F630028O10Rik', 'C1qc', 'Tmpo', 'Lyz2', 'Mrc1', 'Slbp', 'Ctsc', 'Lbr', 'Gatm', 'F13a1', 'Arhgap45', 'C1qb', 'Ccdc82', 'Nav2', 'Tm6sf1', 'Smchd1', 'Nucks1', 'G2e3', 'Apobec3', 'Tubb5', 'Myo5a', 'Dab2', 'Smc1a', 'Git2', 'Arl6ip1', 'Mhrt', 'Gm31251', 'Cpeb3', 'Vegfa', 'Atcayos', 'Prune2', 'Sorbs1', 'Sorbs2', 'Pde4d', 'Rnf207', 'Lmo7', 'Coro6', 'Ryr2', 'Ppip5k2', 'Ivns1abp', 'Myh7b', 'Rbm24', 'Smtn', 'Trim63', 'Dmd', 'D830005E20Rik', 'Rbm20', 'D830024N08Rik', 'Nav2', 'Ppargc1a', 'Ccdc141', 'Pde4dip', 'Dmpk', 'Tnnt2', 'Pde7a', 'Cacna1c', 'Slc4a3', 'Clasp1', 'Gja1', 'Mybpc3', 'Nexn', 'Pcdh7', 'Gm36827', 'Lgals4', 'Obscn', 'Cacnb2', 'Camk2d', 'Palld', 'Neat1', 'Ggnbp1', 'Ttn', 'Tacc2', 'Mapt', 'Pfkfb2', 'Nnt', 'Slc27a1', 'Trim7', 'Ank3', 'Carns1', 'Clip1', 'Ralgapa2', 'Vldlr', 'Ppargc1b', 'Alpk3', '5430431A17Rik', 'Ank2', 'Acacb', 'Myom2', 'Art1', 'Lrrc2', 'Cdh2', 'Enah', 'Dtna', 'Pkp2', 'Lrrfip2', 'Myzap', 'Dot1l', 'Kidins220', 'Tbc1d4', 'Magi2', 'Mlxipl', 'Gpcpd1', 'Hk2', 'Corin', 'Ctnna3', 'Kbtbd12', 'Asph', 'Mlip', 'Retreg1', 'Fblim1', 'Pacsin3', 'Mical3', 'Speg', 'Ldb3', 'Tango2', 'Nebl', 'Agl', 'Grip2', 'Lmod2', 'Ppip5k1', 'Phkg1', 'Cenpa', 'Fhl2', 'Car14', 'Lrtm1', 'D830005E20Rik', 'Rnf207', 'Mhrt', 'Gm31251', 'Atcayos', 'Rbm24', 'Cpeb3', 'Cacnb2', 'Pde7a', 'Myh7b', 'D830024N08Rik', 'Cacna1c', 'Trim63', 'Ppargc1a', 'Mlip', 'Rbm20', 'Pfkfb2', 'Ctnna3', 'Ppip5k2', 'Pkp2', 'Lgals4', 'Vegfa', 'Trim55', 'Ccdc141', 'Mlxipl', 'Coro6', 'Nav2', 'Pde4d', 'Ppp1r3a', 'Sgcd', 'Gm47101', 'Mfn1', 'Nabp1', 'Ank3', 'Dtna', 'Alpk3', 'Ggnbp1', 'Nnt', 'Gm36827', 'Sorbs2', 'Slc8a1', 'Sox6', 'Lmo7', 'Palld', 'Plin4', 'Nexn', 'Lmod2', 'Dmd', 'Gja1', 'Corin', 'Akap6', 'Smtn', 'Slc4a3', 'Art1', 'Carns1', 'Sorbs1', 'Grb14', 'Speg', 'Trim7', 'Prune2', 'Kcnip2', 'Cdh2', 'Ryr2', 'Gm47547', 'Phkg1', '3222401L13Rik', 'Enah', 'Vldlr', 'Lrrc10', 'C130080G10Rik', 'Nebl', 'Pcdh7', 'Agl', 'Dmpk', 'Myzap', 'Clip1', 'Tnnt2', 'Pde4dip', 'Rrad', 'Lrrc2', 'Ppip5k1', 'Neat1', 'Pcgf5', 'Nceh1', 'Obscn', 'Ndufaf4', 'Svil', 'Ank2', 'Tacc2', 'Clasp1', 'Ivns1abp', 'Ttn', 'Asph', '2010111I01Rik', 'Malat1', 'Camk2d', 'Mybpc3', 'Tnnc1', 'Fhl2', 'Ptgds', 'Mb', 'Myl2', 'Actc1', 'Fabp3', 'Cox6a2', 'Slc25a4', 'Myl3', 'Ckm', 'Tpm1', 'Myh6', 'Tnni3', 'Cox7a1', 'Hspb7', 'Des', 'Atp5a1', 'Atp5b', 'Atp5g3', 'Ckmt2', 'Oxct1', 'Cox6c', 'Cryab', 'Cox5a', 'Atp5g1', 'Atp5k', 'Pgam2', 'Cox4i1', 'Ndufa4', 'Chchd10', 'Tcap', 'Atp2a2', 'Pln', 'Cox7c', 'Ttn', 'Tuba4a', 'Cox8b', 'Uqcrq', 'Idh2', 'Ndufa5', 'Uqcr11', 'Myh7', 'Atp5h', 'Actn2', 'Atp5f1', 'Ldhb', 'S100a1', 'Csrp3', 'Aldoa', 'Atp5e', 'Atp5o.1', 'Uqcrfs1', 'Mdh2', 'Sod2', 'Ndufb9', 'Ankrd1', 'Slc25a3', 'Ndufc1', 'Pdha1', 'Mdh1', 'Ech1', 'Cox7b', 'Ndufa13', 'Atp5c1', 'Hrc', 'Atp5j2', 'Xirp2', 'Ndufs2', 'Uqcrh', 'Uqcrb', '2010107E04Rik', 'Atp5j', 'Ndufa2', 'Ndufs6', 'Uqcrc1', 'Acadl', 'Etfb', 'Uqcrc2', 'Usmg5', 'Ndufb8', 'Ndufb7', 'Ndufa1', 'Got1', 'Srl', 'Aco2', 'Dsp', 'Cox6b1', 'Cox5b', 'Cyc1', 'Uqcr10', 'Jph2', 'Acaa2', 'Acadvl', 'Ndufs7', 'Eno3', 'Ndufab1', 'Etfa', 'Trp53inp2', 'Cmya5', 'Zfp106', 'Acta1', 'Dynll2', 'Myl4', 'Myl7', 'Sln', 'Dkk3', 'Stard10', 'Mybphl', 'Sbk3', 'Bmp10', 'Myl1', 'Clu', 'Nppa', 'Myl9', 'Kcnj3', 'Pam', 'Tbx5', 'Gpx3', 'Nudt4', 'Smpx', 'Atp2a2', 'Cryab', 'Ankrd1', 'Chchd10', 'Tnni3', 'Kcnk3', 'Cox6a2', 'Actc1', 'Myh6', 'Cox7c', 'Tpm1', 'Ndufa1', 'Cox6c', 'Cox7a1', 'Atp5b', 'Ndufa4', 'Cox8b', 'Mb', 'Slc25a4', 'Uqcr11', 'Atp5g1', 'Mdh1', 'Uqcrh', 'Atp5a1', 'Csrp3', 'Cycs', 'Atp5e', 'Mtus2', 'Tmod1', 'S100a1', 'Atp5j', 'Cox5a', 'Des', 'Chchd2', 'Cox4i1', 'Atp5g3', 'Cox6b1', 'Mylk3', 'Tcap', 'Atp5k', 'Uqcrq', 'Atp5f1', 'Corin', 'Cox7b', 'Mdh2', 'Fndc5', 'Aldoa', 'Ehd4', 'Atp5d', 'Atp5h', 'Doc2g', 'Uqcrb', 'Uqcrfs1', 'Nppb', 'Hspb7', 'Uqcr10', 'Atp5l', 'Srl', 'Ndufa5', 'Myoz2', 'Angpt1', 'Ndufc1', 'Ndufa11', 'Ndufb9', 'Map1lc3a', 'Ndufa13', 'Gapdh', 'Cox5b', 'Ndufb4', '2010107E04Rik', 'Chrm2', 'Ndufb11', 'Vdac1', 'Atp5j2', 'Cox7a2', 'Atp5c1', 'Eif1', 'Ndufb10', 'Aes', 'Usmg5', 'Ctgf', 'Slc8a1'],
    'PRC2_consistent'  : ['IGF2BP2', 'PLEK', 'KIF21A', 'ID1', 'AFDN', 'SPATS2L', 'CTBP2', 'LMNA', 'RETN', 'BIN1', 'PTK2', 'CST7', 'AHNAK', 'ANXA5', 'CKAP4', 'GOLM1', 'MAT1A', 'SHTN1', 'KCTD12', 'PRKAR2B', 'AP3B2', 'ANXA1', 'HLA-B', 'SMARCA1', 'SPART', 'YPEL5', 'GLRX', 'ANO6', 'FNBP1', 'ACSF2', 'FLNB', 'ALOX5AP', 'SGSH', 'CD38', 'TMEM63A', 'LGALS9', 'ARHGAP25', 'GIMAP2', 'CD2', 'CAMK4', 'CD1A', 'LGALS3BP', 'CD28', 'GNA15', 'UBA7', 'CD1C', 'FYB1', 'TRAC', 'CD1E', 'MAGEA4'],
    'm6a_story'        : ['NFKB', 'IKZF3', 'IRF4', 'MTAP', 'MYC', 'FDFT1', 'BCL2', 'PTEN', 'ASB2', 'RARA', 'FTO', 'METTL3', 'ALKBH', 'YTHDF', 'HNRNP', 'HNRNPC', 'PTEN', 'SP1', 'CEBPA', 'DHCR7', 'HMGCS1', 'SOCS1', 'SOCS3', 'NOTCH1', 'MSI2', 'BRD4', ],
    'Cristina'         : ['KDM6B', 'SPI1', 'PU1', 'CD44', 'CEBPA', 'CEBPB', 'CEBPD', 'CEBPE', 'CEBPG', 'CEBPZ', 'LTF', 'TK1', 'DNMT1'],
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
        'JQ1', 'MXD1', 'MTAP', 'WTAP', 'CHEK1', 'RBM39', 'CCR9', 'CCR7', 'CD19', 'CDK6', 'BRD1', 'BRD2', 'BRD3', 'BRD4', 'APCDD1,' "FOXO1", 'PRPF39','MSI2', 'NFKIA','SERPINB1','CD69','CCL25', 'DCAF15',
        'TNFSF9', 'TNFRSF9', 'CD38', 'TRBC1', 'TRBC2', 'ZBTB7B', 'RUNX3','ZBP1', 'FBXW7', 'KMT2D', 'IGF2','TAL1','PTEN', 'IDH1','PIK3CA','PIK3R1','STAT5B','JAK1','JAK2','DNM2','IL7R','NGN3','CASC3', 'IDO',
        'MEF2C', 'TCF3','PBX1','HLF','ITGAM','AURKB','TLX3','STAG2','SOX11','XBP1','ID1','ID3','XBP1','CD34','LMO2', 'JAK3', 'RIGI','ADAR','ADARB1','ADARB2','IGHV','FOXM1','AURKA',
        'RELA','TERT','ARNT','HDAC1','HDAC2','HDAC6','HDAC10','NPM1','KDM1B'
        ],
    'Tibo' : ["PSMG1", "PSMG2", "PSMG3", "PSMG4", "POMP", "PSMB8", "PSMB10", "PSMB9", "HSPB1", "PSMA4", "PSMA6", "PSMA5", "PSMD4", "PSME1", "PSME2", "PSMD10", "PSME3", "PSMF1", "RAD23A", "PSMB5", "ADRM1", "PSMC3", "PSMB3", "PSMD8", "PSMB6", "PSMA7", "PSMD13", "PSMA3", "PSMB1", "PSMA2", "PSMB7", "PSMA1", "PSMB4", "PSMC1", "PSMC2", "PSMD3", "PSMD7", "PSMB2", "VCP", "TXNL1", "UBQLN1", "PSMD9", "ZFAND2A", "PSMD14", "PSMD1", "RAD23B", "PSMC5", "PSMC4", "PSMD2", "UBR1", "PSMC6", "PSMD6", "PSMD5", "UBE3A", "USP14", "PSMD11", "PSME4", "PSMD12", "PSMA8", "UCHL5", "PAAF1", "UBE3C", "ECPAS", "PRICKLE1"]
	}


# This function helps define gene sets to be highlighted.
# A set of genes is selected based on the string used as input.
# First, the string is tested in the custom set of gene sets in gene_sets dictionary above.
# If that fails, it will use the string to search for a public gene set on Msigdb
# If that fails, it will default to interpreting the input string as a set with a single gene name (the input string)

from gseapy import Msigdb
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
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                      names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                      dtype={0: str})
    print(' -- KTC_SpliceMap : Done')
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

