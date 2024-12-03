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
# This function helps define gene sets to be searched.
# A set of genes is selected based on the string used as input.
# First, the string is tested in the custom set of gene sets in gene_sets dictionary below.
# If that fails, it will use the string to search for a public gene set on Msigdb
# If that fails, it will default to interpreting the input string as a set with a single gene name (the input string)

gene_sets = {
    'None'             : [],
    'm6a readers'      : ['HNRNPC', 'YTHDF1', 'YTHDF2', 'YTHDF3', 'YTHDC1', 'YTHDC2', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'HNRNPA2B1'],
    'splicing factors' : ['SFRS7', 'CUGBP1', 'DAZAP1', 'CUGBP2', 'FMR1', 'A2BP1', 'RBFOX2', 'HNRNPA0', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPC', 'HNRNPC', 'HNRNPD', 'HNRNPD', 'HNRPDL', 'PCBP1', 'PCBP2', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPH3', 'PTBP1', 'HNRNPK', 'HNRNPK', 'HNRNPL', 'HNRPLL', 'HNRNPM', 'FUS', 'HNRNPU', 'TRA2A', 'TRA2B', 'ELAVL2', 'ELAVL4', 'ELAVL1', 'KHSRP', 'MBNL1', 'NOVA1', 'NOVA2', 'PTBP2', 'SFPQ', 'RBM25', 'RBM4', 'KHDRBS1', 'SF3B1', 'SFRS2', 'SF1', 'SFRS1', 'KHDRBS2', 'KHDRBS3', 'SFRS3', 'SFRS9', 'SFRS13A', 'SFRS5', 'SFRS11', 'SFRS6', 'SFRS4', 'TARDBP', 'TIA1', 'TIAL1', 'YBX1', 'ZRANB2', 'ELAVL3', 'RBM5', 'SYNCRIP', 'HNRNPA3', 'QKI', 'RBMX', 'SRRM1', 'ESRP1', 'ESRP2'], # From SpliceAid
    'EpiFactors'       : ['A1CF', 'ACINU', 'ACTB', 'ACTL6A', 'ACTL6B', 'ACTR3B', 'ACTR5', 'ACTR6', 'ACTR8', 'ADNP', 'AEBP2', 'AICDA', 'AIRE', 'ALKBH1', 'ALKBH1', 'ALKBH4', 'ALKBH5', 'ANKRD32', 'ANP32A', 'ANP32B', 'ANP32E', 'APBB1', 'APEX1', 'APOBEC1', 'APOBEC2', 'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G', 'APOBEC3H', 'ARID1A', 'ARID1B', 'ARID2', 'ARID4A', 'ARID4B', 'ARNTL', 'ARRB1', 'ASF1A', 'ASF1B', 'ASH1L', 'ASH2L', 'ASXL1', 'ASXL2', 'ASXL3', 'ATAD2', 'ATAD2B', 'ATF2', 'ATF7IP', 'ATM', 'ATN1', 'ATR', 'ATRX', 'ATXN7', 'ATXN7L3', 'AURKA', 'AURKB', 'AURKC', 'BABAM1', 'BAHD1', 'BANP', 'BAP1', 'BARD1', 'BAZ1A', 'BAZ1B', 'BAZ2A', 'BAZ2B', 'BCOR', 'BCORL1', 'BMI1', 'BPTF', 'BRCA1', 'BRCA2', 'BRCC3', 'BRD1', 'BRD2', 'BRD3', 'BRD4', 'BRD7', 'BRD8', 'BRD9', 'BRDT', 'BRE', 'BRMS1', 'BRMS1L', 'BRPF1', 'BRPF3', 'BRWD1', 'BRWD3', 'BUB1', 'C11orf30', 'C14orf169', 'C17orf49', 'CARM1', 'CBLL1', 'CBX1', 'CBX2', 'CBX3', 'CBX4', 'CBX5', 'CBX6', 'CBX7', 'CBX8', 'CCDC101', 'CDC6', 'CDC73', 'CDK1', 'CDK17', 'CDK2', 'CDK3', 'CDK5', 'CDK7', 'CDK9', 'CDY1', 'CDY1B', 'CDY2A', 'CDY2B', 'CDYL', 'CDYL2', 'CECR2', 'CELF1', 'CELF2', 'CELF3', 'CELF4', 'CELF5', 'CELF6', 'CENPC', 'CHAF1A', 'CHAF1B', 'CHD1', 'CHD1L', 'CHD2', 'CHD3', 'CHD4', 'CHD5', 'CHD6', 'CHD7', 'CHD8', 'CHD9', 'CHEK1', 'CHRAC1', 'CHTOP', 'CHUK', 'CIR1', 'CIT', 'CLNS1A', 'CLOCK', 'CRB2', 'CREBBP', 'CSNK2A1', 'CSRP2BP', 'CTBP1', 'CTBP2', 'CTCF', 'CTCFL', 'CTR9', 'CUL1', 'CUL2', 'CUL3', 'CUL4A', 'CUL4B', 'CUL5', 'CXXC1', 'DAPK3', 'DAXX', 'DDB1', 'DDB2', 'DDX17', 'DDX21', 'DDX5', 'DDX50', 'DEK', 'DHX9', 'DMAP1', 'DNAJC1', 'DNAJC2', 'DND1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DNMT3L', 'DNTTIP2', 'DOT1L', 'DPF1', 'DPF2', 'DPF3', 'DPPA3', 'DPY30', 'DR1', 'DTX3L', 'DZIP3', 'E2F6', 'EED', 'EEF1AKMT3', 'EEF1AKMT4', 'EEF1AKNMT', 'EHMT1', 'EHMT2', 'EID1', 'EID2', 'EID2B', 'EIF4A3', 'ELP2', 'ELP3', 'ELP4', 'ELP5', 'ELP6', 'ENY2', 'EP300', 'EP400', 'EPC1', 'EPC2', 'ERBB4', 'ERCC6', 'EXOSC1', 'EXOSC2', 'EXOSC3', 'EXOSC4', 'EXOSC5', 'EXOSC6', 'EXOSC7', 'EXOSC8', 'EXOSC9', 'EYA1', 'EYA2', 'EYA3', 'EYA4', 'EZH1', 'EZH2', 'FAM175A', 'FAM175B', 'FBL', 'FBRS', 'FBRSL1', 'FOXA1', 'FOXO1', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXP4', 'FTO', 'GADD45A', 'GADD45B', 'GADD45G', 'GATAD1', 'GATAD2A', 'GATAD2B', 'GFI1', 'GFI1B', 'GLYR1', 'GSE1', 'GSG2', 'GTF2I', 'GTF3C4', 'HAT1', 'HCFC1', 'HCFC2', 'HDAC1', 'HDAC10', 'HDAC11', 'HDAC2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6', 'HDAC7', 'HDAC8', 'HDAC9', 'HDGF', 'HDGFL2', 'HELLS', 'HIF1AN', 'HINFP', 'HIRA', 'HIRIP3', 'HJURP', 'HLCS', 'HLTF', 'HMG20A', 'HMG20B', 'HMGB1', 'HMGN1', 'HMGN2', 'HMGN3', 'HMGN4', 'HMGN5', 'HNRNPU', 'HNRPL', 'HNRPM', 'HP1BP3', 'HR', 'HSFX3', 'HSPA1A', 'HSPA1A', 'HSPA1B', 'HSPA1B', 'HUWE1', 'IKBKAP', 'IKZF1', 'IKZF3', 'ING1', 'ING2', 'ING3', 'ING4', 'ING5', 'INO80', 'INO80B', 'INO80C', 'INO80D', 'INO80E', 'JADE1', 'JADE2', 'JADE3', 'JAK2', 'JARID2', 'JDP2', 'JMJD1C', 'JMJD6', 'KANSL1', 'KANSL2', 'KANSL3', 'KAT2A', 'KAT2B', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'KAT8', 'KDM1A', 'KDM1B', 'KDM2A', 'KDM2B', 'KDM3A', 'KDM3B', 'KDM4A', 'KDM4B', 'KDM4C', 'KDM4D', 'KDM4E', 'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'KDM6A', 'KDM6B', 'KDM7A', 'KDM8', 'KEAP1', 'KHDRBS1', 'KLF18', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'KMT2E', 'L3MBTL1', 'L3MBTL2', 'L3MBTL3', 'L3MBTL4', 'LAS1L', 'LBR', 'LEO1', 'LRWD1', 'MAGOH', 'MAP3K7', 'MAPKAPK3', 'MASTL', 'MAX', 'MAZ', 'MBD1', 'MBD2', 'MBD3', 'MBD4', 'MBD5', 'MBD6', 'MBIP', 'MBNL1', 'MBNL3', 'MBTD1', 'MCRS1', 'MDC1', 'MEAF6', 'MECP2', 'MEN1', 'METTL11B', 'METTL14', 'METTL16', 'METTL21A', 'METTL3', 'METTL4', 'MGA', 'MGEA5', 'MINA', 'MIS18A', 'MIS18BP1', 'MLLT1', 'MLLT10', 'MLLT6', 'MORF4L1', 'MORF4L2', 'MOV10', 'MPHOSPH8', 'MPND', 'MRGBP', 'MSH6', 'MSL1', 'MSL2', 'MSL3', 'MST1', 'MTA1', 'MTA2', 'MTA3', 'MTF2', 'MUM1', 'MYBBP1A', 'MYO1C', 'MYSM1', 'NAA60', 'NAP1L1', 'NAP1L2', 'NAP1L4', 'NASP', 'NAT10', 'NAT10', 'NBN', 'NCL', 'NCOA1', 'NCOA2', 'NCOA3', 'NCOA6', 'NCOR1', 'NCOR2', 'NEK6', 'NEK9', 'NFRKB', 'NFYB', 'NFYC', 'NIPBL', 'NOC2L', 'NPAS2', 'NPM1', 'NPM2', 'NSD1', 'NSL1', 'NSRP1', 'NSUN2', 'NSUN6', 'NTMT1', 'NUP98', 'OGT', 'OIP5', 'PADI1', 'PADI2', 'PADI3', 'PADI4', 'PAF1', 'PAGR1', 'PAK2', 'PARG', 'PARP1', 'PARP2', 'PARP3', 'PAXIP1', 'PBK', 'PBRM1', 'PCGF1', 'PCGF2', 'PCGF3', 'PCGF5', 'PCGF6', 'PCNA', 'PDP1', 'PELP1', 'PHC1', 'PHC2', 'PHC3', 'PHF1', 'PHF10', 'PHF12', 'PHF13', 'PHF14', 'PHF19', 'PHF2', 'PHF20', 'PHF20L1', 'PHF21A', 'PHF8', 'PHIP', 'PIWIL4', 'PKM', 'PKN1', 'POGZ', 'POLE3', 'PPARGC1A', 'PPM1G', 'PPP2CA', 'PPP4C', 'PPP4R2', 'PQBP1', 'PRDM1', 'PRDM11', 'PRDM12', 'PRDM13', 'PRDM14', 'PRDM16', 'PRDM2', 'PRDM4', 'PRDM5', 'PRDM6', 'PRDM7', 'PRDM8', 'PRDM9', 'PRKAA1', 'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKAG3', 'PRKCA', 'PRKCB', 'PRKCD', 'PRKDC', 'PRMT1', 'PRMT2', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9', 'PRPF31', 'PRR14', 'PSIP1', 'PTBP1', 'PTBP1', 'PUF60', 'RAD51', 'RAD54B', 'RAD54L', 'RAD54L2', 'RAG1', 'RAG2', 'RAI1', 'RARA', 'RB1', 'RBBP4', 'RBBP5', 'RBBP7', 'RBFOX1', 'RBM11', 'RBM15', 'RBM15B', 'RBM17', 'RBM24', 'RBM25', 'RBM4', 'RBM5', 'RBM7', 'RBM8A', 'RBMY1A1', 'RBX1', 'RCC1', 'RCOR1', 'RCOR3', 'REST', 'RFOX1', 'RING1', 'RLIM', 'RMI1', 'RNF168', 'RNF2', 'RNF20', 'RNF40', 'RNF8', 'RNPS1', 'RPS6KA3', 'RPS6KA4', 'RPS6KA5', 'RPUSD3', 'RRP8', 'RSF1', 'RSRC1', 'RUVBL1', 'RUVBL2', 'RYBP', 'SAFB', 'SAP130', 'SAP18', 'SAP25', 'SAP30', 'SAP30L', 'SATB1', 'SATB2', 'SCMH1', 'SCML2', 'SCML4', 'SENP1', 'SENP3', 'SET', 'SETD1A', 'SETD1B', 'SETD2', 'SETD3', 'SETD5', 'SETD6', 'SETD7', 'SETD8', 'SETDB1', 'SETDB2', 'SETMAR', 'SF3B1', 'SF3B3', 'SFMBT1', 'SFMBT2', 'SFPQ', 'SFSWAP', 'SHPRH', 'SIN3A', 'SIN3B', 'SIRT1', 'SIRT2', 'SIRT6', 'SIRT7', 'SKP1', 'SLU7', 'SMARCA1', 'SMARCA2', 'SMARCA4', 'SMARCA5', 'SMARCAD1', 'SMARCAL1', 'SMARCB1', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2', 'SMARCD3', 'SMARCE1', 'SMEK1', 'SMEK2', 'SMYD1', 'SMYD2', 'SMYD3', 'SMYD4', 'SNAI2', 'SP1', 'SP100', 'SP140', 'SPEN', 'SPOP', 'SRCAP', 'SRRM4', 'SRSF1', 'SRSF10', 'SRSF12', 'SRSF3', 'SRSF6', 'SS18L1', 'SS18L2', 'SSRP1', 'STK4', 'SUDS3', 'SUPT16H', 'SUPT3H', 'SUPT6H', 'SUPT7L', 'SUV39H1', 'SUV39H2', 'SUV420H1', 'SUV420H2', 'SUZ12', 'SYNCRIP', 'TADA1', 'TADA2A', 'TADA2B', 'TADA3', 'TAF1', 'TAF10', 'TAF12', 'TAF1L', 'TAF2', 'TAF3', 'TAF4', 'TAF5', 'TAF5L', 'TAF6', 'TAF6L', 'TAF7', 'TAF8', 'TAF9', 'TAF9B', 'TBL1XR1', 'TDG', 'TDRD3', 'TDRD7', 'TDRKH', 'TET1', 'TET2', 'TET3', 'TEX10', 'TFDP1', 'TFPT', 'THRAP3', 'TLE1', 'TLE2', 'TLE4', 'TLK1', 'TLK2', 'TNP1', 'TNP2', 'TONSL', 'TOP2A', 'TOP2B', 'TP53', 'TP53BP1', 'TRA2B', 'TRIM16', 'TRIM24', 'TRIM27', 'TRIM28', 'TRIM33', 'TRRAP', 'TRUB2', 'TSSK6', 'TTK', 'TYW5', 'U2AF2', 'UBE2A', 'UBE2B', 'UBE2D1', 'UBE2D3', 'UBE2E1', 'UBE2H', 'UBE2N', 'UBE2T', 'UBN1', 'UBR2', 'UBR5', 'UBR7', 'UCHL5', 'UHRF1', 'UHRF2', 'UIMC1', 'USP11', 'USP12', 'USP15', 'USP16', 'USP17L2', 'USP21', 'USP22', 'USP3', 'USP36', 'USP44', 'USP46', 'USP49', 'USP7', 'UTY', 'VDR', 'VIRMA', 'VPS72', 'VRK1', 'WAC', 'WDR5', 'WDR77', 'WDR82', 'WHSC1', 'WHSC1L1', 'WSB2', 'WTAP', 'YAF2', 'YEATS2', 'YEATS4', 'YTHDC1', 'YWHAB', 'YWHAE', 'YWHAZ', 'YY1', 'ZBTB16', 'ZBTB33', 'ZBTB7A', 'ZBTB7C', 'ZC3H13', 'ZCWPW1', 'ZFP57', 'ZGPAT', 'ZHX1', 'ZMYM2', 'ZMYM3', 'ZMYND11', 'ZMYND8', 'ZNF217', 'ZNF516', 'ZNF532', 'ZNF541', 'ZNF592', 'ZNF687', 'ZNF711', 'ZNHIT1', 'ZRANB3', 'ZZZ3'],
    'm6a_writers'      : ['METTL3', 'METTL14', 'METTL16', 'KIAA1429','RBM15', 'WTAP'],
    'm6a_erasers'      : ['FTO', 'ALKBH5'],
    'm6a_readers'      : ['YTHS', 'EIF3', 'HNRNPC', 'HNRNPA2B1', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3'],
    'm6a_re_wr_er'     : ['METTL3', 'METTL14', 'METTL16', 'KIAA1429','RBM15', 'WTAP', 'FTO', 'ALKBH5', 'YTHS', 'EIF3', 'HNRNPC', 'HNRNPA2B1', 'YTHDF1', 'YTHDF2', 'YTHDC1', 'YTHDC2', 'TYSND1', 'SND1', 'PRRC2A', 'LRPPRC', 'FMR1','FMR1NB', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3'],
    'PRC2'             : ['EZH1', 'EZH2', 'EED', 'SUZ12', 'RBBP4', 'RBBP7', 'JARID2', 'PCGF1', 'PCGF2', 'RING1', 'BMI1'],
    'Laura'            : ['NAMPT', 'NAPRT', 'IDO', 'DHFR', 'NNMAT1', 'NNMAT2', 'NNMAT3', 'QPRT', 'MAT2A'],
    'Kevin'            : ['PRPF8', 'SRRM1', 'SRRM2', 'ACIN1', 'RNPS1', 'CLK1', 'CLK2', 'CLK3', 'CLK4'],
    'Master'           : [
            #PRC1/2
            'RING1', 'BMI1', 'PCGF1', 'PCGF2', 'EED', 'SUZ12', 'EZH1', 'EZH2', 'RBBP4', 'RBBP7', 'JARID2', 
            #Chromatin stuff
            'CTCF', 'ARID1A', 'PBRM1', 'ARID1B', 'SMARCA2', 'SMARCA4', 'SMARCC1', 'SMARCC2', 'SMARCD1', 'SMARCD2', 'SMARCD3', 'ACTL6A', 'ACTL6B', 'SMARCB1', 'SA1', 'SA2', 'SMC1A', 'SMC3', 'RAD21', 
            #NOTCH1 upregulated
            'MYC', 'HEY1', 'NFKB1', 'NFKB2', 'CCND1', 'CDKN1A', 'CDKN1B', 'NOTCH1', #Proteins of interest in lab below
            #Laura
            'NAMPT', 'NAPRT', 'NNMAT1', 'NNMAT2', 'NNMAT3', 'IDO1', 'QPRT', 'DHFR', 'MAT2A', 'HES1', 'BRD4', 'DTX1', 'FOXO3', 'SESN1', 'U2AF1', 'U2AF2',
            #Involved in monocyte interactions with vascular endothelial repair. From Wtikowski, 2020
            'PECAM1', 'CD44', 'ITGA4', 'CX3CR1', 'TNFSF10', 'CSF1R ', 
            #m6a proteins
            'METTL3', 'METTL14', 'METTL16', 'WTAP', 'RBM15', 'RBM15B', 'KIAA429', 'ZC3HI3', 'CBLL1', 'VIRMA', 'FTO', 'ALKBH5', 'HNRNPC', 'YTHDF1', 'YTHDF2', 'YTHDF3', 'YTHDC1', 'YTHDC2', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'HNRNPA2B1',
            #Frequently mutated in AML
            'DNMT3A', 'TET2', 'ASXL1','TP53', 'KRAS', 'IDH', 
            #Deaminase stuff
            'APOBEC1', 'APOBEC2', 'APOBEC3A','APOBEC3B','APOBEC3C','APOBEC3D','APOBEC3F','APOBEC3G','APOBEC3H','APOBEC4', 'AID', 
            #Mentioned alongside CCR9 in Mansour's paper, expressed in MOLT4 but not healthy cells
            'SD06', 'OR10R2', 'SKAJ', 'DNAH17', 'LZTFL1', 'TCF19', 'MAD2L1', 'UBASH3A', 'LCK', 'ITGAE',
            #Gene is in depmap as having the lowest median score among all cancers for T-ALL
            'AHCY', 'CHERP','CTPS1','DCDC2','H2BC11','NPTXR','RHOH','SUCLG1','TAS1R2','TONSL','TPT1','ZAP70',
            #Alterations in these genes across 15 subtypes of T-ALL
            'BCL11B', 'KMT2A', 'MLLT10','HOXA9','NKX2','CUL1','CDKN2A','PHF6','LEF1','PSIP1','NOL4L','KMT2E','KAT6A','DDX39B','MYL1','CD99','E2F1','GEMIN8','FBXW7', 
            #Miscellaneous
            'SF3B1', 'BRD9', 'NT5C2', 'USP7', 'BCL2L11', 'MYB', 'RUNX1', 'IKZF3', 'KAT2B', 'PTPN4', 'CYLB', 'ITGB1', 'CCND3', 'PSM3IP', 'CREBBP', 'PSMB9', 'SRSF2', 'SRSF3', 'SRSF6', 'SRSF7',
            'SRSF11', 'NUP85', 'PSMD4', 'UPF1', 'KDM6B', 'IDH2','INTS3','SPI1', 'TCF7','DHCR7','HMGCS1','WT1','FLT3','KIT','MDA5','PKR','NCOR2', 'FMR1','MATR3','EP300','LMO1','IRF4','KDM4C', 'IDO',
            'JQ1', 'MXD1', 'MTAP', 'CHEK1', 'RBM39', 'CCR9', 'CCR7', 'CD19', 'CDK6', 'BRD1', 'BRD2', 'BRD3', 'APCDD1', 'FOXO1', 'PRPF39','MSI2', 'NFKIA','SERPINB1','CD69','CCL25', 'DCAF15',
            'TNFSF9', 'TNFRSF9', 'CD38', 'TRBC1', 'TRBC2', 'ZBTB7B', 'RUNX3','ZBP1', 'KMT2D', 'IGF2','TAL1','PTEN', 'IDH1','PIK3CA','PIK3R1','STAT5B','JAK1','JAK2','DNM2','IL7R','NGN3','CASC3', 
            'MEF2C', 'TCF3','PBX1','HLF','ITGAM','AURKB','TLX3','STAG2','SOX11','ID1','ID3','XBP1','CD34','LMO2', 'JAK3', 'RIGI','ADAR','ADARB1','ADARB2','IGHV','FOXM1','AURKA',
            'RELA','TERT','ARNT','HDAC1','HDAC2','HDAC6','HDAC10','NPM1','KDM1B'
        ]
    }


from gseapy import Msigdb
def KTC_GetGeneSet(name_gene_set, db_version='2024.1.Hs'):
    try:
        # First, try to find the gene set in the `gene_sets` dictionary
        print('Made it Here')
        gene_set = gene_sets[name_gene_set]
        print(f'\n -- NOTE -- Using custom set of genes from gene_sets dictionary: {name_gene_set}')
    except KeyError:
        print('Made it there')
        # If not found, try to fetch it from msigdb, prioritizing 'HALLMARK' over others
        categories = ['h.all', 'c2.cp', 'c5.bp', 'c6.all', 'c7.all']
        msigdb = Msigdb()
        gmt_cache = {}
        gene_set = None  # Initialize gene_set in case no match is found

        for category in categories:
            try:
                if category not in gmt_cache:
                    gmt_cache[category] = msigdb.get_gmt(category=category, dbver=db_version)
                gmt = gmt_cache[category]

                if gmt and name_gene_set in gmt:
                    print(f'\n -- NOTE -- Using gene set from msigdb: {name_gene_set} from {category}')
                    gene_set = gmt[name_gene_set]
                    break  # Exit loop once a match is found
            except Exception as e:
                print(f"Error fetching or searching category {category}: {e}")

        if gene_set is None:
            # If no match is found, use string as a gene name for a set with a single gene
            print(f'\n -- NOTE -- Gene set: {name_gene_set} not found in gene_sets or using Msigdb. Defaulting to interpreting {name_gene_set} as a set with a single gene')
            gene_set = [name_gene_set.upper()]
    return gene_set


#%%
def lol():
    return 2
