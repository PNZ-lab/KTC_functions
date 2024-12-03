# KTC_functions.py
This script serves as a centralized hub for functions and data that other scripts share.

To use this script from another, primary script it needs to be accessible to the primary script. There are several ways to ensure this but having KTC_functions.py in your PYTHONPATH or simply having KTC_functions.py in the same directory as the primary script will serve.
## KTC_pos_to_gene(path_GTF, list_pos)
This function takes as input a path to a GTF file as well as a list of tuples for positions (e.g. [(chr1,12345),(chr2,23456)]).

It then returns gene names for the given posiitions in the given gtf file.

## KTC_GetGeneSet(name_gene_set)
This script takes as input a string and attempts to convert it to a set of genes. It will return a list of genes when it succeeds, trying in order to:

1. Find a value in the gene_sets dictionary in KTC_functions. These are pre-defined lists of genes (e.g. 'Laura' and 'splicing factors' would succeed)
2. Search Msigdb for the string provided to return a set of genes (e.g. 'HALLMARK_MYC_TARGETS_V1' would succeed)
3. If the two attempts above fail, the function will assume that the input string is a single gene (e.g. 'MYC' will return the gene list: ['MYC'])
