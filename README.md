# KTC_functions.py
This script serves as a centralized hub for functions and data that other scripts share.

To use this script from another, primary script it needs to be accessible to the primary script. There are several ways to ensure this but having KTC_functions.py in your PYTHONPATH or simply having KTC_functions.py in the same directory as the primary script will serve.
## KTC_pos_to_gene(path_GTF, list_pos)
This function takes as input a path to a GTF file as well as a list of tuples for positions (e.g. [(chr1,12345),(chr2,23456)]).

It then returns gene names for the given posiitions in the given gtf file.

## KTC_GetGeneSet(name_gene_set)
This script takes as input a string and attempts to convert it to a set of genes. It will return a list of genes when it succeeds, trying in order to:

1. Check if the input is a list. If so, it will interpret it as a set of genes. e.g:
   - KTC_GetGeneSet(['MYC', 'TAL1', 'RBM39'])
2. Locate a predefined set of genes by a name in the dictionary gene_sets in KTC_functions.py, e.g:
   - KTC_GetGeneSet('Laura')
3. Search Msigdb for a set of genes with that name. Default database is human (2024.1.Hs)), e.g:
   - KTC_GetGeneSet('HALLMARK_MYC_TARGETS_V1') # Find names here: https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
   - KTC_GetGeneSet('HALLMARK_MYC_TARGETS_V1', db_version='2024.1.Hs') # Funtionally identical to the above
   - KTC_GetGeneSet('HALLMARK_MYC_TARGETS_V1', db_version='2024.1.Mm') # Searching the mouse equivalent
4. If all of the above fail, it defaults to interpret the string inout as a single gene, e.g:
   - KTC_GetGeneSet('MYC')
    
