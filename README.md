# KTC_functions.py
This script serves as a centralized hub for functions and data that other scripts share.

To use this script from another, primary script it needs to be accessible to the primary script. There are several ways to ensure this but having KTC_functions.py in your PYTHONPATH or simply having KTC_functions.py in the same directory as the primary script will serve.
## KTC_pos_to_gene(path_GTF, list_pos)
This function takes as input a path to a GTF file as well as a list of tuples for positions (e.g. [(chr1,12345),(chr2,23456)]).

It then returns gene names for the given positions in the given gtf file.

## KTC_GetGeneSet(name_gene_set)
This script takes as input a string and attempts to convert it to a set of genes. It will return a list of genes when it succeeds, trying in order to:

1. Check if the input is a list. If so, it will interpret it as a set of genes. e.g:
   - KTC_GetGeneSet(['MYC', 'TAL1', 'RBM39'])
2. Locate a predefined set of genes by a name in the dictionary gene_sets in KTC_functions.py, e.g:
   - KTC_GetGeneSet('Laura')
3. Search Msigdb for a set of genes with that name [1]. Default database is human (2024.1.Hs)), e.g:
   - KTC_GetGeneSet('HALLMARK_MYC_TARGETS_V1')
   - KTC_GetGeneSet('HALLMARK_MYC_TARGETS_V1', db_version='2024.1.Hs') # Funtionally identical to the above
   - KTC_GetGeneSet('HALLMARK_MYC_TARGETS_V1', db_version='2024.1.Mm') # Searching the mouse equivalent
4. If all of the above fail, it defaults to interpreting the string input as a single gene, e.g:
   - KTC_GetGeneSet('MYC')

[1] : Find names of gene sets here: https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
    
## KTC_Splice_Mapper(gtf_file, rmats_file, event_id, splicing_event_type)
This script takes as input a gtf file and the output of KTC_rmats_compiler together with the event ID and the splicing event type for an event in the rMATS file. It will return a graph that highlights the positions denoted in the rMATS file on a map of the locations of introns and exons for the relevant gene.

Example:
```
KTC_Splice_Mapper('/Users/kasperthorhaugechristensen/Downloads/Homo_sapiens.GRCh38.110.chr.gtf', '/Users/kasperthorhaugechristensen/Downloads/KO1_rMATS_compiled.tsv', 80549, 'SE')
```
Yields:

![image](https://github.com/user-attachments/assets/56b57ca9-7e6e-46ba-a2b7-6d8d0985f771)


If several maps are to be generated you can perform the parsing of the gtf a single time using KTC_preParse_gtf and then use that as input in KTC_Splice_Mapper to save time (parsing takes ~15s). KTC_Splice_Mapper will register that it is being fed a parsed gtf file in place of an unparsed gtf file and skip the parsing.

Example:
```
preparsed_gtf = KTC_preParse_gtf('/Users/kasperthorhaugechristensen/Downloads/Mus_musculus.GRCm39.110.chr.gtf')
KTC_Splice_Mapper(preparsed_gtf, '/Volumes/cmgg_pnlab/Kasper/Analyses/Joao/2025_TF_analysis/CD2_v_Vav_rMATS_compiled.tsv', 53351, 'SE')
KTC_Splice_Mapper(preparsed_gtf, '/Volumes/cmgg_pnlab/Kasper/Analyses/Joao/2025_TF_analysis/CD2_v_Vav_rMATS_compiled.tsv', 53351, 'SE')
```
Yields:
![image](https://github.com/user-attachments/assets/0cd318a9-f5e7-4830-8dd3-718a20abdbf2)
![image](https://github.com/user-attachments/assets/0855e3e6-283c-42f5-bcd9-79e0831fb016)
