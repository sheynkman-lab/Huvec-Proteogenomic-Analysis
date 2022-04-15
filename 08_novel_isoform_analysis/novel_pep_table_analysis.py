# file to read in the novel peptides excel sheet add the ec significance 
# also deetrmine how many of the novel peps lead to single mapping vs indirect mapping 
#%%

import pandas as pd 
import os
from huvec_analysis import huvec_config # parameters for plotting

# all intermediate tables go into their own folder 
table_dir = 'metamorpheus_table'
if not os.path.exists(table_dir):
    os.makedirs(table_dir)

# read in the novel pep file 
novel_pep_file_path = '../08_novel_isoform_analysis/novel_peptides_220119.xlsx'

# EC gene list 
human_ec_gene_path = f'{huvec_config.REFERENCE_DIRECTORY}/human_ec_genes_from_karen.txt'
human_EC_genes =pd.read_table(human_ec_gene_path, header=None)[0].to_list() 
novel_peps_table = pd.read_excel(novel_pep_file_path)
#%%
# split the gene name on the : 
novel_peps_table['gene_name'] = novel_peps_table['Gene'].str.split(':').str[1]
novel_peps_table['ec_priority'] = novel_peps_table['gene_name'].isin(human_EC_genes)*1
# find how many of the pb accessions are multi mapping 
sub = '|'
novel_peps_table['unique_isoform'] = novel_peps_table['Accession'].str.find(sub)
novel_peps_table = novel_peps_table[['Filename', 'Scan Number', 'Base Sequence', 'Full Sequence',
       'Accession', 'Gene', 'Score', 'Precursor Charge',
       'PSM Count (unambiguous, <0.01 q-value)', 'Mass Diff (ppm)', 'gene_name',
       'ec_priority', 'unique_isoform']]
novel_peps_table.to_excel('stats/novel_peps_filtered.xlsx')
# %%
