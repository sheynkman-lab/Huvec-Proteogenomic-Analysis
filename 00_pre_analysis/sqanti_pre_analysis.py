# Pre-analysis on the SQANTI input tables to get the configure and make an intermediate table 
#%%
#establish where the file directories are establishe directories
# current_dir = os.path.dirname(os.path.realpath(__file__))
# parent_dir = os.path.dirname(current_dir)
# sys.path.append(parent_dir)
import pandas as pd 
import numpy as np
import os
from huvec_analysis import huvec_config # parameters for plotting

# folder to hold the intermediate SQANTI table 
sqanti_dir = 'sqanti_info'
if not os.path.exists(sqanti_dir):
    os.makedirs(sqanti_dir)

# file paths to the data 
sqanti_info_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/transcriptome_summary/sqanti_isoform_info.tsv'
human_ec_gene_path = f'{huvec_config.REFERENCE_DIRECTORY}/human_ec_genes_from_karen.txt'
pandey_genes_path = f'{huvec_config.REFERENCE_DIRECTORY}/pandey_upreg_genes.txt'
gencode_human_pc_trans_path = f'{huvec_config.REFERENCE_DIRECTORY}/gencode.v38.pc_transcripts.fa'
gencode_human_all_trans_path = f'{huvec_config.REFERENCE_DIRECTORY}/gencode.v38.transcripts.fa'
protein_classification_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/protein_classification/huvec_unfiltered.protein_classification.tsv'
refined_metadata_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/protein_gene_rename/huvec_orf_refined_gene_update.tsv'

# Read in HUVEC isoform information (SQANTI) into pandas table and add in annotations.
sqanti_info = pd.read_table(sqanti_info_path)
sqanti_info['log2cpm'] = np.log2(sqanti_info['cpm'] + 1)
# add the fractional abudance of the isoform 
gene_cpms = sqanti_info.groupby('gene')['cpm'].sum().reset_index()
gene_cpms.columns = ['gene', 'total_gene_cpm'] 
# merge
sqanti_info = pd.merge(sqanti_info, gene_cpms, how = 'inner', on = 'gene')

sqanti_info['fractional_abundance'] = sqanti_info['cpm']/sqanti_info['total_gene_cpm']
protein_classification = pd.read_table(protein_classification_path)
refined_data = pd.read_table(refined_metadata_path, usecols= ['pb_accs', 'base_acc'])

# Map the base accession on the pb accession in order to get the protein classification
base_acc_map = {}
for _, row in refined_data.iterrows():
    for acc in row['pb_accs'].split('|'):
        base_acc_map[acc] = row['base_acc']
sqanti_info['base_acc'] = sqanti_info['pb_acc'].map(base_acc_map)
sqanti_info = sqanti_info.dropna(subset=['base_acc'])
sqanti_info = sqanti_info.merge(protein_classification[['pb', 'protein_classification_base']], left_on= 'base_acc', right_on= 'pb', how = 'left')

# Add a column indicating if the gene is endothelial-associated (from Hirschi lab)
human_ec_genes = pd.read_table(human_ec_gene_path, header = None)[0].to_list()
sqanti_info['ec_priority'] = sqanti_info['gene'].isin(human_ec_genes) * 1

# Add a column to indicate if the Pandey upregulated genes are also within our dataset 
# pandey_upreg_genes = pd.read_table(pandey_genes_path, header = None)[0].to_list()
# sqanti_info['pandey_upreg_gene'] = sqanti_info['gene'].isin(pandey_upreg_genes) * 1

# Output the intermediate table 
sqanti_info = sqanti_info[sqanti_info.cat !='ISM']
# sqanti_info = sqanti_info[['pb_acc', 'len', 'cat', 'gene', 'transcript', 'cat2', 'fl_cts', 'cpm', 'log2cpm', 'fractional_abundance', 'protein_classification_base', 'ec_priority', 'pandey_upreg_gene']]
sqanti_info = sqanti_info[['pb_acc', 'len', 'cat', 'gene', 'transcript', 'cat2', 'fl_cts', 'cpm', 'log2cpm', 'fractional_abundance', 'protein_classification_base', 'ec_priority']]

sqanti_info.to_csv(f'{sqanti_dir}/sqanti_info.tsv', sep = '\t', index= None)

#%%
