# Novel Peptides

# Novel peptides found in Mass spectrometry, but are not found in the GENCODE or UNIPROT databases

# This analysis pulls the spectra info for what scans map to the novel peptide

from huvec_analysis import huvec_config
import pandas as pd
import os

novel_peptides_ifile = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/novel_peptides/huvec_hybrid.pacbio_novel_peptides.tsv'
peptides_ifile = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllPeptides.huvec.hybrid.psmtsv'
abundance_ifile = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/hybrid_protein_database/huvec_refined_high_confidence.tsv'
human_ec_gene_ifile = f'{huvec_config.REFERENCE_DIRECTORY}/human_ec_genes_from_karen.txt'
if not os.path.exists('./stats'):
    os.mkdir('./stats')

novel_peptides = pd.read_table(novel_peptides_ifile)
peptides = pd.read_table(peptides_ifile)
abundance = pd.read_table(abundance_ifile)
human_EC_genes =pd.read_table(human_ec_gene_ifile, header=None)[0].to_list() 

novel_peptides = novel_peptides.merge(peptides, left_on='seq', right_on='Base Sequence', how = 'left')
novel_peptides = novel_peptides.merge(abundance, left_on='acc', right_on='base_acc', how = 'left')
novel_peptides['is_ec_gene'] = novel_peptides['pr_gene'].isin(human_EC_genes) * 1

novel_peptides.to_csv('./stats/pacbio_hybrid_novel_peptides.tsv', sep='\t', index=False)
