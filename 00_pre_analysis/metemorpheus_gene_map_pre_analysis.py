# pre anslysis code for the metamorpheus results on trypsin HCD run 
#%%
import os, sys 
import re
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
import pandas as pd 
import numpy as np
from Bio import SeqIO
from huvec_analysis import huvec_config # parameters for plotting

# all intermediate tables go into their own folder 
table_dir = 'non_calib_MM_tables'
if not os.path.exists(table_dir):
    os.makedirs(table_dir)

#######################
### File locations ###
######################
# file locations for the MM data from both the GENCODE and the PacBio MM searched dataset
# Peptide files 
sample_name = 'huvec'
gene_isoname_file = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/reference_tables/gene_isoname.tsv'
pb_gene_file = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/protein_classification/{sample_name}_genes.tsv'
uniprot_fasta_file = uniprot_fasta_file = f'{huvec_config.REFERENCE_DIRECTORY}/uniprot_reviewed_canonical_and_isoform.fasta'
# Peptides files 
gencode_peps_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/gencode/search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv'
pacbio_hybrid_peps_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllPeptides.huvec.hybrid.psmtsv'
uniprot_peps_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/uniprot/search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv'

# PSM file locations 
gencode_psms_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/gencode/search_results/Task1SearchTask/AllPSMs.psmtsv'
pacbio_psms_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllPSMs.psmtsv'
uniprot_psms_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/uniprot/search_results/Task1SearchTask/AllPSMs.psmtsv'

# Protein groups file location 
# MM searches with the different databases 
gencode_protein_grps_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/gencode/search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv'
pacbio_protein_grps_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllQuantifiedProteinGroups.huvec.hybrid.tsv'
uniprot_protein_grps_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/uniprot/search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv'

# EC gene list 
human_ec_gene_path = f'{huvec_config.REFERENCE_DIRECTORY}/human_ec_genes_from_karen.txt'
human_EC_genes =pd.read_table(human_ec_gene_path, header=None)[0].to_list()
# Function defintions for reading in the tables and establishing the gene maps 
def read_metamorpheus_file(metamorpheus_file, gene_map):
    metamorpheus_table = pd.read_table(metamorpheus_file, usecols= ['File Name', 'Gene Name', 'Protein Accession', 'Base Sequence', 'Full Sequence', 'PSM Count (unambiguous, <0.01 q-value)', 'Decoy/Contaminant/Target','Score', 'QValue', 'Scan Number', 'Precursor Charge', 'Matched Ion Mass-To-Charge Ratios', 'Matched Ion Intensities', 'Total Ion Current', 'Mass Diff (ppm)'])
    metamorpheus_table = metamorpheus_table[metamorpheus_table['Decoy/Contaminant/Target']=='T']
    metamorpheus_table = metamorpheus_table[metamorpheus_table['QValue']<=0.01]
    metamorpheus_table['accs'] = metamorpheus_table['Protein Accession'].str.split('|')
    metamorpheus_table['genes'] = metamorpheus_table['accs'].apply(lambda accs: [gene_map[x] for x in accs])
    metamorpheus_table['gene_name'] = metamorpheus_table['genes'].str.get(0)
    metamorpheus_table['ec_priority'] =  metamorpheus_table['gene_name'].isin(human_EC_genes)*1
    return metamorpheus_table


def read_peptide_file(peptides_file, gene_map):
    return read_metamorpheus_file(peptides_file, gene_map)

def read_psm_file(psm_file, gene_map):
    return read_metamorpheus_file(psm_file, gene_map)

def get_genes_in_protein_group(protein_group, gene_map):
    genes = []
    for acc in protein_group.split('|'):
        if acc in gene_map.keys():
            genes.append(gene_map[acc])
        return genes

def read_protein_group_file(protein_groups_filename, gene_map):
    protein_groups = pd.read_table(protein_groups_filename, index_col=False)
          #  usecols=['Protein Accession', 'Gene', 'Number of Proteins in Group', 'Unique Peptides', 'Shared Peptides','Number of Peptides', 'Protein Decoy/Contaminant/Target', 'Protein QValue'])
    protein_groups = protein_groups[protein_groups['Protein QValue']<= 0.01]
    protein_groups = protein_groups[protein_groups['Protein Decoy/Contaminant/Target'] =='T']
    protein_groups['accs'] = protein_groups['Protein Accession'].str.split('|')
    protein_groups['genes'] = protein_groups['Protein Accession'].apply(lambda pgroup: get_genes_in_protein_group(pgroup, gene_map))
    protein_groups['gene_name'] = protein_groups['genes'].str.get(0)
    # Add a 1 if the gene is considered EC priotity # TODO why is this failing
    #protein_groups['ec_priority'] = protein_groups['genes'].isin(human_EC_genes) * 1
    return protein_groups


#***************************************************
# gene maps
#***************************************************
gencode_gene_map = pd.read_table(gene_isoname_file, names=['gene','isoname'])
gencode_gene_map = pd.Series(gencode_gene_map.gene.values, index=gencode_gene_map.isoname).to_dict()
pacbio_gene_map = pd.read_table(pb_gene_file)
pacbio_gene_map = pd.Series(pacbio_gene_map.pr_gene.values, index=pacbio_gene_map.pb).to_dict()
hybrid_gene_map = {**pacbio_gene_map, **gencode_gene_map}

uniprot_gene_map = {}
for record in SeqIO.parse(uniprot_fasta_file, format='fasta'):
    isoform = record.id.split('|')[1]
    gene_regex = re.search(r'GN=(.*)', record.description)
    if gene_regex is not None:
        gene = gene_regex.group(1).split(' ')[0]
        uniprot_gene_map[isoform] = gene
    else:
        uniprot_gene_map[isoform] = None

# Read in the peptide files 
pacbio_hybrid_peps = read_peptide_file(pacbio_hybrid_peps_path, hybrid_gene_map)
gencode_peptides = read_peptide_file(gencode_peps_path, gencode_gene_map)
uniprot_peptides = read_peptide_file(uniprot_peps_path, uniprot_gene_map)

# Read in the PSMs file 
pacbio_hybrid_psm = read_psm_file(pacbio_psms_path, hybrid_gene_map)
gencode_psm = read_psm_file(gencode_psms_path, gencode_gene_map)
uniprot_psm = read_psm_file(uniprot_psms_path, uniprot_gene_map)

# Read in the protein groups file 
pacbio_hybrid_protein_grp = read_protein_group_file(pacbio_protein_grps_path, hybrid_gene_map)
gencode_protein_grps = read_protein_group_file(gencode_protein_grps_path, gencode_gene_map)
uniprot_protein_grps = read_protein_group_file(uniprot_protein_grps_path, uniprot_gene_map)

# output these files to CSVs for easy access
pacbio_hybrid_peps.to_csv('220330_metamorpheus_table/AllPeptides.PacBioHybrid.tsv', sep = '\t', index= None)
gencode_peptides.to_csv('220330_metamorpheus_table/AllPeptides.GENCODE.tsv', sep = '\t', index= None )
uniprot_peptides.to_csv('220330_metamorpheus_table/AllPeptides.Uniprot.tsv', sep = '\t', index= None)

pacbio_hybrid_psm.to_csv('220330_metamorpheus_table/AllPSMs.PacBioHybrid.tsv', sep = '\t', index=None)
gencode_psm.to_csv('220330_metamorpheus_table/AllPSMs.GENCODE.tsv', sep = '\t', index=None)
uniprot_psm.to_csv('220330_metamorpheus_table/AllPSMs.UniProt.tsv', sep= '\t', index=None)

pacbio_hybrid_protein_grp.to_csv('220330_metamorpheus_table/AllProteinGroups.PacBioHybrid.tsv', sep = '\t', index=None)
gencode_protein_grps.to_csv('220330_metamorpheus_table/AllProteinGrouos.GENCODE.tsv', sep = '\t', index= None)
uniprot_protein_grps.to_csv('220330_metamorpheus_table/AllProteinGroups.UniProt.tsv', sep = '\t', index= None)
#%%

