# for the trypsin hcd only data to get the pecam isoforms 
# Code to analyze the novel peptides that have been detected from the MM searches 
# Performs MM search with the pb derived dataset aganist GENCODE and Uniprot dataset 
# Novel pep analysis needs performed aganist all of the enzymatic digests 
#%%
# read in the needed packages 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from Bio import SeqIO
import os
from huvec_analysis import huvec_config # paramters for plotting

# Establish the right fonts for matplotlib throughout
matplotlib.rc('font', **huvec_config.font)

# All the plots go into a folder 
plot_dir = 'plot'
if not os.path.exists(plot_dir):
   os.makedirs(plot_dir)
# all of the stats will go into a directory
stats_dir = 'stats'
if not os.path.exists(stats_dir):
   os.makedirs(stats_dir)

# Read in the data paths
pacbio_hybrid_peps_path = '../00_pre_analysis/metamorpheus_table/AllPeptides.PacBioHybrid.tsv'

# fasta sequences 
gencode_fasta_path = f'{huvec_config.PIPELINE_RESULTS_DIRECTORY}/gencode_db/gencode_protein.fasta'
uniprot_fasta_path = f'{huvec_config.REFERENCE_DIRECTORY}/uniprot_reviewed_canonical_and_isoform.fasta'

def read_peps_file(filename):
    peps_file = pd.read_table(filename)
    peps_file = peps_file[['Scan Number', 'File Name', 'Gene Name', 'Protein Accession', 'Base Sequence', 'Score', 'QValue', 'gene_name', 'ec_priority']]
    peps_file.columns = ['scan_num', 'file_name', 'gene', 'acc', 'seq', 'score', 'qval', 'gene_name', 'ec_priority']
    return peps_file

# Read in the trypsin peps, filter and keep the columns we need
pacbio_hybrid_peps = read_peps_file(pacbio_hybrid_peps_path)

# Filter for peptide identifications mapped up to PB accessions (AllPeptides file from MM)
# For the pacbio hybrid databases
pacbio_pb_peps = pacbio_hybrid_peps[pacbio_hybrid_peps['acc'].str.startswith('PB')]
pacbio_sample_peps = pacbio_pb_peps['seq'].to_list()

# Import all of the sequencs from gencode to a big string (represents the full concatenated proteome)
# (for search against potential novel peptides)
gencode_seq = [str(rec.seq) for rec in SeqIO.parse(gencode_fasta_path, 'fasta')]
gencode_all_sequences = ','.join(gencode_seq)

# Import all of the sequencs from uniprot to a big string 
uniprot_seq = [str(rec.seq) for rec in SeqIO.parse(uniprot_fasta_path, 'fasta')]
uniprot_all_sequences = ','.join(uniprot_seq)


# Find the genes that are ec relevant for all of the digests

# Finding and writing out the pacbio-specific novel peptides, for each protease digest
#find the novel peptides 
novel_peps = set()
novel_peps_to_gencode = set()
novel_peps_to_uniprot = set()

for pep in pacbio_sample_peps:
    # some peptides are indistinguishable (have I/L), take first one
    if '|' in pep:
        pep = pep.split('|')[0]
    # is the base peptide sequence in the ref database
    if pep not in gencode_all_sequences:
        novel_peps_to_gencode.add(pep)
    if pep not in uniprot_all_sequences:
        novel_peps_to_uniprot.add(pep)
novel_peps = novel_peps_to_gencode.intersection(novel_peps_to_uniprot)

#write out the novel pacbio accession and gencode for each of the novel peptides 
peps_novel = pacbio_pb_peps[pacbio_pb_peps['seq'].isin(novel_peps)]
peps_novel.to_csv('stats/pacbio_novel_peptides.tsv', sep='\t', index=None)

peps_novel_to_gencode = pacbio_pb_peps[pacbio_pb_peps['seq'].isin(novel_peps_to_gencode)]
peps_novel_to_gencode.to_csv('stats/pacbio_novel_peptides_to_gencode.tsv', sep='\t', index=None)


peps_novel_to_uniprot= pacbio_pb_peps[pacbio_pb_peps['seq'].isin(novel_peps_to_uniprot)]
peps_novel_to_uniprot.to_csv('stats/pacbio_novel_peptides_to_uniprot.tsv', sep='\t', index=None)

#%%
import numpy as np
import matplotlib.pyplot as plt
y = np.array([110, 25, 4, 1])
mylabels = ["PacBio-Hybrid", "AspN", "Chmotrypsin", "GluC"]

plt.pie(y, labels = mylabels, autopct='%1.0f%%')
plt.title("Number of novel peptides by digest")
plt.savefig('plot/novel_peps_recovered_by_digest.pdf', bbox_inches= 'tight')
# %%
