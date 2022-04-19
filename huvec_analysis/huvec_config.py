#configuration file with the different dictionaries for integration across all datasets
#%%
#dict for the SQANTI color schemes used without 
sqanti_colors = {
    "FSM": "#6BAED6",
    #"ISM": "#FC8D59",
    "NIC": "#78C679",
    "NNC": "#EE6A50",
    "Genic Genomic":"#969696",
    "Antisense": "#66C2A4",
    "Fusion":"goldenrod1",
    "Intergenic":"darksalmon",
    "Genic Intron":"#41B6C4",
    "All": "#808080"
}

sqanti_protein_colors = {
    'pFSM' : '#1974b5',
    'pNIC' : '#4a943e',
    'pNNC' : '#d94c2e' 
}

#colors used for the different databases
database_colors = {
    'GENCODE' : 'darkblue',
    'UniProt' : 'darkgreen',
    'PacBio' : 'darkred', 
    'PacBio Rescue Resolve' : 'red'   
}

#color scheme used for the classidfication on the short transcripts
transcript_short_classification = {
    'novel_not_in_catalog': 'NNC', 
    'full-splice_match': 'FSM', 
    'novel_in_catalog' : 'NIC',
    'incomplete-splice_match' : 'ISM', 
    'fusion' : 'fusion',
}

# matplotlib font parameters
font = {'family' : 'sans-serif',
        'sans-serif':['Arial'],
        'weight' : 'normal',
        'size'   : 14}

#info related to the pipeline directoties and where the output can be stored too 
#the pipeline  results direcotry
PIPELINE_RESULTS_DIRECTORY='/Volumes/sheynkman/projects/huvec_proteogenomics/huvec'
#PIPELINE_RESULTS_DIRECTORY=  '/Users/gloriasheynkman/Documents/research_drive/projects/huvec_proteogenomics/HUVEC-Proteogenomics/data/huvec'
# PIPELINE_RESULTS_DIRECTORY= '/Users/madison/Documents/Sheynkman_lab/Huvec-project/huvec-proteogenomics/HUVEC_manuscript_analysis/huvec'

#this holds the gtf files and all of the info from gencode 
REFERENCE_DIRECTORY= '/Volumes/sheynkman/projects/huvec_proteogenomics/input'
MZML_DIRECTORY = '/Volumes/sheynkman/ms/ms_data/210831_huvec_fractions_tryp_HCDonly/mzml'
# REFERENCE_DIRECTORY= '/Users/madison/Documents/Sheynkman_lab/Huvec-project/huvec-proteogenomics/HUVEC_manuscript_analysis/input'
#REFERENCE_DIRECTORY= '/Users/gloriasheynkman/Documents/research_drive/projects/huvec_proteogenomics/HUVEC-Proteogenomics/input'
EXPERIMENT_NAME='huvec'
# %%
