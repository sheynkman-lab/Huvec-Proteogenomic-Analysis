DIGEST=42 # Trypsin
INPUT_FASTA=/Volumes/sheynkman/projects/huvec_proteogenomics/huvec/hybrid_protein_database/huvec_hybrid.fasta

rpg -i $INPUT_FASTA -o hybrid_theoretical_peptides.tsv -f tsv -e $DIGEST