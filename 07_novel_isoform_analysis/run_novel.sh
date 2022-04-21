echo "Generating novel peptide table"
python generate_novel_peptide_dataset.py

echo "Generating novel peptide excel data with spectra plots"
python novel_visualization.py

echo "Adding EC gene info"
python novel_pep_table_analysis.py