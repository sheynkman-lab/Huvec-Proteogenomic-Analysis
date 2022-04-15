# HUVEC Proteogenomic pipeline analysis

## Setting up
To set up the project for analysis first go to huvec_analysis/huvec_config.py and update the 
PIPELINE_RESULTS_DIRECTORY and REFERENCE_DIRECTORY to point to the appropriate data directories

Then run the following :

```
pip install -e .
```

This will create a path to the current directory, which will allow analysis notebooks to pull from 
the huvec_analysis package

## Running the Analysis
To perform the analysis, first download the results of the Long-Read Proteogenomics pipeline

Open `huvec_analysis/huvec_config.py` and update the variables `PIPELINE_RESULTS_DIRECTORY` and `REFERENCE_DIRECTORY` to point to the ouptut of the huvec pipeline and reference tables, respectively. 

Open the python notebooks of interest and run interactively to generate plots and statistics of each analysis