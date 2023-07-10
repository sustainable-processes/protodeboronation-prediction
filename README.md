[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://protodeboronation-prediction.streamlit.app/)
<a target="_blank" href="https://colab.research.google.com/github/sustainable-processes/protodeboronation-prediction/blob/main/PDB-prediction.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>


Visit our interactive website [here!](https://protodeboronation-prediction.streamlit.app/)

# Protodeboronation prediction & DFT utilities

This repository contains the code for the paper: "Quantitative in silico prediction of the rate of protodeboronation by a mechanistic DFT-aided algorithm"

## Data
This folder contains all relevant data to reproduce the results & predictions presented in the paper above, including the DFT generated output files (containing initial geometry, optimised geometry, and the energy of the optimised geometry), pKa and pKaH values, and an overview of the mechanistic pathways active for all of the 100 molecules under consideration (determined with hueristics). Note that the DFT energies listed in the csv files in data/Cox-molecules came from the DFT output files which can be found in data/Cox-molecules/DFT (and equivalently for Novel-molecules).

## Models
For the sake of reproducibility, we have also included a save of all the linear models trained for this work. Do note that training a linear model on a few dozen datapoints takes less than a second, so this may as well be done in real time.

## PDB_predictions
The PDB_predictions.py contains all functions required to run the protodeboronation prediction workflow. For an example of how to use these functions, please consider the PDB_predictions.ipynb file.

## DFT_utilities
The DFT_utilities.py contains useful functions for performing DFT calculations. For an example of how to use these functions, please consider the DFT_utilities.ipynb file. These functions include:
 - Automatically generate SLURM files
 - Extract coordinates from a .log output file
 - Extract the optimised energy from a .log output file

For any questions please feel free to open an issue.

## Figure_generation
The figure_generation.ipynb file contains the code used to generate the figures found in the paper and SI.

## Citing

If you find this project useful, we encourage you to

* Star this repository :star: 
* Cite our [paper](https://doi.org/10.1021/acs.jpca.2c08250).
```
@article{wigh_quantitative_2023,
	title = {Quantitative {In} {Silico} {Prediction} of the {Rate} of {Protodeboronation} by a {Mechanistic} {Density} {Functional} {Theory}-{Aided} {Algorithm}},
	issn = {1089-5639},
	url = {https://doi.org/10.1021/acs.jpca.2c08250},
	doi = {10.1021/acs.jpca.2c08250},
	urldate = {2023-03-19},
	journal = {The Journal of Physical Chemistry A},
	author = {Wigh, Daniel S. and Tissot, Matthieu and Pasau, Patrick and Goodman, Jonathan M. and Lapkin, Alexei A.},
	month = mar,
	year = {2023},
}
```

