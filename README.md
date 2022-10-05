# Protodeboronation prediction & DFT utilities

This repository contains the code for the paper: "Quantitative in silico prediction of the rate of protodeboronation by a mechanistic DFT-aided algorithm"

## Data
This folder contains all relevant data to reproduce the results & predictions presented in the paper above, including DFT energy calculations, pKa and pKaH values, and an overview of the mechanistic pathways active for all of the 100 molecules under consideration (determined with hueristics).

## Models
For the sake of reproducibility, we have also included a save of all the linear models trained for this work. Do note that training a linear model on a few dozen datapoints takes less than a second, so this may as well be done in real time

## PDB_predictions
The PDB_predictions.py contains all functions required to run the protodeboronation prediction workflow. For an example of how to use these functions, please consider the PDB_predictions.ipynb file.

## DFT_utilities
The DFT_utilities.py contains useful functions for performing DFT calculations. For an example of how to use these functions, please consider the DFT_utilities.ipynb file. These functions include:
 - Automatically generate SLURM files
 - Extract coordinates from a .log output file
 - Extract the optimised energy from a .log output file

For any questions please feel free to open an issue.
