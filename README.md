# Combining satellite-derived standing deadwood estimates and vegetation modelling to predict climate change impacts on Norway spruce forests in Germany

## Summary
Recent summer droughts (2018–2020) caused unprecedented Norway spruce (Picea abies) mortality in Central Europe, producing large amounts of standing deadwood. Using satellite-derived standing deadwood data, we developed a weather-driven empirical model (GLM) of spruce mortality and integrated it into the dynamic vegetation model LPJ-GUESS. The DGVM reproduces observed mortality patterns and biomass declines in Germany. Future projections (2021–2070) under RCP2.6 and RCP8.5 indicate recurring drought-induced mortality, substantial carbon losses, and increased calamity timber, while harvests of target-diameter timber decline. Without incorporating weather-driven mortality, LPJ-GUESS underestimates drought impacts. These results highlight risks for carbon storage in Norway spruce forests and support a transition to more diverse mixed forests in Germany.

## Overview
This repository contains the code for reproducing the analyses and simulations presented in the publication. Below is an overview of the repository structure and contents:

### Files and Directories

- **`01_model_calibration_satmort.R`**  
  R script for calibrating the GLM predicting weather-related Norway spruce mortality.

- **`02_model_evaluation_satmort.R`**  
  R script for evaluating the performance of the GLM.

 - **`climatic_predictors_2019_2020.zip`**  
  A ZIP file containing the climatic predictor variables and the tree mortality data for the years 2019 and 2020.

- **`lpj-guess_european_applications_tree_species_code_drought_mortality.zip`**  
  A ZIP file containing the LPJ-GUESS source code with the implemented empirical-based drought mortality algorithms.

### Instructions

#### Model Calibration and Evaluation:
- Open `01_model_calibration_satmort.R` and follow the instructions to calibrate the GLM using the provided dataset.
- Use `02_model_evaluation_satmort.R` to evaluate the model performance and reproduce the results in the manuscript.

#### LPJ-GUESS Integration:
- Extract the contents of `lpj-guess_european_applications_tree_species_code_drought_mortality.zip`.
- Follow the documentation within the extracted folder to compile and run LPJ-GUESS with the implemented drought mortality algorithm.

#### Execturing LPJ-GUESS:
- The newly implemented drought mortality feature is compatible exclusively with the cfx input module! To use it, run LPJ-GUESS by providing the instruction file as an argument with the cfx input module. This is because the calculation of predictor variables is partially handled within the cfxinput.cpp module. Execute LPJ-GUESS as follows: ./guess -input cfx instructionfile.ins


## Publication
This repository accompanies the following publication:  


## Citation
If you use this repository, please cite the publication as follows:

