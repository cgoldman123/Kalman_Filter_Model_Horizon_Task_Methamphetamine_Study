# Statistical Analysis

This repository contains the R script and data used for statistical analysis in the manuscript, "Individuals with Methamphetamine Use Disorder Show Reduced Directed Exploration and Learning Rates Independent of an Aversive Interoceptive State Induction"

## File Descriptions

### 1. `horizon_MUD_data.csv`
This is the data file containing participants' fitted parameter estimates and psychological measurements used in statistical analyses.

### 2. `horizon_MUD_data_dictionary.xlsx`
This is a dictionary describing each column of horizon_MUD_data.csv.

### 3. `horizon_MUD_code.Rmd`
This is the R script used to perform statistical analyses.


# Behavioral Modeling Pipeline

This repository contains MATLAB scripts for processing behavioral data, merging task results, and fitting computational models for hierarchical analysis. The pipeline supports tasks run under loaded and unloaded conditions and uses Kalman Filter-based modeling.

## File Descriptions

### 1. `main_script`
This is the primary script that orchestrates the entire pipeline. It:
- Initializes the environment (`rng`, `dbstop`, and paths depending on the operating system).
- Reads the `group_list.csv` file to identify participants and task configurations.
- Merges behavioral data across participants using the `merge_horizon_data` function.
- Saves the merged data to `./modeling_results` with a timestamped filename (e.g., `beh_<timestamp>_HC_iMUD.csv`).
- Calls the `fit_model` function to fit computational models to the behavioral data and saves the results (e.g., `HC_MUD_fits_<timestamp>.csv`).

### 2. `merge_horizon_data`
This function merges behavioral data across participants and task conditions. Key features:
- Accepts paths, subject IDs, group data, and task conditions (`loaded` or `unloaded`) as inputs.
- Calls `compile_data` for each subject and condition to create a consistent format for analysis.
- Filters valid rows based on game length (80 trials) and cumulative trial counts (600).
- Outputs a combined data table and a mapping of subject IDs to indices for later use.

### 3. `compile_data`
This function prepares individual subject data for merging. Key features:
- Identifies task counterbalance conditions (`loaded` or `unloaded`) from the group data.
- Reads individual task data files (`*_events.tsv`) for each subject and condition.
- Creates a structured table for each game, including trial details, rewards, and reaction times.
- Outputs a table formatted for downstream merging.

### 4. `fit_model`
This function performs hierarchical Bayesian modeling on the merged behavioral data. Key features:
- Loads and preprocesses the formatted behavioral data file.
- Structures data for JAGS analysis, including:
  - Forced choice decisions
  - Rewards
  - Game lengths
  - Uncertainty conditions
- Utilizes JAGS via `matjags` for model fitting, with parameters such as decision noise, spatial bias, information bonus, and learning rates.
- Outputs:
  - Fitted parameters for each participant.
  - Variance estimates for key model parameters.
  - Simulated participant choices for validation.
- Saves results to a timestamped CSV file (e.g., `HC_MUD_fits_<timestamp>.csv`).

### 5. `model_Kalman_Filter`
This file defines the Kalman Filter-based model for hierarchical Bayesian inference. Key features:
- Includes hyperpriors for learning rates, information bonuses, decision noise, and side bias.
- Models trial-by-trial inference, updating learning rates and expected rewards for each bandit.
- Computes choice probabilities based on differences in expected rewards, information bonuses, and biases.
- Simulates choices for validation against observed data.

#### Key Components:
- **Hyperpriors**: Define priors for starting and asymptotic learning rates, decision noise, information bonuses, and side biases.
- **Subject-Level Modeling**: Captures individual differences in learning and decision-making.
- **Inference Model**: Implements a Kalman filter for trial-by-trial updates of expected rewards and learning rates.
- **Choice Probabilities**: Computes probabilities of choosing each bandit using a softmax function.

## Dependencies
- MATLAB
- JAGS (Just Another Gibbs Sampler)
- `matjags` MATLAB interface for JAGS

## Contact
For any issues, contact Carter Goldman (cgoldman@laureateinstitute.org)
