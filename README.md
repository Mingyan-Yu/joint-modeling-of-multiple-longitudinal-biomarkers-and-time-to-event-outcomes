# joint-modeling-of-multiple-longitudinal-biomarkers-and-time-to-event-outcomes
Repository for code corresponding to the 'Joint Modeling of Multiple Longitudinal Biomarkers and Survival Outcomes via Threshold Regression: Variability as a Predictor' paper.

Instructions:
Simulation_data_generation.R: use this R file to simulate longitudinal biomarker data and time-to-event outcomes 

Simulation_joint_model_fit.R: use this R file to fit our proposed joint model with the data simulated from Simulation_data_generation.R
Simulation_joint_model_fit_standardize_time.R: use this R file to fit our proposed joint model with the data simulated from Simulation_data_generation.R, but with the measurement time of longitudinal biomarkers standardized to facilitate convergence

Simulation_TSIMThR_fit.R: use this R file to fit a two-stage altervative approach with the data simulated from Simulation_data_generation.R (this is for comparison with our proposed joint model)
Simulation_TSMEThR_fit.R: use this R file to fit another two-stage altervative approach with the data simulated from Simulation_data_generation.R (this is for comparison with our proposed joint model)
