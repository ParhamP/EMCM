clear
addpath('../../')
% Load two EEG time-series data
alpha_eeg = load('data/alpha_eeg.mat').first_eeg;
gamma_eeg = load('data/gamma_eeg.mat').second_eeg;


%create an EMCM object
emcm = EMCM();

alpha_number = 1;
gamma_number = 2;

% Assign the time-series to the object
emcm.set_timeseries(alpha_number, alpha_eeg, 'name', 'Alpha');
emcm.set_timeseries(gamma_number, gamma_eeg, 'name', 'Gamma');

% Visualize the time-series from within the object
emcm.visualize_time_series();

% Set the delay embedding parameters
delay_lag = 1;
embedding_dimension = 60;

% First, generate attractors from time-series data using traditional time-delay-embedding
emcm.generate_single_attractor(alpha_number, delay_lag, embedding_dimension);
emcm.generate_single_attractor(gamma_number, delay_lag, embedding_dimension);

% Visualize the attractors in 3D
emcm.visualize_attractors_3d();

% Apply PCA to get Eigen-Time_Delay Attractors.
% The threshold method determines the number of components to keep.
% The threhold method of 'one' would retain all the components with eigenvalue greater than 1. 

emcm.apply_pca_and_set_r(alpha_number, 'threshold_method', 'one');
emcm.apply_pca_and_set_r(gamma_number, 'threshold_method', 'one');


% Visualize the attractors in 3D
emcm.visualize_attractors_3d();

% Access the complexity measures for eigen attractors (i.e., their number of retained components)
[alpha_complexity, gamma_complexity] = emcm.complexities('barPlot', true);

% Calculate Cross Mapping scores to infer shared dynamics between the two manifolds.
% The 'weigh_by_eigens' option will weigh each dimension's correlation by its eigenvalue.
[alpha_precits_gamma_score, gamma_predicts_alpha_score] = emcm.ccm('weigh_by_eigens', true, 'barPlot', true);
