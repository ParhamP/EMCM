# Eigen Manifold Cross Mapping (EMCM) Toolbox

<!-- <p align="center">
<img src="https://github.com/ParhamP/DENMD/blob/main/assets/DENMD-logos.jpeg?raw=true" width="305" height="305">
</p> -->


EMCM is a MATLAB toolbox for generating high-dimensional manifolds from time-series data using Eigen-Time-Delay Embedding and measuring their shared dynamics using Cross Mapping techniques. 

## Overview

1. Perform Eigen Time-Delay (ETD) embedding of a time series to generate high-dimensional manifolds (Brunton, 2017)
2. Use the eigenvalue statistics to assay complexity and capture meaningful dimensions of activity.
3. Infer shared dynamics between ETD Manifolds Using Sugihara's Convergent Cross Mapping (Sugihara, 2012).

## Download

1. In the command line:
```
git clone https://github.com/ParhamP/EMCM.git
```

2. Or simply "Download Zip" from Github's "Code" tab

## Install

In MATLAB:

```
addpath(genpath('AbsolutePathToToolbox'))
```

## Example

The MATLAB file for this example can be found in "examples/eeg_example/" folder.

```Matlab
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
```
<p align="center">
<img src="assets/timeseries.png?raw=true" width="560" height="420">
</p>

```Matlab
% Set the delay embedding parameters
delay_lag = 1;
embedding_dimension = 60;

% First, generate attractors from time-series data using traditional time-delay-embedding
emcm.generate_single_attractor(alpha_number, delay_lag, embedding_dimension);
emcm.generate_single_attractor(gamma_number, delay_lag, embedding_dimension);

% Visualize the attractors in 3D
emcm.visualize_attractors_3d();
```
<p align="center">
<img src="assets/pre_etd_attractors.png?raw=true" width="427" height="627">
</p>

```Matlab
% Apply PCA to get Eigen-Time_Delay Attractors.
% The threshold method determines the number of components to keep.
% The threhold method of 'one' would retain all the components with eigenvalue greater than 1. 

c1.apply_pca_and_set_r(alpha_number, 'threshold_method', 'one');
c1.apply_pca_and_set_r(gamma_number, 'threshold_method', 'one');


% Visualize the attractors in 3D
emcm.visualize_attractors_3d();
```
<p align="center">
<img src="assets/attractors.png?raw=true" width="427" height="627">
</p>

```Matlab
% Calculate Convergent Cross Mapping (CCM) scores to infer shared dynamics between the two manifolds.
% The 'weigh_by_eigens' option will weigh each dimension's correlation by its eigenvalue.
emcm.ccm('weigh_by_eigens', weigh_dimensions, 'barPlot', true);
```
<p align="center">
<img src="assets/barplot.png?raw=true" width="560" height="420">
</p>


## Collaborators

- Parham Pourdavood
- Michael Jacob, MD, PhD


## References
1. Brunton, S. L., Brunton, B. W., Proctor, J. L., Kaiser, E., & Nathan Kutz, J. (2017). Chaos as an intermittently forced linear system. Nature Communications, 8(1).
2. Sugihara, G., May, R., Ye, H., Hsieh, C. H., Deyle, E., Fogarty, M., & Munch, S. (2012). Detecting causality in complex ecosystems. Science, 338(6106), 496–500.

