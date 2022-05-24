# Delay Embedding and Nonlinear Mutual Dynamics (DENMD) Toolbox

<p align="center">
<img src="https://github.com/ParhamP/DENMD/blob/main/assets/DENMD-logos.jpeg?raw=true" width="305" height="305">
</p>


A toolbox for generating higher-dimensional attractors of time-series data using Time-Delay Embedding methods and measuring nonlinear mutual dynamics between them using cross mapping techniques.

## Overview

We've employed dynamical systems approaches to compile a toolbox that:
1. Performs Eigen Time-Delay (ETD) embedding of a time series to assay complexity and determine meaningful dimensions of activity (Broomhead, 1989).
2. Parses linear and nonlinear dynamics (Brunton, 2017).
3. Assays causal asymmetry in reciprocal interactions between correlated signals (Sugihara, 2012).

## Download

1. In the command line:
```
git clone https://github.com/ParhamP/DENMD.git
```

2. Or simply "Download Zip" from Github's "Code" tab

## Install

In MATLAB:

```
addpath(genpath('AbsolutePathToToolbox'))
```

## Example

The MATLAB file for this example can be found in "examples/lorenz/" folder.

```Matlab
% Generate two time-series data
[x1, y1, z1] = generate_lorenz(45, 10, 8/3, [0 1 1.05], [0 20], 0.000001);
[x2, y2, z2] = generate_lorenz(28, 10, 8/3, [0 1 1.05], [0 30], 0.000001);

%create a DENMD object
denmd = DENMD();

y1_attractor_num = 1;
y2_attractor_num = 2;

% Assign the time-series to the object
denmd.set_timeseries(y1_attractor_num, y1, 'name', 'y1');
denmd.set_timeseries(y2_attractor_num, y2, 'name', 'y2');

% Equalize the size of time-series in case they're different
denmd.equalize_timeseries_edges();

% Visualize the time-series from within the object
denmd.visualize_time_series();
```
<p align="center">
<img src="assets/vis_time_series_wo_title.png?raw=true" width="560" height="420">
</p>

```Matlab
% Set delay embedding parameters
delay_lag = 1;
embedding_dimension = 30;
use_SVD = 1;
SVD_num_components = 4;

% Generate attractors from time-series data using Eigen Time-Delay Embedding
denmd.generate_single_attractor(y1_attractor_num, delay_lag, embedding_dimension, use_SVD, SVD_num_components);
denmd.generate_single_attractor(y2_attractor_num, delay_lag, embedding_dimension, use_SVD, SVD_num_components);

% Visualize the attractors in 3D
denmd.visualize_attractors_3d();
```
<p align="center">
<img src="assets/vis_lorenz_attractors_wo_title.png?raw=true" width="627" height="623">
</p>

```Matlab
% Calculate Convergent Cross Mapping (CCM) scores to infer causality between the two systems
denmd.ccm('barPlot', true);
```
<p align="center">
<img src="assets/ccm_bar_plot.png?raw=true" width="560" height="420">
</p>

```Matlab
% Parse linear and nonlinear dynamics
std_from_mean = 2;
hit_range = 100;
nonlinear_threshold = denmd.estimate_nonlinear_threshold(y1_attractor_num, std_from_mean);
denmd.visualize_havok_attractor(y1_attractor_num, nonlinear_threshold, hit_range)
```
<p align="center">
<img src="assets/havok.png?raw=true" width="627" height="294">
</p>


## Collaborators

- Parham Pourdavood, BA
- Michael Jacob, MD, PhD

## Acknowledgements

This work was funded in part by [Human Energy](https://humanenergy.io/) and a grant from the department of Veterans Affairs.

## References
1. Broomhead D. S. and Jones Roger (1989). Time-series analysisProc. R. Soc. Lond. A423103–121.
2. Brunton, S. L., Brunton, B. W., Proctor, J. L., Kaiser, E., & Nathan Kutz, J. (2017). Chaos as an intermittently forced linear system. Nature Communications, 8(1).
3. Sugihara, G., May, R., Ye, H., Hsieh, C. H., Deyle, E., Fogarty, M., & Munch, S. (2012). Detecting causality in complex ecosystems. Science, 338(6106), 496–500.

