# Delay Embedding and Nonlinear Mutual Dynamics (DENMD) Toolbox

<p align="center">
<img src="https://github.com/ParhamP/DENMD/blob/main/assets/DENMD-logos.jpeg?raw=true" width="305" height="305">
</p>


A toolbox for generating higher-dimensional manifolds of time-series data using delay embedding methods and measuring nonlinear mutual dynamics between them using cross mapping techniques.

## Overview

We've employed dynamical systems approaches to compile a toolbox that:
1. Performs dynamical embedding of a time series to assay complexity and determine meaningful dimensions of activity
2. Parses linear and nonlinear dynamics
3. Assays causal asymmetry in reciprocal interactions between correlated signals

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

The MATLAB file for this example can be found in "examples/lorenz folder".

```Matlab
% Generate two time-series data
[x1, y1, z1] = generate_lorenz(45, 10, 8/3, [0 1 1.05], [0 20], 0.000001);
[x2, y2, z2] = generate_lorenz(28, 10, 8/3, [0 1 1.05], [0 30], 0.000001);

%create a DENMD object
denmd = DENMD();

% Assign the time-series to the object
denmd.set_timeseries(1, y1, 'name', 'y1');
denmd.set_timeseries(2, y2, 'name', 'y2');

% Visualize the time-series from within the object
denmd.visualize_time_series();
```
<p align="center">
<img src="https://github.com/ParhamP/DENMD/blob/main/assets/vis_time_series_1.png?raw=true" width="477" height="398">
</p>
