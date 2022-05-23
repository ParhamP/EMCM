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
denmd.equalize_timeseries_edges()

% Visualize the time-series from within the object
denmd.visualize_time_series();

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

% Calculate Convergent Cross Mapping (CCM) scores to infer causality between the two systems
denmd.ccm('barPlot', true);

% Parse linear and nonlinear dynamics
std_from_mean = 2;
hit_range = 100;
nonlinear_threshold = denmd.estimate_nonlinear_threshold(y1_attractor_num, std_from_mean);
denmd.visualize_havok_attractor(y1_attractor_num, nonlinear_threshold, hit_range)
