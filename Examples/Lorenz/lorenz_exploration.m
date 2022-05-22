addpath('../../')
[x1, y1, z1] = generate_lorenz(45, 10, 8/3, [0 1 1.05], [0 20], 0.000001);
[x2, y2, z2] = generate_lorenz(28, 10, 8/3, [0 1 1.05], [0 30], 0.000001);
%size(y1)
%size(y2)
%y2 = y2(1:length(y1));
%size(y1)
%size(y2)
c = DENMD();
c.set_timeseries(1, y1);
c.set_timeseries(2, y2);
c.equalize_timeseries_edges()
c.generate_single_attractor(1, 1, 30, 1, 4);
c.generate_single_attractor(2, 1, 30, 1, 4);
c.zscore_attractors()
c.visualize_attractors_3d()
c.visualize_single_attractor_3d(1)

nonlinear_thresh_1 = c.estimate_nonlinear_threshold(1, 2);
nonlinear_thresh_2 = c.estimate_nonlinear_threshold(2, 2);
hit_range = 100;
c.visualize_havok_attractor(1, nonlinear_thresh_1, hit_range)

c.generate_linear_dynamics(1, nonlinear_thresh_1, hit_range);
c.generate_nonlinear_dynamics(1, nonlinear_thresh_1, hit_range);
c.generate_linear_dynamics(2, nonlinear_thresh_2, hit_range);
c.generate_nonlinear_dynamics(2, nonlinear_thresh_2, hit_range);
c.visualize_linear_to_nonlinear(2);

scores = c.ccm()
%{
[a, b, cc, d, e, f, g, h] = c.all_linear_and_nonlinear_ccm_scores();
%scores = c.ccm();
sc1 = scores(1);
sc2 = scores(2);
figure; bar([sc1 sc2 a b cc d e f g h]); xticklabels({'X1->X2', 'X2->X1', 'SL', 'SNL', 'OL', 'ONL', 'SL', 'SNL', 'OL', 'ONL'});
sc3 = c.linear_correlation();
%c.visualize_predicted_attractors()
%c.visualize_time_series()
c.visualize_nonlinear_forcing(1)
c.visualize_nonlinear_forcing(2)

x = categorical({'Linear', 'CCM (X1->X2)', 'CCM (X2->X1)'});
x = reordercats(x, {'Linear', 'CCM (X1->X2)', 'CCM (X2->X1)'});
y = [sc3, sc1, sc2];

figure; bar(x, y); title('CCM vs Pearson Correlation'); ylabel('Scores');


%}

