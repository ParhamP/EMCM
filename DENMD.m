classdef DENMD < handle
    % DENMD Class for visualizing and comparing attractors.
    % This class is able to get two time-series, embed them in higher
    % dimensional attractors, and apply cross-mapping to infer causality
    % between them. There are functions for visualization, as well as tools
    % for estimating relevant parameters. 
    properties %(Access=private)
        x1;
        x2;
        Xm;
        Ym;
        Xm_linear;
        Ym_linear;
        Xm_nonlinear;
        Ym_nonlinear;
        Xm_linear_points;
        Ym_linear_points;
        Xm_nonlinear_points;
        Ym_nonlinear_points;
        E;
        E1;
        E2;
        r1;
        r2;
        ts_length;
        x1_pred;
        x2_pred;
    end
    
    methods
      function obj = DENMD()
          % Constructs an DENMD object.
          %
          % OBJ = AttractorAnalysis()
          return
      end
      
      function set_timeseries(obj, num, x)
          % SET_TIMESERIES sets the NUMth time-series X.
          %
          % SET_TIMESERIES(NUM, X)
          if num == 1
              obj.x1 = x;
          elseif num == 2
              obj.x2 = x;
          else
              error("invalid timeseries number.");
          end
      end
      
      function timeseries=get_timeseries(obj, num)
          % GET_TIMESERIES Get the NUMth timeseries.
          %
          %TIMESERIES = GET_TIMESERIES(NUM)
          if num == 1
              timeseries = obj.x1;
          elseif num == 2
              timeseries = obj.x2;
          end
      end
      
      function set_attractor(obj, num, attractor)
          % SET_ATTRACTOR Sets the NUMth attractor to ATTRACTOR.
          %
          % SET_ATTRACTOR(NUM, ATTRACTOR)
          if num == 1
              obj.Xm = attractor;
          elseif num == 2
              obj.Ym = attractor;
          end
      end
      
      function attractor=get_attractor(obj, num)
          % GET_ATTRACTOR Return the NUMth attractor.
          %
          % ATTRACTOR = GET_ATTRACTOR(NUM)
          if num == 1
              attractor = obj.Xm;
          elseif num == 2
              attractor = obj.Ym;
          end
      end
      
      function res = num_linear_points(obj, num)
          if num == 1
              linear_points = obj.Xm_linear_points;
          elseif num == 2
              linear_points = obj.Ym_linear_points;
          end
          res = length(linear_points);
      end
      
      function res = num_nonlinear_points(obj, num)
          if num == 1
              nonlinear_points = obj.Xm_nonlinear_points;
          elseif num == 2
              nonlinear_points = obj.Ym_nonlinear_points;
          end
          res = length(nonlinear_points);
      end
      
      function res = num_total_points(obj, num)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          res = length(arr);
      end
      
      function equalize_edges(obj)
          % EQUALIZE_EDGES Trimes the larger attractor to have the same
          % length as the shorter attractor, in case they're not the same
          % length.
          %
          % EQUALIZE_EDGES()
          X1 = obj.get_attractor(1);
          X2 = obj.get_attractor(2);
          if length(X1) > length(X2)
              X1 = X1(1:length(X2), :);
          else
              X2 = X2(1:length(X1), :);
          end
          obj.set_attractor(1, X1);
          obj.set_attractor(2, X2);
      end
      
      function equalize_timeseries_edges(obj)
          X1 = obj.get_timeseries(1);
          X2 = obj.get_timeseries(2);
          if length(X1) > length(X2)
              X1 = X1(1:length(X2));
          else
              X2 = X2(1:length(X1));
          end
          obj.set_timeseries(1, X1);
          obj.set_timeseries(2, X2);
      end
      
      function sizes(obj, include_linear_and_nonlinear)
          % SIZES Shows the sizes for the attractors. If
          % INCLUDE_LINEAR_AND_NONLINEAR is 1, the sizes of the linear and
          % nonlinear attractors are shown, 0 otherwise. 
          %
          %SIZES()
          %
          %SIZES(INCLUDE_LINEAR_AND_NONLINEAR)
          disp('');
          fprintf('Xm: (%d, %d)\n', size(obj.Xm));
          if nargin > 1 && include_linear_and_nonlinear == 1
              fprintf('Linear Xm: (%d, %d)\n', size(obj.Xm_linear));
              fprintf('NonLinear Xm: (%d, %d)\n', size(obj.Xm_nonlinear));
          end
          
          disp('----------------------');
          
          fprintf('Ym: (%d, %d)\n', size(obj.Ym));
          if nargin > 1 && include_linear_and_nonlinear == 1
              fprintf('Linear Ym: (%d, %d)\n', size(obj.Ym_linear));
              fprintf('NonLinear Ym: (%d, %d)\n', size(obj.Ym_nonlinear));
          end
      end
      
      function visualize_time_series(obj)
          % VISUALIZE_TIME_SERIES Plots the object's time-series. 
          %
          % VISUALIZE_TIME_SERIES()
          figure
          sgtitle('Time Series');
          subplot(2,1,1)
          plot(obj.x1);
          title('x1');
          
          subplot(2,1,2)
          plot(obj.x2);
          title('x2');
      end
      
      function visualize_delayed_timeseries(obj, num, tau, E)
          % VISUALIZE_DELAYED_TIMESERIES Creates a plot with NUMth
          % timeseries and its delayed version, had it been embedded with
          % delay of TAU and dimension of E.
          %
          % VISUALIZE_DELAYED_TIMESERIES(NUM, TAU, E)
          if num == 1
              arr = obj.x1;
              arr = obj.embeder(arr, tau, E);
          elseif num == 2
              arr = obj.x2;
              arr = obj.embeder(arr, tau, E);
          end
          my_point = 1900;
          figure;
          subplot(3,1,1)
          plot(arr(:, 1)); hold on;
          plot(my_point, arr(my_point, 1), 'o', 'Color', 'r');
          title('X(t)', 'FontSize', 20);
          
          subplot(3,1,2)
          plot(arr(:, 2)); hold on;
          plot(my_point - tau, arr(my_point - tau, 2), 'o', 'Color', 'r');
          title('X(t - \tau)', 'FontSize', 20);
          
          subplot(3,1,3)
          plot(arr(:, 3)); hold on;
          plot(my_point - 2 * tau, arr(my_point - 2 * tau, 3), 'o', 'Color', 'r');
          title('X(t - 2\tau)', 'FontSize', 20);
          
          figure;
          
          plot3(arr(:, 1), arr(:, 2), arr(:, 3));
          xlabel('X(t)', 'FontSize', 20)
          ylabel('X(t - \tau)', 'FontSize', 20);
          zlabel('X(t - 2\tau)', 'FontSize', 20);
      end
      
      function Xm = embeder(~, x, tau, E)
          % EMBEDER Embedes X in E dimension with lag of TAU.
          %
          % EMBEDER(X, TAU, E)
          tsize = length(x);
          t_iter = tsize-(tau*(E-1)) - 1;
          Xm = zeros(t_iter, E);
          for ii=1:t_iter
              end_val = ii+tau*(E-1)+1;
              Xm(ii, :) = x(ii:tau:end_val - 1);
          end
      end
      
      function upsample_linear_attractor(obj, num, scale)
          % UPSAMPLE_LINEAR_ATTRACTOR Upsamples NUMth linear attractor with
          %  scale of SCALE.
          if num == 1
              arr = obj.Xm_linear;
          elseif num == 2
              arr = obj.Ym_linear;
          end
          arr_size = size(arr);
          cols = arr_size(2);
          rows = arr_size(1);
          arr_interp = zeros(rows * scale, cols);
          for c=1:cols
              arr_interp(:, c) = interp(arr(:, c), scale);
          end
                   
          if num == 1
              obj.Xm_linear = arr_interp;
          elseif num ==2
              obj.Ym_linear = arr_interp;
          end 
      end
      
      function upsample_nonlinear_attractor(obj, num, scale)
          % UPSAMPLE_NONLINEAR_ATTRACTOR Upsamples NUMth nonlinear
          %  attractor with scale of SCALE.
          if num == 1
              arr = obj.Xm_nonlinear;
          elseif num == 2
              arr = obj.Ym_nonlinear;
          end
          arr_size = size(arr);
          cols = arr_size(2);
          rows = arr_size(1);
          arr_interp = zeros(rows * scale, cols);
          for c=1:cols
              arr_interp(:, c) = interp(arr(:, c), scale);
          end
                   
          if num == 1
              obj.Xm_nonlinear = arr_interp;
          elseif num ==2
              obj.Ym_nonlinear = arr_interp;
          end 
      end
      
      function upsample_attractor(obj, num, scale)
          % UPSAMPLE_ATTRACTOR Upsamples NUMth attractor with scale of
          % SCALE.
          if num == 1
                  arr = obj.Xm;
          elseif num == 2
                arr = obj.Ym;
          end
          
          arr_size = size(arr);
          cols = arr_size(2);
          rows = arr_size(1);
          arr_interp = zeros(rows * scale, cols);
          for c=1:cols
              arr_interp(:, c) = interp(arr(:, c), scale);
          end
                   
          if num == 1
              obj.Xm = arr_interp;
          elseif num == 2
                obj.Ym = arr_interp;
          end 
      end
      
      function downsample_attractor(obj, num, linear_or_not, scale)
          % DOWNSAMPLE_ATTRACTOR Downsamples NUMth attractor with scale
          % of SCALE. LINEAR_OR_NOT is an array of size 2 indicating the
          % linearity and nonlinearity of the first and second attractor
          % using 1 for linear, 2 for nonlinear, and 0 for standard.
          %
          % DOWNSAMPLE_ATTRACTOR(NUM, LINEAR_OR_NOT, SCALE)
          if num == 1
              if linear_or_not(1) == 1
                  arr1 = obj.Xm_linear;
              elseif linear_or_not(1) == 2
                  arr1 = obj.Xm_nonlinear;
              else
                  arr1 = obj.Xm;
              end
              if linear_or_not(2) == 1
                  arr2 = obj.Ym_linear;
              elseif linear_or_not(2) == 2
                  arr2 = obj.Ym_nonlinear;
              else
                  arr2 = obj.Ym;
              end
          elseif num == 2
              if linear_or_not(2) == 1
                  arr1 = obj.Ym_linear;
              elseif linear_or_not(2) == 2
                  arr1 = obj.Ym_nonlinear;
              else
                  arr1 = obj.Ym;
              end
              if linear_or_not(1) == 1
                  arr2 = obj.Xm_linear;
              elseif linear_or_not(1) == 2
                  arr2 = obj.Xm_nonlinear;
              else
                  arr2 = obj.Xm;
              end
          end
          
          if nargin > 3
               arr1_downsampled = downsample(arr1, scale);
          else
              ratio = length(arr1) / length(arr2);
              ratio_floored = floor(ratio);

              arr1_downsampled = downsample(arr1, ratio_floored);
              arr1_downsampled = arr1_downsampled(1:length(arr2), :);
          end
          
          if num == 1
              if linear_or_not(1) == 1
                  obj.Xm_linear = arr1_downsampled;
              elseif linear_or_not(1) == 2
                  obj.Xm_nonlinear = arr1_downsampled;
              else
                  obj.Xm = arr1_downsampled;
              end
          elseif num == 2
              if linear_or_not(2) == 1
                  obj.Ym_linear = arr1_downsampled;
              elseif linear_or_not(2) == 2
                  obj.Ym_nonlinear = arr1_downsampled;
              else
                  obj.Ym = arr1_downsampled;
              end
          end 
      end
      
      function est_dimension = estimate_dimension_using_fnn(obj, num, tau, ...
                                                            max_dim)
          % ESTIMATE_DIMENSION_USING_FNN Estimates the optimal dimension
          % for the NUMth timeseries using false nearest neighbor
          % algorithm. It embeds it with TAU and uses dimensions up to
          % MAX_DIM.
          %
          %ESTIMATE_DIMENSION_USING_FNN(NUM, TAU, MAX_DIM)
          
          rtol = 15.0;
          atol = 2.0;
          min_dim = 1;
          if nargin < 3
            max_dim = 50;
          end
          dim_range = min_dim:max_dim;
          if num == 1
              orig_attractor = obj.get_attractor(num);
              temp_r = obj.r1;
          elseif num == 2
              orig_attractor = obj.get_attractor(num);
              temp_r = obj.r2;
          end
          
          x_ts = obj.get_timeseries(num);
          
          fnn_scores = zeros(1, max_dim);
          
         for d=dim_range
            obj.generate_single_attractor(num, tau, d, 0, nan);
            x_d = obj.get_attractor(num);
            obj.generate_single_attractor(num, tau, d + 1, 0, nan);
            x_d_1 = obj.get_attractor(num);
            [~, ~, f3] = single_dimension_fnn(x_ts, x_d, x_d_1, rtol, atol);
            fnn_scores(d) = f3;
         end
        [~, i] = min(fnn_scores); 
        est_dimension = dim_range(i);
        obj.set_attractor(num, orig_attractor);
        if num == 1
            obj.r1 = temp_r;
        elseif num == 2
            obj.r2 = temp_r;
        end
      end
      
      function estimate_period(obj, num)
          % ESTIMATE_PERIOD Estimates the period of NUMth attractor by
          % plotting its fourier transform peaks. 
          %
          %ESTIMATE_PERIOD(NUM)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end

        y1 = arr(:, 1);

        t = 1:1:length(y1); % Time Vector
        s = y1; % signal Vector
        L = numel(t); % Signal Length
        Ts = mean(diff(t)); % Sampling Interval
        Fs = 1/Ts; % Sampling Frequency
        Fn = Fs/2; % Nyquist Frequency
        sm = s - mean(s); 
        FTs = fft(sm)/L; % Normalised Fourier Transform
        Fv = linspace(0, 1, fix(L/2)+1)*Fn; % Frequency Vector
        Iv = 1:numel(Fv); % Index Vector

        [~,locs] = findpeaks(abs(FTs(Iv)));
        Perd = 1./Fv(locs);
        figure; plot(Perd);
      end
      
      
      function r_est = estimate_r_obs(obj, num, tau, E)
          if num == 1
              x = obj.x1;
              X = obj.embeder(x, tau, E);
                [~,S,~] = svd(X,'econ');
                sigs = diag(S);
                beta = size(X,2)/size(X,1);
                thresh = optimal_SVHT_coef(beta,0) * median(sigs);
                r_est = length(sigs(sigs>thresh));
          elseif num == 2
              y = obj.x2;
              Y = obj.embeder(y, tau, E);
                [~,S,~] = svd(Y,'econ');
                sigs = diag(S);
                beta = size(Y,2)/size(Y,1);
                thresh = optimal_SVHT_coef(beta,0) * median(sigs);
                r_est = length(sigs(sigs>thresh));
          end
      end
      
      function r_est = estimate_r(obj, num, tau, E, thresh)
          if num == 1
              x = obj.x1;
              X = obj.embeder(x, tau, E);
              [~,~,latent,~,~] = pca(X);
          elseif num == 2
              y = obj.x2;
              Y = obj.embeder(y, tau, E);
              [~,~,latent,~,~] = pca(Y);
          end
          total_var_sum = sum(latent);
          for i=1:length(latent)
              sofar_var_sum = sum(latent(1:i));
              var_ratio = sofar_var_sum / total_var_sum;
              if var_ratio >=thresh
                  r_est = i;
                  return
              end
          end
          r_est = length(latent);
      end
      
      function generate_single_attractor(obj, num, tau, E, use_svd, r)
          %addpath('./HAVOK/utils');
          % GENERATE_SINGLE_ATTRACTOR Generates the NUMth attractor with
          % tau of TAU and dimension of E. Binary USE_SVD indicates whether
          % or not to use SVD on the constructed attractor. If svd is
          % applied, the method picks R components of U. 
          if num == 1
              obj.E1 = E;
              obj.r1 = r;
              x = obj.x1;
              X = obj.embeder(x, tau, E);
              if use_svd==1
                [Ux,~,~] = svd(X,'econ');
                X = Ux(1:end,1:r);
                obj.E1 = r;
              end
              obj.Xm = X;
          elseif num == 2
              obj.E2 = E;
              obj.r2 = r;
              y = obj.x2;
              Y = obj.embeder(y, tau, E);
              if use_svd==1
                [Uy,~,~] = svd(Y,'econ');
                Y = Uy(1:end,1:r);
                obj.E2 = r;
              end
              obj.Ym = Y;
          end
      end
      
      function thresh = estimate_nonlinear_threshold(obj, num, ...
                                                     std_from_mean)
          % ESTIMATE_NONLINEAR_THRESHOLD Estimates the threshold for the
          % bursty activity of the nonlinear forcing of NUMth attractor,
          % using standard deviation from the mean STD_FROM_MEAN.
          %
          % THRESH = ESTIMATE_NONLINEAR_THRESHOLD(NUM< STD_FROM_MEAN)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          n = arr(:, end);
          thresh = mean(n)+ std_from_mean * std(n);
      end
      
      function visualize_linear_to_nonlinear(obj, num)
          if num == 1
              arr = obj.Xm;
              self_linear_points = obj.Xm_linear_points;
              self_nonlinear_points = obj.Xm_nonlinear_points;
              other_linear_points = obj.Ym_linear_points;
              other_nonlinear_points = obj.Ym_nonlinear_points;
          elseif num == 2
              arr = obj.Ym;
              self_linear_points = obj.Ym_linear_points;
              self_nonlinear_points = obj.Ym_nonlinear_points;
              other_linear_points = obj.Xm_linear_points;
              other_nonlinear_points = obj.Xm_nonlinear_points;
          end
          [linear_linear_intersection, ~] = intersect(self_linear_points, other_linear_points);
          [nonlinear_nonlinear_intersection, ~] = intersect(self_nonlinear_points, other_nonlinear_points);
          [linear_nonlinear_intersection, ~] = intersect(other_linear_points, self_nonlinear_points);
          [nonlinear_linear_intersection, ~] = intersect(other_nonlinear_points, self_linear_points);
          f = figure;
          f.Position = [300, 50, f.Position(3) * 2.2, f.Position(4) * 1.4];
          plot3(arr(:, 1), arr(:, 2), arr(:, 3), 'k'); hold on;
          scatter3(arr(self_nonlinear_points, 1), arr(self_nonlinear_points, 2), arr(self_nonlinear_points, 3), 'r', 'filled');
          legend('self Linear','self NonLinear');
          f = figure;
          f.Position = [300, 50, f.Position(3) * 2.2, f.Position(4) * 1.4];
          plot3(arr(:, 1), arr(:, 2), arr(:, 3), 'k'); hold on;
          scatter3(arr(other_linear_points, 1), arr(other_linear_points, 2), arr(other_linear_points, 3), 'g', 'filled');
          hold on;
          scatter3(arr(other_nonlinear_points, 1), arr(other_nonlinear_points, 2), arr(other_nonlinear_points, 3), 'y', 'filled');
          legend('self linear', 'other linear', 'other NonLinear');
          
      end
    
      function generate_linear_dynamics(obj, num, thresh, hit_range)
          % GENERATE_LINEAR_DYNAMICS Identifies and stores the linear
          % regions of the NUMth attractor using bursty threshold of
          % THRESH and burst hit range of HIT_RANGE.
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          
          [startvals, endvals] = obj.find_startvals_endvals(num, thresh,...
                                                            hit_range);
            
            linear_points = [];                                            
            if ~isempty(startvals)
                linear_points = [linear_points 1:startvals(1) - 1];
            else
                linear_points = [linear_points 1:length(arr)];
            end
            %count_points = length(linear_points);
            for k=1:length(startvals)
                %new_points = arr(endvals(k) + 1:startvals(k+1),1);
                %count_points = count_points + length(new_points);
                if k ~= length(startvals)
                    linear_points = [linear_points endvals(k) + 1:startvals(k+1) - 1];
                else
                    linear_points = [linear_points endvals(k) + 1:length(arr)];
                end 
            end
            
            all_points_arr = arr(linear_points, :);
            
            if num == 1
              obj.Xm_linear = all_points_arr;
              obj.Xm_linear_points = linear_points;
            elseif num == 2
              obj.Ym_linear = all_points_arr;
              obj.Ym_linear_points = linear_points;
            end
      end
      
      function generate_nonlinear_dynamics(obj, num, thresh, hit_range)
          % GENERATE_NONLINEAR_DYNAMICS Identifies and stores the nonlinear
          % regions of the NUMth attractor using bursty threshold of
          % THRESH and burst hit range of HIT_RANGE.
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          
          [startvals, endvals] = obj.find_startvals_endvals(num, thresh,...
                                                            hit_range);
                                                        
          if isempty(startvals) || startvals(1) == endvals(1)
               if num == 1
                obj.Xm_nonlinear = [];
                obj.Xm_nonlinear_points = [];
              elseif num == 2
                obj.Ym_nonlinear = [];
                obj.Ym_nonlinear_points = [];
              end
          end
            
          nonlinear_points = [];
            
          %count_points = 0;
          for k=1:length(startvals)
              %new_points = arr(startvals(k):endvals(k),1);
              %count_points = count_points + length(new_points);
              nonlinear_points = [nonlinear_points startvals(k):endvals(k)];
          end
            
          all_points_arr = arr(nonlinear_points, :);
            
          if num == 1
            obj.Xm_nonlinear = all_points_arr;
            obj.Xm_nonlinear_points = nonlinear_points;
          elseif num == 2
            obj.Ym_nonlinear = all_points_arr;
            obj.Ym_nonlinear_points = nonlinear_points;
          end
      end
      
      function [startvals, endvals] = find_startvals_endvals(obj, num, ...
                                                             thresh, ...
                                                             hit_range)
          % FIND_STARTVALS_ENDVALS Finds start values and end values for
          % linear and nonlinear regions of the NUMth attractor with thresh
          % of THRESH and hit range of HIT_RANGE
          %
          % [STARTVALS, ENDVALS] = FIND_STARTVALS_ENDVALS(NUM, THRESH,
          % HIT_RANGE)
          if num == 1
              arr = obj.Xm;
              r = obj.r1;
          elseif num == 2
              arr = obj.Ym;
              r = obj.r2;
          end
          
            arrL = 1:length(arr);
            inds = abs(arr(arrL, r)) > thresh;
            startvals = [];
            endvals = [];
            finders=(find(inds));
            if isempty(finders)
                %startvals = [startvals; 1; length(arr)];
                %endvals = [endvals; 1 ; length(arr)];
                return
            end
            start = finders(1);
            clear interval hits endval newhit
            %startvals = [startvals; start];
            while start <= length(arrL) - hit_range
                endmax = start+ hit_range;
                interval = start + 1:endmax;
                hits = find(inds(interval));
                if isempty(hits)
                    endval = start;
                    newhit = find(inds(endval+1:end));
                    if isempty(newhit)
                        %startvals = [startvals; length(arr)];
                        %endvals = [endvals; length(arr)];
                        break
                    end
                    %endvals = [endvals; endval];
                    start = endval+newhit(1); %check?
                    continue
                end
                startvals = [startvals; start];
                
                endval = start+hits(end);
                endvals = [endvals; endval];
                newhit = find(inds(endval+1:end));
                if isempty(newhit)
                    %startvals = [startvals; length(arr)];
                    %endvals = [endvals; length(arr)];
                    break
                end
                start = endval+newhit(1); %check?
            end
      end
      
      function visualize_havok_attractor_experimental(obj, num, my_title, my_dest)
          f = figure('visible','off');
          if num == 1
              arr = obj.Xm;
              nonlinear_points = obj.Xm_nonlinear_points;
          elseif num == 2
              arr = obj.Ym;
              nonlinear_points = obj.Ym_nonlinear_points;
          else
              error('wrong num');
          end
          
          scatter3(arr(:, 1), arr(:, 2), arr(:, 3), 'k', 'filled'); hold on;
          %plot3(arr(:, 1), arr(:, 2), arr(:, 3), 'k'); hold on;
          scatter3(arr(nonlinear_points, 1), arr(nonlinear_points, 2), arr(nonlinear_points, 3), 'r', 'filled');
          axis tight
          axis off
          view(270,0)
          
          title(my_title);
          saveas(f, my_dest)
      end
      
      function visualize_havok_attractor(obj, num, thresh, hit_range, ...
                                          custom_f, my_title, my_dest)
          % VISUALIZE_HAVOK_ATTRACTOR Plots the NUMth attractor with its
          % linear and nonlinear regions highlighted in black and red
          % respectively. The HAVOK's nonlinear threshold is THRESH and hit
          % range is set to HIT_RANGE. If the arguments MY_TITLE and
          % MY_DEST are provided, the figure will be saved in MY_DEST with
          % title MY_TITLE.
          %
          %VISUALIZE_HAVOK_ATTRACTOR(NUM, THRESH, HIT_RANGE)
          %
          % VISUALIZE_HAVOK_ATTRACTOR(NUM, THRESH, HIT_RANGE, MY_TITLE,
          % MY_DEST)
          if nargin < 5
              my_title = '';
              my_dest = '';
          end
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          [startvals, endvals] = obj.find_startvals_endvals(num, thresh, hit_range);
          if nargin == 4
              f = figure;
              f.Position = [f.Position(1), f.Position(2), ...
              f.Position(3) * 1.6, f.Position(4) * 1];
%               f.Position = [300, 50, 560 * 1.6, f.Position(4)];

          elseif nargin == 5
              f = custom_f;
          elseif nargin > 5
              f = figure('visible','off');  
              f.Position = [f.Position(1), f.Position(2), ...
              f.Position(3) * 2, f.Position(4) * 2];
          end
          %{
            if nargin > 4
                f = figure('visible','off');  
            else
                f = figure;
            end
          %}

            
            
            numhits = length(startvals);
            
            for k=1:numhits
                plot3(arr(startvals(k):endvals(k),1),arr(startvals(k):endvals(k),2),arr(startvals(k):endvals(k),3),'r','LineWidth',1.5), hold on
            end
            
            for k=1:numhits-1
                plot3(arr(endvals(k):startvals(k+1),1),arr(endvals(k):startvals(k+1),2),arr(endvals(k):startvals(k+1),3),'Color',[.25 .25 .25],'LineWidth',1.5), hold on
            end
           axis tight
           axis off
%            view(270,0)
           %title('Complete Attractor')
           
           if nargin > 5
              title(my_title);
              saveas(f, my_dest)
           end
      end
      
      function visualize_single_attractor_3d(obj, num, my_title, dest)
          % VISUALIZE_SINGLE_ATTRACTOR Plots the NUMth attractor in 3D. If
          % MY_TITLE and DEST are given, the figure will be saved in DEST
          % with title MY_TITLE.
          %
          % VISUALIZE_SINGLE_ATTRACTOR(NUM)
          %
          % VISUALIZE_SINGLE_ATTRACTOR(NUM, MY_TITLE, DEST)
           if nargin > 2
                f = figure('visible','off');
                f.Position = [300, 50, 560 * 1.6, f.Position(4)];
           else
                f = figure;
                f.Position = [300, 50, 560 * 1.6, f.Position(4)];
           end
          
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          else
              error('wrong num');
          end
          
          plot3(arr(:, 1), arr(:, 2), arr(:, 3), 'k');
          axis tight
          axis off
%           view(270,0)
          if nargin > 2
              title(my_title);
              saveas(f, dest)
          end
      end
      
      function save_timeseries_plot(obj, num, my_title, dest)
          % SAVE_TIMESERIES_PLOT Saves the NUMth time-series in DEST with
          % title MY_TITLE.
          %
          % SAVE_TIMESERIES_PLOT(NUM, MY_TITLE, DEST)
          f = figure('visible','off');    
          
          if num == 1
              arr = obj.x1;
          elseif num == 2
              arr = obj.x2;
          else
              error('wrong num');
          end
          
          plot(arr, 'k');

          title(my_title);
          saveas(f, dest)
          
      end
      
      function visualize_attractors_3d(obj, l_arr)
          % VISUALIZE_ATTRACTORS_3D Plots the two attractors. if L_ARR is
          % provided, it plots either the linear or the nonlinear
          % attractors. L_ARR is an array of size two that indicates
          % whether each attractor is regular, linear, or nonlinear with
          % numbers 0, 1, 2, respectively.
          %
          % VISUALIZE_ATTRACTORS_3D()
          %
          % VISUALIZE_ATTRACTORS_3D(L_ARR)
          f = figure;
          f.Position = [300, 50, 560 * 1.6, 560 * 1.6];
          sgtitle('Attractors');
          subplot(2,1,1)
          X = obj.Xm;
          Y = obj.Ym;
          
          if nargin > 1
              if l_arr(1) == 1
                  X = obj.Xm_linear;
              elseif l_arr(1) == 2
                  X = obj.Xm_nonlinear;
              end
          end
          
          plot3(X(:, 1), X(:, 2), X(:, 3), 'k');
          hold on;
          
          if nargin > 1
              if l_arr(2) == 1
                  Y = obj.Ym_linear;
              elseif l_arr(2) == 2
                  Y = obj.Ym_nonlinear;
              end
          end

          subplot(2,1,2)
          plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'k');
      end
      
      function similar_finds = find_similar_dynamics_points(obj, num, ...
                                                            thresh)
          % FIND_SIMILAR_DYNAMICS_POINTS Finds the points in the prediction
          % made by the NUMth attractor and the other attractor that are
          % close to each other with a threshold THRESH.
          %
          % SIMILAR_FINDS = FIND_SIMILAR_DYNAMICS_POINTS(NUM, THRESH)
          if num == 1
              arr2 = obj.Ym;
              preds = obj.x2_pred;
          elseif num == 2
              arr2 = obj.Xm;
              preds = obj.x1_pred;
          end
          dists = diag(pdist2(preds, arr2));
          similar_finds = find(dists > thresh);
      end
      
      function visualize_predicted_attractors(obj)
          % VISUALIZE_PREDICTED_ATTRACTORS Plots the prediction points the
          % attractors make of each other on top of the original attractors
          % they would be compared to in CCM.
          %
          % VISUALIZE_PREDICTED_ATTRACTORS()
          f = figure;
          f.Position = [100, 100, 560 * 2, 560 * 2];
          sgtitle('Attractors');
          subplot(2,1,1)
          scatter3(obj.x1_pred(:, 1), obj.x1_pred(:, 2), obj.x1_pred(:, 3), 'k');
          hold on;
          plot3(obj.Xm(:, 1), obj.Xm(:, 2), obj.Xm(:, 3), 'r', 'LineWidth', 2);
          subplot(2,1,2)
          scatter3(obj.x2_pred(:, 1), obj.x2_pred(:, 2), obj.x2_pred(:, 3), 'k');
          hold on;
          plot3(obj.Ym(:, 1), obj.Ym(:, 2), obj.Ym(:, 3), 'r', 'LineWidth', 2);
      end
      
      function visualize_nonlinear_forcing(obj, num)
          % VISUALIZE_NONLINEAR_FORCING Plots the nonlinear forcing
          % component of the NUMth SVD'd attractor.
          %
          % VISUALIZE_NONLINEAR_FORCING(num)
          if num == 1
              arr = obj.Xm;
              r = obj.r1;
          elseif num == 2
              arr = obj.Ym;
              r = obj.r2;
          end
          figure; plot(arr(:, r));
      end
      
      function zscore_attractors(obj)
          obj.Xm = reshape(zscore(obj.Xm(:)), size(obj.Xm));
          obj.Ym = reshape(zscore(obj.Ym(:)), size(obj.Ym));
      end
      
      function estimate_tau(obj, num, E, max_tau)
          % ESTIMATE_TAU Picks the optimal TAU from taus in the range of
          % one to MAX_TAU that would result in an embedding in E
          % dimensinos (uses mutual information).
          %
          % ESTIMATE_TAU(NUM, E, MAX_TAU)
          if nargin < 4
              range = 1:1:10;
          else
              range = 1:1:max_tau;
          end
          
          
          
          if num == 1
              orig_attractor = obj.get_attractor(num);
              orig_r = obj.r1;
          elseif num == 2
              orig_attractor = obj.get_attractor(num);
              orig_r = obj.r2;
          end

          mi_scores = [];
          for i=range
              obj.generate_single_attractor(num, i, E, 0, nan);
              arr = obj.get_attractor(num);
              mi_score = mi(arr(:, 1), arr(:, 2));
              mi_scores = [mi_scores mi_score];
          end
          
          figure;plot(range, mi_scores);
        
          obj.set_attractor(num, orig_attractor);
          if num == 1
              obj.r1 = orig_r;
          elseif num == 2
              obj.r2 = orig_r;
          end
      end
      
      function visualize_linear_component(obj, num, comp_num)
          % VISUALIZE_LINEAR_FORCING Plots the linear forcing
          % component of the NUMth SVD'd attractor.
          %
          % VISUALIZE_LINEAR_FORCING(num)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          figure; plot(arr(:, comp_num));
      end
      
      function [x1_res, x2_res] = lagged(~, x1, x2, lag)
          if lag >= 0
              x1_lagged = x1(1 + lag:end);
              x2_corrected = x2(1:end-lag);
              x1_res = x1_lagged;
              x2_res = x2_corrected;
          else
              lag = abs(lag);
              x2_lagged = x2(1 + lag:end);
              x1_corrected = x1(1:end-lag);
              x1_res = x1_corrected;
              x2_res = x2_lagged;
          end
      end
      
      function [X1_res, X2_res] = lagged_attractor(~, X1, X2, lag)
          if lag >= 0
              X1_lagged = X1(1 + lag:end, :);
              X2_corrected = X2(1:end-lag, :);
              X1_res = X1_lagged;
              X2_res = X2_corrected;
          else
              lag = abs(lag);
              X2_lagged = X2(1 + lag:end, :);
              X1_corrected = X1(1:end-lag, :);
              X1_res = X1_corrected;
              X2_res = X2_lagged;
          end
      end

      function W = exp_weight(~, X)
          norm = X(:, 1) + 0.0001;
          numer = exp(bsxfun(@rdivide, -X, norm));
          denom = sum(numer, 2);
          W = bsxfun(@rdivide, numer, denom);
      end
      
      function corr = linear_correlation(obj, lag)
          % LINEAR_CORRELATION Calculates the Pearson correlation between
          % the two time-series.
          %
          % LINEAR_CORRELATION()
          t1 = obj.x1;
          t2 = obj.x2;
          if nargin > 1
              [t1, t2] = obj.lagged(t1, t2, lag);
          end
          corr_ = corrcoef(t1, t2);
          corr = corr_(1, 2);
      end
      
      function [x1x2, x2x1] = granger_scores(obj)
          y = cat(2, obj.x1, obj.x2)';
          ARorder = 2;
          cmp = 2;
          normalize = 1;
          Aff=lsGC_analysis(y,cmp,ARorder,normalize);
          x1x2 = Aff(2, 1);
          x2x1 = Aff(1, 2);
      end
      
      function [F, c_v] = granger(obj, num_lags)
          t1 = obj.x1;
          t2 = obj.x2;
          [F, c_v] = granger_cause(t1, t2, 0.05, num_lags);
      end
      
      function SugiCorr = ccm_obsolete(obj, num_neighbs, linear_or_not, lag)
          % CCM Calculates Convergent Cross Mapping (CCM) scores using
          % NUM_NEIGHBS number of neighbors. LINEAR_OR_NOT is an array of
          % size two that indicates whether the first or the second
          % attractors are regular, linear, or nonlinear with 0, 1, or 2
          % respectively. NUM_NEIGHBS and LINEAR_OR_NOT are optionals.
          % NUM_NEIGHBS defaults ot E + 2 and the two attractors are
          % regular by default.
          %
          % CCM()
          %
          % CCM(NUM_NEIGHBS)
          %
          % CCM(NUM_NEIGHBS, LINEAR_OR_NOT)
          if nargin < 3
              X = obj.Xm;
              Y = obj.Ym;
              lag = 0;
          else
              if linear_or_not(1) == 1
                  X = obj.Xm_linear;
              elseif linear_or_not(1) == 2
                  X = obj.Xm_nonlinear;
              else
                  X = obj.Xm;
              end

              if linear_or_not(2) == 1
                  Y = obj.Ym_linear;
              elseif linear_or_not(2) == 2
                  Y = obj.Ym_nonlinear;
              else
                  Y = obj.Ym;
              end
          end

          if nargin == 1
              num_neighbs = obj.E1; % it's r if havok is on.
              lag = 0;
          end
          
          if nargin == 4
              [X, Y] = obj.lagged_attractor(obj.Xm, obj.Ym, lag);
          end
          
           
          [n1,d1]=knnsearch(X,X,'k',num_neighbs+2);
          [n2,d2]=knnsearch(Y,Y,'k',num_neighbs+2);
          
          n1 = n1(:, 2:end); %NNX
          d1 = d1(:, 2:end);
          n2 = n2(:, 2:end); %NNY
          d2 = d2(:, 2:end);

          X_size = size(X);
          dim = X_size(2);


          x1_p = zeros(size(X));
          x2_p = zeros(size(Y));
          
          W1 = obj.exp_weight(d1);
          W2 = obj.exp_weight(d2);
          
          for j=1:dim
              a1 = reshape(X(n2, j), size(n2));
              a2 = reshape(Y(n1, j), size(n1));
              
              b1 = a1 .* W2;
              b2 = a2 .* W1;
              c1 = sum(b1, 2);
              c2 = sum(b2, 2);

              x1_p(:, j) = c1;
              x2_p(:, j) = c2;
          end
          
          obj.x1_pred = x1_p;
          obj.x2_pred = x2_p;
          
          %dists1 = diag(pdist2(x1_p, X));
          %dists2 = diag(pdist2(x2_p, Y));

          sc1 = zeros(dim, 1);
          sc2 = zeros(dim, 1);

          for ii=1:dim
              p1 = x1_p(:,ii);
              p2 = x2_p(:,ii);

              corr1 = corrcoef(p1, X(:, ii));
              corr2 = corrcoef(p2, Y(:, ii));
              sc1(ii) = corr1(1, 2);
              sc2(ii) = corr2(1, 2);
          end

          SugiCorr = zeros(2, 1);

          SugiCorr(1, 1) = mean(sc1);
          SugiCorr(2, 1) = mean(sc2);
      end
      
      
       function SugiCorr = ccm(obj, linear_or_not, lag)
          % CCM Calculates Convergent Cross Mapping (CCM) scores using
          % NUM_NEIGHBS number of neighbors. LINEAR_OR_NOT is an array of
          % size two that indicates whether the first or the second
          % attractors are regular, linear, or nonlinear with 0, 1, or 2
          % respectively. NUM_NEIGHBS and LINEAR_OR_NOT are optionals.
          % NUM_NEIGHBS defaults ot E + 2 and the two attractors are
          % regular by default.
          %
          % CCM()
          %
          % CCM(NUM_NEIGHBS)
          %
          % CCM(NUM_NEIGHBS, LINEAR_OR_NOT)
          if nargin < 2
              X = obj.Xm;
              Y = obj.Ym;
              lag = 0;
          else
              if linear_or_not(1) == 1
                  X = obj.Xm_linear;
              elseif linear_or_not(1) == 2
                  X = obj.Xm_nonlinear;
              else
                  X = obj.Xm;
              end

              if linear_or_not(2) == 1
                  Y = obj.Ym_linear;
              elseif linear_or_not(2) == 2
                  Y = obj.Ym_nonlinear;
              else
                  Y = obj.Ym;
              end
          end

          if nargin == 1
              %num_neighbs = obj.E1; % it's r if havok is on.
              lag = 0;
          end
          
          if nargin == 3
              [X, Y] = obj.lagged_attractor(obj.Xm, obj.Ym, lag);
          end
          
          X_size = size(X);
          Y_size = size(Y);
          dim1 = X_size(2);
          dim2 = Y_size(2);
          
           
          [n1,d1]=knnsearch(X,X,'k',dim1+2);
          [n2,d2]=knnsearch(Y,Y,'k',dim2+2);
          
          n1 = n1(:, 2:end); %NNX
          d1 = d1(:, 2:end);
          n2 = n2(:, 2:end); %NNY
          d2 = d2(:, 2:end);


          x1_p = zeros(size(X));
          x2_p = zeros(size(Y));
          
          W1 = obj.exp_weight(d1);
          W2 = obj.exp_weight(d2);
          
          for j1=1:dim1
              a1 = reshape(X(n2, j1), size(n2));
              b1 = a1 .* W2;
              c1 = sum(b1, 2);
              x1_p(:, j1) = c1;
          end
          
          for j2=1:dim2
              a2 = reshape(Y(n1, j2), size(n1));
              b2 = a2 .* W1;
              c2 = sum(b2, 2);
              x2_p(:, j2) = c2;
          end
          
          %{
          for j=1:dim1
              a1 = reshape(X(n2, j), size(n2));
              a2 = reshape(Y(n1, j), size(n1));
              
              b1 = a1 .* W2;
              b2 = a2 .* W1;
              c1 = sum(b1, 2);
              c2 = sum(b2, 2);

              x1_p(:, j) = c1;
              x2_p(:, j) = c2;
          end
          %}
          
          obj.x1_pred = x1_p;
          obj.x2_pred = x2_p;
          
          %dists1 = diag(pdist2(x1_p, X));
          %dists2 = diag(pdist2(x2_p, Y));

          sc1 = zeros(dim1, 1);
          sc2 = zeros(dim2, 1);
          
          for ii1=1:dim1
              p1 = x1_p(:,ii1);
              corr1 = corrcoef(p1, X(:, ii1));
              sc1(ii1) = corr1(1, 2);
          end
          
          for ii2=1:dim2
               p2 = x2_p(:,ii2);
               corr2 = corrcoef(p2, Y(:, ii2));
               sc2(ii2) = corr2(1, 2);
          end
          
          %{

          for ii=1:dim
              p1 = x1_p(:,ii);
              p2 = x2_p(:,ii);

              corr1 = corrcoef(p1, X(:, ii));
              corr2 = corrcoef(p2, Y(:, ii));
              sc1(ii) = corr1(1, 2);
              sc2(ii) = corr2(1, 2);
          end
          
          %}

          SugiCorr = zeros(2, 1);

          SugiCorr(1, 1) = mean(sc1);
          SugiCorr(2, 1) = mean(sc2);
      end
      
      function [sc1_l_mean, sc1_nl_mean, sc2_l_mean, sc2_nl_mean] = ...
              linear_and_nonlinear_ccm_scores(obj)
          X = obj.Xm;
          Y = obj.Ym;
          X_size = size(X);
          dim = X_size(2);
          x1_p = obj.x1_pred;
          x2_p = obj.x2_pred;
          sc1_l = zeros(dim, 1);
          sc1_nl = zeros(dim, 1);
          sc2_l = zeros(dim, 1);
          sc2_nl = zeros(dim, 1);
          Xm_linear_indices = obj.Xm_linear_points;
          Xm_nonlinear_indices = obj.Xm_nonlinear_points;
          Ym_linear_indices = obj.Ym_linear_points;
          Ym_nonlinear_indices = obj.Ym_nonlinear_points;
          
          for ii=1:dim
              p1_l = x1_p(Xm_linear_indices, ii);
              p1_nl = x1_p(Xm_nonlinear_indices, ii);
              p2_l = x2_p(Ym_linear_indices, ii);
              p2_nl = x2_p(Ym_nonlinear_indices, ii);
              
              corr1_l = corrcoef(p1_l, X(Xm_linear_indices, ii));
              corr1_nl = corrcoef(p1_nl, X(Xm_nonlinear_indices, ii));
              corr2_l = corrcoef(p2_l, Y(Ym_linear_indices, ii));
              corr2_nl = corrcoef(p2_nl, Y(Ym_nonlinear_indices, ii));
              
              sc1_l(ii) = corr1_l(1, 2);
              sc1_nl(ii) = corr1_nl(1, 2);
              sc2_l(ii) = corr2_l(1, 2);
              sc2_nl(ii) = corr2_nl(1, 2);
          end
          
          sc1_l_mean = mean(sc1_l);
          sc1_nl_mean = mean(sc1_nl);
          sc2_l_mean = mean(sc2_l);
          sc2_nl_mean = mean(sc2_nl);
      end
      
      function [sc1_self_l_mean, sc1_self_nl_mean, ...
               sc1_other_l_mean, sc1_other_nl_mean,...
               sc2_self_l_mean, sc2_self_nl_mean, ...
               sc2_other_l_mean, sc2_other_nl_mean] = ...
              all_linear_and_nonlinear_ccm_scores(obj)
          X = obj.Xm;
          Y = obj.Ym;
          X_size = size(X);
          Y_size = size(Y);
          dim1 = X_size(2);
          dim2 = Y_size(2);
          x1_p = obj.x1_pred;
          x2_p = obj.x2_pred;
          sc1_self_l = zeros(dim1, 1);
          sc1_self_nl = zeros(dim1, 1);
          sc1_other_l = zeros(dim1, 1);
          sc1_other_nl = zeros(dim1, 1);
          sc2_self_l = zeros(dim2, 1);
          sc2_self_nl = zeros(dim2, 1);
          sc2_other_l = zeros(dim2, 1);
          sc2_other_nl = zeros(dim2, 1);
          Xm_linear_indices = obj.Xm_linear_points;
          Xm_nonlinear_indices = obj.Xm_nonlinear_points;
          Ym_linear_indices = obj.Ym_linear_points;
          Ym_nonlinear_indices = obj.Ym_nonlinear_points;
          
           for ii=1:dim1
              p1_self_l = x1_p(Xm_linear_indices, ii);
              p1_self_nl = x1_p(Xm_nonlinear_indices, ii);
              corr1_self_l = corrcoef(p1_self_l, X(Xm_linear_indices, ii));
              corr1_self_nl = corrcoef(p1_self_nl, X(Xm_nonlinear_indices, ii));
              sc1_self_l(ii) = corr1_self_l(1, 2);
              sc1_self_nl(ii) = corr1_self_nl(1, 2);
              
              p1_other_l = x1_p(Ym_linear_indices, ii);
              p1_other_nl = x1_p(Ym_nonlinear_indices, ii);
              corr1_other_l = corrcoef(p1_other_l, X(Ym_linear_indices, ii));
              corr1_other_nl = corrcoef(p1_other_nl, X(Ym_nonlinear_indices, ii));
              sc1_other_l(ii) = corr1_other_l(1, 2);
              sc1_other_nl(ii) = corr1_other_nl(1, 2);
           end
           
           for ii = 1:dim2
              p2_self_l = x2_p(Ym_linear_indices, ii);
              p2_self_nl = x2_p(Ym_nonlinear_indices, ii);
              corr2_self_l = corrcoef(p2_self_l, Y(Ym_linear_indices, ii));
              corr2_self_nl = corrcoef(p2_self_nl, Y(Ym_nonlinear_indices, ii));
              sc2_self_l(ii) = corr2_self_l(1, 2);
              sc2_self_nl(ii) = corr2_self_nl(1, 2);
              
              p2_other_l = x2_p(Xm_linear_indices, ii);
              p2_other_nl = x2_p(Xm_nonlinear_indices, ii);
              corr2_other_l = corrcoef(p2_other_l, Y(Xm_linear_indices, ii));
              corr2_other_nl = corrcoef(p2_other_nl, Y(Xm_nonlinear_indices, ii));
              sc2_other_l(ii) = corr2_other_l(1, 2);
              sc2_other_nl(ii) = corr2_other_nl(1, 2);
           end
          
          sc1_self_l_mean = mean(sc1_self_l);
          sc1_self_nl_mean = mean(sc1_self_nl);
          sc1_other_l_mean = mean(sc1_other_l);
          sc1_other_nl_mean = mean(sc1_other_nl);
          sc2_self_l_mean = mean(sc2_self_l);
          sc2_self_nl_mean = mean(sc2_self_nl);
          sc2_other_l_mean = mean(sc2_other_l);
          sc2_other_nl_mean = mean(sc2_other_nl);
      end
    end
end
