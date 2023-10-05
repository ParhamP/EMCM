classdef EMCM < handle
    % EMCM (Eigen Manifold Cross Mapping) Class for eigen time delaying
    % time series as high dimensional manifolds and assessing how well
    % they predict each other.
    % This class is able to get two time-series, embed them in higher
    % dimensional attractors, and apply cross-mapping. There are functions
    % for visualization, as well as tools for estimating relevant
    % parameters. 
    properties %(Access=private)
        x1;
        x2;
        x1_name;
        x2_name;
        Xm;
        Ym;
        Xm_latents;
        Ym_latents;
        Xm_eigenvalues;
        Ym_eigenvalues;
        Xm_thresh;
        Ym_thresh;
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
        n1;
        n2;
        y1;
        y2;
        corrs_1;
        corrs_2;
    end
    
    methods
      function obj = EMCM()
          % Constructs an AttractorAnalysis object.
          %
          % OBJ = AttractorAnalysis()
          return
      end
      
      function set_timeseries(obj, num, x, varargin)
          % SET_TIMESERIES sets the NUMth time-series X.
          %
          % SET_TIMESERIES(NUM, X)
          p = inputParser;
          paramName = 'name';
          defaultName = "x1";
          addParameter(p,paramName,defaultName);
          parse(p, varargin{:});
          if num == 1
              obj.x1 = x;
              obj.x1_name = p.Results.name;
          elseif num == 2
              obj.x2 = x;
              obj.x2_name = p.Results.name;
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
          min_val = min([min(obj.x1), min(obj.x2)]) - 0.5;
          max_val = max([max(obj.x1), max(obj.x2)]) + 0.5;
          f = figure;
%           sgtitle('Time-Series');
          subplot(2,1,1)
          plot(obj.x1, 'k');
          ylim([min_val max_val]);
          title(obj.x1_name);
          
          subplot(2,1,2)
          plot(obj.x2, 'k');
          ylim([min_val, max_val]);
          title(obj.x2_name);
          f.Position = [354   358   890   535];
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
      
       function [Xm_neighbors, Ym_neighbors] = ...
                                      get_nearest_neighbors_3d(obj, index)
          Xm = obj.Xm;
          Ym = obj.Ym;
          n1 = obj.n1;
          n2 = obj.n2;
          Xm_nearests = n1(index, :);
          Ym_nearests = n2(index, :);
          
          Xm_neighb_1 = Xm(Xm_nearests(1), :);
          Xm_neighb_2 = Xm(Xm_nearests(2), :);
          Xm_neighb_3 = Xm(Xm_nearests(3), :);
          Xm_neighb_4 = Xm(Xm_nearests(4), :);
          
          Ym_neighb_1 = Ym(Ym_nearests(1), :);
          Ym_neighb_2 = Ym(Ym_nearests(2), :);
          Ym_neighb_3 = Ym(Ym_nearests(3), :);
          Ym_neighb_4 = Ym(Ym_nearests(4), :);
          
          Xm_neighbors = {Xm_neighb_1, Xm_neighb_2, Xm_neighb_3, ...
                          Xm_neighb_4};
          Ym_neighbors = {Ym_neighb_1, Ym_neighb_2, Ym_neighb_3, ...
                          Ym_neighb_4};
       end
      
      function [Xm_neighbors, Ym_neighbors] = ...
                                      get_nearest_neighbors_2d(obj, index)
          Xm = obj.Xm;
          Ym = obj.Ym;
          n1 = obj.n1;
          n2 = obj.n2;
          Xm_nearests = n1(index, :);
          Ym_nearests = n2(index, :);
          
          Xm_neighb_1 = Xm(Xm_nearests(1), :);
          Xm_neighb_2 = Xm(Xm_nearests(2), :);
          Xm_neighb_3 = Xm(Xm_nearests(3), :);
          
          Ym_neighb_1 = Ym(Ym_nearests(1), :);
          Ym_neighb_2 = Ym(Ym_nearests(2), :);
          Ym_neighb_3 = Ym(Ym_nearests(3), :);
          
          Xm_neighbors = {Xm_neighb_1, Xm_neighb_2, Xm_neighb_3};
          Ym_neighbors = {Ym_neighb_1, Ym_neighb_2, Ym_neighb_3};
      end

      function visualize_neighbors(obj, num, point_number, dest)
          if num == 1
              neighb_indices1 = obj.n1;
              neighb_indices2 = obj.n2;
              arr1 = obj.Xm;
              arr2 = obj.Ym;
              arr1_pred_arr2 = obj.x2_pred;
              arr2_pred_arr1 = obj.x1_pred;
          elseif num == 2
              neighb_indices1 = obj.n2;
              neighb_indices2 = obj.n1;
              arr1 = obj.Ym;
              arr2 = obj.Xm;
              arr1_pred_arr2 = obj.x1_pred;
              arr2_pred_arr1 = obj.x2_pred;
          end
          
%           current_neighb_indices = neighb_indices(point_number, :);
          
          f = figure('visible','off');
          f.Position = [100, 100, 560 * 2, 560 * 2];
%           sgtitle('Attractors');
          subplot(2,1,1)
          
          plot3(arr1(:, 1), arr1(:, 2), arr1(:, 3)); hold on;
          
          current_neighb_indices_1 = neighb_indices1(point_number, :);
          current_pred_indices_1 = neighb_indices2(point_number, :);
          current_arr2_pred_arr1 = arr2_pred_arr1(point_number, :);

          
%           pattern1 = arr1(current_neighb_indices, :);
          m_size = 5;
          
          plot3(arr1(point_number, 1), arr1(point_number, 2), arr1(point_number, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;

          
%           for i=1:length(current_neighb_indices_1)
%               n = current_neighb_indices_1(i);
%               plot3(arr1(n, 1), arr1(n, 2), arr1(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
%           end
          
          for i=1:length(current_pred_indices_1)
              n = current_pred_indices_1(i);
              plot3(arr1(n, 1), arr1(n, 2), arr1(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor',"#EDB120",'MarkerEdgeColor',"#EDB120"); hold on;
          end
          
          plot3(current_arr2_pred_arr1(1), current_arr2_pred_arr1(2), current_arr2_pred_arr1(3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
          
          plot3([arr1(point_number, 1) current_arr2_pred_arr1(1)], [arr1(point_number, 2) current_arr2_pred_arr1(2)], [arr1(point_number, 3) current_arr2_pred_arr1(3)], 'LineWidth', 3); hold on;
%           plot3(pattern1(:, 1), pattern1(:, 2), pattern1(:, 3)); xlim([-0.2 0.2]); ...
%               ylim([-0.2 0.2]); zlim([-0.2 0.2]); hold on;
          title('Alpha')
          subplot(2,1,2)
          
          plot3(arr2(:, 1), arr2(:, 2), arr2(:, 3)); hold on;
          
          current_neighb_indices_2 = neighb_indices2(point_number, :);
          current_pred_indices_2 = neighb_indices1(point_number, :);
          current_arr1_pred_arr2 = arr1_pred_arr2(point_number, :);
          
%           pattern2 = arr2(current_neighb_indices, :);

          plot3(arr2(point_number, 1), arr2(point_number, 2), arr2(point_number, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
          
%           for i=1:length(current_neighb_indices_2)
%               n = current_neighb_indices_2(i);
%               plot3(arr2(n, 1), arr2(n, 2), arr2(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
%           end
          
          for i=1:length(current_pred_indices_2)
              n = current_pred_indices_2(i);
              plot3(arr2(n, 1), arr2(n, 2), arr2(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor',"#EDB120",'MarkerEdgeColor',"#EDB120"); hold on;
          end
          
          plot3(current_arr1_pred_arr2(1), current_arr1_pred_arr2(2), current_arr1_pred_arr2(3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
          
          plot3([arr2(point_number, 1) current_arr1_pred_arr2(1)], [arr2(point_number, 2) current_arr1_pred_arr2(2)], [arr2(point_number, 3) current_arr1_pred_arr2(3)], 'LineWidth', 3, 'Color', "#77AC30"); hold on;

          f.Position = [283    66   991   886];
          title('Gamma')
%           disp('')
%           plot3(pattern2(:, 1), pattern2(:, 2), pattern2(:, 3)); xlim([-0.2 0.2]); ...
%               ylim([-0.2 0.2]); zlim([-0.2 0.2]); hold on;
%           scatter3(pattern2(:, 1), pattern2(:, 2), pattern2(:, 3)); xlim([-0.2 0.2]); ...
%               ylim([-0.2 0.2]); zlim([-0.2 0.2]); hold on;
          saveas(gcf, [dest num2str(point_number) '.png']);
      end
      
       function visualize_neighbors_separate(obj, num, point_number, dest)
          if num == 1
              neighb_indices1 = obj.n1;
              neighb_indices2 = obj.n2;
              arr1 = obj.Xm;
              arr2 = obj.Ym;
              arr1_pred_arr2 = obj.x2_pred;
              arr2_pred_arr1 = obj.x1_pred;
          elseif num == 2
              neighb_indices1 = obj.n2;
              neighb_indices2 = obj.n1;
              arr1 = obj.Ym;
              arr2 = obj.Xm;
              arr1_pred_arr2 = obj.x1_pred;
              arr2_pred_arr1 = obj.x2_pred;
          end
          
%           current_neighb_indices = neighb_indices(point_number, :);
          
          f = figure('visible','off');
%           f = figure;
          f.Position = [100, 100, 560 * 2, 560 * 2];
          sgtitle(num2str(point_number));
          subplot(2,1,1)
          
          plot3(arr1(:, 1), arr1(:, 2), arr1(:, 3)); hold on;
          
          current_neighb_indices_1 = neighb_indices1(point_number, :);
          current_pred_indices_1 = neighb_indices2(point_number, :);
          current_arr2_pred_arr1 = arr2_pred_arr1(point_number, :);

          
%           pattern1 = arr1(current_neighb_indices, :);
          m_size = 5;
          
          plot3(arr1(point_number, 1), arr1(point_number, 2), arr1(point_number, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;

          
          for i=1:length(current_neighb_indices_1)
              n = current_neighb_indices_1(i);
              plot3(arr1(n, 1), arr1(n, 2), arr1(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
          end
          
%           for i=1:length(current_pred_indices_1)
%               n = current_pred_indices_1(i);
%               plot3(arr1(n, 1), arr1(n, 2), arr1(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor',"#EDB120",'MarkerEdgeColor',"#EDB120"); hold on;
%           end
          
%           plot3(current_arr2_pred_arr1(1), current_arr2_pred_arr1(2), current_arr2_pred_arr1(3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
          
%           plot3([arr1(point_number, 1) current_arr2_pred_arr1(1)], [arr1(point_number, 2) current_arr2_pred_arr1(2)], [arr1(point_number, 3) current_arr2_pred_arr1(3)], 'LineWidth', 3); hold on;

          title('Alpha')
          subplot(2,1,2)
          
          plot3(arr2(:, 1), arr2(:, 2), arr2(:, 3), 'Color', '#e35919'); hold on;
          
          current_neighb_indices_2 = neighb_indices2(point_number, :);
          current_pred_indices_2 = neighb_indices1(point_number, :);
          current_arr1_pred_arr2 = arr1_pred_arr2(point_number, :);
          

          plot3(arr2(point_number, 1), arr2(point_number, 2), arr2(point_number, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on;
          
%           for i=1:length(current_neighb_indices_2)
%               n = current_neighb_indices_2(i);
%               plot3(arr2(n, 1), arr2(n, 2), arr2(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor',"#4DBEEE",'MarkerEdgeColor',"#4DBEEE"); hold on;
%           end
          
          for i=1:length(current_pred_indices_2)
              n = current_pred_indices_2(i);
              plot3(arr2(n, 1), arr2(n, 2), arr2(n, 3), 'o','MarkerSize',m_size, 'MarkerFaceColor',"#EDB120",'MarkerEdgeColor',"#EDB120"); hold on;
          end
          
          plot3(current_arr1_pred_arr2(1), current_arr1_pred_arr2(2), current_arr1_pred_arr2(3), 'o','MarkerSize',m_size, 'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on;
          
          plot3([arr2(point_number, 1) current_arr1_pred_arr2(1)], [arr2(point_number, 2) current_arr1_pred_arr2(2)], [arr2(point_number, 3) current_arr1_pred_arr2(3)], 'LineWidth', 3, 'Color', "#77AC30"); hold on;

          f.Position = [580    77   621   886];
          title('Gamma')
%           disp('')
          f.Renderer = 'painters';
          if ~strcmp(dest, '')
            saveas(gcf, [dest num2str(point_number) '.svg']);
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
      
      function [est_dimension, fnn_scores] = estimate_dimension_using_fnn(obj, num, tau, ...
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
      
      function attractor_subset = get_subset_attractor(obj, num, num_dims)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          attractor_subset = arr(:, 1:num_dims);
      end
      
      function est_dimension = estimate_dimension_using_fnn_special_case(obj, num, ...
                                                            max_dim)
          % ESTIMATE_DIMENSION_USING_FNN Estimates the optimal dimension
          % for the NUMth timeseries using false nearest neighbor
          % algorithm. It embeds it with TAU and uses dimensions up to
          % MAX_DIM.
          %
          %ESTIMATE_DIMENSION_USING_FNN(NUM, TAU, MAX_DIM)
          
          tau = 1;
          rtol = 15.0;
          atol = 2.0;
          min_dim = 1;
          dim_range = min_dim:max_dim;
          
          x_ts = obj.get_timeseries(num);
          
          %c = AttractorAnalysis();
          %c.set_timeseries(num, x_ts);
          
          fnn_scores = zeros(1, max_dim);
          
          obj.generate_single_attractor(num, tau, max_dim + 1, 0, nan);
          
          for d=dim_range
            x_d = obj.get_subset_attractor(num, d);
            x_d_1 = obj.get_subset_attractor(num, d + 1);
            [~, ~, f3] = single_dimension_fnn(x_ts, x_d, x_d_1, rtol, atol);
            fnn_scores(d) = f3;
          end
          [~, i] = min(fnn_scores); 
          est_dimension = dim_range(i);
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
              obj.Xm = X;
          elseif num == 2
              y = obj.x2;
              Y = obj.embeder(y, tau, E);
              [~,~,latent,~,~] = pca(Y);
              obj.Ym = Y;
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
      
      function generate_single_attractor(obj, num, tau, E)
          %addpath('./HAVOK/utils');
          % GENERATE_SINGLE_ATTRACTOR Generates the NUMth attractor with
          % tau of TAU and dimension of E. Binary USE_SVD indicates whether
          % or not to use SVD on the constructed attractor. If svd is
          % applied, the method picks R components of U. 
          if num == 1
              obj.E1 = E;
              x = obj.x1;
              X = obj.embeder(x, tau, E);
              obj.Xm = X;
          elseif num == 2
              obj.E2 = E;
              y = obj.x2;
              Y = obj.embeder(y, tau, E);
              obj.Ym = Y;
          end
      end
      
      function SVD_on_already_attractor(obj, num, r)
          if num == 1
              obj.r1 = r;
              X = obj.Xm;
              [Ux,~,~] = svd(X,'econ');
              X = Ux(1:end,1:r);
              obj.Xm = X;
          elseif num == 2
              obj.r2 = r;
              Y = obj.Ym;
              [Uy,~,~] = svd(Y,'econ');
              Y = Uy(1:end,1:r);
              obj.Ym = Y;
          end
      end
      
      function thresh = find_optimal_threshold(obj, arr, latents)
          beta = size(arr,2)/size(arr,1);
          thresh = optimal_SVHT_coef(beta,0) * median(latents);
          %r = length(sigs(sigs>thresh));
      end
      
      function [r_est, thresh] = apply_pca_and_set_r(obj, num, varargin)
          p = inputParser;
          thresh_method_param = 'threshold_method';
          thresh_method_param_default = 'all';
          addParameter(p, thresh_method_param, thresh_method_param_default);
          parse(p, varargin{:});
          thresh_method = p.Results.threshold_method;
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          arr = zscore(arr);
          [~,U,latent,~] = pca(arr);
          eigenvalues = latent;
          if strcmp(thresh_method, 'gavish')
              thresh = obj.find_optimal_threshold(U, latent);
          elseif strcmp(thresh_method, 'one')
              thresh = 1;
          elseif strcmp(thresh_method, '90_percent')
              total_var_sum = sum(latent);
               for i=1:length(latent)
                   sofar_var_sum = sum(latent(1:i));
                   var_ratio = sofar_var_sum / total_var_sum;
                   if var_ratio >= 0.9
                       thresh = latent(i);
                       break;
                   end
               end
          elseif strcmp(thresh_method, 'all')
              thresh = -Inf;
          else
              error('Wrong unknown threshold method name.');
              %thresh = -Inf;
          end
          r_est = length(latent(latent>thresh));
          arr = U(1:end, 1:r_est);
          latent = latent(1:r_est);
          if num == 1
              obj.r1 = r_est;
              obj.Xm_latents = latent;
              obj.Xm_eigenvalues = eigenvalues;
              obj.Xm_thresh = thresh;
              obj.Xm = arr;
          elseif num == 2
              obj.r2 = r_est;
              obj.Ym_latents = latent;
              obj.Ym_eigenvalues = eigenvalues;
              obj.Ym_thresh = thresh;
              obj.Ym = arr;
          end
      end
      
      function y = linear_reconstruction(obj, num)
          if num == 1
              V = obj.Xm;
              r = obj.r1;
          else
              V = obj.Ym;
              r = obj.r2;
          end
          dt = .001; %.001
          x = V(1:end-1,1:r);
          xprime = V(2:end,1:r);

          Xi = xprime\x;
          B = Xi(1:r-1,r);
          A = Xi(1:r-1,1:r-1);
          sys = ss(A,B,eye(r-1),0*B,dt);

          cut_num = 100;
          L = 1:length(x);
          [y,~] = lsim(sys,x(L,r),dt*(L-1),x(1,1:r-1));
          y = zscore(y);
          y = y(cut_num:end, :);
%           figure; plot3(y(:, 1), y(:, 2), y(:, 3));
%           view(270,0)
          %figure; plot(y(:, 1), y(:, 2));
          if num == 1
              obj.Xm = V(cut_num+1:end, :);
              obj.y1 = y;
          else
              obj.Ym = V(cut_num+1:end, :);
              obj.y2 = y;
          end
      end
      
      
      
      function shuffle_columns_attractor(obj, num)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          [m,n] = size(arr);
          b = arr;
          for i = 1:n
              idx = randperm(m);
              b(idx, i) = arr(:, i) ;
          end
          if num == 1
              obj.Xm = b;
          elseif num == 2
              obj.Ym = b;
          end
      end
      
      function shuffle_rows_attractor(obj, num)
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          end
          [m,n] = size(arr);
          b = arr;
          for i = 1:m
              idx = randperm(n);
              b(i, idx) = arr(i, :) ;
          end
          if num == 1
              obj.Xm = b;
          elseif num == 2
              obj.Ym = b;
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
    
      function [all_points_arr, linear_points] = generate_linear_dynamics(obj, num, thresh, hit_range)
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
      
      function [all_points_arr, nonlinear_points] = generate_nonlinear_dynamics(obj, num, thresh, hit_range)
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
              f.Position(3) * 1, f.Position(4) * 1];
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
           view(270,0)
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
%            if nargin > 2
%                 f = figure('visible','off');  
%            else
%                 f = figure;
%            end
           f = figure;
          
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          else
              error('wrong num');
          end
          
          plot3(arr(:, 1), arr(:, 2), arr(:, 3), 'LineWidth', 1, 'Color', 'k');
           axis tight
           set(gca,'xtick',[])
           set(gca,'xticklabel',[])
           set(gca,'ytick',[])
           set(gca,'yticklabel',[])
           set(gca,'ztick',[])
           set(gca,'zticklabel',[])
%            view([360 360 360])
           %xlim([0 2]);
           %ylim([0 2]);
           %zlim([0 2]);
           axis off
           view(53,33)
%            view(125, 20)
          if nargin > 2
              title(my_title);
              f.Renderer = 'painters';
              saveas(gcf, dest);
%               close(gcf);
          end
      end
      
           function visualize_single_attractor_2d(obj, num, my_title, dest)
          % VISUALIZE_SINGLE_ATTRACTOR Plots the NUMth attractor in 3D. If
          % MY_TITLE and DEST are given, the figure will be saved in DEST
          % with title MY_TITLE.
          %
          % VISUALIZE_SINGLE_ATTRACTOR(NUM)
          %
          % VISUALIZE_SINGLE_ATTRACTOR(NUM, MY_TITLE, DEST)
%            if nargin > 2
%                 f = figure('visible','off');  
%            else
%                 f = figure;
%            end
           figure;
          
          if num == 1
              arr = obj.Xm;
          elseif num == 2
              arr = obj.Ym;
          else
              error('wrong num');
          end
          
          plot(arr(:, 1), arr(:, 2));
           axis tight
%           axis off
%           view(270,0)
          if nargin > 2
              title(my_title);
              saveas(gcf, dest);
              close(gcf);
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
%           sgtitle('Attractors');
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
          title(obj.x1_name);
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
          title(obj.x2_name);
          f.Position = [514    61   611   896];
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
      
      function visualize_predicted_attractors(obj, d1, d2, d3)
          % VISUALIZE_PREDICTED_ATTRACTORS Plots the prediction points the
          % attractors make of each other on top of the original attractors
          % they would be compared to in CCM.
          %
          % VISUALIZE_PREDICTED_ATTRACTORS()
          if nargin < 2
              d1 = 1;
              d2 = 2;
              d3 = 3;
          end
          %obj.Xm = zscore(obj.Xm);
          %obj.Ym = zscore(obj.Ym);
          %obj.x1_pred = zscore(obj.x1_pred);
          %obj.x2_pred = zscore(obj.x2_pred);
          f = figure;
          f.Position = [100, 100, 560 * 2, 560 * 2];
          sgtitle('Attractors');
          subplot(2,1,1)
          scatter3(obj.x1_pred(:, d1), obj.x1_pred(:, d2), obj.x1_pred(:, d3), 'k');
          set(gca,'box','off')
          hold on;
          plot3(obj.Xm(:, d1), obj.Xm(:, d2), obj.Xm(:, d3), 'r', 'LineWidth', 2);
          subplot(2,1,2)
          scatter3(obj.x2_pred(:, d1), obj.x2_pred(:, d2), obj.x2_pred(:, d3), 'k');
          set(gca,'box','off')
          hold on;
          plot3(obj.Ym(:, d1), obj.Ym(:, d2), obj.Ym(:, d3), 'r', 'LineWidth', 2);
          set(gca,'box','off')
      end
      
      function visualize_predicted_attractors_seperated(obj, d1, d2, d3, dest)
          if strcmp(dest, '')
              dest = '';
          end
          % VISUALIZE_PREDICTED_ATTRACTORS Plots the prediction points the
          % attractors make of each other on top of the original attractors
          % they would be compared to in CCM.
          %
          % VISUALIZE_PREDICTED_ATTRACTORS()
          if nargin < 2
              d1 = 1;
              d2 = 2;
              d3 = 3;
          end
          %obj.Xm = zscore(obj.Xm);
          %obj.Ym = zscore(obj.Ym);
          %obj.x1_pred = zscore(obj.x1_pred);
          %obj.x2_pred = zscore(obj.x2_pred);
          f = figure;
          f.Position = [100, 100, 560 * 2, 560 * 2];
          sgtitle('Attractors');
          subplot(2,2,1)
          plot3(obj.Xm(:, d1), obj.Xm(:, d2), obj.Xm(:, d3), '.-');
          subplot(2,2,2)
          plot3(obj.x1_pred(:, d1), obj.x1_pred(:, d2), obj.x1_pred(:, d3), '.-'); %, '.-'
          %set(gca,'box','off')
          %hold on;
          subplot(2,2,3)
          plot3(obj.Ym(:, d1), obj.Ym(:, d2), obj.Ym(:, d3), '.-', 'Color', '#e35919');
          subplot(2,2,4)
          plot3(obj.x2_pred(:, d1), obj.x2_pred(:, d2), obj.x2_pred(:, d3), '.-',  'Color', '#e35919');
          %set(gca,'box','off')
          %hold on;
          %set(gca,'box','off')
          f.Renderer = 'painters';
          if ~strcmp(dest, '')
            saveas(gcf, [dest 'attractors_preds_sidebyside' '.svg']);
          end
      end
      
      function visualize_predicted_attractors_2d(obj, d1, d2)
          % VISUALIZE_PREDICTED_ATTRACTORS Plots the prediction points the
          % attractors make of each other on top of the original attractors
          % they would be compared to in CCM.
          %
          % VISUALIZE_PREDICTED_ATTRACTORS()
          if nargin < 2
              d1 = 1;
              d2 = 2;
          end
          f = figure;
          f.Position = [100, 100, 560 * 2, 560 * 2];
          sgtitle('Attractors');
          subplot(2,1,1)
          scatter(obj.x1_pred(:, d1), obj.x1_pred(:, d2), 'k');
          hold on;
          plot(obj.Xm(:, d1), obj.Xm(:, d2), 'r', 'LineWidth', 2);
          subplot(2,1,2)
          scatter(obj.x2_pred(:, d1), obj.x2_pred(:, d2), 'k');
          hold on;
          plot(obj.Ym(:, d1), obj.Ym(:, d2), 'r', 'LineWidth', 2);
      end
      
      function visualize_nonlinear_forcing(obj, num)
          % VISUALIZE_NONLINEAR_FORCING Plots the nonlinear forcing
          % component of the NUMth SVD'd attractor.
          %
          % VISUALIZE_NONLINEAR_FORCING(num)
          if num == 1
              arr = obj.Xm;
              if obj.r1
                r = obj.r1;
              else
                  r = size(arr, 2);
              end
          elseif num == 2
              arr = obj.Ym;
              if obj.r1
                r = obj.r2;
              else
                  r = size(arr, 2);
              end
          end
          figure; plot(arr(:, r)); ylim([-5 5]);
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
          if lag == 0
              X1_res = X1;
              X2_res = X2;
              return;
          end
          if lag > 0
              X1_lagged = X1(1 + lag:end, :);
              X2_corrected = X2(1:end-lag, :);
              X1_res = X1_lagged;
              X2_res = X2_corrected;
          elseif lag < 0
              lag = abs(lag);
              X2_lagged = X2(1 + lag:end, :);
              X1_corrected = X1(1:end-lag, :);
              X1_res = X1_corrected;
              X2_res = X2_lagged;
          end
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
      
      function apply_random_coordinatess(obj)
          X = obj.Xm;
          m = size(X, 2);
          R = normrnd(0, 1, m, m);
          X = R * X';
          X = X';
          obj.Xm = X;
          
          Y = obj.Ym;
          m = size(Y, 2);
          R = normrnd(0, 1, m, m);
          Y = R * Y';
          Y = Y';
          obj.Ym = Y;
          
      end
      
      
       function SugiCorr = ccm(obj, varargin)
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
          
          
          X = obj.Xm;
          Y = obj.Ym;
          
          X_size = size(X);
          Y_size = size(Y);
          dim1 = X_size(2);
          dim2 = Y_size(2);
          
          p = inputParser;
          bar_param = 'barPlot';
          bar_param_default = false;
          lag_param = 'lag';
          lag_param_default = 0;
          use_dimension_weights = 'weigh_by_eigens';
          use_dimension_weights_default = false;
          num_neighbs_param = 'num_neighbors';
          num_neighbs_default = [dim1 + 2, dim2 + 2];
          
          addParameter(p, bar_param, bar_param_default);
          addParameter(p, lag_param, lag_param_default);
          addParameter(p, use_dimension_weights, ...
                       use_dimension_weights_default);
          addParameter(p, num_neighbs_param, num_neighbs_default);

          parse(p, varargin{:});
          draw_bar_plot = p.Results.barPlot;
          lag = p.Results.lag;
          use_dimension_weights = p.Results.weigh_by_eigens;
          num_neighbs = p.Results.num_neighbors;
          
          [X, Y] = obj.lagged_attractor(X, Y, lag);
          
          num_neighbs_1 = num_neighbs(1);
          num_neighbs_2 = num_neighbs(2);
          
          [n1,d1]=knnsearch(X,X,'k', num_neighbs_1, 'Distance', 'euclidean'); % each point has dim1 + 1 timestamps associated with it
          [n2,d2]=knnsearch(Y,Y,'k', num_neighbs_2, 'Distance', 'euclidean');
          
          %[n1,d1]=knnsearch(X,X,'k', dim1 + 2);
          %[n2,d2]=knnsearch(Y,Y,'k', dim2 + 2);
          
          n1 = n1(:, 2:end); %NNX
          d1 = d1(:, 2:end);
          n2 = n2(:, 2:end); %NNY
          d2 = d2(:, 2:end);
          
          obj.n1 = n1;
          obj.n2 = n2;

          x1_p = zeros(size(X));
          x2_p = zeros(size(Y));
          
          W1 = exp_weight(d1);
          W2 = exp_weight(d2);
          
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
          
          %x1_p = zscore(x1_p);
          %x2_p = zscore(x2_p);
          
          obj.x1_pred = x1_p;
          obj.x2_pred = x2_p;
          

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

          SugiCorr = zeros(2, 1);
          
          obj.corrs_1 = sc1;
          obj.corrs_2 = sc2;
          
          if use_dimension_weights
              weights_X = obj.Xm_latents;
              weights_Y = obj.Ym_latents;
              SugiCorr(1, 1) = wmean(sc1, weights_X);
              SugiCorr(2, 1) = wmean(sc2, weights_Y);
          else
              SugiCorr(1, 1) = mean(sc1);
              SugiCorr(2, 1) = mean(sc2); 
          end
          
          if draw_bar_plot == true
              figure;
              categs = {[obj.x1_name ' predicts ' obj.x2_name],...
                  [obj.x2_name ' predicts ' obj.x1_name]};
              X = categorical(categs);
              X = reordercats(X, categs);
              Y = [mean(sc2) mean(sc1)];
              bar(X,Y, 'k');
              title("Geometric Cross-Parameter Coupling");
          end
          
       end
       
       function errs = calculate_prediction_errors(obj, num)
           if num == 1
               arr = obj.Xm;
               pred = obj.x1_pred;
           else
               arr = obj.Ym;
               pred = obj.x2_pred;
           end
           num_dims = size(arr, 2);
           errs = NaN(1, num_dims);
           for i = 1:num_dims
               dim_arr = arr(:, i);
               dim_pred = pred(:, i);
               err = immse(dim_arr, dim_pred);
               errs(i) = err;
           end
       end
      
       function visualize_predictions_dim_by_dim(obj, num, dest)
           if strcmp(dest, '')
               dest = '';
           end
           if num == 1
               arr = obj.Xm;
               arr_size = size(arr);
               dim = arr_size(2);
               x_p = obj.x1_pred;
           else
               arr = obj.Ym;
               arr_size = size(arr);
               dim = arr_size(2);
               x_p = obj.x2_pred;
           end
           f = figure;
           num_rows = ceil(dim / 3);
           tlo = tiledlayout(num_rows,3,'TileSpacing','none','Padding','none');
           for ii1=1:dim
               p = x_p(:,ii1);
               nexttile, plot(arr(:, ii1)); hold on; plot(p);
               ylim([-9 9]);
               cur_dim = arr(:, ii1);
               rw_res = vratiotest(cur_dim);
               if rw_res == 1
                   random_walk_bool = 'false';
               else
                   random_walk_bool = 'true';
               end
               l_exp = lyapunovExponent(cur_dim);
%                nexttile, scatter(p, arr(:, ii1));
%                ylim([-5 5]);
               %err = immse(p, arr(:, ii1)); 
               corr = corrcoef(p, arr(:, ii1));
               corr = corr(1, 2);
%                dim = [.2 .5 .3 .3];
%                str = ['Pearson Corr = ' num2str(corr)];
%                annotation('textbox',dim,'String',str,'FitBoxToText','on');
               t = title(['Pearson Corr: ', num2str(corr)]);
               set(t,'position',get(t,'position')-[0 2.4 0])
%                legend('Original', 'Prediction');
               %subplot(num_rows, 3, ii1);
%               axis off
           end
           set(tlo.Children,'XTick',[], 'YTick', []); % all in one
           lg = legend('Original', 'Prediction');
           lg.Layout.Tile = 'East';
           f.Position = [1          41        1680         933];
           if num == 1
               sgtitle("X2's Prediction of X1")
               if ~strcmp(dest, '')
                   f.Renderer = 'painters';
                   saveas(f, [dest '2_1' '.svg']);
               end
           else
               sgtitle("X1's Prediction of X2")
               if ~strcmp(dest, '')
                   f.Renderer = 'painters';
                   saveas(f, [dest '1_2' '.svg']);
               end
           end
           
       end
       
       function visualize_lyapunovs_dim_by_dim(obj, num, dest)
           if strcmp(dest, '')
                dest = '';
           end
           if num == 1
               arr = obj.Xm;
               arr_size = size(arr);
               dim = arr_size(2);
               x_p = obj.x1_pred;
           else
               arr = obj.Ym;
               arr_size = size(arr);
               dim = arr_size(2);
               x_p = obj.x2_pred;
           end
           f = figure;
           num_rows = ceil(dim / 3);
           tlo = tiledlayout(num_rows,3,'TileSpacing','none','Padding','none');
           for ii1=1:dim
               p = x_p(:,ii1);
               nexttile, plot(arr(:, ii1));
               ylim([-8.5 8.5]);
               cur_dim = arr(:, ii1);
               rw_res = vratiotest(cur_dim);
               if rw_res == 1
                   random_walk_bool = 'false';
               else
                   random_walk_bool = 'true';
               end
               l_exp = lyapunovExponent(cur_dim);
%                nexttile, scatter(p, arr(:, ii1));
%                ylim([-5 5]);
               %err = immse(p, arr(:, ii1)); 
               corr = corrcoef(p, arr(:, ii1));
               corr = corr(1, 2);
%                dim = [.2 .5 .3 .3];
%                str = ['Pearson Corr = ' num2str(corr)];
%                annotation('textbox',dim,'String',str,'FitBoxToText','on');
               t = title(['Lyapunov: ', num2str(l_exp)]);
               set(t,'position',get(t,'position')-[0 2.3 0])
%                legend('Original', 'Prediction');
               %subplot(num_rows, 3, ii1);
%               axis off
           end
           set(tlo.Children,'XTick',[], 'YTick', []); % all in one
%            lg = legend('Original', 'Prediction');
%            lg.Layout.Tile = 'East'
           if num == 1
               sgtitle("X1")
           else
               sgtitle("X2")
           end
           f.Position = [1          41        1680         933];
           if ~strcmp(dest, '')
               f.Renderer = 'painters';
               saveas(f, [dest num2str(num) '.svg']);
           end
       end
       
       function visualize_prediction_corrs_dim_by_dim(obj, num)
           if num == 1
               arr = obj.Xm;
               arr_size = size(arr);
               dim = arr_size(2);
               x_p = obj.x1_pred;
           else
               arr = obj.Ym;
               arr_size = size(arr);
               dim = arr_size(2);
               x_p = obj.x2_pred;
           end
           figure;
           num_rows = ceil(dim / 4);
           %tlo = tiledlayout(num_rows,3,'TileSpacing','none','Padding','none');
           for ii1=1:dim
               p = x_p(:,ii1);
               %nexttile, 
               subplot(4, num_rows, ii1);
               scatter(p, arr(:, ii1)); axis square;
               xlim([-5 5]);
               ylim([-5 5]);
           end
           %set(tlo.Children,'XTick',[], 'YTick', []); % all in one
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
              all_linear_and_nonlinear_ccm_scores(obj, varargin)
          p = inputParser;

          use_dimension_weights = 'weigh_by_eigens';
          use_dimension_weights_default = false;
            addParameter(p, use_dimension_weights, ...
                       use_dimension_weights_default);

          parse(p, varargin{:});
          use_dimension_weights = p.Results.weigh_by_eigens;
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
           
          if use_dimension_weights
              weights_X = obj.Xm_latents;
              weights_Y = obj.Ym_latents;
              sc1_self_l_mean = wmean(sc1_self_l, weights_X);
              sc1_self_nl_mean = wmean(sc1_self_nl, weights_X);
              sc1_other_l_mean = wmean(sc1_other_l, weights_X);
              sc1_other_nl_mean = wmean(sc1_other_nl, weights_X);
              sc2_self_l_mean = wmean(sc2_self_l, weights_Y);
              sc2_self_nl_mean = wmean(sc2_self_nl, weights_Y);
              sc2_other_l_mean = wmean(sc2_other_l, weights_Y);
              sc2_other_nl_mean = wmean(sc2_other_nl, weights_Y);
          else
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
end
