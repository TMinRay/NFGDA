clear; clc;
close all;
tic
rng(101)

% I/O options
fig_dir = './figures';
nf_pred_dir = './nf_preds';
v06_dir = '../../V06';
clusters_dir = './clusters';
root_forecast_dir = './forecasts';

% Initialize GF event
event_id = 'real_case'; %'sim_exp_len' 'sim_const_len' 'real_case'
% case_names = {'KABX20210708_00'};
case_names = {'KABX20200705_21'};

rois = containers.Map('KeyType','char','ValueType','any');
% rois('KGGW20210722_23') = {[93.93, 113.93, 113.93, 93.93],[220, 220, 230, 230]}; %true 64.577 mi (103.93 km), 225.43 degs

asos = containers.Map('KeyType','char','ValueType','any');
asos('KABX20210708_00') = {["AEG","ABQ","BRG"], ...
    [2.632095137220614, -0.5085441810964234;...
    18.967467286260703, -12.011131240458427;...
    -1.0801345232802038, -56.01381877018714]};

% Tracking options
points_struct = struct('pos', [], 'displ',[], 'dirn', NaN, 'del_t', NaN, ...
    'depth', 0, 'id','', 'prev', NaN, 'update_time', NaN, 'del_dirn', NaN,...
    'cluster', NaN, 'vel_cart', []);
cluster_struct = struct('id', NaN, 'centroid',[],'area',NaN, 'num_anchors', 0, ...
    'vel_cart',[NaN;NaN]);
exp_gf_speed_range = [4, 32]; % Hwang thesis: GF propagation speed is 5-20 mps (18-72 kmph)
max_del_dirn = pi/16; %pi/4;

% Forecast options
anchor_t_index = 4;% 12;

for case_id = 1:length(case_names)

    case_name = case_names{case_id};
    data_path = fullfile(nf_pred_dir, case_name);
    v06_path = fullfile(v06_dir, case_name);

    clusters_path = fullfile(clusters_dir, case_name);
    if ~isfolder(clusters_path)
        mkdir(clusters_path);
    end

    cured_clusters_path = fullfile(clusters_dir, [case_name '_cured_gc']);
    if ~isfolder(cured_clusters_path)
        mkdir(cured_clusters_path);
    end

    if rois.isKey(case_name)
        roi = rois(case_name);
        R_roi = roi{1};
        azdeg_roi = roi{2};
    else
        R_roi = [];
        azdeg_roi = [];
    end

    if asos.isKey(case_name)
        stations_info = asos(case_name);
        station_ids = stations_info{1};
        station_locs = stations_info{2};
    else
        station_ids = [];
        station_locs = [];
    end

    % Initialize gust front events
    intialize_gf_event

    % TODO(pjatau) move outside loop
    exp_del_t = max(ts(2:end) - ts(1:end-1));
    max_age = 1.05* exp_del_t; % 1*1; 2;

    % Initialize tracks
    % track_history = {};
    % indexed_points = cell(1, length(det_history{1}));
    % for i = 1:length(det_history{1})
    %     point = points_struct;
    %     point.pos = det_history{1}(:,i);
    %     point.cluster = cluster_points_history{1}(i);
    %     point.update_time = ts(1);
    %     indexed_points{i} = point;
    % end

    % indexed_points = repmat(points_struct, 1, size(det_history{1},2));
    % point.update_time = ts(1);
    % for i = 1:size(det_history{1},2)
    %     indexed_points(i).pos = det_history{1}(:,i);
    %     indexed_points(i).cluster = cluster_points_history{1}(i);
    % end

    % track_history{1} = indexed_points;

    % Point tracking
    % [track_history, clusters_history] = one_to_one_point_correspondence(ts, det_history,cluster_points_history, track_history, max_age, max_del_dirn, points_struct, exp_gf_speed_range, clusters_history);
    [track_history, clusters_history] = one_to_one_point_correspondence(ts, points_in_scan, max_age, max_del_dirn, points_struct, exp_gf_speed_range, clusters_history);

    % Figures for initial detections
    title_suffix = sprintf("max. angle deviation: %4.1f degrees.", rad2deg(max_del_dirn));
    fig_dir_tracks = fullfile(fig_dir,'hit_miss', case_name);
    plot_tracks(track_history, gt_history, xi2, yi2, ts, 0, fig_axis, false, ...
        true, fig_dir_tracks,'hit-miss',title_suffix, clusters_history,...
        R_roi, azdeg_roi, ppi_descs, station_ids, station_locs);

    % Figures for point tracks
    title_suffix = sprintf("max. angle deviation: %4.1f degrees.", rad2deg(max_del_dirn));
    fig_dir_tracks = fullfile(fig_dir,'tracks', case_name);
    plot_tracks(track_history, gt_history, xi2, yi2, ts, 0, fig_axis, true, ...
        true, fig_dir_tracks,'tracks',title_suffix, clusters_history, R_roi, ...
        azdeg_roi, ppi_descs, station_ids, station_locs);

    close all;

% % False alarm mitigation
% [track_history, clusters_history] = cure_tracks(track_history, clusters_history, ts);
% 
% % Figures for detections after false alarm mitigation
% title_suffix = sprintf("max. angle deviation: %4.1f degrees.", rad2deg(max_del_dirn));
% fig_dir_tracks = fullfile(fig_dir,'hit_miss_cured', case_name);
% plot_tracks(track_history, gt_history, xi2, yi2, ts, 0, fig_axis, false, ...
%     true, fig_dir_tracks,'hit-miss',title_suffix, clusters_history, ...
%     R_roi, azdeg_roi, ppi_descs, station_ids, station_locs);
% 
% close all;

%     % Export detections for PLD, Precision analysis
%     lb_x = -100; lb_y = -100;
%     del_x = 0.5; del_y = 0.5;
% 
%     for it = 1:length(ts)
% 
%         tmp_file = fullfile(clusters_path,[sprintf('%02d', it+1) '.mat']);
%         load(tmp_file);
% 
%         figure
%         subplot(1,2,1);
%         pcolor(groups); shading flat;
% 
%         curr_time = ts(it);
%         indexed_points = track_history{it};
%         indexed_clusters = clusters_history{it};
%         evalbox = gt_history{it};
% 
%         groups = map_points_to_2d_grid(indexed_points,lb_x, lb_y, del_x, del_y);
% 
%         if isempty(indexed_clusters)
%             areas = [];
%         else
%             areas = [indexed_clusters{:}];
%             areas = {areas.area};
%         end
% 
%         subplot(1,2,2);
%         pcolor(groups); shading flat;
% 
%         cured_clusters_file = fullfile(cured_clusters_path,[sprintf('%02d', it+1) '.mat']);
%         save(cured_clusters_file,'groups','areas','evalbox','xi2','yi2');
%     end

    % Forecasting
    close all;

    [tracks_future, gt_future, ts_future, f_clusters_history] = forecast(ts, anchor_t_index, track_history, gt_history, points_struct, clusters_history);

    points_forecast = tracks_future{1};
    if isempty(points_forecast)
        continue
    end

    forecast_id = sprintf('%s_anchor_%d', case_name, anchor_t_index+1);
    % forecast_dir = fullfile(root_forecast_dir, forecast_id);
    % if ~isfolder(forecast_dir)
    %     mkdir(forecast_dir);
    % end

    % plot forecast
    fig_dir_tracks = fullfile(fig_dir,'forecast_results',forecast_id);
    title_suffix = sprintf("max. angle deviation: %4.1f degrees.", rad2deg(max_del_dirn));
    ppi_descs_future = ppi_descs(end-length(ts_future)+1:end);

    plot_tracks(tracks_future, gt_future, xi2, yi2, ts_future, 0, fig_axis,...
        false, true, fig_dir_tracks,'hit-miss',title_suffix, f_clusters_history,...
        R_roi, azdeg_roi,ppi_descs_future, station_ids, station_locs);

    close all;
end
toc