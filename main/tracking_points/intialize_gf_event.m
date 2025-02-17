
% Simulated GF
if strcmp(event_id, 'sim_exp_len')

    max_t = 5;
    ts = 0:max_t;

    num_pts = 10;
    speed = 2;
    thetas = linspace(deg2rad(30), deg2rad(60), num_pts);
    [pos, vel] = intialize_points(2,thetas, speed);
    [det_history, gt_history] = initialize_tracks(pos, vel, ts);
    xi2 = []; yi2 = [];

    fig_axis = [0, 12, 0, 12];

elseif strcmp(event_id, 'sim_const_len')

    max_t = 5;
    ts = 0:max_t;

    num_pts = 10;
    speed = 2;
    x_pos = zeros(1, num_pts);
    y_pos = linspace(0,5, num_pts);
    pos = [x_pos; y_pos];
    vel = [speed*ones(1, num_pts); pi/2*ones(1, num_pts)];
    [det_history, gt_history] = initialize_tracks(pos, vel, ts);
    xi2 = []; yi2 = [];

    fig_axis = [-2, 12, -2, 7];

else % event_id is 'real_case'

    % Real GF
    nf_preds = {dir(fullfile(data_path,'nf_pred*')).name};
    nf_preds = sort(nf_preds);

    % ppi_names = {dir(v06_path).name};
    % ppi_names = ppi_names(3:end);
    % ppi_names = sort(ppi_names);

    det_history = {};
    gt_history = {};
    points_in_scan = {};

    cluster_points_history = {};
    clusters_history = {};
    tstamp0 = datetime(nf_preds{1}(16:23), 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC') + ...
            timeofday(datetime(nf_preds{1}(25:30), 'InputFormat', 'HHmmss'));
    ts = [];
    ppi_descs = {};
    alias_idx = -1;
    for i = 1:length(nf_preds)
        ppi_id = nf_preds{i};
        if ppi_id == "." || ppi_id == ".."
            continue
        end

        outfile_pos = strsplit(ppi_id,'.');
        ppi_id_no_ext = outfile_pos{1};
        ppi_name = ppi_id(12:end);
        tstamp = datetime(ppi_name(5:12), 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC') + ...
                timeofday(datetime(ppi_name(14:19), 'InputFormat', 'HHmmss'));
        radar_id = ppi_name(1:4);
        ppi_desc = [radar_id ', ' char(tstamp, 'MM/dd/yyyy, HH:mm:ss z')];
        ppi_descs{end+1} = ppi_desc;
        disp(ppi_desc);

        % time_hr = get_time_hour_UTC(ppi_name);
        % ts = [ts time_hr];
        ts = [ts hours(tstamp-tstamp0)];

        % % find beginning of 24 to 0 UTC aliasing
        % if length(ts) >= 2 && (ts(end) - ts(end-1) < -12)
        %     alias_idx = length(ts);
        % end

        load(fullfile(data_path,ppi_id));

        % cluster gf points
        radius = 4; % 8; % km
        [groups, areas, centroids] = cluster_points(xi2, yi2, nfout, radius,false);

        save(fullfile(clusters_path,ppi_id),'groups','areas','centroids','evalbox',...
            'xi2','yi2');
        
        mask = nfout;
        x = xi2(mask)';
        y = yi2(mask)';
        gp = groups(mask)';

        init_points = points_struct;
        init_points.update_time = ts(i);
        points_in_scan{i} = repmat(init_points, 1, length(x));
        for p = 1:length(x)
            points_in_scan{i}(p).pos = [x(p);y(p)];
            points_in_scan{i}(p).cluster = gp(p);
        end


        % det_history{i} = [x; y];
        gt_history{i} = evalbox;

        clusters_history{i} = repmat(cluster_struct, 1, length(areas));
        % indexed_clusters = {};
        for i_area = 1:length(areas)
            clusters_history{i}(i_area).area = areas(i_area);
            clusters_history{i}(i_area).centroid = centroids(:,i_area);
            clusters_history{i}(i_area).id = i_area;
            % curr_cluster = cluster_struct;
            % curr_cluster.area = areas(i_area);
            % curr_cluster.centroid = centroids(:,i_area);
            % curr_cluster.id = i_area;
            % indexed_clusters{i_area} = curr_cluster;
        end

        % cluster_points_history{i} = groups(mask)';
        % clusters_history{i} = indexed_clusters;
    end

    % % dealias time steps
    % if alias_idx >= 1
    %     ts(alias_idx:end) = ts(alias_idx:end) + 24;
    % end

    fig_axis = [-100 100 -100 100];
end
