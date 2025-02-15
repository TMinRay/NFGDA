function [] = plot_tracks(track_history, gt_history, x_grid, y_grid, ts, ...
    min_depth, fig_axis, show_tracks, show_most_recent, fig_dir, cmap_type,...
    title_suffix, clusters_history, R_roi, azdeg_roi, ppi_descs, ...
    station_ids, station_locs)

def_color_order = get(gca,'colororder'); close;
len_col_order = length(def_color_order);
c_max_depth = 10;

marker = "o";
marker_size = 80;
marker_alpha = 0.4;

% region of interest
if ~isempty(R_roi)
    x_roi = R_roi .* sind(azdeg_roi);
    y_roi = R_roi .* cosd(azdeg_roi);
    x_roi(end + 1) = x_roi(1);
    y_roi(end + 1) = y_roi(1);
    roi = inpolygon(x_grid, y_grid,x_roi,y_roi);
end

% title_str = sprintf('Min. path depth: %02d.', min_depth);
% title is set below


if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

for it = 1:length(ts)
    indexed_points = track_history{it};
    indexed_clusters = clusters_history{it};
    gt_region = gt_history{it};

    CC = bwconncomp(gt_region);
    num_boxes = CC.NumObjects;

    title_str = ppi_descs{it};

    fig = figure;
    set(fig,'Position',[100 100 500 480]);
    axis(fig_axis);
    grid on;
    title(title_str);
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);

    % add colormap in case there are no tracks
    if strcmp(cmap_type, 'depth')
        colorbar;
        caxis([0 c_max_depth]);
        colormap jet; hold on;
    end

    plot_pos = {};
    most_rec_pos = [];
    plot_depth = {};
    most_rec_depth = [];
    most_rec_cluster = [];

    for i = 1:length(indexed_points)
        point = indexed_points(i);
        depth_cond = point.depth < min_depth;
%         depth_cond = point.num_anchors < min_depth;
        if point.update_time < ts(it)
            continue
        end

        cluster_id = point.cluster;
        cluster_info = indexed_clusters(cluster_id);
        depth_cond = cluster_info.num_anchors < 2;

%         if depth_cond
%             continue
%         end

        track_pos = [];
        track_depth = [];
        most_rec_pos = [most_rec_pos point.pos];
        most_rec_depth = [most_rec_depth point.depth];
        most_rec_cluster = [most_rec_cluster point.cluster];

        %         track_vel = [];
        %         track_dirn = {};
        %         track_del_dirn = {};

        while isfield(point, 'pos')
            if point.depth < -1
                break
            end
            track_pos = [track_pos point.pos];
            track_depth = [track_depth point.depth];
            %             curr_vel = point.displ./point.del_t;
            %             track_vel = [track_vel curr_vel];
            %             track_dirn{end + 1} = point.dirn * 180/pi;
            %             track_del_dirn{end+1} = point.del_dirn * 180/pi;
            point = point.prev;
        end

        plot_pos{end+1} = track_pos;
        plot_depth{end+1} = track_depth;
    end

    % figure

    if show_most_recent

        if ~isempty(most_rec_pos)
            if strcmp(cmap_type, 'depth')
                c = most_rec_depth;
                c(c > c_max_depth) = 10;
            elseif strcmp(cmap_type,'clusters')
                c = zeros(length(most_rec_cluster),3);
                c_idx = 1 + mod(most_rec_cluster-1,len_col_order);
                c = def_color_order(c_idx,:);
            elseif strcmp(cmap_type, 'hit-miss')
                marker = 'o';
                marker_size = 8;
                marker_alpha = 0.8;

                c = zeros(length(most_rec_pos),3);
                c(:,1) = 1;     % initialize all misses
                for i_box = 1:num_boxes
                    idx_box = CC.PixelIdxList{i_box};
                    x_box = x_grid(idx_box)';
                    y_box = y_grid(idx_box)';

                    for idx_c = 1:length(most_rec_depth)
                        for kk = 1:length(x_box)
                            diff1 = most_rec_pos(1,idx_c) - x_box(kk);
                            diff2 = most_rec_pos(2,idx_c) - y_box(kk);
                            if abs(diff1) < 0.5 && abs(diff2) < 0.5
                                c(idx_c,1) = 0;
                            end
                        end
                    end
                end
            else
                marker_alpha = 0.1;
                c = zeros(length(most_rec_pos),3); % 1->7, 8 -> 14
                for jj = 1:length(most_rec_pos)
                    c(jj,:) = def_color_order(1 + mod(jj-1,len_col_order),:);
                end
            end

            scatter(most_rec_pos(1,:), most_rec_pos(2,:),marker_size, c, 'filled','Marker', marker, 'MarkerFaceAlpha',marker_alpha); hold on;
        end
        
    else
        for ii = 1:length(plot_pos)
            track_pos = plot_pos{ii};
            track_depth = plot_depth{ii};

            c = track_depth;
            c(c > c_max_depth) = 10;

            scatter(track_pos(1,:), track_pos(2,:),marker_size, c,'filled','Marker', marker,'MarkerFaceAlpha',marker_alpha); hold on;
        end
    end

    if show_tracks
        for ii = 1:length(plot_pos)
            track_pos = plot_pos{ii};
            curr_color = def_color_order(1 + mod(ii-1,len_col_order),:);
            plot(track_pos(1,:), track_pos(2,:), 'Color',curr_color); hold on;
        end
    end

    % ground truth region
    if ~isempty(gt_region)
        contour(x_grid,y_grid,gt_region,'y-','linewidth',1); hold on;
    end

    % plot region of interest
    if ~isempty(R_roi)
        contour(x_grid, y_grid, double(roi),'Color','m','Linewidth',1); hold on;
    end

    if ~isempty(station_ids)
        for i = 1:length(station_ids)
            station_id = station_ids(i);
            station_loc = station_locs(i,:);
            scatter(station_loc(1), station_loc(2), 20,'green', 'filled','Marker','d', 'MarkerEdgeColor','k');
            text(station_loc(1), station_loc(2),station_id);
        end
    end

    hold off;

    axis(fig_axis);
    grid on;
    title(title_str, 'FontSize',13);
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], ...
        'LineWidth', 2, 'XTick', -100:20:100,'YTick', -100:20:100);
    outfile_track = fullfile(fig_dir, ['track_' sprintf('%03d', it) '.png']);
    print(outfile_track, '-dpng');
    %     close all;
end
end
