function [forecast_history, gt_future, ts_future, forecast_clusters_history] = forecast(ts, anchor_t_index, track_history, gt_history, points_struct, clusters_history)

% curr_t = 2;
curr_t = ts(anchor_t_index);
ts_future = ts(ts > curr_t);
gt_future = gt_history(ts > curr_t);

% get indexed points for the current time
indexed_points = track_history{ts == curr_t};
remove_idx = false(size(indexed_points));
for i = 1:length(indexed_points)
    point = indexed_points{i};
    if point.update_time < curr_t || isempty(point.displ) || isnan(point.del_t)
        remove_idx(i) = true;
    end
end
indexed_points(remove_idx) = [];

indexed_clusters = clusters_history{ts == curr_t};
forecast_history = {};
forecast_clusters_history = {};

% TODO(pjatau) em
% figure
% tmp_pos = [];
% for i = 1:length(indexed_points)
%     point = indexed_points{i};
%     if point.update_time < curr_t
%         continue
%     end
%     tmp_pos = [tmp_pos point.pos];
% end
% 
% scatter(tmp_pos(1,:), tmp_pos(2,:)); 
% axis([-90 90 -90 90]);


last_t = curr_t;
for t_f = ts_future
    for i = 1:length(indexed_points)
        point = indexed_points{i};
        if point.update_time < last_t || isempty(point.displ) || isnan(point.del_t)
            continue
        end

        % calculate future position
        cluster_id = point.cluster;
        vel_cart = indexed_clusters{cluster_id}.vel_cart;

%         vel_cart = point.displ / point.del_t;
        diff_t = t_f - last_t;
        new_pos = point.pos + vel_cart*diff_t;

        % create future point. 
        next_point = points_struct;
        next_point.pos = new_pos;
        next_point.update_time = t_f;
        next_point.prev = point;
        next_point.depth = point.depth + 1;
        next_point.displ = next_point.pos - point.pos;
        next_point.dirn = atan2(next_point.displ(2),next_point.displ(1));
        next_point.del_dirn = calc_small_angle_diff(next_point.dirn, point.dirn);
        next_point.del_t = next_point.update_time - point.update_time;
        next_point.cluster = point.cluster;

        % add to chain
        indexed_points{i} = next_point;
    end
    forecast_history{end+1} = indexed_points;
    forecast_clusters_history{end+1} = indexed_clusters;
    last_t = t_f;
end

end


