function [forecast_history, gt_future, ts_future, forecast_clusters_history] = forecast(ts, anchor_t_index, track_history, gt_history, points_struct, clusters_history)

% curr_t = 2;
curr_t = ts(anchor_t_index);
ts_future = ts(ts > curr_t);
gt_future = gt_history(ts > curr_t);

% get indexed points for the current time
indexed_points = track_history{ts == curr_t};
indexed_clusters = clusters_history{ts == curr_t};
% remove_idx = false(size(indexed_points));
remove_idx = find([indexed_points.update_time]< curr_t | isnan([indexed_points.del_t]));
% remove_idx = find([indexed_points.update_time]< curr_t | isnan([indexed_points.del_t])...
%                     | ismember([indexed_points.cluster], find([]));
% for i = 1:length(indexed_points)
%     point = indexed_points{i};
%     if point.update_time < curr_t || isempty(point.displ) || isnan(point.del_t)
%         remove_idx(i) = true;
%     end
% end
indexed_points(remove_idx) = [];

forecast_history = {};
forecast_clusters_history = {};

last_t = curr_t;

for t_f = ts_future
    diff_t = t_f - curr_t;
    vel_cart = [indexed_clusters.vel_cart];
    point_pos = [indexed_points.pos];
    point_vel = vel_cart(:,[indexed_points.cluster]);
    displ = point_vel*diff_t;
    new_pos = point_pos + displ;
    next_points = indexed_points;

    for i = 1:length(indexed_points)
        next_points(i).prev = indexed_points(i);
        next_points(i).pos = new_pos(:,i);
        next_points(i).update_time = t_f;
        next_points(i).depth = next_points(i).depth + 1;

        next_points(i).displ = displ(:,i);
        next_points(i).dirn = atan2(displ(2,i),displ(1,i));
        next_points(i).del_dirn = calc_small_angle_diff(next_points(i).dirn, next_points(i).prev.dirn);
        next_points(i).del_t = next_points(i).update_time - next_points(i).prev.update_time;


%         point = indexed_points{i};
%         if point.update_time < last_t || isempty(point.displ) || isnan(point.del_t)
%             continue
%         end

%         % calculate future position
%         cluster_id = point.cluster;
%         vel_cart = indexed_clusters{cluster_id}.vel_cart;

% %         vel_cart = point.displ / point.del_t;
%         diff_t = t_f - last_t;
%         new_pos = point.pos + vel_cart*diff_t;

%         % create future point. 
%         next_point = points_struct;
%         next_point.pos = new_pos;
%         next_point.update_time = t_f;
%         next_point.prev = point;
%         next_point.depth = point.depth + 1;
%         next_point.displ = next_point.pos - point.pos;
%         next_point.dirn = atan2(next_point.displ(2),next_point.displ(1));
%         next_point.del_dirn = calc_small_angle_diff(next_point.dirn, point.dirn);
%         next_point.del_t = next_point.update_time - point.update_time;
%         % next_point.cluster = point.cluster;

%         % add to chain
%         % indexed_points{i} = next_point;
    end

    indexed_points = next_points;
    forecast_history{end+1} = indexed_points;
    forecast_clusters_history{end+1} = indexed_clusters;
    curr_t = t_f;
end

end


