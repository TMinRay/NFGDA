function [track_history] = cure_tracks(track_history, clusters_history, ts)

for it = 1:length(ts)
    curr_time = ts(it);
    indexed_points = track_history{it};
    indexed_clusters = clusters_history{it};
    bs = length(indexed_points);
    if isempty(indexed_points)
        continue
    end
    [indexed_points.update_time]~= curr_time;
    remove_idx = find([indexed_points.update_time] ~= curr_time | isnan([indexed_points.del_t]));
    indexed_points(remove_idx) = [];

    vel_cart = [indexed_clusters.vel_cart];
    point_vel = vel_cart(:,[indexed_points.cluster]);
    point_norm = sqrt(sum(point_vel.^2,1));
    num_anchors = [indexed_clusters.num_anchors];
    point_num_anchors = num_anchors([indexed_points.cluster]);

    remove_idx_points = find( point_num_anchors<2 | point_norm<2 | isnan(point_norm) );

    % remove_idx_points = false(size(indexed_points));
    % remove_idx_clusters = false(size(indexed_clusters));

    % for i_point = 1:length(indexed_points)
    %     point = indexed_points{i_point};
    %     if point.update_time ~= curr_time
    %         remove_idx_points(i_point) = true;
    %         continue
    %     end
    %     cluster_id = point.cluster;
    %     cluster_info = indexed_clusters{cluster_id};
    %     num_anchors = cluster_info.num_anchors;
    %     vel_cart = cluster_info.vel_cart;

    %     if num_anchors < 2 || (~isempty(vel_cart) && norm(vel_cart) <= 2)
    %         remove_idx_points(i_point) = true;
    %         remove_idx_clusters(cluster_id) = true;
    %     end
    % end

    indexed_points(remove_idx_points) = [];
    track_history{it} = indexed_points;
    fprintf('Cure tracks frame %d :: %d -> %.3f \n', it, bs, length(indexed_points));

    % for i_cluster = 1:length(remove_idx_clusters)
    %     if remove_idx_clusters(i_cluster)
    %         indexed_clusters{i_cluster}.id = -1;
    %     end
    % end
    % clusters_history{it} = indexed_clusters;
end

end