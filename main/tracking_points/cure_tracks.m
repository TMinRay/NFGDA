function [track_history, clusters_history] = cure_tracks(track_history, clusters_history, ts)

for it = 1:length(ts)
    curr_time = ts(it);
    indexed_points = track_history{it};
    indexed_clusters = clusters_history{it};

    if isempty(indexed_points)
        continue
    end

    remove_idx_points = false(size(indexed_points));
    remove_idx_clusters = false(size(indexed_clusters));

    for i_point = 1:length(indexed_points)
        point = indexed_points{i_point};
        if point.update_time ~= curr_time
            remove_idx_points(i_point) = true;
            continue
        end
        cluster_id = point.cluster;
        cluster_info = indexed_clusters{cluster_id};
        num_anchors = cluster_info.num_anchors;
        vel_cart = cluster_info.vel_cart;

        if num_anchors < 2 || (~isempty(vel_cart) && norm(vel_cart) <= 2)
            remove_idx_points(i_point) = true;
            remove_idx_clusters(cluster_id) = true;
        end
    end

    indexed_points(remove_idx_points) = [];
    track_history{it} = indexed_points;

    for i_cluster = 1:length(remove_idx_clusters)
        if remove_idx_clusters(i_cluster)
            indexed_clusters{i_cluster}.id = -1;
        end
    end
    clusters_history{it} = indexed_clusters;
end

end