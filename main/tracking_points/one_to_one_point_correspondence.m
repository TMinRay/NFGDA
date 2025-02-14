function [track_history, clusters_history] = one_to_one_point_correspondence(ts, det_history, cluster_points_history, track_history,max_age, max_del_dirn,points_struct, exp_gf_speed_range_kmph, clusters_history)

indexed_points = track_history{1};

for it = 2:length(ts)
    curr_time = ts(it);

    % calculate time difference
    curr_del_t = ts(it) - ts(it-1);
    min_dist = 0.95*exp_gf_speed_range_kmph(1)*curr_del_t;
    max_dist = 1.05*exp_gf_speed_range_kmph(2) * curr_del_t;
    
    curr_pos = det_history{it};
    curr_points_cluster = cluster_points_history{it};
    curr_clusters = clusters_history{it};

    % for one-to-one point correspondence
    next_points = cell(size(indexed_points));
    next_points_dist = Inf(size(indexed_points));
    new_tracks = {};

    for i_curr = 1:length(curr_pos)
        point = points_struct;
        point.pos = curr_pos(:, i_curr);
        point.cluster = curr_points_cluster(:, i_curr);
        point.update_time = curr_time;

        % correspond to previous points
        best_dist = inf;
        best_prev = points_struct;
        best_i_prev = -1;
        best_time = -inf;

        for i_prev = 1:length(indexed_points)
            prev_point = indexed_points{i_prev};

            if prev_point.update_time == curr_time
                continue
            end

            % distance relative to previous point's actual or projected
            % position
            displ = point.pos - prev_point.pos;
            dist = norm(displ);
            dirn = atan2(displ(2),displ(1));
            dist_range_cond = dist >= min_dist && dist < max_dist;

            if ~isempty(prev_point.vel_cart) && ~isnan(prev_point.update_time)
%                 vel_prev = prev_point.displ / prev_point.del_t;
                vel_prev = prev_point.vel_cart;
                del_t = point.update_time - prev_point.update_time;
                displ = point.pos - (prev_point.pos + vel_prev*del_t);
                dist = norm(displ);
                dist_range_cond = dist < max_dist;
            end

            dirn_diff = calc_small_angle_diff(dirn, prev_point.dirn);

            % conditions for corresponding points.
            eps = 0.1;
            dist_cond = dist < best_dist || (abs(dist - best_dist) <= eps && best_time < prev_point.update_time);

            if isnan(dirn_diff)
                is_match = dist_range_cond && dist_cond;
            else
                is_match = dist_range_cond && dist_cond && abs(dirn_diff) < max_del_dirn;
            end

            if is_match
                best_dist = dist;
                best_prev = prev_point;
                best_i_prev = i_prev;
                best_time = prev_point.update_time;
            end

        end

        % TODO(pjatau) Add description. Seems to be coasting.
        if best_i_prev < 0 || best_dist >= next_points_dist(best_i_prev)
            new_tracks{end+1} = point;
        else
            curr_point = next_points{best_i_prev};
            if ~isempty(curr_point)
                new_tracks{end+1} = curr_point;
            end
            next_points_dist(best_i_prev) = best_dist;
            next_points{best_i_prev} = point;
        end

    end

%     indexed_points = [indexed_points new_tracks];
    
    % remove old elements
    remove_idx = [];
    cum_cluster_vel = containers.Map('KeyType','char','ValueType','any');
    cum_valid_paths = containers.Map('KeyType','char','ValueType','any');
    cluster_vel_count = containers.Map('KeyType','char','ValueType','any');
    for j = 1:length(indexed_points)
        % TODO(pjatau) update index the remaining next points
        if j <= length(next_points) && ~isempty(next_points{j})
            point = next_points{j};
            point.prev = indexed_points{j};
            point.depth = indexed_points{j}.depth + 1;
            point.displ = point.pos - indexed_points{j}.pos;
            point.dirn = atan2(point.displ(2),point.displ(1));
            point.del_dirn = calc_small_angle_diff(point.dirn, indexed_points{j}.dirn);
            point.del_t = point.update_time - indexed_points{j}.update_time;
            point.vel_cart = point.displ / point.del_t;

            % Calculate velocity of each cluster
            % TODO(pjatau) weight by depth
            vel_cart = point.displ/point.del_t;
            cluster_id_num = point.cluster;
            cluster_id = num2str(cluster_id_num);
            valid_path = point.depth >= 2; % TODO(pjatau) use false alarm threshold
            if ~cum_cluster_vel.isKey(cluster_id)
                cum_cluster_vel(cluster_id) = vel_cart;
                cluster_vel_count(cluster_id) = 1;
                cum_valid_paths(cluster_id) = valid_path;
            else
                cum_cluster_vel(cluster_id) = cum_cluster_vel(cluster_id) + vel_cart;
                cluster_vel_count(cluster_id) = cluster_vel_count(cluster_id) + 1;
                cum_valid_paths(cluster_id) = cum_valid_paths(cluster_id) + valid_path;
            end

            indexed_points{j} = point;
        end

        if curr_time - indexed_points{j}.update_time > max_age
            remove_idx = [remove_idx j];
        end
    end
    indexed_points(remove_idx) = [];

    % Average velocity of clusters
%     avg_cluster_vel = containers.Map('KeyType','char','ValueType','any');
%     avg_cluster_depth = containers.Map('KeyType','char','ValueType','any');
    cluster_ids = unique(curr_points_cluster);

    for c_id = cluster_ids
        c_id_num = c_id;
        c_id = num2str(c_id_num);
        %         disp(cluster_vel_count(c_id));
        if cum_cluster_vel.isKey(c_id)
            curr_clusters{c_id_num}.vel_cart = cum_cluster_vel(c_id)/cluster_vel_count(c_id);
            curr_clusters{c_id_num}.num_anchors = cum_valid_paths(c_id);
%             avg_cluster_depth(c_id) = cum_valid_paths(c_id)/curr_clusters{c_id_num}.area;
            %         disp(c_id);
            %         disp(avg_cluster_vel(c_id));
        end

    end
    

%     % TODO(pjatau) em
%     fig = figure;
%     set(fig,'Position',[100 100 500 480]);
% 
% %     colorbar;
% %     caxis([0 5]);
% %     colormap jet; hold on;
% 
%     for i_key = 1:length(curr_clusters)
%         curr_key_num = curr_clusters{i_key}.id;
%         curr_key = num2str(curr_key_num);
%         mask = curr_points_cluster == curr_key_num;
% 
%         tmp_cluster_pos = curr_pos(:,mask);
%         tmp_pos = curr_clusters{curr_key_num}.centroid;
%         tmp_vel = curr_clusters{curr_key_num}.vel_cart;
% 
%         if ~isempty(tmp_vel)
% %             tmp_depth = avg_cluster_depth(curr_key) * ones(1,length(tmp_cluster_pos));
% 
%             scatter(tmp_cluster_pos(1,:),tmp_cluster_pos(2,:),80, 'filled', 'MarkerFaceAlpha',0.4); hold on;
%             quiver(tmp_pos(1), tmp_pos(2), tmp_vel(1), tmp_vel(2),"LineWidth",1.1); hold on;
%             text(tmp_pos(1), tmp_pos(2),sprintf("%4.2f", norm(tmp_vel)),FontSize=14);
%         end
%     end
%     hold off;
%     axis([-100,100,-100,100]);
%     grid on;
%     set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
%     print(['tmp/' sprintf('%03d', it) '.png'], '-dpng')


    % Assign cluster velocity to new tracks
    for i_new = 1:length(new_tracks)
        curr_point = new_tracks{i_new};
        assert(isempty(curr_point.vel_cart));
        c_id_num = curr_point.cluster;
        new_vel_cart = curr_clusters{c_id_num}.vel_cart;

        if ~isempty(new_vel_cart)
            new_tracks{i_new}.vel_cart = new_vel_cart;
        end
    end

    indexed_points = [indexed_points new_tracks];

    track_history{it} = indexed_points;
    clusters_history{it} = curr_clusters;
end
end
