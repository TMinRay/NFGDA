function [groups] = map_points_to_2d_grid(indexed_points, lb_x, lb_y, del_x, del_y)

groups = zeros(abs(lb_x)*2/del_x + 1, abs(lb_y)*2/del_y + 1);

if isempty(indexed_points)
    disp("There are no points");
end

for i_point = 1:length(indexed_points)
    point = indexed_points{i_point};

    curr_x = point.pos(1);
    curr_y = point.pos(2);
    ind_c = int32(1 + (curr_x - lb_x)/del_x);
    ind_r = int32(1+ (curr_y - lb_y)/del_y);
    groups(ind_r, ind_c) = point.cluster;
end

end