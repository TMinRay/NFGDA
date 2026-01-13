function [p, best_rot, best_x] = polyfit_rotation(x, y, n)
% find rotation that maximizes the range of x
points = [x; y];
[best_rot, new_points] = binary_search_rotation(0, pi, points, 5);

% polynomial fit at best rotation
if best_rot == 0
    p = polyfit(x,y,n);
else
    p = polyfit(new_points(1,:), new_points(2,:), n);
    best_x = new_points(1,:);
end
end
