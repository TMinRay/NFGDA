function [mid, new_points] =  binary_search_rotation(lo, hi, points, max_iter)
eps = pi/128;

t = 1;
while lo < hi && t <= max_iter
    mid = lo + (hi - lo)/2;
    t = t + 1;

    rot_mat = get_rotation_matrix(mid);
    new_points = rot_mat * points;
    f_mid = range(new_points(1,:));

    rot_mat_ahead = get_rotation_matrix(mid + eps);
    new_points_ahead = rot_mat_ahead * points;
    f_ahead = range(new_points_ahead(1,:));

    if f_mid < f_ahead
        lo = mid;
    elseif f_mid > f_ahead
        hi = mid;
    else
        error("Need a larger delta x for differentiation.");
    end

end

end
