function points =  polyval_rotation(p,x, angle_rot)
x = sort(x);
y_pred = polyval(p, x);
points = [x; y_pred];
if angle_rot ~= 0
    rot_mat = get_rotation_matrix(-angle_rot);
    points = rot_mat * points;
end
end
