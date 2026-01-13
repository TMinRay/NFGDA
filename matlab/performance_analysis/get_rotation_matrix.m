function mat = get_rotation_matrix(rot_angle_rad)
mat = [cos(rot_angle_rad) -sin(rot_angle_rad); sin(rot_angle_rad) cos(rot_angle_rad)];
end