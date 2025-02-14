function [diff] = calc_small_angle_diff(dirn1, dirn2)
diff = dirn1 - dirn2;
if diff < -pi
    diff = diff + 2*pi;
elseif diff > pi
    diff = diff - 2*pi;
end
end
