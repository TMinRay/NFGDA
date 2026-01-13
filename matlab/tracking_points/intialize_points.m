function [pos, vel_pol] = intialize_points(radius, thetas, speed)
pos = [radius*sin(thetas); radius*cos(thetas)];
vel_pol = [speed*ones(size(thetas)); thetas];
end