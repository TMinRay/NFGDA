function L = arclength(x,y)
dX = gradient(x);
dY = gradient(y);
dL = hypot(dX,dY);
cum_L = cumtrapz(dL); % integrate arc segments
L = cum_L(end);
end