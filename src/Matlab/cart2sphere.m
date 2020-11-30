function [the,phi,r] = cart2sphere(x,y,z)

r = sqrt(x.^2 + y.^2 + z.^2);
the = acos(z./r);
phi = atan2(y,x);
% ind = (phi<0);
% phi(ind) = phi(ind) + 2*pi;
