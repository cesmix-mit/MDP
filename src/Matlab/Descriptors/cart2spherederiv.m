function [the,phi,r,thex,they,thez,phix,phiy,phiz,rx,ry,rz] = cart2spherederiv(x,y,z)

r = sqrt(x.^2 + y.^2 + z.^2);
the = acos(z./r);
phi = atan2(y,x);

r2 = r.*r;
rxy = sqrt(x.^2 + y.^2);
rxy2 = rxy.*rxy;
rr2 = rxy.*r2;

rx = x./r;
ry = y./r;
rz = z./r;
thex = x.*z./rr2;
they = y.*z./rr2;
thez = -rxy./r2;
phix = -y./rxy2;
phiy = x./rxy2;
phiz = 0*z;

% n = 1000;
% dx = 1e-8;
% the = pi*rand(n,1);
% phi = -pi + 2*pi*rand(n,1);
% r = 0.01 + 0.99*rand(n,1);
% [x,y,z] = sphere2cart(the,phi,r);
% 
% [the,phi,r,thex,they,thez,phix,phiy,phiz,rx,ry,rz] = cart2spherederiv(x,y,z);
% 
% e = zeros(n,9);
% [the2,phi2,r2] = cart2sphere(x+dx,y,z);
% [the1,phi1,r1] = cart2sphere(x-dx,y,z);
% e(:,1) = (the2-the1)/(2*dx) - thex;
% e(:,2) = (phi2-phi1)/(2*dx) - phix;
% e(:,3) = (r2-r1)/(2*dx) - rx;
% 
% [the2,phi2,r2] = cart2sphere(x,y+dx,z);
% [the1,phi1,r1] = cart2sphere(x,y-dx,z);
% e(:,4) = (the2-the1)/(2*dx) - they;
% e(:,5) = (phi2-phi1)/(2*dx) - phiy;
% e(:,6) = (r2-r1)/(2*dx) - ry;
% 
% [the2,phi2,r2] = cart2sphere(x,y,z+dx);
% [the1,phi1,r1] = cart2sphere(x,y,z-dx);
% e(:,7) = (the2-the1)/(2*dx) - thez;
% e(:,8) = (phi2-phi1)/(2*dx) - phiz;
% e(:,9) = (r2-r1)/(2*dx) - rz;


