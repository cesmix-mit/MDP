function [x,y,z] = sphere2cart(the,phi,r)

x = r.*sin(the).*cos(phi);
y = r.*sin(the).*sin(phi);
z = r.*cos(the);





