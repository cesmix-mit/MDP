function [v,w] = boundingbox(r, pbc, a, b, c)
% pbc periodic boundary conditions
% a = (a1, a2, a3), b = (b1,b2,b3), c=(c1,c2,c3) are three vectors defining a parallelepiped 
% r is the distance from the original parallelepiped to the bounding parallelepiped 
% v are 8 vertices defining the original parallelepiped 
% w are 8 vertices defining the bounding parallelepiped 

if nargin<5
    [v, w] = boundingparallelogram(pbc, a, b, r);
else    
    [v, w] = boundingparallelepiped(pbc, a, b, c, r);
end


