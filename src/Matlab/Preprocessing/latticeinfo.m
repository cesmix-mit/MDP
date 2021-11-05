%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [lx, ly, lz, xy, xz, yz, la, lb, lc, cosa, cosb, cosc, V] = latticeinfo(a,b,c)
% (a, b, c): lattice vectors of a parallelepiped with triclinic symmetry. 

if (a(2)>0) || (a(3)>0) || (b(3)>0)
    error("Lattice vectors are invalid");
end

lx = a(1);
ly = b(2);
lz = c(3);
xy = b(1);
xz = c(1);
yz = c(2);

la = norm(a);
lb = norm(b);
lc = norm(c);
cosa = (xy*xz + ly*yz)/(lb*lc);
cosb = xz/lc;
cosc = xy/lb;

V = lx*ly*lz;





