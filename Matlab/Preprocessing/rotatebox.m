%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [a, b, c] = rotatebox(A, B, C)
% https://docs.lammps.org/Howto_triclinic.html
% LAMMPS also allows simulations to be performed in triclinic (non-orthogonal) 
% simulation boxes shaped as a parallelepiped with triclinic symmetry. 
% The parallelepiped has its ?origin? at (xlo,ylo,zlo) and is defined by 3 edge 
% vectors starting from the origin given by a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); 
% c = (xz,yz,zhi-zlo). xy,xz,yz can be 0.0 or positive or negative values and are 
% called ?tilt factors? because they are the amount of displacement applied to faces 
% of an originally orthogonal box to transform it into the parallelepiped. 
% In LAMMPS the triclinic simulation box edge vectors a, b, and c cannot be 
% arbitrary vectors. As indicated, a must lie on the positive x axis. b must 
% lie in the xy plane, with strictly positive y component. c may have any orientation 
% with strictly positive z component. The requirement that a, b, and c have strictly 
% positive x, y, and z components, respectively, ensures that a, b, and c form a 
% complete right-handed basis. These restrictions impose no loss of generality, since 
% it is possible to rotate/invert any set of 3 crystal basis vectors so that they 
% conform to the restrictions.

% For example, assume that the 3 vectors A,B,C are the edge vectors of a general 
% parallelepiped, where there is no restriction on A,B,C other than they form a 
% complete right-handed basis i.e. A x B . C > 0. The equivalent LAMMPS a,b,c are 
% a linear rotation of A, B, and C and can be computed as follows:

Anorm = norm(A);
Bnorm = norm(B);
Cnorm = norm(C);
Ahat = A/Anorm;

ax = Anorm;
bx = dot(B,Ahat);
by = sqrt(Bnorm^2 - bx^2); %norm(cross(Ahat,B));
cx = dot(C,Ahat);
cy = (dot(B, C) - bx*cx)/by;
cz = sqrt(Cnorm^2 - cx^2 - cy^2);

a = [ax 0 0];
b = [bx by 0];
c = [cx cy cz];



