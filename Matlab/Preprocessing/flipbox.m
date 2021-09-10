%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [a, b, c, flipxy, flipxz, flipyz] = flipbox(A, B, C, pbc)
% https://docs.lammps.org/Howto_triclinic.html
% To avoid extremely tilted boxes (which would be computationally inefficient), 
% LAMMPS normally requires that no tilt factor can skew the box more than half 
% the distance of the parallel box length, which is the first dimension in the 
% tilt factor (x for xz). This is required both when the simulation box is 
% created, e.g. via the create_box or read_data commands, as well as when the 
% box shape changes dynamically during a simulation, e.g. via the fix deform 
% or fix npt commands.

% For example, if xlo = 2 and xhi = 12, then the x box length is 10 and the 
% xy tilt factor must be between -5 and 5. Similarly, both xz and yz must be 
% between -(xhi-xlo)/2 and +(yhi-ylo)/2. Note that this is not a limitation, 
% since if the maximum tilt factor is 5 (as in this example), then configurations 
% with tilt = ?, -15, -5, 5, 15, 25, ? are geometrically all equivalent. If the 
% box tilt exceeds this limit during a dynamics run (e.g. via the fix deform 
% command), then the box is ?flipped? to an equivalent shape with a tilt factor 
% within the bounds, so the run can continue. 

% One exception to this rule is if the first dimension in the tilt factor 
% (x for xy) is non-periodic. In that case, the limits on the tilt factor are
% not enforced, since flipping the box in that dimension does not change the 
% atom positions due to non-periodicity. In this mode, if you tilt the system 
% to extreme angles, the simulation will simply become inefficient, due to the 
% highly skewed simulation box.
 
DELTAFLIP = 0.1;

xprd = A(1);
yprd = B(2);
%zprd = C(3);
xy = B(1);
xz = C(1);
yz = C(2);

xtiltmax = (0.5+DELTAFLIP)*xprd;
ytiltmax = (0.5+DELTAFLIP)*yprd;

flipxy = 0;
flipxz = 0;
flipyz = 0;

if (pbc(2)) 
    if (yz < -ytiltmax) 
        %yz = yz + yprd;
        xz = xz + xy;
        flipyz = 1;
    elseif (yz >= ytiltmax) 
        %yz = yz - yprd;
        xz = xz - xy;
        flipyz = -1;
    end
end

if pbc(1)
    if (xz < -xtiltmax) 
        %xz = xz + xprd;
        flipxz = 1;
    elseif (xz >= xtiltmax) 
        %xz = xz - xprd;
        flipxz = -1;
    end
    if (xy < -xtiltmax) 
        %xy = xy + xprd;
        flipxy = 1;
    elseif (xy >= xtiltmax) 
        %xy = xy - xprd;
        flipxy = -1;
    end
end
a = A;
b = B + flipxy*A;  
c = C + flipxz*A + flipyz*B;

