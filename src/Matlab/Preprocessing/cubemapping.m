%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [B2C, C2B] = cubemapping(a, b, c)
% a, b, c are three vectors describing the simulation box
% a = [a1, a2, a3]
% b = [b1, b2, b3]
% c = [c1, c2, c3]

if nargin<3
    C2B = [a(:) b(:)]; % transform the unit square to the simulation box (a, b)  
    B2C = inv(C2B); % transform the simulation box (a, b) to the unit square 
else
    C2B = [a(:) b(:) c(:)]; % transform the unit cube to the simulation box (a, b, c)  
    B2C = inv(C2B); % transform the simulation box (a, b, c) to the unit cube
end

return;

end

