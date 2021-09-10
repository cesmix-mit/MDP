%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [x, R] = transformcoords(a,b,c,A,B,C,X)

% 
% P = [a(:) b(:) c(:)];
% V = a(1)*b(2)*c(3);
% L1 = cross(B,C);
% L2 = cross(C,A);
% L3 = cross(A,B);
% 
% L1 = (L1(:))';
% L2 = (L2(:))';
% L3 = (L3(:))';
% L = [L1; L2; L3];
% 
% R = P*L/V;
% 

% if norm(R'*R - eye(3))>1e-10 || norm(R*R' - eye(3))>1e-10     
%     error("Rotation matrix not orthogonal");
% end

% X -> Lambda by inv(A, B, C) 
% Lambda = lambda
% lambda -> x by (a, b, c)
R = [a(:) b(:) c(:)]*inv([A(:) B(:) C(:)]);

x = R*reshape(X, 3, []);




