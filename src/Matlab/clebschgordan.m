function C = clebschgordan(j1,m1,j2,m2,j,m,fac)
% CLEBSCHGORDAN computes Clebsch-Gordan coefficients.
%
% Inputs:
%   j1 - A scalar giving the first total angular momentum.
%   m1 - A scalar giving the projection of the first total angular
%        momentum.
%   j2 - A scalar giving the second total angular momentum.
%   m2 - A scalar giving the projection of the second total angular
%        momentum.
%   j  - A scalar giving the coupled total angular momentum.
%   m  - A scalar giving the projection of the coupled total angular
%        momentum.
%
% Outputs:
%   C - A scalar giving the required Clebsch-Gordan coefficient.

if m1+m2-m ~= 0 || j < abs(j1-j2) || j > j1+j2
    
    C = 0;
    return    
end

% factorial table
% fac = factable;

%---compute valid k values for summation---%
k = max([0,j2-j-m1,j1-j+m2]):min([j1+j2-j,j1-m1,j2+m2]);

%---compute coefficient---%
C = sqrt((2*j+1)*fac(1+j+j1-j2)*fac(1+j+j2-j1)*fac(1+j1+j2-j)/fac(1+j1+j2+j+1))*...
    sqrt(fac(1+j+m)*fac(1+j-m)*fac(1+j1+m1)*fac(1+j1-m1)*fac(1+j2+m2)*fac(1+j2-m2))*...
    sum(((-1).^k)./(fac(1+k).*fac(1+j1+j2-j-k).*fac(1+j1-m1-k).*fac(1+j2+m2-k).*fac(1+j-j2+m1+k).*fac(1+j-j1-m2+k)));

