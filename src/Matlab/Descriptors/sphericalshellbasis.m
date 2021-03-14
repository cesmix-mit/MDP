function [Sr, Si] = sphericalshellbasis(x,y,z,K,L)
% (x,y,z): Cartesian coordinates
% K : the number of radial functions
% L : the degree of spherical harmonics 
%   Ylmr    : real part of the spherical shell harmonics
%   Ylmi    : imag part of the spherical shell harmonics

x = x(:);
y = y(:);
z = z(:);
nx = length(x);

% Spherical coordinates
[the,phi,r] = cart2sphere(x,y,z);

% Spherical harmonics
[Ylmr, Ylmi] = sphericalharmonics(the,phi,L);

% roots of spherical Bessel functions
N = L+1;
kind = 2;
n = (1:N)' - 1/2;
x0 = besselzero(n, K, kind);

% Spherical Bessel functions
g = zeros(nx,K,N);
for n=1:N
    for k = 1:K
        r0 = x0(n,k)*r;
        g(:,k,n) = -sqrt(pi./(2*r0)).*bessely(n-1/2,r0);
    end
end

Sr = cell(N,1);
Si = cell(N,1);

% spherical shell harmonics for l = 0
Sr{1} = zeros(nx,K);
Si{1} = zeros(nx,K);    
for k = 1:K
    Sr{1}(:,k) = sqrt(1/(4*pi))*g(:,k,1);
    Si{1}(:,k) = 0;
end       

% spherical shell harmonics for l > 0
for l=1:L   
    Sr{l+1} = zeros(nx,K,l+1);
    Si{l+1} = zeros(nx,K,l+1);    
    for m = 1:(l+1)    
        for k = 1:K                
            Sr{l+1}(:,k,m) = g(:,k,l+1).*Ylmr{l+1}(:,m);
            Si{l+1}(:,k,m) = g(:,k,l+1).*Ylmi{l+1}(:,m);            
        end
    end       
end

