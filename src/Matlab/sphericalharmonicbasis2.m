function [Ylm, Plm] = sphericalharmonicbasis2(x,y,z,L)
% (x,y,z): Cartesian coordinates
% L : the degree of spherical harmonics 
%   Ylmr    : real part of the spherical harmonics
%   Ylmi    : imag part of the spherical harmonics

x = x(:);
y = y(:);
z = z(:);
nx = length(x);

% Spherical coordinates
[the,phi,~] = cart2sphere(x,y,z);

% Spherical harmonics
Ylm = sphericalharmonics2(the,phi,L);

N = L+1;
Plm = cell(N,1);
% spherical shell harmonics for l > 0
for l=1:L   
    Plm{l+1} = zeros(nx,l);    
    for m = 0:l                                
        if m>0
            c = 1/((-1)^m);                            
            Plm{l+1}(:,m) =  c*conj(Ylm{l+1}(:,m+1));            
        end            
    end       
end

