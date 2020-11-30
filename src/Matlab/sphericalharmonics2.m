function Ylm = sphericalharmonics2(the,phi,L)
% (the,phi): spherical coordinates
% L : the degree of spherical harmonics 
%   Ylmr    : real part of the spherical harmonics
%   Ylmi    : imag part of the spherical harmonics

fac = factable;

the = the(:);
phi = phi(:);
np = length(phi);

% allocate memory
Ylm = cell(L+1,1);

% precompute
costhe = cos(the);

% l = 0
Ylm{1} = sqrt(1/(4*pi))*ones(np,1);

% l = 1
for l = 1:L
    n = 2*l+1;    
    m = 0:l;
    C = sqrt(n*fac(l - m + 1)./(4*pi*fac(l + m + 1)));

    P = legendre(l,costhe)';
    Ylm{l+1} = zeros([np l+1]);
    for m=0:l
        Ylm{l+1}(:,m+1) = C(m+1)*(P(:,m+1).*exp(sqrt(-1)*m*phi));
    end
end

