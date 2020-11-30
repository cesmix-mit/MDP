function [Ylmr, Ylmi] = sphericalharmonics(the,phi,L)
% (the,phi): spherical coordinates
% L : the degree of spherical harmonics 
%   Ylmr    : real part of the spherical harmonics
%   Ylmi    : imag part of the spherical harmonics

fac = factable;

ns = size(the);
the = the(:);
phi = phi(:);
np = length(phi);

% temporary arrays to compute the recurrence 
tmp = zeros([np L]);
tmp(:,1) = 1;

% allocate memory
P = zeros([np (L+1)]);
Ylmr = cell(L+1,1);
Ylmi = cell(L+1,1);

% l = 0
Ylmr{1} = sqrt(1/(4*pi))*ones(np,1);
Ylmi{1} = zeros(np,1);

% precompute
costhe = cos(the);
a = -sin(the); %-sqrt(1-costhe.^2);

% l = 1
l = 1;
n = 2*l+1;    
m = 0:l;
C = sqrt(n*fac(l - m + 1)./(4*pi*fac(l + m + 1))); 
P(:,l+1) = a;
P(:,l) = costhe;    

Ylmr{2} = zeros([np 2]);
Ylmi{2} = zeros([np 2]);    
Ylmr{2}(:,1) = (C(1)*cos(0))*P(:,1);    % m=0
Ylmr{2}(:,2) = C(2)*(cos(phi).*P(:,2)); % m=1
Ylmi{2}(:,1) = (C(1)*sin(0))*P(:,1);    % m=0
Ylmi{2}(:,2) = C(2)*(sin(phi).*P(:,2)); % m=1
if ns(1)>1 && ns(2)>1
    Ylmr{2} = reshape(Ylmr{2},[ns (l+1)]);
    Ylmi{2} = reshape(Ylmi{2},[ns (l+1)]);
end    

for l=2:L
    Ylmr{l+1} = zeros([np l+1]);
    Ylmi{l+1} = zeros([np l+1]);    
    
    n = 2*l+1;    
    m = 0:l;
    C = sqrt(n*fac(l - m + 1)./(4*pi*fac(l + m + 1))); 
    
    Pll = P(:,l);
    tmp(:,l) = Pll;
    P(:,l+1) = (2*(l-1)+1)*(a.*Pll); % P_{l+1}^{l+1} = -(2*l+1)*sqrt(1-x.^2)*P_l^l
    P(:,l) = (2*(l-1)+1)*(costhe.*Pll);   % P_{l+1}^{l} = (2*l+1)*x*P_l^l
    for m = 1:(l-1)
        % (l-m+1)*P_{l+1}^m = (2*l+1)*(x.*P_l^m - (l+m)*P_{l-1}^m
        Pll = P(:,m);
        P(:,m) = ((2*(l-1)+1)*(costhe.*Pll) - (l+m-2)*tmp(:,m))/(l-m+1);
        tmp(:,m) = Pll;
    end        
        
    for m=1:l+1
        Ylmr{l+1}(:,m) = C(m)*(cos((m-1)*phi).*P(:,m)); % real
        Ylmi{l+1}(:,m) = C(m)*(sin((m-1)*phi).*P(:,m)); % imag
    end    
    
    if ns(1)>1 && ns(2)>1
        Ylmr{l+1} = reshape(Ylmr{l+1},[ns (l+1)]);
        Ylmi{l+1} = reshape(Ylmi{l+1},[ns (l+1)]);
    end    
end

