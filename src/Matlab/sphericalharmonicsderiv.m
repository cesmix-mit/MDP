function [Ylmr, Ylmi, YlmrThe, YlmiThe, YlmrPhi, YlmiPhi] = sphericalharmonicsderiv(the,phi,L)
% (the,phi): spherical coordinates
% L : the degree of spherical harmonics 
%   Ylmr    : real part of the spherical harmonics
%   Ylmi    : imag part of the spherical harmonics

fac = factable;

the = the(:);
phi = phi(:);
np = length(phi);

% temporary arrays to compute the recurrence 
P = zeros([np (L+1)]); % l
tmp = zeros([np L]);   % l - 1
tmp(:,1) = 1;
dP = zeros([np (L+1)]); % l
dtmp = zeros([np L]);   % l - 1
dtmp(:,1) = 0;

% allocate memory
Ylmr = cell(L+1,1);
Ylmi = cell(L+1,1);
YlmrThe = cell(L+1,1);
YlmiThe = cell(L+1,1);
YlmrPhi = cell(L+1,1);
YlmiPhi = cell(L+1,1);

% l = 0
Ylmr{1} = sqrt(1/(4*pi))*ones(np,1);
Ylmi{1} = zeros(np,1);
YlmrThe{1} = zeros(np,1);
YlmiThe{1} = zeros(np,1);
YlmrPhi{1} = zeros(np,1);
YlmiPhi{1} = zeros(np,1);

% precompute
costhe = cos(the);
a = -sin(the);
dcosthe = -sin(the);
da = -cos(the);

% l = 1
l = 1;
n = 2*l+1;    
m = 0:l;
C = sqrt(n*fac(l - m + 1)./(4*pi*fac(l + m + 1))); 
P(:,l) = costhe;    
P(:,l+1) = a;
dP(:,l) = dcosthe;    
dP(:,l+1) = da;

Ylmr{2} = zeros([np 2]);
Ylmi{2} = zeros([np 2]);    
Ylmr{2}(:,1) = (C(1)*cos(0))*P(:,1);    % m=0
Ylmr{2}(:,2) = C(2)*(cos(phi).*P(:,2)); % m=1
Ylmi{2}(:,1) = (C(1)*sin(0))*P(:,1);    % m=0
Ylmi{2}(:,2) = C(2)*(sin(phi).*P(:,2)); % m=1

YlmrThe{2} = zeros([np 2]);
YlmiThe{2} = zeros([np 2]);    
YlmrPhi{2} = zeros([np 2]);
YlmiPhi{2} = zeros([np 2]);    
YlmrThe{2}(:,1) = (C(1)*cos(0))*dcosthe;    % m=0
YlmrThe{2}(:,2) = C(2)*(cos(phi).*da);      % m=1
YlmiThe{2}(:,1) = (C(1)*sin(0))*dcosthe;    % m=0
YlmiThe{2}(:,2) = C(2)*(sin(phi).*da);      % m=1
YlmrPhi{2}(:,1) = 0;                        % m=0
YlmrPhi{2}(:,2) = -C(2)*(sin(phi).*P(:,2)); % m=1
YlmiPhi{2}(:,1) = 0;                        % m=0
YlmiPhi{2}(:,2) = C(2)*(cos(phi).*P(:,2));  % m=1

for l=2:L
    Ylmr{l+1} = zeros([np l+1]);
    Ylmi{l+1} = zeros([np l+1]);    
    YlmrThe{l+1} = zeros([np l+1]);
    YlmiThe{l+1} = zeros([np l+1]);    
    YlmrPhi{l+1} = zeros([np l+1]);
    YlmiPhi{l+1} = zeros([np l+1]);    
    
    n = 2*l+1;    
    m = 0:l;
    C = sqrt(n*fac(l - m + 1)./(4*pi*fac(l + m + 1))); 
    
    Pll = P(:,l);    
    tmp(:,l) = Pll;    
    P(:,l) = (2*(l-1)+1)*(costhe.*Pll);   % P_{l+1}^{l} = (2*l+1)*x*P_l^l
    P(:,l+1) = (2*(l-1)+1)*(a.*Pll); % P_{l+1}^{l+1} = -(2*l+1)*sqrt(1-x.^2)*P_l^l
    dPll = dP(:,l);
    dtmp(:,l) = dPll;
    dP(:,l) = (2*(l-1)+1)*(costhe.*dPll + dcosthe.*Pll);   % P_{l+1}^{l} = (2*l+1)*x*P_l^l
    dP(:,l+1) = (2*(l-1)+1)*(a.*dPll + da.*Pll); % P_{l+1}^{l+1} = -(2*l+1)*sqrt(1-x.^2)*P_l^l
    for m = 1:(l-1)
        % (l-m+1)*P_{l+1}^m = (2*l+1)*(x.*P_l^m - (l+m)*P_{l-1}^m
        Pll = P(:,m);
        P(:,m) = ((2*(l-1)+1)*(costhe.*Pll) - (l+m-2)*tmp(:,m))/(l-m+1);
        tmp(:,m) = Pll;
        dPll = dP(:,m);
        dP(:,m) = ((2*(l-1)+1)*(costhe.*dPll+dcosthe.*Pll) - (l+m-2)*dtmp(:,m))/(l-m+1);
        dtmp(:,m) = dPll;
    end        
        
    for m=1:l+1
        Ylmr{l+1}(:,m) = C(m)*(cos((m-1)*phi).*P(:,m)); % real
        Ylmi{l+1}(:,m) = C(m)*(sin((m-1)*phi).*P(:,m)); % imag
        YlmrThe{l+1}(:,m) = C(m)*(cos((m-1)*phi).*dP(:,m)); % real
        YlmiThe{l+1}(:,m) = C(m)*(sin((m-1)*phi).*dP(:,m)); % imag
        YlmrPhi{l+1}(:,m) = -((m-1)*C(m))*(sin((m-1)*phi).*P(:,m)); % real
        YlmiPhi{l+1}(:,m) =  ((m-1)*C(m))*(cos((m-1)*phi).*P(:,m)); % imag        
    end        
end

