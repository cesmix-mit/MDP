function [Sr, Si, Srx, Sry, Srz, Six, Siy, Siz] = ssharmonicsderiv(x,y,z,K,L)
% (x,y,z): Cartesian coordinates
% K : the number of radial functions
% L : the degree of spherical harmonics 
%   Ylmr    : real part of the spherical shell harmonics
%   Ylmi    : imag part of the spherical shell harmonics

N = L+1;

x = x(:);
y = y(:);
z = z(:);
nx = length(x);

% Spherical coordinates
[the,phi,r,Thex,They,Thez,Phix,Phiy,Phiz,Rx,Ry,Rz] = cart2spherederiv(x,y,z);

% Spherical harmonics
[Ylmr, Ylmi, YlmrThe, YlmiThe, YlmrPhi, YlmiPhi] = sphericalharmonicsderiv(the,phi,L);

% spherical Bessel functions
[g, gR] = sphericalbesselderiv(r,L,K);
gx = zeros(nx,K,N);
gy = zeros(nx,K,N);
gz = zeros(nx,K,N);
for l=1:N   
    for k = 1:K
        gx(:,k,l) = gR(:,k,l).*Rx;
        gy(:,k,l) = gR(:,k,l).*Ry;
        gz(:,k,l) = gR(:,k,l).*Rz;
    end
end

Sr = cell(N,1);
Si = cell(N,1);
Srx = cell(N,1);
Sry = cell(N,1);
Srz = cell(N,1);
Six = cell(N,1);
Siy = cell(N,1);
Siz = cell(N,1);

% spherical shell harmonics for l = 0
Sr{1} = zeros(nx,K);
Si{1} = zeros(nx,K);    
Srx{1} = zeros(nx,K);
Six{1} = zeros(nx,K);    
Sry{1} = zeros(nx,K);
Siy{1} = zeros(nx,K);    
Srz{1} = zeros(nx,K);
Siz{1} = zeros(nx,K);    
for k = 1:K
    Sr{1}(:,k) = sqrt(1/(4*pi))*g(:,k,1);
    Si{1}(:,k) = 0;
    Srx{1}(:,k) = sqrt(1/(4*pi))*(gx(:,k,1));
    Sry{1}(:,k) = sqrt(1/(4*pi))*(gy(:,k,1));
    Srz{1}(:,k) = sqrt(1/(4*pi))*(gz(:,k,1));
end       

% spherical shell harmonics for l > 0
for l=1:L   
    Sr{l+1} = zeros(nx,K,l+1);
    Si{l+1} = zeros(nx,K,l+1);    
    Srx{l+1} = zeros(nx,K,l+1);
    Six{l+1} = zeros(nx,K,l+1);    
    Sry{l+1} = zeros(nx,K,l+1);
    Siy{l+1} = zeros(nx,K,l+1);    
    Srz{l+1} = zeros(nx,K,l+1);
    Siz{l+1} = zeros(nx,K,l+1);    
    for m = 1:(l+1)    
        Ylmrx = (YlmrThe{l+1}(:,m).*Thex + YlmrPhi{l+1}(:,m).*Phix);
        Ylmry = (YlmrThe{l+1}(:,m).*They + YlmrPhi{l+1}(:,m).*Phiy);
        Ylmrz = (YlmrThe{l+1}(:,m).*Thez + YlmrPhi{l+1}(:,m).*Phiz);
        Ylmix = (YlmiThe{l+1}(:,m).*Thex + YlmiPhi{l+1}(:,m).*Phix);
        Ylmiy = (YlmiThe{l+1}(:,m).*They + YlmiPhi{l+1}(:,m).*Phiy);
        Ylmiz = (YlmiThe{l+1}(:,m).*Thez + YlmiPhi{l+1}(:,m).*Phiz);
        for k = 1:K                
            Sr{l+1}(:,k,m) = g(:,k,l+1).*Ylmr{l+1}(:,m);
            Si{l+1}(:,k,m) = g(:,k,l+1).*Ylmi{l+1}(:,m);            
            Srx{l+1}(:,k,m) = (gx(:,k,l+1)).*Ylmr{l+1}(:,m) + g(:,k,l+1).*Ylmrx;
            Sry{l+1}(:,k,m) = (gy(:,k,l+1)).*Ylmr{l+1}(:,m) + g(:,k,l+1).*Ylmry;
            Srz{l+1}(:,k,m) = (gz(:,k,l+1)).*Ylmr{l+1}(:,m) + g(:,k,l+1).*Ylmrz;
            Six{l+1}(:,k,m) = (gx(:,k,l+1)).*Ylmi{l+1}(:,m) + g(:,k,l+1).*Ylmix;
            Siy{l+1}(:,k,m) = (gy(:,k,l+1)).*Ylmi{l+1}(:,m) + g(:,k,l+1).*Ylmiy;
            Siz{l+1}(:,k,m) = (gz(:,k,l+1)).*Ylmi{l+1}(:,m) + g(:,k,l+1).*Ylmiz;
        end
    end       
end

