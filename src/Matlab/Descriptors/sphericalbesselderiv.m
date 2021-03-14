function [g, dg] = sphericalbesselderiv(r,L,K)

% syms n x 
% y(1) = -cos(x)/x;
% y(2) = -cos(x)/(x*x) - sin(x)/x;
% for n = 2:10
%     y(n+1) = ((2*n-1)/x)*y(n) - y(n-1);
% end

% roots of spherical Bessel functions
N = L+1;
kind = 2;
n = (1:N)' - 1/2;
x0 = besselzero(n, K, kind);

r = r(:);
nx = length(r);

% Spherical Bessel functions
g = zeros(nx,K,N);
dg = zeros(nx,K,N);
f = zeros(nx,K,N);
df = zeros(nx,K,N);

n = 1;
for k = 1:K
    x = x0(n,k)*r;
    g(:,k,n) = cos(x)./x;
    dg(:,k,n) = x0(n,k)*(-cos(x)./x.^2 - sin(x)./x);
end

n = 2;
for k = 1:K
    x = x0(n,k)*r;
    g(:,k,n) = cos(x)./(x.^2) + sin(x)./x;
    dg(:,k,n) = x0(n,k)*(cos(x)./x - (2*cos(x))./x.^3 - (2*sin(x))./x.^2);
end

for n=3:N        
    for k = 1:K
        x = x0(n,k)*r;
        f(:,k,1) = -cos(x)./x;
        f(:,k,2) = -cos(x)./(x.^2) - sin(x)./x;
        df(:,k,1) = x0(n,k)*(cos(x)./x.^2 + sin(x)./x);
        df(:,k,2) = x0(n,k)*((2*cos(x))./x.^3 - cos(x)./x  + (2*sin(x))./x.^2);
    end    
    for j = 3:n
        for k = 1:K
            x = x0(n,k)*r;
            f(:,k,j) = ((2*j-3)./x).*f(:,k,j-1) - f(:,k,j-2);        
            df(:,k,j) = ((2*j-3)./x).*df(:,k,j-1) - x0(n,k)*((2*j-3)./x.^2).*f(:,k,j-1) - df(:,k,j-2);        
        end
    end
    g(:,:,n) = -f(:,:,n);
    dg(:,:,n) = -df(:,:,n);
end

end




