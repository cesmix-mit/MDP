function g = sphericalbessel2(r,L,K)

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
f = zeros(nx,K,N);

n = 1;
for k = 1:K
    x = x0(n,k)*r;
    g(:,k,n) = cos(x)./x;
end

n = 2;
for k = 1:K
    x = x0(n,k)*r;
    g(:,k,n) = cos(x)./(x.^2) + sin(x)./x;
end

for n=3:N        
    for k = 1:K
        x = x0(n,k)*r;
        f(:,k,1) = -cos(x)./x;
        f(:,k,2) = -cos(x)./(x.^2) - sin(x)./x;
    end    
    for j = 3:n
        for k = 1:K
            x = x0(n,k)*r;
            f(:,k,j) = ((2*j-3)./x).*f(:,k,j-1) - f(:,k,j-2);        
        end
    end
    g(:,:,n) = -f(:,:,n);
end

end




