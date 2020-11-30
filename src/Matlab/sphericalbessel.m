function g = sphericalbessel(r,L,K)

% roots of spherical Bessel functions
N = L+1;
kind = 2;
n = (1:N)' - 1/2;
x0 = besselzero(n, K, kind);

r = r(:);
nx = length(r);

% Spherical Bessel functions
g = zeros(nx,K,N);
for n=1:N
    for k = 1:K
        x = x0(n,k)*r;
        g(:,k,n) = -sqrt(pi./(2*x)).*bessely(n-1/2,x);
    end
end

