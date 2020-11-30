function [p,px,py,pz] = sspowerderiv(ar, ai, arx, ary, arz, aix, aiy, aiz)
% SSPOWER computes the power spectrum of spherical shell harmonics
%   ar    : real part of the spherical shell harmonic sums
%   ai    : imag part of the spherical shell harmonic sums
%   p     : the power spectrum

L = length(ar);
K = size(ar{1},1);
K2 = K*(K+1)/2;
p = zeros(K2*L,1);

nx = size(arx{1},1);
px = zeros(nx,K2*L);
py = zeros(nx,K2*L);
pz = zeros(nx,K2*L);

indk = getindk(K);
for k = 1:K2
    j = indk(k,1);
    i = indk(k,2);
    for l = 1:L    
        %p(l + (k-1)*L) = (sum(ar{l}(i,:).*ar{l}(j,:)+ai{l}(i,:).*ai{l}(j,:)));
        n = l + (k-1)*L;
        for m = 1:(2*l-1)
            p(n) = p(n) + ar{l}(i,m)*ar{l}(j,m)+ai{l}(i,m)*ai{l}(j,m);
            px(:,n) = px(:,n) + 2*(ar{l}(i,m)*arx{l}(:,j,m)+ai{l}(i,m)*aix{l}(:,j,m));
            py(:,n) = py(:,n) + 2*(ar{l}(i,m)*ary{l}(:,j,m)+ai{l}(i,m)*aiy{l}(:,j,m));
            pz(:,n) = pz(:,n) + 2*(ar{l}(i,m)*arz{l}(:,j,m)+ai{l}(i,m)*aiz{l}(:,j,m));
        end
    end    
end


