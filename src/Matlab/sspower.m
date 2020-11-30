function p = sspower(ar, ai)
% SSPOWER computes the power spectrum of spherical shell harmonics
%   ar    : real part of the spherical shell harmonic sums
%   ai    : imag part of the spherical shell harmonic sums
%   p     : the power spectrum

L = length(ar);
K = size(ar{1},1);
K2 = K*(K+1)/2;
p = zeros(K2*L,1);

indk = getindk(K);
for k = 1:K2
    j = indk(k,1);
    i = indk(k,2);
    for l = 1:L    
        p(l + (k-1)*L) = (sum(ar{l}(i,:).*ar{l}(j,:)+ai{l}(i,:).*ai{l}(j,:)));
    end    
end


% for l = 1:L
%     count = 0;
%     for i = 1:K
%         for j = 1:i
%             count = count + 1;
%             p(count+(l-1)*M) = (sum(ar{l}(i,:).*ar{l}(j,:)+ai{l}(i,:).*ai{l}(j,:)));
%         end
%     end
% end

