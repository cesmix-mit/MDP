function p = shpower2(a)
% SHPOWER2 computes the power spectrum of spherical harmonics
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums
%   p     : the power spectrum

N = length(a);
p = zeros(N,1);
for l = 1:N
    p(l) = a{l}*(a{l}'); 
end

