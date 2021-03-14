function p = shpower(ar, ai)
% SHPOWER computes the power spectrum of spherical harmonics
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums
%   p     : the power spectrum

N = length(ar);
p = zeros(N,1);
for l = 1:N
    p(l) = sum(ar{l}.*ar{l}+ai{l}.*ai{l});
end
