function [b, indl, ar, ai, p] = shspectrum(x,y,z,L)
% SHSUM computes spherical harmonic spectra
% (x,y,z) : Cartesian coordinates
%    L    : degree of spherical harmonic basis
%   p     : the power spectrum
%   b     : the bispectrum
%  indl   : indices where the bispectrum components are nonzero and unique
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums

% compute spherical harmonic sums
[ar, ai] = shsum(x,y,z,L);

% compute the bispectrum of spherical harmonics
b = shbispectrum(ar, ai);

% find non-zero unique bispectrum components
indl = uniquebispectrum(b);

if nargout >= 5
    % compute the power spectrum of spherical harmonics
    p = shpower(ar, ai);
end



