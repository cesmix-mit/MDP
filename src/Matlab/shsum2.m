function a = shsum2(x,y,z,L)
% SHSUM2 computes spherical harmonic sums
% (x,y,z) : Cartesian coordinates
%    L    : degree of spherical harmonic basis
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums

[Ylm, Plm] = sphericalharmonicbasis2(x,y,z,L);

N = L + 1;
a = cell(N,1);

for l = 1:N
    a{l} = reshape(sum(Ylm{l},1), [1 l]);    
    if l>1
        b = reshape(sum(Plm{l},1), [1 l-1]); % (-1)^(m+1)*ar        
        a{l} = [b(:,end:-1:1) a{l}];        
    end
end


