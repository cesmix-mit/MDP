function [ar, ai]= shsum(x,y,z,L)
% SHSUM computes spherical harmonic sums
% (x,y,z) : Cartesian coordinates
%    L    : degree of spherical harmonic basis
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums

% Spherical coordinates
[the,phi,~] = cart2sphere(x,y,z);

[Ylmr, Ylmi] = sphericalharmonics(the,phi,L);

N = L + 1;
ar = cell(N,1);
ai = cell(N,1);

for l = 1:N
    ar{l} = reshape(sum(Ylmr{l},1), [1 l]);
    ai{l} = reshape(sum(Ylmi{l},1), [1 l]);    
    if l>1
        n = size(ar{l},1);
        br = zeros(n,l-1);
        bi = zeros(n,l-1);
        for i = 2:l
            br(:,i-1) = ((-1)^(i-1))*ar{l}(:,i);
            bi(:,i-1) = ((-1)^i)*ai{l}(:,i);                
        end
        ar{l} = [br(:,end:-1:1) ar{l}];
        ai{l} = [bi(:,end:-1:1) ai{l}];
    end
end


