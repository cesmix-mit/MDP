function [ar, ai]= sssum(x,y,z,K,L)
% SHSUM computes spherical shell harmonic sums
% (x,y,z) : Cartesian coordinates
%    K    : the number of radial functions
%    L    : degree of spherical shell harmonic basis
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums

[Sr, Si] = sphericalshellbasis(x,y,z,K,L);

N = L + 1;
ar = cell(N,1);
ai = cell(N,1);

for l = 1:N
    ar{l} = reshape(sum(Sr{l},1), [K l]);
    ai{l} = reshape(sum(Si{l},1), [K l]);    
    if l>1        
        br = zeros(K,l-1);
        bi = zeros(K,l-1);
        for i = 2:l
            br(:,i-1) = ((-1)^(i-1))*ar{l}(:,i);
            bi(:,i-1) = ((-1)^i)*ai{l}(:,i);                
        end
        ar{l} = [br(:,end:-1:1) ar{l}];
        ai{l} = [bi(:,end:-1:1) ai{l}];
    end    
end


