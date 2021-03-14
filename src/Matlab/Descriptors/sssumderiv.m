function [ar, ai, arx, ary, arz, aix, aiy, aiz]= sssumderiv(x,y,z,K,L)
% SHSUM computes spherical shell harmonic sums
% (x,y,z) : Cartesian coordinates
%    K    : the number of radial functions
%    L    : degree of spherical shell harmonic basis
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums

[Sr, Si, arx, ary, arz, aix, aiy, aiz] = ssharmonicsderiv(x,y,z,K,L);

nx = length(x(:));
N = L + 1;
ar = cell(N,1);
ai = cell(N,1);

for l = 1:N
    ar{l} = reshape(sum(Sr{l},1), [K l]);
    ai{l} = reshape(sum(Si{l},1), [K l]);    
    if l>1        
        br = zeros(K,l-1);
        bi = zeros(K,l-1);
        brx = zeros(nx,K,l-1);
        bix = zeros(nx,K,l-1);
        bry = zeros(nx,K,l-1);
        biy = zeros(nx,K,l-1);
        brz = zeros(nx,K,l-1);
        biz = zeros(nx,K,l-1);
        for i = 2:l
            br(:,i-1) = ((-1)^(i-1))*ar{l}(:,i);
            bi(:,i-1) = ((-1)^i)*ai{l}(:,i);                
            brx(:,:,i-1) = ((-1)^(i-1))*arx{l}(:,:,i);
            bix(:,:,i-1) = ((-1)^i)*aix{l}(:,:,i);                
            bry(:,:,i-1) = ((-1)^(i-1))*ary{l}(:,:,i);
            biy(:,:,i-1) = ((-1)^i)*aiy{l}(:,:,i);                
            brz(:,:,i-1) = ((-1)^(i-1))*arz{l}(:,:,i);
            biz(:,:,i-1) = ((-1)^i)*aiz{l}(:,:,i);                
        end
        ar{l} = [br(:,end:-1:1) ar{l}];
        ai{l} = [bi(:,end:-1:1) ai{l}];
        arx{l} = cat(3, brx(:,:,end:-1:1), arx{l});
        aix{l} = cat(3, bix(:,:,end:-1:1), aix{l});
        ary{l} = cat(3, bry(:,:,end:-1:1), ary{l});
        aiy{l} = cat(3, biy(:,:,end:-1:1), aiy{l});
        arz{l} = cat(3, brz(:,:,end:-1:1), arz{l});
        aiz{l} = cat(3, biz(:,:,end:-1:1), aiz{l});
    end    
end


