function [b, bx, by, bz] = ssbispectrumderiv(ar, ai, arx, ary, arz, aix, aiy, aiz, cg, indl, indm, rowm)
% SHSUM computes spherical harmonic spectra
%   ar    : real part of the spherical shell harmonic sums
%   ai    : imag part of the spherical shell harmonic sums
%   b     : the bispectrum of the spherical shell 
%   cg    : Clebsch-Gordan coefficients
%  indl   : indices where the bispectrum components are nonzero and unique
%  indm   : tuple indices (m2, m1, m)
%  rowm   : the number of Clebsch-Gordan coefficients for each tuple (l2, l1, l) 


I = size(indl,1);
K = size(ar{1},1);
indk = getindk(K);

K2 = K*(K+1)/2;
nx = size(arx{1},1); % number of neighbors

b = zeros(I*K2,1);
bx = zeros(nx,K2*I);
by = zeros(nx,K2*I);
bz = zeros(nx,K2*I);
for k = 1:K2
    k2 = indk(k,1);
    k1 = indk(k,2);    
    for i = 1:I
        n = i+(k-1)*I; % index for the nth basis function
        l2 = indl(i,1);
        l1 = indl(i,2);
        l = indl(i,3);     
        tmp = 0;
        nm = rowm(i+1)-rowm(i);
        for j = 1:nm % perform summations over (m,m1,m2) tuple indices
            c = cg(rowm(i)+j);
            m2 = indm(rowm(i)+j,1);
            m1 = indm(rowm(i)+j,2);
            m = indm(rowm(i)+j,3);
            
            a1 = ar{l+1}(k1,m+l+1);
            b1 = ai{l+1}(k1,m+l+1);
            a2 = ar{l1+1}(k2,m1+l1+1);
            b2 = ai{l1+1}(k2,m1+l1+1);
            a3 = ar{l2+1}(k2,m2+l2+1);
            b3 = ai{l2+1}(k2,m2+l2+1);      
            tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                 
            
            t1 = arx{l+1}(:,k1,m+l+1)*a2*a3 + a1*arx{l1+1}(:,k2,m1+l1+1)*a3 + a1*a2*arx{l2+1}(:,k2,m2+l2+1);                
            t2 = arx{l1+1}(:,k2,m1+l1+1)*b1*b3 + a2*aix{l+1}(:,k1,m+l+1)*b3 + a2*b1*aix{l2+1}(:,k2,m2+l2+1);
            t3 = arx{l2+1}(:,k2,m2+l2+1)*b1*b2 + a3*aix{l+1}(:,k1,m+l+1)*b2 + a3*b1*aix{l1+1}(:,k2,m1+l1+1);
            t4 = arx{l+1}(:,k1,m+l+1)*b2*b3 + a1*aix{l1+1}(:,k2,m1+l1+1)*b3 + a1*b2*aix{l2+1}(:,k2,m2+l2+1);
            bx(:,n) = bx(:,n) + c*(t1 + t2 + t3 - t4);
            
            t1 = ary{l+1}(:,k1,m+l+1)*a2*a3 + a1*ary{l1+1}(:,k2,m1+l1+1)*a3 + a1*a2*ary{l2+1}(:,k2,m2+l2+1);                
            t2 = ary{l1+1}(:,k2,m1+l1+1)*b1*b3 + a2*aiy{l+1}(:,k1,m+l+1)*b3 + a2*b1*aiy{l2+1}(:,k2,m2+l2+1);
            t3 = ary{l2+1}(:,k2,m2+l2+1)*b1*b2 + a3*aiy{l+1}(:,k1,m+l+1)*b2 + a3*b1*aiy{l1+1}(:,k2,m1+l1+1);
            t4 = ary{l+1}(:,k1,m+l+1)*b2*b3 + a1*aiy{l1+1}(:,k2,m1+l1+1)*b3 + a1*b2*aiy{l2+1}(:,k2,m2+l2+1);
            by(:,n) = by(:,n) + c*(t1 + t2 + t3 - t4);                        
            
            t1 = arz{l+1}(:,k1,m+l+1)*a2*a3 + a1*arz{l1+1}(:,k2,m1+l1+1)*a3 + a1*a2*arz{l2+1}(:,k2,m2+l2+1);                
            t2 = arz{l1+1}(:,k2,m1+l1+1)*b1*b3 + a2*aiz{l+1}(:,k1,m+l+1)*b3 + a2*b1*aiz{l2+1}(:,k2,m2+l2+1);
            t3 = arz{l2+1}(:,k2,m2+l2+1)*b1*b2 + a3*aiz{l+1}(:,k1,m+l+1)*b2 + a3*b1*aiz{l1+1}(:,k2,m1+l1+1);
            t4 = arz{l+1}(:,k1,m+l+1)*b2*b3 + a1*aiz{l1+1}(:,k2,m1+l1+1)*b3 + a1*b2*aiz{l2+1}(:,k2,m2+l2+1);
            bz(:,n) = bz(:,n) + c*(t1 + t2 + t3 - t4);                        
        end                  
        b(n) = tmp;               
    end        
end


