function b = ssbispectrum(ar, ai, cg, indl, indm, rowm)
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

%imagi = sqrt(-1);

K2 = K*(K+1)/2;
b = zeros(I*K2,1);
for k = 1:K2
    k2 = indk(k,1);
    k1 = indk(k,2);    
    for i = 1:I
        l2 = indl(i,1);
        l1 = indl(i,2);
        l = indl(i,3);     
        tmp = 0;
        %tmr = 0;
        %tmi = 0;
        nm = rowm(i+1)-rowm(i);
        for j = 1:nm
            c = cg(rowm(i)+j);
            m2 = indm(rowm(i)+j,1);
            m1 = indm(rowm(i)+j,2);
            m = indm(rowm(i)+j,3);
%             tmp = tmp + (ar{l+1}(k1,m+l+1) - imagi*ai{l+1}(k1,m+l+1))*c...
%                        *(ar{l1+1}(k2,m1+l1+1) + imagi*ai{l1+1}(k2,m1+l1+1))...
%                        *(ar{l2+1}(k2,m2+l2+1) + imagi*ai{l2+1}(k2,m2+l2+1));                                   
            a1 = ar{l+1}(k1,m+l+1);
            b1 = ai{l+1}(k1,m+l+1);
            a2 = ar{l1+1}(k2,m1+l1+1);
            b2 = ai{l1+1}(k2,m1+l1+1);
            a3 = ar{l2+1}(k2,m2+l2+1);
            b3 = ai{l2+1}(k2,m2+l2+1);      
            %tmp = tmp + c*(a1*(a2*a3 - b2*b3) + b1*(a2*b3 + a3*b2));
            tmp = tmp + c*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                 
            %tmi = tmi + c*(a1*a2*b3 + a1*a3*b2 + b1*b2*b3 - a2*a3*b1);                   
        end                  
        b(i+(k-1)*I) = tmp;               
    end        
end

% b = zeros(K*(K+1)/2,I);
% k = 0;
% for k1 = 1:K
%     for k2 = 1:k1    
%         k = k + 1;
%         for i = 1:I
%             l2 = indl(i,1);
%             l1 = indl(i,2);
%             l = indl(i,3);     
%             tmp = 0;
%             nm = rowm(i+1)-rowm(i);
%             for j = 1:nm
%                 c = cg(rowm(i)+j);
%                 m2 = indm(rowm(i)+j,1);
%                 m1 = indm(rowm(i)+j,2);
%                 m = indm(rowm(i)+j,3);
%                 tmp = tmp + (ar{l+1}(k1,m+l+1) - imagi*ai{l+1}(k1,m+l+1))*c...
%                            *(ar{l1+1}(k2,m1+l1+1) + imagi*ai{l1+1}(k2,m1+l1+1))...
%                            *(ar{l2+1}(k2,m2+l2+1) + imagi*ai{l2+1}(k2,m2+l2+1));                
%             end                        
%             b(k,i) = tmp;               
%         end        
%     end    
% end
% 
% 
