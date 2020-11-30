function [c,indm,rowm] = cgcoefficients(ind)
% cgcoefficients computes nonzero Clebsch-Gordan coefficients
% ind  : tuple indices (l2, l1, l) 
% c    : Clebsch-Gordan coefficients
% indm : tuple indices (m2, m1, m)
% rowm : the number of Clebsch-Gordan coefficients for each tuple (l2, l1, l) 

N = size(ind,1);

nc = 0;
for n = 1:N
    l2 = ind(n,1);
    l1 = ind(n,2);
    l = ind(n,3);
    for m = -l:1:l
        for m1 = -l1:1:l1
            for m2 = -l2:1:l2
                if ~(m1+m2-m ~= 0 || l < abs(l1-l2) || l > l1+l2)                        
                    nc = nc + 1;
                end
            end
        end
    end
end

% factorial table
fac = factable;

c = zeros(nc,1);
indm = zeros(nc,3);
nc = 1;
rowm = zeros(N,1);
for n = 1:N
    l2 = ind(n,1);
    l1 = ind(n,2);
    l = ind(n,3);
    count = 0;
    for m = -l:1:l
        for m1 = -l1:1:l1
            for m2 = -l2:1:l2
                if ~(m1+m2-m ~= 0 || l < abs(l1-l2) || l > l1+l2)                                            
                    %---compute valid k values for summation---%
                    k = max([0,l2-l-m1,l1-l+m2]):min([l1+l2-l,l1-m1,l2+m2]);

                    %---compute coefficient---%
                    c(nc) = sqrt((2*l+1)*fac(1+l+l1-l2)*fac(1+l+l2-l1)*fac(1+l1+l2-l)/fac(1+l1+l2+l+1))*...
                            sqrt(fac(1+l+m)*fac(1+l-m)*fac(1+l1+m1)*fac(1+l1-m1)*fac(1+l2+m2)*fac(1+l2-m2))*...
                            sum(((-1).^k)./(fac(1+k).*fac(1+l1+l2-l-k).*fac(1+l1-m1-k).*fac(1+l2+m2-k).*fac(1+l-l2+m1+k).*fac(1+l-l1-m2+k)));
                    
                    indm(nc,:) = [m2 m1 m];   
                    nc = nc + 1;
                    count = count + 1;
                end                        
            end
        end
    end
    rowm(n) = count;
end
rowm = cumsum([0; rowm]);
