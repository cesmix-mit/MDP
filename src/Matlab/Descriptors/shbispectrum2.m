function b = shbispectrum2(a)
% SHBISPECTRUM2 computes the bispectrum of spherical harmonics
%   ar    : real part of the spherical harmonic sums
%   ai    : imag part of the spherical harmonic sums
%   b     : the bispectrum

N = length(a);
L = N-1;

% factorial table
fac = factable;

b = zeros(N,N,N);
for l = 0:L
    for l1 = 0:L
        for l2 = 0:L
            tmp = 0;
            for m = -l:1:l
                for m1 = -l1:1:l1
                    for m2 = -l2:1:l2
                        cg = clebschgordan(l2,m2,l1,m1,l,m,fac);
                        if cg~=0
                            tmp = tmp + conj(a{l+1}(m+l+1))*cg...
                                       *(a{l1+1}(m1+l1+1))...
                                       *(a{l2+1}(m2+l2+1));
                        end                                
                    end
                end
            end
            if abs(tmp)>=1e-11
                if abs(real(tmp)) >= abs(imag(tmp)) 
                    b(l2+1,l1+1,l+1) = real(tmp);
                else
                    b(l2+1,l1+1,l+1) = imag(tmp);            
                end
            end
        end
    end    
end
