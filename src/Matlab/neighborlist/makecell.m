function [eta1, eta2, eta3] = makecell(smin, smax, ds)

dim = length(smin);

if dim==2
    ne1 = floor(1/ds(1));
    ne2 = floor(1/ds(2));
    eta1 = [smin(1) linspace(0,1,ne1+1) smax(1)];
    eta2 = [smin(2) linspace(0,1,ne2+1) smax(2)];  
    eta3 = 0;
elseif dim==3
    ne1 = floor(1/ds(1));
    ne2 = floor(1/ds(2));
    ne3 = floor(1/ds(3));
    eta1 = [smin(1) linspace(0,1,ne1+1) smax(1)];
    eta2 = [smin(2) linspace(0,1,ne2+1) smax(2)];                        
    eta3 = [smin(3) linspace(0,1,ne3+1) smax(3)];    
end


