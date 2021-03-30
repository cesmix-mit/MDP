function Hy = HD_Hessian_Call(func,x)
    hx          = hyperdual(x);
    n           = size(x,1);
    Hy          = zeros(n,n);
    for j=1:n
        for k=1:n
            hxjk        = hx;
            hxjk.dx1(j)	= 1;
            hxjk.dx2(k) = 1;
            hyjk      	= func(hxjk);
            Hy(j,k)     = hyjk.dx1x2;
            Hy(k,j)     = hyjk.dx1x2;
        end
    end
    
end