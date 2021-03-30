function Dy = HD_Jacobian_Call(func,x)
    hx          = hyperdual(x);
    n           = size(x,1);
    Dy          = zeros(1,n);
    for j=1:n
        hxj         = hx;
        hxj.dx1(j)  = 1;
        hyj         = func(hxj);
        Dy(1,j)     = hyj.dx1;
    end
    
end