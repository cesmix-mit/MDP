function Hyperdual_Validation()

    clear;
    close all;
    clc;
    
    number      = 5;
    % set x, set n, set func
    switch number
        case 1
            n           = 3;
            x           = randn(n,1);
            x(2)        = abs(x(2));    % f1 requires that x(2) is positive.
            func       	= @(x) f1(x);
        case 2
            n           = 1;
            x           = randn(n,1);
            x(1)        = abs(x(1));    % f2 requires that x(1) is positive.
            func       	= @(x) f2(x);
        case 3
            n           = 9;
            x           = randn(n,1);
            func       	= @(x) f3(x);
        case 4
            n           = 2;
            x           = randn(n,1);
            func       	= @(x) f4(x);
        case 5
            n           = 3;
            x           = randn(n,1);
            x(2)        = abs(x(2));    % f1 requires that x(2) is positive.
            func       	= @(x) f5(x);
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Jacobian
    Dy          = HD_Jacobian_Call(func,x);
    x
    Dy
    % validate Jacobian
    Dz          = zeros(1,n);
    for j=1:n
        hxj         = x;
        hxj(j)      = hxj(j) + 1e-50 * sqrt(-1);
        hyj         = func(hxj);
        Dz(1,j)     = 1e50 * imag(hyj);
    end
    
    disp('Jacobian')
    err_Dy  = norm(Dy-Dz,'inf')/norm(Dz,'inf') %#ok<NOPRT,NASGU>
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Hessian
    Hy          = HD_Hessian_Call(func,x);
    
    % validate Hessian
    Hz          = zeros(n,n);
    for j=1:n
        hxj         = x;
        hxj(j)      = hxj(j) + 1e-50 * sqrt(-1);
        hyj         = HD_Jacobian_Call(func,hxj);
        Hz(j,:)     = 1e50 * imag(hyj);
    end
    
    disp('Hessian')
    err_Hy  = norm(Hy-Hz,'inf')/norm(Hz,'inf') %#ok<NOPRT,NASGU>
    
end