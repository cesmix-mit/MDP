function y = f5(x) % 3
    % x(2) must be positive
    f   = @(x,y) x.^2 .* (y.^2).^1.5;
    y = f( f1(x) , f4(x(1:2)) );
end