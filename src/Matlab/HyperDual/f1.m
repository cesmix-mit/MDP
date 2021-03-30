function y = f1(x)
    y   = x(1) + x(2).^2 .* x(3) - x(1)./x(3) + x(2).^x(1);
    % x(2) must be positive
end