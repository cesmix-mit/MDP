function y = f4(x)
    a = 1; b = 100;
    y = (a - x(1)).^2 + b .* (x(2)-x(1).^2).^2;
end