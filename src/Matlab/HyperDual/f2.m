function y = f2(x)
    y   = x; % x must be positive
    for N=1:10
        y   = 0.5 * (y - x./y);
    end
end