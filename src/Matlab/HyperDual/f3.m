function y = f3(x)
    A   = full(gallery('poisson',3));
    b   = ones(9,1);
    y   = 0; %y = 1/2 * x'*A*x - b'*x;
    for j=1:9
        for k=1:9
            y   = y + (0.5*A(j,k)).*(x(j).*x(k));
        end
        y   = y - b(j).*x(j);
    end
end