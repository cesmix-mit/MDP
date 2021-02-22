function [xi, ilist, n, k] = createilist(xi, ximages, ws, B2C)

dim = size(xi,1);
n = size(xi,2);
m = size(ximages,2);
ilist = (1:n);
k = 0;

if dim==2
    for i = 1:n
        for j = 2:m
            xj = xi(:,i) + ximages(:,j);
            xs = B2C*xj;        
            if (ws(1,1) <= xs(1)) && (xs(1) <= ws(1,3)) &&  (ws(2,1) <= xs(2)) && (xs(2) <= ws(2,3))
                k = k + 1;
                xi(:,n+k) = xj;
                ilist(n+k) = i;
            end
        end    
    end
else
    for i = 1:n
        for j = 2:m
            xj = xi(:,i) + ximages(:,j);
            xs = B2C*xj;        
            if (ws(1,1) <= xs(1)) && (xs(1) <= ws(1,7)) &&  (ws(2,1) <= xs(2)) && (xs(2) <= ws(2,7)) &&  (ws(3,1) <= xs(3)) && (xs(3) <= ws(3,7))
                k = k + 1;
                xi(:,n+k) = xj;
                ilist(n+k) = i;
            end
        end    
    end    
end



