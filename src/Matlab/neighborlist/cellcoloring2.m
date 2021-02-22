function cellcolor = cellcoloring2(nc,bsize)

dim = length(nc);
if dim==2
    nc1 = nc(1);
    nc2 = nc(2);
    bs1 = bsize(1);
    bs2 = bsize(2);
    nb1 = ceil(nc1/bs1); % number of blocks
    nb2 = ceil(nc2/bs2); % number of blocks
        
    cellcolor = zeros(nc1,nc2);
    bid = reshape(0:1:((bs1*bs2)-1),bs1, bs2);
    
    for i = 1:nb1
        if i==1, i1 = 1; else, i1 = i1 + bs1; end
        i2 = min(i1 + bs1 - 1, nc1);            
        for j = 1:nb2
            if j==1, j1 = 1; else, j1 = j1 + bs2; end
            j2 = min(j1 + bs2 - 1, nc2);          
            blockcolor = bid((1:(i2-i1+1)),(1:(j2-j1+1)));
            cellcolor(i1:i2,j1:j2) = blockcolor;            
        end
    end    
elseif dim==3
    nc1 = nc(1);
    nc2 = nc(2);
    nc3 = nc(3);
    bs1 = bsize(1);
    bs2 = bsize(2);
    bs3 = bsize(3);
    nb1 = ceil(nc1/bs1); % number of blocks
    nb2 = ceil(nc2/bs2); % number of blocks
    nb3 = ceil(nc3/bs3); % number of blocks
    
    cellcolor = zeros(nc1,nc2,nc3);
    bid = reshape(0:1:((bs1*bs2*bs3)-1),bs1, bs2, bs3);
        
    for i = 1:nb1
        if i==1, i1 = 1; else, i1 = i1 + bs1; end
        i2 = min(i1 + bs1 - 1, nc1);            
        for j = 1:nb2
            if j==1, j1 = 1; else, j1 = j1 + bs2; end
            j2 = min(j1 + bs2 - 1, nc2);          
            for k = 1:nb3
                if k==1, k1 = 1; else, k1 = k1 + bs3; end
                k2 = min(k1 + bs3 - 1, nc3);                        
                blockcolor = bid((1:(i2-i1+1)),(1:(j2-j1+1)),(1:(k2-k1+1)));
                cellcolor(i1:i2,j1:j2,k1:k2) = blockcolor;            
            end
        end
    end        
end

end



