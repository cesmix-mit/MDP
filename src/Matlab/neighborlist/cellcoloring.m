function cellcolor = cellcoloring(nc,bsize)

dim = length(nc);
if dim==2
    nc1 = nc(1);
    nc2 = nc(2);
    bs1 = bsize(1);
    bs2 = bsize(2);
    bid = reshape(0:1:((bs1*bs2)-1),bs1, bs2);
    cellcolor = zeros(nc1,nc2);
    
    nb1 = ceil(nc1/bs1); % number of blocks
    id1 = bs1*ones(1,nb1+1);
    id1(1) = 0;
    id1(end) = nc1 - (nb1-1)*bs1;
    for i=2:(nb1+1)
        id1(i) = id1(i) + id1(i-1);
    end
    id1 = id1 + 1;
    
    nb2 = ceil(nc2/bs2); % number of blocks
    id2 = bs2*ones(1,nb2+1);
    id2(1) = 0;
    id2(end) = nc2 - (nb2-1)*bs2;
    for i=2:(nb2+1)
        id2(i) = id2(i) + id2(i-1);
    end
    id2 = id2 + 1;
    
    for i = 1:nb1
        ni = id1(i+1) - id1(i);
        ii = id1(i):(id1(i+1)-1);
        for j = 1:nb2
            nj = id2(j+1) - id2(j);
            jj = id2(j):(id2(j+1)-1);
            blockcolor = bid((1:ni),(1:nj));
            cellcolor(ii,jj) = blockcolor;
        end
    end
elseif dim==3
    nc1 = nc(1);
    nc2 = nc(2);
    nc3 = nc(3);
    bs1 = bsize(1);
    bs2 = bsize(2);
    bs3 = bsize(3);
    
    bid = reshape(0:1:((bs1*bs2*bs3)-1),bs1, bs2, bs3);
    cellcolor = zeros(nc1,nc2,nc3);
    
    nb1 = ceil(nc1/bs1); % number of blocks
    id1 = bs1*ones(1,nb1+1);
    id1(1) = 0;
    id1(end) = nc1 - (nb1-1)*bs1;
    for i=2:(nb1+1)
        id1(i) = id1(i) + id1(i-1);
    end
    id1 = id1 + 1;
    
    nb2 = ceil(nc2/bs2); % number of blocks
    id2 = bs2*ones(1,nb2+1);
    id2(1) = 0;
    id2(end) = nc2 - (nb2-1)*bs2;
    for i=2:(nb2+1)
        id2(i) = id2(i) + id2(i-1);
    end
    id2 = id2 + 1;
    
    nb3 = ceil(nc3/bs3); % number of blocks
    id3 = bs3*ones(1,nb3+1);
    id3(1) = 0;
    id3(end) = nc3 - (nb3-1)*bs3;
    for i=2:(nb3+1)
        id3(i) = id3(i) + id3(i-1);
    end
    id3 = id3 + 1;
        
    for i = 1:nb1
        ni = id1(i+1) - id1(i);
        ii = id1(i):(id1(i+1)-1);
        for j = 1:nb2
            nj = id2(j+1) - id2(j);
            jj = id2(j):(id2(j+1)-1);
            for k = 1:nb3
                nk = id3(k+1) - id3(k);
                kk = id3(k):(id3(k+1)-1);            
                blockcolor = bid((1:ni),(1:nj),(1:nk));
                cellcolor(ii,jj,kk) = blockcolor;
            end
        end
    end    
end

end



