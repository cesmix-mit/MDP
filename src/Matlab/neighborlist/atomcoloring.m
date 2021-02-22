function [clist,atomcolor] = atomcoloring(xi, P2S, cellcolor, eta1, eta2, eta3, inum, gnum)

natom = inum+gnum;
dim = size(xi,1);
clist = zeros(dim,natom);

colormax = max(cellcolor(:))+1;
atomcolor = zeros(colormax,natom);

if dim==2
    nc1 = length(eta1)-1;
    nc2 = length(eta2)-1;    
    
    for i=1:natom
        xt = P2S*xi(:,i);
        for j1 = 1:nc1
            if (eta1(j1) <= xt(1)) && (xt(1)<= eta1(j1+1))
                break;
            end
        end
        for j2 = 1:nc2
            if (eta2(j2) <= xt(2)) && (xt(2)<= eta2(j2+1))
                break;
            end
        end
        clist(1,i) = j1;
        clist(2,i) = j2;
        color = cellcolor(j1,j2);
        atomcolor(color+1,i) = 1;
    end
else
    nc1 = length(eta1)-1;
    nc2 = length(eta2)-1;
    nc3 = length(eta3)-1;
    
    for i=1:(inum+gnum)
        xt = P2S*xi(:,i);
        for j1 = 1:nc1
            if (eta1(j1) <= xt(1)) && (xt(1)<= eta1(j1+1))
                break;
            end
        end
        for j2 = 1:nc2
            if (eta2(j2) <= xt(2)) && (xt(2)<= eta2(j2+1))
                break;
            end
        end
        for j3 = 1:nc3
            if (eta3(j3) <= xt(3)) && (xt(3)<= eta3(j3+1))
                break;
            end
        end
        clist(1,i) = j1;
        clist(2,i) = j2;
        clist(3,i) = j3;
        color = cellcolor(j1,j2,j3);
        atomcolor(color+1,i) = 1;
    end    
end


