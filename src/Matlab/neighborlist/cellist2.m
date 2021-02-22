function [clist, c2ilist, c2inum, cellnum] = cellist2(xi, ws, P2S, inum, gnum)

dim = size(xi,1);
if dim==2
    ne1 = floor(1/abs(ws(1,1)));
    ne2 = floor(1/abs(ws(2,1)));
    eta1 = [ws(1,1) linspace(0,1,ne1+1) ws(1,3)];
    eta2 = [ws(2,1) linspace(0,1,ne2+1) ws(2,3)];    
    clist = zeros(1,inum+gnum);
    nc1 = length(eta1)-1;
    nc2 = length(eta2)-1;
    cellnum = [nc1 nc2];
    c2ilist = zeros(10,nc1*nc2);
    c2inum = zeros(1,nc1*nc2);
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
        l = j1 + (j2-1)*nc1;
        clist(i) = l;        
        c2inum(l) = c2inum(l) + 1;
        c2ilist(c2inum(l),l) = i;    
    end
else
    ne1 = floor(1/abs(ws(1,1)));
    ne2 = floor(1/abs(ws(2,1)));
    ne3 = floor(1/abs(ws(3,1)));
    eta1 = [ws(1,1) linspace(0,1,ne1+1) ws(1,7)];
    eta2 = [ws(2,1) linspace(0,1,ne2+1) ws(2,7)];    
    eta3 = [ws(3,1) linspace(0,1,ne3+1) ws(3,7)];    
    clist = zeros(1,inum+gnum);
    nc1 = length(eta1)-1;
    nc2 = length(eta2)-1;
    nc3 = length(eta3)-1;
    cellnum = [nc1 nc2 nc3];
    c2ilist = zeros(10,nc1*nc2*nc3);
    c2inum = zeros(1,nc1*nc2*nc3);
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
        l = j1 + (j2-1)*nc1 + (j3-1)*nc1*nc2;      
        clist(i) = l;  % cell l
%         clist(1,i) = j1;
%         clist(2,i) = j2;
%         clist(3,i) = j3;
        c2inum(l) = c2inum(l) + 1;
        c2ilist(c2inum(l),l) = i;    
    end    
end


