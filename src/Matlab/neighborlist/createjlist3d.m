function jlist = createjlist3d(xi, ilist, clist, c2ilist, c2inum, r, inum)

jlist = cell(inum,1);
rsq = r*r;
for i = 1:inum    
    jlist{i} = [];
    xt = xi(:,i);
    j1 = clist(1,i);
    j2 = clist(2,i);
    j3 = clist(3,i);
    for i1 = -1:1:1
        k1 = j1 + i1;
        for i2 = -1:1:1
            k2 = j2 + i2;
            for i3 = -1:1:1
                k3 = j3 + i3;            
%                 [j1 j2 j3]
%                 [k1 k2 k3]                
                m = c2inum(k1,k2,k3);            
                for l = 1:m
                    j = c2ilist(l,k1,k2,k3);
                    rij = (xt(1)-xi(1,j)).^2 + (xt(2)-xi(2,j)).^2 + (xt(3)-xi(3,j)).^2;
                    if rij<=rsq
                        jlist{i} = [jlist{i} j];
                    end
                end
            end
        end
    end
    jlist{i} = setdiff(jlist{i},i);
end

