function jlist = createvlist(xi, atomtype, ilist, clist, c2ilist, c2inum, ellipsoid, inum, dim)

jlist = cell(inum,1);

if dim==2
    for i = 1:inum  % loop over each i atom 
        jlist{i} = [];
        t = atomtype(i); % element of atom i 
        A = ellipsoid(:,:,t); % ellipsoid for element t        
        xt = xi(:,i);
        j1 = clist(1,i);  % cell (j1, j2) contains atom i
        j2 = clist(2,i);
        for i1 = -1:1:1   % loop over neighboring cells of cell (j1, j2)
            k1 = j1 + i1;
            for i2 = -1:1:1
                k2 = j2 + i2;
                % cell (k1, k2) is a neighbor of cell (j1, j2)
                m = c2inum(k1,k2);            
                for l = 1:m  % loop over each atom j of cell (k1, k2)
                    j = c2ilist(l,k1,k2); % atom j
                    % distance between atom i and atom j 
                    rij = (xt-xi(:,j))'*(A*(xt-xi(:,j)));
                    if rij<=1
                        % add atom j into the list
                        jlist{i} = [jlist{i} j];
                    end
                end
            end
        end
        jlist{i} = setdiff(jlist{i},i);
    end
else
    for i = 1:inum    
        jlist{i} = [];
        t = atomtype(i); % element of atom i 
        A = ellipsoid(:,t); % ellipsoid for element t        
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
                    m = c2inum(k1,k2,k3);            
                    for l = 1:m
                        j = c2ilist(l,k1,k2,k3);
                        rij = (xt-xi(:,j))'*(A*(xt-xi(:,j)));
                        if rij<=1
                            jlist{i} = [jlist{i} j];
                        end
                    end
                end
            end
        end
        jlist{i} = setdiff(jlist{i},i);
    end
end


