function jlist = createjlist2d(xi, ilist, clist, c2ilist, c2inum, r, inum)

jlist = cell(inum,1);
rsq = r*r;
for i = 1:inum    % loop over each i atom 
    jlist{i} = [];
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
                rij = (xt(1)-xi(1,j)).^2 + (xt(2)-xi(2,j)).^2;
                if rij<=rsq
                    jlist{i} = [jlist{i} ilist(j)];
                end
            end
        end
    end
    jlist{i} = setdiff(jlist{i},i);
end

