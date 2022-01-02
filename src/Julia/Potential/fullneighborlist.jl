function fullneighborlist(x, a, b, c, pbc, rcutmax)

latticemethod = 0    
if (pbc[1] == 1) & (rcutmax >= 1.0*a[1])       
    latticemethod = 1
end
if (pbc[2] == 1) & (rcutmax >= 1.0*b[2])       
    latticemethod = 1
end
if (length(c)>2) & (length(pbc) > 2)
    if (pbc[3] == 1) & (rcutmax >= 1.0*c[3])       
        latticemethod = 1
    end
end

if (latticemethod == 1)
    y, alist, neighlist, neighnum = latticecoords(x, rcutmax, pbc, a, b, c)
    return y, alist, neighlist, neighnum
end

# offset simulation box for periodic boundary conditions
boxoffset = [rcutmax; rcutmax; rcutmax]; 

# bounding box
B2R, R2B, v, w, vr, wr, pimages = prebox(pbc[:], boxoffset, a[:], b[:], c[:]);

y, alist = atomlist(x, pimages, wr, B2R);

neighi, neighj, neighnum = neighborlist(y, rcutmax);

n = size(x,2);
neighi = neighi[neighi .<= n];
ne = length(neighi);
neighlist = neighj[1:ne];
neighnum = neighnum[1:n];
neighnum = [0; cumsum(neighnum)];

return y, alist, neighlist, neighnum

end

