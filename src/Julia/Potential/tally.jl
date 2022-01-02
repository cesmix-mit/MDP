function tallysingle(ei, fi, ai)

# energies per atom i
e = accumarray(ai, ei);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fi,1);

# forces
f = zeros(dim,N);
for d = 1:dim
    # tally forces on each dimension
    f[d,:] = reshape(accumarray(ai, fi[d,:]),(1, N)); 
end

return e, f

end

function tallysinglevirial(fij, xij, ai, factor)

    dim = size(fij,1);
    N = size(fij,2);
    v = zeros(6,N);

    vij = xij[1,:].*fij[1,:]
    vi = accumarray(ai, vij);
    v[1,:] = vi;

    vij = xij[2,:].*fij[2,:]
    vi = accumarray(ai, vij);
    v[2,:] = vi;

    vij = xij[1,:].*fij[2,:]
    vi = accumarray(ai, vij);
    v[6,:] = vi;    

    if dim>2
        vij = xij[3,:].*fij[3,:]
        vi = accumarray(ai, vij);
        v[3,:] = vi;

        vij = xij[2,:].*fij[3,:]
        vi = accumarray(ai, vij);
        v[4,:] = vi;

        vij = xij[1,:].*fij[3,:]
        vi = accumarray(ai, vij);
        v[5,:] = vi;
    end

    v = v*factor;

    return v
end

function tallypair(eij, fij, ai, aj)

# energies per atom i
e = accumarray(ai, eij);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fij,1);

# forces
# f = accumarray2d(ai, fij) - accumarray2d(aj, fij);
f = zeros(dim,N);
for d = 1:dim
    # tally forces on each dimension
    f[d,:] = reshape(accumarray(ai, fij[d,:]),(1, N)) - reshape(accumarray(aj, fij[d,:]),(1, N)); 
end

return e, f

end

function tallypairvirial(fij, xij, ai, aj, factor)

    dim = size(fij,1);
    #N = size(fij,2);

    vij = xij[1,:].*fij[1,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij);

    v = zeros(6,length(vi));
    v[1,:] = vi;

    vij = xij[2,:].*fij[2,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij);
    v[2,:] = vi;

    vij = xij[1,:].*fij[2,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij);
    v[6,:] = vi;    

    if dim>2
        vij = xij[3,:].*fij[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij);
        v[3,:] = vi;

        vij = xij[2,:].*fij[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij);
        v[4,:] = vi;

        vij = xij[1,:].*fij[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij);
        v[5,:] = vi;
    end

    v = v*factor;

    return v
end

function tallytriplet(eijk, fij, fik, ai, aj, ak)

# energies per atom i
e = accumarray(ai, eijk);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fij,1);

# forces
# f = accumarray2d(ai, fij) - accumarray2d(aj, fij);
# f = f + accumarray2d(ai, fik) - accumarray2d(ak, fik);
f = zeros(dim,N);
for d = 1:dim
    # tally forces on each dimension
    f[d,:] = reshape(accumarray(ai, fij[d,:]),(1, N)) - reshape(accumarray(aj, fij[d,:]),(1, N)); 
    f[d,:] = reshape(accumarray(ai, fik[d,:]),(1, N)) - reshape(accumarray(ak, fik[d,:]),(1, N)); 
end

return e, f

end

function tallytripletvirial(fij, fik, xij, xik, ai, aj, ak, factor)

    dim = size(fij,1);
    #N = size(fij,2);
    #v = zeros(6,N);

    vij = (xij[1,:].*fij[1,:] + xik[1,:].*fik[1,:])
    vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij)

    v = zeros(6,length(vi));
    v[1,:] = vi;

    vij = (xij[2,:].*fij[2,:] + xik[2,:].*fik[2,:])
    vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij)
    v[2,:] = vi;

    vij = xij[1,:].*fij[2,:] + xik[1,:].*fik[2,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij)
    v[6,:] = vi;    

    if dim>2
        vij = xij[3,:].*fij[3,:] + xik[3,:].*fik[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij)
        v[3,:] = vi;

        vij = xij[2,:].*fij[3,:] + xik[2,:].*fik[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij)
        v[4,:] = vi;

        vij = xij[1,:].*fij[3,:] + xik[1,:].*fik[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij)
        v[5,:] = vi;
    end

    v = v*factor;

    return v
end

function tallyquadlet(eijkl, fij, fik, fil, ai, aj, ak, al)

# energies per atom i
e = accumarray(ai, eijkl);
#etot = sum(e); # potential energy

N = length(e);
dim = size(fij,1);

# forces
# f = accumarray2d(ai, fij) - accumarray2d(aj, fij);
# f = f + accumarray2d(ai, fik) - accumarray2d(ak, fik);
# f = f + accumarray2d(ai, fil) - accumarray2d(al, fil);
f = zeros(dim,N);
for d = 1:dim
    # tally forces on each dimension
    f[d,:] = reshape(accumarray(ai, fij[d,:]),(1, N)) - reshape(accumarray(aj, fij[d,:]),(1, N)); 
    f[d,:] = reshape(accumarray(ai, fik[d,:]),(1, N)) - reshape(accumarray(ak, fik[d,:]),(1, N)); 
    f[d,:] = reshape(accumarray(ai, fil[d,:]),(1, N)) - reshape(accumarray(al, fil[d,:]),(1, N)); 
end

return e, f

end

function tallyquadletvirial(fij, fik, fil, xij, xik, xil, ai, aj, ak, al, factor)

    dim = size(fij,1);
    #N = size(fij,2);
    #v = zeros(6,N);

    vij = xij[1,:].*fij[1,:] + xik[1,:].*fik[1,:] + xil[1,:].*fil[1,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij) + accumarray(al, vij)

    v = zeros(6,length(vi));
    v[1,:] = vi;

    vij = xij[2,:].*fij[2,:] + xik[2,:].*fik[2,:] + xil[2,:].*fil[2,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij) + accumarray(al, vij)
    v[2,:] = vi;

    vij = xij[1,:].*fij[2,:] + xik[1,:].*fik[2,:] + xil[1,:].*fil[2,:]
    vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij) + accumarray(al, vij)
    v[6,:] = vi;    

    if dim>2
        vij = xij[3,:].*fij[3,:] + xik[3,:].*fik[3,:] + xil[3,:].*fil[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij) + accumarray(al, vij)
        v[3,:] = vi;

        vij = xij[2,:].*fij[3,:] + xik[2,:].*fik[3,:] + xil[2,:].*fil[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij) + accumarray(al, vij)
        v[4,:] = vi;

        vij = xij[1,:].*fij[3,:] + xik[1,:].*fik[3,:] + xil[1,:].*fil[3,:]
        vi = accumarray(ai, vij) + accumarray(aj, vij) + accumarray(ak, vij) + accumarray(al, vij)
        v[5,:] = vi;
    end

    v = v*factor;

    return v
end



















