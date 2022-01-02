function numlattices(rcutmax, pbc, a, b, c)

    d = (length(c) < 3) ? 2 : 3

    m = 0
    n = 0
    p = 0

    if pbc[1] == 1
        m = Int64(round(rcutmax/a[1]))
    end
    if pbc[2] == 1
        n = Int64(round(rcutmax/b[2]))
    end
    if d==3
        if pbc[3] == 1
            p = Int64(round(rcutmax/c[3]))
        end
    end

    return m, n, p
end


function referenceorigins(m, n, p, d)

xref = zeros(d, (2*m+1)*(2*n+1)*(2*p+1))
for i = 1:(2*p+1) # z
    for j = 1:(2*n+1) # y
        for k = 1:(2*m+1) # x
            t = k + (2*m+1)*(j-1) + (2*m+1)*(2*n+1)*(i-1)
            xref[1, t] = k-m-1;  
            xref[2, t] = j-n-1;  
            if d==3 
                xref[3, t] = i-p-1;  
            end
        end
    end
end

return xref

end

function latticeorigins(a, b, c, m, n, p)

    d = (length(c) < 3) ? 2 : 3
 
    xo = referenceorigins(m, n, p, d)

    if d==2
        B2R, R2B = domainmaps(a, b);    
    else
        B2R, R2B = domainmaps(a, b, c);    
    end

    latorigins = R2B*xo; 

    return latorigins
end

function latticecoords(x, rcutmax, pbc, a, b, c)

    m, n, p = numlattices(rcutmax, pbc, a, b, c)

    latorigins = latticeorigins(a, b, c, m, n, p)
    ind = (m+1) + (2*m+1)*(n) + (2*m+1)*(2*n+1)*(p)

    d = size(x,1)
    nx = size(x,2)
    nl = size(latorigins,2)  # number of lattices

    y = zeros(d, nx*nl)
    for j = 1:nx
        for k = 1:d
            y[k,j] = x[k,j] 
        end
    end
    q = nx;
    for i = 1:nl
        if i != ind            
            for j = 1:nx
                q = q + 1
                for k = 1:d
                    y[k,q] = latorigins[k,i] + x[k,j]
                end
            end
        end
    end

    alist = zeros(Int64,nx*nl)
    for i = 1:nl
        for j = 1:nx
            alist[j + nx*(i-1)] = j
        end
    end

    # get 
    neighi, neighj, neighnum = neighborlist(y, rcutmax);

    neighi = neighi[neighi .<= nx];
    ne = length(neighi);

    # list of neighbors for each atom i in x
    neighlist = neighj[1:ne];

    # number of neighbors for each atom i in x
    neighnum = neighnum[1:nx];
    neighnum = [0; cumsum(neighnum)];

    return y, alist, neighlist, neighnum
end

