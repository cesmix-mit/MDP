function empiricaldescriptors(x, q, t, a, b, c, pbc, rcutmax, eta, kappa, emdescriptors)

    # compute full neighbor list
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    if !isempty(q)
        q = q[:,alist];
    end
    
    dim, N = size(x);
    d = []
    dd = []    
    dv = []
    for i = 1:length(emdescriptors)        
        
        style, bondtype, bondedatomtypes, potentialfunc, mu, rcut = getpotential(emdescriptors[i]);    
        potentialfunc = getfield(Main, Symbol(potentialfunc));
        
        # single potentials
        if (style == "single") & (bondtype == "nonbonded")            
            ilist = Array(1:N);                        
         elseif (style == "single") & (bondtype == "bonded")     
            typei = bondedatomtypes[1];
            ilist = findatomtype(Array(1:N), t, typei);      
        end
        if (style == "single")
            xi, qi, ai, ti = neighsingles(x, q, t, ilist);        
            ei, fi = potentialfunc(xi, qi, ti, mu, eta, kappa);
            ea, fa = tallysingle(ei, fi, ai);   
            d = [d; sum(ea)];         
            if length(dd)==0
                dd = fa;
            else
                dd = cat(dd, fa, dims=3)
            end    
            va = tallysinglevirial(fi, xi, ai, -1.0)
            if length(dv)==0
                dv = sum(va, dims=2)
            else
                dv = [dv sum(va, dims=2)]
            end            
        end
    
        # pair potentials
        if (style == "pair") & (bondtype == "nonbonded")            
            ilist = Int64.(Array(1:N));                        
            pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);
         elseif (style == "pair") & (bondtype == "bonded")     
            typei = bondedatomtypes[1];
            typej = bondedatomtypes[2];            
    
            ilist = findatomtype(Array(1:N), t, typei);      
            pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, typej, neighlist, neighnum, rcut*rcut);                                                
        end
        if (style == "pair")
            xij, ai, aj, ti, tj, qi, qj = neighpairs(y, q, pairlist, pairnum, t, ilist, alist);            
            eij, fij = potentialfunc(xij, qi, qj, ti, tj, mu, eta, kappa);     
            # tm0 = 1.0*fij   
            # dim, na = size(xij)
            # for i1 = 1:dim
            #     for j1 = 1:na
            #         tm1 = 1.0*xij
            #         tm1[i1,j1] = tm1[i1,j1] - 1e-6
            #         tm2, ~ = potentialfunc(tm1, qi, qj, ti, tj, mu, eta, kappa);     
            #         tm0[i1,j1,:] = (eij[j1,:] - tm2[j1,:])/1e-6
            #     end
            # end
            # display(abs.(fij-tm0))
            # error("here")
            for m = 1:size(eij,2)
                ea, fa = tallypair(eij[:,m], fij[:,:,m], ai, aj);            
                d = [d; sum(ea)];
                if length(dd)==0
                    dd = fa;
                else
                    dd = cat(dd, fa, dims=3)
                end    
                va = tallypairvirial(fij[:,:,m], xij, ai, aj, -1.0/2)
                if length(dv)==0
                    dv = sum(va, dims=2)
                else
                    dv = [dv sum(va, dims=2)]
                end                
            end
        end
    
        # triplet potentials
        if (style == "triplet") & (bondtype == "nonbonded")            
            ilist = Array(1:N);                     
            pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);        
         elseif (style == "triplet") & (bondtype == "bonded")     
            typei = bondedatomtypes[1];
            typej = bondedatomtypes[2];            
            typek = bondedatomtypes[3];            
    
            ilist = findatomtype(Array(1:N), t, typei);      
            pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej; typek], neighlist, neighnum, rcut*rcut);                                                
        end
        if (style == "triplet")
            tripletlist, tripletnum = neightripletlist(pairlist, pairnum, ilist);
            xij, xik, ai, aj, ak, ti, tj, tk, qi, qj, qk = neightriplets(y, q, tripletlist, tripletnum, t, ilist, alist);
            eij, fij, fik = potentialfunc(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);                                
            # tmj = 1.0*fik  
            # tmk = 1.0*fik   
            # dim, na = size(xij)
            # na = 10;
            # for i1 = 1:dim
            #     for j1 = 1:na
            #         tm1 = 1.0*xij
            #         tm1[i1,j1] = tm1[i1,j1] - 1e-6
            #         tm2, ~ = potentialfunc(tm1, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);                                
            #         tmj[i1,j1,:] = (eij[j1,:] - tm2[j1,:])/1e-6
            #         tm1 = 1.0*xik
            #         tm1[i1,j1] = tm1[i1,j1] - 1e-6
            #         tm2, ~ = potentialfunc(xij, tm1, qi, qj, qk, ti, tj, tk, mu, eta, kappa);                                
            #         tmk[i1,j1,:] = (eij[j1,:] - tm2[j1,:])/1e-6
            #     end
            # end
            # display(abs.(fij[:,1:na,:]-tmj[:,1:na,:]))
            # display(abs.(fik[:,1:na,:]-tmk[:,1:na,:]))
            # error("here")
            for m = 1:size(eij,2)
                ea, fa = tallytriplet(eij[:,m], fij[:,:,m], fik[:,:,m], ai, aj, ak);                                             
                d = [d; sum(ea)];
                if length(dd)==0
                    dd = fa;
                else
                    dd = cat(dd, fa, dims=3)
                end    
                va = tallytripletvirial(fij[:,:,m], fik[:,:,m], xij, xik, ai, aj, ak, -1.0/3)
                if length(dv)==0
                    dv = sum(va, dims=2)
                else
                    dv = [dv sum(va, dims=2)]
                end                
            end            
        end
        
        # quadlet potentials
        if (style == "quadlet") & (bondtype == "nonbonded")            
            ilist = Array(1:N);                        
            pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcut*rcut);        
         elseif (style == "quadlet") & (bondtype == "bonded")     
            typei = bondedatomtypes[1];
            typej = bondedatomtypes[2];            
            typek = bondedatomtypes[3];      
            typel = bondedatomtypes[4];      
    
            ilist = findatomtype(Array(1:N), t, typei);      
            pairlist, pairnum = neighpairlisttypej(y, ilist, alist, t, [typej typek typel], neighlist, neighnum, rcut*rcut);                                                
        end
        if (style == "quadlet")
            quadletlist, quadletnum = neighquadletlist(pairlist, pairnumsum, ilist);        
            xij, xik, xil, ai, aj, ak, al, ti, tj, tk, tl, qi, qj, qk, ql = 
                        neighquadlets(x, q, quadletlist, quadletnum, t, ilist, alist);
            eij, fij, fik, fil = potentialfunc(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa);            
            for m = 1:size(eij,2)
                ea, fa = tallyquadlet(eij[:,m], fij[:,:,m], fik[:,:,m], fil[:,:,m], ai, aj, ak, al);                            
                d = [d; sum(ea)];
                if length(dd)==0
                    dd = fa;
                else
                    dd = cat(dd, fa, dims=3)
                end    
                va = tallyquadletvirial(fij[:,:,m], fik[:,:,m], fil[:,:,m], xij, xik, xil, ai, aj, ak, al, -1.0/4)     
                if length(dv)==0
                    dv = sum(va, dims=2)
                else
                    dv = [dv sum(va, dims=2)]
                end                
            end
        end
            
    end    
        
    return d, dd, dv
    
end
    
    
    
    
    
    