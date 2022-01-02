using Preprocessing, Potential


function linearsolve(Ae, Af, As, be, bf, bs)

    cefs = zeros(0)
    if (length(Ae)>0) & (length(Af)>0) & (length(As)>0) & (length(be)>0) & (length(bf)>0) & (length(bs)>0) 
        cefs = cat(cat(Ae, Af, dims=1), As, dims=1)\[be[:]; bf[:]; bs[:]]
    end

    cef = zeros(0)
    if (length(Ae)>0) & (length(Af)>0) & (length(be)>0) & (length(bf)>0) 
        cef = cat(Ae, Af, dims=1)\[be[:]; bf[:]]                
    end

    ce = zeros(0)
    if (length(Ae)>0) & (length(be)>0)
        ce = Ae\be
    end

    cf = zeros(0)
    if (length(Af)>0) & (length(bf)>0)
        cf = Af\bf                
    end

    return ce, cf, cef, cefs
end

function lsqsolve(Ae, Af, As, be, bf, bs)

    cefs = zeros(0)
    cef = zeros(0)
    ce = zeros(0)
    cf = zeros(0)

    sumAe = sum(abs.(Ae[:])) > 1e-15
    sumAf = sum(abs.(Af[:])) > 1e-15
    sumAs = sum(abs.(As[:])) > 1e-15
    sumbe = sum(abs.(be[:])) > 1e-15
    sumbf = sum(abs.(bf[:])) > 1e-15
    sumbs = sum(abs.(bs[:])) > 1e-15

    if (sumAe) & (sumAf) & (sumAs) & (sumbe) & (sumbf) & (sumbs) 
        A = Ae+Af+As
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end
        cefs = (A)\(be+bf+bs)
    end
    if (sumAe) & (sumAf) & (sumbe) & (sumbf) 
        A = Ae+Af
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end
        cef = A\(be+bf)
    end
    if (sumAe) & (sumbe) 
        A = 1.0*Ae
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end
        ce = A\be
    end
    if (sumAf) & (sumbf) 
        A = 1.0*Af
        for i = 1:size(A,1)
            A[i,i] = A[i,i]*(1 + 1e-15)
        end        
        if abs(A[1,1])<=1e-10
            A[1,1] = 1e-10
        end
        cf = A\bf
    end

    return ce, cf, cef, cefs
end

function insertarray(ce, cei, i)
    if (length(cei)>0)
        ce[:,i] = cei
    end
    return ce
end

function catarray(Ae, be, Aei, bei)
    if (length(Aei) > 0) & (length(bei) > 0)
        Ae = cat(Ae, Aei, dims=1)
        be = [be; bei]
    end
    return Ae, be
end

function calmae(Aei, bei, cei)
    mae = -1.0
    if (length(bei) > 0) & (length(cei)>0) & (length(Aei)>0)       
        mae = sum(abs.(bei - Aei*cei))/length(bei)       
    end
    return mae
end

function calerrors(Aei, bei, cei)
    mae = -1.0
    rmse = - 1.0
    rsq = -1.0
    ssr = -1.0
    avg = -1.0
    if (length(bei) > 0) & (length(cei)>0) & (length(Aei)>0)       
        res = bei - Aei*cei
        mae = sum(abs.(res))/length(bei)       
        ssr = sum(res.^2)
        mse = ssr /length(bei);
        rmse = sqrt(mse) 
        avg = sum(bei)/length(bei);
        rsq = 1.0 - ssr/sum((bei .- avg).^2);
    end
    return mae, rmse, rsq, ssr, avg
end

function applyweight(Aei, bei, we)
    if (length(bei) > 0) & (length(Aei) > 0) & (length(we) > 0) 
        Aei = we[1]*Aei
        bei = we[1]*bei
    end
    return Aei, bei 
end

function applytranspose(Aei, bei, m, n)    
    Cei = zeros(m,n)
    if (length(bei) > 0) & (length(Aei) > 0) 
        Cei = Aei'*bei        
    end
    return Cei
end

function linearfit(data, emdescriptors, mldescriptors, eta, kappa, pbc=nothing, a=nothing, b=nothing, c=nothing, method=nothing, normalizeenergy=nothing)
    
    if normalizeenergy === nothing
        normalizeenergy = 0 # divide energy by number of atoms
    end
    if method === nothing
        method = "lsq"
    end

    n = length(data)
    emae = -ones(Float64,n,4)
    fmae = -ones(Float64,n,4)
    smae = -ones(Float64,n,4)
    ce = 1.0
    cf = 1.0
    cef = 1.0
    cefs = 1.0
    Ae = 1.0
    Af = 1.0
    As = 1.0
    be = 1.0
    bf = 1.0
    bs = 1.0
    De = 1.0
    Df = 1.0
    Ds = 1.0
    de = 1.0
    df = 1.0
    ds = 1.0
    Le = 1.0
    Lf = 1.0
    Ls = 1.0
    he = 1.0
    hf = 1.0
    hs = 1.0
    for i = 1:n
        # training data set
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        training = Array(1:config.nconfigs)      

        Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.emldescriptors(config, training, eta, kappa, normalizeenergy, emdescriptors, mldescriptors)
        
        m = size(Aei,2) # total number of descriptors
        if i == 1            
            ce = ce*zeros(Float64,m,n)
            cf = cf*zeros(Float64,m,n)
            cef = cef*zeros(Float64,m,n)
            cefs = cefs*zeros(Float64,m,n)
            De = De*zeros(m,m,n)
            Df = Df*zeros(m,m,n)
            Ds = Ds*zeros(m,m,n)
            de = de*zeros(m,n)
            df = df*zeros(m,n)
            ds = ds*zeros(m,n)
            Le = Le*zeros(m,m)
            Lf = Lf*zeros(m,m)
            Ls = Ls*zeros(m,m)
            he = he*zeros(m)
            hf = hf*zeros(m)
            hs = hs*zeros(m)
        end
        
        De[:,:,i] = applytranspose(Aei, Aei, m, m)    
        Df[:,:,i] = applytranspose(Afi, Afi, m, m)    
        Ds[:,:,i] = applytranspose(Asi, Asi, m, m)    
        de[:,i] = applytranspose(Aei, bei, m, 1)    
        df[:,i] = applytranspose(Afi, bfi, m, 1)    
        ds[:,i] = applytranspose(Asi, bsi, m, 1)    

        Aei, bei = applyweight(Aei, bei, config.we)
        Afi, bfi = applyweight(Afi, bfi, config.wf)
        Asi, bsi = applyweight(Asi, bsi, config.ws)
        
        Le += applytranspose(Aei, Aei, m, m)    
        Lf += applytranspose(Afi, Afi, m, m)    
        Ls += applytranspose(Asi, Asi, m, m)    
        he += applytranspose(Aei, bei, m, 1)    
        hf += applytranspose(Afi, bfi, m, 1)    
        hs += applytranspose(Asi, bsi, m, 1)    

        if length(data)==1
            if method == "lsq"
                cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
            else
                cei, cfi, cefi, cefsi = linearsolve(Aei, Afi, Asi, bei, bfi, bsi)                        
            end            
            ce = insertarray(ce, cei, i)
            cf = insertarray(cf, cfi, i)
            cef = insertarray(cef, cefi, i)
            cefs = insertarray(cefs, cefsi, i)

            emae[i,1] = calmae(Aei, bei, cei)
            emae[i,2] = calmae(Aei, bei, cfi)
            emae[i,3] = calmae(Aei, bei, cefi)
            emae[i,4] = calmae(Aei, bei, cefsi)
            fmae[i,1] = calmae(Afi, bfi, cei)
            fmae[i,2] = calmae(Afi, bfi, cfi)
            fmae[i,3] = calmae(Afi, bfi, cefi)
            fmae[i,4] = calmae(Afi, bfi, cefsi)
            smae[i,1] = calmae(Asi, bsi, cei)
            smae[i,2] = calmae(Asi, bsi, cfi)
            smae[i,3] = calmae(Asi, bsi, cefi)
            smae[i,4] = calmae(Asi, bsi, cefsi)
            
            return ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
        end

        if (i == 1) | (method == "lsq")
            Ae = Aei
            Af = Afi
            As = Asi
            be = bei
            bf = bfi
            bs = bsi
        else         
            Ae, be = catarray(Ae, be, Aei, bei)
            Af, bf = catarray(Af, bf, Afi, bfi)
            As, bs = catarray(As, bs, Asi, bsi)
        end
        
        if (method == "lsq")
            cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
        else
            cei, cfi, cefi, cefsi = linearsolve(Ae, Af, As, be, bf, bs)
        end

        ce = insertarray(ce, cei, i)
        cf = insertarray(cf, cfi, i)
        cef = insertarray(cef, cefi, i)
        cefs = insertarray(cefs, cefsi, i)

        emae[i,1] = calmae(Ae, be, cei)
        emae[i,2] = calmae(Ae, be, cfi)
        emae[i,3] = calmae(Ae, be, cefi)
        emae[i,4] = calmae(Ae, be, cefsi)
        fmae[i,1] = calmae(Af, bf, cei)
        fmae[i,2] = calmae(Af, bf, cfi)
        fmae[i,3] = calmae(Af, bf, cefi)
        fmae[i,4] = calmae(Af, bf, cefsi)
        smae[i,1] = calmae(As, bs, cei)
        smae[i,2] = calmae(As, bs, cfi)
        smae[i,3] = calmae(As, bs, cefi)
        smae[i,4] = calmae(As, bs, cefsi)
    end

    return ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
end

function calmean(emae, szbe)
    n,p,m = size(emae)
    em = -ones(p,m)       
    if (sum(szbe) > 0) & (sum(emae[:]) > 0)
        for k = 1:p
            for j = 1:m
                em[k,j] = 0.0
                for i = 1:n            
                    em[k,j] = em[k,j] + emae[i,k,j]*szbe[i]
                end
                em[k,j] =  em[k,j]/sum(szbe)
            end
        end
    end
    return em
end

function validation(data, emdescriptors, mldescriptors, ce, cf, cef, cefs, eta, kappa, pbc=nothing, a=nothing, b=nothing, c=nothing)    

    normalizeenergy = 0
    m = size(ce,2)
    n = length(data)
    emae = -ones(Float64,n,4,m)
    fmae = -ones(Float64,n,4,m)
    smae = -ones(Float64,n,4,m)
    szbe = zeros(Int64, n)
    szbf = zeros(Int64, n)
    szbs = zeros(Int64, n)
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        valid = Array(1:config.nconfigs)      

        Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.emldescriptors(config, valid, eta, kappa, normalizeenergy, emdescriptors, mldescriptors)
        szbe[i] = length(bei)        
        szbf[i] = length(bfi)        
        szbs[i] = length(bsi)        

        for j = 1:m
            cei = ce[:,j]
            cfi = cf[:,j]
            cefi = cef[:,j]
            cefsi = cefs[:,j] 
            emae[i,1,j] = calmae(Aei, bei, cei)
            emae[i,2,j] = calmae(Aei, bei, cfi)
            emae[i,3,j] = calmae(Aei, bei, cefi)
            emae[i,4,j] = calmae(Aei, bei, cefsi)
            fmae[i,1,j] = calmae(Afi, bfi, cei)
            fmae[i,2,j] = calmae(Afi, bfi, cfi)
            fmae[i,3,j] = calmae(Afi, bfi, cefi)
            fmae[i,4,j] = calmae(Afi, bfi, cefsi)
            smae[i,1,j] = calmae(Asi, bsi, cei)
            smae[i,2,j] = calmae(Asi, bsi, cfi)
            smae[i,3,j] = calmae(Asi, bsi, cefi)
            smae[i,4,j] = calmae(Asi, bsi, cefsi)
        end
    end

    eme = calmean(emae, szbe)
    fme = calmean(fmae, szbf)
    sme = calmean(smae, szbs)

    return eme, fme, sme, emae, fmae, smae, szbe, szbf, szbs
end

function linearfit2(data, descriptors, potential, optim, pbc=nothing, a=nothing, b=nothing, c=nothing)
        
    data = deletenothing(data)
    descriptors = deletenothing(descriptors)
    potential = deletenothing(potential)
    
    if pbc === nothing        
        pbc = [1; 1; 1] # periodic boundary conditions
    end

    lossfunc, method, normalizeenergy, normalizestress, eta, kappa, 
        weightouter, etaspace, kappalist = getoptim(optim)

    n = length(data)
    emae = -ones(Float64,n,4)
    fmae = -ones(Float64,n,4)
    smae = -ones(Float64,n,4)
    ce = 1.0
    cf = 1.0
    cef = 1.0
    cefs = 1.0
    Ae = 1.0
    Af = 1.0
    As = 1.0
    be = 1.0
    bf = 1.0
    bs = 1.0
    De = 1.0
    Df = 1.0
    Ds = 1.0
    de = 1.0
    df = 1.0
    ds = 1.0
    Le = 1.0
    Lf = 1.0
    Ls = 1.0
    he = 1.0
    hf = 1.0
    hs = 1.0
    for i = 1:n
        # training data set
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        training = Array(1:config.nconfigs)      

        Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.emldescriptors2(config, training, eta, kappa, 
                normalizeenergy, normalizestress, descriptors, potential)

        m = size(Aei,2) # total number of descriptors
        if i == 1            
            ce = ce*zeros(Float64,m,n)
            cf = cf*zeros(Float64,m,n)
            cef = cef*zeros(Float64,m,n)
            cefs = cefs*zeros(Float64,m,n)
            De = De*zeros(m,m,n)
            Df = Df*zeros(m,m,n)
            Ds = Ds*zeros(m,m,n)
            de = de*zeros(m,n)
            df = df*zeros(m,n)
            ds = ds*zeros(m,n)
            Le = Le*zeros(m,m)
            Lf = Lf*zeros(m,m)
            Ls = Ls*zeros(m,m)
            he = he*zeros(m)
            hf = hf*zeros(m)
            hs = hs*zeros(m)
        end
        
        De[:,:,i] = applytranspose(Aei, Aei, m, m)    
        Df[:,:,i] = applytranspose(Afi, Afi, m, m)    
        Ds[:,:,i] = applytranspose(Asi, Asi, m, m)    
        de[:,i] = applytranspose(Aei, bei, m, 1)    
        df[:,i] = applytranspose(Afi, bfi, m, 1)    
        ds[:,i] = applytranspose(Asi, bsi, m, 1)    

        Aei, bei = applyweight(Aei, bei, config.we)
        Afi, bfi = applyweight(Afi, bfi, config.wf)
        Asi, bsi = applyweight(Asi, bsi, config.ws)
                
        Le += applytranspose(Aei, Aei, m, m)    
        Lf += applytranspose(Afi, Afi, m, m)    
        Ls += applytranspose(Asi, Asi, m, m)    
        he += applytranspose(Aei, bei, m, 1)    
        hf += applytranspose(Afi, bfi, m, 1)    
        hs += applytranspose(Asi, bsi, m, 1)            
        
        if length(data)==1
            if optim.method == "lsq"
                cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
            else
                cei, cfi, cefi, cefsi = linearsolve(Aei, Afi, Asi, bei, bfi, bsi)                        
            end            
            ce = insertarray(ce, cei, i)
            cf = insertarray(cf, cfi, i)
            cef = insertarray(cef, cefi, i)
            cefs = insertarray(cefs, cefsi, i)

            emae[i,1] = calmae(Aei, bei, cei)
            emae[i,2] = calmae(Aei, bei, cfi)
            emae[i,3] = calmae(Aei, bei, cefi)
            emae[i,4] = calmae(Aei, bei, cefsi)
            fmae[i,1] = calmae(Afi, bfi, cei)
            fmae[i,2] = calmae(Afi, bfi, cfi)
            fmae[i,3] = calmae(Afi, bfi, cefi)
            fmae[i,4] = calmae(Afi, bfi, cefsi)
            smae[i,1] = calmae(Asi, bsi, cei)
            smae[i,2] = calmae(Asi, bsi, cfi)
            smae[i,3] = calmae(Asi, bsi, cefi)
            smae[i,4] = calmae(Asi, bsi, cefsi)
            
            return ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
        end

        if (i == 1) | (optim.method == "lsq")
            Ae = 1.0*Aei
            Af = 1.0*Afi
            As = 1.0*Asi
            be = 1.0*bei
            bf = 1.0*bfi
            bs = 1.0*bsi
        else         
            Ae, be = catarray(Ae, be, Aei, bei)
            Af, bf = catarray(Af, bf, Afi, bfi)
            As, bs = catarray(As, bs, Asi, bsi)                        
        end                        

        if (optim.method == "lsq")
            cei, cfi, cefi, cefsi = lsqsolve(Le, Lf, Ls, he, hf, hs)
        else
            cei, cfi, cefi, cefsi = linearsolve(Ae, Af, As, be, bf, bs)
        end

        # if i==n
        #     display(i)
        #     A = Le+Lf+Ls
        #     b = he+hf+hs                
        #     fileID = open("Ajl","w"); write(fileID,A[:]); close(fileID);
        #     fileID = open("bjl","w"); write(fileID,b[:]); close(fileID);
        #     display(A)    
        #     display(b)    
        #     coeff = (0.5*(A+A'))\b
        #     display(coeff)
        #     # display(size(Ae))
        #     # display(size(Af))
        #     # display(size(As))
        #     coeff = cat(cat(Ae, Af, dims=1), As, dims=1)\[be[:]; bf[:]; bs[:]]
        #     # A = Ae+Af+As
        #     # b = be+bf+bs                        
        #     # coeff = A\b
        #     display(coeff)
        #     error("here")
        # end

        ce = insertarray(ce, cei, i)
        cf = insertarray(cf, cfi, i)
        cef = insertarray(cef, cefi, i)
        cefs = insertarray(cefs, cefsi, i)

        emae[i,1] = calmae(Ae, be, cei)
        emae[i,2] = calmae(Ae, be, cfi)
        emae[i,3] = calmae(Ae, be, cefi)
        emae[i,4] = calmae(Ae, be, cefsi)
        fmae[i,1] = calmae(Af, bf, cei)
        fmae[i,2] = calmae(Af, bf, cfi)
        fmae[i,3] = calmae(Af, bf, cefi)
        fmae[i,4] = calmae(Af, bf, cefsi)
        smae[i,1] = calmae(As, bs, cei)
        smae[i,2] = calmae(As, bs, cfi)
        smae[i,3] = calmae(As, bs, cefi)
        smae[i,4] = calmae(As, bs, cefsi)
    end

    if lossfunc == "energy" 
        coeff = ce[:,end]
    elseif lossfunc == "force" 
        coeff = cf[:,end]
    elseif lossfunc == "energyforce" 
        coeff = cef[:,end]
    elseif lossfunc == "energyforcestress" 
        coeff = cefs[:,end]
    end    
    
    return coeff, ce, cf, cef, cefs, emae, fmae, smae, De, Df, Ds, de, df, ds
end


# function calmae(emae, szbe)
#     n,p,m = size(emae)
#     em = -ones(p,m)       
#     if (sum(szbe) > 0) & (sum(emae[:]) > 0)
#         for k = 1:p
#             for j = 1:m
#                 em[k,j] = 0.0
#                 for i = 1:n            
#                     em[k,j] = em[k,j] + emae[i,k,j]*szbe[i]
#                 end
#                 em[k,j] =  em[k,j]/sum(szbe)
#             end
#         end
#     end
#     return em
# end
