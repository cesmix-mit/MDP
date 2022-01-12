function getdescriptors(descriptors)

    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   

    n = length(descriptors) 
    emdescriptors = Array{Any}(nothing, n)
    mldescriptors = nothing;

    n = 0
    for i = 1:length(descriptors)    
        if typeof(descriptors[i]) == Potential.PotentialStruct
            n = n + 1
            emdescriptors[n] = descriptors[i]                         
        else
            mldescriptors = descriptors[i]                         
        end
    end

    if n>0
        if length(descriptors) > n
            for i = length(descriptors):-1:(n+1)    
                deleteat!(emdescriptors, i)            
            end   
        end     
    else
        emdescriptors = []
    end

    return emdescriptors, mldescriptors
end

function maximumcutoff(descriptors)

    emdescriptors, mldescriptor = getdescriptors(descriptors)
    
    rcutmax = 0.0
    for i = 1:length(emdescriptors)                
        rcutmax = max(rcutmax, emdescriptors[i].rcut)        
    end
    if mldescriptor !== nothing
        if lowercase(mldescriptor.name) == "snap"
            rcutmax = max(mldescriptor.rcutmax,rcutmax)
        end
    end

    return rcutmax 
end

function emlpotential(x, q, t, alist, neighlist, neighnum, eta, kappa, descriptors)

    empotential, mlpotential = getdescriptors(descriptors)
    etot, e, f, v, vtot = empiricalpotential(x, q, t, alist, neighlist, neighnum, eta, kappa, empotential)

    if mlpotential !== nothing
        eatom, fatom, vatom  = snappotential(x, t, alist, neighlist, neighnum, mlpotential)
        e = e + eatom
        f = f + fatom
        v = v + vatom 
        etot = sum(e);
        vtot = sum(v,dims=2);
    end

    return etot, e, f, v, vtot
end

