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

function setcutoff(emdescriptors, mldescriptors, eta)
    if length(eta)>0
        for i = 1:length(emdescriptors)    
            if i==1
                emdescriptors[i].rcut = eta[1] # long-range empirical potentials
            else
                emdescriptors[i].rcut = eta[1] # short-range empirical potentials
            end
        end
        if mldescriptors !== nothing
            if lowercase(mldescriptors.name) == "snap"
                if length(emdescriptors) == 0    
                    mldescriptors.rcutfac = eta[1]
                else
                    mldescriptors.rcutfac = eta[1]                
                end
            end
        end
    end
    return emdescriptors, mldescriptors
end

function maespace(traindata, validdata, emdescriptors, mldescriptors, etaspace, kappa, pbc=nothing, a=nothing, b=nothing, c=nothing, method=nothing, normalizeenergy=nothing)

    etamin = etaspace[:,1]
    etamax = etaspace[:,2]
    N = Int64.(etaspace[:,3])
    eta = tensornodes(etamin, etamax, N)
    
    # read data 
    trainconfig = Array{Any}(undef, length(traindata))
    for i = 1:length(traindata)
        trainconfig[i], indices = Preprocessing.readconfigdata(traindata[i], pbc, a, b, c)                  
    end
    validconfig = Array{Any}(undef, length(validdata))
    for i = 1:length(validdata)
        validconfig[i], indices = Preprocessing.readconfigdata(validdata[i], pbc, a, b, c)                  
    end

    neta = size(eta,1)
    eme = -1.0
    fme = -1.0
    sme = -1.0
    for n = 1:neta
        emdescriptors, mldescriptors = setcutoff(emdescriptors, mldescriptors, eta[n,:])
        ce, cf, cef, cefs = linearfit(trainconfig, emdescriptors, mldescriptors, eta[n,:], kappa, pbc, a, b, c, method, normalizeenergy)
        emae, fmae, smae = validation(validconfig, emdescriptors, mldescriptors, ce, cf, cef, cefs, eta[n,:], kappa, pbc, a, b, c)        

        if n == 1
            sz = size(emae)
            eme = -ones(neta, sz[1], sz[2])
            fme = -ones(neta, sz[1], sz[2])
            sme = -ones(neta, sz[1], sz[2])            
        end
        eme[n,:,:] = emae
        fme[n,:,:] = fmae
        sme[n,:,:] = smae
        
        emae = round.(emae; digits=9)
        fmae = round.(fmae; digits=9)
        smae = round.(smae; digits=9)
        print("            Energy Fitting      Force Fitting      Energy+Force Fitting      Energy+Force+Stress Fitting \n");
        print("Energy MAE:   " * string(emae[1]) * "         " * string(emae[2]) * "         " * string(emae[3]) * "           " * string(emae[4]) * "\n");
        print("Force MAE:    " * string(fmae[1]) * "         " * string(fmae[2]) * "         " * string(fmae[3]) * "           " * string(fmae[4]) * "\n");
        print("Stress MAE:   " * string(smae[1]) * "         " * string(smae[2]) * "         " * string(smae[3]) * "           " * string(smae[4]) * "\n");        
    end

    return eta, eme, fme, sme
end

function optimize(traindata, testdata, descriptors, etaspace, kappa, weightouter, lossfunc, pbc=nothing, a=nothing, b=nothing, c=nothing, method=nothing, normalizeenergy=nothing)

    # getdescriptors
    emdescriptors, mldescriptors = getdescriptors(descriptors)

    # compute MAE on the eta space grid
    etapts, emae, fmae, smae = maespace(traindata, testdata, emdescriptors, mldescriptors, 
                etaspace, kappa, pbc, a, b, c, method, normalizeenergy)    
    
    # compute the loss function on the eta space grid
    if lossfunc == "energy" 
        f = weightouter[1]*emae[:,1,end] 
    elseif lossfunc == "force" 
        f = weightouter[2]*emae[:,2,end] 
    elseif lossfunc == "energyforce" 
        f = weightouter[1]*emae[:,3,end] + weightouter[2]*fmae[:,3,end]
    elseif lossfunc == "energyforcestress" 
        f = weightouter[1]*emae[:,4,end] + weightouter[2]*fmae[:,4,end] + weightouter[3]*smae[:,4,end]
    else
        error("lossfunc is not valid. It must be energy, force, energyforce, or energyforcestress")
    end    

    # interpolate the loss function on the eta space grid
    pc = mpolyinterp(etapts, f, etaspace)

    # optimize the loss function on the eta space
    eta, fmin, iter  = gradientdescent(pc, etaspace)

    # set cut-off radius 
    emdescriptors, mldescriptors = setcutoff(emdescriptors, mldescriptors, eta)

    # compute linear coefficients
    ce, cf, cef, cefs = linearfit(traindata, emdescriptors, mldescriptors, eta, kappa, pbc, a, b, c, method, normalizeenergy)

    if lossfunc == "energy" 
        coeff = ce
    elseif lossfunc == "force" 
        coeff = cf
    elseif lossfunc == "energyforce" 
        coeff = cef
    elseif lossfunc == "energyforcestress" 
        coeff = cefs
    end    

    return eta, coeff, fmin, iter, pc, f, etapts, emae, fmae, smae
end

function validate(data, descriptors, coeff, eta, kappa, pbc=nothing, a=nothing, b=nothing, c=nothing)    

    # getdescriptors
    emdescriptors, mldescriptors = getdescriptors(descriptors)

    # set cut-off radius 
    # if eta !== nothing
    #     emdescriptors, mldescriptors = setcutoff(emdescriptors, mldescriptors, eta)
    # end

    normalizeenergy = 0
    n = length(data)
    emae = -ones(Float64,n)
    fmae = -ones(Float64,n)
    smae = -ones(Float64,n)
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

        emae[i] = calmae(Aei, bei, coeff)
        fmae[i] = calmae(Afi, bfi, coeff)
        smae[i] = calmae(Asi, bsi, coeff)
    end

    eme = sum(emae.*szbe)/sum(szbe)
    fme = sum(fmae.*szbf)/sum(szbf)
    sme = sum(smae.*szbs)/sum(szbs)

    return eme, fme, sme, emae, fmae, smae, szbe, szbf, szbs
end

function setcutoff2(descriptors, eta)
    if length(eta)>0
        for i = 1:length(descriptors)    
            if typeof(descriptors[i]) == Potential.PotentialStruct
                descriptors[i].rcut = eta[1] 
            elseif typeof(descriptors[i]) == Potential.SnaStruct
                descriptors[i].rcutfac = eta[1]
                descriptors[3].rcutmax = eta[1]
            end            
        end
    end
    return descriptors
end

function validate2(data, descriptors, potential, optim, coeff, pbc=nothing, a=nothing, b=nothing, c=nothing)    
    
    data = deletenothing(data)
    descriptors = deletenothing(descriptors)
    potential = deletenothing(potential)

    lossfunc, method, normalizeenergy, normalizestress, eta, kappa, 
        weightouter, etaspace, kappalist = getoptim(optim)
    
    n = length(data)
    emae = -ones(Float64,n)
    fmae = -ones(Float64,n)
    smae = -ones(Float64,n)
    ermse = -ones(Float64,n)
    frmse = -ones(Float64,n)
    srmse = -ones(Float64,n)
    ersq = -ones(Float64,n)
    frsq = -ones(Float64,n)
    srsq = -ones(Float64,n)
    essr = -ones(Float64,n)
    fssr = -ones(Float64,n)
    sssr = -ones(Float64,n)
    eavg = -ones(Float64,n)
    favg = -ones(Float64,n)
    savg = -ones(Float64,n)
    szbe = zeros(Int64, n)
    szbf = zeros(Int64, n)
    szbs = zeros(Int64, n)
    be = []
    bf = []
    bs = []
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, indices = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        valid = Array(1:config.nconfigs)      

        Aei, Afi, Asi, bei, bfi, bsi, nfi = Potential.emldescriptors2(config, valid, eta, kappa, 
                        normalizeenergy, normalizestress, descriptors, potential)

        szbe[i] = length(bei)        
        szbf[i] = length(bfi)        
        szbs[i] = length(bsi)        

        emae[i], ermse[i], ersq[i], essr[i], eavg[i] = calerrors(Aei, bei, coeff)
        fmae[i], frmse[i], frsq[i], fssr[i], favg[i] = calerrors(Afi, bfi, coeff)
        smae[i], srmse[i], srsq[i], sssr[i], savg[i] = calerrors(Asi, bsi, coeff)

        if szbe[i] > 0
            be = [be; bei]
        end 
        if szbf[i] > 0
            bf = [bf; bfi]
        end 
        if szbs[i] > 0
            bs = [bs; bsi]
        end 
    end

    energyerrors = -ones((n+1),3)    
    forceerrors = -ones((n+1),3)    
    stresserrors = -ones((n+1),3)    

    if sum(szbe) > 0 
        ssr = sum(essr)
        avg = sum(eavg.*szbe)/sum(szbe)          
        rsq = 1.0 - ssr/sum((be .- avg).^2);
        energyerrors[1,1] = sum(emae.*szbe)/sum(szbe)  # MAE
        energyerrors[1,2] = sqrt(sum((ermse.^2).*szbe)/sum(szbe)) # RMSE 
        energyerrors[1,3] = rsq
        energyerrors[2:end,:] = [emae ermse ersq]        
    end

    if sum(szbf) > 0 
        ssr = sum(fssr)
        avg = sum(favg.*szbf)/sum(szbf)          
        rsq = 1.0 - ssr/sum((bf .- avg).^2);
        forceerrors[1,1] = sum(fmae.*szbf)/sum(szbf)
        forceerrors[1,2] =  sqrt(sum((frmse.^2).*szbf)/sum(szbf))    
        forceerrors[1,3] = rsq 
        forceerrors[2:end,:] = [fmae frmse frsq]   
    end    
    
    if sum(szbs) > 0 
        ssr = sum(sssr)
        avg = sum(savg.*szbs)/sum(szbs)          
        rsq = 1.0 - ssr/sum((bs .- avg).^2);
        stresserrors[1,1] = sum(smae.*szbs)/sum(szbs)
        stresserrors[1,2] = sqrt(sum((srmse.^2).*szbs)/sum(szbs))        
        stresserrors[1,3] = rsq 
        stresserrors[2:end,:] = [smae srmse srsq]   
    end

    return energyerrors, forceerrors, stresserrors
end

function maespace2(traindata, validdata, descriptors, potential, optim, pbc=nothing, a=nothing, b=nothing, c=nothing)

    # generate the eta space grid
    etaspace = optim.etaspace
    etamin = etaspace[:,1]
    etamax = etaspace[:,2]
    N = Int64.(etaspace[:,3])
    eta = tensornodes(etamin, etamax, N)
        
    # read data 
    trainconfig = Array{Any}(undef, length(traindata))
    for i = 1:length(traindata)
        trainconfig[i], indices = Preprocessing.readconfigdata(traindata[i], pbc, a, b, c)                  
    end
    validconfig = Array{Any}(undef, length(validdata))
    for i = 1:length(validdata)
        validconfig[i], indices = Preprocessing.readconfigdata(validdata[i], pbc, a, b, c)                  
    end

    neta = size(eta,1)
    emae = -ones(neta)
    fmae = -ones(neta)
    smae = -ones(neta)
    for n = 1:neta
        optim.eta = eta[n,:]
        descriptors = setcutoff2(descriptors, eta[n,:])      
        coeff,~ = linearfit2(trainconfig, descriptors, potential, optim, pbc, a, b, c)
        eerr, ferr, serr = validate2(validconfig, descriptors, potential, optim, coeff, pbc, a, b, c)    
        emae[n] = eerr[1]
        fmae[n] = ferr[1]
        smae[n] = serr[1]                
        a1 = num2string(emae[n], 16); 
        a2 = num2string(fmae[n], 16); 
        a3 = num2string(smae[n], 16);
        mystr = "$n  "
        for j = 1:length(eta[n,:])
            a = num2string(eta[n,j], 16)
            mystr = mystr * a * "  "
        end
        print(mystr * a1 * "  " * a2 * "  " * a3 * " \n")
    end

    return eta, emae, fmae, smae
end

function optimize2(traindata, testdata, descriptors, potential, optim, pbc=nothing, a=nothing, b=nothing, c=nothing)

    traindata = deletenothing(traindata)
    testdata = deletenothing(testdata)
    descriptors = deletenothing(descriptors)
    potential = deletenothing(potential)

    # compute errors on the eta space grid
    etapts, emae, fmae, smae = maespace2(traindata, testdata, descriptors, potential, optim, pbc, a, b, c)    

    lossfunc = optim.lossfunc
    weightouter = optim.weightouter
    etaspace = optim.etaspace

    # compute the loss function on the eta space grid
    if lossfunc == "energy" 
        f = weightouter[1]*emae
    elseif lossfunc == "force" 
        f = weightouter[2]*emae
    elseif lossfunc == "energyforce" 
        f = weightouter[1]*emae + weightouter[2]*fmae
    elseif lossfunc == "energyforcestress" 
        f = weightouter[1]*emae + weightouter[2]*fmae + weightouter[3]*smae
    else
        error("lossfunc is not valid. It must be energy, force, energyforce, or energyforcestress")
    end    

    # interpolate the loss function on the eta space grid
    pc = mpolyinterp(etapts, f, etaspace)

    # optimize the loss function on the eta space
    eta, fmin, iter  = gradientdescent(pc, etaspace)
    optim.eta = eta

    # set cut-off radius
    descriptors = setcutoff2(descriptors, eta)        
    
    # compute coefficient 
    coeff,~ = linearfit2(traindata, descriptors, potential, optim, pbc, a, b, c)

    return eta, coeff, fmin, iter, pc, f, etapts, emae, fmae, smae
end

