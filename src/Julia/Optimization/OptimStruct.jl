mutable struct OptimStruct

    lossfunc::String 
    method::String     
    normalizeenergy::Bool
    normalizestress::Bool        
    eta::Vector{Float64}
    kappa::Vector{Int64}
    weightouter::Vector{Float64}
    etaspace::Matrix{Float64}
    kappalist::Matrix{Int64}    

end

function setoptim(lossfunc=nothing, method=nothing, normalizeenergy=nothing, normalizestress=nothing,
     eta=nothing, kappa=nothing, weightouter=nothing, etaspace=nothing, kappalist=nothing)

    if lossfunc === nothing
        lossfunc = "energyforce"
    end
    if method === nothing
        method = "lsq"
    end
    if normalizeenergy === nothing
        normalizeenergy = false
    end
    if normalizestress === nothing
        normalizestress = false
    end
    if eta === nothing
        eta = []
    end
    if kappa === nothing
        kappa = []
    end
    if weightouter === nothing
        weightouter = [1.0, 1.0, 1.0]
    end
    if etaspace === nothing
        etaspace = reshape([0.0 1.0],(1,2))
    end
    if kappalist === nothing
        kappalist = reshape([0],(1,1))
    end

    optim = OptimStruct(lossfunc, method, normalizeenergy, normalizestress, eta, kappa, weightouter, etaspace, kappalist)

    return optim
end

function getoptim(optim::OptimStruct)

    lossfunc = optim.lossfunc
    method = optim.method    
    normalizeenergy = optim.normalizeenergy 
    normalizestress = optim.normalizestress
    eta = optim.eta
    kappa = optim.kappa
    weightouter = optim.weightouter
    etaspace = optim.etaspace
    kappalist = optim.kappalist 

    return lossfunc, method, normalizeenergy, normalizestress, eta, kappa, weightouter, etaspace, kappalist
end

