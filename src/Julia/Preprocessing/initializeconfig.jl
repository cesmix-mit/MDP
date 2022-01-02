#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

mutable struct CONFIGStruct

    dim::Int64;
    ncx::Int64;
    ncv::Int64;
    ncq::Int64;
    nce::Int64;
    ncf::Int64;
    ncs::Int64;
    nci::Int64;
    nct::Int64;
    ncg::Int64;
    ncz::Int64;
    ncm::Int64;
    nco::Int64;
    ncl::Int64;
    ncp::Int64;
    nconfigs::Int64;
    natom::Array{Int64,2};   # number of atoms per each configuration
    natomall::Int64;# number of atoms for all configurations
    
    training::Array{Int64,1};
    validation::Array{Int64,1};
    testing::Array{Int64,1}; 

    # simulation box for each configuration
    a::Array{Float64,2}; # the 1st principal vector of the simulation box
    b::Array{Float64,2}; # the 2nd principal vector of the simulation box
    c::Array{Float64,2}; # the 3rd principal vector of the simulation box
    lattice::Array{Float64,2};  # energies for all configurations
    pbc::Array{Int64,2};  # energies for all configurations
    e::Array{Float64,2};  # energies for all configurations
    stress::Array{Float64,2};  # energies for all configurations

    we::Array{Float64,2};  # energies weights for all configurations
    wf::Array{Float64,2};  # energies weights for all configurations
    ws::Array{Float64,2};  # energies weights for all configurations

    t::Array{Int64,2};   # atom types for all configurations
    Z::Array{Int64,2};   # atom types for all configurations
    move::Array{Int64,2};   # atom types for all configurations
    tags::Array{Int64,2};   # atom types for all configurations
    group::Array{Int64,2};   # atom types for all configurations
    mass::Array{Float64,2}; # atom positions for all configurations
    x::Array{Float64,2}; # atom positions for all configurations
    q::Array{Float64,2}; # atom charges for all configurations
    v::Array{Float64,2}; # atom velocities for all configurations
    f::Array{Float64,2}; # atom forces for all configurations    

    CONFIGStruct() = new();
end

function initializeconfig(app)

    config = CONFIGStruct();     
    nconfigs = 0;            # number of configurations
    dim = 0
    nci = 0;
    nct = 0;
    ncg = 0;
    ncx = 0;
    ncv = 0;
    ncf = 0;
    ncq = 0;
    nce = 0;
    ncs = 0;
    ncz = 0;
    ncm = 0;
    nco = 0;
    ncl = 0;
    ncp = 0;

    config.nconfigs = nconfigs;
    config.dim = dim;
    config.ncx = ncx;
    config.ncv = ncv;
    config.ncq = ncq;
    config.nce = nce;
    config.ncf = ncf;
    config.ncs = ncs;
    config.nci = nci;
    config.nct = nct;
    config.ncg = ncg;
    config.ncz = ncz;
    config.ncm = ncm;
    config.nco = nco;
    config.ncl = ncl;
    config.ncp = ncp;
    
    config.natom = ones(1, nconfigs);   # number of atoms per each configuration
    config.natomall = sum(config.natom);# number of atoms for all configurations

    config.training = Array(1:nconfigs)  # training configurations
    config.validation = zeros(Int64,nconfigs)   # validation configurations
    config.testing = zeros(Int64,nconfigs)      # testing configurations

    # simulation box for each configuration
    config.a = zeros(dim, nconfigs); # the 1st principal vector of the simulation box
    config.b = zeros(dim, nconfigs); # the 2nd principal vector of the simulation box
    config.c = zeros(dim, nconfigs); # the 3rd principal vector of the simulation box
    config.lattice = zeros(ncl, nconfigs); # stresses for all configurations
    config.pbc = zeros(ncp, nconfigs); # periodic boundary conditions
    config.e = zeros(nce, nconfigs);  # energies for all configurations
    config.stress = zeros(ncs, nconfigs); # stresses for all configurations

    config.we = ones(1, nconfigs);      # energy weight per each configuration
    config.wf = ones(1, nconfigs);      # force weight per each configuration
    config.ws = ones(1, nconfigs);      # stress weight per each configuration

    config.Z = zeros(ncz, config.natomall);    # atom numbers for all configurations
    config.mass = zeros(ncm, config.natomall); # atom masses for all configurations
    config.move = zeros(nco, config.natomall);   # atom move masks for all configurations
    config.tags = zeros(nci, config.natomall); # atom tags for all configurations
    config.t = zeros(nct, config.natomall);   # atom types for all configurations
    config.group = zeros(ncg, config.natomall);   # atom groups for all configurations
    config.x = zeros(ncx, config.natomall); # atom positions for all configurations
    config.q = zeros(ncq, config.natomall); # atom charges for all configurations
    config.v = zeros(ncv, config.natomall); # atom velocities for all configurations
    config.f = zeros(ncf, config.natomall); # atom forces for all configurations
    
    return config;
end

function catconfig(config1, config2)

    config = initializeconfig(0)
    config.nconfigs = config1.nconfigs + config2.nconfigs
    config.dim = config1.dim;
    config.ncx = config1.ncx;
    config.ncv = config1.ncv;
    config.ncq = config1.ncq;
    config.nce = config1.nce;
    config.ncf = config1.ncf;
    config.ncs = config1.ncs;
    config.nci = config1.nci;
    config.nct = config1.nct;
    config.ncg = config1.ncg;
    config.ncz = config1.ncz;
    config.ncm = config1.ncm;
    config.nco = config1.nco;
    config.ncl = config1.ncl;
    config.ncp = config1.ncp;
    
    config.natom = cat(config1.natom, config2.natom, dims=2)
    config.natomall = sum(config.natom);# number of atoms for all configurations

    config.a = cat(config1.a, config2.a, dims=2)
    config.b = cat(config1.b, config2.b, dims=2)
    config.c = cat(config1.c, config2.c, dims=2)
    config.lattice = cat(config1.lattice, config2.lattice, dims=2)
    config.pbc = cat(config1.pbc, config2.pbc, dims=2)
    config.e = cat(config1.e, config2.e, dims=2)
    config.stress = cat(config1.stress, config2.stress, dims=2)

    config.we = cat(config1.we, config2.we, dims=2)
    config.wf = cat(config1.wf, config2.wf, dims=2)
    config.ws = cat(config1.ws, config2.ws, dims=2)

    config.Z = cat(config1.Z, config2.Z, dims=2)
    config.mass = cat(config1.mass, config2.mass, dims=2)
    config.move = cat(config1.move, config2.move, dims=2)
    config.tags = cat(config1.tags, config2.tags, dims=2)
    config.t = cat(config1.t, config2.t, dims=2)
    config.group = cat(config1.group, config2.group, dims=2)
    config.x = cat(config1.x, config2.x, dims=2)
    config.q = cat(config1.q, config2.q, dims=2)
    config.v = cat(config1.v, config2.v, dims=2)
    config.f = cat(config1.f, config2.f, dims=2)
    
    return config
end

function extractarray(a, indices)

    if size(a,2) > 1
        a = a[:,indices]
    end

    return a
end

function extractconfig(config, indices)

    natom = config.natom[:]
    config.natom = config.natom[:,indices]
    config.natomall = sum(config.natom);    
    config.nconfigs = length(config.natom)

    config.a = extractarray(config.a, indices)
    config.b = extractarray(config.b, indices)
    config.c = extractarray(config.c, indices)
    config.e = extractarray(config.e, indices)
    config.pbc = extractarray(config.pbc, indices)
    config.lattice = extractarray(config.lattice, indices)    
    config.stress = extractarray(config.stress, indices)

    config.we = extractarray(config.we, indices)
    config.wf = extractarray(config.wf, indices)
    config.ws = extractarray(config.ws, indices)

    natom = [0; cumsum(natom)]
    ind = []
    for i = 1:length(indices)
        j = indices[i]
        if i == 1
            ind = Array((natom[j]+1):natom[j+1]);
        else
            ind = [ind; Array((natom[j]+1):natom[j+1])];
        end
    end

    config.Z = extractarray(config.Z, ind)
    config.mass = extractarray(config.mass, ind)
    config.move = extractarray(config.move, ind)
    config.tags = extractarray(config.tags, ind)
    config.t = extractarray(config.t, ind)
    config.group = extractarray(config.group, ind)
    config.x = extractarray(config.x, ind)
    config.q = extractarray(config.q, ind)
    config.v = extractarray(config.v, ind)
    config.f = extractarray(config.f, ind)
    
    return config
end


using Random

function traingconfigs(nconfigs, training_percentage, validation_percentage, testing_percentage, randomize)

    if randomize==1
        ind = randperm(nconfigs)
    else
        ind = Array(1:nconfigs)
    end

    ntrain = Int64(round(training_percentage*nconfigs/100.0))    
    nvalid = Int64(round(validation_percentage*nconfigs/100.0))
    ntest = Int64(round(testing_percentage*nconfigs/100.0))

    training = ind[1:ntrain]
    validation = ind[(ntrain+1):(ntrain+nvalid)]
    testing = ind[(ntrain+nvalid+1):min(ntrain+nvalid+ntest, nconfigs)]

    return training, validation, testing
end

function setconfig(config, training_percentage, validation_percentage, testing_percentage, randomize)

    config.training, config.validation, config.testing = 
        traingconfigs(config.nconfigs, training_percentage, validation_percentage, testing_percentage, randomize)

    return config
end

