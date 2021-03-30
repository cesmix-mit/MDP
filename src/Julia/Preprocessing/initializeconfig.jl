mutable struct CONFIGStruct

    dim::Int64;
    ncx::Int64;
    ncv::Int64;
    ncq::Int64;
    nce::Int64;
    ncf::Int64;
    nconfigs::Int64;
    natom::Array{Int64,2};   # number of atoms per each configuration
    natomall::Int64;# number of atoms for all configurations
    
    # simulation box for each configuration
    a::Array{Float64,2}; # the 1st principal vector of the simulation box
    b::Array{Float64,2}; # the 2nd principal vector of the simulation box
    c::Array{Float64,2}; # the 3rd principal vector of the simulation box
    e::Array{Float64,2};  # energies for all configurations
    
    t::Array{Int64,2};   # atom types for all configurations
    x::Array{Float64,2}; # atom positions for all configurations
    q::Array{Float64,2}; # atom charges for all configurations
    v::Array{Float64,2}; # atom velocities for all configurations
    f::Array{Float64,2}; # atom forces for all configurations    

    CONFIGStruct() = new();
end

function initializeconfig(app)

    config = CONFIGStruct();    
    nconfigs = app.nconfigs;            # number of configurations
    dim = app.dim;
    ncx = app.ncx;
    ncv = app.ncv;
    ncf = app.ncf;
    ncq = app.ncq;
    nce = app.nce;

    config.dim = dim;
    config.ncx = ncx;
    config.ncv = ncv;
    config.ncq = ncq;
    config.nce = nce;
    config.ncf = ncf;
    config.nconfigs = nconfigs;
    config.natom = ones(1, nconfigs);   # number of atoms per each configuration
    config.natomall = sum(config.natom);# number of atoms for all configurations

    # simulation box for each configuration
    config.a = zeros(dim, nconfigs); # the 1st principal vector of the simulation box
    config.b = zeros(dim, nconfigs); # the 2nd principal vector of the simulation box
    config.c = zeros(dim, nconfigs); # the 3rd principal vector of the simulation box
    config.e = ones(nce, nconfigs);  # energies for all configurations

    config.t = zeros(1, config.natomall);   # atom types for all configurations
    config.x = zeros(ncx, config.natomall); # atom positions for all configurations
    config.q = zeros(ncq, config.natomall); # atom charges for all configurations
    config.v = zeros(ncv, config.natomall); # atom velocities for all configurations
    config.f = zeros(ncf, config.natomall); # atom forces for all configurations

    return config;
end
