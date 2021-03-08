mutable struct AppStruct
    codename::String; # MDP
    version::String;  # Set version
    appname::String;  # application name
    platform::String; # CPU or GPU
    cpucompiler::String; # Path to CPU compiler
    mpicompiler::String; # Path to MPI compiler
    gpucompiler::String; # Path to GPU compiler
    mpirun::String;      # Path to MPI run command and MPI run options
    cpuflags::String;    # options for CPU compiler
    gpuflags::String;    # options for GGU compiler
    model::String;# used to indicate potential model
    modelfile::String;# model file name
    modelnumber::Int64;
    potentialfile::String;# model file name    
    preprocessmode::Int64; # preprocessing mode
    mpiprocs::Int64; # number of MPI ranks
    
    nd::Int64; # physical dimension
    dim::Int64; # physical dimension
    natomtype::Int64; # number of atom types
    ncq::Int64; # number of compoments of q
    ncmu::Int64; # number of compoments of eta
    nceta::Int64; # number of compoments of eta
    nckappa::Int64;# number of compoments of kappa
    ncmu1a::Int64;# number of compoments mu1a
    ncmu1b::Int64;# number of compoments mu1b
    ncmu2a::Int64;# number of compoments mu2a
    ncmu2b::Int64;# number of compoments mu2b
    ncmu2c::Int64;# number of compoments mu2c
    ncmu3a::Int64;# number of compoments mu3a
    ncmu3b::Int64;# number of compoments mu3b
    ncmu3c::Int64;# number of compoments mu3c
    ncmu4a::Int64;# number of compoments mu4a
    ncmu4b::Int64;# number of compoments mu4b
    
    flag::Array{Int64,2};   # flag parameters
    problem::Array{Int64,2};# problem parameters
    factor::Array{Float64,2};  # factors
    physicsparam::Array{Float64,2}; # physical parameters
    solversparam::Array{Float64,2}; # solvers parameters
    mu::Array{Float64,2}; # physical parameters
    eta::Array{Float64,2}; # physical parameters
    kappa::Array{Int64,2}; # physical parameters
    
    AppStruct() = new();
end


function initializeapp(version)
    app = AppStruct();

    app.codename = "Exasim";
    app.version = version;
    app.appname = "app";
    app.platform = "cpu";
    app.cpucompiler = "g++";
    app.mpicompiler = "mpicxx";
    app.gpucompiler = "nvcc";
    app.mpirun = "mpirun";
    app.cpuflags = "-O2 -ldl -lm -lblas -llapack";
    app.gpuflags = "-lcudart -lcublas";
    app.model="ModelD";
    app.modelfile = "";    
    app.modelnumber = 0;
    app.potentialfile = "";
    app.preprocessmode = 1;
    app.mpiprocs = 1;

    app.nd = 3;
    app.dim = 3;
    app.natomtype = 1;
    app.ncq = 1;
    app.ncmu = 1;
    app.nceta = 1;
    app.nckappa = 1;
    app.ncmu1a = 1;
    app.ncmu1b = 1;
    app.ncmu2a = 1;
    app.ncmu2b = 1;
    app.ncmu2c = 1;
    app.ncmu3a = 1;
    app.ncmu3b = 1;
    app.ncmu3c = 1;
    app.ncmu4a = 1;
    app.ncmu4b = 1;
    
    app.flag = [0 0];
    app.problem = [0 0];
    app.factor = [0.0 0.0];  # factors
    app.physicsparam = [0.0 0.0]; # physical parameters
    app.solversparam = [0.0 0.0]; # solvers parameters
    app.mu = [0.0 0.0]; # physical parameters
    app.eta = [0.0 0.0]; # solvers parameters
    app.kappa = [0 0];
    
    return app;
end

