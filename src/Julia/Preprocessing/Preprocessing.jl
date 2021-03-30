__precompile__()

module Preprocessing

using Revise

export initializeapp, initializeconfig, initializemdp, preprocessing, readconfig

include("initializeapp.jl");
include("initializeconfig.jl");
include("initializemdp.jl");
include("readconfig.jl");
include("writeapp.jl");
include("writeconfig.jl");
include("checkconfig.jl");
include("readLAMMPS.jl");
include("cubemapping.jl");
include("boxperiodicimages.jl");

function preprocessing(app)

# read configurations from data file
config = readconfig(app, app.configmode);     

# number of configurations
app.nconfigs = config.nconfigs;          

# nonbonded single potentials
if  length(app.pot1a[:]) > 0
    app.np1a = length(app.pot1a);  # number of nonbonded single potentials 
    app.ncmu1a = length(app.mu1a); # length of mu1a
else
    app.np1a = 0;
    app.ncmu1a = 0;
end

# bonded single potentials
if  length(app.pot1b[:]) > 0
    app.np1b = length(app.pot1b);  # number of nonbonded single potentials 
    app.ncmu1b = length(app.mu1b); # length of mu1a
else
    app.np1b = 0;
    app.ncmu1b = 0;
end

# nonbonded pair potentials
if  length(app.pot2a[:]) > 0
    app.np2a = length(app.pot2a);  # number of nonbonded single potentials 
    app.ncmu2a = length(app.mu2a); # length of mu1a
else
    app.np2a = 0;
    app.ncmu2a = 0;
end

# bonded pair potentials
if  length(app.pot2b[:]) > 0
    app.np2b = length(app.pot2b);  # number of nonbonded single potentials 
    app.ncmu2b = length(app.mu2b); # length of mu1a
else
    app.np2b = 0;
    app.ncmu2b = 0;
end

# two-body bond order potentials
if  length(app.pot2c[:]) > 0
    app.np2c = length(app.pot2c);  # number of nonbonded single potentials 
    app.ncmu2c = length(app.mu2c); # length of mu1a
else
    app.np2c = 0;
    app.ncmu2c = 0;
end

# nonbonded triplet potentials
if  length(app.pot3a[:]) > 0
    app.np3a = length(app.pot3a);  # number of nonbonded single potentials 
    app.ncmu3a = length(app.mu3a); # length of mu1a
else
    app.np3a = 0;
    app.ncmu3a = 0;
end

# bonded triplet potentials
if  length(app.pot3b[:]) > 0
    app.np3b = length(app.pot3b);  # number of nonbonded single potentials 
    app.ncmu3b = length(app.mu3b); # length of mu1a
else
    app.np3b = 0;
    app.ncmu3b = 0;
end

# three-body bond order potentials
if  length(app.pot3c[:]) > 0
    app.np3c = length(app.pot3c);  # number of nonbonded single potentials 
    app.ncmu3c = length(app.mu3c); # length of mu1a
else
    app.np3c = 0;
    app.ncmu3c = 0;
end

# nonbonded quadruplet potentials
if  length(app.pot4a[:]) > 0
    app.np4a = length(app.pot4a);  # number of nonbonded single potentials 
    app.ncmu4a = length(app.mu4a); # length of mu1a
else
    app.np4a = 0;
    app.ncmu4a = 0;
end

# bonded quadruplet potentials
if  length(app.pot4b[:]) > 0
    app.np4b = length(app.pot4b);  # number of nonbonded single potentials 
    app.ncmu4b = length(app.mu4b); # length of mu1a
else
    app.np4b = 0;
    app.ncmu4b = 0;
end

app.flag = [app.descriptor app.spectrum app.training app.runMD app.potentialform app.neighpair app.energycal app.forcecal app.stresscal app.neighcell app.decomposition app.chemtype app.dftdata];

rcut = [reshape([app.rcutml],1,1) app.rcut2a app.rcut2b app.rcut2c app.rcut3a app.rcut3b app.rcut3c app.rcut4a app.rcut4b];
rcutmax = maximum(rcut);     
app.rcutsqmax = rcutmax^2;        
app.boxoffset = [rcutmax rcutmax rcutmax];
app.simparam = [app.time app.dt];        

app.ndims = zeros(20,1);
app.ndims[1] = app.dim;
app.ndims[2] = app.L;
app.ndims[3] = app.K;
app.ndims[4] = app.ntimesteps;
app.ndims[5] = app.nab;
app.ndims[6] = app.mpiprocs;  
app.ndims[7] = app.backend;
app.ndims[8] = app.nconfigs;
app.ndims[9] = app.natomtype;
app.ndims[10] = app.nmoletype; 
app.ndims[11] = app.maxnumneighbors;

m = 1;
for i = 1:config.nconfigs
    B2C, C2B = cubemapping(config.a[:,i], config.b[:,i], config.c[:,i]);
    ximages = boxperiodicimages(app.pbc, config.a[:,i], config.b[:,i], config.c[:,i]);    
    n = config.natom[i];
    config.x[:,m:(m+n-1)] = checkconfig(config.x[:,m:(m+n-1)], ximages, B2C, C2B);    
    m = m + n;
end

writeapp(app, app.appname * "app.bin");
writeconfig(config, app.appname * "config.bin");
               
mv(app.appname * "app.bin", app.sourcepath * "C++/Main/" * app.appname * "app.bin", force = true);
mv(app.appname * "config.bin", app.sourcepath * "C++/Main/" * app.appname * "config.bin", force = true);

return app, config

end

end
