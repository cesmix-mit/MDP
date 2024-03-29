# External packages
using Revise, DelimitedFiles, SymPy

version = "Version0.1"
cdir = pwd();
push!(LOAD_PATH, cdir);

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Julia", cdir);
sourcepath = cdir[1:(ii[1]-1)];
include(sourcepath * "/Installation/setpath.jl");

# MDP packages
using Preprocessing, Gencode, Postprocessing

app = Preprocessing.initializeapp(sourcepath,version);


#using Gencode
# create pde structure and mesh structure
#app = Gencode.initializeapp(version);

app.dim = 3;
app.ncq = 1;
app.nceta = 10;
app.nckappa = 2;
app.ncmu1a = 1;
app.ncmu1b = 1;
app.ncmu2a = 10;
app.ncmu2b = 10;
app.ncmu2c = 10;
app.ncmu3a = 10;
app.ncmu3b = 10;
app.ncmu3c = 10;
app.ncmu4a = 4;
app.ncmu4b = 4;
app.natomtype = 2;
app.potentialfile = "potential.jl";

include(app.potentialfile);  # include the potential file
Gencode.gencode(app); # Generate C++ code
