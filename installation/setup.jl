# Add MDP to Julia search path
push!(LOAD_PATH, MDPpath * "src/Julia");
push!(LOAD_PATH, MDPpath * "src/Julia/Gencode");
push!(LOAD_PATH, MDPpath * "src/Julia/Preprocessing");
push!(LOAD_PATH, MDPpath * "src/Julia/Postprocessing");
push!(LOAD_PATH, MDPpath * "src/Julia/Potential");
push!(LOAD_PATH, MDPpath * "src/Julia/Optimization");
push!(LOAD_PATH, MDPpath * "src/Julia/Integration");

# Set Julia's PATH enviroment variable so that MDP can call external programs
ENV["PATH"] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";

using Preprocessing, Potential, Optimization, Integration 

# load Snap library 
Potential.loadsnap(MDPpath)

# load integration library 
Integration.loadintegration(MDPpath)

# initialize arrays with nothing 
traindata = Array{Any}(nothing, 100)
testdata = Array{Any}(nothing, 100)
validdata = Array{Any}(nothing, 100)
descriptors = Array{Any}(nothing, 100)
potentials  = Array{Any}(nothing, 100)

# initialize app 
app = Preprocessing.initializeapp(cdir,MDPpath);
