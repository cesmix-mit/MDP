#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function mdp(app)

# create exec folder if it does not exist
bindir = "exec";
cd(app.sourcepath);
if isdir(bindir) == 0
    mkdir(bindir);    
end
cd(app.currentdir);

# preprocess input data
app, config = Preprocessing.preprocessing(app); 

# generate code
Gencode.gencode(app); # Generate C++ code

# compile C++ source code
# cd(app.sourcepath * "C++/Main");
# Gencode.compile(app);
cd(app.sourcepath * bindir);
if isfile("CMakeCache.txt")    
    rm("CMakeCache.txt", force=true);
end
if app.platform == "gpu"
    if isfile("libcpuCore.a") && isfile("libgpuCore.a") && isfile("gpuMDP")
        run(Gencode.string2cmd("cmake -D MDP_POTENTIALS=ON -D MDP_CUDA=ON ../Installation"));
    else
        run(Gencode.string2cmd("cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON ../Installation"));
    end
elseif app.platform == "cpu"
    if isfile("libcpuCore.a") && isfile("cpuMDP")
        run(Gencode.string2cmd("cmake -D MDP_POTENTIALS=ON ../Installation"));
    else
        run(Gencode.string2cmd("cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation"));
    end
end
run(Gencode.string2cmd("cmake --build ."));

# run code
#Gencode.runcode(app);
if app.platform == "cpu"
    run(Gencode.string2cmd("./cpuMDP " * app.appname * " out"));
else app.platform == "gpu"
    run(Gencode.string2cmd("./gpuMDP " * app.appname * " out"));
end
cd(app.currentdir);

if app.training > 0
    filename = app.sourcepath * "C++/Main/coefficients.bin";
    tmp = reinterpret(Float64,read(filename));
    app.coeff = reshape(tmp, 1, length(tmp));
end

return app, config

end
