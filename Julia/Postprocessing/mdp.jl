#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function mdp(app)

# preprocess input data
app, config = Preprocessing.preprocessing(app); 

# generate code
Gencode.gencode(app); # Generate C++ code

# compile code
cd(app.sourcepath * "C++/Main");
Gencode.compile(app);

# run code
#Gencode.runcode(app);
if app.platform == "cpu"
    run(Gencode.string2cmd("./cpuMDP " * app.appname * " out"));
else app.platform == "gpu"
    run(Gencode.string2cmd("./gpuMDP " * app.appname * " out"));
end
cd(app.currentdir);

#cd(app.currentdir);
if app.training > 0
    filename = app.sourcepath * "C++/Main/coefficients.bin";
    tmp = reinterpret(Float64,read(filename));
    app.coeff = reshape(tmp, 1, length(tmp));
end

return app, config

end
