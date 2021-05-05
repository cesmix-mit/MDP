function [app, config] = mdp(app)

% preprocess input data
[app,config] = preprocessing(app); 

% generate code
gencode(app);                      

% compile code
cd(char(app.sourcepath + "C++/Main"));
compile(app);

% run code
if app.platform == "cpu"
    eval("!./cpuMDP " + app.appname + " out");
elseif app.platform == "gpu"
    eval("!./gpuMDP " + app.appname + " out");    
end

cd(char(app.currentdir));

if app.training > 0
    filename = app.sourcepath + "C++/Main/coefficients.bin";
    fileID = fopen(filename,'r');
    app.coeff = fread(fileID,'double');
    fclose(fileID);
end


