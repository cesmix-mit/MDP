function [app, config] = mdp(app)

% preprocess input data
[app,config] = preprocessing(app); 

% generate code
gencode(app);                      

% compile code
cd(app.sourcepath + "C++/Main");
compile(app);

% run code
eval("!./cpuMDP " + app.appname + " out");

cd(app.currentdir);

if app.training > 0
    filename = app.sourcepath + "C++/Main/coefficients.bin";
    fileID = fopen(filename,'r');
    app.coeff = fread(fileID,'double');
    fclose(fileID);
end


