function [app, config] = mdp(app)

% preprocess input data
[app,config] = preprocessing(app); 

% generate code
gencode(app);                      

% compile code
cd(app.sourcepath + "C++/Main");
compile(app.cpucompiler, app.gpucompiler);

% run code
eval("!./cpuMDP " + app.appname + " out");
cd(app.currentdir);


