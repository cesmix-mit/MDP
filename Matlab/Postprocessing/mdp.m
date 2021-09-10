%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [app, config] = mdp(app)

% create exec folder if it does not exist
bindir = "exec";
cd(char(app.sourcepath));
if exist(bindir, "file") == 0
    mkdir(char(bindir));    
end
cd(char(app.currentdir));

% preprocess input data
[app,config] = preprocessing(app); 

% generate code
gencode(app);                      

%compile(app);
%buildpotlib2(app);

% compile C++ source code
cd(char(app.sourcepath + bindir));
if exist("CMakeCache.txt", "file")
    delete(char("CMakeCache.txt"));
end
if app.platform == "gpu"
    if app.buildexec
        eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON ../Installation");
    else
        if exist("libcpuCore.a", "file") && exist("libgpuCore.a", "file") && exist("gpuMDP", "file")
            eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CUDA=ON ../Installation");
        else
            eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON ../Installation");
        end
    end
elseif app.platform == "cpu"
    if app.buildexec
        eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation");
    else
        if exist("libcpuCore.a", "file") && exist("cpuMDP", "file")
           eval("!cmake -D MDP_POTENTIALS=ON ../Installation");
        else
            eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation");
        end
    end
end
eval("!cmake --build .");

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





















% cdir = pwd(); ii = strfind(cdir, "Applications");
% sourcepath = cdir(1:(ii-1));
% run(sourcepath + "Installation/setpath.m");
% 
% app = initializeapp(sourcepath);
% app.appname = "snapTA06A";
% app.currentdir = cdir; % current directory
% app.cpucompiler = "/usr/local/Cellar/llvm/11.1.0/bin/clang++";
% app.cpuflags = "-Xclang -load -Xclang /usr/local/lib/ClangEnzyme-11.dylib";
% app.cpumacros = "-std=c++11 -D _ENZYME";
% app.cpuflags = "";
% app.cpumacros = "-std=c++11";
% app.cpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
% app.cpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
% app.cpumacros = "-std=c++11 -D _ENZYME";
% app.gpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
% app.gpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
% app.gpumacros = "--cuda-gpu-arch=sm_60 -std=c++11 -D _ENZYME";
% app.buildcorelib = 1;
% app.globalfreq = 10;        % frequency to print and save default global outputs (pe, ke, ce, te, temp, press)
% app.peratomfreq = 20;       % frequency to save default peratom outputs (id, t, x, f)
% app.globaloutputs = ["stresses"];    % additional global outputs (stresses, com, vcm)
% app.peratomoutputs = ["velocity", "mass", "image"];   % additional peratom outputs (velocity, mass, virial, image)
