%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function buildcorelib(app)

disp("Building core libraries");

cpucompiler = app.cpucompiler;
gpucompiler = app.gpucompiler;
gpumacros = app.gpumacros;

if exist("../Core/cpuCore.a", "file") == 0
    comstr = cpucompiler + " -fPIC -ffast-math -O3 -c ../Core/cpuCore.cpp -o ../Core/cpuCore.o";
    eval(char("!" + comstr));
    comstr = "ar rvs ../Core/cpuCore.a ../Core/cpuCore.o";
    eval(char("!" + comstr));
    if ismac==1
        comstr = cpucompiler + " --shared ../Core/cpuCore.o -o ../Core/libcpuCore.dylib";
    else
        comstr = cpucompiler + " --shared ../Core/cpuCore.o -o ../Core/libcpuCore.so";
    end
    eval(char("!" + comstr));
else
    disp("cpuCore.a already exists. Delete it if you want to rebuild core libaries.");
end

if exist("../Core/gpuCore.a", "file") == 0
    if ~isempty(char(gpucompiler))
        comstr = gpucompiler + " " + gpumacros + " -D_FORCE_INLINES -O3 -c -fPIC -w ../Core/gpuCore.cu -o ../Core/gpuCore.o";
        eval(char("!" + comstr));
        comstr = "ar rvs ../Core/gpuCore.a ../Core/gpuCore.o";
        eval(char("!" + comstr));
        if ismac==1
            comstr = gpucompiler + " --shared ../Core/gpuCore.o -o ../Core/libgpuCore.dylib";
        else
            comstr = gpucompiler + " --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so";
        end
        eval(char("!" + comstr));
    end
else
    disp("gpuCore.a already exists. Delete it if you want to rebuild core libaries.");
end


