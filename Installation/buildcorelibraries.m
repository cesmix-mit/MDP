function [cpucompiler, gpucompiler, gpumacros, gpuflags] = buildcorelibraries(cpucompiler, gpucompiler, cuda_arch)

% g++ -fPIC -ffast-math -O3 -c cpuCore.cpp -o cpuCore.o
% ar rvs cpuCore.a cpuCore.o
% nvcc -std=c++11  --compiler-options '-fPIC'  -arch=sm_60  -D_FORCE_INLINES -O3 -c -w gpuCore.cu -o gpuCore.o
% ar rvs gpuCore.a gpuCore.o

if nargin<1
    cpucompiler = "";
end
if nargin<2
    gpucompiler = "";
end
if nargin<3
    cuda_arch = "";
end

disp("Building core libraries ...");

[cpucompiler, gpucompiler, gpumacros, gpuflags] = getcompilers(cpucompiler, gpucompiler, cuda_arch);

if exist("../C++/Core/cpuCore.o", "file")
    delete(char("../C++/Core/cpuCore.o"));
end
if exist("../C++/Core/cpuCore.a", "file")
    delete(char("../C++/Core/cpuCore.a"));
end
if exist("../C++/Core/cpuCore.dylib", "file")
    delete(char("../C++/Core/cpuCore.dylib"));
end
if exist("../C++/Core/cpuCore.so", "file")
    delete(char("../C++/Core/cpuCore.so"));
end
if exist("../C++/Core/gpuCore.o", "file")
    delete(char("../C++/Core/gpuCore.o"));
end
if exist("../C++/Core/gpuCore.a", "file")
    delete(char("../C++/Core/gpuCore.a"));
end
if exist("../C++/Core/gpuCore.dylib", "file")
    delete(char("../C++/Core/gpuCore.dylib"));
end
if exist("../C++/Core/gpuCore.so", "file")
    delete(char("../C++/Core/gpuCore.so"));
end

comstr = cpucompiler + " -fPIC -ffast-math -O3 -c ../C++/Core/cpuCore.cpp -o ../C++/Core/cpuCore.o";
eval(char("!" + comstr));
comstr = "ar rvs ../C++/Core/cpuCore.a ../C++/Core/cpuCore.o";
eval(char("!" + comstr));
if ismac==1
    comstr = cpucompiler + " --shared ../C++/Core/cpuCore.o -o ../C++/Core/cpuCore.dylib";
else
    comstr = cpucompiler + " --shared ../C++/Core/cpuCore.o -o ../C++/Core/cpuCore.so";
end
eval(char("!" + comstr));

if (isempty(char(gpucompiler))==0)
    comstr = gpucompiler + " -std=c++11 " + gpuflags + gpumacros + " -D_FORCE_INLINES -O3 -c -w ../C++/Core/gpuCore.cu -o ../C++/Core/gpuCore.o";
    eval(char("!" + comstr));
    comstr = "ar rvs ../C++/Core/gpuCore.a ../C++/Core/gpuCore.o";
    eval(char("!" + comstr));
    if ismac==1
        comstr = gpucompiler + gpuflags + "--shared ../C++/Core/gpuCore.o -o ../C++/Core/gpuCore.dylib";
    else
        comstr = gpucompiler + gpuflags + "--shared ../C++/Core/gpuCore.o -o ../C++/Core/gpuCore.so";
    end
    eval(char("!" + comstr));
end
