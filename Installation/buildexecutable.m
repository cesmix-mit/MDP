function [cpustr, gpustr] = buildexecutable(cpucompiler, gpucompiler, cuda_arch, enzyme)

% g++  -fPIC -ffast-math -O3 -c ../Potentials/cpuEmpiricalPotentials.cpp -o ../Potentials/cpuEmpiricalPotentials.o
% g++ --shared ../Potentials/cpuEmpiricalPotentials.o -o ../Potentials/cpuEmpiricalPotentials.so
% g++ -std=c++11 main.cpp -o cpuMDP ../Core/cpuCore.a -ffast-math -O3 -lblas -llapack  ../Potentials/cpuEmpiricalPotentials.so  
% nvcc -std=c++11  --compiler-options '-fPIC'  -arch=sm_70  -O3 -D_FORCE_INLINES -c -w ../Potentials/gpuEmpiricalPotentials.cu -o ../Potentials/gpuEmpiricalPotentials.o 
% nvcc --shared ../Potentials/gpuEmpiricalPotentials.o -o ../Potentials/gpuEmpiricalPotentials.so
% g++ -std=c++11  -D _CUDA main.cpp -o gpuMDP ../Core/cpuCore.a ../Core/gpuCore.a -ffast-math -O3 -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack  ../Potentials/cpuEmpiricalPotentials.so   ../Potentials/gpuEmpiricalPotentials.so 
% g++ -std=c++11  -D _CUDA -I/usr/tce/packages/cuda/cuda-10.1.243/include main.cpp -o gpuMDP ../Core/cpuCore.a ../Core/gpuCore.a -ffast-math -O3 -fno-unroll-loops -fPIC -lblas -llapack /usr/tce/packages/cuda/cuda-10.1.243/lib64/libcudart.so /usr/tce/packages/cuda/cuda-10.1.243/lib64/libcublas.so ../Potentials/cpuEmpiricalPotentials.so   ../Potentials/gpuEmpiricalPotentials.so

if nargin<1
    cpucompiler = "";
end
if nargin<2
    gpucompiler = "";
end
if nargin<3
    cuda_arch = "";
end
if nargin<4
    enzyme = "";
end

disp("Building executable ...");

[cpucompiler, gpucompiler, gpumacros, gpuflags] = getcompilers(cpucompiler, gpucompiler, cuda_arch);
if (nargin<4) || (isempty(char(enzyme))==1)
    enzyme = "";
    enzymemacros = "";
else    
    enzymemacros = "-D _ENZYME";
end

cdir = pwd(); 
cd(char("../C++/Main"));

cplstr = cpucompiler + " -fPIC -ffast-math -O3 -c ../Potentials/cpuEmpiricalPotentials.cpp -o ../Potentials/cpuEmpiricalPotentials.o";
eval(char("!" + cplstr));    
if ismac==1
    comstr = cpucompiler + " --shared ../Potentials/cpuEmpiricalPotentials.o -o ../Potentials/cpuEmpiricalPotentials.dylib";
    potlib = " ../Potentials/cpuEmpiricalPotentials.dylib ";
else
    comstr = cpucompiler + " --shared ../Potentials/cpuEmpiricalPotentials.o -o ../Potentials/cpuEmpiricalPotentials.so";
    potlib = " ../Potentials/cpuEmpiricalPotentials.so ";
end
eval(char("!" + comstr));
%cpustr = cpucompiler + " -std=c++11 " + enzymemacros + " main.cpp -o cpuMDP ../Core/cpuCore.a ../Potentials/cpuEmpiricalPotentials.o -ffast-math -O3 -lblas -llapack " + enzyme;
cpustr = cpucompiler + " -std=c++11 " + enzymemacros + " main.cpp -o cpuMDP ../Core/cpuCore.a -ffast-math -O3 -lblas -llapack " + potlib + " " + enzyme;
% if ismac 
%     cpustr = cpucompiler + " -std=c++11 " + enzymemacros + " main.cpp -o cpuMDP -ffast-math -O3 -lblas -llapack ../Core/cpuCore.dylib " + potlib + enzyme;
% else
%     cpustr = cpucompiler + " -std=c++11 " + enzymemacros + " main.cpp -o cpuMDP -ffast-math -O3 -lblas -llapack ../Core/cpuCore.so " + potlib + enzyme;
% end
eval(char("!" + cpustr));

if ~isempty(char(gpucompiler))
    cplstr = gpucompiler + " -std=c++11 " + gpuflags + gpumacros  + enzymemacros + " -O3 -D_FORCE_INLINES -c -w ../Potentials/gpuEmpiricalPotentials.cu -o ../Potentials/gpuEmpiricalPotentials.o " + enzyme; 
    eval(char("!" + cplstr));
    if ismac==1
        comstr = gpucompiler + gpuflags + "--shared ../Potentials/gpuEmpiricalPotentials.o -o ../Potentials/gpuEmpiricalPotentials.dylib";
        gpupotlib = " ../Potentials/gpuEmpiricalPotentials.dylib ";
    else
        comstr = gpucompiler + gpuflags + "--shared ../Potentials/gpuEmpiricalPotentials.o -o ../Potentials/gpuEmpiricalPotentials.so";
        gpupotlib = " ../Potentials/gpuEmpiricalPotentials.so ";
    end
    eval(char("!" + comstr));
    %gpustr = cpucompiler + " " + gpumacros + " " + enzymemacros + " -D _CUDA main.cpp -o gpuMDP ../Core/cpuCore.a ../Core/gpuCore.a ../Potentials/gpuEmpiricalPotentials.o -ffast-math -O3 -fno-vectorize -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack " + enzyme;
    gpustr = cpucompiler + " -std=c++11 " + enzymemacros + " -D _CUDA main.cpp -o gpuMDP ../Core/cpuCore.a ../Core/gpuCore.a -ffast-math -O3 -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack " + potlib + " " + gpupotlib + " " + enzyme;
%     if ismac
%         gpustr = cpucompiler + " -std=c++11 " + enzymemacros + " -D _CUDA main.cpp -o gpuMDP -ffast-math -O3 -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack ../Core/cpuCore.dylib ../Core/gpuCore.dylib" + potlib + " " + gpupotlib + " " + enzyme;
%     else
%         gpustr = cpucompiler + " -std=c++11 " + enzymemacros + " -D _CUDA main.cpp -o gpuMDP -ffast-math -O3 -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack ../Core/cpuCore.so ../Core/gpuCore.so" + potlib + " " + gpupotlib + " " + enzyme;
%     end
    eval(char("!" + gpustr));
else
    gpustr = "";
end

cd(cdir);
