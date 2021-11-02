function buildpotlib(app)

[cpucompiler, gpucompiler] = getcompilers(app.cpucompiler, app.gpucompiler, app.cuda_arch);

disp("Build potential library for CPU platform using " + cpucompiler + " ...");

cplstr = cpucompiler + " -fPIC -ffast-math -O3 -c ../Potentials/cpuEmpiricalPotentials.cpp -o ../Potentials/cpuEmpiricalPotentials.o";
eval(char("!" + cplstr));    
if ismac==1
    comstr = cpucompiler + " --shared ../Potentials/cpuEmpiricalPotentials.o -o ../Potentials/cpuEmpiricalPotentials.dylib";
else
    comstr = cpucompiler + " --shared ../Potentials/cpuEmpiricalPotentials.o -o ../Potentials/cpuEmpiricalPotentials.so";
end
eval(char("!" + comstr));

[~,foundnvcc] = findexec("nvcc");
if (isempty(char(gpucompiler))==0) && (foundnvcc==1)   
    disp("Build potential library for GPU platform using " + gpucompiler + " ...");
    cuda_arch = app.cuda_arch;
    enzyme = app.enzyme;
    if contains(gpucompiler,"nvcc")
         gpuarch = " -arch=";
    else
         gpuarch = " --cuda-gpu-arch=";
    end  
    if (isempty(char(cuda_arch))==1)
        gpumacros = gpuarch + "sm_60 ";
    else
        gpumacros = gpuarch + cuda_arch + " ";
    end
    if contains(gpucompiler,"nvcc")
        gpuflags = " --compiler-options '-fPIC' ";
    else
        gpuflags = " -fPIC ";
    end
    if (isempty(char(enzyme))==1)
        enzyme = "";
        enzymemacros = "";
    else    
        enzymemacros = "-D _ENZYME";
    end
    cplstr = gpucompiler + " -std=c++11 " + gpuflags + gpumacros + enzymemacros + " -O3 -D_FORCE_INLINES -c -w ../Potentials/gpuEmpiricalPotentials.cu -o ../Potentials/gpuEmpiricalPotentials.o " + enzyme; 
    eval(char("!" + cplstr));
    if ismac==1
        comstr = gpucompiler + gpuflags + " --shared ../Potentials/gpuEmpiricalPotentials.o -o ../Potentials/gpuEmpiricalPotentials.dylib";
    else
        comstr = gpucompiler + gpuflags + " --shared ../Potentials/gpuEmpiricalPotentials.o -o ../Potentials/gpuEmpiricalPotentials.so";
    end
    eval(char("!" + comstr));
end


