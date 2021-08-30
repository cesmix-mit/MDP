function [cpucompiler, gpucompiler, gpumacros, gpuflags] = getcompilers(cpucompiler, gpucompiler, cuda_arch)

if (nargin<1) || (isempty(char(cpucompiler))==1)
    [clang,found] = findexec("clang++");
    if found==1
        cpucompiler = clang;
    else
        [gcc,found] = findexec("g++");
        if found==1
            cpucompiler = gcc;
        end
    end    
    if found==0
        error("MDP could not find C++ compiler. Please set cpucompiler to the path of a C++ compiler.");
    end
else
    [cpucompiler,found] = findexec(cpucompiler);
    if found==0
        [clang,found] = findexec("clang++");
        if found==1
            cpucompiler = clang;
        else
            [gcc,found] = findexec("g++");
            if found==1
                cpucompiler = gcc;
            end
        end    
        if found==0
            error("MDP could not find a C++ compiler. Please install it and set cpucompiler to the path of a C++ compiler.");
        end
    end
end

if (nargin<2) || (isempty(char(gpucompiler))==1)
    [clang,foundgpu] = findexec("clang++");
    if foundgpu==1
        gpucompiler = clang;
    else
        [nvcc,foundgpu] = findexec("nvcc");
        if foundgpu==1
            gpucompiler = nvcc;
        end
    end    
    [~,foundnvcc] = findexec("nvcc");
    if (foundgpu==0) || (foundnvcc==0)
        disp("CUDA Toolkit is not found on your system. MDP's CUDA code will not be compiled.");
        gpucompiler = ""; 
        gpumacros = "";
        gpuflags = "";
    end
else
    [gpucompiler,found] = findexec(gpucompiler);
    if found==0
        [clang,foundgpu] = findexec("clang++");
        if foundgpu==1
            gpucompiler = clang;
        else
            [nvcc,foundgpu] = findexec("nvcc");
            if foundgpu==1
                gpucompiler = nvcc;
            end
        end    
        [~,foundnvcc] = findexec("nvcc");
        if (foundgpu==0) || (foundnvcc==0)
            disp("CUDA Toolkit is not found on your system. MDP's CUDA code will not be compiled.");   
            gpucompiler = ""; 
            gpumacros = "";
            gpuflags = "";
        end        
    end    
end

if (isempty(char(gpucompiler))==0)
    if (nargin<3) || (isempty(char(cuda_arch))==1)
        if contains(gpucompiler,"nvcc")
            gpumacros = " -arch=sm_60 ";
        else
            gpumacros = " --cuda-gpu-arch=sm_60 ";
        end
    else
        if contains(gpucompiler,"nvcc")
            gpumacros = " -arch=" + cuda_arch + " ";
        else
            gpumacros = " --cuda-gpu-arch=" + cuda_arch + " ";
        end
    end
end

if (isempty(char(gpucompiler))==0)
    if contains(gpucompiler,"nvcc")
        gpuflags = " --compiler-options '-fPIC' ";
    else
        gpuflags = " -fPIC ";
    end
end

