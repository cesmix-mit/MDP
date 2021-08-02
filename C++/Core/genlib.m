function genlib(cpucompiler, gpucompiler)

if nargin<1
    cpucompiler = "g++";
end
if nargin<2
    gpucompiler = [];
end

delete('cpuCore.a');
delete('libcpuCore.so');
compilerstr = cell(6,1);
compilerstr{1} = cpucompiler + " -fPIC -O3 -c cpuCore.cpp";
compilerstr{2} = cpucompiler + " --shared cpuCore.o -o libcpuCore.so";
compilerstr{3} = "ar rvs cpuCore.a cpuCore.o";
for i = 1:3
    eval(char("!" + compilerstr{i}));
end
delete('cpuCore.o');

if ~isempty(gpucompiler)
    delete('gpuCore.a');
    delete('libgpuCore.so');
    compilerstr{4} = gpucompiler + " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w gpuCore.cu";
    compilerstr{5} = cpucompiler + " --shared gpuCore.o -o libgpuCore.so";
    compilerstr{6} = "ar rvs gpuCore.a gpuCore.o";
    for i = 4:6
        eval(char("!" + compilerstr{i}));
    end
    delete('gpuCore.o');
end


% !g++ -fPIC -O3 -c cpuCore.cpp
% !g++ --shared cpuCore.o -o libcpuCore.so
% !ar rvs cpuCore.a cpuCore.o        
% 
% !nvcc -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w gpuCore.cu
% !g++ --shared gpuCore.o -o libgpuCore.so
% !ar rvs gpuCore.a gpuCore.o        


