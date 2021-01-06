function genlib(cpucompiler::String, gpucompiler::String)
# Examples: genlib("g++","");
#           genlib("g++","nvcc");

if length(cpucompiler)==0
    error("cpucompiler is empty");
end

rm("cpuCore.a");
rm("libcpuCore.so");
str = `$cpucompiler -fPIC -O3 -c cpuCore.cpp`;
run(str);
str = `$cpucompiler --shared cpuCore.o -o libcpuCore.so`;
run(str);
str = `ar rvs cpuCore.a cpuCore.o`;
run(str);
rm("cpuCore.o");

if length(gpucompiler)>0
    rm("gpuCore.a");
    rm("libgpuCore.so");
    str = `$gpucompiler -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w gpuCore.cu`;
    run(str);
    str = `$cpucompiler --shared gpuCore.o -o libgpuCore.so`;
    run(str);
    str = `ar rvs gpuCore.a gpuCore.o`;
    run(str);
    rm("gpuCore.o");
end

end

