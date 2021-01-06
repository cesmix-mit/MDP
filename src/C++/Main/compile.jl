function compile(cpucompiler::String, gpucompiler::String)
# Examples: compile("g++","");
#           compile("g++","nvcc");

if length(cpucompiler)==0
    error("cpucompiler is empty");
end

str = `$cpucompiler -fPIC -O3 -c ../Core/cpuCore.cpp -o ../Core/cpuCore.o`;
run(str);
str = `$cpucompiler --shared ../Core/cpuCore.o -o ../Core/libcpuCore.so`;
run(str);
str = `ar rvs ../Core/cpuCore.a ../Core/cpuCore.o`;
run(str);
str = `$cpucompiler -std=c++11 main.cpp -o cpuSerialGOLFF ../Core/cpuCore.a -O3`;
run(str);
rm("../Core/cpuCore.o");

if length(gpucompiler)>0
    str = `$gpucompiler -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w ../Core/gpuCore.cu -o ../Core/gpuCore.o`;
    run(str);
    str = `$cpucompiler --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so`;
    run(str);
    str = `ar rvs ../Core/gpuCore.a ../Core/gpuCore.o`;
    run(str);
    str = `$cpucompiler -std=c++11 -D _CUDA main.cpp -o gpuSerialGOLFF ../Core/gpuCore.a ../Core/cpuCore.a -O3 -lcudart -lcublas`;
    run(str);
    rm("../Core/gpuCore.o");
end

end

