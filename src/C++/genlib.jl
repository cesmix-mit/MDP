function genlib(cpucompiler::String, gpucompiler::String)
# Examples: genlib("g++","");
#           genlib("g++","nvcc");

if length(cpucompiler)==0
    error("cpucompiler is empty");
end

str = `$cpucompiler -fPIC -O3 -c opuGolff.cpp`;
run(str);
str = `$cpucompiler --shared cpuGolff.o -o libcpuGolff.so`;
run(str);
str = `ar rvs opuGolff.a opuGolff.o`;
run(str);

end

