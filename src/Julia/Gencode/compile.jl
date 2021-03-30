function compile(app)

print("compile code...\n");

codename = app.codename;
cpucompiler = app.cpucompiler;
mpicompiler = app.mpicompiler;
gpucompiler = app.gpucompiler;
cpuflags = app.cpuflags;
gpuflags = app.gpuflags;
sourcepath = app.sourcepath;
currentdir = app.currentdir;    

ii = findlast(sourcepath, currentdir);
tmp = findall("/", currentdir[(ii[end]):end]);
up = length(tmp);

codedir = "";
for i = 1:up
    codedir = codedir * "../";
end
maindir = codedir * "C++/Main/";
coredir = codedir * "C++/Core/";

compilerstr = Array{String, 1}(undef, 12);

compilerstr[1] = cpucompiler * " -fPIC -std=c++11 -ffast-math -O3 -c " * coredir * "cpuCore.cpp -o " * coredir * "cpuCore.o"; 
compilerstr[2] = cpucompiler * " --shared " * coredir * "cpuCore.o -o " * coredir * "libcpuCore.so";
compilerstr[3] = "ar rvs " * coredir * "cpuCore.a " * coredir * "cpuCore.o";

str1 = cpucompiler * " -std=c++11 " * maindir * "main.cpp -o " * maindir * "cpu" * codename * " " * coredir * "cpuCore.a ";
str2 = "-ffast-math -O3 " * app.cpuflags;
compilerstr[4] = str1 * str2;

for i = 1:4
    run(Gencode.string2cmd(compilerstr[i]));
end

# if length(gpucompiler)>0
#     str = `$gpucompiler -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w ../Core/gpuCore.cu -o ../Core/gpuCore.o`;
#     run(str);
#     str = `$cpucompiler --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so`;
#     run(str);
#     str = `ar rvs ../Core/gpuCore.a ../Core/gpuCore.o`;
#     run(str);
#     str = `$cpucompiler -std=c++11 -D _CUDA main.cpp -o gpuSerialGOLFF ../Core/gpuCore.a ../Core/cpuCore.a -O3 -lblas -lcudart -lcublas`;
#     run(str);
#     rm("../Core/gpuCore.o");
# end

end

