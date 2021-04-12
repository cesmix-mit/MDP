function compilerstr = compile(app)

disp("Compiling C++ source code");

cpucompiler = app.cpucompiler;
cpuflags = app.cpuflags;
gpucompiler = app.gpucompiler;
gpuflags = app.gpuflags;

if nargin<1
    cpucompiler = "g++";
end
if nargin<2
    gpucompiler = [];
end

compilerstr = cell(8,1);
for i = 1:8
    compilerstr{i} = "";
end

compilerstr{1} = cpucompiler + " -std=c++11 main.cpp -o cpuMDP -ffast-math -O3 " + cpuflags;
for i = 1:1
    eval(char("!" + compilerstr{i}));
end

% delete('../Core/cpuCore.a');
% delete('../Core/libcpuCore.so');
% compilerstr{1} = cpucompiler + " -fPIC -std=c++11 -ffast-math -O3 -c ../Core/cpuCore.cpp -o ../Core/cpuCore.o";
% compilerstr{2} = cpucompiler + " --shared ../Core/cpuCore.o -o ../Core/libcpuCore.so";
% compilerstr{3} = "ar rvs ../Core/cpuCore.a ../Core/cpuCore.o";
% compilerstr{4} = cpucompiler + " -std=c++11 main.cpp -o cpuMDP ../Core/cpuCore.a -ffast-math -O3 -lm" + cpuflags;
% for i = 1:4
%     eval(char("!" + compilerstr{i}));
% end
% delete('../Core/cpuCore.o');

if ~isempty(gpucompiler)
    delete('../Core/gpuCore.a');
    delete('../Core/libgpuCore.so');
    compilerstr{5} = gpucompiler + " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w ../Core/gpuCore.cu -o ../Core/gpuCore.o";
    compilerstr{6} = cpucompiler + " --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so";
    compilerstr{7} = "ar rvs ../Core/gpuCore.a ../Core/gpuCore.o";
    compilerstr{8} = cpucompiler + " -std=c++11 -D _CUDA main.cpp -o gpuMDP ../Core/gpuCore.a ../Core/cpuCore.a -O3 -lblas -lcudart -lcublas" + gpuflags;
    for i = 5:8
        eval(char("!" + compilerstr{i}));
    end
    delete('../Core/gpuCore.o');
end

% g++ -fPIC -O3 -c ../Core/cpuCore.cpp -o ../Core/cpuCore.o
% g++ --shared ../Core/cpuCore.o -o ../Core/libcpuCore.so
% ar rvs ../Core/cpuCore.a ../Core/cpuCore.o
% g++ -std=c++11 main.cpp -o cpuSerialGOLFF ../Core/cpuCore.a -O3

% nvcc -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w ../Core/gpuCore.cu -o ../Core/gpuCore.o
% g++ --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so
% ar rvs ../Core/gpuCore.a ../Core/gpuCore.o
% g++ -std=c++11 -D _CUDA main.cpp -o gpuSerialGOLFF ../Core/gpuCore.a ../Core/cpuCore.a -O3 -lcudart -lcublas


