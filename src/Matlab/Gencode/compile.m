function compilerstr = compile(app)

disp("Compiling C++ source code");

cpucompiler = app.cpucompiler;
gpucompiler = app.gpucompiler;
cpuflags = app.cpuflags;
gpuflags = app.gpuflags;
cpumacros = app.cpumacros;
gpumacros = app.gpumacros;

compilerstr = cell(4,1);
for i = 1:4
    compilerstr{i} = "";
end

% cpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
% cpumacros = "-std=c++11 -D _ENZYME -D _DEBUG";
% cpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so"
compilerstr{1} = cpucompiler + " " + cpumacros + " main.cpp -o cpuMDP -ffast-math -O3 -lblas -llapack " + cpuflags;
eval(char("!" + compilerstr{1}));

if ~isempty(char(gpucompiler))
    % gpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
    % gpumacros = "--cuda-gpu-arch=sm_60 -std=c++11 -D _ENZYME -D _DEBUG";
    % gpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so"
    compilerstr{2} = gpucompiler + " " + gpumacros + " -D_FORCE_INLINES -O3 -c -fPIC -w ../Core/gpuCore.cu -o ../Core/gpuCore.o";
    compilerstr{3} = gpucompiler + " " + gpumacros + " -D_FORCE_INLINES -O3 -c -fPIC -w ../Potentials/gpuEmpiricalPotentials.cu -o ../Core/gpuEmpiricalPotentials.o " + gpuflags; 
    compilerstr{4} = cpucompiler + " " + gpumacros + " -D _CUDA main.cpp -o gpuMDP ../Core/gpuCore.o ../Core/gpuEmpiricalPotentials.o -ffast-math -O3 -fno-vectorize -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack " + gpuflags;
    for i = 2:4
        eval(char("!" + compilerstr{i}));
    end
end















% compilerstr = cell(8,1);
% for i = 1:8
%     compilerstr{i} = "";
% end
% 
% if isempty(strfind(app.cpuflags, "Enzyme")) == 0
%     compilerstr{1} = cpucompiler + " -std=c++11 -D _ENZYME -D _DEBUG main.cpp -o cpuMDP -ffast-math -O3 " + cpuflags;
% else
%     compilerstr{1} = cpucompiler + " -std=c++11 -D _DEBUG main.cpp -o cpuMDP -ffast-math -O3 " + cpuflags;
%     %compilerstr{1} = cpucompiler + " -std=c++11 -D _CUDA main.cpp -o cpuMDP -ffast-math -O3 -lcudart -lcublas " + cpuflags;
% end
% for i = 1:1
%     eval(char("!" + compilerstr{i}));
% end
% 
% if ~isempty(gpucompiler)
%     %compilerstr{5} = gpucompiler + " -arch=sm_60 -std=c++11 -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w ../Core/gpuCore.cu -o ../Core/gpuCore.o";
%     compilerstr{5} = gpucompiler + " --cuda-gpu-arch=sm_60 -std=c++11 -D_FORCE_INLINES -O3 -c -fPIC -w ../Core/gpuCore.cu -o ../Core/gpuCore.o";
%     compilerstr{6} = gpucompiler + " --cuda-gpu-arch=sm_60 -std=c++11 -D_FORCE_INLINES -O3 -c -fPIC -w ../Potentials/gpuEmpiricalPotentials.cu -o ../Core/gpuEmpiricalPotentials.o" + " -Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
%     %compilerstr{8} = cpucompiler + " --cuda-gpu-arch=sm_60 -std=c++11 -D _DEBUG -D _CUDA main.cpp -o gpuMDP ../Core/gpuCore.o ../Core/gpuEmpiricalPotentials.o -ffast-math -O3 -lcudart -lcublas "  + cpuflags;
%      compilerstr{8} = cpucompiler + " --cuda-gpu-arch=sm_60 -std=c++11 -D _DEBUG -D _CUDA main.cpp -o gpuMDP ../Core/gpuCore.o ../Core/gpuEmpiricalPotentials.o -ffast-math -O2 -fno-vectorize -fno-unroll-loops -fPIC -lcudart -lcublas " + cpuflags;
%     for i = 5:8
%         eval(char("!" + compilerstr{i}));
%     end
% end
% 
% 
% delete('../Core/cpuCore.a');
% delete('../Core/libcpuCore.so');
% compilerstr{1} = cpucompiler + " -fPIC -std=c++11 -ffast-math -O3 -c ../Core/cpuCore.cpp -o ../Core/cpuCore.o";
% compilerstr{2} = cpucompiler + " --shared ../Core/cpuCore.o -o ../Core/libcpuCore.so";
% compilerstr{3} = "ar rvs ../Core/cpuCore.a ../Core/cpuCore.o";
% compilerstr{4} = cpucompiler + " -std=c++11 main.cpp -o cpuMDP ../Core/cpuCore.a -ffast-math -O3 " + cpuflags;
% for i = 1:4
%     eval(char("!" + compilerstr{i}));
% end
% delete('../Core/cpuCore.o');


%compilerstr{6} = cpucompiler + " --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so";
%compilerstr{7} = "ar rvs ../Core/gpuCore.a ../Core/gpuCore.o";
%compilerstr{8} = cpucompiler + " -std=c++11 -D _CUDA main.cpp -o gpuMDP ../Core/gpuCore.a ../Core/cpuCore.a -O3 -lblas -lcudart -lcublas" + gpuflags;
% g++ -fPIC -O3 -c ../Core/cpuCore.cpp -o ../Core/cpuCore.o
% g++ --shared ../Core/cpuCore.o -o ../Core/libcpuCore.so
% ar rvs ../Core/cpuCore.a ../Core/cpuCore.o
% g++ -std=c++11 main.cpp -o cpuSerialGOLFF ../Core/cpuCore.a -O3

% nvcc -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w ../Core/gpuCore.cu -o ../Core/gpuCore.o
% g++ --shared ../Core/gpuCore.o -o ../Core/libgpuCore.so
% ar rvs ../Core/gpuCore.a ../Core/gpuCore.o
% g++ -std=c++11 -D _CUDA main.cpp -o gpuSerialGOLFF ../Core/gpuCore.a ../Core/cpuCore.a -O3 -lcudart -lcublas


