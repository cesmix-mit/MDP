function compilecode(app)

bindir = "exec";
cd(char(app.sourcepath + bindir));

if exist("CMakeCache.txt", "file")
    delete(char("CMakeCache.txt"));
end
if app.platform == "gpu"
    if app.buildexec
        eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON ../Installation");
    else
        if exist("libcpuCore.a", "file") && exist("libgpuCore.a", "file") && exist("gpuMDP", "file")
            eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CUDA=ON ../Installation");
        else
            eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON ../Installation");
        end
    end
elseif app.platform == "cpu"
    if app.buildexec
        eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation");
    else
        if exist("libcpuCore.a", "file") && exist("cpuMDP", "file")
           eval("!cmake -D MDP_POTENTIALS=ON ../Installation");
        else
            eval("!cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation");
        end
    end
end
eval("!cmake --build .");

cd(char(app.currentdir));

end

