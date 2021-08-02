#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function single(filename::String, u, xi, qi, ti, ai, mu, gen)

foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

if gen==0
    stropu = "template <typename T> void " * opufile;

    tmp = "(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, int dim, int ncq, int nmu, int ng)\n";
    
    stropu = stropu * tmp * "{\n";
    stropu = stropu * "}\n\n";
    
    tmp = "template void " * opufile;
    tmp = tmp * "(double *, double *, double *, int *, int *, double *, int, int, int, int);\n";
    tmp = tmp * "template void " * opufile;
    tmp = tmp * "(float *, float *, float *, int *, int *, float *, int, int, int, int);\n";
    
    stropu = stropu * tmp;
    strcpu = replace(stropu, "opu" => "cpu");
    strgpu = replace(stropu, "opu" => "gpu");
    
    ioopu = open(foldername * "/" * opufile * ".cpp", "w");
    write(ioopu, stropu);
    close(ioopu);
    
    iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
    write(iocpu, strcpu);
    close(iocpu);
    
    iogpu = open(foldername * "/" * gpufile * ".cu", "w");
    write(iogpu, strgpu);
    close(iogpu);
else    
    stropu = "template <typename T> void " * opufile;
    strgpu = "template <typename T>  __global__  void kernel" * gpufile;

    tmp = "(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, int dim, int ncq, int nmu, int ng)\n";

    stropu = stropu * tmp * "{\n";
    stropu = stropu * "\tfor (int i = 0; i <ng; i++) {\n";

    strgpu = strgpu * tmp * "{\n";
    strgpu = strgpu * "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
    strgpu = strgpu * "\twhile (i<ng) {\n";

    str = "";
    str = varsassign(str, "mu", length(mu), 0);
    str = varsassign(str, "xi", length(xi), 2);
    str = varsassign(str, "qi", length(qi), 2);
    str = varsassign(str, "ti", length(ti), 2);
    str = varsassign(str, "ai", length(ai), 2);
    str = sympyassign(str, u);

    stropu = stropu * str * "\t}\n" * "}\n\n";
    tmp = "template void " * opufile;
    tmp = tmp * "(double *, double *, double *, int *, int *, double *, int, int, int, int);\n";
    tmp = tmp * "template void " * opufile;
    tmp = tmp * "(float *, float *, float *, int *, int *, float *, int, int, int, int);\n";
    stropu = stropu * tmp;

    ioopu = open(foldername * "/" * opufile * ".cpp", "w");
    write(ioopu, stropu);
    close(ioopu);

    strgpu = strgpu * str * "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu * "\t}\n" * "}\n\n";
    tmp = "template <typename T> void " * gpufile;
    tmp = tmp * "(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, int dim, int ncq, int nmu, int ng)\n";
    tmp = tmp * "{\n";
    tmp = tmp * "\tint blockDim = 256;\n";
    tmp = tmp * "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp * "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp * "\tkernel" * gpufile * "<<<gridDim, blockDim>>>(u, xi, qi, ti, ai, mu, dim, ncq, nmu, ng);\n";
    tmp = tmp * "}\n\n";
    tmp = tmp * "template void " * gpufile;
    tmp = tmp * "(double *, double *, double *, int *, int *, double *, int, int, int, int);\n";
    tmp = tmp * "template void " * gpufile;
    tmp = tmp * "(float *, float *, float *, int *, int *, float *, int, int, int, int);\n";
    strgpu = strgpu * tmp;

    iogpu = open(foldername * "/" * gpufile * ".cu", "w");
    write(iogpu, strgpu);
    close(iogpu);

    strcpu = replace(stropu, "opu" => "cpu");
    strcpu = replace(strcpu, "for (int i = 0; i <ng; i++) {" => "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

    iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
    write(iocpu, strcpu);
    close(iocpu);
end

end
