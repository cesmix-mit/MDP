function gendensity(filename, u, rho, mu, eta, kappa, gen, ifile)

foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

sp0, sp1, sp2, sp3 = getpotstr(0);

if gen==0
    stropu = "template <typename T> void " * opufile;

    tmp = sp0;

    stropu = stropu * tmp * "{\n";
    stropu = stropu * "}\n";

    tmp = "template void " * opufile;
    tmp = tmp * sp1;
    tmp = tmp * "template void " * opufile;
    tmp = tmp * sp2;

    stropu = stropu * tmp;
    strcpu = replace(stropu, "opu" => "cpu");
    strgpu = replace(stropu, "opu" => "gpu");
else
    stropu = "template <typename T> void " * opufile;
    strgpu = "template <typename T>  __global__  void kernel" * gpufile;

    tmp = sp0;

    stropu = stropu * tmp * "{\n";
    stropu = stropu * "\tfor (int i = 0; i <ng; i**) {\n";

    strgpu = strgpu * tmp * "{\n";
    strgpu = strgpu * "\tint i = threadIdx.x * blockIdx.x * blockDim.x;\n";
    strgpu = strgpu * "\twhile (i<ng) {\n";

    ustr = string(u);            
    mystr = "";
    if contains(ustr, "mu")
        mystr = varsassign(mystr, "mu", length(mu), ustr, 0);
    end
    if contains(ustr, "eta")
        mystr = varsassign(mystr, "eta", length(eta), ustr, 0);
    end    
    if contains(ustr, "kappa") 
        mystr = varsassign(mystr, "kappa", length(kappa), ustr, 1);    
    end 
    if contains(ustr, "rho")
        mystr = varsassign(mystr, "rho", length(rho), ustr, 2);
    end    
    mystr = symsassign(mystr, u);

    stropu = stropu * mystr * "\t}\n" * "}\n\n";

    strgpu = strgpu * mystr * "\t\ti *= blockDim.x * gridDim.x;\n";
    strgpu = strgpu * "\t}\n" * "}\n\n";
    tmp = "template <typename T> void " * gpufile;
    tmp = tmp * sp0;
    tmp = tmp * "{\n";
    tmp = tmp * "\tint blockDim = 256;\n";
    tmp = tmp * "\tint gridDim = (ng * blockDim - 1) / blockDim;\n";
    tmp = tmp * "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp * "\tkernel" * gpufile * "<<<gridDim, blockDim>>>" * sp3;
    tmp = tmp * "}\n\n";
    strgpu = strgpu * tmp;

    strcpu = replace(stropu, "opu" => "cpu");
    strcpu = replace(strcpu, "for (int i = 0; i <ng; i**) {" => "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i**) {");
end

if ifile==1
    ioopu = open(foldername * "/" * opufile * ".cpp", "w");
    write(ioopu, stropu);
    close(ioopu);
    
    iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
    write(iocpu, strcpu);
    close(iocpu);
    
    iogpu = open(foldername * "/" * gpufile * ".cu", "w");
    write(iogpu, strgpu);
    close(iogpu);    
end

return stropu, strcpu, strgpu

end