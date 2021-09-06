#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function genpairgrad(foldername,filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile)

#foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

sp0, sp1, sp2, sp3 = getpotstr(2);

if gen==0
    stropu = "template <typename T> void " * opufile;

    tmp = sp0;

    stropu = stropu * tmp * "{\n";
    stropu = stropu * "}\n";
    stropu = replace(stropu, "(T *u, T *xij," => "(T *u, T *u_xij, T *xij,");        

    tmp = "template void " * opufile;
    tmp = tmp * sp1;
    tmp = tmp * "template void " * opufile;
    tmp = tmp * sp2;
    tmp = replace(tmp, "(double *," => "(double *, double *,");            
    tmp = replace(tmp, "(float *," => "(float *, float *,");                

    stropu = stropu * tmp;
    strcpu = replace(stropu, "opu" => "cpu");
    strgpu = replace(stropu, "opu" => "gpu");    
else
    stropu = "template <typename T> void " * opufile;
    strgpu = "template <typename T>  __global__  void kernel" * gpufile;

    # here
    tmp = sp0;
    st2 = strgpu * sp0 * "{\n";    

    stropu = stropu * tmp * "{\n";
    stropu = stropu * "\tfor (int i = 0; i <ng; i++) {\n";

    strgpu = strgpu * tmp * "{\n";
    strgpu = strgpu * "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
    strgpu = strgpu * "\twhile (i<ng) {\n";
        
    ustr = string(u);         
    mystr = "";
    if (contains(ustr, "mu")) 
        mystr = varsassign(mystr, "mu", length(mu), ustr, 0);
    end
    if (contains(ustr, "eta")) 
        mystr = varsassign(mystr, "eta", length(eta), ustr, 0);
    end    
    if (contains(ustr, "kappa")) 
        mystr = varsassign(mystr, "kappa", length(kappa), ustr, 1);    
    end    
    if (contains(ustr, "xij"))
        mystr = varsassign(mystr, "xij", length(xij), ustr, 2);
    end
    if (contains(ustr, "qi"))
        mystr = varsassign(mystr, "qi", length(qi), ustr, 2);
    end
    if (contains(ustr, "qj"))
        mystr = varsassign(mystr, "qj", length(qj), ustr, 2);
    end
    if (contains(ustr, "ti"))
        mystr = varsassign(mystr, "ti", length(ti), ustr, 3);
    end
    if (contains(ustr, "tj"))
        mystr = varsassign(mystr, "tj", length(tj), ustr, 3);
    end
    if (contains(ustr, "ai"))
        mystr = varsassign(mystr, "ai", length(ai), ustr, 3);
    end
    if (contains(ustr, "aj"))
        mystr = varsassign(mystr, "aj", length(aj), ustr, 3);
    end    
    dim = length(xij);
    m = length(u); 
    dudx = [];   
    for j=1:m
        dudxij = [SymPy.symbols("dudx$i") for i=1:dim];
        for i=1:dim
            dudxij[i] = sympy.diff(u[m], xij[i]);
        end
        dudx = [dudx; u[m]; dudxij[:]];
    end
    mystr = derivassign(mystr, dudx, "dudx");   
    for j=1:m
        mystr = replace(mystr, "dudx[" * string((j-1)*(dim+1)) * " + i*" * string(length(dudx)) * "]" => "u[" * string(j-1) * " + i*" * string(m) * "]");       
        for d = 1:dim
            mystr = replace(mystr, "dudx[" * string((j-1)*(dim+1)+d) * " + i*" * string(length(dudx)) * "]" => "u_xij[" * string(d-1) * " + i*" * string(dim) * "]");       
        end
    end

    stropu = stropu * mystr * "\t}\n" * "}\n";
    stropu = replace(stropu, "(T *u, T *xij" => "(T *u, T *u_xij, T *xij");

    strgpu = strgpu * mystr * "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu * "\t}\n" * "}\n\n";
    strgpu = replace(strgpu, "(T *u, T *xij" => "(T *u, T *u_xij, T *xij");

    tmp = "template <typename T> void " * gpufile;
    tmp = tmp * sp0;
    tmp = tmp * "{\n";
    tmp = tmp * "\tint blockDim = 256;\n";
    tmp = tmp * "\tint gridDim = (ng * blockDim - 1) / blockDim;\n";
    tmp = tmp * "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp * "\tkernel" * gpufile * "<<<gridDim, blockDim>>>" * sp3;
    tmp = tmp * "}\n";
    tmp = replace(tmp, "(T *u" => "(T *u, T *u_xij");    
    tmp = replace(tmp, "(u," => "(u, u_xij,");  
    strgpu = strgpu * tmp;

    strcpu = replace(stropu, "opu" => "cpu");
    strcpu = replace(strcpu, "for (int i = 0; i <ng; i++) {" => "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");
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
