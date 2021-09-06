#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function genquadrupletgrad(foldername,filename, u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile)

#foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

sp0, sp1, sp2, sp3 = getpotstr(4);

if gen==0
    stropu = "template <typename T> void " * opufile;

    tmp = sp0;
    stropu = stropu * tmp * "{\n";
    stropu = stropu * "}\n";
    stropu = replace(stropu, "(T *u, T *xij," => "(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij,");    

    tmp = "template void " * opufile;
    tmp = tmp * sp1;
    tmp = tmp * "template void " * opufile;
    tmp = tmp * sp2;
    tmp = replace(tmp, "(double *," => "(double *, double *, double *, double *,");            
    tmp = replace(tmp, "(float *," => "(float *, float *, float *, float *,");                

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
    if (contains(ustr, "xik"))
        mystr = varsassign(mystr, "xik", length(xik), ustr, 2);
    end
    if (contains(ustr, "xil"))
        mystr = varsassign(mystr, "xil", length(xil), ustr, 2);
    end
    if (contains(ustr, "qi"))
        mystr = varsassign(mystr, "qi", length(qi), ustr, 2);
    end
    if (contains(ustr, "qj"))
        mystr = varsassign(mystr, "qj", length(qj), ustr, 2);
    end
    if (contains(ustr, "qk"))
        mystr = varsassign(mystr, "qk", length(qk), ustr, 2);
    end
    if (contains(ustr, "ql"))
        mystr = varsassign(mystr, "ql", length(ql), ustr, 2);
    end
    if (contains(ustr, "ti"))
        mystr = varsassign(mystr, "ti", length(ti), ustr, 3);
    end
    if (contains(ustr, "tj"))
        mystr = varsassign(mystr, "tj", length(tj), ustr, 3);
    end
    if (contains(ustr, "tk"))
        mystr = varsassign(mystr, "tk", length(tk), ustr, 3);
    end
    if (contains(ustr, "tl"))
        mystr = varsassign(mystr, "tl", length(tl), ustr, 3);
    end
    if (contains(ustr, "ai"))
        mystr = varsassign(mystr, "ai", length(ai), ustr, 3);
    end
    if (contains(ustr, "aj"))
        mystr = varsassign(mystr, "aj", length(aj), ustr, 3);
    end
    if (contains(ustr, "ak"))
        mystr = varsassign(mystr, "ak", length(ak), ustr, 3);
    end
    if (contains(ustr, "al"))
        mystr = varsassign(mystr, "al", length(al), ustr, 3);
    end

    dim = length(xij);
    dudxij = [SymPy.symbols("dudxij$i") for i=1:dim];
    dudxik = [SymPy.symbols("dudxik$i") for i=1:dim];
    dudxil = [SymPy.symbols("dudxil$i") for i=1:dim];
    for i=1:dim
        dudxij[i] = sympy.diff(u[1], xij[i]);
        dudxik[i] = sympy.diff(u[1], xik[i]);
        dudxil[i] = sympy.diff(u[1], xil[i]);
    end
    dudx = [u[1]; dudxij[:]; dudxik[:]; dudxil[:]];
    mystr = derivassign(mystr, dudx, "dudx");   

    mystr = replace(mystr, "dudx[0 + i*" * string(length(dudx)) * "]" => "u[i]");       
    for d = 1:dim
        mystr = replace(mystr, "dudx[" * string(d) * " + i*" * string(length(dudx)) * "]" => "u_xij[" * string(d-1) * " + i*" * string(dim) * "]");     
        mystr = replace(mystr, "dudx[" * string(dim+d) * " + i*" * string(length(dudx)) * "]" => "u_xik[" * string(d-1) * " + i*" * string(dim) * "]");                                 
        mystr = replace(mystr, "dudx[" * string(2*dim+d) * " + i*" * string(length(dudx)) * "]" => "u_xil[" * string(d-1) * " + i*" * string(dim) * "]");                                 
    end    

    stropu = stropu * mystr * "\t}\n" * "}\n";
    stropu = replace(stropu, "(T *u, T *xij" => "(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij");

    strgpu = strgpu * mystr * "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu * "\t}\n" * "}\n\n";
    strgpu = replace(strgpu, "(T *u, T *xij" => "(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij");   

    tmp = "template <typename T> void " * gpufile;
    tmp = tmp * sp0;
    tmp = tmp * "{\n";
    tmp = tmp * "\tint blockDim = 256;\n";
    tmp = tmp * "\tint gridDim = (ng * blockDim - 1) / blockDim;\n";
    tmp = tmp * "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp * "\tkernel" * gpufile * "<<<gridDim, blockDim>>>" * sp3;
    tmp = tmp * "}\n";
    tmp = replace(tmp, "(T *u" => "(T *u, T *u_xij, T *u_xik, T *u_xil");    
    tmp = replace(tmp, "(u," => "(u, u_xij, u_xik, u_xil,");  
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
