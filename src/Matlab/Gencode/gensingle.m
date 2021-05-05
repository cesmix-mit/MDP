function [stropu, strcpu, strgpu] = gensingle(filename, u, xi, qi, ti, ai, mu, eta, kappa, gen, ifile)

foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

[sp0, sp1, sp2, sp3] = getpotstr(1);

if gen==0
    stropu = "template <typename T> void " + opufile;

    tmp = sp0;

    stropu = stropu + tmp + "{\n";
    stropu = stropu + "}\n";

    tmp = "template void " + opufile;
    tmp = tmp + sp1;
    tmp = tmp + "template void " + opufile;
    tmp = tmp + sp2;

    stropu = stropu + tmp;
    strcpu = strrep(stropu, "opu", "cpu");
    strgpu = strrep(stropu, "opu", "gpu");
    
    tmp = strgpu;    
    tmp = strrep(tmp, "(T *u, T *xi,", "Gradient(T *u, T *du, T *u_xi, T *xi,");            
    tmp = strrep(tmp, "(double *,", "Gradient(double *, double *, double *,");            
    tmp = strrep(tmp, "(float *,", "Gradient(float *, float *, float *,");            
    strgpu = strgpu + "\n" + tmp;
else
    stropu = "template <typename T> void " + opufile;
    strgpu = "template <typename T>  __global__  void kernel" + gpufile;

    % here
    tmp = sp0;
    st2 = strgpu + sp0 + "{\n";    

    stropu = stropu + tmp + "{\n";
    stropu = stropu + "\tfor (int i = 0; i <ng; i++) {\n";

    strgpu = strgpu + tmp + "{\n";
    strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
    strgpu = strgpu + "\twhile (i<ng) {\n";

    ustr = string(u);            
    mystr = "";
    if any(contains(ustr, "mu")) 
        mystr = varsassign(mystr, "mu", length(mu), ustr, 0);
    end
    if any(contains(ustr, "eta")) 
        mystr = varsassign(mystr, "eta", length(eta), ustr, 0);
    end    
    if any(contains(ustr, "kappa")) 
        mystr = varsassign(mystr, "kappa", length(kappa), ustr, 1);    
    end 
    if any(contains(ustr, "xi"))
        mystr = varsassign(mystr, "xi", length(xi), ustr, 2);
    end
    if any(contains(ustr, "qi"))
        mystr = varsassign(mystr, "qi", length(qi), ustr, 2);
    end
    if any(contains(ustr, "ti"))
        mystr = varsassign(mystr, "ti", length(ti), ustr, 3);
    end
    if any(contains(ustr, "ai"))
        mystr = varsassign(mystr, "ai", length(ai), ustr, 3);
    end    
    mystr = symsassign(mystr, u);

    stropu = stropu + mystr + "\t}\n" + "}\n";

    strgpu = strgpu + mystr + "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu + "\t}\n" + "}\n\n";
    
    % here
    st1 = strgpu;    
    
    tmp = "template <typename T> void " + gpufile;
    tmp = tmp + sp0;
    tmp = tmp + "{\n";
    tmp = tmp + "\tint blockDim = 256;\n";
    tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>" + sp3;
    tmp = tmp + "}\n";
    strgpu = strgpu + tmp;

    % here
    st1 = strrep(st1, "__global__  void kernelgpu", "__device__  void devicegpu");     
    st1 = strrep(st1, "T *", "T *__restrict__ ");      
    st1 = strrep(st1, "int *", "int *__restrict__ ");      
    st2 = strrep(st2, "(T *u, T *xij,", "Gradient(T *u, T *du, T *u_xi, T *xi,");                    
    st2 = strrep(st2, "*", "*__restrict__ ");      
    st3 = "\t__enzyme_autodiff((void*)devicegpu" + filename +  "<T>, \n";
    st3 = st3 +  "\t\tenzyme_dup, u, du, \n";
    st3 = st3 +  "\t\tenzyme_dup, xi, u_xi, \n";
    st3 = st3 +  "\t\tenzyme_const, qi, \n";
    st3 = st3 +  "\t\tenzyme_const, ti, \n";
    st3 = st3 +  "\t\tenzyme_const, ai, \n";
    st3 = st3 +  "\t\tenzyme_const, mu, \n";
    st3 = st3 +  "\t\tenzyme_const, eta, \n";
    st3 = st3 +  "\t\tenzyme_const, kappa, \n";
    st3 = st3 +  "\t\tdim, ncq, nmu, neta, nkappa, ng); \n";      
    st2 = st2 + st3 + "}\n";
    st4 = tmp;
    st4 = strrep(st4, "(T *u, T *xij,", "Gradient(T *u, T *du, T *u_xi, T *xi,");   
    st4 = strrep(st4, "u, xij,", "u, du, u_xi, xi,");   
    st4 = strrep(st4, "<<<", "Gradient<<<");   
    st2 = st2 + "\n" + st4;
    strgpu = strgpu + "\n" + st1 + "\n" + st2;    
    
    strcpu = strrep(stropu, 'opu', "cpu");
    strcpu = strrep(strcpu, "for (int i = 0; i <ng; i++) {", "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");
end

if ifile==1
    fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
    fprintf(fid, char(stropu));
    fclose(fid);

    fid = fopen(foldername + "/" + gpufile + ".cu", 'w');
    fprintf(fid, char(strgpu));
    fclose(fid);

    fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
    fprintf(fid, char(strcpu));
    fclose(fid);    
end

end
