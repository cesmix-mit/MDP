function [stropu, strcpu, strgpu] = genpair(filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile)

foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

[sp0, sp1, sp2, sp3] = getpotstr(2);

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
else
    stropu = "template <typename T> void " + opufile;
    strgpu = "template <typename T>  __global__  void kernel" + gpufile;

    tmp = sp0;

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
    if any(contains(ustr, "xij"))
        mystr = varsassign(mystr, "xij", length(xij), ustr, 2);
    end
    if any(contains(ustr, "qi"))
        mystr = varsassign(mystr, "qi", length(qi), ustr, 2);
    end
    if any(contains(ustr, "qj"))
        mystr = varsassign(mystr, "qj", length(qj), ustr, 2);
    end
    if any(contains(ustr, "ti"))
        mystr = varsassign(mystr, "ti", length(ti), ustr, 3);
    end
    if any(contains(ustr, "tj"))
        mystr = varsassign(mystr, "tj", length(tj), ustr, 3);
    end
    if any(contains(ustr, "ai"))
        mystr = varsassign(mystr, "ai", length(ai), ustr, 3);
    end
    if any(contains(ustr, "aj"))
        mystr = varsassign(mystr, "aj", length(aj), ustr, 3);
    end    
    mystr = symsassign(mystr, u);

    stropu = stropu + mystr + "\t}\n" + "}\n";

    strgpu = strgpu + mystr + "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu + "\t}\n" + "}\n\n";
    tmp = "template <typename T> void " + gpufile;
    tmp = tmp + sp0;
    tmp = tmp + "{\n";
    tmp = tmp + "\tint blockDim = 256;\n";
    tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>" + sp3;
    tmp = tmp + "}\n";
    strgpu = strgpu + tmp;

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
