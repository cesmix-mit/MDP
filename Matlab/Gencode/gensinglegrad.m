%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [stropu, strcpu, strgpu] = gensinglegrad(filename, u, xi, qi, ti, ai, mu, eta, kappa, gen, ifile)

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
    stropu = strrep(stropu, "(T *u, T *xi,", "(T *u, T *u_xi, T *xi,");            
    
    tmp = "template void " + opufile;
    tmp = tmp + sp1;
    tmp = tmp + "template void " + opufile;
    tmp = tmp + sp2;    
    tmp = strrep(tmp, "(double *,", "(double *, double *,");            
    tmp = strrep(tmp, "(float *,", "(float *, float *,");                
    
    stropu = stropu + tmp;
    strcpu = strrep(stropu, "opu", "cpu");
    strgpu = strrep(stropu, "opu", "gpu");       
    %strgpu = strgpu + "\n" + tmp;
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
    %mystr = symsassign(mystr, u);
    dudxi = 0*xi;
    for i = 1:length(xi)
        if (i==1)
            dudxi(i) = diff(u, 'xi1');            
        elseif (i==2)
            dudxi(i) = diff(u, 'xi2');            
        elseif (i==3)
            dudxi(i) = diff(u, 'xi3');                
        end        
    end
    dudx = [u; dudxi(:)];
    mystr = derivassign(mystr, dudx, 'dudx');
    mystr = strrep(mystr, "dudx[0 + i*" + num2str(length(dudx)) + "]", "u[i]");       
    for d = 1:dim
        mystr = strrep(mystr, "dudx[" + num2str(d) + " + i*" + num2str(length(dudx)) + "]", ...
                              "u_xi[" + num2str(d-1) + " + i*" + num2str(dim) + "]");       
    end
    
    stropu = stropu + mystr + "\t}\n" + "}\n";
    stropu = strrep(stropu, "(T *u", "(T *u, T *u_xi");
    
    strgpu = strgpu + mystr + "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu + "\t}\n" + "}\n\n";
    strgpu = strrep(strgpu, "(T *u", "(T *u, T *u_xi");
    
    tmp = "template <typename T> void " + gpufile;
    tmp = tmp + sp0;
    tmp = tmp + "{\n";
    tmp = tmp + "\tint blockDim = 256;\n";
    tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>" + sp3;
    tmp = tmp + "}\n";
    tmp = strrep(tmp, "(T *u", "(T *u, T *u_xi");    
    tmp = strrep(tmp, "(u,", "(u, u_xi,");  
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
