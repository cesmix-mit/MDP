%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [stropu, strcpu, strgpu] = genquadrupletgrad(filename, u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile)

foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

[sp0, sp1, sp2, sp3] = getpotstr(4);
dim = length(xij);

if gen==0
    stropu = "template <typename T> void " + opufile;

    tmp = sp0;

    stropu = stropu + tmp + "{\n";
    stropu = stropu + "}\n";
    stropu = strrep(stropu, "(T *u, T *xij,", "(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij,");    
    
    tmp = "template void " + opufile;
    tmp = tmp + sp1;
    tmp = tmp + "template void " + opufile;
    tmp = tmp + sp2;
    tmp = strrep(tmp, "(double *,", "(double *, double *, double *, double *,");            
    tmp = strrep(tmp, "(float *,", "(float *, float *, float *, float *,");                
    
    stropu = stropu + tmp;
    strcpu = strrep(stropu, "opu", "cpu");
    strgpu = strrep(stropu, "opu", "gpu");        
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
    if any(contains(ustr, "xij"))
        mystr = varsassign(mystr, "xij", length(xij), ustr, 2);
    end
    if any(contains(ustr, "xik"))
        mystr = varsassign(mystr, "xik", length(xik), ustr, 2);
    end
    if any(contains(ustr, "xil"))
        mystr = varsassign(mystr, "xil", length(xil), ustr, 2);
    end
    if any(contains(ustr, "qi"))
        mystr = varsassign(mystr, "qi", length(qi), ustr, 2);
    end
    if any(contains(ustr, "qj"))
        mystr = varsassign(mystr, "qj", length(qj), ustr, 2);
    end
    if any(contains(ustr, "qk"))
        mystr = varsassign(mystr, "qk", length(qk), ustr, 2);
    end
    if any(contains(ustr, "ql"))
        mystr = varsassign(mystr, "ql", length(ql), ustr, 2);
    end
    if any(contains(ustr, "ti"))
        mystr = varsassign(mystr, "ti", length(ti), ustr, 3);
    end
    if any(contains(ustr, "tj"))
        mystr = varsassign(mystr, "tj", length(tj), ustr, 3);
    end
    if any(contains(ustr, "tk"))
        mystr = varsassign(mystr, "tk", length(tk), ustr, 3);
    end
    if any(contains(ustr, "tl"))
        mystr = varsassign(mystr, "tl", length(tl), ustr, 3);
    end
    if any(contains(ustr, "ai"))
        mystr = varsassign(mystr, "ai", length(ai), ustr, 3);
    end
    if any(contains(ustr, "aj"))
        mystr = varsassign(mystr, "aj", length(aj), ustr, 3);
    end
    if any(contains(ustr, "ak"))
        mystr = varsassign(mystr, "ak", length(ak), ustr, 3);
    end
    if any(contains(ustr, "al"))
        mystr = varsassign(mystr, "al", length(al), ustr, 3);
    end
    %mystr = symsassign(mystr, u);
    dudxij = 0*xij; dudxik = 0*xik; dudxil = 0*xil; 
    for i = 1:length(xij)
        if (i==1)
            dudxij(i) = diff(u, 'xij1');          
            dudxik(i) = diff(u, 'xik1');  
            dudxil(i) = diff(u, 'xil1');  
        elseif (i==2)
            dudxij(i) = diff(u, 'xij2');          
            dudxik(i) = diff(u, 'xik2');     
            dudxil(i) = diff(u, 'xil2');  
        elseif (i==3)
            dudxij(i) = diff(u, 'xij3');        
            dudxik(i) = diff(u, 'xik3');     
            dudxil(i) = diff(u, 'xil3');  
        end        
    end
    dudx = [u; dudxij(:); dudxik(:); dudxil(:)];
    mystr = derivassign(mystr, dudx, 'dudx');
    mystr = strrep(mystr, "dudx[0 + i*" + num2str(length(dudx)) + "]", "u[i]");       
    for d = 1:dim
        mystr = strrep(mystr, "dudx[" + num2str(d) + " + i*" + num2str(length(dudx)) + "]", ...
                              "u_xij[" + num2str(d-1) + " + i*" + num2str(dim) + "]");     
        mystr = strrep(mystr, "dudx[" + num2str(dim+d) + " + i*" + num2str(length(dudx)) + "]", ...
                              "u_xik[" + num2str(d-1) + " + i*" + num2str(dim) + "]");                                 
        mystr = strrep(mystr, "dudx[" + num2str(2*dim+d) + " + i*" + num2str(length(dudx)) + "]", ...
                              "u_xil[" + num2str(d-1) + " + i*" + num2str(dim) + "]");                                                           
    end        
    
    stropu = stropu + mystr + "\t}\n" + "}\n\n";
    stropu = strrep(stropu, "(T *u, T *xij", "(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij");
    
    strgpu = strgpu + mystr + "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu + "\t}\n" + "}\n\n";
    strgpu = strrep(strgpu, "(T *u, T *xij", "(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij");   
    
    tmp = "template <typename T> void " + gpufile;
    tmp = tmp + sp0;
    tmp = tmp + "{\n";
    tmp = tmp + "\tint blockDim = 256;\n";
    tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>" + sp3;
    tmp = tmp + "}\n\n";
    tmp = strrep(tmp, "(T *u", "(T *u, T *u_xij, T *u_xik, T *u_xil");    
    tmp = strrep(tmp, "(u,", "(u, u_xij, u_xik, u_xil,");  
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
