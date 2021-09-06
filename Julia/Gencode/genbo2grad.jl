#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function genbo2grad(foldername, filename, u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, natomtype)

# if (u[1]==0)
#     m = 0;
# else    
#     n,m = size(u);
#     if (n==0) || (n == 1)
#         m = 0;    
#     elseif (n == natomtype+1)
#         if u[n,1] == 0
#             m = 0;
#         end    
#     else
#         error("Number of functions for the two-body bond order potential must be equal to 1 plus the number of atom types");
#     end
# end

if length(u)>0
    if (u[1]==0)
        m = 0;
        n = 0;
    else
        n,m = size(u);
        # if (n==0) || (n == 1)
        #     m = 0;    
        # elseif (n == natomtype+1)
        #     if u[end,1] == 0
        #         m = 0;
        #     end    
        # else
        #     error("Number of functions for the two-body bond order potential must be equal to 1 plus the number of atom types");
        # end
    end
else
    n = 0;
    m = 0;
end

gen = 1;
ifile = 0;
stropu = "";
strcpu = "";
strgpu = "";
potnum = 0; 
fn = "PaircDensityGradient";

for i = 1:m
    ui = u[:,i];    
    if (ui[n] != 0)
        potnum = potnum + 1;        
        stropu1, strcpu1, strgpu1 = genpairgrad(foldername, filename*string(potnum), ui[1:(n-1)], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
        stropu1 = modifystr(stropu1, natomtype);
        strcpu1 = modifystr(strcpu1, natomtype);
        strgpu1 = modifystr(strgpu1, natomtype);
        strgpu1 = replace(strgpu1, "\t\ti += blockDim.x * gridDim.x;\n\t\t}" => "\t\t}\n\t\ti += blockDim.x * gridDim.x;");          
        stropu2, strcpu2, strgpu2 = gendensitygrad(foldername, fn*string(potnum), ui[n], rho, mu, eta, kappa, gen, ifile);
        stropu = stropu * stropu1 * "\n" * stropu2 * "\n"; 
        strcpu = strcpu * strcpu1 * "\n" * strcpu2 * "\n"; 
        strgpu = strgpu * strgpu1 * "\n" * strgpu2 * "\n"; 
    end    
end

if potnum==0
    ifile = 0;
    gen = 0;
    stropu, strcpu, strgpu = genpairgrad(foldername, filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    tmpopu, tmpcpu, tmpgpu = gendensitygrad(foldername, fn, u, mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  * "\n" * tmpopu * "\n";
    strcpu = strcpu  * "\n" * tmpcpu * "\n";
    strgpu = strgpu  * "\n" * tmpgpu * "\n";    
    cppfiles(foldername, filename, stropu, strcpu, strgpu, 1);
else
    sp0, sp1, sp2, sp3 = getpotstr(2);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    tmpopu = modstr11(tmpopu);
    tmpcpu = modstr11(tmpcpu);
    tmpgpu = modstr11(tmpgpu);
    stropu = stropu * tmpopu * "\n";
    strcpu = strcpu * tmpcpu * "\n";
    strgpu = strgpu * tmpgpu * "\n";    
    
    sp0, sp1, sp2, sp3 = getpotstr(0);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(fn, potnum, sp0, sp1, sp2, sp3);
    tmpopu = modstr12(tmpopu);
    tmpcpu = modstr12(tmpcpu);
    tmpgpu = modstr12(tmpgpu);    
    stropu = stropu * tmpopu;
    strcpu = strcpu * tmpcpu;
    strgpu = strgpu * tmpgpu;    
        
    cppfiles(foldername, filename, stropu, strcpu, strgpu, 0);        
end    

end

function modifystr(oldstr, natomtype)
    
newstr = oldstr;
newstr = replace(newstr, "u_xij[" => "\tu_xij[");
for i = 1:natomtype
    s = "u[" * string(i-1) * " * i*" * string(natomtype) * "]";    
    if contains(newstr, s)   
        news = "if (tj[i] == " * string(i) * ") \n\t\t\t" * "u[i]";
        if (i>1)
            news = "}\n \t\t" * news;
        end
        newstr = replace(newstr, s => news);    
    end
end
newstr = replace(newstr, "\t}\n}" => "\t\t}\n\t}\n}");    

return newstr 

end

function modstr11(tmpgpu)

tmpgpu = replace(tmpgpu, "(T *u," => "(T *u => T *u_xij,");   
tmpgpu = replace(tmpgpu, "(u," => "(u, u_xij,");   
tmpgpu = replace(tmpgpu, "(double *" => "(double *, double *");   
tmpgpu = replace(tmpgpu, "(float *" => "(float *, float *");   

return tmpgpu

end
        
function modstr12(tmpgpu)

tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_rho,");   
tmpgpu = replace(tmpgpu, "(u," => "(u, u_rho,");   
tmpgpu = replace(tmpgpu, "(double *" => "(double *, double *");   
tmpgpu = replace(tmpgpu, "(float *" => "(float *, float *");   

return tmpgpu

end
    