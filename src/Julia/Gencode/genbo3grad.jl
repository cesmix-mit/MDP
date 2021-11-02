#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function genbo3grad(foldername, filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, natomtype)

# if (u[1]==0)
#     m = 0;
# else
#     n,m = size(u);
#     if (n==0) || (n == 1)
#         m = 0;    
#     elseif (n == natomtype+2)
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
fp = "TripletcPairGradient";
fn = "TripletcDensityGradient";

for i = 1:m
    ui = u[:,i];    
    if ui[n] != 0
        potnum = potnum + 1;        
        stropu0, strcpu0, strgpu0 = gentripletgrad(foldername, filename * string(potnum), ui[1:(n-2)], xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);        
        stropu1, strcpu1, strgpu1 = genpairgrad(foldername, fp * string(potnum), ui[n-1], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);        
        stropu2, strcpu2, strgpu2 = gendensitygrad(foldername, fn * string(potnum), ui[end], rho, mu, eta, kappa, gen, ifile);
        stropu = stropu * stropu0 * "\n" * stropu1 * "\n" * stropu2 * "\n"; 
        strcpu = strcpu * strcpu0 * "\n" * strcpu1 * "\n" * strcpu2 * "\n"; 
        strgpu = strgpu * strgpu0 * "\n" * strgpu1 * "\n" * strgpu2 * "\n"; 
    end    
end

if potnum==0
    ifile = 0;
    gen = 0;
    stropu, strcpu, strgpu = gentripletgrad(foldername, filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
    tmqopu, tmqcpu, tmqgpu = genpairgrad(foldername, fp, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    tmpopu, tmpcpu, tmpgpu = gendensitygrad(foldername, fn, u, mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  * "\n" * tmqopu * "\n" * tmpopu * "\n";
    strcpu = strcpu  * "\n" * tmqcpu * "\n" * tmpcpu * "\n";
    strgpu = strgpu  * "\n" * tmqgpu * "\n" * tmpgpu * "\n";    
    cppfiles(foldername, filename, stropu, strcpu, strgpu, 1);
else
    sp0, sp1, sp2, sp3 = getpotstr(3);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    tmpopu = modstr0(tmpopu);
    tmpcpu = modstr0(tmpcpu);
    tmpgpu = modstr0(tmpgpu);        
    stropu = stropu * tmpopu * "\n";
    strcpu = strcpu * tmpcpu * "\n";
    strgpu = strgpu * tmpgpu * "\n";    

    sp0, sp1, sp2, sp3 = getpotstr(2);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(fp, potnum, sp0, sp1, sp2, sp3);
    tmpopu = modstr1(tmpopu);
    tmpcpu = modstr1(tmpcpu);
    tmpgpu = modstr1(tmpgpu);
    stropu = stropu * tmpopu * "\n";
    strcpu = strcpu * tmpcpu * "\n";
    strgpu = strgpu * tmpgpu * "\n";    

    sp0, sp1, sp2, sp3 = getpotstr(0);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(fn, potnum, sp0, sp1, sp2, sp3);
    tmpopu = modstr2(tmpopu);
    tmpcpu = modstr2(tmpcpu);
    tmpgpu = modstr2(tmpgpu);        
    stropu = stropu * tmpopu;
    strcpu = strcpu * tmpcpu;
    strgpu = strgpu * tmpgpu;    

    cppfiles(foldername, filename, stropu, strcpu, strgpu, 0);
end    

end

function modstr0(tmpgpu)

tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_xij, T *u_xik,");   
tmpgpu = replace(tmpgpu, "(u," => "(u, u_xij, u_xik,");   
tmpgpu = replace(tmpgpu, "(double *" => "(double *, double *, double *");   
tmpgpu = replace(tmpgpu, "(float *" => "(float *, float *, float *");   
return tmpgpu

end

function modstr1(tmpgpu)

tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_xij,");   
tmpgpu = replace(tmpgpu, "(u," => "(u, u_xij,");   
tmpgpu = replace(tmpgpu, "(double *" => "(double *, double *");   
tmpgpu = replace(tmpgpu, "(float *" => "(float *, float *");   
return tmpgpu

end

function modstr2(tmpgpu)

tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_rho,");   
tmpgpu = replace(tmpgpu, "(u," => "(u, u_rho,");   
tmpgpu = replace(tmpgpu, "(double *" => "(double *, double *");   
tmpgpu = replace(tmpgpu, "(float *" => "(float *, float *");   
return tmpgpu

end


# function modifystr3(oldstr, natomtype)
    
# newstr = oldstr;
# for i = 1:natomtype
#     s = "u[" * string(i-1) * " * i*" * string(natomtype) * "]";    
#     if contains(newstr, s)   
#         news = "if (tk[i] == " * string(i) * ") \n\t\t\t" * "u[i]";
#         newstr = replace(newstr, s, news);    
#     end
# end

# return newstr

# end


