#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function genpotentialgrad(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot)

n = length(u);
gen = 1;
ifile = 0;
stropu = "";
strcpu = "";
strgpu = "";
potnum = 0; 

for i = 1:n
    ui = u[i];
    if ui != 0        
        potnum = potnum + 1;
        fn = filename * string(potnum);
        if pot==1
            tmpopu, tmpcpu, tmpgpu = gensinglegrad(foldername, fn, ui, xi, qi, ti, ai, mu, eta, kappa, gen, ifile);
        elseif pot==2
            tmpopu, tmpcpu, tmpgpu = genpairgrad(foldername, fn, ui, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
        elseif pot==3
            tmpopu, tmpcpu, tmpgpu = gentripletgrad(foldername, fn, ui, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
        elseif pot==4    
            tmpopu, tmpcpu, tmpgpu = genquadrupletgrad(foldername, fn, ui, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile);
        end
        tmpopu = modifystr1(tmpopu);
        tmpcpu = modifystr1(tmpcpu);
        tmpgpu = modifystr1(tmpgpu);
        stropu = stropu * tmpopu * "\n";
        strcpu = strcpu * tmpcpu * "\n";
        strgpu = strgpu * tmpgpu * "\n";        
    end    
end    

if potnum==0
    nopotentialgrad(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot);
else
    sp0, sp1, sp2, sp3 = getpotstr(pot);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    # stropu = stropu * tmpopu;
    # strcpu = strcpu * tmpcpu;
    # strgpu = strgpu * tmpgpu;        
    if pot==1
        tmpopu = replace(tmpopu, "(T *u," => "(T *u, T *u_xi,");   
        tmpopu = replace(tmpopu, "(u," => "(u, u_xi,");   
        tmpopu = replace(tmpopu, "Gradient(double *" => "Gradient(double *, double *");   
        tmpopu = replace(tmpopu, "Gradient(float *" => "Gradient(float *, float *");              
        tmpcpu = replace(tmpcpu, "(T *u," => "(T *u, T *u_xi,");   
        tmpcpu = replace(tmpcpu, "(u," => "(u, u_xi,");   
        tmpcpu = replace(tmpcpu, "Gradient(double *" => "Gradient(double *, double *");   
        tmpcpu = replace(tmpcpu, "Gradient(float *" => "Gradient(float *, float *");              
        tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_xi,");   
        tmpgpu = replace(tmpgpu, "(u," => "(u, u_xi,");   
        tmpgpu = replace(tmpgpu, "Gradient(double *" => "Gradient(double *, double *");   
        tmpgpu = replace(tmpgpu, "Gradient(float *" => "Gradient(float *, float *");      
        # tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xi,");   
        # tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xi,");   
        # tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*");   
        # tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*");      
    elseif pot==2
        tmpopu = replace(tmpopu, "(T *u," => "(T *u, T *u_xij,");
        tmpopu = replace(tmpopu, "(u," => "(u, u_xij,");    
        tmpopu = replace(tmpopu, "Gradient(double *," => "Gradient(double *, double *,");        
        tmpopu = replace(tmpopu, "Gradient(float *," => "Gradient(float *, float *,");        
        tmpcpu = replace(tmpcpu, "(T *u," => "(T *u, T *u_xij,");
        tmpcpu = replace(tmpcpu, "(u," => "(u, u_xij,");    
        tmpcpu = replace(tmpcpu, "Gradient(double *," => "Gradient(double *, double *,");        
        tmpcpu = replace(tmpcpu, "Gradient(float *," => "Gradient(float *, float *,");        
        tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_xij,");
        tmpgpu = replace(tmpgpu, "(u," => "(u, u_xij,");    
        tmpgpu = replace(tmpgpu, "Gradient(double *," => "Gradient(double *, double *,");        
        tmpgpu = replace(tmpgpu, "Gradient(float *," => "Gradient(float *, float *,");                
        # tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xij,");   
        # tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xij,");   
        # tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*");   
        # tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*");   
    elseif pot==3
        tmpopu = replace(tmpopu, "(T *u," => "(T *u, T *u_xij, T *u_xik,");
        tmpopu = replace(tmpopu, "(u," => "(u, u_xij, u_xik,");    
        tmpopu = replace(tmpopu, "Gradient(double *," => "Gradient(double *, double *, double *,");        
        tmpopu = replace(tmpopu, "Gradient(float *," => "Gradient(float *, float *, float *,");                
        tmpcpu = replace(tmpcpu, "(T *u," => "(T *u, T *u_xij, T *u_xik,");
        tmpcpu = replace(tmpcpu, "(u," => "(u, u_xij, u_xik,");    
        tmpcpu = replace(tmpcpu, "Gradient(double *," => "Gradient(double *, double *, double *,");        
        tmpcpu = replace(tmpcpu, "Gradient(float *," => "Gradient(float *, float *, float *,");                
        tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_xij, T *u_xik,");
        tmpgpu = replace(tmpgpu, "(u," => "(u, u_xij, u_xik,");    
        tmpgpu = replace(tmpgpu, "Gradient(double *," => "Gradient(double *, double *, double *,");        
        tmpgpu = replace(tmpgpu, "Gradient(float *," => "Gradient(float *, float *, float *,");        
        # tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xij, T *u_xik,");   
        # tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xij, u_xik,");   
        # tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*, double*");   
        # tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*, float*");      
    elseif pot==4       
        tmpopu = replace(tmpopu, "(T *u," => "(T *u, T *u_xij, T *u_xik, T *u_xil,");
        tmpopu = replace(tmpopu, "(u," => "(u, u_xij, u_xik, u_xil,");    
        tmpopu = replace(tmpopu, "Gradient(double *," => "Gradient(double *, double *, double *, double *,");        
        tmpopu = replace(tmpopu, "Gradient(float *," => "Gradient(float *, float *, float *, float *,");                                
        tmpcpu = replace(tmpcpu, "(T *u," => "(T *u, T *u_xij, T *u_xik, T *u_xil,");
        tmpcpu = replace(tmpcpu, "(u," => "(u, u_xij, u_xik, u_xil,");    
        tmpcpu = replace(tmpcpu, "Gradient(double *," => "Gradient(double *, double *, double *, double *,");        
        tmpcpu = replace(tmpcpu, "Gradient(float *," => "Gradient(float *, float *, float *, float *,");                                
        tmpgpu = replace(tmpgpu, "(T *u," => "(T *u, T *u_xij, T *u_xik, T *u_xil,");
        tmpgpu = replace(tmpgpu, "(u," => "(u, u_xij, u_xik, u_xil,");    
        tmpgpu = replace(tmpgpu, "Gradient(double *," => "Gradient(double *, double *, double *, double *,");        
        tmpgpu = replace(tmpgpu, "Gradient(float *," => "Gradient(float *, float *, float *, float *,");                        
        # tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xij, T *u_xik, T *u_xil,");   
        # tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xij, u_xik, u_xil,");   
        # tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*, double*, double*");   
        # tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*, float*, float*");      
    end
    stropu = stropu * "\n" * tmpopu;  
    strcpu = strcpu * "\n" * tmpcpu;  
    strgpu = strgpu * "\n" * tmpgpu;  

    cppfiles(foldername, filename, stropu, strcpu, strgpu, 0);
end    

end

# function modifystr1(oldstr)
    
# newstr = oldstr;
# s = "u[" * string(0) * " + i*" * string(1) * "]";    
# if contains(newstr, s)   
#     news = "u[i]";
#     newstr = replace(newstr, s => news);    
# end
   
# return newstr 

# end
