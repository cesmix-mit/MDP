function genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot)

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
            tmpopu, tmpcpu, tmpgpu = gensingle(foldername, fn, ui, xi, qi, ti, ai, mu, eta, kappa, gen, ifile);
        elseif pot==2
            tmpopu, tmpcpu, tmpgpu = genpair(foldername, fn, ui, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
        elseif pot==3
            tmpopu, tmpcpu, tmpgpu = gentriplet(foldername, fn, ui, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
        elseif pot==4    
            tmpopu, tmpcpu, tmpgpu = genquadruplet(foldername, fn, ui, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile);
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
    nopotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot);
else
    sp0, sp1, sp2, sp3 = getpotstr(pot);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu * tmpopu;
    strcpu = strcpu * tmpcpu;
    strgpu = strgpu * tmpgpu;    

    # here
    if pot==1
        tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xi,");   
        tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xi,");   
        tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*");   
        tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*");      
    elseif pot==2
        tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xij,");   
        tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xij,");   
        tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*");   
        tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*");   
    elseif pot==3
        tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xij, T *u_xik,");   
        tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xij, u_xik,");   
        tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*, double*");   
        tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*, float*");      
    elseif pot==4       
        tmpgpu = replace(tmpgpu, "(T *u," => "Gradient(T *u, T *du, T *u_xij, T *u_xik, T *u_xil,");   
        tmpgpu = replace(tmpgpu, "(u," => "Gradient(u, du, u_xij, u_xik, u_xil,");   
        tmpgpu = replace(tmpgpu, "(double *" => "Gradient(double *, double *, double*, double*, double*");   
        tmpgpu = replace(tmpgpu, "(float *" => "Gradient(float *, float *, float*, float*, float*");      
    end
    strgpu = strgpu * "\n" * tmpgpu;  

    cppfiles(foldername, filename, stropu, strcpu, strgpu, 0);
end    

end

function modifystr1(oldstr)
    
newstr = oldstr;
s = "u[" * string(0) * " + i*" * string(1) * "]";    
if contains(newstr, s)   
    news = "u[i]";
    newstr = replace(newstr, s => news);    
end
   
return newstr 

end
    