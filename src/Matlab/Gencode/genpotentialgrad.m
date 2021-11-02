%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function genpotentialgrad(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot)

n = length(u);
gen = 1;
ifile = 0;
stropu = "";
strcpu = "";
strgpu = "";
potnum = 0; 

for i = 1:n
    ui = u(i);
    if ui ~= 0        
        potnum = potnum+1;
        fn = filename + num2str(potnum);        
        if pot==1
            [tmpopu, tmpcpu, tmpgpu] = gensinglegrad(fn, ui, xi, qi, ti, ai, mu, eta, kappa, gen, ifile);
        elseif pot==2
            [tmpopu, tmpcpu, tmpgpu] = genpairgrad(fn, ui, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);                        
        elseif pot==3
            [tmpopu, tmpcpu, tmpgpu] = gentripletgrad(fn, ui, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
        elseif pot==4    
            [tmpopu, tmpcpu, tmpgpu] = genquadrupletgrad(fn, ui, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile);
        end
        tmpopu = modifystr(tmpopu);
        tmpcpu = modifystr(tmpcpu);
        tmpgpu = modifystr(tmpgpu);
        stropu = stropu + tmpopu + "\n";
        strcpu = strcpu + tmpcpu + "\n";
        strgpu = strgpu + tmpgpu + "\n";        
    end    
end    

if potnum==0
    nopotentialgrad(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot);    
else
    [sp0, sp1, sp2, sp3] = getpotstr(pot);
    sp0 = strrep(sp0, "int ng)", "int ng, int potnum)");
    sp1 = strrep(sp1, "int);", "int, int);");
    sp2 = strrep(sp2, "int);", "int, int);");
    [tmpopu, tmpcpu, tmpgpu] = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);    
    if pot==1        
        tmpopu = strrep(tmpopu, "(T *u,", "(T *u, T *u_xi,");   
        tmpopu = strrep(tmpopu, "(u,", "(u, u_xi,");   
        tmpopu = strrep(tmpopu, "Gradient(double *", "Gradient(double *, double *");   
        tmpopu = strrep(tmpopu, "Gradient(float *", "Gradient(float *, float *");              
        tmpcpu = strrep(tmpcpu, "(T *u,", "(T *u, T *u_xi,");   
        tmpcpu = strrep(tmpcpu, "(u,", "(u, u_xi,");   
        tmpcpu = strrep(tmpcpu, "Gradient(double *", "Gradient(double *, double *");   
        tmpcpu = strrep(tmpcpu, "Gradient(float *", "Gradient(float *, float *");              
        tmpgpu = strrep(tmpgpu, "(T *u,", "(T *u, T *u_xi,");   
        tmpgpu = strrep(tmpgpu, "(u,", "(u, u_xi,");   
        tmpgpu = strrep(tmpgpu, "Gradient(double *", "Gradient(double *, double *");   
        tmpgpu = strrep(tmpgpu, "Gradient(float *", "Gradient(float *, float *");      
    elseif pot==2
        tmpopu = strrep(tmpopu, "(T *u,", "(T *u, T *u_xij,");
        tmpopu = strrep(tmpopu, "(u,", "(u, u_xij,");    
        tmpopu = strrep(tmpopu, "Gradient(double *,", "Gradient(double *, double *,");        
        tmpopu = strrep(tmpopu, "Gradient(float *,", "Gradient(float *, float *,");        
        tmpcpu = strrep(tmpcpu, "(T *u,", "(T *u, T *u_xij,");
        tmpcpu = strrep(tmpcpu, "(u,", "(u, u_xij,");    
        tmpcpu = strrep(tmpcpu, "Gradient(double *,", "Gradient(double *, double *,");        
        tmpcpu = strrep(tmpcpu, "Gradient(float *,", "Gradient(float *, float *,");        
        tmpgpu = strrep(tmpgpu, "(T *u,", "(T *u, T *u_xij,");
        tmpgpu = strrep(tmpgpu, "(u,", "(u, u_xij,");    
        tmpgpu = strrep(tmpgpu, "Gradient(double *,", "Gradient(double *, double *,");        
        tmpgpu = strrep(tmpgpu, "Gradient(float *,", "Gradient(float *, float *,");                
    elseif pot==3
        tmpopu = strrep(tmpopu, "(T *u,", "(T *u, T *u_xij, T *u_xik,");
        tmpopu = strrep(tmpopu, "(u,", "(u, u_xij, u_xik,");    
        tmpopu = strrep(tmpopu, "Gradient(double *,", "Gradient(double *, double *, double *,");        
        tmpopu = strrep(tmpopu, "Gradient(float *,", "Gradient(float *, float *, float *,");                
        tmpcpu = strrep(tmpcpu, "(T *u,", "(T *u, T *u_xij, T *u_xik,");
        tmpcpu = strrep(tmpcpu, "(u,", "(u, u_xij, u_xik,");    
        tmpcpu = strrep(tmpcpu, "Gradient(double *,", "Gradient(double *, double *, double *,");        
        tmpcpu = strrep(tmpcpu, "Gradient(float *,", "Gradient(float *, float *, float *,");                
        tmpgpu = strrep(tmpgpu, "(T *u,", "(T *u, T *u_xij, T *u_xik,");
        tmpgpu = strrep(tmpgpu, "(u,", "(u, u_xij, u_xik,");    
        tmpgpu = strrep(tmpgpu, "Gradient(double *,", "Gradient(double *, double *, double *,");        
        tmpgpu = strrep(tmpgpu, "Gradient(float *,", "Gradient(float *, float *, float *,");        
    elseif pot==4       
        tmpopu = strrep(tmpopu, "(T *u,", "(T *u, T *u_xij, T *u_xik, T *u_xil,");
        tmpopu = strrep(tmpopu, "(u,", "(u, u_xij, u_xik, u_xil,");    
        tmpopu = strrep(tmpopu, "Gradient(double *,", "Gradient(double *, double *, double *, double *,");        
        tmpopu = strrep(tmpopu, "Gradient(float *,", "Gradient(float *, float *, float *, float *,");                                
        tmpcpu = strrep(tmpcpu, "(T *u,", "(T *u, T *u_xij, T *u_xik, T *u_xil,");
        tmpcpu = strrep(tmpcpu, "(u,", "(u, u_xij, u_xik, u_xil,");    
        tmpcpu = strrep(tmpcpu, "Gradient(double *,", "Gradient(double *, double *, double *, double *,");        
        tmpcpu = strrep(tmpcpu, "Gradient(float *,", "Gradient(float *, float *, float *, float *,");                                
        tmpgpu = strrep(tmpgpu, "(T *u,", "(T *u, T *u_xij, T *u_xik, T *u_xil,");
        tmpgpu = strrep(tmpgpu, "(u,", "(u, u_xij, u_xik, u_xil,");    
        tmpgpu = strrep(tmpgpu, "Gradient(double *,", "Gradient(double *, double *, double *, double *,");        
        tmpgpu = strrep(tmpgpu, "Gradient(float *,", "Gradient(float *, float *, float *, float *,");                        
    end
    stropu = stropu + "\n" + tmpopu;  
    strcpu = strcpu + "\n" + tmpcpu;  
    strgpu = strgpu + "\n" + tmpgpu;  
    
    cppfiles(filename, stropu, strcpu, strgpu, 0);
end    

end

function newstr = modifystr(oldstr)
    
newstr = oldstr;
s = "u[" + num2str(0) + " + i*" + num2str(1) + "]";    
if contains(newstr, s)   
    news = "u[i]";
    newstr = strrep(newstr, s, news);    
end

end
