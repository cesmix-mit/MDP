function genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot)

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
            [tmpopu, tmpcpu, tmpgpu] = gensingle(fn, ui, xi, qi, ti, ai, mu, eta, kappa, gen, ifile);
        elseif pot==2
            [tmpopu, tmpcpu, tmpgpu] = genpair(fn, ui, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
        elseif pot==3
            [tmpopu, tmpcpu, tmpgpu] = gentriplet(fn, ui, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
        elseif pot==4    
            [tmpopu, tmpcpu, tmpgpu] = genquadruplet(fn, ui, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile);
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
    nopotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot);
else
    [sp0, sp1, sp2, sp3] = getpotstr(pot);
    sp0 = strrep(sp0, "int ng)", "int ng, int potnum)");
    sp1 = strrep(sp1, "int);", "int, int);");
    sp2 = strrep(sp2, "int);", "int, int);");
    [tmpopu, tmpcpu, tmpgpu] = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu + tmpopu;
    strcpu = strcpu + tmpcpu;
    strgpu = strgpu + tmpgpu;    
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