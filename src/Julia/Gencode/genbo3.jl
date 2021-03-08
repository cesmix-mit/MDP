function genbo3(filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, natomtype)

n,m = size(u);
if (n==0) || (n == 1)
    m = 0;    
elseif (n == natomtype+2)
    if u[n,1] == 0
        m = 0;
    end    
else
    error("Number of functions for the two-body bond order potential must be equal to 1 plus the number of atom types");
end
    
gen = 1;
ifile = 0;
stropu = "";
strcpu = "";
strgpu = "";
potnum = 0; 
fp = "TripletcPair";
fn = "TripletcDensity";

for i = 1:m
    ui = u[:,i];    
    if ui[n] != 0
        potnum = potnum + 1;        
        stropu0, strcpu0, strgpu0 = gentriplet(filename * string(potnum), ui[1:(n-2)], xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);        
        stropu0 = modifystr3(stropu0, natomtype);
        strcpu0 = modifystr3(strcpu0, natomtype);
        strgpu0 = modifystr3(strgpu0, natomtype);
        stropu1, strcpu1, strgpu1 = genpair(fp * string(potnum), ui[n-1], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);        
        stropu2, strcpu2, strgpu2 = gendensity(fn * string(potnum), ui[end], rho, mu, eta, kappa, gen, ifile);
        stropu = stropu * stropu0 * "\n" * stropu1 * "\n" * stropu2 * "\n"; 
        strcpu = strcpu * strcpu0 * "\n" * strcpu1 * "\n" * strcpu2 * "\n"; 
        strgpu = strgpu * strgpu0 * "\n" * strgpu1 * "\n" * strgpu2 * "\n"; 
    end    
end

if potnum==0
    ifile = 0;
    gen = 0;
    stropu, strcpu, strgpu = gentriplet(filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
    tmqopu, tmqcpu, tmqgpu = genpair(fp, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    tmpopu, tmpcpu, tmpgpu = gendensity(fn, u, mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  * "\n" * tmqopu * "\n" * tmpopu * "\n";
    strcpu = strcpu  * "\n" * tmqcpu * "\n" * tmpcpu * "\n";
    strgpu = strgpu  * "\n" * tmqgpu * "\n" * tmpgpu * "\n";    
    cppfiles(filename, stropu, strcpu, strgpu, 1);
else
    sp0, sp1, sp2, sp3 = getpotstr(3);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu * tmpopu * "\n";
    strcpu = strcpu * tmpcpu * "\n";
    strgpu = strgpu * tmpgpu * "\n";    
    
    sp0, sp1, sp2, sp3 = getpotstr(2);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(fp, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu * tmpopu * "\n";
    strcpu = strcpu * tmpcpu * "\n";
    strgpu = strgpu * tmpgpu * "\n";    
    
    sp0, sp1, sp2, sp3 = getpotstr(0);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(fn, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu * tmpopu;
    strcpu = strcpu * tmpcpu;
    strgpu = strgpu * tmpgpu;    
    
    cppfiles(filename, stropu, strcpu, strgpu, 0);
end    

end

function modifystr3(oldstr, natomtype)
    
newstr = oldstr;
for i = 1:natomtype
    s = "u[" * string(i-1) * " * i*" * string(natomtype) * "]";    
    if contains(newstr, s)   
        news = "if (tk[i] == " * string(i) * ") \n\t\t\t" * "u[i]";
        newstr = replace(newstr, s, news);    
    end
end

return newstr

end


