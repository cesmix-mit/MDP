function genbo2(filename, u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, natomtype)

n,m = size(u);
if (n==0) || (n == 1)
    m = 0;    
elseif (n == natomtype+1)
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
fn = "PaircDensity";

for i = 1:m
    ui = u[:,i];    
    if (ui[n] != 0)
        potnum = potnum + 1;        
        stropu1, strcpu1, strgpu1 = genpair(filename*string(potnum), ui[1:(n-1)], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
        stropu1 = modifystr2(stropu1, natomtype);
        strcpu1 = modifystr2(strcpu1, natomtype);
        strgpu1 = modifystr2(strgpu1, natomtype);
        stropu2, strcpu2, strgpu2 = gendensity(fn*string(potnum), ui[n], rho, mu, eta, kappa, gen, ifile);
        stropu = stropu * stropu1 * "\n" * stropu2 * "\n"; 
        strcpu = strcpu * strcpu1 * "\n" * strcpu2 * "\n"; 
        strgpu = strgpu * strgpu1 * "\n" * strgpu2 * "\n"; 
    end    
end

if potnum==0
    ifile = 0;
    gen = 0;
    stropu, strcpu, strgpu = genpair(filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    tmpopu, tmpcpu, tmpgpu = gendensity(fn, u, mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  * "\n" * tmpopu * "\n";
    strcpu = strcpu  * "\n" * tmpcpu * "\n";
    strgpu = strgpu  * "\n" * tmpgpu * "\n";    
    cppfiles(filename, stropu, strcpu, strgpu, 1);
else
    sp0, sp1, sp2, sp3 = getpotstr(2);
    sp0 = replace(sp0, "int ng)" => "int ng, int potnum)");
    sp1 = replace(sp1, "int);" => "int, int);");
    sp2 = replace(sp2, "int);" => "int, int);");
    tmpopu, tmpcpu, tmpgpu = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
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

function modifystr2(oldstr, natomtype)
    
newstr = oldstr;
for i = 1:natomtype
    s = "u[" * string(i-1) * " * i*" * string(natomtype) * "]";    
    if contains(newstr, s)   
        news = "if (tj[i] == " * string(i) * ") \n\t\t\t" * "u[i]";
        newstr = replace(newstr, s => news);    
    end
end

return newstr 

end


