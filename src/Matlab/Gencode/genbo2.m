function genbo2(filename, u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, natomtype)

[n,m] = size(u);
if (n==0) || (n == 1)
    m = 0;    
elseif (n == natomtype+1)
    if u(end,1) == 0
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
    ui = u(:,i);    
    if ui(end) ~= 0
        potnum = potnum+1;        
        [stropu1, strcpu1, strgpu1] = genpair(filename+num2str(potnum), ui(1:end-1), xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
        stropu1 = modifystr(stropu1, natomtype);
        strcpu1 = modifystr(strcpu1, natomtype);
        strgpu1 = modifystr(strgpu1, natomtype);
        [stropu2, strcpu2, strgpu2] = gendensity(fn+num2str(potnum), ui(end), rho, mu, eta, kappa, gen, ifile);
        stropu = stropu + stropu1 + "\n" + stropu2 + "\n"; 
        strcpu = strcpu + strcpu1 + "\n" + strcpu2 + "\n"; 
        strgpu = strgpu + strgpu1 + "\n" + strgpu2 + "\n"; 
    end    
end

if potnum==0
    ifile = 0;
    gen = 0;
    [stropu, strcpu, strgpu] = genpair(filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    [tmpopu, tmpcpu, tmpgpu] = gendensity(fn, u, mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  + "\n" + tmpopu + "\n";
    strcpu = strcpu  + "\n" + tmpcpu + "\n";
    strgpu = strgpu  + "\n" + tmpgpu + "\n";    
    cppfiles(filename, stropu, strcpu, strgpu, 1);
else
    [sp0, sp1, sp2, sp3] = getpotstr(2);
    sp0 = strrep(sp0, "int ng)", "int ng, int potnum)");
    sp1 = strrep(sp1, "int);", "int, int);");
    sp2 = strrep(sp2, "int);", "int, int);");
    [tmpopu, tmpcpu, tmpgpu] = genpotnum(filename, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu + tmpopu + "\n";
    strcpu = strcpu + tmpcpu + "\n";
    strgpu = strgpu + tmpgpu + "\n";    
    
    [sp0, sp1, sp2, sp3] = getpotstr(0);
    sp0 = strrep(sp0, "int ng)", "int ng, int potnum)");
    sp1 = strrep(sp1, "int);", "int, int);");
    sp2 = strrep(sp2, "int);", "int, int);");
    [tmpopu, tmpcpu, tmpgpu] = genpotnum(fn, potnum, sp0, sp1, sp2, sp3);
    stropu = stropu + tmpopu;
    strcpu = strcpu + tmpcpu;
    strgpu = strgpu + tmpgpu;    
    
    cppfiles(filename, stropu, strcpu, strgpu, 0);
end    

end

function newstr = modifystr(oldstr, natomtype)
    
newstr = oldstr;
for i = 1:natomtype
    s = "u[" + num2str(i-1) + " + i*" + num2str(natomtype) + "]";    
    if contains(newstr, s)   
        news = "if (tj[i] == " + num2str(i) + ") \n\t\t\t" + "u[i]";
        newstr = strrep(newstr, s, news);    
    end
end

end


