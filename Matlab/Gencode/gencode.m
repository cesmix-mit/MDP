%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function gencode(app)

disp("Generate C++ code ...");
if ~exist(char("app"), 'dir')
    mkdir(char("app"));
end

[xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al] = syminit(app);
eta = sym('eta',[app.nceta 1]);
kappa = sym('eta',[app.nckappa 1]);

potentialfile = str2func(app.potentialfile);
pot = potentialfile();

potform = 1;
filename = "Singlea";
mu = sym('mu',[app.ncmu1a 1]);
if isfield(pot, char(filename))    
    u = pot.Singlea(xi, qi, ti, mu, eta, kappa);     
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Singleb";
mu = sym('mu',[app.ncmu1b 1]);
if isfield(pot, char(filename))        
    u = pot.Singleb(xi, qi, ti, mu, eta, kappa);    
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);        
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

potform = 2;
filename = "Paira";
mu = sym('mu',[app.ncmu2a 1]);
if isfield(pot, char(filename))    
    u = pot.Paira(xij, qi, qj, ti, tj, mu, eta, kappa);
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Pairb";
mu = sym('mu',[app.ncmu2b 1]);
if isfield(pot, char(filename))    
    u = pot.Pairb(xij, qi, qj, ti, tj, mu, eta, kappa);
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Pairc";
mu = sym('mu',[app.ncmu2c 1]);
rho = sym('rho',[1 1]);
if isfield(pot, char(filename))    
    u = pot.Pairc(xij, qi, qj, ti, tj, rho, mu, eta, kappa);    
    genbo2(filename, u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, app.natomtype);       
    genbo2grad(filename + "Gradient", u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, app.natomtype);           
else
    ifile = 0;
    gen = 0;    
    [stropu, strcpu, strgpu] = genpair(filename, [], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    [tmpopu, tmpcpu, tmpgpu] = gendensity("PaircDensity", [], mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  + "\n" + tmpopu + "\n";
    strcpu = strcpu  + "\n" + tmpcpu + "\n";
    strgpu = strgpu  + "\n" + tmpgpu + "\n";    
    cppfiles(filename, stropu, strcpu, strgpu, 1);
    genbo2grad(filename + "Gradient", [], xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, app.natomtype);       
end

potform = 3;
filename = "Tripleta";
mu = sym('mu',[app.ncmu3a 1]);
if isfield(pot, char(filename))    
    u = pot.Tripleta(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Tripletb";
mu = sym('mu',[app.ncmu3b 1]);
if isfield(pot, char(filename))    
    u = pot.Tripletb(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Tripletc";
mu = sym('mu',[app.ncmu3c 1]);
rho = sym('rho',[1 1]);
if isfield(pot, char(filename))    
    u = pot.Tripletc(xij, xik, qi, qj, qk, ti, tj, tk, rho, mu, eta, kappa);
    genbo3(filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, app.natomtype);
    genbo3grad(filename+"Gradient", u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, app.natomtype);
else
    ifile = 0;
    gen = 0;
    fp = "TripletcPair";
    fn = "TripletcDensity";    
    [stropu, strcpu, strgpu] = gentriplet(filename, [], xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
    [tmqopu, tmqcpu, tmqgpu] = genpair(fp, [], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    [tmpopu, tmpcpu, tmpgpu] = gendensity(fn, [], mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  + "\n" + tmqopu + "\n" + tmpopu + "\n";
    strcpu = strcpu  + "\n" + tmqcpu + "\n" + tmpcpu + "\n";
    strgpu = strgpu  + "\n" + tmqgpu + "\n" + tmpgpu + "\n";    
    cppfiles(filename, stropu, strcpu, strgpu, 1);    
    genbo3grad(filename+"Gradient", [], xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, app.natomtype);
end

potform = 4;
filename = "Quadrupleta";
mu = sym('mu',[app.ncmu4a 1]);
if isfield(pot, char(filename))    
    u = pot.Quadrupleta(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu,  eta, kappa);
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Quadrupletb";
mu = sym('mu',[app.ncmu4b 1]);
if isfield(pot, char(filename))    
    u = pot.Quadrupletb(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa);
    genpotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(filename + "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

movefile('app/*', char(app.sourcepath + "C++/Potentials"));
rmdir('app');

end
