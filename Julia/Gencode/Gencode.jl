#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module Gencode

using Revise, SymPy

#export syminit, gencode, compilecode
export syminit, gencode, compile, runcode, tring2cmd

#include("initializeapp.jl");
include("syminit.jl");
include("varsassign.jl");
include("symsassign.jl");
include("genpotnum.jl");
include("getpotstr.jl");
include("cppfiles.jl");
include("contains.jl");
include("gendensity.jl");
include("gensingle.jl");
include("genpair.jl");
include("gentriplet.jl");
include("genquadruplet.jl");
include("genbo2.jl");
include("genbo3.jl");
include("genpotential.jl");
include("nopotential.jl");
include("string2cmd.jl");
include("compile.jl");
include("runcode.jl");

include("derivassign.jl");
include("genpotentialgrad.jl");
include("nopotentialgrad.jl");
include("gensinglegrad.jl");
include("genpairgrad.jl");
include("gentripletgrad.jl");
include("genquadrupletgrad.jl");
include("gendensitygrad.jl");
include("genbo2grad.jl");
include("genbo3grad.jl");

function gencode(app)

print("generate code...\n");
if !isdir("app")
    mkdir("app");
end

potential = getfield(Main, Symbol("Main"))    

xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al = syminit(app);
eta = [SymPy.symbols("eta$i") for i=1:app.nceta];
kappa = [SymPy.symbols("kappa$i") for i=1:app.nckappa];

foldername = app.sourcepath * "C++/Potentials";
#foldername = app.sourcepath * "Julia/Gencode/app";

potform = 1;
filename = "Singlea";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu1a];
if isdefined(potential, Symbol(filename)) 
    u = potential.Singlea(xi, qi, ti, mu, eta, kappa);     
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Singleb";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu1b];
if isdefined(potential, Symbol(filename)) 
    u = potential.Singleb(xi, qi, ti, mu, eta, kappa);     
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

potform = 2;
filename = "Paira";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu2a];
if isdefined(potential, Symbol(filename)) 
    u = potential.Paira(xij, qi, qj, ti, tj, mu, eta, kappa);     
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Pairb";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu2b];
if isdefined(potential, Symbol(filename)) 
    u = potential.Pairb(xij, qi, qj, ti, tj, mu, eta, kappa);         
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Pairc";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu2c];
rho = [SymPy.symbols("rho$i") for i=1:1];
if isdefined(potential, Symbol(filename)) 
    u = potential.Pairc(xij, qi, qj, ti, tj, rho, mu, eta, kappa);     
    genbo2(foldername, filename, u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, app.natomtype);    
    genbo2grad(foldername, filename * "Gradient", u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, app.natomtype);       
    #function genbo2grad(foldername, filename, u, xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, natomtype)
else
    ifile = 0;
    gen = 0;    
    stropu, strcpu, strgpu = genpair(foldername, filename, [], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    tmpopu, tmpcpu, tmpgpu = gendensity(foldername, "PaircDensity", [], mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  * "\n" * tmpopu * "\n";
    strcpu = strcpu  * "\n" * tmpcpu * "\n";
    strgpu = strgpu  * "\n" * tmpgpu * "\n";    
    cppfiles(foldername, filename, stropu, strcpu, strgpu, 1);    
    genbo2grad(foldername, filename * "Gradient", [], xij, qi, qj, ti, tj, ai, aj, rho, mu, eta, kappa, app.natomtype);       
end

potform = 3;
filename = "Tripleta";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu3a];
if isdefined(potential, Symbol(filename)) 
    u = potential.Tripleta(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Tripletb";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu3b];
if isdefined(potential, Symbol(filename)) 
    u = potential.Tripletb(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa);
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Tripletc";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu3c];
rho = [SymPy.symbols("rho$i") for i=1:1];
if isdefined(potential, Symbol(filename)) 
    u = potential.Tripletc(xij, xik, qi, qj, qk, ti, tj, tk, rho, mu, eta, kappa);
    genbo3(foldername, filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, app.natomtype);
    genbo3grad(foldername, filename * "Gradient", u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, app.natomtype);
else
    ifile = 0;
    gen = 0;
    fp = "TripletcPair";
    fn = "TripletcDensity";    
    stropu, strcpu, strgpu = gentriplet(foldername, filename, [], xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
    tmqopu, tmqcpu, tmqgpu = genpair(foldername, fp, [], xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
    tmpopu, tmpcpu, tmpgpu = gendensity(foldername, fn, [], mu, mu, eta, kappa, gen, ifile);    
    stropu = stropu  * "\n" * tmqopu * "\n" * tmpopu * "\n";
    strcpu = strcpu  * "\n" * tmqcpu * "\n" * tmpcpu * "\n";
    strgpu = strgpu  * "\n" * tmqgpu * "\n" * tmpgpu * "\n";    
    cppfiles(foldername, filename, stropu, strcpu, strgpu, 1);       
    genbo3grad(foldername, filename * "Gradient", [], xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, eta, kappa, app.natomtype); 
end

potform = 4;
filename = "Quadrupleta";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu4a];
if isdefined(potential, Symbol(filename)) 
    u = potential.Quadrupleta(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu,  eta, kappa);
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

filename = "Quadrupletb";
mu = [SymPy.symbols("mu$i") for i=1:app.ncmu4b];
if isdefined(potential, Symbol(filename)) 
    u = potential.Quadrupletb(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu,  eta, kappa);
    genpotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
else
    nopotential(foldername, filename, [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
    genpotentialgrad(foldername, filename * "Gradient", [], xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, potform);    
end

end

end
