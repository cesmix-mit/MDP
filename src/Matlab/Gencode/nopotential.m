%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function nopotential(filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot)

ifile = 0;
gen = 0;
if pot==0
    [stropu, strcpu, strgpu] = gendensity(filename, u, mu, mu, eta, kappa, gen, ifile);
elseif pot==1    
    [stropu, strcpu, strgpu] = gensingle(filename, u, xi, qi, ti, ai, mu, eta, kappa, gen, ifile);
elseif pot==2
    [stropu, strcpu, strgpu] = genpair(filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
elseif pot==3
    [stropu, strcpu, strgpu] = gentriplet(filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
elseif pot==4    
    [stropu, strcpu, strgpu] = genquadruplet(filename, u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile);
end    
cppfiles(filename, stropu, strcpu, strgpu, 1);
