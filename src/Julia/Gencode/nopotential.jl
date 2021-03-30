function nopotential(foldername, filename, u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, pot)

ifile = 0;
gen = 0;
if pot==0
    stropu, strcpu, strgpu = gendensity(foldername,filename, u, mu, mu, eta, kappa, gen, ifile);
elseif pot==1    
    stropu, strcpu, strgpu = gensingle(foldername, filename, u, xi, qi, ti, ai, mu, eta, kappa, gen, ifile);
elseif pot==2
    stropu, strcpu, strgpu = genpair(foldername, filename, u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, gen, ifile);
elseif pot==3
    stropu, strcpu, strgpu = gentriplet(foldername, filename, u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, gen, ifile);
elseif pot==4    
    stropu, strcpu, strgpu = genquadruplet(foldername, filename, u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, gen, ifile);
end    
cppfiles(foldername, filename, stropu, strcpu, strgpu, 1);

end
