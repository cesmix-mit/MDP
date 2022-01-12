function scaledbessel(xij, qi, qj, ti, tj, mu, eta, kappa)

dim, N = size(xij);
rcut = eta[1]
alpha = eta[2]
nbf = length(mu)

dij = sqrt.(xij[1,:].*xij[1,:] .+ xij[2,:].*xij[2,:] .+ xij[3,:].*xij[3,:]);    
dr = xij./reshape(dij,(1,N));
x =  (1.0 .- exp.(-alpha*dij/rcut))/(1-exp(-alpha));
dx = -((alpha/rcut)*exp.(-(alpha*dij/rcut)))/(exp(-alpha) - 1);

eij = zeros(N);
fij = zeros(dim, N);
for i = 1:nbf    
    a = i*pi;
    eij = eij + mu[i]*(sqrt(2/rcut)/i)*(sin.(a*x)./dij);
    b = (sqrt(2/rcut)/i)*(sin.(a*x)./(dij.^2) .- (a*cos.(a*x).*dx)./dij);
    fij = fij + -mu[i]*dr.*reshape(b,(1,N))
end

return eij, fij 

end



