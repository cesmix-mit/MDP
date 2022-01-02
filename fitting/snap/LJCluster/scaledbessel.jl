function scaledbessel(xij, qi, qj, ti, tj, mu, eta, kappa)

dim, N = size(xij);
rcut = eta[1]
alpha = eta[2]
nbf = kappa[1]

dij = sqrt.(xij[1,:].*xij[1,:] .+ xij[2,:].*xij[2,:] .+ xij[3,:].*xij[3,:]);    
dr = xij./reshape(dij,(1,N));
x =  (1.0 .- exp.(-alpha*dij/rcut))/(1-exp(-alpha));
dx = -((alpha/rcut)*exp.(-(alpha*dij/rcut)))/(exp(-alpha) - 1);

rbf = zeros(N, nbf);
drbf = zeros(dim, N, nbf);
for i = 1:nbf    
    a = i*pi;
    rbf[:,i] = (sqrt(2/rcut)/i)*(sin.(a*x)./dij);
    b = (sqrt(2/rcut)/i)*(sin.(a*x)./(dij.^2) .- (a*cos.(a*x).*dx)./dij);
    drbf[:,:,i] = -dr.*reshape(b,(1,N))
    # for j = 1:3
    #     drbf[j,:,i] = (xij[i,:].*b)./dij;
    # end
end

return rbf, drbf

end


