function Singlea(xi, qi, ti, mu, eta, kappa)
    u = 0.0;    
    return u;
end

# one-body bonded potentials 
function Singleb(xi, qi, ti, mu, eta, kappa)
    u = 0.0;    
    return u;

# two-body non-bonded potentials 
function Paira(xij, qi, qj, ti, tj, mu, eta, kappa)
    #u = 0.0;    
    r2 = xij[1]*xij[1] + xij[2]*xij[2] + xij[3]*xij[3];    
    r4 = r2*r2;
    r6 = r2*r4;
    r12 = r6*r6;
    A = mu[1];
    B = mu[2];
    u = (A/r12 - B/r6); 
    return u;
end

# two-body bonded potentials 
function Pairb(xij, qi, qj, ti, tj, mu, eta, kappa)
    u = 0.0;        
    return u;
end

# two-body bond order potentials (EAM-like potentials)
function Pairc(xij, qi, qj, ti, tj, rho, mu, eta, kappa)
    u = 0.0;    
    return u;
end

# three-body non-bonded potentials 
function Tripleta(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa)
    u = 0.0;        
    return u; 
end

# three-body bonded potentials 
function Tripletb(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa)
    u = 0.0;            
    return u;
end

# three-body bond porder potentials (Terssof-like potentials)
function Tripletc(xij, xik, qi, qj, qk, ti, tj, tk, rho, mu, eta, kappa)
    u = 0.0;
    return u;
end

# four-body non-bonded potentials
function Quadrupleta(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa)
    u = 0.0;    
    return u;
end

# four-body bonded potentials 
function Quadrupletb(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa)
    u = 0.0;    
    return u;
end

