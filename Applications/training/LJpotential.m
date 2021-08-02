function pot = LJpotential
pot.Singlea = @Singlea;
pot.Singleb = @Singleb;
pot.Paira = @Paira;
pot.Pairb = @Pairb;
pot.Pairc = @Pairc;
pot.Tripleta = @Tripleta;
pot.Tripletb = @Tripletb;
pot.Tripletc = @Tripletc;
pot.Quadrupleta = @Quadrupleta;
pot.Quadrupletb = @Quadrupletb;
end

% one-body non-bonded potentials 
function u = Singlea(xi, qi, ti, mu, eta, kappa)
    u = 0.0;    
end

% one-body bonded potentials 
function u = Singleb(xi, qi, ti, mu, eta, kappa)
    u = 0.0;    
end

% two-body non-bonded potentials 
function u = Paira(xij, qi, qj, ti, tj, mu, eta, kappa)
    %u = 0.0;    
    % xij = xj - xi, where xi is the position vector of atom i and xj is
    % the position vector of atom j
    r2 = xij(1)*xij(1) + xij(2)*xij(2) + xij(3)*xij(3);    
    %r = sqrt(r2);
    r4 = r2*r2;
    r6 = r2*r4;
    r12 = r6*r6;
    
    A = mu(1);
    B = mu(2);
    u = (A/r12 - B/r6);        
%     u(2) = 1/r4;
%     u(3) = 1/r2;
%     u(4) = 1/sqrt(r2);
end

% two-body bonded potentials 
function u = Pairb(xij, qi, qj, ti, tj, mu, eta, kappa)
    u = 0.0;    
end

% two-body bond order potentials (EAM-like potentials)
function u = Pairc(xij, qi, qj, ti, tj, rho, mu, eta, kappa)
    u = 0.0;
end

% three-body non-bonded potentials 
function u = Tripleta(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa)

    u = 0.0;         
end

% three-body bonded potentials 
function u = Tripletb(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa)
    u = 0.0;            
end

% three-body bond porder potentials (Terssof-like potentials)
function u = Tripletc(xij, xik, qi, qj, qk, ti, tj, tk, rho, mu, eta, kappa)
    u = 0.0;        
end

% four-body non-bonded potentials
function u = Quadrupleta(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa)
    u = 0.0;    
end

% four-body bonded potentials 
function u = Quadrupletb(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa)
    u = 0.0;        
end

