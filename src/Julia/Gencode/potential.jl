#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function Singlea(xi, qi, ti, mu, eta, kappa)
    u = 0.0;    
    return u;
end

# one-body bonded potentials 
function Singleb(xi, qi, ti, mu, eta, kappa)
    u = 0.0;    
    return u;
end

# two-body non-bonded potentials 
function Paira(xij, qi, qj, ti, tj, mu, eta, kappa)
    #u = 0.0;    
    r2 = xij[1]*xij[1] + xij[2]*xij[2] + xij[3]*xij[3];    
    r = sqrt(r2);
    r3 = r2*r; r4 = r2*r2; r5 = r4*r; r6 = r2*r4;
    r7 = r6*r; r12 = r6*r6;
    
    # shitfted Lennard-Jones potential   
    R = eta[1];     # cut-off radius for shitfted Lennard-Jones potential   
    R2 = R*R;
    R4 = R2*R2;
    R6 = R4*R2;
    R12 = R6*R6;
    A = mu[1];
    B = mu[2];
    uLJ = (A/r12 - B/r6) - (A/R12 - B/R6);    
    
    # shielded Coulomb potential (cf. ReaxFF)
    Rcut = eta[2]; # cut-off radius for Coulomb potential   
    tap0 = 1; 
    tap4 = -35/(Rcut^4);
    tap5 = 84/(Rcut^5);
    tap6 = -70/(Rcut^6);
    tap7 = 20/(Rcut^7);
    tap = tap0 + tap4*r4 + tap5*r5 + tap6*r6 + tap7*r7;
    gamma = mu[3];
    C  = mu[4];
    rs = (r3 + (1/gamma)^3)^(1/3);
    uC = tap*C*qi[1]*qj[1]/rs; 
    
    # shielded Morse potential (cf. ReaxFF)
    Rcut = eta[3]; # cut-off radius for shielded Morse potential   
    tap0 = 1; 
    tap4 = -35/(Rcut^4);
    tap5 = 84/(Rcut^5);
    tap6 = -70/(Rcut^6);
    tap7 = 20/(Rcut^7);
    tap = tap0 + tap4*r4 + tap5*r5 + tap6*r6 + tap7*r7;
    gamma = mu[5];
    D  = mu[6];
    a  = mu[7];
    re  = mu[8];    
    p = mu[9];
    rs = (r^p + (1/gamma)^p)^(1.0/p);
    tm = a*(1 - rs/re);
    uM = tap*D*(exp(tm) - 2*exp(0.5*tm)); 
    
    u = [uLJ uC uM];

    return u;
end

# two-body bonded potentials 
function Pairb(xij, qi, qj, ti, tj, mu, eta, kappa)
    #u = 0.0;    
    r2 = xij[1]*xij[1] + xij[2]*xij[2] + xij[3]*xij[3];    
    r = sqrt(r2);
    
    # bonded Morse potential 
    D  = mu[1];
    a  = mu[2];
    rc = mu[3];            
    tm = 2*a*(rc - r);
    uM = D*(exp(tm) - 2*exp(tm));     
    
    # bonded harmonic potential 
    K  = mu[4];
    r0  = mu[5];
    rc  = mu[6];                
    uH = K*((r-r0)^2 - (rc-r0)^2);
    
    u = [uM uH];

    return u;
end

# two-body bond order potentials (EAM-like potentials)
function Pairc(xij, qi, qj, ti, tj, rho, mu, eta, kappa)
    #u = 0.0;
    r2 = xij[1]*xij[1] + xij[2]*xij[2] + xij[3]*xij[3];    
    r = sqrt(r2);
    
    # EAM Potential 
    a1 = mu[1];
    m1 = mu[2];    
    a2 = mu[3];
    m2 = mu[4];    
    rhoaa = (a1/r)^m1;
    rhoab = (a2/r)^m2;
    Fa = -sqrt(rho[1]);        
    ua = [rhoaa; rhoab; Fa];
    
    u = [ua ua];

    return u;
end

# three-body non-bonded potentials 
function Tripleta(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa)
    u = 0.0;        
    return u; 
end

# three-body bonded potentials 
function Tripletb(xij, xik, qi, qj, qk, ti, tj, tk, mu, eta, kappa)
    #u = 0.0;    
        
    xdot = xij[1]*xik[1] + xij[2]*xik[2] + xij[3]*xik[3];
    rij2 = xij[1]*xij[1] + xij[2]*xij[2] + xij[3]*xij[3];    
    rik2 = xik[1]*xik[1] + xik[2]*xik[2] + xik[3]*xik[3];    
    rij = sqrt(rij2);  # distance between atom i and atom j  
    rik = sqrt(rik2);  # distance between atom i and atom k
    
    # cos(theta), where theta is the angle between xij and xik
    costhe = xdot/(rij*rik);
        
    rjk2 = rij2 + rik2 - 2*rij*rik*costhe;
    rjk = sqrt(rjk2); # distance between atom j and atom k
    
    u1 = mu[1]*(costhe - mu[2])^2;   
    u2 = u1 + mu[3]*(mu[4] - rjk)^2;
    u = [u1 u2];

    return u;
end

# three-body bond porder potentials (Terssof-like potentials)
function Tripletc(xij, xik, qi, qj, qk, ti, tj, tk, rho, mu, eta, kappa)
    #u = 0.0;
    
    xdot = xij[1]*xik[1] + xij[2]*xik[2] + xij[3]*xik[3];
    rij2 = xij[1]*xij[1] + xij[2]*xij[2] + xij[3]*xij[3];    
    rik2 = xik[1]*xik[1] + xik[2]*xik[2] + xik[3]*xik[3];    
    rij = sqrt(rij2);  # distance between atom i and atom j  
    rik = sqrt(rik2);  # distance between atom i and atom k
    
    # cos(theta), where theta is the angle between xij and xik
    costhe = xdot/(rij*rik);
   
    B = mu[1];
    lambda2 = mu[2];
    lambda3 = mu[3];
    gamma = mu[4];
    c = mu[5];
    d = mu[6];
    cos0 = mu[7];
    m = mu[8];
    n = mu[9];
    beta = mu[10];
    
    fA = -B*exp(-lambda2*rij);
    fB = gamma*(1 + c^2/d^2 - c^2/(d^2 + (costhe-cos0)^2))*exp(lambda3^m*(rij - rik)^m);
    fC = gamma*(1 + c^2/d^2 - c^2/(d^2 + (costhe-cos0)^2))*exp(lambda3^m*(rij - rik)^m);
    fD = (1 + beta^n*rho[1]^n)^(-1.0/(2*n));
    u = [fB; fA; fD]; 

    u = reshape(u,(3,1));
    
    return u;
end

# four-body non-bonded potentials
function Quadrupleta(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa)
    u = 0.0;    
    return u;
end

# four-body bonded potentials 
function Quadrupletb(xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, mu, eta, kappa)
    #u = 0.0;    
    
    # calculate the distance d from atom i to a plane defined by 3 atoms (j, k, l)
    xjk = xik - xij;
    xjl = xil - xij;
    n1 = xjk[2]*xjl[3] - xjk[3]*xjl[2];
    n2 = xjk[3]*xjl[1] - xjk[1]*xjl[3];
    n3 = xjk[1]*xjl[2] - xjk[2]*xjl[1];
    nn = sqrt(n1*n1 + n2*n2 + n3*n3);
    d = (n1*xij[1] + n2*xij[2] + n3*xij[3])/nn;    
    
    u = mu[1]*d^2 + mu[2]*d^4;
    return u;
end

