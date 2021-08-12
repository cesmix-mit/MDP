function pot = ZBLpotential
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

function [e,ep,epp] = ezbl(r,zi,zj)
  pzbl = 0.23;
  a0 = 0.46850;
  c1 = 0.02817;
  c2 = 0.28022;
  c3 = 0.50986;
  c4 = 0.18175;
  d1 = 0.20162;
  d2 = 0.40290;
  d3 = 0.94229;
  d4 = 3.19980;
  qqr2e = 14.399645;
  qelectron = 1.0;

  ainv = (power(zi,pzbl) + power(zj,pzbl))/(a0);
  d1a = d1*ainv;
  d2a = d2*ainv;
  d3a = d3*ainv;
  d4a = d4*ainv;
  zze = zi*zj*qqr2e*qelectron*qelectron;
  
  %[d1a d2a d3a d4a] %-> 2.3090    4.6140   10.7912   36.6444
  
  e1 = exp(-d1a*r);
  e2 = exp(-d2a*r);
  e3 = exp(-d3a*r);
  e4 = exp(-d4a*r);
  
  rinv = 1/r;
  sum = c1*e1;
  sum = sum + c2*e2;
  sum = sum + c3*e3;
  sum = sum + c4*e4;

  sum_p = -c1*d1a*e1;
  sum_p = sum_p - c2*d2a*e2;
  sum_p = sum_p - c3*d3a*e3;
  sum_p = sum_p - c4*d4a*e4;
  
  sum_pp = c1*e1*d1a*d1a;
  sum_pp = sum_pp + c2*e2*d2a*d2a;
  sum_pp = sum_pp + c3*e3*d3a*d3a;
  sum_pp = sum_pp + c4*e4*d4a*d4a;
  
  e = zze*sum*rinv;      
  ep = zze*(sum_p - sum*rinv)*rinv;
  epp = zze*(sum_pp - 2.0*sum_p*rinv + 2.0*sum*rinv*rinv)*rinv;
end

% two-body non-bonded potentials 
function u = Paira(xij, qi, qj, ti, tj, mu, eta, kappa)
    
    r2 = xij(1)*xij(1) + xij(2)*xij(2) + xij(3)*xij(3);    
    r = sqrt(r2);    

    cut_inner = 4;
    cut_outer = 4.8;
    zi = 73;
    zj = 73;    
%     cut_inner = mu(1);
%     cut_outer = mu(2);
%     zi = mu(3);
%     zj = mu(3);
    
    % ZBL potential
    u = ezbl(r,zi,zj);
    
    tc = cut_outer - cut_inner;
    [fc, fcp, fcpp] = ezbl(cut_outer, zi, zj);
    swa = (-3.0*fcp + tc*fcpp)/(tc*tc);
    swb = ( 2.0*fcp - tc*fcpp)/(tc*tc*tc);
    swc = -fc + (tc/2.0)*fcp - (tc*tc/12.0)*fcpp;

    %sw1 = swa;
    %sw2 = swb;
    sw3 = swa/3.0;
    sw4 = swb/4.0;
    sw5 = swc;
    %[sw3 sw4 sw5] % ->  [0.0456   -0.0343   -0.0163] 
    
    t = r - cut_inner;
    sinner = sw5;
    souter = sw5 + t*t*t * (sw3 + sw4*t);
    
    % switching function
    a = 0.5 + 0.5*tanh(1e4*t); % a = 0 when t<0, and a = 1 when t>0 
    s = a*souter + (1-a)*sinner;
    u = u + s;
      
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

