function zbl(xij, qi, qj, ti, tj, mu, eta, kappa)

    ng = length(ti)
    u = zeros(ng)
    u_xij = zeros(3,ng)
    for i = 1:ng    
        xij1 = xij[1,i];
        xij2 = xij[2,i];
        xij3 = xij[3,i];    
        t2 = xij1*xij1;
        t3 = xij2*xij2;
        t4 = xij3*xij3;
        t5 = t2+t3+t4;
        t6 = sqrt(t5);
        t7 = t6*1.0E2;
        t8 = t7-4.0E2;
        t9 = tanh(t8);
        t10 = t6-4.0;
        t11 = t10*t10;
        t12 = 1.0/sqrt(t5);
        t13 = t6*3.428291787990583E-2;
        t14 = t13-1.827782636402304E-1;
        t15 = t9*(1.0/2.0);
        t16 = t15+1.0/2.0;
        t21 = t6*2.308969885292557;
        t17 = exp(-t21);
        t23 = t6*4.614046060829141;
        t18 = exp(-t23);
        t25 = t6*3.664438963872198E1;
        t19 = exp(-t25);
        t27 = t6*1.079118754693148E1;
        t20 = exp(-t27);
        t22 = t17*2.16164490013485E3;
        t24 = t18*2.15028801532051E4;
        t26 = t19*1.394671496625875E4;
        t28 = t20*3.91244681854013E4;
        t29 = t22+t24+t26+t28;
        t30 = t9*t9;
        t31 = t30-1.0;
        t32 = t10*t11*t14;
        t33 = t32+1.625187736600612E-2;
        t34 = 1.0/t5^(3.0/2.0);
        u[i] = t9*8.125938683003058E-3+t12*t29-t16*t33-8.125938683003058E-3;
        u_xij[1,i] = -t12*(t12*t17*xij1*4.991172977107606E3+t12*t18*xij1*9.921527946737711E4+t12*t19*xij1*5.11068857403781E5+t12*t20*xij1*4.221994738626192E5)-t16*(t10*t11*t12*xij1*3.428291787990583E-2+t11*t12*t14*xij1*3.0)-t12*t31*xij1*8.125938683003058E-1-t29*t34*xij1+t12*t31*t33*xij1*5.0E1;
        u_xij[2,i] = -t12*(t12*t17*xij2*4.991172977107606E3+t12*t18*xij2*9.921527946737711E4+t12*t19*xij2*5.11068857403781E5+t12*t20*xij2*4.221994738626192E5)-t16*(t10*t11*t12*xij2*3.428291787990583E-2+t11*t12*t14*xij2*3.0)-t12*t31*xij2*8.125938683003058E-1-t29*t34*xij2+t12*t31*t33*xij2*5.0E1;
        u_xij[3,i] = -t12*(t12*t17*xij3*4.991172977107606E3+t12*t18*xij3*9.921527946737711E4+t12*t19*xij3*5.11068857403781E5+t12*t20*xij3*4.221994738626192E5)-t16*(t10*t11*t12*xij3*3.428291787990583E-2+t11*t12*t14*xij3*3.0)-t12*t31*xij3*8.125938683003058E-1-t29*t34*xij3+t12*t31*t33*xij3*5.0E1;    
    end
    u = 0.5*u
    u_xij = 0.5*u_xij

    return u, u_xij
end


# function ezbl(r,zi,zj)
#     pzbl = 0.23;
#     a0 = 0.46850;
#     c1 = 0.02817;
#     c2 = 0.28022;
#     c3 = 0.50986;
#     c4 = 0.18175;
#     d1 = 0.20162;
#     d2 = 0.40290;
#     d3 = 0.94229;
#     d4 = 3.19980;
#     qqr2e = 14.399645;
#     qelectron = 1.0;
  
#     ainv = (sympy.Pow(zi,pzbl) + sympy.Pow(zj,pzbl))/(a0);
#     d1a = d1*ainv;
#     d2a = d2*ainv;
#     d3a = d3*ainv;
#     d4a = d4*ainv;
#     zze = zi*zj*qqr2e*qelectron*qelectron;
    
#     #[d1a d2a d3a d4a] #-> 2.3090    4.6140   10.7912   36.6444
    
#     e1 = exp(-d1a*r);
#     e2 = exp(-d2a*r);
#     e3 = exp(-d3a*r);
#     e4 = exp(-d4a*r);
    
#     rinv = 1/r;
#     sum = c1*e1;
#     sum = sum + c2*e2;
#     sum = sum + c3*e3;
#     sum = sum + c4*e4;
  
#     sum_p = -c1*d1a*e1;
#     sum_p = sum_p - c2*d2a*e2;
#     sum_p = sum_p - c3*d3a*e3;
#     sum_p = sum_p - c4*d4a*e4;
    
#     sum_pp = c1*e1*d1a*d1a;
#     sum_pp = sum_pp + c2*e2*d2a*d2a;
#     sum_pp = sum_pp + c3*e3*d3a*d3a;
#     sum_pp = sum_pp + c4*e4*d4a*d4a;
    
#     e = zze*sum*rinv;      
#     ep = zze*(sum_p - sum*rinv)*rinv;
#     epp = zze*(sum_pp - 2.0*sum_p*rinv + 2.0*sum*rinv*rinv)*rinv;
  
#     return e, ep, epp
#   end
  
#   function zbl(xij, qi, qj, ti, tj, mu, eta, kappa)
      
#       r2 = xij[1,:].*xij[1,:] + xij[2,:]*xij[2,:] + xij[3,:]*xij[3,:];    
#       r = sqrt.(r2);    
  
#       cut_inner = 4;
#       cut_outer = 4.8;
#       zi = 73;
#       zj = 73;    
      
#       # ZBL potential
#       u = ezbl(r,zi,zj)[1];
      
#       tc = cut_outer - cut_inner;
#       fc, fcp, fcpp = ezbl(cut_outer, zi, zj);
#       swa = (-3.0*fcp + tc*fcpp)/(tc*tc);
#       swb = ( 2.0*fcp - tc*fcpp)/(tc*tc*tc);
#       swc = -fc + (tc/2.0)*fcp - (tc*tc/12.0)*fcpp;
  
#       #sw1 = swa;
#       #sw2 = swb;
#       sw3 = swa/3.0;
#       sw4 = swb/4.0;
#       sw5 = swc;
#       #[sw3 sw4 sw5] # ->  [0.0456   -0.0343   -0.0163] 
      
#       t = r - cut_inner;
#       sinner = sw5;
#       souter = sw5 + t*t*t * (sw3 + sw4*t);
      
#       # switching function
#       a = 0.5 + 0.5*tanh(1e2*t); # a = 0 when t<0, and a = 1 when t>0 
#       s = a*souter + (1-a)*sinner;
#       u = u + s;      
#       return u;
#   end
  
  