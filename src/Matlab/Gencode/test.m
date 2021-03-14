app.nd = 3;
app.ncq = 1;
app.nceta = 10;
app.nckappa = 2;
app.ncmu1a = 1;
app.ncmu1b = 1;
app.ncmu2a = 10;
app.ncmu2b = 10;
app.ncmu2c = 10;
app.ncmu3a = 10;
app.ncmu3b = 10;
app.ncmu3c = 10;
app.ncmu4a = 4;
app.ncmu4b = 4;
app.ncmu4c = 4;
app.natomtype = 2;
app.potentialfile = "potential0";
app.kappa = [app.nd app.natomtype];
gencode(app);




% [xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al] = syminit(app);
% 
% gen = 1; ifile = 1;
% 
% pot = 1;
% mu = sym('mu',[app.ncmu1a 1]);
% u = [sin(xi(1))*sin(xi(2))*sin(xi(3)) sin(xi(2))*mu(1)]; 
% genpotential("Singlea", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, pot);
% 
% pot = 2;
% mu = sym('mu',[app.ncmu2a 1]);
% u = [sin(xij(1))*sin(xij(2))*sin(xij(3)) + cos(qi(1)) + mu(1)*mu(1)  cos(xij(2))*sin(mu(1))]; 
% genpotential("Paira", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, pot);
% 
% pot = 3;
% mu = sym('mu',[app.ncmu3a 1]);
% u = [sin(xij(1))*sin(xij(2))*sin(xij(3)) + cos(qi(1)) + mu(1)*mu(1)  cos(xik(2))*sin(mu(2))]; 
% genpotential("Tripleta", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, pot);
% 
% pot = 4;
% mu = sym('mu',[app.ncmu4a 1]);
% u = [sin(xij(1))*sin(xij(2))*sin(xij(3)) + cos(qi(1)) + mu(1)*mu(1)  cos(xik(2))*sin(mu(2)) tan(xil(2))*cos(mu(3))]; 
% genpotential("Quadrupleta", u, xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, pot);
% 
% natomtype=2;
% mu = sym('mu',[3 1]);
% rho = sym('rho', [1 1]);
% u = [sin(xij(1))*sin(xij(2))*sin(xij(3)) + cos(qi(1)) + mu(1)*mu(1)  cos(xij(2))*sin(mu(2)) tan(rho(1))*cos(mu(3))]; 
% u = [u(:) u(:)];
% %u = [0; 0; 0];
% genbo2("Pairc", u, xij, qi, qj, ti, tj, ai, aj, rho, mu, natomtype);
% 
% u = [sin(xij(1))*sin(xij(2))*sin(xij(3))  cos(qi(1))+mu(1)*mu(1)  cos(xij(2))*sin(mu(2)) tan(rho(1))*cos(mu(3))]; 
% u = [u(:) u(:) u(:)];
% %u = [0; 0; 0];
% genbo3("Tripletc", u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, rho, mu, natomtype);
% 
