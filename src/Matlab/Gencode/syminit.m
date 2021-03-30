function [xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al] = syminit(app)
%, mu1a, mu1b, mu2a, mu2b, mu2c, mu3a, mu3b, mu3c, mu4a, mu4b, mu4c

nd = app.dim;
ncq = app.ncq;
% ncmu1a = app.ncmu1a;
% ncmu1b = app.ncmu1b;
% ncmu2a = app.ncmu2a;
% ncmu2b = app.ncmu2b;
% ncmu2c = app.ncmu2c;
% ncmu3a = app.ncmu3a;
% ncmu3b = app.ncmu3b;
% ncmu3c = app.ncmu3c;
% ncmu4a = app.ncmu4a;
% ncmu4b = app.ncmu4b;
% ncmu4c = app.ncmu4c;

xij = sym('xij',[nd 1]);
xik = sym('xik',[nd 1]);
xil = sym('xil',[nd 1]);

xi = sym('xi',[nd 1]); 
xj = sym('xj',[nd 1]); 
xk = sym('xk',[nd 1]); 
xl = sym('xl',[nd 1]); 

qi = sym('qi',[ncq 1]); 
qj = sym('qj',[ncq 1]); 
qk = sym('qk',[ncq 1]); 
ql = sym('ql',[ncq 1]); 

ti = sym('ti'); 
tj = sym('tj'); 
tk = sym('tk'); 
tl = sym('tl'); 

ai = sym('ai'); 
aj = sym('aj'); 
ak = sym('ak'); 
al = sym('al'); 

% mu1a = sym('mu1a',[ncmu1a 1]);
% mu1b = sym('mu1b',[ncmu1b 1]);
% 
% mu2a = sym('mu2a',[ncmu2a 1]);
% mu2b = sym('mu2b',[ncmu2b 1]);
% mu2c = sym('mu2c',[ncmu2c 1]);
% 
% mu3a = sym('mu3a',[ncmu3a 1]);
% mu3b = sym('mu3b',[ncmu3b 1]);
% mu3c = sym('mu3c',[ncmu3c 1]);
% 
% mu4a = sym('mu4a',[ncmu4a 1]);
% mu4b = sym('mu4b',[ncmu4b 1]);
% mu4c = sym('mu4c',[ncmu4c 1]);

end

