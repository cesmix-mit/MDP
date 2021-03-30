function syminit(app)

nd = app.dim;
ncq = app.ncq;
# ncmu1a = app.ncmu1a;
# ncmu1b = app.ncmu1b;
# ncmu2a = app.ncmu2a;
# ncmu2b = app.ncmu2b;
# ncmu2c = app.ncmu2c;
# ncmu3a = app.ncmu3a;
# ncmu3b = app.ncmu3b;
# ncmu3c = app.ncmu3c;
# ncmu4a = app.ncmu4a;
# ncmu4b = app.ncmu4b;
# ncmu4c = app.ncmu4c;

xij = [SymPy.symbols("xij$i") for i=1:nd];
xik = [SymPy.symbols("xik$i") for i=1:nd];
xil = [SymPy.symbols("xil$i") for i=1:nd];

xi = [SymPy.symbols("xi$i") for i=1:nd];
xj = [SymPy.symbols("xj$i") for i=1:nd];
xk = [SymPy.symbols("xk$i") for i=1:nd];
xl = [SymPy.symbols("xl$i") for i=1:nd];

qi = [SymPy.symbols("qi$i") for i=1:ncq];
qj = [SymPy.symbols("qj$i") for i=1:ncq];
qk = [SymPy.symbols("qk$i") for i=1:ncq];
ql = [SymPy.symbols("ql$i") for i=1:ncq];

ti = SymPy.symbols("ti");
tj = SymPy.symbols("tj");
tk = SymPy.symbols("tk");
tl = SymPy.symbols("tl");

ai = SymPy.symbols("ai");
aj = SymPy.symbols("aj");
ak = SymPy.symbols("ak");
al = SymPy.symbols("al");

# mu1a = [SymPy.symbols("mu1a$i") for i=1:ncmu1a];
# mu1b = [SymPy.symbols("mu1b$i") for i=1:ncmu1a];

# mu2a = [SymPy.symbols("mu2a$i") for i=1:ncmu2a];
# mu2b = [SymPy.symbols("mu2b$i") for i=1:ncmu2b];
# mu2c = [SymPy.symbols("mu2c$i") for i=1:ncmu2c];

# mu3a = [SymPy.symbols("mu3a$i") for i=1:ncmu3a];
# mu3b = [SymPy.symbols("mu3b$i") for i=1:ncmu3b];
# mu3c = [SymPy.symbols("mu3c$i") for i=1:ncmu3c];

# mu4a = [SymPy.symbols("mu4a$i") for i=1:ncmu4a];
# mu4b = [SymPy.symbols("mu4b$i") for i=1:ncmu4b];
# mu4c = [SymPy.symbols("mu4c$i") for i=1:ncmu4c];

return xij, xik, xil, xi, xj, xk, xl, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al
#, mu1a, mu1b, mu2a, mu2b, mu2c, mu3a, mu3b, mu3c, mu4a, mu4b, mu4c

end

