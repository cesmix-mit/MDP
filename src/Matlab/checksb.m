nx = 1000;
x = linspace(0.02,1,nx);
x = x(:);
N = 4;
K = 10;
g = sphericalbessel2(x,N,K);
for n = 1:(N+1)
    figure(n);clf;plot(x,g(:,:,n),'LineWidth',1); axis([0 1 -0.4 2]);
    set(gca,'FontSize',16);
    xlabel('$r/r_{\max}$', 'interpreter', 'latex', 'FontSize', 22);
    ylabel('$y_{\ell}(r)$', 'interpreter', 'latex', 'FontSize', 22);
    legend('$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','interpreter', 'latex', 'FontSize', 18);
end

x = rand;
dx = 1e-6;
[g,dg] = sphericalbesselderiv(x,N,K);
gp = sphericalbessel(x+dx,N,K);
gm = sphericalbessel(x-dx,N,K);
fdg = (gp-gm)/(2*dx);
fprintf("Maximum derivative error: %g\n", max(abs(dg(:)-fdg(:)))./max(abs(dg(:))));



