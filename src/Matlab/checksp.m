R = linspace(0.05,1,1000);
y = 0;
z = 0;
L = 1;
K = L+1;

x = 0.5;
[~, indl] = shspectrum(x,y,z,L);
[cg,indm,rowm] = cgcoefficients(indl);
[ar, ai] = sssum(x,y,z,K,L);
p = sspower(ar, ai);
b = ssbispectrum(ar, ai, cg, indl, indm, rowm);

p = zeros(length(p),length(R));
b = zeros(length(b),length(R));
for i = 1:length(R)
    x = R(i);
    [ar, ai] = sssum(x,y,z,K,L);
    p(:,i) = sspower(ar, ai);
    b(:,i) = ssbispectrum(ar, ai, cg, indl, indm, rowm);    
end

figure(1); clf;  plot(R,p(1:1:end,:),'LineWidth',1);
axis([0.05 1 -0.2 1]);
set(gca,'FontSize',20);
xlabel('$r/r_{\max}$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$p_{k_1k_2l}(r)$', 'interpreter', 'latex', 'FontSize', 24);

figure(2); clf;  plot(R,b(1:1:end,:),'LineWidth',1);
axis([0.05 1 -1 1]);
set(gca,'FontSize',20);
xlabel('$r/r_{\max}$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$b_{k_1k_2l_1l_2l}(r)$', 'interpreter', 'latex', 'FontSize', 24);

maxb = max(b,[],2);
indp = find(maxb>0.2);
indm = setdiff(1:size(b,1),indp);

figure(3); clf;  plot(R,b(indp,:),'LineWidth',1);
axis([0.05 1 -0.2 1]);
set(gca,'FontSize',20);
xlabel('$r/r_{\max}$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$b_{k_1k_2l_1l_2l}(r)$', 'interpreter', 'latex', 'FontSize', 24);

% figure(4); clf;  plot(R,b(indm,:),'LineWidth',1);
% axis([0.05 1 -1 0.2]);
% set(gca,'FontSize',20);
% xlabel('$r/r_{\max}$', 'interpreter', 'latex', 'FontSize', 24);
% ylabel('$b_{k_1k_2l_1l_2l}(r)$', 'interpreter', 'latex', 'FontSize', 24);

