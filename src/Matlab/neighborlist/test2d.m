a = [1 0.0];
b = [0.5 1]*2;
r = [0.1 0.1]*2.0;
pbc = [1 1];
n = 100;
dim = length(a);
xx = rand(2,n);

[P2S, S2P] = cubemapping(a, b);
x = S2P*xx;
[v, w] = boundingbox(r, [1 1], a, b);
vs = P2S*v; ws = P2S*w;
ximages = boxperiodicimages(pbc, a, b);
%[eta1, eta2, eta3] = makecell(ws(:,1), ws(:,3), [r(1)/norm(a); r(2)/norm(b)]);
[eta1, eta2, eta3] = makecell(ws(:,1), ws(:,3), abs(ws(:,1)));
cellcolor = cellcoloring2([length(eta1)-1 length(eta2)-1],[3 3]);

[xi, ilist, inum, gnum] = createilist(x, ximages, ws, P2S);
[clist, atomcolor] = atomcoloring(xi, P2S, cellcolor, eta1, eta2, eta3, inum, gnum);
n = inum;

[cellist, c2ilist, c2inum, cellnum] = cellist2(xi, ws, P2S, inum, gnum);
[cellcounts, cell2atom] = createcell2atom(cellist, cellnum);
counts = [0 cumsum(c2inum)];
tm = counts-cellcounts;
max(abs(tm(:)))

figure(1); clf;
hold on; axis equal;
plotboundingbox(v);
plotboundingbox(w);
nc = size(cellcolor);
for i = 1:nc(1)    
    for j = 1:nc(2)
        xy = S2P*[eta1(i) eta1(i+1); eta2(j) eta2(j+1)];
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
        text(mean(xy(1,:)),mean(xy(2,:)),num2str(cellcolor(i,j)),txtpars{:});
    end
end
for i = 1:size(atomcolor,1)
    for j = 1:size(atomcolor,2)
        if atomcolor(i,j)==1
            txtpars={'fontname','times','fontsize',18,'horizontala','center','BackgroundColor',[1,1,1]};
            cell = clist(:,j);
            text(xi(1,j),xi(2,j),num2str(cellcolor(cell(1),cell(2))),txtpars{:});
        end
    end
end
plot(xi(1,1:n),xi(2,1:n),'ob','LineWidth',1);
plot(xi(1,n+1:end),xi(2,n+1:end),'or','LineWidth',1);
for i = 1:length(eta1)
    xy = S2P*[eta1(i) eta1(i); eta2(1) eta2(end)];
    plot(xy(1,:),xy(2,:),'-k');
end
for i = 1:length(eta2)        
    xy = S2P*[eta1(1) eta1(end); eta2(i) eta2(i)];
    plot(xy(1,:),xy(2,:),'-k');
end

return;


x = S2P*xx;
[xi, ilist, inum, gnum] = createilist(x, ximages, ws, P2S);
jlist = createjlist(xi, ilist, r, inum);
[clist, c2ilist, c2inum] = cellist(xi, ws, P2S, inum, gnum);
%jlist2 = createjlist2d(xi, ilist, clist, c2ilist, c2inum, r, inum);
atomtype = ones(1,n);
A = [1/r^2 0; 0 1/r^2];
R = chol(A); Rinv = inv(R);
jlist2 = createvlist(xi, atomtype, ilist, clist, c2ilist, c2inum, A, inum, dim);
for i = 1:inum  
    if max(abs(jlist2{i}-jlist{i})) > 0
        error("jlist is wrong");
    end
end
B = A*4;
Q = chol(B); Qinv = inv(Q);
neighlist = neighborlist(xi, atomtype, jlist, B, inum);

nc = size(c2inum);
ncell = prod(nc);
natom = inum+gnum;
cl = clist(1,:) + nc(1)*(clist(2,:)-1);
c2in = zeros(1,ncell);
for k = 1:ncell
    i2 = 0;
    for i = 1:natom
        if cl(i) == k
            i2 = i2+1;
        end
    end
    c2in(k) = i2;
end
if max(abs(c2in(:)-c2inum(:))) > 0
    error("c2in is wrong");
end

c2incum = cumsum([0 c2in]);

c2il = zeros(1,ncell);
for k = 1:ncell
    i1 = c2incum(k);
    i2 = 0;
    for i = 1:natom
        if cl(i) == k
            i2 = i2+1;
            c2il(i1+i2) = i;            
        end
    end
end
c2ili = c2ilist(c2ilist>0);
% if max(abs(c2il(:)-c2ili(:))) > 0
%     error("c2il is wrong");
% end

% initialize the list containting the number of atoms in each cell
c2in = 0*c2in;
for i=1:natom
    k = cl(i);    
    c2in(k) = c2in(k) + 1; 
    i1 = c2incum(k);
    i2 = c2in(k);
    c2il(i1+i2) = i;
end
% if max(abs(c2il(:)-c2ili(:))) > 0
%     error("c2il is wrong");
% end


%     
%     for (int k=0; k<ncell; k++)   // cell k
%     {
%         int i1 = c2inumsum[k];    // 
%         int i2 = 0;
%         for (int i=0; i<natom; i++)   // atom i
%             if (clist[i]==k) {    // if cell k contains atom i
%                 c2ilist[i1+i2] = i;
%                 i2 += 1;
%             }
%     }    


figure(1); clf;
hold on;
plotline(v(:,1), v(:,2));
plotline(v(:,2), v(:,3));
plotline(v(:,3), v(:,4));
plotline(v(:,4), v(:,1));
plot(x(1,:),x(2,:),'ob','LineWidth',1);
v1 = v; v1(1,:) = v1(1,:)+1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)+1,x(2,:),'or');
v1 = v; v1(1,:) = v1(1,:)-1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)-1,x(2,:),'or');
v1 = v; v1(1,:) = v1(1,:)+1; v1(2,:) = v1(2,:)+1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)+1,x(2,:)+1,'or');
v1 = v; v1(1,:) = v1(1,:)+0; v1(2,:) = v1(2,:)+1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)+0,x(2,:)+1,'or');
v1 = v; v1(1,:) = v1(1,:)-1; v1(2,:) = v1(2,:)+1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)-1,x(2,:)+1,'or');
v1 = v; v1(1,:) = v1(1,:)+1; v1(2,:) = v1(2,:)-1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)+1,x(2,:)-1,'or');
v1 = v; v1(1,:) = v1(1,:)+0; v1(2,:) = v1(2,:)-1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)+0,x(2,:)-1,'or');
v1 = v; v1(1,:) = v1(1,:)-1; v1(2,:) = v1(2,:)-1; 
plotline(v1(:,1), v1(:,2));
plotline(v1(:,2), v1(:,3));
plotline(v1(:,3), v1(:,4));
plotline(v1(:,4), v1(:,1));
plot(x(1,:)-1,x(2,:)-1,'or');
% m = [50 100];
% for i = 1:2
%     plot(Rinv(1,1)*cos(t)+Rinv(1,2)*sin(t)+xi(1,m(i)),Rinv(2,1)*cos(t)+Rinv(2,2)*sin(t)+xi(2,m(i)),'-m');
%     plot(Qinv(1,1)*cos(t)+Qinv(1,2)*sin(t)+xi(1,m(i)),Qinv(2,1)*cos(t)+Qinv(2,2)*sin(t)+xi(2,m(i)),'-k');
% end
axis equal;
axis off;

% plotline(w(:,1), w(:,2));
% plotline(w(:,2), w(:,3));
% plotline(w(:,3), w(:,4));
% plotline(w(:,4), w(:,1));

figure(2); clf;
hold on;
plot(xi(1,1:n),xi(2,1:n),'ob','LineWidth',1);
plot(xi(1,n+1:end),xi(2,n+1:end),'or','LineWidth',1);
plotline(v(:,1), v(:,2));
plotline(v(:,2), v(:,3));
plotline(v(:,3), v(:,4));
plotline(v(:,4), v(:,1));
plotline(w(:,1), w(:,2));
plotline(w(:,2), w(:,3));
plotline(w(:,3), w(:,4));
plotline(w(:,4), w(:,1));
t = linspace(0,2*pi,200);
m = [7 50 100 138];
% for i = 1:length(m)
%     plot(Rinv(1,1)*cos(t)+Rinv(1,2)*sin(t)+xi(1,m(i)),Rinv(2,1)*cos(t)+Rinv(2,2)*sin(t)+xi(2,m(i)),'-m');
% end
% for i = 1:length(m)    
%     plot(Qinv(1,1)*cos(t)+Qinv(1,2)*sin(t)+xi(1,m(i)),Qinv(2,1)*cos(t)+Qinv(2,2)*sin(t)+xi(2,m(i)),'-k');
% end
ne1 = floor(1./abs(ws(1,1)));
ne2 = floor(1./abs(ws(2,1)));
eta1 = [ws(1,1) linspace(0,1,ne1+1) ws(1,2)];
eta2 = [ws(2,1) linspace(0,1,ne2+1) ws(2,3)];    
for i = 1:length(eta1)
    xy = S2P*[eta1(i) eta1(i); eta2(1) eta2(end)];
    plot(xy(1,:),xy(2,:),'-k');
end
for i = 1:length(eta2)        
    xy = S2P*[eta1(1) eta1(end); eta2(i) eta2(i)];
    plot(xy(1,:),xy(2,:),'-k');
end
axis equal;

xs = P2S*xi;
figure(3); clf;
hold on;
plot(xs(1,1:n),xs(2,1:n),'ob');
plot(xs(1,n+1:end),xs(2,n+1:end),'or');
plotline(vs(:,1), vs(:,2));
plotline(vs(:,2), vs(:,3));
plotline(vs(:,3), vs(:,4));
plotline(vs(:,4), vs(:,1));
plotline(ws(:,1), ws(:,2));
plotline(ws(:,2), ws(:,3));
plotline(ws(:,3), ws(:,4));
plotline(ws(:,4), ws(:,1));
% for i = 1:length(m)
%     tmp = [Rinv(1,1)*cos(t)+Rinv(1,2)*sin(t)+xi(1,m(i));
%            Rinv(2,1)*cos(t)+Rinv(2,2)*sin(t)+xi(2,m(i))];
%     xy = P2S*tmp;
%     plot(xy(1,:),xy(2,:),'-m');
% end
ne1 = floor(1./abs(ws(1,1)));
ne2 = floor(1./abs(ws(2,1)));
eta1 = [ws(1,1) linspace(0,1,ne1+1) ws(1,2)];
eta2 = [ws(2,1) linspace(0,1,ne2+1) ws(2,3)];    
for i = 1:length(eta1)
    plot([eta1(i) eta1(i)],[eta2(1) eta2(end)],'-g');
end
for i = 1:length(eta2)
    plot([eta1(1) eta1(end)],[eta2(i) eta2(i)],'-g');
end

axis equal;

% a = [1 0 0];
% b = [0.2 1 0];
% c = [0.2 0.2 1];
% r = 0.2;
% [w1] = boundingparallelepiped(a, b, c, r);

% c = a+b;
% a1 = a(1); a2 = a(2);
% b1 = b(1); b2 = b(2);
% 
% [P2S, S2P] = parallelogram2unitsquare(a, b);
% 
% r = 0.1;
% t = linspace(0,2*pi,200);
% 
% [ps,~] = squaremesh(12,12,0,1);
% pp = S2P*ps;
% %figure(1); clf;
% hold on;
% %plot(pp(1,:), pp(2,:), 'o');
% plot([0 a(1)], [0 a(2)], 'k-', 'LineWidth', 2); 
% plot([0 b(1)], [0 b(2)], 'k-', 'LineWidth', 2);
% plot([a(1) c(1)], [a(2) c(2)], 'k-', 'LineWidth', 2);
% plot([b(1) c(1)], [b(2) c(2)], 'k-', 'LineWidth', 2);
% plot(r*cos(t),r*sin(t),'-');
% axis equal;
% 
% % a2*x - a1*y = 0 -> y = a2/a1 * x (1st principal line)
% % a2*x - a1*y = c a parallel line to 1st principal line
% % The distance between the two lines is c/norm(a) 
% c = r*norm(a);
% 
% % b2*x - b1*y = 0 -> x = b1/b2 * y (2nd principal line)
% % b2*x - b1*y = -d a parallel line to 2nd principal line
% % The distance between the two lines is d/norm(b)
% d = r*norm(b);
% 
% % intersection of y = (a2/a2)*x - c/a1 and x = (b1/b2)*y - d/b2
% % x = (b1/b2)*((a2/a1)*x - c/a1) - d/b2
% % (1 - (b1/b2)*((a2/a1)*x = (b1/b2)*(-c/a1) - d/b2
% % s = (a(2)/a(1))*(b(1)/b(2));
% % x1 = -((b(1)/b(2))*(c/a(1)) + d/b(2))/(1-s);
% % y1 = (a(2)/a(1))*x1 - c/a(1);
% 
% % intersection of y = (a2/a2)*x - c/a1 + (i-1)*c/a1 and x = (b1/b2)*y - d/b2
% % intersection of y = (a2/a2)*x - (c-i*c)/a1 and x = (b1/b2)*y - d/b2
% m = 11;
% x1 = zeros(1,m); y1 = x1;
% s = (a(2)/a(1))*(b(1)/b(2));
% for i = 1:m
%     ci = c - (i-1)*c;
%     x1(i) = -((b(1)/b(2))*(ci/a(1)) + d/b(2))/(1-s);
%     y1(i) = (a(2)/a(1))*x1(i) - ci/a(1);
% end
% plot(x1,y1,'-o');
% 
% % intersection of y = (a2/a2)*x - c/a1 and x = (b1/b2)*y - d/b2 + (i-1)*d/b2
% % intersection of y = (a2/a2)*x - (c-i*c)/a1 and x = (b1/b2)*y - d/b2
% x2 = zeros(1,m); y2 = x2;
% for i = 1:m
%     di = d - (i-1)*d;
%     x2(i) = -((b(1)/b(2))*(c/a(1)) + di/b(2))/(1-s);
%     y2(i) = (a(2)/a(1))*x2(i) - c/a(1);
% end
% plot(x2,y2,'-o');
% 
% dx1 = x1(2) - x1(1);
% dy1 = y1(2) - y1(1);
% dx2 = x2(2) - x2(1);
% dy2 = y2(2) - y2(1);
% c = a+b; % corner opposite of the origin 
% d = c + [dx1 dy1] + [dx2 dy2];
% plot(d(1),d(2),'-o');
% 
% ds = sqrt(dx1.^2 + dy1.^2);
% na = floor(norm(a)/ds(1));
% nb = floor(norm(b)/ds(1));
% rs = max(norm(a)/na, norm(b)/nb);
% 
% % for i = 1:(m+1)
% %     x = linspace(x1(i),x1(i)+1.3,10);    
% %     y = (a(2)/a(1))*x - c + (i-1)*c;
% %     plot(x,y,'-');
% % end
% 
% % y = linspace(y1(1),1,10);
% % x = (b(1)/b(2))*y - d;
% % plot(x,y,'-');
% % for i = 1:10
% %     x = (b(1)/b(2))*y - d + i*d;
% %     plot(x,y,'-');
% % end
% 
% 
% %plot(0.376422+r*cos(t),0.2386395+r*sin(t),'-');
% 
% % p1 = p;
% % p1 = -0.1 + 1.2*p;
% % figure(1); clf;
% % plot(p1(1,:), p1(2,:), 'o');
% % 
% % 
% % p2 = S2P*p1;
% % figure(2); clf;
% % plot(p1(1,:), p1(2,:), 'o');
% % hold on;
% % plot(p2(1,:), p2(2,:), '*');
% % plot(0.1*cos(linspace(0,2*pi,200)),0.1*sin(linspace(0,2*pi,200)),'-');
% % plot(p2(1,1:12), p2(2,1:12), '-');
% % axis equal;
% % 
