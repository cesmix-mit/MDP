% a = [1 0 0];
% b = [0 1 0];
% c = [0 0 1];
% r = 0.1;
a = config.a(:,1);
b = config.b(:,1);
c = config.c(:,1);
r = 8.5;
pbc = [1 1 1];
n = 256;
dim = length(a);
%xx = rand(3,n);
x = config.x(:,n+1:2*n);

[B2C, C2B] = cubemapping(a, b, c);
[v, w] = boundingbox(r, pbc, a, b, c);
ximages = boxperiodicimages(pbc, a, b, c);
vc = B2C*v;
wc = B2C*w;
x = checkconfig(x, ximages, B2C, C2B);

% xx = B2C*x;
% pp = B2C*ximages;
% for i = 1:n
%     xt = xx(:,i);
%     if (0<=xt(1)) && (xt(1)<=1) && (0<=xt(2)) && (xt(2)<=1) && (0<=xt(3)) && (xt(3)<=1)
%     else        
%         xp = xt + pp;        
%         for j = 1:27
%             xq = xp(:,j);
%             if (0<=xq(1)) && (xq(1)<=1) && (0<=xq(2)) && (xq(2)<=1) && (0<=xq(3)) && (xq(3)<=1)
%                 xx(:,i) = xq;                
%                 break;
%             end
%         end
%     end
% end
% x = C2B*xx;

[xi, ilist, inum, gnum] = createilist(x, ximages, wc, B2C);
jlist = createjlist(xi, ilist, r, inum);
[clist, c2ilist, c2inum, cellnum] = cellist(xi, wc, B2C, inum, gnum);
jlist2 = createjlist3d(xi, ilist, clist, c2ilist, c2inum, r, inum);
for i = 1:inum  
    if max(abs(sort(jlist2{i})-sort(jlist{i}))) > 0
        a = sort(jlist2{i})
        b = sort(jlist{i})        
        error("jlist is wrong");
    end
end

epsilon = 0.01029849;
sigma = 3.4;
A = 4*epsilon*sigma^12;
B = 4*epsilon*sigma^6;
u = zeros(1,inum);
ux = zeros(3,inum);
rc2 = r*r;
rc6 = rc2*rc2*rc2;
rc8 = rc6*rc2;
rc12 = rc6*rc6;
rc14 = rc12*rc2;
for i = 1:inum
    ja = jlist{i};    
    xij = xi(:,ja) - xi(:,i);      
    r2 = sum(xij.^2,1);
    r6 = r2.*r2.*r2; 
    r8 = r6.*r2;
    r12 = r6.*r6;
    r14 = r12.*r2;    
    u(i) = 0.5*sum(A./r12 - B./r6);
    xij1 = xij(1,:);
    xij2 = xij(2,:);
    xij3 = xij(3,:);
    ux(1,i) = sum((6*B*xij1)./r8 - (12*A*xij1)./r14);
    ux(2,i) = sum((6*B*xij2)./r8 - (12*A*xij2)./r14);
    ux(3,i) = sum((6*B*xij3)./r8 - (12*A*xij3)./r14);
    % (6*B*x1)/(x1^2 + x2^2 + x3^2)^4 - (12*A*x1)/(x1^2 + x2^2 + x3^2)^7
%     for j = 1:length(ja)        
%         ux(:,i) = ux(:,i) + (6*B*xij(:,j)/r8(j) - 12*A*xij(:,j)/r14(j));       
%     end
end
utot = sum(u);


% a = [1 0.2];
% b = [0.2 1];
% r = 0.2;
% pbc = [1 1];
% n = 100;
% dim = length(a);
% xx = rand(2,n);
% 
% [B2C, S2P] = parallelogram2unitsquare(a, b);
% [v, w] = boundingparallelogram(a, b, r);
% ximages = CreatePeriodicImages(a, b, pbc);
% vs = B2C*v;
% ws = B2C*w;
% 
% x = S2P*xx;
% [xi, ilist, inum, gnum] = createilist(x, ximages, ws, B2C);
% jlist = createjlist(xi, ilist, r, inum);
% [clist, c2ilist, c2inum] = cellist(xi, ws, B2C, inum, gnum);
% jlist2 = createjlist2d(xi, ilist, clist, c2ilist, c2inum, r, inum);
% for i = 1:inum  
%     jlist2{i}-jlist{i}
% end
% 
% figure(1); clf;
% hold on;
% plotline(v(:,1), v(:,2));
% plotline(v(:,2), v(:,3));
% plotline(v(:,3), v(:,4));
% plotline(v(:,4), v(:,1));
% plotline(w(:,1), w(:,2));
% plotline(w(:,2), w(:,3));
% plotline(w(:,3), w(:,4));
% plotline(w(:,4), w(:,1));
% plot(x(1,:),x(2,:),'o');
% 
% figure(2); clf;
% hold on;
% plot(xi(1,1:n),xi(2,1:n),'ob');
% plot(xi(1,n+1:end),xi(2,n+1:end),'or');
% plotline(v(:,1), v(:,2));
% plotline(v(:,2), v(:,3));
% plotline(v(:,3), v(:,4));
% plotline(v(:,4), v(:,1));
% plotline(w(:,1), w(:,2));
% plotline(w(:,2), w(:,3));
% plotline(w(:,3), w(:,4));
% plotline(w(:,4), w(:,1));
% t = linspace(0,2*pi,200);
% m = 50;
% plot(r*cos(t)+xi(1,m),r*sin(t)+xi(2,m),'-');
% jlist2{m}
% axis equal;
% 
% xs = B2C*xi;
% figure(3); clf;
% hold on;
% plot(xs(1,1:n),xs(2,1:n),'ob');
% plot(xs(1,n+1:end),xs(2,n+1:end),'or');
% plotline(vs(:,1), vs(:,2));
% plotline(vs(:,2), vs(:,3));
% plotline(vs(:,3), vs(:,4));
% plotline(vs(:,4), vs(:,1));
% plotline(ws(:,1), ws(:,2));
% plotline(ws(:,2), ws(:,3));
% plotline(ws(:,3), ws(:,4));
% plotline(ws(:,4), ws(:,1));
% 
% ne1 = floor(1./abs(ws(1,1)));
% ne2 = floor(1./abs(ws(2,1)));
% eta1 = [ws(1,1) linspace(0,1,ne1+1) ws(1,2)];
% eta2 = [ws(2,1) linspace(0,1,ne2+1) ws(2,3)];    
% for i = 1:length(eta1)
%     plot([eta1(i) eta1(i)],[eta2(1) eta2(end)],'-');
% end
% for i = 1:length(eta2)
%     plot([eta1(1) eta1(end)],[eta2(i) eta2(i)],'-');
% end
% 
% axis equal;
% 
% % a = [1 0 0];
% % b = [0.2 1 0];
% % c = [0.2 0.2 1];
% % r = 0.2;
% % [w1] = boundingparallelepiped(a, b, c, r);
% 
% % c = a+b;
% % a1 = a(1); a2 = a(2);
% % b1 = b(1); b2 = b(2);
% % 
% % [B2C, S2P] = parallelogram2unitsquare(a, b);
% % 
% % r = 0.1;
% % t = linspace(0,2*pi,200);
% % 
% % [ps,~] = squaremesh(12,12,0,1);
% % pp = S2P*ps;
% % %figure(1); clf;
% % hold on;
% % %plot(pp(1,:), pp(2,:), 'o');
% % plot([0 a(1)], [0 a(2)], 'k-', 'LineWidth', 2); 
% % plot([0 b(1)], [0 b(2)], 'k-', 'LineWidth', 2);
% % plot([a(1) c(1)], [a(2) c(2)], 'k-', 'LineWidth', 2);
% % plot([b(1) c(1)], [b(2) c(2)], 'k-', 'LineWidth', 2);
% % plot(r*cos(t),r*sin(t),'-');
% % axis equal;
% % 
% % % a2*x - a1*y = 0 -> y = a2/a1 * x (1st principal line)
% % % a2*x - a1*y = c a parallel line to 1st principal line
% % % The distance between the two lines is c/norm(a) 
% % c = r*norm(a);
% % 
% % % b2*x - b1*y = 0 -> x = b1/b2 * y (2nd principal line)
% % % b2*x - b1*y = -d a parallel line to 2nd principal line
% % % The distance between the two lines is d/norm(b)
% % d = r*norm(b);
% % 
% % % intersection of y = (a2/a2)*x - c/a1 and x = (b1/b2)*y - d/b2
% % % x = (b1/b2)*((a2/a1)*x - c/a1) - d/b2
% % % (1 - (b1/b2)*((a2/a1)*x = (b1/b2)*(-c/a1) - d/b2
% % % s = (a(2)/a(1))*(b(1)/b(2));
% % % x1 = -((b(1)/b(2))*(c/a(1)) + d/b(2))/(1-s);
% % % y1 = (a(2)/a(1))*x1 - c/a(1);
% % 
% % % intersection of y = (a2/a2)*x - c/a1 + (i-1)*c/a1 and x = (b1/b2)*y - d/b2
% % % intersection of y = (a2/a2)*x - (c-i*c)/a1 and x = (b1/b2)*y - d/b2
% % m = 11;
% % x1 = zeros(1,m); y1 = x1;
% % s = (a(2)/a(1))*(b(1)/b(2));
% % for i = 1:m
% %     ci = c - (i-1)*c;
% %     x1(i) = -((b(1)/b(2))*(ci/a(1)) + d/b(2))/(1-s);
% %     y1(i) = (a(2)/a(1))*x1(i) - ci/a(1);
% % end
% % plot(x1,y1,'-o');
% % 
% % % intersection of y = (a2/a2)*x - c/a1 and x = (b1/b2)*y - d/b2 + (i-1)*d/b2
% % % intersection of y = (a2/a2)*x - (c-i*c)/a1 and x = (b1/b2)*y - d/b2
% % x2 = zeros(1,m); y2 = x2;
% % for i = 1:m
% %     di = d - (i-1)*d;
% %     x2(i) = -((b(1)/b(2))*(c/a(1)) + di/b(2))/(1-s);
% %     y2(i) = (a(2)/a(1))*x2(i) - c/a(1);
% % end
% % plot(x2,y2,'-o');
% % 
% % dx1 = x1(2) - x1(1);
% % dy1 = y1(2) - y1(1);
% % dx2 = x2(2) - x2(1);
% % dy2 = y2(2) - y2(1);
% % c = a+b; % corner opposite of the origin 
% % d = c + [dx1 dy1] + [dx2 dy2];
% % plot(d(1),d(2),'-o');
% % 
% % ds = sqrt(dx1.^2 + dy1.^2);
% % na = floor(norm(a)/ds(1));
% % nb = floor(norm(b)/ds(1));
% % rs = max(norm(a)/na, norm(b)/nb);
% % 
% % % for i = 1:(m+1)
% % %     x = linspace(x1(i),x1(i)+1.3,10);    
% % %     y = (a(2)/a(1))*x - c + (i-1)*c;
% % %     plot(x,y,'-');
% % % end
% % 
% % % y = linspace(y1(1),1,10);
% % % x = (b(1)/b(2))*y - d;
% % % plot(x,y,'-');
% % % for i = 1:10
% % %     x = (b(1)/b(2))*y - d + i*d;
% % %     plot(x,y,'-');
% % % end
% % 
% % 
% % %plot(0.376422+r*cos(t),0.2386395+r*sin(t),'-');
% % 
% % % p1 = p;
% % % p1 = -0.1 + 1.2*p;
% % % figure(1); clf;
% % % plot(p1(1,:), p1(2,:), 'o');
% % % 
% % % 
% % % p2 = S2P*p1;
% % % figure(2); clf;
% % % plot(p1(1,:), p1(2,:), 'o');
% % % hold on;
% % % plot(p2(1,:), p2(2,:), '*');
% % % plot(0.1*cos(linspace(0,2*pi,200)),0.1*sin(linspace(0,2*pi,200)),'-');
% % % plot(p2(1,1:12), p2(2,1:12), '-');
% % % axis equal;
% % % 