function [v,w] = boundingparallelepiped(pbc, a, b, c, r)
% a = (a1, a2, a3), b = (b1,b2,b3), c=(c1,c2,c3) are three vectors defining a parallelepiped 
% r is the distance from the original parallelepiped to the bounding parallelepiped 
% w are 8 vertices defining the bounding parallelepiped 

if length(r)==1
    r = [r r r];
end

a = a(:); b = b(:); c = c(:);

a1 = a(1); a2 = a(2); a3 = a(3);
b1 = b(1); b2 = b(2); b3 = b(3);
c1 = c(1); c2 = c(2); c3 = c(3);
norma = norm(a);
normb = norm(b);
normc = norm(c);

% vertices of the parallelogram defined by a, b, and c
v1 = [0; 0; 0];
v2 = a;
v3 = a+b;
v4 = b;
v5 = c;
v6 = a+c;
v7 = a+b+c;
v8 = b+c;

% the 1st plane defined by a and b
p1 = [a2*b3-a3*b2 a3*b1-a1*b3 a1*b2-a2*b1];
normp1 = sqrt(p1(1)*p1(1) + p1(2)*p1(2) + p1(3)*p1(3));

% Since the distance between (0,0,0) and the 1st plane p1(1)*x + p1(2)*y + p1(3)*z = d1
% is equal to r(3), we have
d1 = -sign(p1(3))*r(3)*normp1;

% the 2nd plane defined by b and c
p2 = [b2*c3-b3*c2 b3*c1-b1*c3 b1*c2-b2*c1];
normp2 = sqrt(p2(1)*p2(1) + p2(2)*p2(2) + p2(3)*p2(3));

% Since the distance between (0,0,0) and the 2nd plane p2(1)*x + p2(2)*y + p2(3)*z = d2
% is equal to r(1), we have
d2 = -sign(p2(1))*r(1)*normp2;

% the 3rd plane defined by c and a
p3 = [c2*a3-c3*a2 c3*a1-c1*a3 c1*a2-c2*a1];
normp3 = sqrt(p3(1)*p3(1) + p3(2)*p3(2) + p3(3)*p3(3));

% Since the distance between (0,0,0) and the 2nd plane p3(1)*x + p3(2)*y + p3(3)*z = d3
% is equal to r(2), we have
d3 = -sign(p3(2))*r(2)*normp3;

% intersection of the above three planes
w1 = [p1; p2; p3]\[d1; d2; d3];

% find e = (e1,e2,e3) such that e1*a + e2*b + e3*c = w1
e = [a(:) b(:) c(:)]\w1;

% distance between w1 and the point e2*b + e3*c
sa = norm(w1 - e(2)*b - e(3)*c);

% distance between w1 and the point e1*a + e3*c
sb = norm(w1 - e(1)*a - e(3)*c);

% distance between w1 and the point e1*a + e2*b
sc = norm(w1 - e(1)*a - e(2)*b);

% % distance between w1 and the 2nd plane p2 defined by b and c
% sa1 = abs(w1(1)*p2(1)+w1(2)*p2(2)+w1(3)*p2(3))/normp2;
% 
% % distance between w1 and the 3rd plane p3 defined by c and a
% sb1 = abs(w1(1)*p3(1)+w1(2)*p3(2)+w1(3)*p3(3))/normp3;
% 
% % distance between w1 and the 1st plane p3 defined by a and b
% sc1 = abs(w1(1)*p1(1)+w1(2)*p1(2)+w1(3)*p1(3))/normp1;
% 
% [sa sa1 sb sb1 sc sc1]

% % length of the bounding parallelepiped along the 1st axis
% la = norma + 2*sa;
% % length of the bounding parallelepiped along the 2nd axis
% lb = normb + 2*sb;
% % length of the bounding parallelepiped along the 3rd axis
% lc = normc + 2*sc;

% length of the bounding parallelepiped along the 1st axis
la = norma + 2*sa*pbc(1);
% length of the bounding parallelepiped along the 2nd axis
lb = normb + 2*sb*pbc(2);
% length of the bounding parallelepiped along the 3rd axis
lc = normc + 2*sc*pbc(3);

% the 1st vertex of  the bounding parallelepiped
w1 = pbc(1)*e(1)*a + pbc(2)*e(2)*b + pbc(3)*e(3)*c;

% the 2nd vertex of  the bounding parallelepiped
w2 = w1 + la*a/norma;
% the 4th vertex of  the bounding parallelepiped
w4 = w1 + lb*b/normb;
% the 3rd vertex of  the bounding parallelepiped
w3 = w2 + w4 - w1;
% the 5th vertex of  the bounding parallelepiped
w5 = w1 + lc*c/normc;
% the 6th vertex of  the bounding parallelepiped
w6 = w5 + la*a/norma;
% the 8th vertex of  the bounding parallelepiped
w8 = w5 + lb*b/normb;
% the 3rd vertex of  the bounding parallelepiped
w7 = w6 + w8 - w5;

v = [v1 v2 v3 v4 v5 v6 v7 v8];
w = [w1 w2 w3 w4 w5 w6 w7 w8];    

% if (pbc(1)==0) && (pbc(2)==0) && (pbc(3)==0)
%     w = v;
% elseif (pbc(1)==1) && (pbc(2)==0) && (pbc(3)==0)
%     
% elseif (pbc(1)==0) && (pbc(2)==1) && (pbc(3)==0)
%     
% elseif (pbc(1)==0) && (pbc(2)==0) && (pbc(3)==1)
%     
% elseif (pbc(1)==1) && (pbc(2)==1) && (pbc(3)==0)
%     
% elseif (pbc(1)==1) && (pbc(2)==0) && (pbc(3)==1)
%     
% elseif (pbc(1)==0) && (pbc(2)==1) && (pbc(3)==1)
%     
% elseif (pbc(1)==1) && (pbc(2)==1) && (pbc(3)==1)
%     w = [w1 w2 w3 w4 w5 w6 w7 w8];    
% end


figure(1); clf;
hold on;
plotline(v1, v2);
plotline(v2, v3);
plotline(v3, v4);
plotline(v4, v1);
plotline(v5, v6);
plotline(v6, v7);
plotline(v7, v8);
plotline(v8, v5);
plotline(v1, v5);
plotline(v2, v6);
plotline(v3, v7);
plotline(v4, v8);

plotline(w1, w2);
plotline(w2, w3);
plotline(w3, w4);
plotline(w4, w1);
plotline(w5, w6);
plotline(w6, w7);
plotline(w7, w8);
plotline(w8, w5);
plotline(w1, w5);
plotline(w2, w6);
plotline(w3, w7);
plotline(w4, w8);
%plot3(w1(1),w1(2),w1(3),'o','LineWidth',2);
axis equal
view(3);

