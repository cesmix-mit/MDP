function [v, w] = boundingparallelogram(pbc, a, b, r)
% a = (a1, a2) and b = (b1,b2) are two vectors defining a parallelogram 
% r is the distance from the original parallelogram to the bounding parallelogram 
% w are four vertices defining the bounding parallelogram 

if length(r)==1
    r = [r r];
end

a = a(:); b = b(:);

a1 = a(1); a2 = a(2);
b1 = b(1); b2 = b(2);
norma = norm(a);
normb = norm(b);

% vertices of the parallelogram defined by a and b
v1 = [0; 0];
v2 = a;
v3 = a+b;
v4 = b;

% a2*x - a1*y = 0 -> y = a2/a1 * x (1st principal axis)
% a2*x - a1*y = da a parallel line to 1st principal axis
% The distance between the two lines is da/norma. As da/norma = r, we get
da = r(2)*norma;
% a2*x - a1*y = da is a line which is below and parallel to the 1st axis at a distance r

% apply the same formula to the second axis b
db = r(1)*normb;
% b2*x - b1*y = -db is a line which is left and parallel to the 2nd axis at a distance r

% intersection of a2*x - a1*y = da and b2*x - b1*y = -db
w1 = [a2 -a1; b2 -b1]\[da; -db]; % the 1st vertex of the bounding parallelogram

% find e = (e1,e2) such that e1*a + e2*b = w1
e = [a b]\w1;

% distance between w1 and e(1)*a
sb = norm(w1-e(1)*a);
% distance between w1 and e(2)*b
sa = norm(w1-e(2)*b);

% % length of the bounding parallelogram along the 1st axis
% l1 = norm(a) + 2*sa;
% % length of the bounding parallelogram along the 2nd axis
% l2 = norm(b) + 2*sb;

% length of the bounding parallelogram along the 1st axis
l1 = norm(a) + 2*sa*pbc(1);
% length of the bounding parallelogram along the 2nd axis
l2 = norm(b) + 2*sb*pbc(2);

% the 1st vertex of  the bounding parallelepiped
w1 = pbc(1)*e(1)*a + pbc(2)*e(2)*b;
% the 2nd vertex of  the bounding parallelogram
w2 = w1 + l1*a/norm(a);
% the 3rd vertex of  the bounding parallelogram
w3 = v3-w1;
% the 4th vertex of  the bounding parallelogram
w4 = w1 + l2*b/norm(b);

v = [v1 v2 v3 v4];
w = [w1 w2 w3 w4];

% if (pbc(1)==0) && (pbc(2)==0)
%     w = v;
% elseif (pbc(1)==1) && (pbc(2)==0)
%     w1 = e(1)*a;
%     w3 = v3-w1;
%     sa = norm(w1);
%     l1 = norm(a) + 2*sa;
%     w2 = w1 + l1*a/norm(a);    
%     w4 = w1 + b;
%     w = [w1 w2 w3 w4];
% elseif (pbc(1)==0) && (pbc(2)==1)    
%     w1 = e(2)*b;
%     w2 = w1 + a;
%     w3 = v3-w1;
%     sb = norm(w1);
%     l2 = norm(b) + 2*sb;
%     w4 = w1 + l2*b/norm(b);        
%     w = [w1 w2 w3 w4];
% end

figure(1); clf;
hold on;
plotboundingbox(v);
plotboundingbox(w);
% t = linspace(0,2*pi,200);
% plot(r(1)*cos(t),r(2)*sin(t),'b-','LineWidth',1);
% plot(v3(1)+r(1)*cos(t),v3(2)+r(2)*sin(t),'b-','LineWidth',1);
% plot(a1+r(1)*cos(t),a2+r(2)*sin(t),'b-','LineWidth',1);
% plot(b1+r(1)*cos(t),b2+r(2)*sin(t),'b-','LineWidth',1);
axis equal;


% % c = (c1, c2) corner opposite of the origin 
% c1 = a1 + b1; 
% c2 = a2 + b2;
% 
% % a2*x - a1*y = 0 -> y = a2/a1 * x (1st principal line)
% % a2*x - a1*y = c a parallel line to 1st principal line
% % The distance between the two lines is c/norma 
% c = r*norma;
% 
% % b2*x - b1*y = 0 -> x = b1/b2 * y (2nd principal line)
% % b2*x - b1*y = -d a parallel line to 2nd principal line
% % The distance between the two lines is d/normb
% d = r*normb;
% 
% % intersection of y = (a2/a2)*x - c/a1 and x = (b1/b2)*y - d/b2
% % x = (b1/b2)*((a2/a1)*x - c/a1) - d/b2
% % (1 - (b1/b2)*((a2/a1)*x = (b1/b2)*(-c/a1) - d/b2
% s = (a2/a1)*(b1/b2);
% x1 = -((b1/b2)*(c/a1) + d/b2)/(1-s);
% y1 = (a2/a1)*x1 - c/a1;
% 
% w1
% [x1 y1]
% 
% % the corner opposite to (x1, y1)
% x3 = c1 - x1;
% y3 = c2 - y1;
% 
% w3
% [x3 y3]
% 
% % alpha*a + beta*b = (x1,y1)
% alfa = [a(:) b(:)]\[x1; y1];
% dx1 = -alfa(1)*a1;
% dy1 = -alfa(1)*a2;
% dx2 = -alfa(2)*b1;
% dy2 = -alfa(2)*b2;
% ds1 = sqrt(dx1.^2 + dy1.^2);
% ds2 = sqrt(dx2.^2 + dy2.^2);
% 
% [ds1 norm(w1-e(1)*a)]
% 
% e
% alfa
% 
% l1 = norma + 2*ds1; 
% x2 = x1 + l1*a1/norma;
% y2 = y1 + l1*a2/norma;
% 
% l2 = normb + 2*ds2;
% x4 = x1 + l2*b1/normb;
% y4 = y1 + l2*b2/normb;
% 
% t = linspace(0,2*pi,200);
% figure(2); clf;
% hold on;
% plot(r*cos(t),r*sin(t),'-');
% plot(c1+r*cos(t),c2+r*sin(t),'-');
% plot(a1+r*cos(t),a2+r*sin(t),'-');
% plot(b1+r*cos(t),b2+r*sin(t),'-');
% plot([0 a1], [0 a2], 'k-', 'LineWidth', 2); 
% plot([0 b1], [0 b2], 'k-', 'LineWidth', 2);
% plot([a1 c1], [a2 c2], 'k-', 'LineWidth', 2);
% plot([b1 c1], [b2 c2], 'k-', 'LineWidth', 2);
% plot(x1, y1, 'bo', 'LineWidth', 2);
% plot(x2, y2, 'bo', 'LineWidth', 2);
% plot(x3, y3, 'bo', 'LineWidth', 2);
% plot(x4, y4, 'bo', 'LineWidth', 2);
% plot([x1 x2], [y1 y2], 'k-', 'LineWidth', 2); 
% plot([x2 x3], [y2 y3], 'k-', 'LineWidth', 2); 
% plot([x3 x4], [y3 y4], 'k-', 'LineWidth', 2); 
% plot([x4 x1], [y4 y1], 'k-', 'LineWidth', 2); 
% % plot(alfa(1)*a1, alfa(1)*a2, 'bo', 'LineWidth', 2);
% % plot(alfa(2)*b1, alfa(2)*b2, 'bo', 'LineWidth', 2);
% axis equal;
% 
% 
% 
