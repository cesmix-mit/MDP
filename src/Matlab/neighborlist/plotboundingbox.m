function plotboundingbox(v)

v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);
v4 = v(:,4);
plotline(v1, v2);
plotline(v2, v3);
plotline(v3, v4);
plotline(v4, v1);
