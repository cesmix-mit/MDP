function [X,Y,Z] = rotc(R, x, y, z)

tmp = R*([x(:) y(:) z(:)]');
X = tmp(1,:)';
Y = tmp(2,:)';
Z = tmp(3,:)';
