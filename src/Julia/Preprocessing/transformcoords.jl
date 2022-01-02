function transformcoords(a,b,c,A,B,C,X)
    
R = [a[:] b[:] c[:]]*inv([A[:] B[:] C[:]]);

d = size(X,1)
n = size(X,2)
x = R*reshape(X, (d, n));

return x, R

end



