function rotatelattice(A, B, C)

Anorm = norm(A);
Bnorm = norm(B);
Cnorm = norm(C);
Ahat = A/Anorm;

ax = Anorm;
bx = dot(B,Ahat);
by = sqrt(Bnorm^2 - bx^2); #norm(cross(Ahat,B));
cx = dot(C,Ahat);
cy = (dot(B, C) - bx*cx)/by;
cz = sqrt(Cnorm^2 - cx^2 - cy^2);

a = [ax 0 0];
b = [bx by 0];
c = [cx cy cz];

return a, b, c

end




