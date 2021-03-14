
for L = 1:10
    n = 100;
    the = pi*rand(n,1);
    phi = -pi + 2*pi*rand(n,1);
    r = ones(n,1);
    [x,y,z] = sphere2cart(the,phi,r);
%     [The,Phi,~] = cart2sphere(x,y,z);
%     [max(abs(the-The)) max(abs(phi-Phi))]
    
    euler = sort([2*pi*rand pi*rand 2*pi*rand]); 
    R = euler2rotm(euler);
    [X,Y,Z] = rotc(R, x, y, z);
    
    [ar,ai] = shsum(x,y,z,L);
    p = shpower(ar,ai);
    b = shbispectrum(ar,ai);

    a2 = shsum2(X,Y,Z,L);
    p2 = shpower2(a2);
    b2 = shbispectrum2(a2);

    fprintf("Euler's angles: [%g %g %g]\n", euler);
    fprintf("Maximum power spectrum error: %g\n", max(abs(p(:)-p2(:))));
    fprintf("Maximum bispectrum error: %g\n", max(abs(b(:)-b2(:))));
end

L = 10;
the = pi*rand;
phi = -pi + 2*pi*rand;
[Ylmr, Ylmi, YlmrThe, YlmiThe, YlmrPhi, YlmiPhi] = sphericalharmonicsderiv(the,phi,L);

dx = 1e-6;
[Ylmr2, Ylmi2] = sphericalharmonics(the+dx,phi,L);
[Ylmr1, Ylmi1] = sphericalharmonics(the-dx,phi,L);
for l=1:(L+1)
    e1 = (Ylmr2{l}-Ylmr1{l})/(2*dx) - YlmrThe{l};
    e2 = (Ylmi2{l}-Ylmi1{l})/(2*dx) - YlmiThe{l};
    fprintf("Maximum derivative error: %g\n", max(sqrt(e1(:).^2+e2(:).^2)));
end

[Ylmr2, Ylmi2] = sphericalharmonics(the,phi+dx,L);
[Ylmr1, Ylmi1] = sphericalharmonics(the,phi-dx,L);
for l=1:(L+1)
    e1 = (Ylmr2{l}-Ylmr1{l})/(2*dx) - YlmrPhi{l};
    e2 = (Ylmi2{l}-Ylmi1{l})/(2*dx) - YlmiPhi{l};
    fprintf("Maximum derivative error: %g\n", max(sqrt(e1(:).^2+e2(:).^2)));
end


