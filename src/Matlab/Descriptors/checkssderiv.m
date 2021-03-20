L = 3;
K = 3;

x = [-1e-10, 1e-10, 1e-10, -1e-10]/1e10;
y = [-1e-10, -1e-10, 1e-10, 1e-10]/1e10;
z = [1, 1, 1, 1];
[Sr, Si, Srx, Sry, Srz, Six, Siy, Siz] = ssharmonicsderiv(x,y,z,K,L);

for i = 1:L+1
    Srx{i}
    Sry{i}
    Srz{i}
    Six{i}
    Siy{i}
    Siz{i}
end