fn = "alistcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "alistgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "neighnumcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "neighnumgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "neighlistcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "neighlistgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "clistcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "clistgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "c2anumcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "c2anumgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "c2anumsumcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "c2anumsumgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "c2alistcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "c2alistgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "neighnumcpu2.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
fn = "neighnumgpu2.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "neighlistcpu2.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'int');
fclose(fileID);
tm1 = reshape(tm1,app.maxnumneighbors,[]);
tm1 = sort(tm1,1);
fn = "neighlistgpu2.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'int');
fclose(fileID);
tm2 = reshape(tm2,app.maxnumneighbors,[]);
tm2 = sort(tm2,1);
max(abs(tm1(:)-tm2(:)))

fn = "xijcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "xijgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "srcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "srgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "sicpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "sigpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "dcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "dgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "ddcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "ddgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "ccpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "cgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "bcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "bgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "acpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "agpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

xij = reshape(tm1,3,[]);

x = xij(1,:);
y = xij(2,:);
z = xij(3,:);
r = sqrt(x.*x + y.*y + z.*z);
phi = atan2(y,x);                                
r2 = r.*r;
rxy = sqrt(x.*x + y.*y);
rxy2 = rxy.*rxy;
rr2 = rxy.*r2;
Rx = x./r;
Ry = y./r;
Rz = z./r;
Thex = x.*z./rr2;
They = y.*z./rr2;
Thez = -rxy./r2;
Phix = -y./rxy2;
Phiy = x./rxy2;
Phiz = 0.0;        

costhe = z./r; 
l=1; m= 0;
C = sqrt((2*l+1)*factorial(l-m)/(4*pi*factorial(l+m)));
Ylmr1 = C*costhe;

l=1; m= 1;
C = sqrt((2*l+1)*factorial(l-m)/(4*pi*factorial(l+m)));
Ylmr2 = C*cos(phi).*(-sqrt((r-z).*(r+z))./r);

