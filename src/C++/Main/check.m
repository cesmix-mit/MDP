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

fn = "eicpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "eigpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "ocpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "ogpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "ecpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "egpu.bin";
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

fn = "xijcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "xijgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "eijcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "eijgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "fijcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "fijgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "ecpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "egpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))

fn = "fcpu.bin";
fileID = fopen(fn,'r');
tm1 = fread(fileID,'double');
fclose(fileID);
fn = "fgpu.bin";
fileID = fopen(fn,'r');
tm2 = fread(fileID,'double');
fclose(fileID);
max(abs(tm1(:)-tm2(:)))
