/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CONFIGURATION
#define __CONFIGURATION

#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "setunits.cpp"
#include "readdata.cpp"

// template <typename T> void cpuSmallMatrixInverse(T *invA, T *A, int dim)
// {                 
//     if (dim==2) {
//         T detA = A[0]*A[3] - A[1]*A[2];
//         invA[0] = A[3]/detA;
//         invA[1] = -A[1]/detA;
//         invA[2] = -A[2]/detA;
//         invA[3] = A[0]/detA;
//     }
//     else if (dim==3)
//     {
//         T a11 = A[0];
//         T a21 = A[1];
//         T a31 = A[2];
//         T a12 = A[3];
//         T a22 = A[4];
//         T a32 = A[5];
//         T a13 = A[6];
//         T a23 = A[7];
//         T a33 = A[8];        
//         T detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);
// 
//         invA[0] = (a22*a33 - a23*a32)/detA;
//         invA[1] = (a23*a31 - a21*a33)/detA;
//         invA[2] = (a21*a32 - a22*a31)/detA;
//         invA[3] = (a13*a32 - a12*a33)/detA;
//         invA[4] = (a11*a33 - a13*a31)/detA;
//         invA[5] = (a12*a31 - a11*a32)/detA;
//         invA[6] = (a12*a23 - a13*a22)/detA;
//         invA[7] = (a13*a21 - a11*a23)/detA;
//         invA[8] = (a11*a22 - a12*a21)/detA;        
//     }    
// }
// template void cpuSmallMatrixInverse(double*, double*, int);
// template void cpuSmallMatrixInverse(float*, float*, int);
// 
// template <typename T> void cpuSmallMatrixInverse(T *invA, T *A1, T *A2, T *A3)
// {                 
//     T a11 = A1[0];
//     T a21 = A1[1];
//     T a31 = A1[2];
//     T a12 = A2[0];
//     T a22 = A2[1];
//     T a32 = A2[2];
//     T a13 = A3[0];
//     T a23 = A3[1];
//     T a33 = A3[2];        
//     T detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);
// 
//     invA[0] = (a22*a33 - a23*a32)/detA;
//     invA[1] = (a23*a31 - a21*a33)/detA;
//     invA[2] = (a21*a32 - a22*a31)/detA;
//     invA[3] = (a13*a32 - a12*a33)/detA;
//     invA[4] = (a11*a33 - a13*a31)/detA;
//     invA[5] = (a12*a31 - a11*a32)/detA;
//     invA[6] = (a12*a23 - a13*a22)/detA;
//     invA[7] = (a13*a21 - a11*a23)/detA;
//     invA[8] = (a11*a22 - a12*a21)/detA;            
// }
// template void cpuSmallMatrixInverse(double*, double*, double*, double*);
// template void cpuSmallMatrixInverse(float*, float*, float*, float*);
// 
// template <typename T> void BoundingBox3D(T *vc, T *wc, T *v, T *w, T *a, T *b, T *c, T *r, int *pbc)
// {
//     T a1 = a[0]; 
//     T a2 = a[1];
//     T a3 = a[2];
//     T b1 = b[0]; 
//     T b2 = b[1];
//     T b3 = b[2];
//     T c1 = c[0]; 
//     T c2 = c[1];
//     T c3 = c[2];
//     
//     T norma = sqrt(cpuArraySquareSum(a, 3));
//     T normb = sqrt(cpuArraySquareSum(b, 3));
//     T normc = sqrt(cpuArraySquareSum(c, 3));
//     
//     // 8 vertices of the parallelogram defined by a, b, c
//     v[0] = 0.0;
//     v[1] = 0.0;
//     v[2] = 0.0;
//     v[3] = a1;
//     v[4] = a2;
//     v[5] = a3;    
//     v[6] = a1 + b1;
//     v[7] = a2 + b2;
//     v[8] = a3 + b3;    
//     v[9] = b1;
//     v[10] = b2;
//     v[11] = b3;
//     v[12] = c1;
//     v[13] = c2;
//     v[14] = c3;
//     v[15] = a1 + c1;
//     v[16] = a2 + c2;
//     v[17] = a3 + c3;    
//     v[18] = a1 + b1 + c1;
//     v[19] = a2 + b2 + c2;
//     v[20] = a3 + b3 + c3;    
//     v[21] = b1 + c1;
//     v[22] = b2 + c2;
//     v[23] = b3 + c3;
//     
//     // the 1st plane defined by a and b
//     T p1[3] = {a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1};
//     T normp1 = sqrt(cpuArraySquareSum(p1, 3));    
//     // Since the distance between (0,0,0) and the 1st plane p1(1)*x + p1(2)*y + p1(3)*z = d1
//     // is equal to r, we have
//     T d1 = -r[2]*normp1;
//     printf("%g %g %g %g\n", d1, p1[2], r[2], normp1);
//     //[d1 p1(3) r(3) normp1]
//             
//     // the 2nd plane defined by b and c
//     T p2[3] = {b2*c3-b3*c2, b3*c1-b1*c3, b1*c2-b2*c1};
//     T normp2 = sqrt(cpuArraySquareSum(p2, 3));
//     T d2 = -r[0]*normp2;
//     
//     // the 3rd plane defined by c and a
//     T p3[3] = {c2*a3-c3*a2, c3*a1-c1*a3, c1*a2-c2*a1};
//     T normp3 = sqrt(cpuArraySquareSum(p3, 3));
//     T d3 = -r[1]*normp3;
//         
//     printf("%g %g %g\n", d1, d2, d3);
//     
//     // intersection of the above three planes
//     //w1 = [p1; p2; p3]\[d1; d2; d3];
//     T A[9] = {p1[0], p2[0], p3[0], p1[1], p2[1], p3[1], p1[2], p2[2], p3[2]};
//     T invA[9];    
//     cpuSmallMatrixInverse(invA, A, 3);        
//     w[0] = invA[0]*d1 + invA[3]*d2 + invA[6]*d3;
//     w[1] = invA[1]*d1 + invA[4]*d2 + invA[7]*d3;
//     w[2] = invA[2]*d1 + invA[5]*d2 + invA[8]*d3;
//     
//     printf("%g %g %g\n", A[0], A[3], A[6]);
//     printf("%g %g %g\n", A[1], A[4], A[7]);
//     printf("%g %g %g\n", A[2], A[5], A[8]);
//     
//     printf("%g %g %g\n", invA[0], invA[3], invA[6]);
//     printf("%g %g %g\n", invA[1], invA[4], invA[7]);
//     printf("%g %g %g\n", invA[2], invA[5], invA[8]);
//     printf("%g %g %g\n", w[0], w[1], w[2]);
//     
//     // find e = (e1,e2,e3) such that e1*a + e2*b + e3*c = w1
//     //e = [a(:) b(:) c(:)]\w1;    
//     cpuSmallMatrixInverse(invA, a, b, c);
//     T e1 = invA[0]*w[0] + invA[3]*w[1] + invA[6]*w[2];
//     T e2 = invA[1]*w[0] + invA[4]*w[1] + invA[7]*w[2];
//     T e3 = invA[2]*w[0] + invA[5]*w[1] + invA[8]*w[2];
//     
//     T onet = (T) 1.0;
//     T sa[3], sb[3], sc[3];
//     // distance between w1 and the point e2*b + e3*c    
//     cpuArrayAdd3Vectors(sa, w, b, c, onet, -e2, -e3, 3);
//     
//     // distance between w1 and the point e1*a + e3*c
//     cpuArrayAdd3Vectors(sb, w, a, c, onet, -e1, -e3, 3);
// 
//     // distance between w1 and the point e1*a + e2*b
//     cpuArrayAdd3Vectors(sc, w, a, b, onet, -e1, -e2, 3);
//     
//     // length of the bounding parallelepiped along the 1st axis
//     T la = norma + (2.0*pbc[0])*sqrt(cpuArraySquareSum(sa, 3));
//     // length of the bounding parallelepiped along the 2nd axis
//     T lb = normb + (2.0*pbc[1])*sqrt(cpuArraySquareSum(sb, 3));
//     // length of the bounding parallelepiped along the 3rd axis
//     T lc = normc + (2.0*pbc[2])*sqrt(cpuArraySquareSum(sc, 3));
//     
//     // the 1st vertex of the bounding parallelepiped   
//     cpuArrayAdd3Vectors(&w[0], a, b, c, pbc[0]*e1, pbc[1]*e2, pbc[2]*e3, 3);    
//     printf("%g %g %g\n", w[0], w[1], w[2]);
// //     w[0] = pbc[0]*e[0]*a1 + pbc[1]*e[1]*b1 + pbc[c]*e[2]*c1;
// //     w[1] = pbc[0]*e[0]*a2 + pbc[1]*e[1]*b2 + pbc[c]*e[2]*c2;
// //     w[2] = pbc[0]*e[0]*a3 + pbc[1]*e[1]*b3 + pbc[c]*e[2]*c3;    
//     // the 2nd vertex of the bounding parallelepiped: w2 = w1 + la*a/norma;
//     cpuArrayAXPBY(&w[3], &w[0], a, onet, la/norma, 3);
//     // the 4th vertex of the bounding parallelepiped: w4 = w1 + lb*b/normb;
//     cpuArrayAXPBY(&w[9], &w[0], b, onet, lb/normb, 3);
//     // the 3rd vertex of the bounding parallelepiped: w3 = w2 + w4 - w1;
//     cpuArrayAdd3Vectors(&w[6], &w[0], &w[3], &w[9], -onet, onet, onet, 3);    
//     // the 5th vertex of the bounding parallelepiped: w5 = w1 + lc*c/normc;
//     cpuArrayAXPBY(&w[12], &w[0], c, onet, lc/normc, 3);
//     // the 6th vertex of the bounding parallelepiped: w6 = w5 + la*a/norma;
//     cpuArrayAXPBY(&w[15], &w[12], a, onet, la/norma, 3);
//     // the 8th vertex of the bounding parallelepiped: w8 = w5 + lb*b/normb;
//     cpuArrayAXPBY(&w[21], &w[12], b, onet, lb/normb, 3);
//     // the 3rd vertex of the bounding parallelepiped: w7 = w6 + w8 - w5;
//     cpuArrayAdd3Vectors(&w[18], &w[12], &w[15], &w[21], -onet, onet, onet, 3);    
//     
//     for (int i=0; i<8; i++) {
//         vc[3*i+0] = invA[0]*v[3*i+0] + invA[3]*v[3*i+1] + invA[6]*v[3*i+2];
//         vc[3*i+1] = invA[1]*v[3*i+0] + invA[4]*v[3*i+1] + invA[7]*v[3*i+2];
//         vc[3*i+2] = invA[2]*v[3*i+0] + invA[5]*v[3*i+1] + invA[8]*v[3*i+2];
//         wc[3*i+0] = invA[0]*w[3*i+0] + invA[3]*w[3*i+1] + invA[6]*w[3*i+2];
//         wc[3*i+1] = invA[1]*w[3*i+0] + invA[4]*w[3*i+1] + invA[7]*w[3*i+2];
//         wc[3*i+2] = invA[2]*w[3*i+0] + invA[5]*w[3*i+1] + invA[8]*w[3*i+2];        
//     }    
// }
// template void BoundingBox3D(double*, double*, double*, double*, double*, double*, double*, double*, int*);
// template void BoundingBox3D(float*, float*, float*, float*, float*, float*, float*, float*, int*);

// void implReadAppStruct(appstruct &app, commonstruct &common, string filein, Int mpiprocs, Int mpirank, Int backend)
// {
//     string filename = filein + "app.bin";                    
//     
//     // Open file to read
//     ifstream in(filename.c_str(), ios::in | ios::binary);
//     
//     if (!in) 
//         error("Unable to open file " + filename);
//        
//     if (mpirank==0)
//         printf("Read app struct from a binary file...\n");   
//     
//     /* Read data to app structure */            
//     app.lsize = readiarrayfromdouble(in, 1);
//     app.nsize = readiarrayfromdouble(in, app.lsize[0]);
//     app.ndims = readiarrayfromdouble(in, app.nsize[0]);
//     app.flags = readiarrayfromdouble(in, app.nsize[1]);
//     app.bcs = readiarrayfromdouble(in, app.nsize[2]);
//     app.pbc = readiarrayfromdouble(in, app.nsize[3]);
//     readarray(in, &app.boxoffset, app.nsize[4]);  
//     app.atomnumber = readiarrayfromdouble(in, app.nsize[5]);            
//     readarray(in, &app.atommass, app.nsize[6]);    
//     readarray(in, &app.atomcharge, app.nsize[7]);    
//     readarray(in, &app.simulaparam, app.nsize[8]);   
//     readarray(in, &app.solversparam, app.nsize[9]);       
//     readarray(in, &app.eta, app.nsize[10]);       
//     app.kappa = readiarrayfromdouble(in, app.nsize[11]);    
//     readarray(in, &app.muml, app.nsize[12]);       
//     readarray(in, &app.mu1a, app.nsize[13]);       
//     readarray(in, &app.mu1b, app.nsize[14]);       
//     readarray(in, &app.mu2a, app.nsize[15]);       
//     readarray(in, &app.mu2b, app.nsize[16]);       
//     readarray(in, &app.mu2c, app.nsize[17]);       
//     readarray(in, &app.mu3a, app.nsize[18]);       
//     readarray(in, &app.mu3b, app.nsize[19]);       
//     readarray(in, &app.mu3c, app.nsize[20]);       
//     readarray(in, &app.mu4a, app.nsize[21]);       
//     readarray(in, &app.mu4b, app.nsize[22]);       
//     app.pot1a = readiarrayfromdouble(in, app.nsize[23]);    
//     app.pot1b = readiarrayfromdouble(in, app.nsize[24]);    
//     app.pot2a = readiarrayfromdouble(in, app.nsize[25]);    
//     app.pot2b = readiarrayfromdouble(in, app.nsize[26]);    
//     app.pot2c = readiarrayfromdouble(in, app.nsize[27]);    
//     app.pot3a = readiarrayfromdouble(in, app.nsize[28]);    
//     app.pot3b = readiarrayfromdouble(in, app.nsize[29]);    
//     app.pot3c = readiarrayfromdouble(in, app.nsize[30]);    
//     app.pot4a = readiarrayfromdouble(in, app.nsize[31]);    
//     app.pot4b = readiarrayfromdouble(in, app.nsize[32]);      
//     readarray(in, &app.rcutsqml, app.nsize[33]);       
//     readarray(in, &app.rcutsq2a, app.nsize[34]);       
//     readarray(in, &app.rcutsq2b, app.nsize[35]);       
//     readarray(in, &app.rcutsq2c, app.nsize[36]);       
//     readarray(in, &app.rcutsq3a, app.nsize[37]);       
//     readarray(in, &app.rcutsq3b, app.nsize[38]);       
//     readarray(in, &app.rcutsq3c, app.nsize[39]);       
//     readarray(in, &app.rcutsq4a, app.nsize[40]);       
//     readarray(in, &app.rcutsq4b, app.nsize[41]);       
//     readarray(in, &app.rcutsq, app.nsize[42]);       
//     app.atom1b = readiarrayfromdouble(in, app.nsize[43]);    
//     app.atom2b = readiarrayfromdouble(in, app.nsize[44]);    
//     app.atom2c = readiarrayfromdouble(in, app.nsize[45]);    
//     app.atom3b = readiarrayfromdouble(in, app.nsize[46]);    
//     app.atom3c = readiarrayfromdouble(in, app.nsize[47]);    
//     app.atom4b = readiarrayfromdouble(in, app.nsize[48]);        
//     common.traininglist = readiarrayfromdouble(in, app.nsize[50]);        
//     common.validatelist = readiarrayfromdouble(in, app.nsize[51]);        
//     common.trainingnum = app.nsize[50];
//     common.validatenum = app.nsize[51];
//     readarray(in, &app.nveparam, app.nsize[52]);    
//     readarray(in, &app.nvtparam, app.nsize[53]);    
//     
//     int n = 0;
//     for (int i=13; i<=22; i++) 
//         n += app.nsize[i];
//     TemplateMalloc(&app.muep, n, 0); 
//     
//     for (int i=0; i < app.nsize[13]; i++) 
//         app.muep[i] = app.mu1a[i];
//     n = app.nsize[13];    
//     for (int i=0; i < app.nsize[14]; i++) 
//         app.muep[i+n] = app.mu1b[i];
//     n += app.nsize[14];    
//     for (int i=0; i < app.nsize[15]; i++) 
//         app.muep[i+n] = app.mu2a[i];
//     n += app.nsize[15];    
//     for (int i=0; i < app.nsize[16]; i++) 
//         app.muep[i+n] = app.mu2b[i];
//     n += app.nsize[16];    
//     for (int i=0; i < app.nsize[17]; i++) 
//         app.muep[i+n] = app.mu2c[i];
//     n += app.nsize[17];        
//     for (int i=0; i < app.nsize[18]; i++) 
//         app.muep[i+n] = app.mu3a[i];
//     n += app.nsize[18];    
//     for (int i=0; i < app.nsize[19]; i++) 
//         app.muep[i+n] = app.mu3b[i];
//     n += app.nsize[19];    
//     for (int i=0; i < app.nsize[20]; i++) 
//         app.muep[i+n] = app.mu3c[i];
//     n += app.nsize[20];    
//     for (int i=0; i < app.nsize[21]; i++) 
//         app.muep[i+n] = app.mu4a[i];
//     n += app.nsize[21];    
//     for (int i=0; i < app.nsize[22]; i++) 
//         app.muep[i+n] = app.mu4b[i];
//     n += app.nsize[22];    
//     
//     // Close file:
//     in.close();        
// }

void implSetAppStruct(appstruct &app, appstruct &happ, Int backend)
{        
    /* Allocate memory on GPU */            
    TemplateMalloc(&app.nsize, happ.lsize[0], backend); 
    TemplateMalloc(&app.ndims, happ.nsize[0], backend); 
    TemplateMalloc(&app.flags, happ.nsize[1], backend); 
    TemplateMalloc(&app.bcs, happ.nsize[2], backend); 
    TemplateMalloc(&app.pbc, happ.nsize[3], backend); 
    TemplateMalloc(&app.boxoffset, happ.nsize[4], backend);  
    TemplateMalloc(&app.atomnumber, happ.nsize[5], backend); 
    TemplateMalloc(&app.atommass, happ.nsize[6], backend); 
    TemplateMalloc(&app.atomcharge, happ.nsize[7], backend); 
    TemplateMalloc(&app.simulaparam, happ.nsize[8], backend);     
    TemplateMalloc(&app.solversparam, happ.nsize[9], backend);
    TemplateMalloc(&app.eta, happ.nsize[10], backend); 
    TemplateMalloc(&app.kappa, happ.nsize[11], backend); 
    TemplateMalloc(&app.muml, happ.nsize[12], backend); 
    TemplateMalloc(&app.mu1a, happ.nsize[13], backend); 
    TemplateMalloc(&app.mu1b, happ.nsize[14], backend); 
    TemplateMalloc(&app.mu2a, happ.nsize[15], backend); 
    TemplateMalloc(&app.mu2b, happ.nsize[16], backend); 
    TemplateMalloc(&app.mu2c, happ.nsize[17], backend); 
    TemplateMalloc(&app.mu3a, happ.nsize[18], backend); 
    TemplateMalloc(&app.mu3b, happ.nsize[19], backend); 
    TemplateMalloc(&app.mu3c, happ.nsize[20], backend); 
    TemplateMalloc(&app.mu4a, happ.nsize[21], backend); 
    TemplateMalloc(&app.mu4b, happ.nsize[22], backend); 
    TemplateMalloc(&app.pot1a, happ.nsize[23], backend); 
    TemplateMalloc(&app.pot1b, happ.nsize[24], backend); 
    TemplateMalloc(&app.pot2a, happ.nsize[25], backend); 
    TemplateMalloc(&app.pot2b, happ.nsize[26], backend); 
    TemplateMalloc(&app.pot2c, happ.nsize[27], backend); 
    TemplateMalloc(&app.pot3a, happ.nsize[28], backend); 
    TemplateMalloc(&app.pot3b, happ.nsize[29], backend); 
    TemplateMalloc(&app.pot3c, happ.nsize[30], backend); 
    TemplateMalloc(&app.pot4a, happ.nsize[31], backend); 
    TemplateMalloc(&app.pot4b, happ.nsize[32], backend);     
    TemplateMalloc(&app.rcutsqml, happ.nsize[33], backend); 
    TemplateMalloc(&app.rcutsq2a, happ.nsize[34], backend); 
    TemplateMalloc(&app.rcutsq2b, happ.nsize[35], backend); 
    TemplateMalloc(&app.rcutsq2c, happ.nsize[36], backend); 
    TemplateMalloc(&app.rcutsq3a, happ.nsize[37], backend); 
    TemplateMalloc(&app.rcutsq3b, happ.nsize[38], backend); 
    TemplateMalloc(&app.rcutsq3c, happ.nsize[39], backend); 
    TemplateMalloc(&app.rcutsq4a, happ.nsize[40], backend); 
    TemplateMalloc(&app.rcutsq4b, happ.nsize[41], backend); 
    TemplateMalloc(&app.rcutsq, happ.nsize[42], backend);         
    TemplateMalloc(&app.atom1b, happ.nsize[43], backend); 
    TemplateMalloc(&app.atom2b, happ.nsize[44], backend); 
    TemplateMalloc(&app.atom2c, happ.nsize[45], backend); 
    TemplateMalloc(&app.atom3b, happ.nsize[46], backend); 
    TemplateMalloc(&app.atom3c, happ.nsize[47], backend); 
    TemplateMalloc(&app.atom4b, happ.nsize[48], backend);                
    
//     readarray(in, &app.nveparam, app.nsize[52]);    
//     readarray(in, &app.nvtparam, app.nsize[53]);    
//     readarray(in, &app.snaparam, app.nsize[54]);    
//     readarray(in, &app.snapelemradius, app.nsize[55]);    
//     readarray(in, &app.snapelemweight, app.nsize[56]);    
//     readarray(in, &app.snapcoeff, app.nsize[57]);    
//     readarray(in, &app.createvelocity, app.nsize[58]);    
    
    int n = 0;
    for (int i=13; i<=22; i++) 
        n += happ.nsize[i];
    TemplateMalloc(&app.muep, n, backend); 
    
    
    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
        CUDA_CHECK( cudaMemcpy(app.nsize, happ.nsize, happ.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.ndims, happ.ndims, happ.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.flags, happ.flags, happ.nsize[1]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.bcs, happ.bcs, happ.nsize[2]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pbc, happ.pbc, happ.nsize[3]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.boxoffset, happ.boxoffset, happ.nsize[4]*sizeof(dstype), cudaMemcpyHostToDevice ) );                                      
        CUDA_CHECK( cudaMemcpy(app.atomnumber, happ.atomnumber, happ.nsize[5]*sizeof(Int), cudaMemcpyHostToDevice ) );                         
        CUDA_CHECK( cudaMemcpy(app.atommass, happ.atommass, happ.nsize[6]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.atomcharge, happ.atomcharge, happ.nsize[7]*sizeof(dstype), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.simulaparam, happ.simulaparam, happ.nsize[8]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.solversparam, happ.solversparam, happ.nsize[9]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.eta, happ.eta, happ.nsize[10]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.kappa, happ.kappa, happ.nsize[11]*sizeof(Int), cudaMemcpyHostToDevice ) );   
        CUDA_CHECK( cudaMemcpy(app.muml, happ.muml, happ.nsize[12]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu1a, happ.mu1a, happ.nsize[13]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu1b, happ.mu1b, happ.nsize[14]*sizeof(dstype), cudaMemcpyHostToDevice ) );    
        CUDA_CHECK( cudaMemcpy(app.mu2a, happ.mu2a, happ.nsize[15]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu2b, happ.mu2b, happ.nsize[16]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu2c, happ.mu2c, happ.nsize[17]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu3a, happ.mu3a, happ.nsize[18]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu3b, happ.mu3b, happ.nsize[19]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu3c, happ.mu3c, happ.nsize[20]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu4a, happ.mu4a, happ.nsize[21]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.mu4b, happ.mu4b, happ.nsize[22]*sizeof(dstype), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.pot1a, happ.pot1a, happ.nsize[23]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot1b, happ.pot1b, happ.nsize[24]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot2a, happ.pot2a, happ.nsize[25]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot2b, happ.pot2b, happ.nsize[26]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot2c, happ.pot2c, happ.nsize[27]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot3a, happ.pot3a, happ.nsize[28]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot3b, happ.pot3b, happ.nsize[29]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot3c, happ.pot3c, happ.nsize[30]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot4a, happ.pot4a, happ.nsize[31]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.pot4b, happ.pot4b, happ.nsize[32]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.rcutsqml, happ.rcutsqml, happ.nsize[33]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq2a, happ.rcutsq2a, happ.nsize[34]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq2b, happ.rcutsq2b, happ.nsize[35]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq2c, happ.rcutsq2c, happ.nsize[36]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq3a, happ.rcutsq3a, happ.nsize[37]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq3b, happ.rcutsq3b, happ.nsize[38]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq3c, happ.rcutsq3c, happ.nsize[39]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq4a, happ.rcutsq4a, happ.nsize[40]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.rcutsq4b, happ.rcutsq4b, happ.nsize[41]*sizeof(dstype), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.rcutsq, happ.rcutsq, happ.nsize[42]*sizeof(dstype), cudaMemcpyHostToDevice ) );                          
        CUDA_CHECK( cudaMemcpy(app.atom1b, happ.atom1b, happ.nsize[43]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.atom2b, happ.atom2b, happ.nsize[44]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.atom2c, happ.atom2c, happ.nsize[45]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.atom3b, happ.atom3b, happ.nsize[46]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(app.atom3c, happ.atom3c, happ.nsize[47]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(app.atom4b, happ.atom4b, happ.nsize[48]*sizeof(Int), cudaMemcpyHostToDevice ) );                                      
        CUDA_CHECK( cudaMemcpy(app.muep, happ.muep, n*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif        
    }        
}

// void implReadConfigStruct(configstruct &config, string filein, Int mpiprocs, Int mpirank, Int backend)
// {
//     string filename;
//     if (mpiprocs > 1) {
//         Int filenumber = mpirank+1; //file number     
//         filename = filein + "config" + NumberToString(filenumber) + ".bin";                    
//     }
//     else
//         filename = filein + "config" + ".bin";                    
//     
//     // Open file to read
//     ifstream in(filename.c_str(), ios::in | ios::binary);
//     
//     if (!in) 
//         error("Unable to open file " + filename);
//        
//     if (mpirank==0)
//         printf("Read config struct from a binary file...\n");   
//     
//     /* Read data to config structure */            
//     config.lsize = readiarrayfromdouble(in, 1);
//     config.nsize = readiarrayfromdouble(in, config.lsize[0]);
//     config.ndims = readiarrayfromdouble(in, config.nsize[0]);
//     config.natoms = readiarrayfromdouble(in, config.nsize[1]);        
//     readarray(in, &config.a, config.nsize[2]);   
//     readarray(in, &config.b, config.nsize[3]);   
//     readarray(in, &config.c, config.nsize[4]);   
//     readarray(in, &config.e, config.nsize[5]);   
//     config.t = readiarrayfromdouble(in, config.nsize[6]);
//     readarray(in, &config.x, config.nsize[7]);   
//     readarray(in, &config.q, config.nsize[8]);   
//     readarray(in, &config.v, config.nsize[9]);   
//     readarray(in, &config.f, config.nsize[10]);      
//     readarray(in, &config.we, config.nsize[11]);      
//     readarray(in, &config.wf, config.nsize[12]);      
//                 
//     TemplateMalloc(&config.natomssum, config.nsize[1]+1, 0); 
//     cpuCumsum(config.natomssum, config.natoms, config.nsize[1]+1); 
//         
//     // Close file:
//     in.close();            
// }
// 
// void implSetCommonStruct(commonstruct &common, appstruct &app, configstruct &config, string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
// {
//     common.filein = filein;
//     common.fileout = fileout;
//     common.mpiProcs = mpiprocs;
//     common.mpiRank = mpirank;
//     common.backend = backend;
//     
//     common.dim = app.ndims[0];
//     common.L = app.ndims[1];  // order of radial basis functions
//     common.K = app.ndims[2];  // order of spherical harmonics          
//     common.ntimesteps = app.ndims[3]; // number of time steps
//     common.nab = app.ndims[4]; // number of atoms per block
//     common.natomtypes = app.ndims[8]; // number of atom types
//     common.nmoletypes = app.ndims[9]; // number of molecule types
//     common.neighmax = app.ndims[10]; // % maximum number of neighbors allowed
//     
//     common.descriptor = app.flags[0];   // descriptor flag: 0 -> Spherical Harmonics Bessel
//     common.spectrum = app.flags[1];     // spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
//     common.training = app.flags[2];     // 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
//     common.runMD = app.flags[3];        // 0 no MD simulation, 1 -> run MD simulation
//     common.potential = app.flags[4];    // 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential        
//     common.neighpair = app.flags[5];  // 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
//     common.energycal = app.flags[6];    // turns energy calculation on or off
//     common.forcecal = app.flags[7];     // turns force calculation on or off
//     common.stresscal = app.flags[8];   // turns stress calculation on or off
//     common.neighcell = app.flags[9];   // 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
//     common.decomposition = app.flags[10]; // 0 -> force decomposition, 1 -> atom decomposition
//     common.chemtype = app.flags[11];      // 0 -> single atom-type basis functions, 1 -> double atom-type basis functions 
//     common.dftdata = app.flags[12];       // 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
//     common.unitstyle = app.flags[13];       // 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
//     common.ensemblemode = app.flags[14];   
//     
//     // set unit parameters
//     setunits(common, common.unitstyle);
//         
//     common.time = app.simulaparam[0]; // current time    
//     common.dt = app.simulaparam[1];    
//     common.currentstep = (int) common.time/common.dt;
//     common.rcutml = sqrt(app.rcutsqml[0]);
//     
//     common.nflags = app.nsize[1];
//     common.nsimulaparam = app.nsize[8];   
//     common.nsolversparam = app.nsize[9];       
//     common.neta = app.nsize[10];  
//     common.nkappa = app.nsize[11];  
//     common.nmuml = app.nsize[12];  
//     common.nmu1a = app.nsize[13];  
//     common.nmu1b = app.nsize[14];  
//     common.nmu2a = app.nsize[15];  
//     common.nmu2b = app.nsize[16];  
//     common.nmu2c = app.nsize[17];  
//     common.nmu3a = app.nsize[18];  
//     common.nmu3b = app.nsize[19];  
//     common.nmu3c = app.nsize[20];  
//     common.nmu4a = app.nsize[21];  
//     common.nmu4b = app.nsize[22];      
//     common.npot1a = app.nsize[23];  
//     common.npot1b = app.nsize[24];  
//     common.npot2a = app.nsize[25];  
//     common.npot2b = app.nsize[26];  
//     common.npot2c = app.nsize[27];  
//     common.npot3a = app.nsize[28];  
//     common.npot3b = app.nsize[29];  
//     common.npot3c = app.nsize[30];  
//     common.npot4a = app.nsize[31];  
//     common.npot4b = app.nsize[32];                  
//     common.natom1b = app.nsize[43];  
//     common.natom2b = app.nsize[44];  
//     common.natom2c = app.nsize[45];  
//     common.natom3b = app.nsize[46];  
//     common.natom3c = app.nsize[47];  
//     common.natom4b = app.nsize[48];  
//     common.nmu[0] = 0;
//     common.nmu[1] = common.nmu1a;
//     common.nmu[2] = common.nmu[1] + common.nmu1b;
//     common.nmu[3] = common.nmu[2] + common.nmu2a;
//     common.nmu[4] = common.nmu[3] + common.nmu2b;
//     common.nmu[5] = common.nmu[4] + common.nmu2c;
//     common.nmu[6] = common.nmu[5] + common.nmu3a;
//     common.nmu[7] = common.nmu[6] + common.nmu3b;
//     common.nmu[8] = common.nmu[7] + common.nmu3c;
//     common.nmu[9] = common.nmu[8] + common.nmu4a;
//     common.nmu[10] = common.nmu[9] + common.nmu4b;
//     common.Nempot = common.npot1a + common.npot1b;
//     common.Nempot += common.npot2a + common.npot2b + common.npot2c;
//     common.Nempot += common.npot3a + common.npot3b + common.npot3c;
//     common.Nempot += common.npot4a + common.npot4b;
//     
//     common.nconfigs = config.nsize[1]; // number of configurations
//     common.ne = config.nsize[5];
//     common.nt = config.nsize[6];
//     common.nx = config.nsize[7];
//     common.nq = config.nsize[8];
//     common.nv = config.nsize[9];
//     common.nf = config.nsize[10];
//     common.ncq = common.nq/(common.nt);        
//                 
//     common.inum = config.natoms[0];    
//     common.nlocal = common.inum;    
//     common.inummax = config.natoms[0];        
//     for (int i=1; i<common.nconfigs; i++)         
//         if (config.natoms[i] > common.inummax)
//             common.inummax = config.natoms[i];   
//     
// #ifdef HAVE_MPI        
//     MPI_Allreduce(&common.nlocal, &common.natoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
// #else    
//     common.natoms = common.nlocal;
// #endif    
//     common.tdof = (common.natoms-1)*common.dim;    
//             
//     for (int i=0; i<common.dim; i++) {       
//         common.pbc[i] = app.pbc[i];
//         common.boxoffset[i] = app.boxoffset[i];
//     }
//     
//     common.dtarray[0] = common.dt;     // dt
//     common.dtarray[1] = 0.5*common.dt; // dtf
//     common.dtarray[2] = common.dt;     // dtv
//     common.dtarray[3] = 0;             // beginstep 
//     common.dtarray[4] = common.ntimesteps; // endstep
//     common.dtarray[5] = common.currentstep;      
//     
//     if (common.ensemblemode == 1) { // NVE limit ensemble                
//       common.vlimitsq = (app.nveparam[0]/common.dt) * (app.nveparam[0]/common.dt);
//     }
//     
//     if (common.ensemblemode == 2) { // NVT ensemble                
//         common.tarray[0] = app.nvtparam[0]; // temp start
//         common.tarray[1] = app.nvtparam[1]; // temp stop
//         common.tarray[2] = 1/app.nvtparam[2]; // temp frequency    
//         common.tarray[5] = common.tdof; 
//         common.tarray[6] = common.boltz; 
//         common.tarray[7] = (app.nsize[53]>3) ? app.nvtparam[3] : 0.0; // drag factor added to barostat/thermostat         
//         common.tarray[8] = common.mvv2e; 
//         common.mtchain   = (app.nsize[53]>4) ? (int) app.nvtparam[4] : 3;         
//         common.nc_tchain = (app.nsize[53]>5) ? (int) app.nvtparam[5] : 1;                 
//     }
//     
//     TemplateMalloc(&common.pot1a, app.nsize[23], 0); 
//     TemplateMalloc(&common.pot1b, app.nsize[24], 0); 
//     TemplateMalloc(&common.pot2a, app.nsize[25], 0); 
//     TemplateMalloc(&common.pot2b, app.nsize[26], 0); 
//     TemplateMalloc(&common.pot2c, app.nsize[27], 0); 
//     TemplateMalloc(&common.pot3a, app.nsize[28], 0); 
//     TemplateMalloc(&common.pot3b, app.nsize[29], 0); 
//     TemplateMalloc(&common.pot3c, app.nsize[30], 0); 
//     TemplateMalloc(&common.pot4a, app.nsize[31], 0); 
//     TemplateMalloc(&common.pot4b, app.nsize[32], 0);         
//     TemplateMalloc(&common.atom1b, app.nsize[43], 0); 
//     TemplateMalloc(&common.atom2b, app.nsize[44], 0); 
//     TemplateMalloc(&common.atom2c, app.nsize[45], 0); 
//     TemplateMalloc(&common.atom3b, app.nsize[46], 0); 
//     TemplateMalloc(&common.atom3c, app.nsize[47], 0); 
//     TemplateMalloc(&common.atom4b, app.nsize[48], 0);        
//     
//     cpuArrayCopy(common.pot1a, app.pot1a, app.nsize[23]);
//     cpuArrayCopy(common.pot1b, app.pot1b, app.nsize[24]);
//     cpuArrayCopy(common.pot2a, app.pot2a, app.nsize[25]);
//     cpuArrayCopy(common.pot2b, app.pot2b, app.nsize[26]);
//     cpuArrayCopy(common.pot2c, app.pot2c, app.nsize[27]);
//     cpuArrayCopy(common.pot3a, app.pot3a, app.nsize[28]);
//     cpuArrayCopy(common.pot3b, app.pot3b, app.nsize[29]);
//     cpuArrayCopy(common.pot3c, app.pot3c, app.nsize[30]);
//     cpuArrayCopy(common.pot4a, app.pot4a, app.nsize[31]);
//     cpuArrayCopy(common.pot4b, app.pot4b, app.nsize[32]);    
//     cpuArrayCopy(common.atom1b, app.atom1b, app.nsize[43]);
//     cpuArrayCopy(common.atom2b, app.atom2b, app.nsize[44]);
//     cpuArrayCopy(common.atom2c, app.atom2c, app.nsize[45]);
//     cpuArrayCopy(common.atom3b, app.atom3b, app.nsize[46]);
//     cpuArrayCopy(common.atom3c, app.atom3c, app.nsize[47]);
//     cpuArrayCopy(common.atom4b, app.atom4b, app.nsize[48]);    
// }

void implSetAtomBlocks(commonstruct &common) 
{
    Int inum = common.inum;
    Int ns = min(inum, common.nab);  // number of atoms per block
    common.nba = min((Int) round(inum/ns), 16);     // number of blocks    
    Int na = (Int) round(inum/common.nba); // number of atoms per block    
        
    for (int i=0; i<=common.nba; i++)
        common.ablks[i] = i*na;
    common.ablks[common.nba] = inum;
    
    common.nabmax = 0;
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        na = e2 - e1; 
        if (na > common.nabmax)
            common.nabmax = na;
    }    
}

void implGetConfiguration(Int *atomtype, dstype *x, commonstruct &common, configstruct &config, Int ci)
{
    Int inum = config.natoms[ci];    
    Int start = config.natomssum[ci];    
    common.inum = inum;    
    
    if (common.backend <= 1) {               
        cpuArrayCopy(x, &config.x[common.dim*start], common.dim*inum);
        cpuArrayCopy(atomtype, &config.t[start], inum);
    }
    else {
#ifdef HAVE_CUDA                
        CUDA_CHECK( cudaMemcpy(atomtype, &config.t[start], inum*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CUDA_CHECK( cudaMemcpy(x, &config.x[common.dim*start], common.dim*inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif        
    }        
}

void implGetAtomtypes(Int *atomtype, commonstruct &common, configstruct &config, Int ci)
{
    common.inum = config.natoms[ci];    
    Int start = config.natomssum[ci];        
    if (common.backend <= 1) {               
        cpuArrayCopy(atomtype, &config.t[start], common.inum);
    }
    else {
#ifdef HAVE_CUDA                
        CUDA_CHECK( cudaMemcpy(atomtype, &config.t[start], common.inum*sizeof(Int), cudaMemcpyHostToDevice ) );                      
#endif        
    }        
}

void implGetPositions(dstype *x, commonstruct &common, configstruct &config, Int ci)
{
    common.inum = config.natoms[ci];    
    Int start = config.natomssum[ci];      
    if (common.backend <= 1) {
        cpuArrayCopy(x, &config.x[common.dim*start], common.dim*common.inum);    
    }
    else {
#ifdef HAVE_CUDA                        
        CUDA_CHECK( cudaMemcpy(x, &config.x[common.dim*start], common.dim*common.inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif    
    }
}

void implGetVelocities(dstype *v, commonstruct &common, configstruct &config, Int ci)
{
    Int start = config.natomssum[ci];      
    if (common.backend <= 1) {
        cpuArrayCopy(v, &config.v[common.dim*start], common.dim*common.inum);    
    }
    else {
#ifdef HAVE_CUDA                        
        CUDA_CHECK( cudaMemcpy(v, &config.v[common.dim*start], common.dim*common.inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif    
    }
}

void implGetForces(dstype *f, commonstruct &common, configstruct &config, Int ci)
{
    Int start = config.natomssum[ci];        
    if (common.backend <= 1) {
        cpuArrayCopy(f, &config.f[common.dim*start], common.dim*common.inum);    
    }
    else {
#ifdef HAVE_CUDA                                
        CUDA_CHECK( cudaMemcpy(f, &config.f[common.dim*start], common.dim*common.inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif    
    }
}

void implGetCharges(dstype *q, commonstruct &common, configstruct &config, Int ci)
{
    Int start = config.natomssum[ci]; 
    if (common.backend <= 1) {
        cpuArrayCopy(q, &config.q[common.dim*start], common.ncq*common.inum);    
    }
    else {
#ifdef HAVE_CUDA                                
        CUDA_CHECK( cudaMemcpy(q, &config.q[common.dim*start], common.ncq*common.inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif            
    }        
}

void implGetEnergy(dstype *e, commonstruct &common, configstruct &config, Int ci)
{    
    if (common.backend <= 1) {
        cpuArrayCopy(e, &config.e[ci], 1);    
    }
    else {
#ifdef HAVE_CUDA                                
        CUDA_CHECK( cudaMemcpy(e, &config.e[ci], sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif                    
    }
}

void implSetNeighborStruct(neighborstruct &nb, commonstruct &common, configstruct &config, Int ci)
{               
    Int dim = common.dim;    
    Int m = (dim==2) ? dim*4 : dim*8;         
    Int n = (dim==2) ? dim*9 : dim*27;         
        
    Int ne[dim], cellnum[3];
    dstype a[dim], b[dim], c[dim], cellsize[dim], s2rmap[dim*dim];        
    dstype refvertices[m], rbvertices[m], boxvertices[m], bbvertices[m];
    dstype pimages[n];
    
    dstype simbox[dim*dim];    
    for (Int i=0; i<dim; i++) {
        simbox[i] = config.a[i + dim*dim*ci];
        simbox[dim+i] = config.b[i + dim*dim*ci];
        if (dim > 2)
            simbox[2*dim+i] = config.c[i + dim*dim*ci];
    }                    
    cpuSmallMatrixInverse(s2rmap, simbox, dim);       
    for (Int i=0; i<dim; i++) {
        a[i] = simbox[i];
        b[i] = simbox[dim+i];        
        c[i] = (dim > 2) ? simbox[2*dim+i] : 0.0;        
    }                    
        
    if (dim==2) {
        cpuBoundingBox2D(refvertices, rbvertices, boxvertices, bbvertices, 
                a, b, common.boxoffset, common.pbc);
        
        common.pnum = cpuPeriodicImages2D(pimages, a, b, common.pbc);                
    }
    else if (dim==3){
        cpuBoundingBox3D(refvertices, rbvertices, boxvertices, bbvertices, 
                a, b, c, common.boxoffset, common.pbc);        
//         BoundingBox3D(refvertices, rbvertices, boxvertices, bbvertices, 
//                 a, b, c, common.boxoffset, common.pbc);
        
        common.pnum = cpuPeriodicImages3D(pimages, a, b, c, common.pbc);                
    }    
    
//     printArray2D(refvertices, 3, 8, common.backend);  
//     printArray2D(rbvertices, 3, 8, common.backend);  
//     printArray2D(boxvertices, 3, 8, common.backend);  
//     printArray2D(bbvertices, 3, 8, common.backend);  
//     printArray2D(pimages, 3, common.pnum, common.backend);  
//     error("here");
    
    dstype *smin = &rbvertices[0];
    dstype *smax = (dim == 2) ? &rbvertices[dim*2] :  &rbvertices[dim*6];        
    common.cnum = 1;
    for (Int i=0; i<dim; i++) {        
        ne[i] = (Int) round(1.0/(smax[i]-1.0));   
        ne[i] = (ne[i] < 1) ? 1 : ne[i];            
        cellnum[i] = ne[i] + 2;                
        cellsize[i] = 1.0/((dstype) ne[i]);
        common.cnum = common.cnum*cellnum[i];
    }
    if (dim==3) 
        cellnum[2] = 0;    
    
    dstype eta1[cellnum[0]+1], eta2[cellnum[1]+1], eta3[cellnum[2]+1];
    cpuMakeReferenceGrid(eta1, smin[0], smax[0], cellsize[0], cellnum[0]);
    cpuMakeReferenceGrid(eta2, smin[1], smax[1], cellsize[1], cellnum[1]);
    common.rbbvol = (smax[0]-smin[0])*(smax[1]-smin[1]);
    if (dim==3) {
        common.rbbvol = common.rbbvol*(smax[2]-smin[2]);
        cpuMakeReferenceGrid(eta3, smin[2], smax[2], cellsize[2], cellnum[2]);    
    }
    common.anummax = (Int) ceil(1.2*common.rbbvol*common.inum);
            
    nb.freememory(common.backend);        
    TemplateMalloc(&nb.cellnum, dim, common.backend);         
    TemplateMalloc(&nb.cellsize, dim, common.backend);         
    TemplateMalloc(&nb.a, dim, common.backend); 
    TemplateMalloc(&nb.b, dim, common.backend); 
    TemplateMalloc(&nb.c, dim, common.backend);    
    TemplateMalloc(&nb.s2rmap, dim*dim, common.backend);              
    TemplateMalloc(&nb.refvertices, m, common.backend); 
    TemplateMalloc(&nb.rbvertices,  m, common.backend); 
    TemplateMalloc(&nb.boxvertices, m, common.backend); 
    TemplateMalloc(&nb.bbvertices, m, common.backend);    
    TemplateMalloc(&nb.pimages, n, common.backend);    
    TemplateMalloc(&nb.eta1, cellnum[0]+1, common.backend);   
    TemplateMalloc(&nb.eta2, cellnum[1]+1, common.backend);   
    TemplateMalloc(&nb.eta3, cellnum[2]+1, common.backend);  
            
    common.inum = config.natoms[ci];    
    TemplateMalloc(&nb.atomtype, common.inum, common.backend);  
    implGetAtomtypes(nb.atomtype, common, config, ci);                
    
    // need to check size of nb.alist
    TemplateMalloc(&nb.alist, common.anummax, common.backend);  
    TemplateMalloc(&nb.neighlist, common.inum*common.neighmax, common.backend);  
    TemplateMalloc(&nb.neighnum, common.inum, common.backend);      
    //TemplateMalloc(&nb.neighnumsum, common.inum+1, common.backend);      
    //for (Int i=0; i<common.inum*common.neighmax; i++)
    //    nb.neighlist[i] = 0;
            
    if (common.backend <= 1) {
        for (Int i=0; i<dim; i++)
            nb.cellnum[i] = cellnum[i];
        cpuArrayCopy(nb.cellsize, cellsize, dim);
        cpuArrayCopy(nb.a, a, dim);
        cpuArrayCopy(nb.b, b, dim);
        cpuArrayCopy(nb.c, c, dim);
        cpuArrayCopy(nb.s2rmap, s2rmap, dim*dim);
        cpuArrayCopy(nb.refvertices, refvertices, m);
        cpuArrayCopy(nb.rbvertices, rbvertices, m);
        cpuArrayCopy(nb.boxvertices, boxvertices, m);
        cpuArrayCopy(nb.bbvertices, bbvertices, m);
        cpuArrayCopy(nb.pimages, pimages, n);
        cpuArrayCopy(nb.eta1, eta1, cellnum[0]+1);
        cpuArrayCopy(nb.eta2, eta2, cellnum[1]+1);
        cpuArrayCopy(nb.eta3, eta3, cellnum[2]+1);
    }
    else if (common.backend==2) { // GPU
#ifdef HAVE_CUDA        
        CUDA_CHECK( cudaMemcpy(nb.cellnum, cellnum, dim*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.cellsize, cellsize, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.a, a, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.b, b, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.c, c, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.s2rmap, s2rmap, dim*dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.refvertices, refvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.rbvertices, rbvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.boxvertices, boxvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.bbvertices, bbvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.pimages, pimages, n*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.eta1, eta1, (cellnum[0]+1)*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.eta2, eta2, (cellnum[1]+1)*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CUDA_CHECK( cudaMemcpy(nb.eta3, eta3, (cellnum[2]+1)*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif        
    }           
}

void implSetTempStruct(tempstruct & tmp, commonstruct &common, snastruct &sna, shstruct &sh)
{
    // need to check size of tmp.intmem and tmp.tmpmem
    Int anummax= common.anummax;
    Int inum = common.inum;
    Int neighmax = common.neighmax;
    Int pnum = common.pnum;
    
    Int nabmax = inum; // common.nabmax;
    int ijnum = nabmax*neighmax;
        
    Int n1 = pnum*inum + 2*inum + 1;
    Int n2 = 4*anummax;    
    Int n3 = 1 + 3*nabmax + 11*nabmax*neighmax; 
    common.nintmem = max(max(n1, n2), n3);
    TemplateMalloc(&tmp.intmem, common.nintmem, common.backend);      
    
    //printf("%i %i %i %i\n", inum, nabmax, neighmax, pnum);
    
    n1 = ijnum*(7*common.dim+3*common.ncq+1);
    common.ntmpmem = n1;    
    if (common.descriptor == 0) {
        Int nbasis = common.K*(common.L+1)*(common.L+1);        
        n2 = 2*nabmax*nbasis + ijnum*(8*nbasis) + nabmax*common.Nbf;
        n3 = 2*nabmax*nbasis + ijnum*(6*nbasis) + 3*nabmax*neighmax*common.Nbf;    
        common.ntmpmem = max(max(n1, n2), n3);    
    } else if (common.descriptor == 1) {
        int dim = common.dim;
        n2 = ijnum*dim;
        n2 += max(max(sna.idxu_max*ijnum, sna.idxz_max*sna.ndoubles*nabmax), sna.idxu_max*sna.nelements*nabmax); 
        n2 += max(max(sna.idxu_max*ijnum, sna.idxz_max*sna.ndoubles*nabmax), sna.idxu_max*sna.nelements*nabmax); 
        n2 += sna.idxu_max*dim*ijnum;
        n2 += sna.idxu_max*dim*ijnum;    
        n2 += max(2*sna.idxu_max*sna.nelements*nabmax + sna.ncoeff*nabmax + sna.idxb_max*sna.ntriples*nabmax, dim*ijnum);        
        
        n2 = ijnum*dim;
        n2 += 4*sna.idxu_max*sna.nelements*nabmax;
        common.ntmpmem = max(n1, n2);   
    }
    
    printf("Temporary memory for integer array: %i\n", common.nintmem);
    printf("Temporary memory for float array: %i\n", common.ntmpmem);
    
    TemplateMalloc(&tmp.tmpmem, common.ntmpmem, common.backend);          
}

void implSetSysStruct(sysstruct & sys, commonstruct &common, configstruct &config, Int ci)
{   
    TemplateMalloc(&sys.e, common.inummax, common.backend);      
    TemplateMalloc(&sys.f, common.dim*common.inummax, common.backend);      
    TemplateMalloc(&sys.vatom, 6*common.inummax, common.backend);      
    TemplateMalloc(&sys.image, common.dim*common.inummax, common.backend);  
    ArraySetValue(sys.image, 0, common.dim*common.inummax, common.backend);  
    
    //printf("%i %i %i %i\n", common.dim, common.inum, common.inummax, common.anummax);
    if (common.nx>0) {
        // need to check size of sys.x        
        TemplateMalloc(&sys.x, common.dim*max(common.anummax,common.inummax), common.backend);  
        implGetPositions(sys.x, common, config, ci);
    }
    
    TemplateMalloc(&sys.xhold, common.dim*common.inummax, common.backend);  
    ArrayCopy(sys.xhold, sys.x, common.dim*common.inummax, common.backend);  
    
    if (common.nv > 0) {
        TemplateMalloc(&sys.v, common.dim*common.inummax, common.backend);  
        implGetVelocities(sys.v, common, config, ci);
    }
    
    if (common.nq > 0) {
        TemplateMalloc(&sys.q, common.ncq*common.inummax, common.backend);  
        implGetCharges(sys.q, common, config, ci);           
    }
}

void implNeighborList(neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, Int inum)
{
            
    Int dim = common.dim;
    Int cnum = common.cnum;
    Int pnum = common.pnum;
    Int neighmax = common.neighmax;    
    Int ntype = common.natomtypes;
    Int gnum, anum;            
    
    //printf("%i %i %i %i %i %i\n", dim, cnum, pnum, inum, neighmax, ntype);
    
    Int *inside = &tmp.intmem[0]; // inum*pnum
    Int *glistnum = &tmp.intmem[inum*common.pnum]; // inum
    Int *glistnumsum = &tmp.intmem[inum*common.pnum + inum]; // inum+1
    Int *d_sums = &tmp.intmem[inum*common.pnum + 2*inum+1]; // inum+1
    Int *d_incr = &tmp.intmem[inum*common.pnum + 3*inum+2]; // inum+1
    
    if (dim==2) {        
        // include periodic images of atoms I as ghost atoms
        AtomList2D(nb.alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, nb.pimages, 
                 nb.rbvertices, nb.s2rmap, inum, pnum, dim, common.backend);
        
        // number of ghost atoms
        gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
        anum = inum + gnum;
        common.inum = inum;
        common.gnum = gnum;
        common.anum = anum;
        
        if (common.neighcell == 0) { // O(N^2) algorithm to form the neighbor list                    
            FullNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, anum, inum, neighmax, dim, common.backend);                                    
        }
        else { // form neighbor list using cell list
            Int *clist = &tmp.intmem[0]; //anum
            Int *c2alist = &tmp.intmem[anum]; // anum
            Int *c2anum = &tmp.intmem[2*anum]; // anum
            Int *c2anumsum = &tmp.intmem[3*anum]; // anum            
            Int *dtemp = &tmp.intmem[4*anum]; 
            
            CellList2D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, inum, anum, dim, common.backend);  
            Cell2AtomList(c2alist, c2anumsum, c2anum, clist, dtemp, anum, cnum, common.backend);                        
            FullNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist,   
                       c2alist, c2anumsum, nb.cellnum, inum, neighmax, dim, common.backend);                             
        }
    }
    else {        
        AtomList3D(nb.alist, inside, glistnumsum, glistnum, d_sums, d_incr, x, nb.pimages, 
                 nb.rbvertices, nb.s2rmap, inum, pnum, dim, common.backend);
        
        gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
        anum = inum + gnum;
        common.inum = inum;
        common.gnum = gnum;
        common.anum = anum;                               
        
        if (common.neighcell == 0) { // O(N^2) algorithm to form the neighbor list
            FullNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, anum, inum, neighmax, dim, common.backend);                                    
        }
        else { // form neighbor list using cell list
            Int *clist = &tmp.intmem[0]; //anum
            Int *c2alist = &tmp.intmem[anum]; // anum
            Int *c2anum = &tmp.intmem[2*anum]; // anum
            Int *c2anumsum = &tmp.intmem[3*anum]; // anum            
            Int *dtemp = &tmp.intmem[4*anum]; 
            
            CellList3D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, inum, anum, dim, common.backend);  
            Cell2AtomList(c2alist, c2anumsum, c2anum, clist, dtemp, anum, cnum, common.backend);                        
            FullNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist,   
                       c2alist, c2anumsum, nb.cellnum, inum, neighmax, dim, common.backend);                 
        }
    }
    
// #ifdef HAVE_DEBUG                      
//     //printArray2D(nb.neighnum, 1, inum, common.backend);  
//     //printArray2D(nb.neighlist, neighmax, inum, common.backend);              
//     string fn = (common.backend == 2) ? "alistgpu.bin" : "alistcpu.bin";
//     writearray2file(fn, nb.alist, anum, common.backend); 
//     fn = (common.backend == 2) ? "neighnumgpu.bin" : "neighnumcpu.bin";
//     writearray2file(fn, nb.neighnum, inum, common.backend);
//     fn = (common.backend == 2) ? "neighlistgpu.bin" : "neighlistcpu.bin";
//     writearray2file(fn, nb.neighlist, inum*neighmax, common.backend);    
// #endif                        
    
    //error("here");

    if (anum > common.anummax)
        error("Memory allocation for ghost atoms is insffucient");    
}

// void implReadInputFiles(appstruct &app, configstruct &config, commonstruct &common, 
//         string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
// {
//     implReadAppStruct(app, common, filein, mpiprocs, mpirank, backend);    
//     implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);            
//     implSetCommonStruct(common, app, config, filein,  fileout, mpiprocs, mpirank, backend);          
// }

void implSetConfiguration(neighborstruct &nb, tempstruct &tmp, sysstruct &sys, appstruct &app, configstruct &config, commonstruct &common, Int ci)
{
    implSetNeighborStruct(nb, common, config, ci);         
    implSetSysStruct(sys, common, config, ci);       
    implSetAtomBlocks(common);    
}

#endif

