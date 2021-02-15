#ifndef __CPUNEIGHBORLIST
#define __CPUNEIGHBORLIST

// template <typename T> void cpuBoundingBox2D(T *vc, T*wc, T *v, T *w, T *a, T *b, T *r, int *pbc)
// {
//     T a1 = a[0]; 
//     T a2 = a[1];
//     T b1 = b[0]; 
//     T b2 = b[1];
//     T norma = sqrt(cpuArraySquareSum(a, 2));
//     T normb = sqrt(cpuArraySquareSum(b, 2));
//     
//     // vertices of the parallelogram defined by a and b
//     v[0] = 0.0;
//     v[1] = 0.0;
//     v[2] = a1;
//     v[3] = a2;
//     v[4] = a1 + b1;
//     v[5] = a2 + b2;
//     v[6] = b1;
//     v[7] = b2;
//     
//     // Find da such that the distance from the origin (0,0)
//     // to the line a2*x - a1*y = da is equal to r
//     T da = r[1]*norma;
//     // Find db such that the distance from the origin (0,0)
//     // to the line b2*x - b1*y = -db is equal to r
//     T db = r[0]*normb;
//     
//     // intersection of a2*x - a1*y = da and b2*x - b1*y = -db
//     T A[4] = {a2, b2, -a1, -b1};
//     T invA[4];    
//     cpuSmallMatrixInverse(invA, A, 2);    
//     w[0] = invA[0]*da - invA[2]*db;
//     w[1] = invA[1]*da - invA[3]*db;
//         
//     // find e = (e1,e2) such that e1*a + e2*b = w1
//     A[0] = a1;
//     A[1] = a2;
//     A[2] = b1;
//     A[3] = b2;
//     cpuSmallMatrixInverse(invA, A, 2);    
//     T e[2];    
//     e[0] = invA[0]*w[0] + invA[2]*w[1];
//     e[1] = invA[1]*w[0] + invA[3]*w[1];    
//     
//     // distance between w1 and e(1)*a
//     T sb = sqrt((w[0]-e[0]*a[0])*(w[0]-e[0]*a[0]) + (w[1]-e[0]*a[1])*(w[1]-e[0]*a[1]));
//     // distance between w1 and e(2)*a
//     T sa = sqrt((w[0]-e[1]*b[0])*(w[0]-e[1]*b[0]) + (w[1]-e[1]*b[1])*(w[1]-e[1]*b[1]));
//    
//     // length of the bounding parallelogram along the 1st axis
//     T l1 = norma + 2*sa*(pbc[0]);
//     // length of the bounding parallelogram along the 2nd axis
//     T l2 = normb + 2*sb*(pbc[1]);
//     
//     // the 1st vertex of  the bounding parallelepiped
//     w[0] = pbc[0]*e[0]*a1 + pbc[1]*e[1]*b1;
//     w[1] = pbc[0]*e[0]*a2 + pbc[1]*e[1]*b2;
//     
//     // the 2nd vertex of  the bounding parallelogram
//     w[2] = w[0] + l1*a[0]/norma;
//     w[3] = w[1] + l1*a[1]/norma;
//     
//     // the 3rd vertex of  the bounding parallelogram
//     w[4] = v[4] - w[0];
//     w[5] = v[5] - w[1];
//     
//     // the 4th vertex of  the bounding parallelogram
//     w[6] = w[0] + l2*b[0]/normb;
//     w[7] = w[1] + l2*b[1]/normb;    
//     
//     // bounding box in the reference domain
//     for (int i=0; i<4; i++) {
//         vc[2*i+0] = invA[0]*v[2*i+0] + invA[2]*v[2*i+1];
//         vc[2*i+1] = invA[1]*v[2*i+0] + invA[3]*v[2*i+1];
//         wc[2*i+0] = invA[0]*w[2*i+0] + invA[2]*w[2*i+1];
//         wc[2*i+1] = invA[1]*w[2*i+0] + invA[3]*w[2*i+1]; 
//     }        
// }
// template void cpuBoundingBox2D(double*, double*, double*, double*, double*, double*, double*, int*);
// template void cpuBoundingBox2D(float*, float*, float*, float*, float*, float*, float*, int*);
// 
// template <typename T> void cpuBoundingBox3D(T *vc, T *wc, T *v, T *w, T *a, T *b, T *c, T *r, int *pbc)
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
//     // intersection of the above three planes
//     //w1 = [p1; p2; p3]\[d1; d2; d3];
//     T A[9] = {p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]};
//     T invA[9];    
//     cpuSmallMatrixInverse(invA, A, 3);        
//     w[0] = invA[0]*d1 + invA[3]*d2 + invA[6]*d3;
//     w[1] = invA[1]*d1 + invA[4]*d2 + invA[7]*d3;
//     w[2] = invA[2]*d1 + invA[5]*d2 + invA[8]*d3;
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
// template void cpuBoundingBox3D(double*, double*, double*, double*, double*, double*, double*, double*, int*);
// template void cpuBoundingBox3D(float*, float*, float*, float*, float*, float*, float*, float*, int*);
// 
// template <typename T> int cpuPeriodicImages2D(T *pimages, T *a, T *b, int *pbc)
// {
//     T px[2][9];
//     
//     px[0][0] = 0.0;             px[1][0] = 0.0;     
//     px[0][1] = a[0];            px[1][1] = a[1];    
//     px[0][2] =-a[0];            px[1][2] =-a[1];    
//     px[0][3] = b[0];            px[1][3] = b[1];      
//     px[0][4] =-b[0];            px[1][4] =-b[1];      
//     px[0][5] = a[0]+b[0];       px[1][5] = a[1]+b[1]; 
//     px[0][6] = a[0]-b[0];       px[1][6] = a[1]-b[1]; 
//     px[0][7] =-a[0]+b[0];       px[1][7] =-a[1]+b[1]; 
//     px[0][8] =-a[0]-b[0];       px[1][8] =-a[1]-b[1]; 
//     
//     int nimages;
//     int indx[9];
//     if ((pbc[0]==0) && (pbc[1]==0)) {
//          nimages = 1;
//          int ind[] = {0};
//          for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==1) && (pbc[1]==0)) {
//         nimages = 3;
//         int ind[] = {0,1,2};      
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==0) && (pbc[1]==1)) {
//         nimages = 3;
//         int ind[] = {0,3,4};    
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==1) && (pbc[1]==1)) {
//         nimages = 9;
//         int ind[] = {0,1,2,3,4,5,6,7,8};
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else {
//     }
//     
//     for (int i=0; i<nimages; i++) {
//         pimages[2*i+0] = px[0][indx[i]];
//         pimages[2*i+1] = px[1][indx[i]];
//     }    
//     
//     return nimages;
// }
// template int cpuPeriodicImages2D(double*, double*, double*, int*);
// template int cpuPeriodicImages2D(float*, float*, float*, int*);
    
// template <typename T> int cpuPeriodicImages3D(T *pimages, T *a, T *b, T *c, int *pbc)
// {
//     T px[3][27];
//     
//     px[0][0] = 0.0;             px[1][0] = 0.0;             px[2][0] = 0.0;
//     px[0][1] = a[0];            px[1][1] = a[1];            px[2][1] = a[2];
//     px[0][2] =-a[0];            px[1][2] =-a[1];            px[2][2] =-a[2];
//     px[0][3] = b[0];            px[1][3] = b[1];            px[2][3] = b[2];
//     px[0][4] =-b[0];            px[1][4] =-b[1];            px[2][4] =-b[2];
//     px[0][5] = a[0]+b[0];       px[1][5] = a[1]+b[1];       px[2][5] = a[2]+b[2];
//     px[0][6] = a[0]-b[0];       px[1][6] = a[1]-b[1];       px[2][6] = a[2]-b[2];
//     px[0][7] =-a[0]+b[0];       px[1][7] =-a[1]+b[1];       px[2][7] =-a[2]+b[2];
//     px[0][8] =-a[0]-b[0];       px[1][8] =-a[1]-b[1];       px[2][8] =-a[2]-b[2];
//     
//     px[0][9] =  c[0]+0.0;       px[1][9] =  c[1]+0.0;       px[2][9] =  c[2]+0.0;
//     px[0][10] = c[0]+a[0];      px[1][10] = c[1]+a[1];      px[2][10] = c[2]+a[2];
//     px[0][11] = c[0]-a[0];      px[1][11] = c[1]-a[1];      px[2][11] = c[2]-a[2];
//     px[0][12] = c[0]+b[0];      px[1][12] = c[1]+b[1];      px[2][12] = c[2]+b[2];
//     px[0][13] = c[0]-b[0];      px[1][13] = c[1]-b[1];      px[2][13] = c[2]-b[2];
//     px[0][14] = c[0]+a[0]+b[0]; px[1][14] = c[1]+a[1]+b[1]; px[2][14] = c[2]+a[2]+b[2];
//     px[0][15] = c[0]+a[0]-b[0]; px[1][15] = c[1]+a[1]-b[1]; px[2][15] = c[2]+a[2]-b[2];
//     px[0][16] = c[0]-a[0]+b[0]; px[1][16] = c[1]-a[1]+b[1]; px[2][16] = c[2]-a[2]+b[2];
//     px[0][17] = c[0]-a[0]-b[0]; px[1][17] = c[1]-a[1]-b[1]; px[2][17] = c[2]-a[2]-b[2];
// 
//     px[0][18] =-c[0]+0.0;       px[1][18] =-c[1]+0.0;       px[2][18] =-c[2]+0.0;
//     px[0][19] =-c[0]+a[0];      px[1][19] =-c[1]+a[1];      px[2][19] =-c[2]+a[2];
//     px[0][20] =-c[0]-a[0];      px[1][20] =-c[1]-a[1];      px[2][20] =-c[2]-a[2];
//     px[0][21] =-c[0]+b[0];      px[1][21] =-c[1]+b[1];      px[2][21] =-c[2]+b[2];
//     px[0][22] =-c[0]-b[0];      px[1][22] =-c[1]-b[1];      px[2][22] =-c[2]-b[2];
//     px[0][23] =-c[0]+a[0]+b[0]; px[1][23] =-c[1]+a[1]+b[1]; px[2][23] =-c[2]+a[2]+b[2];
//     px[0][24] =-c[0]+a[0]-b[0]; px[1][24] =-c[1]+a[1]-b[1]; px[2][24] =-c[2]+a[2]-b[2];
//     px[0][25] =-c[0]-a[0]+b[0]; px[1][25] =-c[1]-a[1]+b[1]; px[2][25] =-c[2]-a[2]+b[2];
//     px[0][26] =-c[0]-a[0]-b[0]; px[1][26] =-c[1]-a[1]-b[1]; px[2][26] =-c[2]-a[2]-b[2];
//     
//     int nimages;
//     int indx[27];
//     if ((pbc[0]==0) && (pbc[1]==0) && (pbc[2]==0)) {
//          nimages = 1;
//          int ind[] = {0};
//          for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==1) && (pbc[1]==0) && (pbc[2]==0)) {
//         nimages = 3;
//         int ind[] = {0,1,2};      
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==0) && (pbc[1]==1) && (pbc[2]==0)) {
//         nimages = 3;
//         int ind[] = {0,3,4};    
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==0) && (pbc[1]==0) && (pbc[2]==1)) {
//         nimages = 3;
//         int ind[] = {0,9,18};
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==1) && (pbc[1]==1) && (pbc[2]==0)) {
//         nimages = 9;
//         int ind[] = {0,1,2,3,4,5,6,7,8};
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==1) && (pbc[1]==0) && (pbc[2]==1)) {
//         nimages = 9;
//         int ind[] = {0,1,2,9,10,11,18,19,20};
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==0) && (pbc[1]==1) && (pbc[2]==1)) {
//         nimages = 9;
//         int ind[] = {0,3,4,9,12,13,18,21,22};
//         for (int i=0; i<nimages; i++)
//             indx[i] = ind[i];
//     }
//     else if ((pbc[0]==1) && (pbc[1]==1) && (pbc[2]==1)) {
//         nimages = 27;        
//         for (int i=0; i<nimages; i++)
//             indx[i] = i;
//     }
//     else {
//     }
//     
//     for (int i=0; i<nimages; i++) {
//         pimages[3*i+0] = px[0][indx[i]];
//         pimages[3*i+1] = px[1][indx[i]];
//         pimages[3*i+2] = px[2][indx[i]];
//     }    
//     
//     return nimages;
// }
// template int cpuPeriodicImages3D(double*, double*, double*, double*, int*);
// template int cpuPeriodicImages3D(float*, float*, float*, float*, int*);
// 
// template <typename T> void cpuMakeReferenceGrid(T *eta, T smin, T smax, int nc, int pbc)
// {
//     if (pbc==1) { // periodic boundary
//         eta[0] = smin;
//         eta[nc] = smax;
//         for (int i=1; i<nc; i++)                   
//             eta[i] = (i-1)*(((T) 1.0)/((T) (nc-2)));
//     }
//     else {
//         for (int i=0; i<=nc; i++)                   
//             eta[i] = i*(((T) 1.0)/((T) nc));
//     }    
// }
// template void cpuMakeReferenceGrid(double*, double, double, int, int);
// template void cpuMakeReferenceGrid(float*, float, float, int, int);
// 
// /************************************************************************************/
// 
// void cpuGridColoring(int *cellcolor, int *nc, int *bsize, int dim)
// {
//     if (dim==2) {
//         int nc1 = nc[0];
//         int nc2 = nc[1];
//         int bs1 = bsize[0];
//         int bs2 = bsize[1];
//         int nb1 = ceil(((double) nc1)/((double) bs1)); // number of blocks
//         int nb2 = ceil(((double) nc2)/((double) bs2)); // number of blocks
//         
//         int bid[bs1*bs2];
//         for (int i=0; i<bs1*bs2; i++)            
//              bid[i] = i;
//         
//         for (int i=0; i<nb1; i++) {
//             int i1, i2;
//             if (i==1) 
//                 i1 = (i==1) ? 0 : (i1 + bs1 - 1);
//             i2 = i1 + bs1 - 1;    
//             if (i2 >= nc1)
//                 i2 = nc1 - 1;
//             for (int j=0; j<nb2; j++) {
//                 int j1, j2;
//                 if (j == 1) 
//                     j1 = (j==1) ? 0 : (j1 + bs2 - 1);
//                 j2 = j1 + bs2 - 1;    
//                 if (j2 >= nc2)
//                     j2 = nc2 - 1;              
//                 for (int m2=j1; m2<=j2; m2++)
//                     for (int m1=i1; m1<=i2; m1++)                    
//                         cellcolor[m1 + nc[0]*m2] = bid[m1-i1 + (m2-j1)*bs1];
//             }
//         }                        
//     }
//     else if (dim==3) {
//         int nc1 = nc[0];
//         int nc2 = nc[1];
//         int nc3 = nc[2];
//         int bs1 = bsize[0];
//         int bs2 = bsize[1];
//         int bs3 = bsize[2];        
//         int nb1 = ceil(((double) nc1)/((double) bs1)); // number of blocks
//         int nb2 = ceil(((double) nc2)/((double) bs2)); // number of blocks
//         int nb3 = ceil(((double) nc3)/((double) bs3)); // number of blocks
//         
//         int bid[bs1*bs2*bs3];
//         for (int i=0; i<bs1*bs2*bs3; i++)            
//              bid[i] = i;
//         
//         for (int i=0; i<nb1; i++) {
//             int i1, i2;
//             if (i==1) 
//                 i1 = (i==1) ? 0 : (i1 + bs1 - 1);
//             i2 = i1 + bs1 - 1;    
//             if (i2 >= nc1)
//                 i2 = nc1 - 1;
//             for (int j=0; j<nb2; j++) {
//                 int j1, j2;
//                 if (j==1) 
//                     j1 = (j==1) ? 0 : (j1 + bs2 - 1);
//                 j2 = j1 + bs2 - 1;    
//                 if (j2 >= nc2)
//                     j2 = nc2 - 1;                              
//                 for (int k=0; k<nb3; k++) {
//                     int k1, k2;
//                     if (k==1) 
//                         k1 = (k==1) ? 0 : (k1 + bs3 - 1);
//                     k2 = k1 + bs3 - 1;    
//                     if (k2 >= nc3)
//                         k2 = nc3 - 1;          
//                     
//                     for (int m3=k1; m3<=k2; m3++)
//                         for (int m2=j1; m2<=j2; m2++)
//                             for (int m1=i1; m1<=i2; m1++)                    
//                                 cellcolor[m1 + nc[0]*m2 + nc[0]*nc[1]*m3] = bid[m1-i1 + (m2-j1)*bs1 + (m3-k1)*bs1*bs2];                
//                 }
//             }
//         }                                
//     }
// }
// 
// template <typename T> int cpuCreateIlist(int *ilist, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
// {
//     int k = 0;
//     
//     if (dim==2) {
//         for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
//             ilist[i] = i;         // add atom i to the list
//             for (int j=1; j<m; j++) { // loop over each periodic image of atom i
//                 T xj0 = x[i*dim+0] + pimages[j*dim+0];  // periodic image of x      
//                 T xj1 = x[i*dim+1] + pimages[j*dim+1];        
//                 T xc0 = B2C[0]*xj0 + B2C[2]*xj1;        // map it to the unit square      
//                 T xc1 = B2C[1]*xj0 + B2C[3]*xj1;        
//                 /// check if the mapped point is inside the bounding box
//                 if ((wc[0] <= xc0) && (xc0 <= wc[4]) &&  (wc[1] <= xc1) && (xc1 <= wc[5])) {
//                     x[dim*(n+k)+0] = xj0; // add the periodic image as a ghost atom
//                     x[dim*(n+k)+1] = xj1; //
//                     ilist[n+k] = i;       // add atom i to the list
//                     k = k + 1;            // increase the number of ghost atoms by 1    
//                 }
//             }
//         }    
//     }
//     else if (dim==3){
//         for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
//             ilist[i] = i; // add atom i to the list
//             for (int j=1; j<m; j++) { // loop over each periodic image of atom i
//                 T xj0 = x[i*dim+0] + pimages[j*dim+0];    // periodic image of x          
//                 T xj1 = x[i*dim+1] + pimages[j*dim+1];        
//                 T xj2 = x[i*dim+2] + pimages[j*dim+2];        
//                 T xc0 = B2C[0]*xj0 + B2C[3]*xj1 + B2C[6]*xj2;  // map it to the unit square            
//                 T xc1 = B2C[1]*xj0 + B2C[4]*xj1 + B2C[7]*xj2;        
//                 T xc2 = B2C[2]*xj0 + B2C[5]*xj1 + B2C[8]*xj2;   
//                 /// check if the mapped point is inside the bounding box
//                 if ((wc[0] <= xc0) && (xc0 <= wc[18]) && (wc[1] <= xc1) && (xc1 <= wc[19]) && (wc[2] <= xc2) && (xc2 <= wc[20])) {
//                     x[dim*(n+k)+0] = xj0; // add the periodic image as a ghost atom
//                     x[dim*(n+k)+1] = xj1;  
//                     x[dim*(n+k)+2] = xj2;
//                     ilist[n+k] = i;      // add atom i to the list
//                     k = k + 1;           // increase the number of ghost atoms by 1              
//                 }
//             }
//         }            
//     }    
//     return k;
// }
// template int cpuCreateIlist(int*, double*, double*, double*, double*, int, int, int);
// template int cpuCreateIlist(int*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuGhostAtoms(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{    
    if (dim==2) {
        for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
            glistnum[i] = 0; // set the number of ghost atoms for atom i to 0
            for (int j=1; j<m; j++) { // loop over each periodic image of atom i
                T xj0 = x[i*dim+0] + pimages[j*dim+0];  // periodic image of x      
                T xj1 = x[i*dim+1] + pimages[j*dim+1];        
                T xc0 = B2C[0]*xj0 + B2C[2]*xj1;        // map it to the unit square      
                T xc1 = B2C[1]*xj0 + B2C[3]*xj1;                        
                /// check if the mapped point is inside the bounding box
                if ((wc[0] <= xc0) && (xc0 <= wc[4]) &&  (wc[1] <= xc1) && (xc1 <= wc[5])) 
                    glistnum[i] += 1;            // increase the number of ghost atoms by 1                    
            }
        }    
    }
    else if (dim==3){
        for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
            glistnum[i] = 0;  // set the number of ghost atoms for atom i to 0
            for (int j=1; j<m; j++) { // loop over each periodic image of atom i
                T xj0 = x[i*dim+0] + pimages[j*dim+0];    // periodic image of x          
                T xj1 = x[i*dim+1] + pimages[j*dim+1];        
                T xj2 = x[i*dim+2] + pimages[j*dim+2];        
                T xc0 = B2C[0]*xj0 + B2C[3]*xj1 + B2C[6]*xj2;  // map it to the unit square            
                T xc1 = B2C[1]*xj0 + B2C[4]*xj1 + B2C[7]*xj2;        
                T xc2 = B2C[2]*xj0 + B2C[5]*xj1 + B2C[8]*xj2;   
                /// check if the mapped point is inside the bounding box
                if ((wc[0] <= xc0) && (xc0 <= wc[18]) && (wc[1] <= xc1) && (xc1 <= wc[19]) && (wc[2] <= xc2) && (xc2 <= wc[20])) 
                    glistnum[i] += 1;            // increase the number of ghost atoms by 1                    
            }
        }            
    }        
}
template void cpuGhostAtoms(int*, double*, double*, double*, double*, int, int, int);
template void cpuGhostAtoms(int*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuCreateAtomList(int *alist, int *glistnumsum, int *glistnum, int *atomtype,
        T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{           
    // a list contains the starting position of the ghost atom of every atom i
    cpuCumsum(glistnumsum, glistnum, n+1); 
    
    if (dim==2) {
        for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box            
            alist[i] = i;         // add atom i to the list
            int q = n + glistnumsum[i]; // offset the starting position by n
            int k = 0;            
            for (int j=1; j<m; j++) { // loop over each periodic image of atom i
                T xj0 = x[i*dim+0] + pimages[j*dim+0];  // periodic image of x      
                T xj1 = x[i*dim+1] + pimages[j*dim+1];        
                T xc0 = B2C[0]*xj0 + B2C[2]*xj1;        // map it to the unit square      
                T xc1 = B2C[1]*xj0 + B2C[3]*xj1;        
                /// check if the mapped point is inside the bounding box
                if ((wc[0] <= xc0) && (xc0 <= wc[4]) &&  (wc[1] <= xc1) && (xc1 <= wc[5])) {                    
                    x[dim*(q+k)+0] = xj0; // add the periodic image as a ghost atom
                    x[dim*(q+k)+1] = xj1; //
                    atomtype[q+k] = atomtype[i];
                    alist[q+k] = i;       // add atom i to the list
                    k += 1;            // increase the number of ghost atoms for atom i by 1    
                }
            }
        }    
    }
    else if (dim==3){
        for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
            alist[i] = i; // add atom i to the list
            int q = n + glistnumsum[i]; // offset the starting position by n
            int k = 0;            
            for (int j=1; j<m; j++) { // loop over each periodic image of atom i
                T xj0 = x[i*dim+0] + pimages[j*dim+0];    // periodic image of x          
                T xj1 = x[i*dim+1] + pimages[j*dim+1];        
                T xj2 = x[i*dim+2] + pimages[j*dim+2];        
                T xc0 = B2C[0]*xj0 + B2C[3]*xj1 + B2C[6]*xj2;  // map it to the unit square            
                T xc1 = B2C[1]*xj0 + B2C[4]*xj1 + B2C[7]*xj2;        
                T xc2 = B2C[2]*xj0 + B2C[5]*xj1 + B2C[8]*xj2;   
                /// check if the mapped point is inside the bounding box
                if ((wc[0] <= xc0) && (xc0 <= wc[18]) && (wc[1] <= xc1) && (xc1 <= wc[19]) && (wc[2] <= xc2) && (xc2 <= wc[20])) {
                    x[dim*(q+k)+0] = xj0; // add the periodic image as a ghost atom
                    x[dim*(q+k)+1] = xj1;  
                    x[dim*(q+k)+2] = xj2;
                    atomtype[q+k] = atomtype[i];
                    alist[q+k] = i;      // add atom i to the list
                    k += 1;           // increase the number of ghost atoms by 1              
                }
            }
        }            
    }    
}
template void cpuCreateAtomList(int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void cpuCreateAtomList(int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuCellList(int *clist, int *c2anum, T *xi, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int inum, int gnum, int dim)
{    
    int k, n = inum+gnum;
    if (dim==2) {        
        k = nc[0]*nc[1];   // the total number of cells
        for (int i=0; i<k; i++)
            c2anum[i] = 0; // initialize the number of atoms in every cell to zero
        
        for (int i=0; i<n; i++) { // loop over each atom i
            int j1, j2;
            k = dim*i;
            // position of atom i in the unit square/cube
            T xt0 = B2C[0]*xi[k] + B2C[2]*xi[k+1];        
            T xt1 = B2C[1]*xi[k] + B2C[3]*xi[k+1];     
            // identify a cell containing atom i 
            for (j1=0; j1<nc[0]; j1++)
                if ((eta1[j1] <= xt0) && (xt0<= eta1[j1+1]))
                    break;
            for (j2=0; j2<nc[1]; j2++) 
                if ((eta2[j2] <= xt1) && (xt1<= eta2[j2+1]))
                    break;
            clist[i] = j1 + nc[0]*j2; // link that cell to atom i
            c2anum[clist[i]] += 1;    // increase the number of atoms in that cell by 1 (use atomicAdd on GPU to avoid race condition)            
        }                        
    }        
    else if (dim==3) {                
        k = nc[0]*nc[1]*nc[2]; // the total number of cells
        for (int i=0; i<k; i++)
            c2anum[i] = 0;    // initialize the number of atoms in every cell to zero
        
        for (int i=0; i<n; i++) { // loop over each atom i
            int j1, j2, j3;
            k = dim*i;
            // position of atom i in the unit square/cube
            T xt0 = B2C[0]*xi[k+0] + B2C[3]*xi[k+1] + B2C[6]*xi[k+2];        
            T xt1 = B2C[1]*xi[k+0] + B2C[4]*xi[k+1] + B2C[7]*xi[k+2];                                    
            T xt2 = B2C[2]*xi[k+0] + B2C[5]*xi[k+1] + B2C[8]*xi[k+2];        
            // identify a cell containing atom i 
            for (j1=0; j1<nc[0]; j1++)
                if ((eta1[j1] <= xt0) && (xt0<= eta1[j1+1]))
                    break;
            for (j2=0; j2<nc[1]; j2++) 
                if ((eta2[j2] <= xt1) && (xt1<= eta2[j2+1]))
                    break;
            for (j3=0; j3<nc[2]; j3++) 
                if ((eta3[j3] <= xt2) && (xt2<= eta3[j3+1]))
                    break;
            clist[i] = j1 + nc[0]*j2 + nc[0]*nc[1]*j3; // link that cell to atom i
            c2anum[clist[i]] += 1; // use atomicAdd on GPU to avoid race condition
        }                        
    }        
}
template void cpuCellList(int*, int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void cpuCellList(int*, int*, float*, float*, float*, float*, float*, int*, int, int, int);

void cpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *nc, int inum, int gnum, int dim)
{   
    // number of cells
    int ncell = (dim==2) ? nc[0]*nc[1] : nc[0]*nc[1]*nc[2];
    
    // number of atoms
    int natom = inum+gnum;        
    
    // a list contains the starting positions of the first atom of a cell
    cpuCumsum(c2anumsum, c2anum, ncell+1); 
    
    // initialize the list containting the number of atoms in each cell
    for (int k=0; k<ncell; k++)
        c2anum[k] = 0;
    
    for (int i=0; i<natom; i++) { // for each atom i           
        int k = clist[i];         // cell k contains atom i           
        int i1 = c2anumsum[k];    // starting position for cell k
        int i2 = c2anum[k];       // position of the atom i in cell k
        c2alist[i1+i2] = i;       // add atom i to the list of cell k
        c2anum[k] += 1;           // increase the number of atoms in cell k by 1 (using atomicAdd on GPU to avoid race condition)
    }                            
}

// void cpuCellList3(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *nc, int inum, int gnum, int dim)
// {            
//     // number of cells
//     int ncell = (dim==2) ? nc[0]*nc[1] : nc[0]*nc[1]*nc[2];
//     
//     // number of atoms
//     int natom = inum+gnum;        
// 
//     for (int k=0; k<ncell; k++)      // cell k
//     {
//         int i2 = 0;
//         for (int i=0; i<natom; i++)  // atom i
//             if (clist[i]==k)         // if cell k contains atom i
//                i2 += 1;              // increase the number of atoms in cell k by 1
//         c2anum[k] = i2;              // the number of atoms in cell k
//     }    
//     
//     // a list contains the starting positions of the first atom of a cell
//     cpuCumsum(c2anumsum, c2anum, ncell+1); 
//     
//     for (int k=0; k<ncell; k++)       // cell k
//     {
//         int i1 = c2anumsum[k];        // starting position for cell k
//         int i2 = 0;                   // staring position of the first atom in cell k
//         for (int i=0; i<natom; i++)   // atom i
//             if (clist[i]==k) {        // if cell k contains atom i
//                 c2alist[i1+i2] = i;   // add atom i to the list of cell k
//                 i2 += 1;              // increase the number of atoms in cell k by 1
//             }
//     }    
// }
// 

template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *alist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
{
    T onet = (T) 1.0;
    if (dim==2) {
        T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
        T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
        T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
        T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                          
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = alist[ii];
            verletnum[i] = 0;
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i            
            int j = clist[i];         // cell j contains atom i           
            int j1 = j%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    int k = k1 + nc[0]*k2; // cell k
                    int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];  // starting position of the first atom in cell k
                    for (int l=0; l<m ; l++) {
                        j = c2alist[s+l];  // atom j
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi
                        // distance between atom i and atom j 
                        T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                        if (rij <= onet)
                            verletnum[i] += 1;
                    }
                }
            }                
        }        
    }
    else if (dim==3) {
        T A00 = ellipsoid[0]; // ellipsoid for element t   
        T A10 = ellipsoid[1]; // ellipsoid for element t   
        T A20 = ellipsoid[2]; // ellipsoid for element t   
        T A01 = ellipsoid[3]; // ellipsoid for element t   
        T A11 = ellipsoid[4]; // ellipsoid for element t   
        T A21 = ellipsoid[5]; // ellipsoid for element t   
        T A02 = ellipsoid[6]; // ellipsoid for element t   
        T A12 = ellipsoid[7]; // ellipsoid for element t   
        T A22 = ellipsoid[8]; // ellipsoid for element t                                                                   
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = alist[ii];
            verletnum[i] = 0;
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int j = clist[i];         // cell j contains atom i       
            int n = j%(nc[0]*nc[1]);
            int j3 = (j-n)/(nc[0]*nc[1]);            
            int j1 = n%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    for (int i3=-1; i3<=1; i3++) {
                        int k3 = j3 + i3;                    
                        int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                        int m = c2anum[k];     // number of atoms in cell k
                        int s = c2anumsum[k];  // starting position of the first atom in cell k
                        for (int l=0; l<m ; l++) {
                            j = c2alist[s+l];  // atom j
                            T xij0 = x[j*dim] - xi0;  // xj - xi
                            T xij1 = x[j*dim+1] - xi1; // xj - xi
                            T xij2 = x[j*dim+2] - xi2; // xj - xi
                            // distance between atom i and atom j 
                            T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                            if (rij <= onet)
                                verletnum[i] += 1;
                        }
                    }
                }
            }                
        }                
    }        
}
template void cpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int);
template void cpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int);

template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum,  
        int *alist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
{
    // a list contains the starting positions of the first j atom 
    cpuCumsum(verletnumsum, verletnum, inum+1); 
    
    T onet = (T) 1.0;
    if (dim==2) {
        T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
        T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
        T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
        T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                  
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = alist[ii];
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int nstart = verletnumsum[i];     // starting 
            int ninc = 0;                // increment
            int j = clist[i];         // cell j contains atom i           
            int j1 = j%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    int k = k1 + nc[0]*k2; // cell k
                    int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];  // starting position of the first atom in cell k
                    for (int l=0; l<m ; l++) { //loop over each atom j in cell k
                        j = c2alist[s+l];  // atom j
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi
                        // distance between atom i and atom j 
                        T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                        if (rij <= onet) {
                            verletlist[nstart + ninc] = j; // add atom j into the list
                            ninc += 1;
                        }
                    }
                }
            }               
            verletnum[i] = ninc;
        }        
    }
    else if (dim==3) {
        T A00 = ellipsoid[0]; // ellipsoid for element t   
        T A10 = ellipsoid[1]; // ellipsoid for element t   
        T A20 = ellipsoid[2]; // ellipsoid for element t   
        T A01 = ellipsoid[3]; // ellipsoid for element t   
        T A11 = ellipsoid[4]; // ellipsoid for element t   
        T A21 = ellipsoid[5]; // ellipsoid for element t   
        T A02 = ellipsoid[6]; // ellipsoid for element t   
        T A12 = ellipsoid[7]; // ellipsoid for element t   
        T A22 = ellipsoid[8]; // ellipsoid for element t                                                                           
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = alist[ii];
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int nstart = verletnumsum[i];     // starting 
            int ninc = 0;                // increment            
            int j = clist[i];         // cell j contains atom i       
            int n = j%(nc[0]*nc[1]);
            int j3 = (j-n)/(nc[0]*nc[1]);            
            int j1 = n%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    for (int i3=-1; i3<=1; i3++) {
                        int k3 = j3 + i3;                    
                        int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                        int m = c2anum[k];     // number of atoms in cell k
                        int s = c2anumsum[k];  // starting position of the first atom in cell k
                        for (int l=0; l<m ; l++) {
                            j = c2alist[s+l];  // atom j
                            T xij0 = x[j*dim] - xi0;  // xj - xi
                            T xij1 = x[j*dim+1] - xi1; // xj - xi
                            T xij2 = x[j*dim+2] - xi2; // xj - xi
                            // distance between atom i and atom j 
                            T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                            if (rij <= onet) {
                                verletlist[nstart + ninc] = j; // add atom j into the list
                                ninc += 1;
                            }                            
                        }
                    }
                }
            }                
            verletnum[i] = ninc;
        }        
    }        
}
template void cpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
template void cpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);

template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, 
        int *alist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
{
    if (dim==2) {    
        T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
        T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
        T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
        T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                              
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j              
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                if (rij <= 1.0) // if atom j is inside the ellipsoid 
                    ninc += 1;  // increase the number of neighbors by 1               
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        T A00 = ellipsoid[0]; // ellipsoid for element t   
        T A10 = ellipsoid[1]; // ellipsoid for element t   
        T A20 = ellipsoid[2]; // ellipsoid for element t   
        T A01 = ellipsoid[3]; // ellipsoid for element t   
        T A11 = ellipsoid[4]; // ellipsoid for element t   
        T A21 = ellipsoid[5]; // ellipsoid for element t   
        T A02 = ellipsoid[6]; // ellipsoid for element t   
        T A12 = ellipsoid[7]; // ellipsoid for element t   
        T A22 = ellipsoid[8]; // ellipsoid for element t                                                                                   
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j                
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                if (rij <= 1.0) // if atom j is inside the ellipsoid 
                    ninc += 1;  // increase the number of neighbors by 1                               
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int, int);
template void cpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int, int);

template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *alist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
{
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(neighnumsum, neighnum, inum+1); 
        
    if (dim==2) {    
        T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
        T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
        T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
        T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                                          
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j                
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                    neighlist[nstart+ninc] = j;
                    ninc += 1;  // increase the number of neighbors by 1               
                }
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        T A00 = ellipsoid[0]; // ellipsoid for element t   
        T A10 = ellipsoid[1]; // ellipsoid for element t   
        T A20 = ellipsoid[2]; // ellipsoid for element t   
        T A01 = ellipsoid[3]; // ellipsoid for element t   
        T A11 = ellipsoid[4]; // ellipsoid for element t   
        T A21 = ellipsoid[5]; // ellipsoid for element t   
        T A02 = ellipsoid[6]; // ellipsoid for element t   
        T A12 = ellipsoid[7]; // ellipsoid for element t   
        T A22 = ellipsoid[8]; // ellipsoid for element t                                                                                           
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j              
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                    neighlist[nstart+ninc] = j;
                    ninc += 1;  // increase the number of neighbors by 1               
                }
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int);
template void cpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int);

template <typename T> void cpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, 
        int *alist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
{    
    if (dim==2) {    
        T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
        T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
        T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
        T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                                                  
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j        
                if (j >= i) {
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                    if (rij <= 1.0) // if atom j is inside the ellipsoid 
                        ninc += 1;  // increase the number of neighbors by 1            
                }
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        T A00 = ellipsoid[0]; // ellipsoid for element t   
        T A10 = ellipsoid[1]; // ellipsoid for element t   
        T A20 = ellipsoid[2]; // ellipsoid for element t   
        T A01 = ellipsoid[3]; // ellipsoid for element t   
        T A11 = ellipsoid[4]; // ellipsoid for element t   
        T A21 = ellipsoid[5]; // ellipsoid for element t   
        T A02 = ellipsoid[6]; // ellipsoid for element t   
        T A12 = ellipsoid[7]; // ellipsoid for element t   
        T A22 = ellipsoid[8]; // ellipsoid for element t                                                                                                   
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j                
                if (j >= i) {
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    T xij2 = x[j*dim+2] - xi2; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                    if (rij <= 1.0) // if atom j is inside the ellipsoid 
                        ninc += 1;  // increase the number of neighbors by 1        
                }
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuHalfNeighNum(int*, double*, double*, int*, int*, int*, int*, int, int);
template void cpuHalfNeighNum(int*, float*, float*, int*, int*, int*, int*, int, int);

template <typename T> void cpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *alist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
{
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(neighnumsum, neighnum, inum+1); 
        
    if (dim==2) {    
        T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
        T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
        T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
        T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                                                          
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i            
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j     
                if (j >= i) {
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                    if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                        neighlist[nstart+ninc] = j;
                        ninc += 1;  // increase the number of neighbors by 1               
                    }
                }
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        T A00 = ellipsoid[0]; // ellipsoid for element t   
        T A10 = ellipsoid[1]; // ellipsoid for element t   
        T A20 = ellipsoid[2]; // ellipsoid for element t   
        T A01 = ellipsoid[3]; // ellipsoid for element t   
        T A11 = ellipsoid[4]; // ellipsoid for element t   
        T A21 = ellipsoid[5]; // ellipsoid for element t   
        T A02 = ellipsoid[6]; // ellipsoid for element t   
        T A12 = ellipsoid[7]; // ellipsoid for element t   
        T A22 = ellipsoid[8]; // ellipsoid for element t                                                                                                           
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = alist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i            
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j           
                if (j >= i) {
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    T xij2 = x[j*dim+2] - xi2; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                    if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                        neighlist[nstart+ninc] = j;
                        ninc += 1;  // increase the number of neighbors by 1               
                    }
                }
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuHalfNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int);
template void cpuHalfNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int);

template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, int inum, int dim)
{    

    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = m;              
    }
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    if (dim==2) {    
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];       // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i            
            int itype = atomtype[i];
            int istart = neighnumsum[i];               
            int m = anum[ii];        // number of neighbors around i             
            int nstart = anumsum[ii];   
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = neighlist[istart+l];  // atom j  
                int k = nstart + l;                   
                xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
                xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi     
                ai[k]        = i;
                aj[k]        = j; // should be alist[j];         
                ti[k]        = itype;       
                tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
            }
        }
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int itype = atomtype[i];
            int istart = neighnumsum[i];   
            int m = anum[ii];        // number of neighbors around i             
            int nstart = anumsum[ii];   
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = neighlist[istart+l];  // atom j    
                int k = nstart+l;                 
                xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
                xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi                                                
                xij[k*dim+2] = x[j*dim+2] - xi2; // xj - xi      
                ai[k]        = i;
                aj[k]        = j;                
                ti[k]        = itype;       
                tj[k]        = atomtype[j];
            }                        
        }
    }    
}
template void cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
template void cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
   
template <typename T> int cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
       int typei, int inum, int dim)
{    
    // form ilist from olist
    int tnum = 0;
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];           // atom i
        if (atomtype[i] == typei)
            tlist[tnum++] = i;
    }
    
    cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, tlist, 
            neighlist, neighnum, neighnumsum, atomtype, tnum, dim);     
    
    return tnum;
}
template int cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template int cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, int typej, int inum, int dim)
{    
    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int istart = neighnumsum[i];  
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = neighlist[istart+l];  // atom j              
            if (atomtype[j] == typej) 
                anum[ii] += 1;
        }                
    }
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    if (dim==2) {    
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];       // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int itype = atomtype[i];
            int istart = neighnumsum[i];   
            int m = neighnum[i];     // number of neighbors around i                     
            int nstart = anumsum[ii];      
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = neighlist[istart+l];  // atom j  
                int jtype = atomtype[j];
                if (jtype == typej) {                    
                    int k = nstart + anum[ii];                   
                    xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
                    xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi     
                    ai[k]        = i;
                    aj[k]        = j;          
                    ti[k]        = itype;       
                    tj[k]        = jtype;    
                }
            }
        }
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int itype = atomtype[i];
            int istart = neighnumsum[i];   
            int m = neighnum[i];     // number of neighbors around i             
            int nstart = anumsum[ii];   
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = neighlist[istart+l];  // atom j    
                int jtype = atomtype[j];
                if (jtype == typej) {                    
                    int k = nstart + anum[ii];                 
                    xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
                    xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi                                                
                    xij[k*dim+2] = x[j*dim+2] - xi2; // xj - xi      
                    ai[k]        = i;
                    aj[k]        = j;                
                    ti[k]        = itype;       
                    tj[k]        = jtype;    
                }
            }            
        }
    }    
}
template void cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);   

template <typename T> int cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
       int typei, int typej, int inum, int dim)
{    
    // form ilist from olist
    int tnum = 0;
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];           // atom i
        if (atomtype[i] == typei)
            tlist[tnum++] = i;
    }
    
    cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, tlist, 
            neighlist, neighnum, neighnumsum, atomtype, typej, tnum, dim);
    
    return tnum;
}
template int cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
template int cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);


// template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
//         int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
// {    
//     if (dim==2) {    
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//             int i = ilist[ii];  // atom i
//             int itype = atomtype[i];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             int m = neighnum[i];     // number of neighbors around i 
//             int nstart = neighnumsum[i];   
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int k = nstart+l;
//                 int j = neighlist[k];  // atom j     
//                 xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
//                 xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi      
//                 ti[k]       = itype;
//                 tj[k]       = atomtype[j];
//             }
//         }
//     }
//     else if (dim==3) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
//             int i = ilist[ii];  // atom i
//             int itype = atomtype[i];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             T xi2 = x[i*dim+2];      // position of atom i
//             int m = neighnum[i];     // number of neighbors around i 
//             int nstart = neighnumsum[i];               
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int k = nstart+l;
//                 int j = neighlist[k];  // atom j     
//                 xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
//                 xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi                                                
//                 xij[k*dim+2] = x[j*dim+2] - xi2; // xj - xi       
//                 ti[k]       = itype;
//                 tj[k]       = atomtype[j];
//             }            
//         }
//     }
// }
// template void cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// template void cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

// template <typename T> void cpuGetNeighPairs(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
//         int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
// {    
//     if (dim==2) {    
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//             int i = ilist[ii];  // atom i
//             int itype = atomtype[i];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i            
//             int m = neighnum[i];     // number of neighbors around i 
//             int nstart = neighnumsum[i];   
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int k = nstart+l;
//                 int j = neighlist[k];  // atom j     
//                 xi[k*dim]   = xi0;        // position of atom i
//                 xi[k*dim+1] = xi1;       // position of atom i                
//                 xj[k*dim]   = x[j*dim];  // xj 
//                 xj[k*dim+1] = x[j*dim+1]; // xj    
//                 ti[k]       = itype;
//                 tj[k]       = atomtype[j];
//             }
//         }
//     }
//     else if (dim==3) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
//             int i = ilist[ii];  // atom i
//             int itype = atomtype[i];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             T xi2 = x[i*dim+2];      // position of atom i
//             int m = neighnum[i];     // number of neighbors around i 
//             int nstart = neighnumsum[i];               
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int k = nstart+l;
//                 int j = neighlist[k];  // atom j     
//                 xi[k*dim]   = xi0;        // position of atom i
//                 xi[k*dim+1] = xi1;       // position of atom i                
//                 xi[k*dim+2] = xi2;       // position of atom i                
//                 xj[k*dim]   = x[j*dim];  // xj - xi
//                 xj[k*dim+1] = x[j*dim+1]; // xj - xi                                                
//                 xj[k*dim+2] = x[j*dim+2]; // xj - xi       
//                 ti[k]       = itype;
//                 tj[k]       = atomtype[j];
//             }            
//         }
//     }
// }
// template void cpuGetNeighPairs(double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


// ***************************************************************/
template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim)
{
    int dimsq = dim*dim;
    T onet = (T) 1.0;
    if (dim==2) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = ilist[ii];
            verletnum[i] = 0;
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int ti = atomtype[i];      // element of atom i 
            int j = clist[i];         // cell j contains atom i           
            int j1 = j%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    int k = k1 + nc[0]*k2; // cell k
                    int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];  // starting position of the first atom in cell k
                    for (int l=0; l<m ; l++) {
                        j = c2alist[s+l];  // atom j
                        int tj = atomtype[j];      // element of atom j
                        int tij = dimsq*(ti + tj*ntype);
                        T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                        T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                        T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                        T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                          
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi
                        // distance between atom i and atom j 
                        T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                        if (rij <= onet)
                            verletnum[i] += 1;
                    }
                }
            }                
        }        
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = ilist[ii];
            verletnum[i] = 0;
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int ti = atomtype[i];      // element of atom i 
            int j = clist[i];         // cell j contains atom i       
            int n = j%(nc[0]*nc[1]);
            int j3 = (j-n)/(nc[0]*nc[1]);            
            int j1 = n%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    for (int i3=-1; i3<=1; i3++) {
                        int k3 = j3 + i3;                    
                        int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                        int m = c2anum[k];     // number of atoms in cell k
                        int s = c2anumsum[k];  // starting position of the first atom in cell k
                        for (int l=0; l<m ; l++) {
                            j = c2alist[s+l];  // atom j
                            int tj = atomtype[j];      // element of atom j
                            int tij = dimsq*(ti + tj*ntype);
                            T A00 = ellipsoid[tij]; // ellipsoid for element t   
                            T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                            T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                            T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                            T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                            T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                            T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                            T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                            T A22 = ellipsoid[tij+8]; // ellipsoid for element t                               
                            T xij0 = x[j*dim] - xi0;  // xj - xi
                            T xij1 = x[j*dim+1] - xi1; // xj - xi
                            T xij2 = x[j*dim+2] - xi2; // xj - xi
                            // distance between atom i and atom j 
                            T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                            if (rij <= onet)
                                verletnum[i] += 1;
                        }
                    }
                }
            }                
        }                
    }        
}
template void cpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
     int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim)
{
    // a list contains the starting positions of the first j atom 
    cpuCumsum(verletnumsum, verletnum, inum+1); 
    
    int dimsq = dim*dim;
    T onet = (T) 1.0;
    if (dim==2) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = ilist[ii];
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int ti = atomtype[i];      // element of atom i             
            int nstart = verletnumsum[i];     // starting 
            int ninc = 0;                // increment
            int j = clist[i];         // cell j contains atom i           
            int j1 = j%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    int k = k1 + nc[0]*k2; // cell k
                    int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];  // starting position of the first atom in cell k
                    for (int l=0; l<m ; l++) { //loop over each atom j in cell k
                        j = c2alist[s+l];  // atom j
                        int tj = atomtype[j];      // element of atom j
                        int tij = dimsq*(ti + tj*ntype);
                        T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                        T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                        T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                        T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                  
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi
                        // distance between atom i and atom j 
                        T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                        if (rij <= onet) {
                            verletlist[nstart + ninc] = j; // add atom j into the list
                            ninc += 1;
                        }
                    }
                }
            }               
            verletnum[i] = ninc;
        }        
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
            int i = ilist[ii];
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int ti = atomtype[i];      // element of atom i  
            int nstart = verletnumsum[i];     // starting 
            int ninc = 0;                // increment            
            int j = clist[i];         // cell j contains atom i       
            int n = j%(nc[0]*nc[1]);
            int j3 = (j-n)/(nc[0]*nc[1]);            
            int j1 = n%nc[0];
            int j2 = (j-j1)/nc[0];
            for (int i1=-1; i1<=1; i1++) {
                int k1 = j1 + i1;
                for (int i2=-1; i2<=1; i2++) {
                    int k2 = j2 + i2;
                    for (int i3=-1; i3<=1; i3++) {
                        int k3 = j3 + i3;                    
                        int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                        int m = c2anum[k];     // number of atoms in cell k
                        int s = c2anumsum[k];  // starting position of the first atom in cell k
                        for (int l=0; l<m ; l++) {
                            j = c2alist[s+l];  // atom j
                            int tj = atomtype[j];      // element of atom j
                            int tij = dimsq*(ti + tj*ntype);
                            T A00 = ellipsoid[tij]; // ellipsoid for element t   
                            T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                            T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                            T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                            T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                            T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                            T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                            T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                            T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                           
                            T xij0 = x[j*dim] - xi0;  // xj - xi
                            T xij1 = x[j*dim+1] - xi1; // xj - xi
                            T xij2 = x[j*dim+2] - xi2; // xj - xi
                            // distance between atom i and atom j 
                            T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                            if (rij <= onet) {
                                verletlist[nstart + ninc] = j; // add atom j into the list
                                ninc += 1;
                            }                            
                        }
                    }
                }
            }                
            verletnum[i] = ninc;
        }        
    }        
}
template void cpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
{
    int dimsq = dim*dim;
    
    if (dim==2) {    
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int ti = atomtype[i];      // element of atom i           
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j              
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                  
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                if (rij <= 1.0) // if atom j is inside the ellipsoid 
                    ninc += 1;  // increase the number of neighbors by 1               
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int ti = atomtype[i];      // element of atom i           
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j                
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for element t   
                T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                           
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                if (rij <= 1.0) // if atom j is inside the ellipsoid 
                    ninc += 1;  // increase the number of neighbors by 1                               
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void cpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
{
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(neighnumsum, neighnum, inum+1); 
    
    int dimsq = dim*dim;
    
    if (dim==2) {    
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int ti = atomtype[i];      // element of atom i
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j                
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                                  
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                    neighlist[nstart+ninc] = j;
                    ninc += 1;  // increase the number of neighbors by 1               
                }
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int ti = atomtype[i];      // element of atom i
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j              
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for element t   
                T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                                           
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                    neighlist[nstart+ninc] = j;
                    ninc += 1;  // increase the number of neighbors by 1               
                }
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
{
    int dimsq = dim*dim;
    
    if (dim==2) {    
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int ti = atomtype[i];      // element of atom i           
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j        
                if (j >= i) {
                    int tj = atomtype[j];      // element of atom j
                    int tij = dimsq*(ti + tj*ntype);
                    T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                    T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                    T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                    T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                  
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                    if (rij <= 1.0) // if atom j is inside the ellipsoid 
                        ninc += 1;  // increase the number of neighbors by 1            
                }
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int ti = atomtype[i];      // element of atom i           
            int m = verletnum[i]; // number of atoms around i 
            int start = verletnumsum[i];     // starting 
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[start+l];  // atom j                
                if (j >= i) {
                    int tj = atomtype[j];      // element of atom j
                    int tij = dimsq*(ti + tj*ntype);
                    T A00 = ellipsoid[tij]; // ellipsoid for element t   
                    T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                    T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                    T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                    T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                    T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                    T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                    T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                    T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                           
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    T xij2 = x[j*dim+2] - xi2; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                    if (rij <= 1.0) // if atom j is inside the ellipsoid 
                        ninc += 1;  // increase the number of neighbors by 1        
                }
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuHalfNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void cpuHalfNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
{
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(neighnumsum, neighnum, inum+1); 
    
    int dimsq = dim*dim;
    
    if (dim==2) {    
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            int ti = atomtype[i];      // element of atom i
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j     
                if (j >= i) {
                    int tj = atomtype[j];      // element of atom j
                    int tij = dimsq*(ti + tj*ntype);
                    T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                    T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                    T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                    T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                                  
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                    if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                        neighlist[nstart+ninc] = j;
                        ninc += 1;  // increase the number of neighbors by 1               
                    }
                }
            }
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
    else if (dim==3) {
        for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
            int i = ilist[ii];  // atom i
            T xi0 = x[i*dim];        // position of atom i
            T xi1 = x[i*dim+1];      // position of atom i
            T xi2 = x[i*dim+2];      // position of atom i
            int ti = atomtype[i];      // element of atom i
            int m = verletnum[i]; // number of atoms around i 
            int jstart = verletnumsum[i];     // starting 
            int nstart = neighnumsum[i];   
            int ninc = 0;                // increment 
            for (int l=0; l<m ; l++) {   // loop over each atom around atom i
                int j = verletlist[jstart+l];  // atom j           
                if (j >= i) {
                    int tj = atomtype[j];      // element of atom j
                    int tij = dimsq*(ti + tj*ntype);
                    T A00 = ellipsoid[tij]; // ellipsoid for element t   
                    T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                    T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                    T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                    T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                    T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                    T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                    T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                    T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                                           
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    T xij2 = x[j*dim+2] - xi2; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                    if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                        neighlist[nstart+ninc] = j;
                        ninc += 1;  // increase the number of neighbors by 1               
                    }
                }
            }            
            neighnum[i] = ninc; // number of neighbors of atom i
        }
    }
}
template void cpuHalfNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuHalfNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


// template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, 
//         int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {
//     int dimsq = dim*dim;
//     T onet = (T) 1.0;
//     if (dim==2) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//             int i = ilist[ii];
//             verletnum[i] = 0;
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//             int j = clist[i];         // cell j contains atom i           
//             int j1 = j%nc[0];
//             int j2 = (j-j1)/nc[0];
//             for (int i1=-1; i1<=1; i1++) {
//                 int k1 = j1 + i1;
//                 for (int i2=-1; i2<=1; i2++) {
//                     int k2 = j2 + i2;
//                     int k = k1 + nc[0]*k2; // cell k
//                     int m = c2anum[k];     // number of atoms in cell k
//                     int s = c2anumsum[k];  // starting position of the first atom in cell k
//                     for (int l=0; l<m ; l++) {
//                         j = c2alist[s+l];  // atom j
//                         T xij0 = x[j*dim] - xi0;  // xj - xi
//                         T xij1 = x[j*dim+1] - xi1; // xj - xi
//                         // distance between atom i and atom j 
//                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                         if (rij <= onet)
//                             verletnum[i] += 1;
//                     }
//                 }
//             }                
//         }        
//     }
//     else if (dim==3) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//             int i = ilist[ii];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             T xi2 = x[i*dim+2];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
//             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
//             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
//             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
//             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
//             int j = clist[i];         // cell j contains atom i       
//             int n = j%(nc[0]*nc[1]);
//             int j3 = (j-n)/(nc[0]*nc[1]);            
//             int j1 = n%nc[0];
//             int j2 = (j-j1)/nc[0];
//             for (int i1=-1; i1<=1; i1++) {
//                 int k1 = j1 + i1;
//                 for (int i2=-1; i2<=1; i2++) {
//                     int k2 = j2 + i2;
//                     for (int i3=-1; i3<=1; i3++) {
//                         int k3 = j3 + i3;                    
//                         int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
//                         int m = c2anum[k];     // number of atoms in cell k
//                         int s = c2anumsum[k];  // starting position of the first atom in cell k
//                         for (int l=0; l<m ; l++) {
//                             j = c2alist[s+l];  // atom j
//                             T xij0 = x[j*dim] - xi0;  // xj - xi
//                             T xij1 = x[j*dim+1] - xi1; // xj - xi
//                             T xij2 = x[j*dim+2] - xi2; // xj - xi
//                             // distance between atom i and atom j 
//                             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//                             if (rij <= onet)
//                                 verletnum[i] += 1;
//                         }
//                     }
//                 }
//             }                
//         }                
//     }        
// }
// template void cpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int);
// template void cpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int);
// 
// template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
//      int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {
//     // a list contains the starting positions of the first j atom 
//     cpuCumsum(verletnumsum, verletnum, inum+1); 
//     
//     int dimsq = dim*dim;
//     T onet = (T) 1.0;
//     if (dim==2) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//             int i = ilist[ii];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//             int nstart = verletnumsum[i];     // starting 
//             int ninc = 0;                // increment
//             int j = clist[i];         // cell j contains atom i           
//             int j1 = j%nc[0];
//             int j2 = (j-j1)/nc[0];
//             for (int i1=-1; i1<=1; i1++) {
//                 int k1 = j1 + i1;
//                 for (int i2=-1; i2<=1; i2++) {
//                     int k2 = j2 + i2;
//                     int k = k1 + nc[0]*k2; // cell k
//                     int m = c2anum[k];     // number of atoms in cell k
//                     int s = c2anumsum[k];  // starting position of the first atom in cell k
//                     for (int l=0; l<m ; l++) { //loop over each atom j in cell k
//                         j = c2alist[s+l];  // atom j
//                         T xij0 = x[j*dim] - xi0;  // xj - xi
//                         T xij1 = x[j*dim+1] - xi1; // xj - xi
//                         // distance between atom i and atom j 
//                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                         if (rij <= onet) {
//                             verletlist[nstart + ninc] = j; // add atom j into the list
//                             ninc += 1;
//                         }
//                     }
//                 }
//             }               
//             verletnum[i] = ninc;
//         }        
//     }
//     else if (dim==3) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//             int i = ilist[ii];
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             T xi2 = x[i*dim+2];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
//             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
//             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
//             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
//             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
//             int nstart = verletnumsum[i];     // starting 
//             int ninc = 0;                // increment            
//             int j = clist[i];         // cell j contains atom i       
//             int n = j%(nc[0]*nc[1]);
//             int j3 = (j-n)/(nc[0]*nc[1]);            
//             int j1 = n%nc[0];
//             int j2 = (j-j1)/nc[0];
//             for (int i1=-1; i1<=1; i1++) {
//                 int k1 = j1 + i1;
//                 for (int i2=-1; i2<=1; i2++) {
//                     int k2 = j2 + i2;
//                     for (int i3=-1; i3<=1; i3++) {
//                         int k3 = j3 + i3;                    
//                         int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
//                         int m = c2anum[k];     // number of atoms in cell k
//                         int s = c2anumsum[k];  // starting position of the first atom in cell k
//                         for (int l=0; l<m ; l++) {
//                             j = c2alist[s+l];  // atom j
//                             T xij0 = x[j*dim] - xi0;  // xj - xi
//                             T xij1 = x[j*dim+1] - xi1; // xj - xi
//                             T xij2 = x[j*dim+2] - xi2; // xj - xi
//                             // distance between atom i and atom j 
//                             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//                             if (rij <= onet) {
//                                 verletlist[nstart + ninc] = j; // add atom j into the list
//                                 ninc += 1;
//                             }                            
//                         }
//                     }
//                 }
//             }                
//             verletnum[i] = ninc;
//         }        
//     }        
// }
// template void cpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
// template void cpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
// 
// template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
// {
//     int dimsq = dim*dim;
//     
//     if (dim==2) {    
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//             int i = ilist[ii];  // atom i
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t      
//             int m = verletnum[i]; // number of atoms around i 
//             int start = verletnumsum[i];     // starting 
//             int ninc = 0;                // increment 
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int j = verletlist[start+l];  // atom j                
//                 T xij0 = x[j*dim] - xi0;  // xj - xi
//                 T xij1 = x[j*dim+1] - xi1; // xj - xi
//                 // distance between atom i and atom j 
//                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                 if (rij <= 1.0) // if atom j is inside the ellipsoid 
//                     ninc += 1;  // increase the number of neighbors by 1               
//             }
//             neighnum[i] = ninc; // number of neighbors of atom i
//         }
//     }
//     else if (dim==3) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
//             int i = ilist[ii];  // atom i
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             T xi2 = x[i*dim+2];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
//             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
//             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
//             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
//             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
//             int m = verletnum[i]; // number of atoms around i 
//             int start = verletnumsum[i];     // starting 
//             int ninc = 0;                // increment 
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int j = verletlist[start+l];  // atom j                
//                 T xij0 = x[j*dim] - xi0;  // xj - xi
//                 T xij1 = x[j*dim+1] - xi1; // xj - xi
//                 T xij2 = x[j*dim+2] - xi2; // xj - xi
//                 // distance between atom i and atom j 
//                 T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//                 if (rij <= 1.0) // if atom j is inside the ellipsoid 
//                     ninc += 1;  // increase the number of neighbors by 1                               
//             }            
//             neighnum[i] = ninc; // number of neighbors of atom i
//         }
//     }
// }
// template void cpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int);
// template void cpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int);
// 
// template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
//         int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
// {
//     // a list contains the starting positions of the first neighbor 
//     cpuCumsum(neighnumsum, neighnum, inum+1); 
//     
//     int dimsq = dim*dim;
//     
//     if (dim==2) {    
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//             int i = ilist[ii];  // atom i
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t      
//             int m = verletnum[i]; // number of atoms around i 
//             int jstart = verletnumsum[i];     // starting 
//             int nstart = neighnumsum[i];   
//             int ninc = 0;                // increment 
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int j = verletlist[jstart+l];  // atom j                
//                 T xij0 = x[j*dim] - xi0;  // xj - xi
//                 T xij1 = x[j*dim+1] - xi1; // xj - xi
//                 // distance between atom i and atom j 
//                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                 if (rij <= 1.0) { // if atom j is inside the ellipsoid 
//                     neighlist[nstart+ninc] = j;
//                     ninc += 1;  // increase the number of neighbors by 1               
//                 }
//             }
//             neighnum[i] = ninc; // number of neighbors of atom i
//         }
//     }
//     else if (dim==3) {
//         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
//             int i = ilist[ii];  // atom i
//             T xi0 = x[i*dim];        // position of atom i
//             T xi1 = x[i*dim+1];      // position of atom i
//             T xi2 = x[i*dim+2];      // position of atom i
//             int t = atomtype[i];      // element of atom i 
//             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
//             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
//             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
//             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
//             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
//             int m = verletnum[i]; // number of atoms around i 
//             int jstart = verletnumsum[i];     // starting 
//             int nstart = neighnumsum[i];   
//             int ninc = 0;                // increment 
//             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//                 int j = verletlist[jstart+l];  // atom j                
//                 T xij0 = x[j*dim] - xi0;  // xj - xi
//                 T xij1 = x[j*dim+1] - xi1; // xj - xi
//                 T xij2 = x[j*dim+2] - xi2; // xj - xi
//                 // distance between atom i and atom j 
//                 T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//                 if (rij <= 1.0) { // if atom j is inside the ellipsoid 
//                     neighlist[nstart+ninc] = j;
//                     ninc += 1;  // increase the number of neighbors by 1               
//                 }
//             }            
//             neighnum[i] = ninc; // number of neighbors of atom i
//         }
//     }
// }
// template void cpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int);
// template void cpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int);



#endif


