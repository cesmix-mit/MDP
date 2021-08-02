/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CPUARRAYOPERATIONS
#define MDP_CPUARRAYOPERATIONS

void cpuTripletnum(int* output, int* input, int length) 
{	
	for (int ii = 1; ii < length; ++ii)	
		output[ii] = (input[ii]-1)*input[ii]/2;       
}
void cpuQuadrupletnum(int* output, int* input, int length) 
{	
	for (int ii = 1; ii < length; ++ii)	
		output[ii] = (input[ii]-2)*(input[ii]-1)*input[ii]/6;       
}

void cpuCumsum(int* output, int* input, int length) 
{
	output[0] = 0; 
	for (int j = 1; j < length; ++j)	
		output[j] = input[j - 1] + output[j - 1];	
}

void cpuArrayFill(int* output, int start, int length) 
{	
	for (int j = 0; j < length; ++j)	
		output[j] = start + j;
}

// template <typename T> T cpuArrayErrorNorm(T *a, T *b, int n)
// {
//     T e = (a[0]-b[0])*(a[0]-b[0]);
//     for (int i=1; i<n; i++)        
//         e += (a[i]-b[i])*(a[i]-b[i]);    
//     return sqrt(e);
// }

template <typename T> T cpuArrayMaxError(T *a, T *b, int n)
{
    T e = fabs(a[0]-b[0]);
    for (int i=1; i<n; i++)
        if (fabs(a[i]-b[i])>e)
            e = fabs(a[i]-b[i]);    
    return e;
}

// template <typename T> T cpuArrayNorm(T *a, int n)
// {
//     T e = a[0]*a[0];
//     for (int i=1; i<n; i++)        
//         e += a[i]*a[i];    
//     return sqrt(e);
// }

template <typename T> T cpuArrayMin(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

template <typename T> T cpuArrayMax(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

template <typename T> T cpuArraySum(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    return b;
}

template <typename T> T cpuArraySquareSum(T *a, int n)
{
    T b = a[0]*a[0];
    for (int i=1; i<n; i++)
            b = b + a[i]*a[i];    
    return b;
}

template <typename T> void cpuArraySumEveryColumn(T *b, T *a, int m, int n)
{
    // m: number of columns
    for (int j=0; j<m; j++) {
        int k = n*j;
        b[j] = a[0 + k];
        for (int i=1; i<n; i++)
            b[j] += a[i + k];    
    }
}

template <typename T> void cpuArrayMaxEveryColumn(T *b, T *a, int m, int n)
{
    for (int j=0; j<m; j++) {
        int k = n*j;
        b[j] = a[0 + k];
        for (int i=1; i<n; i++)
            b[j] = (b[j] > a[i+k]) ? b[j] : a[i+k];    
    }
}

template <typename T> void cpuArrayMinEveryColumn(T *b, T *a, int m, int n)
{
    for (int j=0; j<m; j++) {
        int k = n*j;
        b[j] = a[0 + k];
        for (int i=1; i<n; i++)
            b[j] = (b[j] > a[i+k]) ? a[i+k] : b[j];    
    }
}

template <typename T> void cpuArraySumEveryRow(T *b, T *a, int m, int n)
{
    // m: number of rows
    for (int j=0; j<m; j++) {        
        b[j] = a[j + m*0];
        for (int i=1; i<n; i++)
            b[j] += a[j + m*i];    
    }
}

template <typename T> void cpuArrayMaxEveryRow(T *b, T *a, int m, int n)
{
    for (int j=0; j<m; j++) {        
        b[j] = a[j + m*0];
        for (int i=1; i<n; i++)
            b[j] = (b[j] >  a[j + m*i]) ? b[j] : a[j + m*i];    
    }
}

template <typename T> void cpuArrayMinEveryRow(T *b, T *a, int m, int n)
{
    for (int j=0; j<m; j++) {        
        b[j] = a[j + m*0];
        for (int i=1; i<n; i++)
            b[j] = (b[j] >  a[j + m*i]) ? a[j + m*i] : b[j];    
    }
}

template <typename T> T cpuArrayMean(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    b = b/((T) n);
    return b;
}


template <typename T> void cpuArrayTranspose(T *A, T *B, int m, int n)
{  
    int k = m*n;
    for (int ii=0; ii<k; ii++) {
        int i = ii%m;        
        int j = (ii-i)/m;                
        A[j+n*i] = B[ii]; // ii = i + m*j
    }
}

template <typename T> void cpuGetArrayAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{  
    for (int ii=0; ii<n; ii++) {
        int i = colind[ii];   
        for (int j=0; j<m; j++)
            A[ii*m+j] = B[i*m+j]; 
    }
}

template <typename T> void cpuArrayTransposeAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{  
    for (int ii=0; ii<n; ii++) {
        int i = colind[ii];   
        for (int j=0; j<m; j++)
            A[ii+j*n] = B[i*m+j];        
    }
}

template <typename T> void cpuGetArrayAtRowIndex(T *A, T *B, int *rowind, int m, int n)
{      
    for (int ii=0; ii<m; ii++) {
        int i = rowind[ii];   
        for (int j=0; j<n; j++)
            A[j*m+ii] = B[j*m+i];        
    }
}

template <typename T> void cpuArrayTransposeAtRowIndex(T *A, T *B, int *rowind, int m, int n)
{  
    for (int ii=0; ii<m; ii++) {
        int i = rowind[ii];   
        for (int j=0; j<n; j++)
            A[ii*n+j] = B[j*m+i];
    }
}

template <typename T> void cpuGetArrayAtIndex(T *y, T *x, int *ind, int n)
{    
    //#pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[i] = x[ind[i]];
}

template <typename T> void cpuPutArrayAtIndex(T *y, T *x, int *ind, int n)
{    
    //#pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] = x[i];
}

template <typename T> void cpuArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{    
    //#pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] += x[i];
}

template <typename T> void cpuArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{    
    //#pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] -= x[i];
}

template <typename T> void cpuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{    
    //#pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] = y[ind[i]] + a*x[i];
}

template <typename T> void cpuArraySetValue(T *y, T a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a;        
}

template <typename T> void cpuArrayMultiplyScalar(T *y, T a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a*y[i];        
}

template <typename T> void cpuArrayAddScalar(T *y, T a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = y[i] + a;        
}

template <typename T> void cpuArrayCopy(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = x[i];        
}

template <typename T> void cpuArrayMinus(T *y, T *x, int n)
{    
    //#pragma omp parallel for    
    for (int i=0; i<n; i++) 
        y[i] = -x[i];        
}

template <typename T> void cpuArrayAbs(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = fabs(x[i]);        
}

template <typename T> void cpuArraySqrt(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = sqrt(x[i]);        
}

template <typename T> void cpuArraySin(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = sin(x[i]);        
}

template <typename T> void cpuArrayCos(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = cos(x[i]);        
}

template <typename T> void cpuArrayTan(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = tan(x[i]);        
}

template <typename T> void cpuArrayAsin(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = asin(x[i]);        
}

template <typename T> void cpuArrayAcos(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = acos(x[i]);        
}

template <typename T> void cpuArrayAtan(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = atan(x[i]);        
}

template <typename T> void cpuArraySinh(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = sinh(x[i]);        
}

template <typename T> void cpuArrayCosh(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = cosh(x[i]);        
}

template <typename T> void cpuArrayTanh(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = tanh(x[i]);        
}

template <typename T> void cpuArrayAsinh(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = asinh(x[i]);        
}

template <typename T> void cpuArrayAcosh(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = acosh(x[i]);        
}

template <typename T> void cpuArrayAtanh(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = atanh(x[i]);        
}

template <typename T> void cpuArrayExp(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = exp(x[i]);        
}

template <typename T> void cpuArrayLog(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = log(x[i]);        
}

template <typename T> void cpuArrayCeil(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = ceil(x[i]);        
}

template <typename T> void cpuArrayFloor(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = floor(x[i]);        
}

template <typename T> void cpuArrayErf(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = erf(x[i]);        
}

template <typename T> void cpuArrayErfc(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = erfc(x[i]);        
}

template <typename T> void cpuArraySquare(T *y, T *x, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = x[i]*x[i];        
}

template <typename T> void cpuArrayPower(T *y, T *x, int p, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) {
        y[i] = x[i];    
        for (int j=1; j<p; j++)
            y[i] = y[i]*x[i];
    }
}

template <typename T> void cpuArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    //#pragma omp parallel for
    for (int i=0; i<n; i++)    
        C[i+n*i] = a*C[i+n*i];                
}

template <typename T> void cpuArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{            
    //#pragma omp parallel for
    for (int i=0; i<n; i++)    
        C[i+n*i] += a*x[i];                
}

template <typename T> void cpuArrayRowAverage(T *y, T *x, int m, int n)
{    
    //#pragma omp parallel 
    for (int j=0; j<n; j++)
    {        
        T avg = 0;
        int i;
        for (i=0; i<m; i++)
            avg = avg + x[i + m*j];
        avg = avg/((T) m);
        for (i=0; i<m; i++)
            y[i + m*j] = avg;         
    }        
}

template <typename T> void cpuArrayRowSum(T *y, T *x, int m, int n)
{    
    //#pragma omp parallel 
    for (int j=0; j<n; j++)
    {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += x[i + m*j]; // y[j] += x[j + n*i]
    }        
}

template <typename T> void cpuArrayRowSquareSum(T *y, T *x, int m, int n)
{    
    //#pragma omp parallel 
    for (int j=0; j<n; j++)
    {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += x[i + m*j]*x[i + m*j];
    }        
}

template <typename T> void cpuArrayDistSquareSum(T *y, T *x1, T *x2, int m, int n)
{    
    //#pragma omp parallel 
    for (int j=0; j<n; j++)
    {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += (x2[i + m*j]-x1[i + m*j])*(x2[i + m*j]-x1[i + m*j]);
    }        
}

template <typename T> void cpuArrayRowkAXPB(T *y, T *x, T a, T b, int m, int n, int k)
{    
    //#pragma omp parallel 
    for (int j=0; j<n; j++)    
        y[k + m*j] = a*x[k + m*j] + b;    
}

template <typename T> void cpuArrayAXPB(T *y, T *x, T a, T b, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a*x[i]+b;        
}

template <typename T> void cpuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        z[i] = a*x[i] + b*y[i];        
}

template <typename T> void cpuArrayAXY(T *s, T *x, T *y, T a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i];        
}

template <typename T> void cpuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{    
    //#pragma omp parallel for
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i]*z[i];        
}

template <typename T> void cpuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{    
    //#pragma omp parallel for    
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i] + b*z[i];        
}

template <typename T> void cpuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{    
    //#pragma omp parallel for        
    for (int i=0; i<n; i++) 
        s[i] = a*x[i] + b*y[i] + c*z[i];        
}

template <typename T> void cpuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    //#pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        un[idx] = u[i+I*j+I*J*k];
    }            
}

template <typename T> void cpuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    //#pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+I*J*k] = un[idx];        
    }            
}

template <typename T> void cpuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    //#pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        C[i+I*j+I*J*s] = 0.0;
        for (int k=0; k<K; k++)
            C[i+I*j+I*J*s] += A[i+I*k+I*K*s]*B[k+K*j+K*J*s];
    }            
}

template <typename T> void cpuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    //#pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        for (int k=0; k<K; k++)
            C[i+I*j+I*J*s] += A[i+I*k+I*K*s]*B[k+K*j+K*J*s];
    }            
}

template <typename T> void cpuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{    
    //#pragma omp parallel for
    for (int i=0; i<nent; i++)         
    {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++)
            ucg[i] += udg[cgent2dgent[rowent2elem[i]+k]]; 
        ucg[i] = ucg[i]/((T) nelem);
    }
}

template <typename T> void cpuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{    
    //#pragma omp parallel for
    for (int i=0; i<nent; i++)         
    {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++) {
            int e = colent2elem[rowent2elem[i]+k];
            for (int j=0; j<npe; j++)
                ucg[i] += udg[j+npe*e]; 
        }
        ucg[i] = ucg[i]/((T) (nelem*npe));
    }
}

template void cpuGetArrayAtIndex(double*, double*, int*, int);
template void cpuPutArrayAtIndex(double*, double*, int*, int);
template void cpuArrayPlusXAtIndex(double*, double*, int*, int);
template void cpuArrayMinusXAtIndex(double*, double*, int*, int);
template void cpuArrayAXPYAtIndex(double*, double*, double, int*, int);

template double cpuArrayMaxError(double*, double *, int);
template double cpuArrayMin(double*, int);
template double cpuArrayMax(double*, int);
template double cpuArraySum(double*, int);
template double cpuArraySquareSum(double*, int);
template double cpuArrayMean(double*, int);
template void cpuArraySetValue(double*, double, int);
template void cpuArrayAddScalar(double*, double, int);
template void cpuArrayMultiplyScalar(double*, double, int);

template void cpuArrayCopy(double*, double*, int);
template void cpuArrayMinus(double*, double*, int);
template void cpuArrayAbs(double*, double*, int);
template void cpuArraySqrt(double*, double*, int);
template void cpuArraySin(double*, double*, int);
template void cpuArrayCos(double*, double*, int);
template void cpuArrayTan(double*, double*, int);
template void cpuArrayAsin(double*, double*, int);
template void cpuArrayAcos(double*, double*, int);
template void cpuArrayAtan(double*, double*, int);
template void cpuArraySinh(double*, double*, int);
template void cpuArrayCosh(double*, double*, int);
template void cpuArrayTanh(double*, double*, int);
template void cpuArrayAsinh(double*, double*, int);
template void cpuArrayAcosh(double*, double*, int);
template void cpuArrayAtanh(double*, double*, int);
template void cpuArrayExp(double*, double*, int);
template void cpuArrayLog(double*, double*, int);
template void cpuArrayCeil(double*, double*, int);
template void cpuArrayFloor(double*, double*, int);
template void cpuArrayErf(double*, double*, int);
template void cpuArrayErfc(double*, double*, int);
template void cpuArraySquare(double*, double*, int);
template void cpuArrayPower(double*, double*, int, int);

template void cpuArrayMultiplyScalarDiagonal(double*, double, int);
template void cpuArrayAddVectorToDiagonal(double*, double*, double, int);
template void cpuArrayRowAverage(double*, double*, int, int);
template void cpuArrayRowSum(double*, double*, int, int);
template void cpuArrayRowSquareSum(double*, double*, int, int);
template void cpuArrayDistSquareSum(double*, double*, double*, int, int);
template void cpuArrayAXPB(double*, double*, double, double, int);
template void cpuArrayRowkAXPB(double*, double*, double, double, int, int, int);
template void cpuArrayAXPBY(double*, double*, double*, double, double, int);
template void cpuArrayAXY(double*, double*, double*, double, int);
template void cpuArrayAXYZ(double*, double*, double*, double*, double, int);
template void cpuArrayAXYPBZ(double*, double*, double*, double*, double, double, int);
template void cpuArrayAdd3Vectors(double*, double*, double*, double*, double, double, double, int);
template void cpuArrayExtract(double*, double*, int, int, int, int, int, int, int, int, int);
template void cpuArrayInsert(double*, double*, int, int, int, int, int, int, int, int, int);
template void cpuArrayGemmBatch(double*, double*, double*, int, int, int, int);
template void cpuArrayGemmBatch1(double*, double*, double*, int, int, int, int);
template void cpuArrayDG2CG(double*, double*, int*, int*, int);
template void cpuArrayDG2CG2(double*, double*, int*, int*, int, int);

template void cpuGetArrayAtIndex(float*, float*, int*, int);
template void cpuPutArrayAtIndex(float*, float*, int*, int);
template void cpuArrayPlusXAtIndex(float*, float*, int*, int);
template void cpuArrayMinusXAtIndex(float*, float*, int*, int);
template void cpuArrayAXPYAtIndex(float*, float*, float, int*, int);

template float cpuArrayMaxError(float*, float*, int);
template float cpuArrayMin(float*, int);
template float cpuArrayMax(float*, int);
template float cpuArraySum(float*, int);
template float cpuArraySquareSum(float*, int);
template float cpuArrayMean(float*, int);
template void cpuArraySetValue(float*, float, int);
template void cpuArrayAddScalar(float*, float, int);
template void cpuArrayMultiplyScalar(float*, float, int);

template void cpuArrayCopy(float*, float*, int);
template void cpuArrayMinus(float*, float*, int);
template void cpuArrayAbs(float*, float*, int);
template void cpuArraySqrt(float*, float*, int);
template void cpuArraySin(float*, float*, int);
template void cpuArrayCos(float*, float*, int);
template void cpuArrayTan(float*, float*, int);
template void cpuArrayAsin(float*, float*, int);
template void cpuArrayAcos(float*, float*, int);
template void cpuArrayAtan(float*, float*, int);
template void cpuArraySinh(float*, float*, int);
template void cpuArrayCosh(float*, float*, int);
template void cpuArrayTanh(float*, float*, int);
template void cpuArrayAsinh(float*, float*, int);
template void cpuArrayAcosh(float*, float*, int);
template void cpuArrayAtanh(float*, float*, int);
template void cpuArrayExp(float*, float*, int);
template void cpuArrayLog(float*, float*, int);
template void cpuArrayCeil(float*, float*, int);
template void cpuArrayFloor(float*, float*, int);
template void cpuArrayErf(float*, float*, int);
template void cpuArrayErfc(float*, float*, int);
template void cpuArraySquare(float*, float*, int);
template void cpuArrayPower(float*, float*, int, int);

template void cpuArrayMultiplyScalarDiagonal(float*, float, int);
template void cpuArrayAddVectorToDiagonal(float*, float*, float, int);
template void cpuArrayRowAverage(float*, float*, int, int);
template void cpuArrayRowSum(float*, float*, int, int);
template void cpuArrayRowSquareSum(float*, float*, int, int);
template void cpuArrayDistSquareSum(float*, float*, float*, int, int);
template void cpuArrayAXPB(float*, float*, float, float, int);
template void cpuArrayRowkAXPB(float*, float*, float, float, int, int, int);
template void cpuArrayAXPBY(float*, float*, float*, float, float, int);
template void cpuArrayAXY(float*, float*, float*, float, int);
template void cpuArrayAXYZ(float*, float*, float*, float*, float, int);
template void cpuArrayAXYPBZ(float*, float*, float*, float*, float, float, int);
template void cpuArrayAdd3Vectors(float*, float*, float*, float*, float, float, float, int);
template void cpuArrayExtract(float*, float*, int, int, int, int, int, int, int, int, int);
template void cpuArrayInsert(float*, float*, int, int, int, int, int, int, int, int, int);
template void cpuArrayGemmBatch(float*, float*, float*, int, int, int, int);
template void cpuArrayGemmBatch1(float*, float*, float*, int, int, int, int);
template void cpuArrayDG2CG(float*, float*, int*, int*, int);
template void cpuArrayDG2CG2(float*, float*, int*, int*, int, int);

template int cpuArrayMin(int*, int);
template int cpuArrayMax(int*, int);
template int cpuArraySum(int*, int);
template int cpuArraySquareSum(int*, int);
template void cpuArrayCopy(int*, int*, int);
template void cpuArraySetValue(int*, int, int);
template void cpuArrayAddScalar(int*, int, int);
template void cpuArrayMultiplyScalar(int*, int, int);
template void cpuArrayAXPB(int*, int*, int, int, int);
template void cpuArrayAXPBY(int*, int*, int*, int, int, int);
template void cpuArrayExtract(int*, int*, int, int, int, int, int, int, int, int, int);
template void cpuArrayInsert(int*, int*, int, int, int, int, int, int, int, int, int);

template void cpuArrayTranspose(double *A, double *B, int m, int n);
template void cpuArrayTranspose(float *A, float *B, int m, int n);
template void cpuArrayTranspose(int *A, int *B, int m, int n);

template void cpuGetArrayAtColumnIndex(double *A, double *B, int *colind, int m, int n);
template void cpuGetArrayAtColumnIndex(float *A, float *B, int *colind, int m, int n);
template void cpuGetArrayAtColumnIndex(int *A, int *B, int *colind, int m, int n);

template void cpuArrayTransposeAtColumnIndex(double *A, double *B, int *colind, int m, int n);
template void cpuArrayTransposeAtColumnIndex(float *A, float *B, int *colind, int m, int n);
template void cpuArrayTransposeAtColumnIndex(int *A, int *B, int *colind, int m, int n);

template void cpuGetArrayAtRowIndex(double *A, double *B, int *colind, int m, int n);
template void cpuGetArrayAtRowIndex(float *A, float *B, int *colind, int m, int n);
template void cpuGetArrayAtRowIndex(int *A, int *B, int *colind, int m, int n);

template void cpuArrayTransposeAtRowIndex(double *A, double *B, int *colind, int m, int n);
template void cpuArrayTransposeAtRowIndex(float *A, float *B, int *colind, int m, int n);
template void cpuArrayTransposeAtRowIndex(int *A, int *B, int *colind, int m, int n);

template void cpuArraySumEveryColumn(int*, int*, int, int);
template void cpuArraySumEveryColumn(float*, float*, int, int);
template void cpuArraySumEveryColumn(double*, double*, int, int);

template void cpuArrayMaxEveryColumn(int*, int*,  int, int);
template void cpuArrayMaxEveryColumn(float*, float*,  int, int);
template void cpuArrayMaxEveryColumn(double*, double*, int, int);

template void cpuArrayMinEveryColumn(int*, int*,  int, int);
template void cpuArrayMinEveryColumn(float*, float*,  int, int);
template void cpuArrayMinEveryColumn(double*, double*, int, int);

template void cpuArraySumEveryRow(int*, int*, int, int);
template void cpuArraySumEveryRow(float*, float*, int, int);
template void cpuArraySumEveryRow(double*, double*, int, int);

template void cpuArrayMaxEveryRow(int*, int*, int, int);
template void cpuArrayMaxEveryRow(float*, float*, int, int);
template void cpuArrayMaxEveryRow(double*, double*, int, int);

template void cpuArrayMinEveryRow(int*, int*, int, int);
template void cpuArrayMinEveryRow(float*, float*, int, int);
template void cpuArrayMinEveryRow(double*, double*, int, int);


#endif


