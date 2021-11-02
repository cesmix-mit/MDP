/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_OMPARRAYOPERATIONS
#define MDP_OMPARRAYOPERATIONS

void ompTripletnum(int* output, int* input, int length) 
{	
    #pragma omp parallel for
	for (int ii = 1; ii < length; ++ii)	
		output[ii] = (input[ii]-1)*input[ii]/2;       
}
void ompQuadrupletnum(int* output, int* input, int length) 
{	
    #pragma omp parallel for
	for (int ii = 1; ii < length; ++ii)	
		output[ii] = (input[ii]-2)*(input[ii]-1)*input[ii]/6;       
}

void ompCumsum(int* output, int* input, int length) 
{
	output[0] = 0; 
	for (int j = 1; j < length; ++j)	
		output[j] = input[j - 1] + output[j - 1];	
}

void ompArrayFill(int* output, int start, int length) 
{	
    #pragma omp parallel for
	for (int j = 0; j < length; ++j)	
		output[j] = start + j;
}

// template <typename T> T ompArrayErrorNorm(T *a, T *b, int n)
// {
//     T e = (a[0]-b[0])*(a[0]-b[0]);
//     for (int i=1; i<n; i++)        
//         e += (a[i]-b[i])*(a[i]-b[i]);    
//     return sqrt(e);
// }

template <typename T> T ompArrayMaxError(T *a, T *b, int n)
{
    T e = fabs(a[0]-b[0]);
    for (int i=1; i<n; i++)
        if (fabs(a[i]-b[i])>e)
            e = fabs(a[i]-b[i]);    
    return e;
}

// template <typename T> T ompArrayNorm(T *a, int n)
// {
//     T e = a[0]*a[0];
//     for (int i=1; i<n; i++)        
//         e += a[i]*a[i];    
//     return sqrt(e);
// }

template <typename T> T ompArrayMin(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

template <typename T> T ompArrayMax(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

template <typename T> T ompArraySum(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    return b;
}

template <typename T> T ompArraySquareSum(T *a, int n)
{
    T b = a[0]*a[0];
    for (int i=1; i<n; i++)
            b = b + a[i]*a[i];    
    return b;
}

template <typename T> T ompArrayMean(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    b = b/((T) n);
    return b;
}

template <typename T> void ompGetArrayAtIndex(T *y, T *x, int *ind, int n)
{    
    #pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[i] = x[ind[i]];
}

template <typename T> void ompPutArrayAtIndex(T *y, T *x, int *ind, int n)
{    
    #pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] = x[i];
}

template <typename T> void ompArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{    
    #pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] += x[i];
}

template <typename T> void ompArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{    
    #pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] -= x[i];
}

template <typename T> void ompArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{    
    #pragma omp parallel for
    for (int i = 0; i<n; i++)    
        y[ind[i]] = y[ind[i]] + a*x[i];
}

template <typename T> void ompArraySetValue(T *y, T a, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a;        
}

template <typename T> void ompArrayMultiplyScalar(T *y, T a, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a*y[i];        
}

template <typename T> void ompArrayAddScalar(T *y, T a, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = y[i] + a;        
}

template <typename T> void ompArrayCopy(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = x[i];        
}

template <typename T> void ompArrayMinus(T *y, T *x, int n)
{    
    #pragma omp parallel for    
    for (int i=0; i<n; i++) 
        y[i] = -x[i];        
}

template <typename T> void ompArrayAbs(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = fabs(x[i]);        
}

template <typename T> void ompArraySqrt(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = sqrt(x[i]);        
}

template <typename T> void ompArraySin(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = sin(x[i]);        
}

template <typename T> void ompArrayCos(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = cos(x[i]);        
}

template <typename T> void ompArrayTan(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = tan(x[i]);        
}

template <typename T> void ompArrayAsin(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = asin(x[i]);        
}

template <typename T> void ompArrayAcos(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = acos(x[i]);        
}

template <typename T> void ompArrayAtan(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = atan(x[i]);        
}

template <typename T> void ompArraySinh(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = sinh(x[i]);        
}

template <typename T> void ompArrayCosh(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = cosh(x[i]);        
}

template <typename T> void ompArrayTanh(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = tanh(x[i]);        
}

template <typename T> void ompArrayAsinh(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = asinh(x[i]);        
}

template <typename T> void ompArrayAcosh(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = acosh(x[i]);        
}

template <typename T> void ompArrayAtanh(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = atanh(x[i]);        
}

template <typename T> void ompArrayExp(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = exp(x[i]);        
}

template <typename T> void ompArrayLog(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = log(x[i]);        
}

template <typename T> void ompArrayCeil(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = ceil(x[i]);        
}

template <typename T> void ompArrayFloor(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = floor(x[i]);        
}

template <typename T> void ompArrayErf(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = erf(x[i]);        
}

template <typename T> void ompArrayErfc(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = erfc(x[i]);        
}

template <typename T> void ompArraySquare(T *y, T *x, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = x[i]*x[i];        
}

template <typename T> void ompArrayPower(T *y, T *x, int p, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        y[i] = x[i];    
        for (int j=1; j<p; j++)
            y[i] = y[i]*x[i];
    }
}

template <typename T> void ompArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    #pragma omp parallel for
    for (int i=0; i<n; i++)    
        C[i+n*i] = a*C[i+n*i];                
}

template <typename T> void ompArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{            
    #pragma omp parallel for
    for (int i=0; i<n; i++)    
        C[i+n*i] += a*x[i];                
}

template <typename T> void ompArrayRowAverage(T *y, T *x, int m, int n)
{    
    #pragma omp parallel 
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

template <typename T> void ompArrayRowkAXPB(T *y, T *x, T a, T b, int m, int n, int k)
{    
    #pragma omp parallel 
    for (int j=0; j<n; j++)    
        y[k + m*j] = a*x[k + m*j] + b;    
}

template <typename T> void ompArrayAXPB(T *y, T *x, T a, T b, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        y[i] = a*x[i]+b;        
}

template <typename T> void ompArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        z[i] = a*x[i] + b*y[i];        
}

template <typename T> void ompArrayAXY(T *s, T *x, T *y, T a, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i];        
}

template <typename T> void ompArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i]*z[i];        
}

template <typename T> void ompArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{    
    #pragma omp parallel for    
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i] + b*z[i];        
}

template <typename T> void ompArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{    
    #pragma omp parallel for        
    for (int i=0; i<n; i++) 
        s[i] = a*x[i] + b*y[i] + c*z[i];        
}

template <typename T> void ompArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        un[idx] = u[i+I*j+I*J*k];
    }            
}

template <typename T> void ompArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    #pragma omp parallel for
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+I*J*k] = un[idx];        
    }            
}

template <typename T> void ompArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    #pragma omp parallel for
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

template <typename T> void ompArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    #pragma omp parallel for
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

template <typename T> void ompArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{    
    #pragma omp parallel for
    for (int i=0; i<nent; i++)         
    {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++)
            ucg[i] += udg[cgent2dgent[rowent2elem[i]+k]]; 
        ucg[i] = ucg[i]/((T) nelem);
    }
}

template <typename T> void ompArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{    
    #pragma omp parallel for
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

template void ompGetArrayAtIndex(double*, double*, int*, int);
template void ompPutArrayAtIndex(double*, double*, int*, int);
template void ompArrayPlusXAtIndex(double*, double*, int*, int);
template void ompArrayMinusXAtIndex(double*, double*, int*, int);
template void ompArrayAXPYAtIndex(double*, double*, double, int*, int);

template double ompArrayMaxError(double*, double *, int);
template double ompArrayMin(double*, int);
template double ompArrayMax(double*, int);
template double ompArraySum(double*, int);
template double ompArraySquareSum(double*, int);
template double ompArrayMean(double*, int);
template void ompArraySetValue(double*, double, int);
template void ompArrayAddScalar(double*, double, int);
template void ompArrayMultiplyScalar(double*, double, int);

template void ompArrayCopy(double*, double*, int);
template void ompArrayMinus(double*, double*, int);
template void ompArrayAbs(double*, double*, int);
template void ompArraySqrt(double*, double*, int);
template void ompArraySin(double*, double*, int);
template void ompArrayCos(double*, double*, int);
template void ompArrayTan(double*, double*, int);
template void ompArrayAsin(double*, double*, int);
template void ompArrayAcos(double*, double*, int);
template void ompArrayAtan(double*, double*, int);
template void ompArraySinh(double*, double*, int);
template void ompArrayCosh(double*, double*, int);
template void ompArrayTanh(double*, double*, int);
template void ompArrayAsinh(double*, double*, int);
template void ompArrayAcosh(double*, double*, int);
template void ompArrayAtanh(double*, double*, int);
template void ompArrayExp(double*, double*, int);
template void ompArrayLog(double*, double*, int);
template void ompArrayCeil(double*, double*, int);
template void ompArrayFloor(double*, double*, int);
template void ompArrayErf(double*, double*, int);
template void ompArrayErfc(double*, double*, int);
template void ompArraySquare(double*, double*, int);
template void ompArrayPower(double*, double*, int, int);

template void ompArrayMultiplyScalarDiagonal(double*, double, int);
template void ompArrayAddVectorToDiagonal(double*, double*, double, int);
template void ompArrayRowAverage(double*, double*, int, int);
template void ompArrayAXPB(double*, double*, double, double, int);
template void ompArrayRowkAXPB(double*, double*, double, double, int, int, int);
template void ompArrayAXPBY(double*, double*, double*, double, double, int);
template void ompArrayAXY(double*, double*, double*, double, int);
template void ompArrayAXYZ(double*, double*, double*, double*, double, int);
template void ompArrayAXYPBZ(double*, double*, double*, double*, double, double, int);
template void ompArrayAdd3Vectors(double*, double*, double*, double*, double, double, double, int);
template void ompArrayExtract(double*, double*, int, int, int, int, int, int, int, int, int);
template void ompArrayInsert(double*, double*, int, int, int, int, int, int, int, int, int);
template void ompArrayGemmBatch(double*, double*, double*, int, int, int, int);
template void ompArrayGemmBatch1(double*, double*, double*, int, int, int, int);
template void ompArrayDG2CG(double*, double*, int*, int*, int);
template void ompArrayDG2CG2(double*, double*, int*, int*, int, int);

template void ompGetArrayAtIndex(float*, float*, int*, int);
template void ompPutArrayAtIndex(float*, float*, int*, int);
template void ompArrayPlusXAtIndex(float*, float*, int*, int);
template void ompArrayMinusXAtIndex(float*, float*, int*, int);
template void ompArrayAXPYAtIndex(float*, float*, float, int*, int);

template float ompArrayMaxError(float*, float*, int);
template float ompArrayMin(float*, int);
template float ompArrayMax(float*, int);
template float ompArraySum(float*, int);
template float ompArraySquareSum(float*, int);
template float ompArrayMean(float*, int);
template void ompArraySetValue(float*, float, int);
template void ompArrayAddScalar(float*, float, int);
template void ompArrayMultiplyScalar(float*, float, int);

template void ompArrayCopy(float*, float*, int);
template void ompArrayMinus(float*, float*, int);
template void ompArrayAbs(float*, float*, int);
template void ompArraySqrt(float*, float*, int);
template void ompArraySin(float*, float*, int);
template void ompArrayCos(float*, float*, int);
template void ompArrayTan(float*, float*, int);
template void ompArrayAsin(float*, float*, int);
template void ompArrayAcos(float*, float*, int);
template void ompArrayAtan(float*, float*, int);
template void ompArraySinh(float*, float*, int);
template void ompArrayCosh(float*, float*, int);
template void ompArrayTanh(float*, float*, int);
template void ompArrayAsinh(float*, float*, int);
template void ompArrayAcosh(float*, float*, int);
template void ompArrayAtanh(float*, float*, int);
template void ompArrayExp(float*, float*, int);
template void ompArrayLog(float*, float*, int);
template void ompArrayCeil(float*, float*, int);
template void ompArrayFloor(float*, float*, int);
template void ompArrayErf(float*, float*, int);
template void ompArrayErfc(float*, float*, int);
template void ompArraySquare(float*, float*, int);
template void ompArrayPower(float*, float*, int, int);

template void ompArrayMultiplyScalarDiagonal(float*, float, int);
template void ompArrayAddVectorToDiagonal(float*, float*, float, int);
template void ompArrayRowAverage(float*, float*, int, int);
template void ompArrayAXPB(float*, float*, float, float, int);
template void ompArrayRowkAXPB(float*, float*, float, float, int, int, int);
template void ompArrayAXPBY(float*, float*, float*, float, float, int);
template void ompArrayAXY(float*, float*, float*, float, int);
template void ompArrayAXYZ(float*, float*, float*, float*, float, int);
template void ompArrayAXYPBZ(float*, float*, float*, float*, float, float, int);
template void ompArrayAdd3Vectors(float*, float*, float*, float*, float, float, float, int);
template void ompArrayExtract(float*, float*, int, int, int, int, int, int, int, int, int);
template void ompArrayInsert(float*, float*, int, int, int, int, int, int, int, int, int);
template void ompArrayGemmBatch(float*, float*, float*, int, int, int, int);
template void ompArrayGemmBatch1(float*, float*, float*, int, int, int, int);
template void ompArrayDG2CG(float*, float*, int*, int*, int);
template void ompArrayDG2CG2(float*, float*, int*, int*, int, int);

template int ompArrayMin(int*, int);
template int ompArrayMax(int*, int);
template int ompArraySum(int*, int);
template int ompArraySquareSum(int*, int);
template void ompArrayCopy(int*, int*, int);
template void ompArraySetValue(int*, int, int);
template void ompArrayAddScalar(int*, int, int);
template void ompArrayMultiplyScalar(int*, int, int);
template void ompArrayAXPB(int*, int*, int, int, int);
template void ompArrayAXPBY(int*, int*, int*, int, int, int);
template void ompArrayExtract(int*, int*, int, int, int, int, int, int, int, int, int);
template void ompArrayInsert(int*, int*, int, int, int, int, int, int, int, int, int);

#endif


