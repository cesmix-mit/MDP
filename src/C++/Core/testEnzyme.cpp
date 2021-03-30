#include <iostream>
#include <math.h>

//!/usr/local/Cellar/llvm/11.1.0/bin/clang++ -std=c++11 -ffast-math -O3 testEnzyme.cpp -Xclang -load -Xclang /usr/local/lib/ClangEnzyme-11.dylib -lm
        
using namespace std;

template <typename... Args>
void __enzyme_autodiff(void*, Args... args);
int enzyme_const, enzyme_dup;

void print2darray(double* a, int m, int n)
{
    //cout.precision(4);
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++)
            cout << scientific << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3darray(double* a, int m, int n, int p)
{
    //cout.precision(8);
    for (int k=0; k<p; k++) {
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++)
                cout << scientific << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

template <typename T> void u1(T *__restrict__ u, T *__restrict__ x, int ncu, int N)
{
    for (int i = 0; i<N; i++) {
        u[0+ncu*i] = x[i]*x[i];
        u[1+ncu*i] = sin(x[i]);
        u[2+ncu*i] = cos(x[i]);
    }        
}
template void u1(double *, double *, int, int);
template void u1(float *, float *, int, int);

template <typename T> void du1(T *__restrict__ u, T *__restrict__ du, 
        T *__restrict__ x, T *__restrict__ dx, int ncu, int N)
{
    __enzyme_autodiff((void*)u1<T>, 
                        enzyme_dup, u, du,
                        enzyme_dup, x, dx,
                        ncu, N);    
}
template void du1(double *, double *, double *, double *,  int, int);
template void du1(float *, float *, float *, float *,  int, int);

template <typename T> void calu(T *__restrict__ u, T *__restrict__ x, int ncu, int ncx, int N)
{
    for (int i = 0; i<N; i++) {
        T x1 = x[0 + ncx*i];
        T x2 = x[1 + ncx*i];
        T x3 = x[2 + ncx*i];
        u[0+ncu*i] = x1*x1 + x2*x2 + x3*x3;
        u[1+ncu*i] = sin(x1+x2+x3);
    }        
}
template void calu(double *, double *, int, int, int);
template void calu(float *, float *, int, int, int);

template <typename T> void caldu(T *__restrict__ u, T *__restrict__ du, 
        T *__restrict__ x, T *__restrict__ dx, int ncu, int ncx, int N)
{
    __enzyme_autodiff((void*)calu<T>, 
                        enzyme_dup, u, du,
                        enzyme_dup, x, dx,
                        ncu, ncx, N);    
}
template void caldu(double *, double *, double *, double *, int, int, int);
template void caldu(float *, float *, float *, float *, int, int, int);

template <typename T> void caluv(T *__restrict__ u, T *__restrict__ v, T *__restrict__ x, int ncx, int N)
{
    for (int i = 0; i<N; i++) {
        T x1 = x[0 + ncx*i];
        T x2 = x[1 + ncx*i];
        T x3 = x[2 + ncx*i];
        u[i] = x1*x1 + x2*x2 + x3*x3;
        v[i] = sin(x1+x2+x3);
    }        
}
template void caluv(double *, double *, double *, int, int);
template void caluv(float *, float *, float *, int, int);


template <typename T> void calduv(T *__restrict__ u, T *__restrict__ du, 
        T *__restrict__ v, T *__restrict__ dv, 
        T *__restrict__ x, T *__restrict__ dx, 
        int ncx, int N)
{
    __enzyme_autodiff((void*)caluv<T>, 
                        enzyme_dup, u, du,
                        enzyme_dup, v, dv,
                        enzyme_dup, x, dx, 
                        ncx, N);    
}
template void calduv(double *, double *, double *, double *,  double *, double *, int, int);
template void calduv(float *, float *, float *, float *, float *, float *, int, int);

// template <typename T> void caluv(T *__restrict__ u, T *__restrict__ x,  
//         T *__restrict__ y,  T *__restrict__ z, int ncx, int N)
// {
//     for (int i = 0; i<N; i++) {
//         T x1 = x[i];
//         T x2 = y[i];
//         T x3 = z[i];
//         u[i] = x1*x1 + x2*x2 + x3*x3;
//         v[i] = sin(x1+x2+x3);
//     }        
// }
// template void caluv(double *, double *, double *, int, int);
// template void caluv(float *, float *, float *, int, int);

#define ncu 3
#define ncx 3
#define N 10

int main() {
    double x[N];
    double u[ncu*N];
    double du[ncu*N];
    double dx[ncu*N];   
    
    for (int i = 0; i<N; i++)        
        x[i] = i*0.1;    
    
    for (int i = 0; i<ncu*N; i++) 
        du[i] = 1.0;
    for (int i = 0; i<ncu*N; i++) 
        dx[i] = 0.0;
       
    du1(u, du, x, dx, ncu, N);
        
    print2darray(x, 1, N);
    print2darray(u, ncu, N);
    print2darray(du, ncu, N);
    print2darray(dx, ncu, N);    
    
//     double x[ncx*N];
//     double u[ncu*N];
//     double du[ncu*N];
//     double dx[ncx*ncu*N];   
//     
//     for (int i = 0; i<N; i++)
//         for (int j = 0; j<ncx; j++)
//             x[j + ncx*i] = j*1.0 + i*0.1;    
//     
//     for (int i = 0; i<ncu*N; i++) 
//         du[i] = 1.0;
//     for (int i = 0; i<ncx*ncu*N; i++) 
//         dx[i] = 0.0;
//        
//     caldu(u, du, x, dx, ncu, ncx, N);
//         
//     print2darray(x, ncx, N);
//     print2darray(u, ncu, N);
//     print2darray(dx, ncx, ncu*N);    

    
//     double x[ncx*N];
//     double u[ncu*N];
//     double v[ncu*N];
//     double du[ncu*N];
//     double dv[ncu*N];
//     double dudx[2*ncx*ncu*N];
//     double dvdx[2*ncx*ncu*N];
//     
//     for (int i = 0; i<N; i++)
//         for (int j = 0; j<ncx; j++)
//             x[j + ncx*i] = j*1.0 + i*0.1;    
//                 
//     for (int i = 0; i<ncu*N; i++) {
//         du[i] = 1.0;
//         dv[i] = 0.0;        
//     }    
//     for (int i = 0; i<ncx*ncu*N; i++) {
//         dudx[i] = 0.0;
//         dvdx[i] = 0.0;        
//     }    
//     
//     calduv(u, du, v, dv, x, dudx, ncx, N);    
//     print2darray(x, ncx, N);
//     print2darray(u, ncu, N);
//     print2darray(dudx, ncx, N);   
//     
//     for (int i = 0; i<ncu*N; i++) {
//         du[i] = 0.0;
//         dv[i] = 1.0;        
//     }    
//     for (int i = 0; i<ncx*ncu*N; i++) {
//         dudx[i] = 0.0;
//         dvdx[i] = 0.0;        
//     }    
//     
//     calduv(u, du, v, dv, x, dvdx, ncx, N);    
//     print2darray(v, ncu, N);
//     print2darray(dvdx, ncx, N);     
//     
    
//     for (int i = 0; i<ncu*N; i++) {
//         du[i] = 1.0;
//         dv[i] = 0.0;
//         dudx[i] = 0.0;
//         dvdx[i] = 0.0;
//     }

//     caldu(u, du, dudx, x, ncu, ncx, N);
//         
//     print2darray(x, ncx, N);
//     print2darray(u, ncu, N);
//     print2darray(dudx, ncx, N);    
//     
//     for (int i = 0; i<ncu*N; i++) {
//         du[i] = 1.0;
//         dv[i] = 0.0;
//     }    
}
