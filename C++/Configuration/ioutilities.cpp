/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __IOUTILITIES
#define __IOUTILITIES

template <typename T> string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

template <typename T> T * copyarray(T *b, Int N)
{
    T *a;
    if (N>0) {        
        a = (T*) malloc (sizeof (T)*N);
        for (Int i=0; i<N; i++)
            a[i] = b[i];
    }    
    else {
        a = NULL;
    }
    return a;
}

template <typename T> void readarray(ifstream &in, T **a, Int N)
{
    if (N>0) {        
        *a = (T*) malloc (sizeof (T)*N);
        in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );                           
    }    
}

template <typename T> void readarraynomalloc(ifstream &in, T **a, Int N)
{
    if (N>0) {        
        in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );                           
    }    
}

template <typename T> T * readarray(ifstream &in, Int N)
{
    T *a;
    if (N>0) {        
        a = (T*) malloc (sizeof (T)*N);
        in.read( reinterpret_cast<char*>( a ), sizeof(T)*N );        
    }    
    else {
        a = NULL;
    }
    return a;
}

Int * readiarrayfromdouble(ifstream &in, Int N)
{
    Int *a;
    if (N>0) {      
        a = (Int*) malloc (sizeof (Int)*N);
        double read;
        for (unsigned i = 0; i < N; i++) {
            in.read( reinterpret_cast<char*>( &read ), sizeof read );
            a[i] = (Int) round(read);
        }        
    }
    else {
        a = NULL;
    }
    return a;
}

#ifdef  HAVE_CUDA

template <typename T> T compareArrayCPUandGPU(T* a, T* gpub, Int N)
{
    T diff = 0;
    T *b;
    b = (T*) malloc (sizeof (T)*N);
    cudaMemcpy(b, gpub, N*sizeof(T), cudaMemcpyDeviceToHost);    
    for (Int i=0; i<N; i++)
        if (fabs(a[i]-b[i])>diff)
            diff = fabs(a[i]-b[i]);        
    free(b);
    return diff;
}

template <typename T> void readarrayManaged(ifstream &in, T **a, Int N)
{
    if (N>0) {        
        cudaTemplateMallocManaged(a, N);
        in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );        
    }    
}

Int * readiarrayfromdoubleManaged(ifstream &in, Int N)
{
    Int *a;
    if (N>0) {                              
        CUDA_CHECK(cudaMallocManaged((void **)&a, sizeof (Int)*N));
    
        double read;
        for (unsigned i = 0; i < N; i++) {
            in.read( reinterpret_cast<char*>( &read ), sizeof read );
            a[i] = (Int) round(read);
        }        
    }    
    return a;
}

template <typename T> void readarrayZeroCopy(ifstream &in, T **a, Int N)
{
    if (N>0) {        
        cudaTemplateHostAlloc(a, N, cudaHostAllocMapped);
        in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );        
    }    
}

Int * readiarrayfromdoubleZeroCopy(ifstream &in, Int N)
{
    Int *a;
    if (N>0) {                              
        unsigned int flags = cudaHostAllocMapped;
        
        CUDA_CHECK(cudaHostAlloc((void **)&a, sizeof (Int)*N, flags));
    
        double read;
        for (unsigned i = 0; i < N; i++) {
            in.read( reinterpret_cast<char*>( &read ), sizeof read );
            a[i] = (Int) round(read);
        }        
    }    
    return a;
}

#endif

template <typename T> void writearray(ofstream &out, T *a, Int N)
{
    if (N>0)       
        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );
}

void writeiarraytodouble(ofstream &out, Int *a, Int N)
{
    if (N>0) {
        double b;
        for (unsigned i = 0; i < N; i++) {
            b = (double) a[i];
            out.write( reinterpret_cast<char*>( &b ), sizeof(double) );
        }
    }
}

void writedouble(ofstream &out, double b)
{        
    out.write( reinterpret_cast<char*>( &b ), sizeof(double) );    
}

bool fileexists(string filename) 
{
    ifstream ifile(filename.c_str());
    return (bool)ifile;
}

template <typename T> void readarrayfromfile(string filename, T **a, Int N)
{
    if (N>0) {
        // Open file to read
        ifstream in(filename.c_str(), ios::in | ios::binary);

        if (!in) {
            error("Unable to open file " + filename);
        }

        if (in) {
            *a = (T*) malloc (sizeof (T)*N);
            in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );        
        }

        in.close();
    }
}

template <typename T> void readarrayfromfile(string filename, T **a, Int N, Int backend)
{
    if (N>0) {
        // Open file to read
        ifstream in(filename.c_str(), ios::in | ios::binary);

        if (!in) {
            error("Unable to open file " + filename);
        }
        
        if (backend==2) { //GPU
#ifdef  HAVE_CUDA                        
            T *a_host;            
            a_host = (T*) malloc (sizeof (T)*N);            
            
            // read data from file
            in.read( reinterpret_cast<char*>(a_host ), sizeof(T)*N );        

            // transfer data from CPU to GPU 
            cudaMemcpy(*a, &a_host[0], N*sizeof(T), cudaMemcpyHostToDevice);    
                        
            free(a_host);
#endif            
        }
        else {
            //*a = (T*) malloc (sizeof (T)*N);
            in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );        
        }
        
        in.close();
    }
}

template <typename T> void writearray2file(string filename, T *a, Int N)
{
    if (N>0) {
        // Open file to read
        ofstream out(filename.c_str(), ios::out | ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }

        out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );

        out.close();
    }
}

template <typename T> void writearray(ofstream &out, T *a, Int N, Int backend)
{
    if (N>0) {        
        if (backend==2) { //GPU
#ifdef  HAVE_CUDA                        
            T *a_host;            
            a_host = (T*) malloc (sizeof (T)*N);            
            
            // transfer data from GPU to CPU to save in a file
            cudaMemcpy(&a_host[0], &a[0], N*sizeof(T), cudaMemcpyDeviceToHost);    
            
            out.write( reinterpret_cast<char*>( &a_host[0] ), sizeof(T) * N );
            
            free(a_host);
#endif            
        }
        else 
            out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );                            
    }
}

template <typename T> void writearray2file(string filename, T *a, Int N, Int backend)
{
    if (N>0) {        
        // Open file to read
        ofstream out(filename.c_str(), ios::out | ios::binary);

        if (!out) {
            error("Unable to open file " + filename);
        }

        if (backend==2) { //GPU
#ifdef  HAVE_CUDA                        
            T *a_host;            
            a_host = (T*) malloc (sizeof (T)*N);            
            
            // transfer data from GPU to CPU to save in a file
            cudaMemcpy(&a_host[0], &a[0], N*sizeof(T), cudaMemcpyDeviceToHost);    
            
            out.write( reinterpret_cast<char*>( &a_host[0] ), sizeof(T) * N );
            
            free(a_host);
#endif            
        }
        else 
            out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );                    
        
        out.close();
    }
}

void writeTimeStepSize2File(string filename, Int timeStep, Int DIRKstage, double time, double dt)
{
    ofstream out(filename.c_str(), ios::out | ios::app);
    if (!out)
        cout <<"Unable to open file" << filename << endl;
    else
        out << NumberToString(timeStep) + "\t" + NumberToString(DIRKstage) + "\t" + NumberToString(time) + "\t" + NumberToString(dt) + "\n";
    out.close();
}

void writeScalarField2File(string filename, double* field, Int ne, Int* ndims)
{
    Int npv = ndims[9];
    
    ofstream out(filename.c_str(), ios::out | ios::binary);
    if (!out)
        cout <<"Unable to open file" << filename << endl;
    else if (out) {
        out.write( reinterpret_cast<char*>( &field[0] ), sizeof(double) * npv*ne );
    }
    out.close();
}

#endif