#include<fftw3.h>
#include<mkl_dfti.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstring>
#include<string>
#include<vector>
#include<chrono>
#include<random>
#include<cmath>
#include<omp.h>


#ifndef TWOPI
#define TWOPI M_PI*2
#endif

#define IN_PARALLEL

#ifndef PARAMETERS_NUM
#define PARAMETERS_NUM 10
#endif

extern int Resolution;                    // MRA scale parameter J
extern int BaseType;                      // 0 for B_Spline, 1 for Daubechies phi
extern int phiGenus;                      // Daubechies Wavelet Genus n
extern int SampRate;                      // Wavelet Phi sampling rate (points / 1)
extern int KernelFunc;                    // window function, 0:shell, 1:sphere, 2:Gaussian
extern double Radius;                     // window radius R in Mpc/h
extern double SimBoxL;                    // simulation box length in Mpc/h
extern int Threads;                       // number of threads that used
extern uint64_t GridLen;                  // side length of MRA frame, == 2^J
extern uint64_t GridNum;                  // number of cubes, == (2^J)^3
extern std::string DIREC;
extern std::string RESOL;
extern std::string RADII;
extern std::string GENUS;
extern std::string DataDirec;
extern std::vector<double> phi;

// kernel function declared
double WindowFunction_Shell(double R, double ki, double kj, double kk);
double WindowFunction_Sphere(double R, double ki, double kj, double kk);
double WindowFunction_Gaussian(double R, double ki, double kj, double kk);
double WindowFunction_Dual_Ring(double R, double theta, double ki, double kj, double kk);


// needed for binary I/O
template<class T> char* as_bytes(T& i)
{
    void* addr = &i;                            // get the address of the first byte
                                                // of memory used to store the object
    return static_cast<char*>(addr);            // treat that memory as bytes
}

// in case the .bin file store in diffrent endianness
template<class T> void readBigEndian(std::ifstream& i, T& a)
{
    for(int n = sizeof(a) - 1; n >= 0 ; --n)
    {
        i.read(((char*) &a) + n, sizeof(char));
    }
}

// Millennium Run galaxy catalog
struct Galaxy
{
    float x, y, z;                              // position in Mpc/h;
    float vx, vy, vz;                           // velocity in km/s;
    float Mag_u, Mag_g, Mag_r, Mag_i, Mag_z;    // total galaxy magnitudes in (AB) standard SDSS filters
    float BulgeMag_u, BulgeMag_g,
           BulgeMag_r, BulgeMag_i, BulgeMag_z;  // bulge magnitude only
    float StellarMass, BulgeMass, ColdGas, 
           HotGas, EjectedMass, BlackHoleMass, Sfr; // all mass in 10^10Msun/h, Sfr in Msun/yr

    int size() {return 23;}                     // number of parameter {x,y,x,vx,vy...} == 23

    // initialize Galaxy from contineous memory <float*> a
    void init(float* a)                         
    {
        x               =   a[0];
        y               =   a[1];
        z               =   a[2];
        vx              =   a[3];
        vy              =   a[4];
        vz              =   a[5];
        Mag_u           =   a[6];
        Mag_g           =   a[7];
        Mag_r           =   a[8];
        Mag_i           =   a[9];
        Mag_z           =   a[10];
        BulgeMag_u      =   a[11];
        BulgeMag_g      =   a[12];
        BulgeMag_r      =   a[13];
        BulgeMag_i      =   a[14];
        BulgeMag_z      =   a[15];
        StellarMass     =   a[16];
        BulgeMass       =   a[17];
        ColdGas         =   a[18];
        HotGas          =   a[19];
        EjectedMass     =   a[20];
        BlackHoleMass   =   a[21];
        Sfr             =   a[22];
    }
};

struct Index
{
    int i;
    int j;
    int k;
};

struct Particle
{
    double x;
    double y;
    double z;
    double weight;
};

struct Offset
{
    double dx;
    double dy;
    double dz;
    Offset(double x, double y, double z)
    {
        dx = x;
        dy = y;
        dz = z;
    }
};

// defined in mracs.cpp
void welcome();
void read_parameter();
double* sfc_offset(std::vector<Particle>& p, Offset v);
double* sfc(std::vector<Particle>& p);
std::vector<double> Daubechies_Phi(const int phiGenus);
std::vector<double> B_Spline(const int n, const int sampleRate);
double* Spectrum1(std::vector<double>& v, double k0, double k1, size_t N_k);
double* Spectrum(std::vector<double>& v, double k0, double k1, size_t N_k);
double* PowerSpectrum(std::vector<double>& v, double k0, double k1, size_t N_k);
double* PowerPhiFunc(const size_t N);
double* wfc(const double Radius, const double theta);
double* convol3d(double* s, double* w);
double* convol_c2r(fftw_complex* sc, double* w);
fftw_complex* sfc_r2c(double* s);
double array_sum(double* w, size_t N);
double inner_product(double* v0, double* v1, size_t N);
void force_kernel_type(int x);
void result_interpret(const double* s, std::vector<Particle>& p0, std::vector<double>& result);
void de_duplicate_push_back(std::vector<Index>& index, const int i, const int j, const int k);
void fill_index_set(const double R, std::vector<Index>& inner_index, std::vector<Index>& cross_index);
void count_in_sphere(const double R, const double SimBoxL, std::vector<Particle>& p, std::vector<Particle>& p0, std::vector<double>& result);
void stupid_count(const double R, std::vector<Particle>& p, std::vector<Particle>& p0, std::vector<double>& stupid_result);

