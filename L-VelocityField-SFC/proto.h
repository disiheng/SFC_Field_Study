
#ifndef PROTO_H
#define PROTO_H

#include "allvars.h"

double density_statistics(fftw_real *grid);
double log_density( float *grid);
void read_parameter_file(char *fname);
void velocity_field(int dim);
void save_grids(int dim);
void get_slabs_file(void);
void load_grid_field(int value);
void load_catalogue(void);
void load_snapshot(void);
void load_halo_catalogue_single(int FileNr);
void density(Particle * P, fftw_real *grid, long long NumPart, int mode, int value);
void mapping(Particle * P, fftw_real *grid, long long FirstPart, long long NumPart, int mode, int value);
void init_grid(void);
void clean_grid(void);

double D1(double astart, double aend);
double D2(double astart, double aend);
double growth(double a);
double growth_int(double a, void *param);
double F_Omega(double a);
double F2_Omega(double a);
double hubble_function(double a);

FILE *my_fopen(char *file, char *mode);

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);

void mymalloc_init(void);
void dump_memory_table(void);
void report_detailed_memory_usage_of_largest_task(size_t *OldHighMarkBytes, const char *label, const char *func, const char *file, int line);

void MSG(  const char* format, time_t t0, ... ) ;
void MSG0( const char* format, time_t t0, int task, ... ) ;
void ERR(  const char* format, ... ) ;
void DBG0( const char* format, int mydbg, time_t t0, int task, ... );
void DBG(  const char* format, int mydbg, time_t t0, ... ) ;

#endif

