
#ifndef GADGET_H
#define GADGET_H

long get_number_of_particles(char *, int Num);
long long get_total_number_of_particles(char *, int Num);
long get_number_of_groups(char *, int, int);
long get_number_of_groups_total(char *, int);

int get_subhalo_file(long *FileMap, int Num, int NFiles, int type );
double get_boxsize(char *filename, int Num);
int get_total_number_of_groups(char *, int );
int get_group_catalogue(char *, int,  int *,int *, int *);
float get_redshift(char *, int Num);
int get_nfiles(char *);
float get_mass(char *filename);

float get_hashcell_size(char *, char *, int , int, int *, int *);
int get_hash_table_size(char *, int , char *, int * );
int get_hash_table(char *, int , char *, int *, int *, int *, int *);

int get_group_coordinates(char *, int, char *, \
                          int *, int *,int, int *, \
                          int *, int, int, int, \
                          float *,  float *,  float *, \
                          float *, float *, float * );

int id_sort_compare_key(const void *a, const void *b);
int id_sort_groups(const void *a, const void *b);
double fof_periodic_wrap(double x, double );
double fof_periodic(double x, double);

typedef int peanokey;

peanokey peano_hilbert_key(int x, int y, int z, int bits);
void peano_hilbert_key_inverse(peanokey key, int bits, int *x, int *y, int *z);
 
int millenniumdir(int Num, int readFile);
int hrdir(int Num, int readFile);

#endif
//define FORMAT2
