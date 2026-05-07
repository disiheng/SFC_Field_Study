
#ifndef UTILS_H
#define UTILS_H

#include <time.h>

void MSG(  const char* format, time_t t0, ... ) ;
void MSG0( const char* format, time_t t0, int task, ... ) ;
void ERR(  const char* format, ... ) ;
void DBG0( const char* format, int mydbg, time_t t0, int task, ... ); 
void DBG(  const char* format, int mydbg, time_t t0, ... ) ;

size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
FILE *my_fopen(char *file, char *mode);
int wrap(int i, int n);
float fwrap(float x, float size);

#endif
