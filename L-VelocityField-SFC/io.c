#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <stdarg.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "allvars.h"
#include "proto.h"
#include "io.h"


void MSG(const char *format, time_t t0, ...)
{
  va_list args;

  printf(" [%7.2f] ", difftime(time(NULL), t0));
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  printf("\n");
  fflush(0);
}

void MSG0(const char *format, time_t t0, int task, ...)
{
  va_list args;

  if(task == 0)
    {
      printf(" [%7.2f] ", difftime(time(NULL), t0));
      va_start(args, format);
      vprintf(format, args);
      va_end(args);
      printf("\n");
      fflush(0);
    }
}

void ERR(const char *format, ...)
{
  va_list args;

  printf(" **ERR** ");
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  printf("\n");
  fflush(0);
}
void DBG(const char *format, int mydbg, time_t t0, ...)
{
  va_list args;

  if(mydbg == 1)
    {
      printf(" --DBG-- [%7.2f] ", difftime(time(NULL), t0));
      va_start(args, format);
      vprintf(format, args);
      va_end(args);
      printf("\n");
      fflush(0);
    }
}

void DBG0(const char *format, int mydbg, time_t t0, int task, ...)
{
  va_list args;

  if(mydbg == 1 && task == 0)
    {
      printf(" --DBG-- [%7.2f] ", difftime(time(NULL), t0));
      va_start(args, format);
      vprintf(format, args);
      va_end(args);
      printf("\n");
      fflush(0);
    }
}



/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

#define MAXREAD 1000000000

  if(size * nmemb > 0)
    {
      if(size * nmemb > MAXREAD)
	{
	  nmemb = size * nmemb;
	  size = 1;

	  while(nmemb > MAXREAD)
	    {
	      if((nwritten = fwrite(ptr, size, MAXREAD, stream)) != MAXREAD)
		{
		  printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
		  fflush(stdout);
		  exit(777);
		}

	      ptr += MAXREAD;
	      nmemb -= MAXREAD;
	    }

	}

      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
	  printf("I/O error (fwrite) has occured: %s\n", strerror(errno));
	  fflush(stdout);
	  exit(777);
	}
    }
  else
    nwritten = 0;

  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(size * nmemb == 0)
    return 0;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) has occured: end of file\n");
      else
	printf("I/O error (fread) has occured: %s\n", strerror(errno));
      fflush(stdout);
      exit(778);
    }
  return nread;
}

FILE *my_fopen(char *file, char *mode)
{

    FILE *fp;
    fp = fopen(file, mode);
    if (fp == NULL) {
            printf(" Problem opening %s \n",file);
        exit(0);
        }
    else
    {
        if (fp)
        {
            return(fp);
        }
        else
        {
                printf(" Couldn't open %s \n",file);
            exit(0);
        }
        }
        //return(fp);
}



