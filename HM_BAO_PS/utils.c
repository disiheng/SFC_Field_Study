#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>

#include "utils.h"

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

/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) has occured.\n");
      fflush(stdout);
      exit(0);
    }
  return nwritten;
}

int wrap(int i, int n)
{
  int ii;
  ii = i;
  if(i >= n)
    ii -= n;
  if(i < 0)
    ii += n;
  return ii;
}

float fwrap(float x, float size)
{
  float xx;
  xx = x;
  if(x >= size)
    xx -= size;
  if(x < 0.0)
    xx += size;
  return xx;
}

/* This catches I/O errors occuring for fread(). In this case we better stop.
 */

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) has occured.\n");
      fflush(stdout);
      exit(0);
    }
  return nread;
}


FILE *my_fopen(char *file, char *mode)
{

  FILE *fp;
  fp = fopen(file, mode);
  if(fp == NULL)
    {
      printf(" Problem opening %s \n", file);
      exit(0);
    }
  else
    {
      if(fp)
	{
	  return (fp);
	}
      else
	{
	  printf(" Couldn't open %s \n", file);
	  exit(0);
	}
    }
}
