#include <stdlib.h>
#include "byte_swap.h"

/* This byte swaps an array of 4-byte ints */
void byte_swap_int(int n, int *ptr)
{
  int i;
  char *ch;
  char temp;

  for(i = 0; i < n; i++)
    {
      ch = (char *) (ptr + i);
      temp = *ch;
      *ch = *(ch + 3);
      *(ch + 3) = temp;
      temp = *(ch + 1);
      *(ch + 1) = *(ch + 2);
      *(ch + 2) = temp;
    }
}

/* This byte swaps an array of 4-byte ints */
void byte_swap_long(int n, long *ptr)
{
  int i;
  char *ch;
  char temp;

  for(i = 0; i < n; i++)
    {
      ch = (char *) (ptr + i);
      temp = *ch;
      *ch = *(ch + 3);
      *(ch + 3) = temp;
      temp = *(ch + 1);
      *(ch + 1) = *(ch + 2);
      *(ch + 2) = temp;
    }
}

/* This byte swaps an array of 4-byte ints */
void byte_swap_llong(int n, unsigned long long *ptr)
{
  int i;
  char *ch;
  char temp;

  for(i = 0; i < n; i++)
    {

      ch = (char *) (ptr + i);
      temp = *ch;
      *ch = *(ch + 7);
      *(ch + 7) = temp;

      temp = *(ch + 1);
      *(ch + 1) = *(ch + 6);
      *(ch + 6) = temp;

      temp = *(ch + 2);
      *(ch + 2) = *(ch + 5);
      *(ch + 5) = temp;

      temp = *(ch + 3);
      *(ch + 3) = *(ch + 4);
      *(ch + 4) = temp;
    }
}

/* This byte swaps an array of 4-byte floats */
void byte_swap_float(int n, float *ptr)
{
  int i;
  char *ch;
  char temp;

  for(i = 0; i < n; i++)
    {
      ch = (char *) (ptr + i);
      temp = *ch;
      *ch = *(ch + 3);
      *(ch + 3) = temp;
      temp = *(ch + 1);
      *(ch + 1) = *(ch + 2);
      *(ch + 2) = temp;
    }
}

/* This byte swaps an array of 8-byte doubles */
void byte_swap_double(int n, double *ptr)
{
  int i;
  char *ch;
  char temp;

  for(i = 0; i < n; i++)
    {
      ch = (char *) (ptr + i);
      temp = *ch;
      *ch = *(ch + 7);
      *(ch + 7) = temp;

      temp = *(ch + 1);
      *(ch + 1) = *(ch + 6);
      *(ch + 6) = temp;

      temp = *(ch + 2);
      *(ch + 2) = *(ch + 5);
      *(ch + 5) = temp;

      temp = *(ch + 3);
      *(ch + 3) = *(ch + 4);
      *(ch + 4) = temp;
    }
}
