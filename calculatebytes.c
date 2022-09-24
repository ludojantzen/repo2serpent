/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatebytes.c                               */
/*                                                                           */
/* Created:       2011/11/11 (JLe)                                           */
/* Last modified: 2018/03/14 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates total number of bytes                             */
/*                                                                           */
/* Comments: - To be called after each memory allocation                     */
/*           - ACE array poistettiin TOTAL_BYTES -laskurista 28.5.2012       */
/*             (2.1.6)                                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateBytes:"

/*****************************************************************************/

void CalculateBytes()
{
  long mem, real, nt;

  /* Reset sizes */

  mem = 0;
  real = 0;

  /* Get number of OpenMP threads */

  nt = (long)RDB[DATA_OMP_MAX_THREADS];

  /* ASCII array */

  mem = mem + (long)RDB[DATA_ASCII_DATA_SIZE]*sizeof(char);
  real = real + (long)RDB[DATA_ASCII_DATA_SIZE]*sizeof(char);

  /* Main data array */

  mem = mem + (long)RDB[DATA_ALLOC_MAIN_SIZE]*sizeof(double);
  real = real + (long)RDB[DATA_REAL_MAIN_SIZE]*sizeof(double);

  /* ACE array */

  real = real + (long)RDB[DATA_REAL_ACE_SIZE]*sizeof(double);

  /* First results array */

  mem = mem + (long)RDB[DATA_ALLOC_RES1_SIZE]*sizeof(double);
  real = real + (long)RDB[DATA_REAL_RES1_SIZE]*sizeof(double);
  
  /* Private array */

  mem = mem + (long)RDB[DATA_ALLOC_PRIVA_SIZE]*nt*sizeof(double);
  real = real + (long)RDB[DATA_REAL_PRIVA_SIZE]*nt*sizeof(double);

  /* Buffer array */

  if ((long)RDB[DATA_OPTI_SHARED_BUF] == YES)
    {
      mem = mem + (long)RDB[DATA_ALLOC_BUF_SIZE]*sizeof(double);
      real = real + (long)RDB[DATA_REAL_BUF_SIZE]*sizeof(double);
    }
  else
    {
      mem = mem + (long)RDB[DATA_ALLOC_BUF_SIZE]*nt*sizeof(double);
      real = real + (long)RDB[DATA_REAL_BUF_SIZE]*nt*sizeof(double);
    }

  /* RES2 array */

  if ((long)RDB[DATA_OPTI_SHARED_RES2] == YES)
    {
      mem = mem + (long)RDB[DATA_ALLOC_RES2_SIZE]*sizeof(double);
      real = real + (long)RDB[DATA_REAL_RES2_SIZE]*sizeof(double);
    }
  else
    {
      mem = mem + (long)RDB[DATA_ALLOC_RES2_SIZE]*nt*sizeof(double);
      real = real + (long)RDB[DATA_REAL_RES2_SIZE]*nt*sizeof(double);
    }

  /* RES3 array */

  mem = mem + (long)RDB[DATA_ALLOC_RES3_SIZE]*sizeof(double);
  real = real + (long)RDB[DATA_REAL_RES3_SIZE]*sizeof(double);
  
  /* Put size */

  WDB[DATA_TOTAL_BYTES] = (double)mem;
  WDB[DATA_REAL_BYTES] = (double)real;
}  

/*****************************************************************************/
