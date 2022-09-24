/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : warn.c                                         */
/*                                                                           */
/* Created:       2010/09/14 (JLe)                                           */
/* Last modified: 2017/03/05 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Prints warning message                                       */
/*                                                                           */
/* Comments: - Tähän yhdeksi argumentiksi pointteri laskuriin?               */
/*           - Noi messaget vois kerätä myös erilliseen fileen               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Warn:"

/*****************************************************************************/

void Warn(char *func, char *format, ...)
{
  va_list argp;
  va_start (argp, format);

  /* Print warning message */

  if (mpitasks > 1)
    fprintf(errp, "\n***** %s (seed = %lu, MPI task = %d)\n", TimeStamp(), 
            parent_seed, mpiid);
  else
    fprintf(errp, "\n***** %s (seed = %lu)\n", TimeStamp(), parent_seed);
  
  fprintf(errp, "Warning message from function %s\n\n", func);
  vfprintf(errp, format, argp);
  fprintf(errp, "\n\n");
}

/*****************************************************************************/
