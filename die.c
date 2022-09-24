/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : die.c                                          */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Terminates run in fatal error                                */
/*                                                                           */
/* Comments: - From Serpent 1.1.8                                            */
/*                                                                           */
/*           - Use this for errors in code, user errors terminate the run    */
/*             with Error()                                                  */
/*                                                                           */
/*           - This function fails if called after FreeMem().                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Die:"

/*****************************************************************************/

int Die(char *func, ...)
{
  long ptr, n;
  char *format;
  va_list argp;
  va_start (argp, func);

  /* Print error message */

  fprintf(errp, "\n***** %s:\n\n", TimeStamp());

  fprintf(errp, " - MPI task         = %d\n", mpiid);
  fprintf(errp, " - OpenMP thread    = %d\n", OMP_THREAD_NUM);
  fprintf(errp, " - RNG parent seed  = %lu\n", parent_seed);
  if ((SEED0 != NULL) && (SEED0[OMP_THREAD_NUM*RNG_SZ] > 0))
    fprintf(errp, " - RNG history seed = %lu\n",
            SEED0[OMP_THREAD_NUM*RNG_SZ]);

  /* Check if simulation is running */

  if ((long)RDB[DATA_SIMULATION_COMPLETED] == NO)
    {
      /* Get history index */

      n = -1;
      ptr = (long)RDB[DATA_PTR_PRIVA_HIS_IDX];
      if (ptr > VALID_PTR)
        if ((long)RDB[DATA_PRIVA_MEM_READY] == YES)
          n = (long)GetPrivateData(ptr, OMP_THREAD_NUM);

      if (n > -1)
        fprintf(errp, " - RNG history idx  = %ld\n\n", n);
      else
        fprintf(errp, "\n");
    }
  else
    fprintf(errp, "\n");

  fprintf(errp, "Fatal error in function %s\n\n", func);

  format = va_arg(argp, char *);
  vfprintf(errp, format, argp);

  fprintf(errp, "\n\n");

  /*
  if ((int)RDB[DATA_RUNNING_MODE] == RUNNING_MODE_XSTEST)
    return 0;
  */

  if ((long)RDB[DATA_TERMINATE_ON_DIE] == YES)
    {
      fprintf(errp, "Simulation aborted.\n\n");

      /* Deinitialize the socket communicator */

      DeinitSocket();

      /* Exit with value -1 to terminate all MPI tasks */

      exit(-1);
    }
  else
    return 0;
}

/*****************************************************************************/
