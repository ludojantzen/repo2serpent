/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initmpi.c                                      */
/*                                                                           */
/* Created:       2010/11/23 (JLe)                                           */
/* Last modified: 2018/08/08 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Initializes parallel calculation using MPI                   */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "InitMPI:"

/*****************************************************************************/

void InitMPI(int argc, char **argv)
{

#ifdef MPI

  int rc;

  /***************************************************************************/

  /****** MPI is available in compilation. Initialise. ***********************/

  rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS)
    {
      fprintf(stdout, "MPI initialization failed.\n");
      exit(-1);
    }

  /* Reset number of tasks and id */

  mpitasks = 1;
  mpiid = 0;

  /* Get number of tasks and id */

  rc = MPI_Comm_size(MPI_COMM_WORLD, &mpitasks);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  my_comm = MPI_COMM_WORLD;

  /***************************************************************************/

#else

  /***************************************************************************/

  /***** Not compiled in MPI mode. Set number of tasks to 1. *****************/

  mpitasks = 1;
  mpiid = 0;

  /***************************************************************************/

#endif
}

/*****************************************************************************/
