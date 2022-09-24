/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectbuf.c                                   */
/*                                                                           */
/* Created:       2010/12/27 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects scoring buffer in reproducible MPI mode             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectBuf:"

/*****************************************************************************/

void CollectBuf()
{
#ifdef MPI

  double *buff;
  long sz;

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "BUF array not ready for access");

  /* Check number of tasks */

  if (mpitasks == 1)
    return;

  /* Check for domain decomposition and MPI reproducibility option */

  if (((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO) &&
      ((long)RDB[DATA_DD_DECOMPOSE] == NO))
    return;

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

  /* Get size of results data block */

  sz = (long)RDB[DATA_ALLOC_BUF_SIZE];

  /* Allocate memory for results */

  if (mpiid == 0)
    buff = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
  else
    buff = NULL;

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Reduce data */

  MPITransfer(BUF, buff, sz, 0, MPI_METH_RED);

  /* Move data to original block */

  if (mpiid == 0)
    {
      /* Copy data */

      memcpy(BUF, buff, sz*sizeof(double));

      /* Free buffer */

      Mem(MEM_FREE, buff);
    }

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Broadcast data to other tasks */

  MPITransfer(BUF, NULL, sz, 0, MPI_METH_BC);

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
}

/*****************************************************************************/
