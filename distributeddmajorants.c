/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : distributeddmajorants.c                        */
/*                                                                           */
/* Created:       2018/03/14 (JLe)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Broadcasts majorants to other MPI tasks and calculates       */
/*              absolute majorant.                                           */
/*                                                                           */
/* Comments: - Used only with domain decomposition in optimization modes     */
/*             2-4. When macroscopic cross sections are not pre-calculated,  */
/*             the maximum atomic densities are combined in                  */
/*             distributeddmazadens.c.                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DistributeDDMajorants:"

/*****************************************************************************/

void DistributeDDMajorants()
{

#ifdef MPI

  long loc0, ne, n, m;
  double *xs, *xs0;

  /* Check number of tasks */

  if (mpitasks == 1)
    return;

  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Check if macroscopic xs are generated or multi-group majorant is used */

  if (((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO) &&
      ((long)RDB[DATA_OPTI_MG_MODE] == NO))
    return;

  fprintf(outp, "Sharing majorants to other domains...\n");

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

  /* Check for multi-group mode */

  if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
    {
      /* Get pointer to energy grid */

      loc0 = (long)RDB[DATA_COARSE_MG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Number of points */

      ne = (long)RDB[loc0 + ENERGY_GRID_NE];

      /* Get pointer to majorant reaction data */

      loc0 = (long)RDB[DATA_PTR_MAJORANT];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Get pointer to data */

      loc0 = (long)RDB[loc0 + REACTION_PTR_MGXS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
    }
  else
    {
     /* Get pointer to energy grid */

      loc0 = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Number of points */

      ne = (long)RDB[loc0 + ENERGY_GRID_NE];

      /* Get pointer to majorant reaction data */

      loc0 = (long)RDB[DATA_PTR_MAJORANT];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Get pointer to data */

      loc0 = (long)RDB[loc0 + REACTION_PTR_MAJORANT_XS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
    }

  /* Allocate memory */

  if (mpiid == 0)
    xs0 = (double *)Mem(MEM_ALLOC, mpitasks*ne, sizeof(double));
  else
    xs0 = NULL;

  xs = (double *)Mem(MEM_ALLOC, mpitasks*ne, sizeof(double));

  /* Loop over data */

  for (n = 0; n < ne; n++)
    xs[mpiid*ne + n] = RDB[loc0 + n];

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Reduce data */

  MPITransfer(xs, xs0, ne*mpitasks, 0, MPI_METH_RED);

  /* Calculate maximum */

  if (mpiid == 0)
    for (m = 1; m < mpitasks; m++)
      for (n = 0; n < ne; n++)
        if (xs0[m*ne + n] > xs[n])
          xs[n] = xs0[m*ne + n];

  /* Synchronise */

  MPI_Barrier(my_comm);

  /* Broadcast data */

  MPITransfer(xs, NULL, ne, 0, MPI_METH_BC);

  /* Put data */

  for (n = 0; n < ne; n++)
    WDB[loc0 + n] = xs[mpiid*ne + n];

  /* Free allocated memory */

  if (mpiid == 0)
    Mem(MEM_FREE, xs0);

  Mem(MEM_FREE, xs);

  fprintf(outp, "OK.\n\n");

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
}

/*****************************************************************************/
