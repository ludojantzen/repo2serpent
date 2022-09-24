/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : distributematerialdata.c                       */
/*                                                                           */
/* Created:       2011/12/23 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Distributes material-wise data (reaction lists, macroscopic  */
/*              cross sections, etc.) calculated by MPI tasks.               */
/*                                                                           */
/* Comments: - Processing of non-burnable materials is not parallelized      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DistributeMaterialData:"

/*****************************************************************************/

void DistributeMaterialData()
{
#ifdef MPI

  long mat, id, ptr, sz;

  /* Check number of MPI tasks */

  if (mpitasks == 1)
    return;

  /* Check if domain decomposition is in use */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    return;

  fprintf(outp, "Waiting for results from other MPI tasks...\n");

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

  /***************************************************************************/

  /***** Collect data ********************************************************/

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check divisor type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check burn-flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /* Get MPI id */

          id = (long)RDB[mat + MATERIAL_MPI_ID];
          CheckValue(FUNCTION_NAME, "id", "", id, 0, mpitasks - 1);

          /* Get pointer to data */

          ptr = (long)RDB[mat + MATERIAL_PTR_DATA_BLOCK];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get block size */

          sz = (long)RDB[mat + MATERIAL_DATA_BLOCK_SIZE];

          /* Synchronise */

          MPI_Barrier(my_comm);

          /* Broadcast data to other tasks */

          MPITransfer(&WDB[ptr], NULL, sz, id, MPI_METH_BC);

          /* Synchronise */

          MPI_Barrier(my_comm);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
}

/*****************************************************************************/
