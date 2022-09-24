/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : distributeddmaxadens.c                         */
/*                                                                           */
/* Created:       2018/03/14 (JLe)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Broadcasts maximum atomic densities from reaction list       */
/*              structures to other MPI tasks and calculates absolute        */
/*              maxima.                                                      */
/*                                                                           */
/* Comments: - Used only with domain decomposition in optimization mode 1.   */
/*             When macroscopic cross sections are pre-calculated or         */
/*             MG-majorant is in use, the majorants are combined in          */
/*             distributeddmajorants.c                                       */
/*                                                                           */
/*           - Tonne voi joutua lisäämään vielä ures-listankin.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DistributeDDMaxAdens:"

/*****************************************************************************/

void DistributeDDMaxAdens()
{

#ifdef MPI

  long mat, loc0, loc1, sz, n, m;
  double *adens, *adens0;

  /* Check number of tasks */

  if (mpitasks == 1)
    return;

  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == NO)
    return;

  /* Check if macroscopic xs are generated or multi-group majorant is used */

  if (((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES) ||
      ((long)RDB[DATA_OPTI_MG_MODE] == YES))
    return;

  fprintf(outp, "Sharing maximum densities to other domains...\n");

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check division */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
        {
          /* Skip material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check that material is handled in all domains */

      if ((long)RDB[mat + MATERIAL_MPI_ID] > -1)
        Die(FUNCTION_NAME, "Not included in domain");

      /* Get pointer to total reaction list (tähän ehkä myös ures-lista?) */

      loc0 = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get list size */

      sz = ListSize(loc1);
      CheckValue(FUNCTION_NAME, "sz", "", sz, 1, 10000);

      /* Allocate memory */

      if (mpiid == 0)
        adens0 = (double *)Mem(MEM_ALLOC, mpitasks*sz, sizeof(double));
      else
        adens0 = NULL;

      adens = (double *)Mem(MEM_ALLOC, mpitasks*sz, sizeof(double));

      /* Loop over data */

      n = 0;
      while (loc1 > VALID_PTR)
        {
          /* Read data into table */

          adens[mpiid*sz + n++] = RDB[loc1 + RLS_DATA_MAX_ADENS];

          /* Next reaction */

          loc1 = NextItem(loc1);
        }

      /* Check count */

      if (n != sz)
        Die(FUNCTION_NAME, "Error in count");

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Reduce data */

      MPITransfer(adens, adens0, sz*mpitasks, 0, MPI_METH_RED);

      /* Calculate maximum data */

      if (mpiid == 0)
        for (m = 1; m < mpitasks; m++)
          for (n = 0; n < sz; n++)
            if (adens0[m*sz + n] > adens[n])
              adens[n] = adens0[m*sz + n];

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Broadcast data */

      MPITransfer(adens, NULL, sz, 0, MPI_METH_BC);

      /* Put data */

      n = 0;

      loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
      while (loc1 > VALID_PTR)
        {
          /* Read data into table */

          WDB[loc1 + RLS_DATA_MAX_ADENS] = adens[n++];

          /* Next reaction */

          loc1 = NextItem(loc1);
        }

      /* Free allocated memory */

      if (mpiid == 0)
        Mem(MEM_FREE, adens0);

      Mem(MEM_FREE, adens);

      /* Next material */

      mat = NextItem(mat);
    }

  fprintf(outp, "OK.\n\n");

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
}

/*****************************************************************************/
