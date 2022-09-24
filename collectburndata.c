/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectburndata.c                              */
/*                                                                           */
/* Created:       2011/12/23 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Collects material compositions from MPI parallelized burnup  */
/*              calculation                                                  */
/*                                                                           */
/* Comments: - This subroutine collects and distributes the compositions of  */
/*             burnable materials after parallelized depletion solver.       */
/*             Reaction lists are formed and material-wise total cross       */
/*             sections calculated in parallel in subroutines called from    */
/*             PrepareTransportCycle().                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectBurnData:"

/*****************************************************************************/

void CollectBurnData()
{
#ifdef MPI

  long mat, iso, nuc, i, id, nnuc, sz;
  double *N;

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

  /***** Check that material compositions are sorted by ZAI ******************/

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn-flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /* First nuclide should be lost (ZAI = -1) */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get ZAI */

          if ((i = (long)RDB[nuc + NUCLIDE_ZAI]) != -1)
            Die(FUNCTION_NAME, "First ZAI should be -1");

          /* Loop over remaining */

          iso = NextItem(iso);
          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Compare */

              if ((long)RDB[nuc + NUCLIDE_ZAI] > i)
                i = (long)RDB[nuc + NUCLIDE_ZAI];
              else
                Die(FUNCTION_NAME, "Composition must be in ascending order");

              /* Next */

              iso = NextItem(iso);
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Collect data ********************************************************/

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn-flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
        {
          /* Get MPI id */

          id = (long)RDB[mat + MATERIAL_MPI_ID];
          CheckValue(FUNCTION_NAME, "id", "", id, 0, mpitasks - 1);

          /* Pointer to composition */

          iso = (long)RDB[mat + MATERIAL_PTR_COMP];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

          /* Get number of nuclides size */

          nnuc = ListSize(iso);

          /* Get array size (nuclides + density + burnup) */

          sz = nnuc + 2;

          /* Allocate memory for temporary array */

          N = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

          /* Check id and read composition */

          if (id == mpiid)
            {
              /* Loop over composition */

              i = 0;

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Put density */

                  if (i < nnuc)
                    N[i++] = RDB[iso + COMPOSITION_ADENS];
                  else
                    Die(FUNCTION_NAME, "Indexing error");

                  /* Next */

                  iso = NextItem(iso);
                }

              /* Check count */

              if (i != nnuc)
                Die(FUNCTION_NAME, "Indexing error");

              /* Put material density and burnup */

              N[i++] = RDB[mat + MATERIAL_ADENS];
              N[i++] = RDB[mat + MATERIAL_BURNUP];

              /* Check count */

              if (i != sz)
                Die(FUNCTION_NAME, "Indexing error");
            }

          /* Synchronise */

          MPI_Barrier(my_comm);

          /* Broadcast data to other tasks */

          MPITransfer(N, NULL, sz, id, MPI_METH_BC);

          /* Synchronise */

          MPI_Barrier(my_comm);

          /* Write data to compositions */

          if (id != mpiid)
            {
              /* Loop over composition */

              i = 0;

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Put density */

                  if (i < nnuc)
                    WDB[iso + COMPOSITION_ADENS] = N[i++];
                  else
                    Die(FUNCTION_NAME, "Indexing error");

                  /* Next */

                  iso = NextItem(iso);
                }

              /* Check count */

              if (i != nnuc)
                Die(FUNCTION_NAME, "Indexing error");

              /* Get material density and burnup */

              WDB[mat + MATERIAL_ADENS] = N[i++];
              WDB[mat + MATERIAL_BURNUP] = N[i++];

              /* Check count */

              if (i != sz)
                Die(FUNCTION_NAME, "Indexing error");
            }

          /* Free temporary array */

          Mem(MEM_FREE, N);
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
