/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sumdivcompositions.c                           */
/*                                                                           */
/* Created:       2012/05/13 (JLe)                                           */
/* Last modified: 2019/03/27 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates average compositions for divided burnable         */
/*              materials.                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SumDivCompositions:"

/*****************************************************************************/

void SumDivCompositions()
{
  long mat0, mat1, iso0, iso1, sz, n;
  double vol, tot, *adens, *adens0, bu;

  /* Avoid compiler warning */

  adens = NULL;
  adens0 = NULL;
  bu = 0.0;
  sz = 0;
  n = 0;

  /***************************************************************************/

  /***** Calculate average compositions **************************************/

  /* Reset compositions */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Check if domain decomposition is in use */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
          ((long)RDB[mat0 + MATERIAL_MPI_ID] > -1) &&
          ((long)RDB[mat0 + MATERIAL_MPI_ID] != mpiid))
        {
          /* Pointer to next */

          mat0 = NextItem(mat0);

          /* Cycle loop */

          continue;
        }

      /* Check div type */

      if ((long)RDB[mat0 + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        {
          /* Loop over composition and reset */

          iso0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso0 > VALID_PTR)
            {
              /* Reset density */

              WDB[iso0 + COMPOSITION_ADENS] = 0.0;

              /* Next nuclides */

              iso0 = NextItem(iso0);
            }

          /* Reset burnup */

          WDB[mat0 + MATERIAL_BURNUP] = 0.0;
        }

      /* Next material */

      mat0 = NextItem(mat0);
    }

  /* Put compositions */

  mat1 = (long)RDB[DATA_PTR_M0];
  while (mat1 > VALID_PTR)
    {
      /* Check if domain decomposition is in use */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
          ((long)RDB[mat1 + MATERIAL_MPI_ID] > -1) &&
          ((long)RDB[mat1 + MATERIAL_MPI_ID] != mpiid))
        {
          /* Pointer to next */

          mat1 = NextItem(mat1);

          /* Cycle loop */

          continue;
        }

      /* Check burn-flag and pointer to parent */

      if (((long)RDB[mat1 + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
          ((mat0 = (long)RDB[mat1 + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR))
        {
          /* Get total volume */

          vol = RDB[mat0 + MATERIAL_VOLUME];

          /* Volume can be zero in reprocessing mode */

          if (vol == 0.0)
            {
              /* Check reprocessing mode */

              if ((long)RDB[DATA_PTR_REP0] < VALID_PTR)
                Warn(FUNCTION_NAME, "Zero volume and no reprocessing");

              /* Check second volume */

              if (RDB[mat1 + MATERIAL_VOLUME] != 0.0)
                Warn(FUNCTION_NAME, "mismatch in volumes");
            }
          else
            {
              /* Check volume */

              CheckValue(FUNCTION_NAME, "vol", "", vol, ZERO, INFTY);

              /* Loop over composition and add to sum */

              iso0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
              CheckPointer(FUNCTION_NAME, "(iso0)", DATA_ARRAY, iso0);

              iso1 = (long)RDB[mat1 + MATERIAL_PTR_COMP];
              CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);

              while (iso0 > VALID_PTR)
                {
                  /* Check second pointer */

                  CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);

                  /* Check nuclide pointers */

                  if ((long)RDB[iso0 + COMPOSITION_PTR_NUCLIDE] !=
                      (long)RDB[iso1 + COMPOSITION_PTR_NUCLIDE])
                    Die(FUNCTION_NAME, "Mismatch in composition");

                  /* Add to atomic density */

                  WDB[iso0 + COMPOSITION_ADENS] =
                    RDB[iso0 + COMPOSITION_ADENS] +
                    WDB[iso1 + COMPOSITION_ADENS]*
                    RDB[mat1 + MATERIAL_VOLUME]/vol;

                  /* Next nuclides */

                  iso0 = NextItem(iso0);
                  iso1 = NextItem(iso1);
                }

              /* Add to burnup */

              WDB[mat0 + MATERIAL_BURNUP] = RDB[mat0 + MATERIAL_BURNUP] +
                RDB[mat1 + MATERIAL_BURNUP]*RDB[mat1 + MATERIAL_VOLUME]/vol;
            }
        }

      /* Next material */

      mat1 = NextItem(mat1);
    }

  /***************************************************************************/

  /***** Collect data from other domains *************************************/

#ifdef MPI

  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Start timers */

      StartTimer(TIMER_MPI_OVERHEAD);
      StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

      /* Loop over materials */

      mat0 = (long)RDB[DATA_PTR_M0];
      while (mat0 > VALID_PTR)
        {
          /* Check div type */

          if ((long)RDB[mat0 + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
            {
              /* Pointer to composition */

              iso0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
              CheckPointer(FUNCTION_NAME, "(iso0)", DATA_ARRAY, iso0);

              /* Get size */

              sz = ListSize(iso0);

              /* Allocate memory for vectors */

              if (mpiid == 0)
                adens0 = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
              else
                adens0 = NULL;

              adens = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

              /* Loop over composition */

              n = 0;
              while (iso0 > VALID_PTR)
                {
                  /* Get atomic density */

                  adens[n++] = RDB[iso0 + COMPOSITION_ADENS];

                  /* Next */

                  iso0 = NextItem(iso0);
                }

              /* Synchronise */

              MPI_Barrier(my_comm);

              /* Reduce data */

              MPITransfer(adens, adens0, sz, 0, MPI_METH_RED);
              MPITransfer(&WDB[mat0 + MATERIAL_BURNUP], &bu, 1, 0,
                          MPI_METH_RED);
              /* Copy */

              if (mpiid == 0)
                {
                  memcpy(adens, adens0, sz*sizeof(double));
                  WDB[mat0 + MATERIAL_BURNUP] = bu;
                }

              /* Synchronise */

              MPI_Barrier(my_comm);

              /* Broadcast data to other tasks */

              MPITransfer(adens, NULL, sz, 0, MPI_METH_BC);
              MPITransfer(&WDB[mat0 + MATERIAL_BURNUP], NULL, 1, 0,
                          MPI_METH_BC);

              /* Loop over composition */

              n = 0;

              iso0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
              while (iso0 > VALID_PTR)
                {
                  /* Put data */

                  WDB[iso0 + COMPOSITION_ADENS] = adens[n++];

                  /* Next */

                  iso0 = NextItem(iso0);
                }

              /* Free allocated memory */

              if (mpiid == 0)
                Mem(MEM_FREE, adens0);

              Mem(MEM_FREE, adens);
            }
          /* Next material */

          mat0 = NextItem(mat0);
        }

      /* Stop timers */

      StopTimer(TIMER_MPI_OVERHEAD);
      StopTimer(TIMER_MPI_OVERHEAD_TOTAL);
    }

#endif

  /***************************************************************************/

  /***** Remaining stuff *****************************************************/

  /* Calculate total atomic density (Lisätty 22.5.2015 / 2.1.24, että */
  /* restart-filesta voi lukea parent-materiaalin niin että tiheys    */
  /* vastaa todellista koostumusta.) */

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Check div type */

      if ((long)RDB[mat0 + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        {
          /* Check if domain decomposition is in use */

          if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
              ((long)RDB[mat0 + MATERIAL_MPI_ID] > -1) &&
              ((long)RDB[mat0 + MATERIAL_MPI_ID] != mpiid))
            Die(FUNCTION_NAME, "Error in ID");

          /* Reset density */

          WDB[mat0 + MATERIAL_ADENS] = 0.0;

          /* Loop over composition and reset */

          iso0 = (long)RDB[mat0 + MATERIAL_PTR_COMP];
          while (iso0 > VALID_PTR)
            {
              /* Add to total */

              WDB[mat0 + MATERIAL_ADENS] = RDB[mat0 + MATERIAL_ADENS] +
                RDB[iso0 + COMPOSITION_ADENS];

              /* Next nuclides */

              iso0 = NextItem(iso0);
            }
        }

      /* Next material */

      mat0 = NextItem(mat0);
    }

  /* Calculate total burnup */

  tot = 0.0;

  mat0 = (long)RDB[DATA_PTR_M0];
  while (mat0 > VALID_PTR)
    {
      /* Add to total */

      if ((long)RDB[mat0 + MATERIAL_DIV_PTR_PARENT] < VALID_PTR)
        tot = tot + RDB[mat0 + MATERIAL_BURNUP]*
          RDB[mat0 + MATERIAL_INI_FISS_MDENS]*RDB[mat0 + MATERIAL_VOLUME];

      /* Next material */

      mat0 = NextItem(mat0);
    }

  /* Put total cumulative real burnup */

  if (RDB[DATA_INI_BURN_FMASS] > 0.0)
    WDB[DATA_BURN_CUM_REAL_BURNUP] = tot/RDB[DATA_INI_BURN_FMASS]/1000.0;

  /***************************************************************************/
}

/*****************************************************************************/
