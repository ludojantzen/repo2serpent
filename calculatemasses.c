/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatemasses.c                              */
/*                                                                           */
/* Created:       2011/05/31 (JLe)                                           */
/* Last modified: 2019/08/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates material and fissile masses etc.                  */
/*                                                                           */
/* Comments: - NOTE: tän kanssa pitää olla varovainen, sillä palamat ym.     */
/*             tehdään alkuperäiseen fissiiliin massaan, eli niitä ei        */
/*             päivitetä ensimmäisen askeleen jälkeen. Serpent 1:ssä on      */
/*             optio sille että tuo massa muuttuu, mutta se on nyt jätetty   */
/*             tästä pois.                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateMasses:"

/*****************************************************************************/

void CalculateMasses()
{
  long mat, iso, nuc, ptr;
  double fmass, mdens, adens, fm0[4], fm[4], sum;

  /* Avoid compiler warning */

  fm[0] = 0.0;
  fm0[0] = 0.0;

  /***************************************************************************/

  /***** Calculate masses ****************************************************/

  /* Reset total masses */

  WDB[DATA_TOT_FMASS] = 0.0;
  WDB[DATA_TOT_BURN_FMASS] = 0.0;

  /* Reset initial masses */

  if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
    {
      /* Reset totals */

      WDB[DATA_INI_FMASS] = 0.0;
      WDB[DATA_INI_BURN_FMASS] = 0.0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Reset material-wise masses */

          WDB[mat + MATERIAL_INI_FMASS] = 0.0;

          /* Next */

          mat = NextItem(mat);
        }
    }

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if domain decomposition is in use */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Get material density */

      mdens = RDB[mat + MATERIAL_MDENS];

      /* Check value */

      if (mdens < ZERO)
        {
          /* Check if used in continuous reprocessing */

          if ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] > VALID_PTR)
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }
          else
            Die(FUNCTION_NAME, "material %s mass density too low (%E)",
                GetText(mat + MATERIAL_PTR_NAME), mdens);
        }

      /* Calculate material mass */

      WDB[mat + MATERIAL_MASS] = RDB[mat + MATERIAL_VOLUME]*mdens;

      /* Reset fissile mass */

      fmass = 0.0;

      /* Loop over composition and calculate mass of fissile isotopes */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Get atomic density and pointer to nuclide */

          adens = RDB[iso + COMPOSITION_ADENS];
          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check fissile flag and add to mass */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
            fmass = fmass + adens*RDB[nuc + NUCLIDE_AW]*
              RDB[mat + MATERIAL_VOLUME]/N_AVOGADRO;

          /* Next isotope */

          iso = NextItem(iso);
        }

      /* Check fissile flag and burnup step */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
          !((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
        {
          /* Exlcude divided materials */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
            {
              /* Add to total fissile mass */

              WDB[DATA_TOT_FMASS] = RDB[DATA_TOT_FMASS] + 1E-3*fmass;

              if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
                WDB[DATA_INI_FMASS] = RDB[DATA_INI_FMASS] + 1E-3*fmass;

              /* Put material-wise mass */

              if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
                {
                  /* Put material-wise mass */

                  WDB[mat + MATERIAL_INI_FMASS] = 1E-3*fmass;

                  /* Add to parent */

                  if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
                      > VALID_PTR)
                    WDB[ptr + MATERIAL_INI_FMASS] =
                      RDB[ptr + MATERIAL_INI_FMASS] + 1E-3*fmass;
                }

              /* Check burn-flag */

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                {
                  /* Add to burnable fissile mass */

                  WDB[DATA_TOT_BURN_FMASS] = RDB[DATA_TOT_BURN_FMASS] +
                    1E-3*fmass;

                  if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
                    WDB[DATA_INI_BURN_FMASS] = RDB[DATA_INI_BURN_FMASS] +
                      1E-3*fmass;
                }
            }
        }

      /* Check burn flag and mass */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
          (RDB[mat + MATERIAL_MASS] <= 0.0) &&
          ((long)RDB[mat + MATERIAL_VOL_COUNT] > 0))
        Error(mat, "Volume or mass of material \"%s\" must be given",
              GetText(mat + MATERIAL_PTR_NAME));

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** FIMA-calculation ****************************************************/

  /* Reset number of actinide atoms */

  if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
    WDB[DATA_FIMA_ACT0] = 0.0;

  WDB[DATA_FIMA_ACT] = 0.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check if domain decomposition is in use */

      if (((long)RDB[DATA_DD_DECOMPOSE] == YES) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] > -1) &&
          ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check burn and fissile flag */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) ||
          !((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Reset sum */

      sum = 0.0;

      /* Loop over composition and calculate number of atoms */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Get atomic density */

          adens = RDB[iso + COMPOSITION_ADENS];

          /* Get pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Include actinides and decay series */

          if ((long)RDB[nuc + NUCLIDE_ZAI] > 810000)
            sum = sum + adens;

          /* Next */

          iso = NextItem(iso);
        }

      /* Put value */

      if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
        {
          WDB[mat + MATERIAL_FIMA_ADENS0] = sum;
          WDB[DATA_FIMA_ACT0] = RDB[DATA_FIMA_ACT0]
            + sum*RDB[mat + MATERIAL_VOLUME]/BARN;
        }
      else
        {
          WDB[mat + MATERIAL_FIMA_ADENS] = sum;
          WDB[DATA_FIMA_ACT] = RDB[DATA_FIMA_ACT]
            + sum*RDB[mat + MATERIAL_VOLUME]/BARN;
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Collect results from other domains **********************************/

#ifdef MPI

  /* Check domain decomposition */

  if ((long)RDB[DATA_DD_DECOMPOSE] == YES)
    {
      /* Start timers */

      StartTimer(TIMER_MPI_OVERHEAD);
      StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check div type */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
            {
              /* Next material */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

           /* Set data */

          fm[0] = RDB[mat + MATERIAL_INI_FMASS];
          fm[1] = RDB[mat + MATERIAL_FIMA_ADENS0];
          fm[2] = RDB[mat + MATERIAL_FIMA_ADENS];

          fm0[0] = 0.0;
          fm0[1] = 0.0;
          fm0[2] = 0.0;

          /* Synchronise */

          MPI_Barrier(my_comm);

          /* Reduce data */

          MPITransfer(fm, fm0, 3, 0, MPI_METH_RED);

          /* Synchronise */

          MPI_Barrier(my_comm);

          /* Broadcast data to other tasks */

          MPITransfer(fm0, NULL, 3, 0, MPI_METH_BC);

          /* Put data */

          if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
            {
              WDB[mat + MATERIAL_INI_FMASS] = fm0[0];
              WDB[mat + MATERIAL_FIMA_ADENS0] = fm0[1];
            }

          WDB[mat + MATERIAL_FIMA_ADENS] = fm0[2];

          /* Next material */

          mat = NextItem(mat);
        }

      /* Set data */

      fm[0] = RDB[DATA_TOT_FMASS];
      fm[1] = RDB[DATA_TOT_BURN_FMASS];
      fm[2] = RDB[DATA_FIMA_ACT];
      fm[3] = RDB[DATA_FIMA_ACT0];

      fm0[0] = 0.0;
      fm0[1] = 0.0;
      fm0[2] = 0.0;
      fm0[3] = 0.0;

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Reduce data */

      MPITransfer(fm, fm0, 4, 0, MPI_METH_RED);

      /* Synchronise */

      MPI_Barrier(my_comm);

      /* Broadcast data to other tasks */

      MPITransfer(fm0, NULL, 4, 0, MPI_METH_BC);

      /* Copy data from vector */

      WDB[DATA_TOT_FMASS] = fm0[0];
      WDB[DATA_TOT_BURN_FMASS] = fm0[1];
      WDB[DATA_FIMA_ACT] = fm0[2];

      if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
        WDB[DATA_FIMA_ACT0] = fm0[3];

      /* Put initial masses */

      if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
        {
          WDB[DATA_INI_FMASS] = WDB[DATA_TOT_FMASS];
          WDB[DATA_INI_BURN_FMASS] = WDB[DATA_TOT_BURN_FMASS];
        }

      /* Stop timers */

      StopTimer(TIMER_MPI_OVERHEAD);
      StopTimer(TIMER_MPI_OVERHEAD_TOTAL);
    }

#endif

  /***************************************************************************/

  /* Reset initial mass flag */

  WDB[DATA_BURN_CALC_INI_MASS] = (double)NO;
}

/*****************************************************************************/
