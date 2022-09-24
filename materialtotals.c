/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : materialtotals.c                               */
/*                                                                           */
/* Created:       2011/01/02 (JLe)                                           */
/* Last modified: 2018/11/13 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates material-wise total cross sections                */
/*                                                                           */
/* Comments: - Rutiinia muutettu radikaalisti 2.11.2011 (2.0.37)             */
/*                                                                           */
/*           - Toi väliaikaisen muistin varaaminen tehdään nyt vähän         */
/*             hölmösti, mutta rutiini on tolla tavalla selkeämpi.           */
/*                                                                           */
/*           - Tosta poistettiin #ifdef-lauseella kommentoitu vaihtoehtoinen */
/*             tapa 14.3.2017 / 1.1.31.                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MaterialTotals:"

/* Use local function to simplify OpenMP implementation */

void MaterialTotals0(long);

/*****************************************************************************/

void MaterialTotals()
{
  long mat;

  /* Check reconstruction */

  if (((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO) &&
      ((long)RDB[DATA_OPTI_MG_MODE] == NO))
    return;

  fprintf(outp, "Calculating macroscopic cross sections:\n\n");

  /***************************************************************************/

  /***** Parallelization by material *****************************************/

  /* Reset thread numbers */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check thread number */

      WDB[mat + MATERIAL_OMP_ID] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Start parallel timer */

  StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private (mat)
#endif
  {
    /* Print */

    PrintProgress(0, 0);

    /* Loop over materials */

    mat = (long)RDB[DATA_PTR_M0];
    while (mat > VALID_PTR)
      {
        /* Test parallel id's */

        if (MyParallelMat(mat, YES) == YES)
          {
            /* Process */

            MaterialTotals0(mat);

            /* Print */

            PrintProgress(mat, 1);
          }

        /* Next material */

        mat = NextItem(mat);
      }
  }

  /* Stop parallel timer */

  StopTimer(TIMER_OMP_PARA);

  /* Print */

  PrintProgress(0, 100);

  /***********************************************************************/
}

/*****************************************************************************/

/*****************************************************************************/

void MaterialTotals0(long mat)
{
  long n, mode, rea, loc0, loc1, loc2, loc3, ptr, pte, sz, i0, m, ne, erg, id;
  long nemax, nuc, nr, idx, fortms, arr;
  double Emin, Emax, adens, *xs, *tot, *nubar, mult, *fisse;

  /* Check div type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    return;

  if ((long)RDB[mat + MATERIAL_MPI_ID] > -1)
    if ((long)RDB[mat + MATERIAL_MPI_ID] != mpiid)
      Die(FUNCTION_NAME, "WTF?");

  /* Get OpenMP id */

  id = (long)RDB[mat + MATERIAL_OMP_ID];

  /***************************************************************************/

  /***** Multi-group total cross sections ************************************/

  /* Check mode */

  if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
    {
      /* Get number of energy groups */

      ne = (long)RDB[DATA_COARSE_MG_NE];
      CheckValue(FUNCTION_NAME, "ne", "", ne, 10, 50000);

      /* Pointer to total xs (tï¿½hï¿½n voi kosahtaa jos on pelkkï¿½ï¿½ gammadataa) */

      loc0 = (long)RDB[mat + MATERIAL_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(loc0 1)", DATA_ARRAY, loc0);

      /* Pointer to data */

      loc0 = (long)RDB[loc0 + REACTION_PTR_MGXS];
      CheckPointer(FUNCTION_NAME, "(loc0 2)", DATA_ARRAY, loc0);

      /* Reset data */

      memset(&WDB[loc0], 0.0, ne*sizeof(double));

      /* Pointer to total list */

      loc1 = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while (NextReaction(loc1, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
        {
          /* Check reaction pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get pointer to xs data */

          loc2 = (long)RDB[rea + REACTION_PTR_MGXS];
          CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

          /* Get pointer for checks */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check poisons */

          if (((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES) &&
              ((long)RDB[nuc + NUCLIDE_ZAI] == 541350))
            if (adens > 0.0)
              Die(FUNCTION_NAME, "Xe-135 density is not zero");

          if (((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES) &&
              ((long)RDB[nuc + NUCLIDE_ZAI] == 621490))
            if (adens > 0.0)
              Die(FUNCTION_NAME, "Sm-149 density is not zero");

          /* Set density to zero if nuclide is not to be included in material */
          /* totals and majorant */

          if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_CRIT_ITER) &&
              ((long)RDB[mat + MATERIAL_PTR_ITER_ISO_LIST] > VALID_PTR))
            adens = 0.0;

          /* Set density to zero if nuclide is not to be included in material */
          /* totals and majorant due to atomic density being supplied         */
          /* separately */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS)
            {
              /* Loop through data interfaces to check if this material is */
              /* part of the list of materials for which this nuclides     */
              /* atomic density is supplied separately */

              if ((arr = (long)RDB[mat + MATERIAL_PTR_DATAIFC_ARR]) > VALID_PTR)
                {
                  /* Loop over data interfaces for this material */

                  while ((ptr = (long)RDB[arr]) > VALID_PTR)
                    {
                      /* ptr is a DATAIFC block, check nuclide */

                      if ((long)RDB[ptr + DATAIFC_PTR_NUCLIDE] == nuc)
                        {
                          /* Exclude from material total */

                          adens = 0.0;
                          /*
                          fprintf(outp, "Excluding %s from precalculated material total XS for %s\n",
                                  GetText(nuc + NUCLIDE_PTR_NAME),
                                  GetText(mat + MATERIAL_PTR_NAME));
                          */
                          /* No need to check further */

                          break;
                        }
                      /* Next data interface */

                      arr++;
                    }
                }
            }

          /* Add to total */

          for (n = 0; n < ne; n++)
            WDB[loc0 + n] = RDB[loc0 + n] + adens*RDB[loc2 + n];
        }
    }

  /***************************************************************************/

  /***** Continuous-energy cross sections ************************************/

  /* Check reconstruction */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO)
    return;

  /* Get maximum unionized grid size */

  nemax = -1;

  /* Check neutron and photon grids */

  if ((erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) > VALID_PTR)
    nemax = (long)RDB[erg + ENERGY_GRID_NE];
  if ((erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID]) > VALID_PTR)
    if ((long)RDB[erg + ENERGY_GRID_NE] > nemax)
      nemax = (long)RDB[erg + ENERGY_GRID_NE];

  /* Allocate memory for temporary cross section array */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == NO)
    xs = WorkArray(DATA_PTR_WORK_PRIVA_GRID1, PRIVA_ARRAY, nemax, id);
  else
    xs = NULL;

  /* Allocate memory for temporary nubar array */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
    nubar = WorkArray(DATA_PTR_WORK_PRIVA_GRID2, PRIVA_ARRAY, nemax, id);
  else
    nubar = NULL;

  /* Allocate memory for temporary fission energy deposition array */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
    fisse = WorkArray(DATA_PTR_WORK_PRIVA_GRID3, PRIVA_ARRAY, nemax, id);
  else
    fisse = NULL;

  /* Set number of modes */

  nr = 17;

  /* Loop over lists */

  for (n = 0; n < nr; n++)
    {
      /* Avoid compiler warning */

      mode = -1;
      loc0 = -1;
      erg = -1;
      fortms = 0;

      /* Get mode */

      if (n == 0)
        {
          mode = MATERIAL_PTR_TMP_MAJORANT_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];

          /* This will be included for CE TMS */

          fortms = 1;
        }
      else if (n == 1)
        {
          mode = MATERIAL_PTR_TOT_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_TOTXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 2)
        {
          mode = MATERIAL_PTR_ELA_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_ELAXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 3)
        {
          mode = MATERIAL_PTR_ABS_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_ABSXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 4)
        {
          mode = MATERIAL_PTR_FISS_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_FISSXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 5)
        {
          mode = MATERIAL_PTR_HEATT_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_HEATTXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 6)
        {
          mode = MATERIAL_PTR_PHOTP_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_PHOTPXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 7)
        {
          mode = MATERIAL_PTR_PROTP_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_PROTPXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 8)
        {
          mode = MATERIAL_PTR_DEUTP_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_DEUTPXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 9)
        {
          mode = MATERIAL_PTR_TRITP_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_TRITPXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 10)
        {
          mode = MATERIAL_PTR_HE3P_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_HE3PXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 11)
        {
          mode = MATERIAL_PTR_HE4P_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_HE4PXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 12)
        {
          mode = MATERIAL_PTR_INLP_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_INLPXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 13)
        {
          mode = MATERIAL_PTR_FISS_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_FISSE];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 14)
        {
          mode = MATERIAL_PTR_FISS_REA_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_NSF];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
        }
      else if (n == 15)
        {
          mode = MATERIAL_PTR_PHOT_TOT_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID];

          /* This will be included for CE TMS */

          fortms = 1;
        }
      else if (n == 16)
        {
          mode = MATERIAL_PTR_PHOT_HEAT_LIST;
          loc0 = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS];
          erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID];

          /* This will be included for CE TMS */

          fortms = 1;
        }
      else
        Die(FUNCTION_NAME, "Overflow");

      /* Only include TMP-majorant and photon cross sections for CE TMS */

      if (((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE) &&
          (fortms == 0))
        continue;

      /* Cycle loop if nuclide has no reactions of this type */

      if ((long)RDB[mat + mode] < VALID_PTR)
        continue;

      /* Check reaction and unionized grid pointer */

      CheckPointer(FUNCTION_NAME, "(loc0 3)", DATA_ARRAY, loc0);
      CheckPointer(FUNCTION_NAME, "(erg1)", DATA_ARRAY, erg);

      /* Get number of energy points */

      ne = (long)RDB[erg + ENERGY_GRID_NE];

      /* Pointer to cross section data */

      ptr = (long)RDB[loc0 + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      tot = &WDB[ptr];

      /* Reset total array */

      memset(tot, 0.0, ne*sizeof(double));

      /* Pointer to list data (pointteri tarkistetaan tuolla ylempï¿½nï¿½) */

      loc1 = (long)RDB[mat + mode];

      /* Reset reaction pointer (rewind list) */

      rea = -1;

      /* Loop over reactions */

      while (NextReaction(loc1, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
        {
          /* Check reaction pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get pointer to nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check type (heating and photon production are special) */

          if ((mode != MATERIAL_PTR_HEATT_REA_LIST) &&
              (mode != MATERIAL_PTR_PHOTP_REA_LIST) &&
              (mode != MATERIAL_PTR_PROTP_REA_LIST) &&
              (mode != MATERIAL_PTR_DEUTP_REA_LIST) &&
              (mode != MATERIAL_PTR_TRITP_REA_LIST) &&
              (mode != MATERIAL_PTR_HE3P_REA_LIST) &&
              (mode != MATERIAL_PTR_HE4P_REA_LIST) &&
              (mode != MATERIAL_PTR_PHOT_HEAT_LIST))
            if (((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_SUM) &&
                ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_PARTIAL))
              Die(FUNCTION_NAME, "Invalid reaction type");

          /* Check poisons */

          if (((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES) &&
              ((long)RDB[nuc + NUCLIDE_ZAI] == 541350))
            if (adens > 0.0)
              Die(FUNCTION_NAME, "Xe-135 density is not zero");

          if (((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES) &&
              ((long)RDB[nuc + NUCLIDE_ZAI] == 621490))
            if (adens > 0.0)
              Die(FUNCTION_NAME, "Sm-149 density is not zero");

          /* Set density to zero if nuclide is not to be included in material */
          /* totals and majorant */

          if (((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_CRIT_ITER) &&
              ((long)RDB[mat + MATERIAL_PTR_ITER_ISO_LIST] > VALID_PTR))
            adens = 0.0;

          /* Set density to zero if nuclide is not to be included in material */
          /* totals and majorant due to atomic density being supplied         */
          /* separately */

          if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS)
            {
              /* Loop through data interfaces to check if this material is */
              /* part of the list of materials for which this nuclides     */
              /* atomic density is supplied separately */

              if ((arr = (long)RDB[mat + MATERIAL_PTR_DATAIFC_ARR]) > VALID_PTR)
                {
                  /* Loop over data interfaces for this material */

                  while ((ptr = (long)RDB[arr]) > VALID_PTR)
                    {
                      /* ptr is a DATAIFC block, check nuclide */

                      if ((long)RDB[ptr + DATAIFC_PTR_NUCLIDE] == nuc)
                        {
                          /* Exclude from material total */

                          adens = 0.0;
                          /*
                          fprintf(outp, "Excluding %s from precalculated material total XS for %s\n",
                                  GetText(nuc + NUCLIDE_PTR_NAME),
                                  GetText(mat + MATERIAL_PTR_NAME));
                          */
                          /* No need to check further */

                          break;
                        }
                      /* Next data interface */

                      arr++;
                    }
                }
            }

          /* Get multiplier */

          mult = RDB[rea + REACTION_WGT_F];
          CheckValue(FUNCTION_NAME, "mult", "", mult, 1.0, 5.0);

          /* Check if cross section is fission neutron or energy production, */
          /* or multiplying scattering. */

          if (loc0 == (long)RDB[mat + MATERIAL_PTR_FISSE])
            {
              /* Replace multiplier with fission energy */
              /* (JLe: muutettu 9.3.2016 / 2.1.26) */

              if (mult != 1.0)
                Die(FUNCTION_NAME, "Error in multiplier");
              else
                {
                  /*
                    mult = RDB[rea + REACTION_Q];
                    CheckValue(FUNCTION_NAME, "mult", "", mult, 150.0, 280.0);
                    mult = mult*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ;
                  */
                  /*
                  mult = RDB[nuc + NUCLIDE_FISSE];
                  CheckValue(FUNCTION_NAME, "mult", "", mult/MEV, 150.0,
                             280.0);
                  */
                  InterpolateFissE(fisse, rea, id);
                }
            }
          else if (loc0 == (long)RDB[mat + MATERIAL_PTR_NSF])
            {
              /* Check multiplier */

              if (mult != 1.0)
                      Die(FUNCTION_NAME, "Error in multiplier");

              /* Get nubar array */

              InterpolateNubar(nubar, rea);
            }
          else if (mode == MATERIAL_PTR_INLP_REA_LIST)
            {
              /* Check value and adjust */

              if (mult < 2.0)
                Die(FUNCTION_NAME, "Invalid multiplication");
              else
                mult = mult - 1.0;
            }

          /* Pointer to data */

          loc2 = (long)RDB[rea + REACTION_PTR_XS];

          /* Number of points and index to first point */

          sz = (long)RDB[rea + REACTION_XS_NE];
          i0 = (long)RDB[rea + REACTION_XS_I0];

          /* Get pointer to nuclide energy grid */

          pte = (long)RDB[rea + REACTION_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

          /* Check if nuclide uses unionized grid */

          if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == YES)
            {
              /* Check that pointers match */

              if (pte != erg)
                Die(FUNCTION_NAME, "Mismatch in grid pointer");

              /* Check if nubar data */

              if (loc0 == (long)RDB[mat + MATERIAL_PTR_NSF])
                {
                  /* Loop over non-zero values and add to sum */

                  for (m = 0; m < sz; m++)
                    tot[i0 + m] = tot[i0 + m] +
                      mult*adens*RDB[loc2 + m]*nubar[i0 + m];
                }

              /* Check if fission energy deposition data */

              else if (loc0 == (long)RDB[mat + MATERIAL_PTR_FISSE])
                {
                  /* Loop over non-zero values and add to sum */

                  for (m = 0; m < sz; m++)
                    tot[i0 + m] = tot[i0 + m] +
                      mult*adens*RDB[loc2 + m]*fisse[i0 + m];
                }
              else
                {
                  /* Loop over non-zero values and add to sum */

                  for (m = 0; m < sz; m++)
                    tot[i0 + m] = tot[i0 + m] + mult*adens*RDB[loc2 + m];
                }
            }
          else
            {
              /* Pointer to grid data */

              loc3 = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(loc3)", DATA_ARRAY, loc3);

              pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

              /* Histogrammityyppisia majoranttivaikutusaloja ei voi
                 interpoloida uudelle gridille samalla tavalla kuin
                 muita vaikutusaloja --> tï¿½ï¿½ pitï¿½ï¿½ tehdï¿½ erikseen */
              if( mode != MATERIAL_PTR_TMP_MAJORANT_LIST){

                /* Reconstruct data to temporary array */

                if (loc0 == (long)RDB[mat + MATERIAL_PTR_HEATTXS])
                  InterpolateData(&RDB[loc3], xs, ne, &RDB[pte + i0],
                                  &RDB[loc2], sz, 0, NULL, NULL, YES);
                else
                  InterpolateData(&RDB[loc3], xs, ne, &RDB[pte + i0],
                                  &RDB[loc2], sz, 0, NULL, NULL, NO);

                /* Add points (NOTE: Tï¿½ï¿½ kï¿½ydï¿½ï¿½n nyt nollasta asti lï¿½pi, */
                /* toi InterpolateData() voisi palauttaa tiedon siitï¿½, */
                /* mikï¿½ ensimmï¿½inen nollasta poikkeava piste on. */
                /* NOTE2: se palauttaa, sitï¿½ ei oo vaan otettu tossa */
                /* seuraavassa vielï¿½ huomioon */

                /* Check if nubar data */

                if (loc0 == (long)RDB[mat + MATERIAL_PTR_NSF])
                  {
                  for (m = 0; m < ne; m++)
                    tot[m] = tot[m] + mult*adens*xs[m]*nubar[m];
                  }

                /* Check if fission energy deposition data */

                else if (loc0 == (long)RDB[mat + MATERIAL_PTR_FISSE])
                  {
                  for (m = 0; m < ne; m++)
                    tot[m] = tot[m] + mult*adens*xs[m]*fisse[m];
                  }
                else
                  {
                    for (m = 0; m < ne; m++)
                      tot[m] = tot[m] + mult*adens*xs[m];
                  }
              }

              /* Lisataan TMP-majorantti, huomioidaan histogrammityyppi */

              else{
                idx = 0;

                for (m = 0; m < ne; m++){

                  /* Find correct interval from nuclide reaction grid */

                  while(RDB[pte + idx + 1] < RDB[loc3 + m + 1])
                    idx++;

                  /* Add contribution */

                  tot[m] = tot[m] + mult*adens*RDB[loc2 + idx];

                }
              }
            }
        }
    }

  /* Free allocated memory */
  /*
  if (xs != NULL)
    Mem(MEM_FREE, xs);
  */
}

/*****************************************************************************/
