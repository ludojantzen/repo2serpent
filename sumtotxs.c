/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sumtotxs.c                                     */
/*                                                                           */
/* Created:       2010/12/29 (JLe)                                           */
/* Last modified: 2018/12/12 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets reaction list for partials and calculates nuclide-wise  */
/*              total and total absorption cross section                     */
/*                                                                           */
/* Comments: - Total absorption is an array of zeros if nuclide has no       */
/*             absorption channels                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SumTotXS:"

/*****************************************************************************/

void SumTotXS(long nuc)
{
  long rea0, erg, sz, loc0, ptr, pte, rea, i0, ne, i, type, ty, mode;
  long ptp, n, add, fiss, mt, sab;
  double max, *E, *xs, tmp[2];
  const double *E0, *xs0;

  /* Get type */

  type = (long)RDB[nuc + NUCLIDE_TYPE];

  /* Set flags */

  fiss = ((long)RDB[DATA_EDEP_MODE] > EDEP_MODE_MT458) &&
    ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE);

  sab = (long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA;

  /* Loop (total photon production is taken from total block */
  /* in readacefile.c) */

  for (n = 0; n < 3; n++)
    {
      /* Avoid compiler warning */

      mode = -1;

      /* Get mode */

      if (n == 0)
        mode = NUCLIDE_PTR_TOTXS;
      else if (n == 1)
        mode = NUCLIDE_PTR_SUM_ABSXS;
      else if (n == 2)
        mode = NUCLIDE_PTR_NFXS;
      else if (n == 3)
        mode = NUCLIDE_PTR_PHOTPRODXS;
      else
        Die(FUNCTION_NAME, "Overflow");

      /* Only total included for photon nuclides */

      if ((type == NUCLIDE_TYPE_PHOTON) && (n > 0))
        break;

      /* Non-fission cross section for non-fission KERMA calculation only */
      /* for fissionable free nuclides in edep > 1 */

      if ((mode == NUCLIDE_PTR_NFXS) && (sab || (!fiss)))
        continue;

      /* Create reaction structures */

      if (mode == NUCLIDE_PTR_PHOTPRODXS)
        {
          rea0 = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);
          WDB[nuc + mode] = (double)rea0;
        }
      else
        rea0 = NewItem(nuc + mode, REACTION_BLOCK_SIZE);

      /* Put nuclide pointers */

      WDB[rea0 + REACTION_PTR_NUCLIDE] = (double)nuc;

      /* Put type */

      if (mode == NUCLIDE_PTR_PHOTPRODXS || mode == NUCLIDE_PTR_NFXS)
        WDB[rea0 + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;
      else
        WDB[rea0 + REACTION_TYPE] = (double)REACTION_TYPE_SUM;

      /* Put mt */

      if (type == NUCLIDE_TYPE_PHOTON)
        WDB[rea0 + REACTION_MT] = 501.0;
      else if (mode == NUCLIDE_PTR_TOTXS)
        WDB[rea0 + REACTION_MT] = 1.0;
      else if (mode == NUCLIDE_PTR_SUM_ABSXS)
        WDB[rea0 + REACTION_MT] = 101.0;
      else if (mode == NUCLIDE_PTR_PHOTPRODXS)
        WDB[rea0 + REACTION_MT] = 202.0;
      else if (mode == NUCLIDE_PTR_NFXS)
        WDB[rea0 + REACTION_MT] = 999.0;
      else
        Die(FUNCTION_NAME, "WTF?");

      /* Put multiplier */

      WDB[rea0 + REACTION_WGT_F] = 1.0;

      /* Allocate memory for previous values */

      AllocValuePair(rea0 + REACTION_PTR_PREV_XS);

      /* Put awr */

      WDB[rea0 + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

      /* Reset maximum */

      max = 0.0;

      /* Get pointer to energy grid data */

      if ((erg = (long)RDB[nuc + NUCLIDE_PTR_EGRID]) < VALID_PTR)
        {
          /*******************************************************************/

          /***** No common grid in use ***************************************/

          /* Reset number of grid points and poiter */

          sz = 0;
          E = NULL;

          /* Loop over reaction modes and unionize grid */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check type */

              if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                {
                  /* Pointer to nuclide-wise grid */

                  erg = (long)RDB[rea + REACTION_PTR_EGRID];
                  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

                  /* Number of points and pointer to data */

                  ne = (long)RDB[erg + ENERGY_GRID_NE];

                  ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Add points to common grid */

                  E = AddPts(E, &sz, &RDB[ptr], ne);

                  /* Add extra point just below and above boundaries (this */
                  /* is needed to account for the cut-off of S(a,b) data). */

                  tmp[0] = 0.999999999*RDB[rea + REACTION_EMIN];
                  tmp[1] = 1.000000001*RDB[rea + REACTION_EMAX];

                  /* Cut values to nuclide limits (this is needed to prevent */
                  /* non-zero at limiting points) */

                  if (tmp[0] < RDB[nuc + NUCLIDE_EMIN])
                    tmp[0] = RDB[nuc + NUCLIDE_EMIN];

                  if (tmp[1] > RDB[nuc + NUCLIDE_EMAX])
                    tmp[1] = RDB[nuc + NUCLIDE_EMAX];

                  /* Add points to array */

                  E = AddPts(E, &sz, tmp, 2);
                }

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Generate grid */

          erg = MakeEnergyGrid(sz, 0, 0, -1, E, EG_INTERP_MODE_LIN);

          /* Put minimum and maximum energy */

          WDB[rea0 + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
          WDB[rea0 + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];

          /* Reset ures energy boundaries */

          WDB[rea0 + REACTION_URES_EMIN] = INFTY;
          WDB[rea0 + REACTION_URES_EMAX] = -INFTY;

          /* Put pointer to grid */

          WDB[rea0 + REACTION_PTR_EGRID] = (double)erg;

          /* Put first point and number of points */

          WDB[rea0 + REACTION_XS_I0] = 0.0;
          WDB[rea0 + REACTION_XS_NE] = (double)sz;

          /* Allocate memory for data */

          loc0 = ReallocMem(DATA_ARRAY, sz);

          /* Put pointer */

          WDB[rea0 + REACTION_PTR_XS] = (double)loc0;

          /* Allocate memory for temporary cross section array */

          xs = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

          /* Loop over reaction modes and reconstruct cross sections */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Reset add flag */

              add = NO;

              /* Get mt */

              mt = (long)RDB[rea + REACTION_MT];

              /* Check type and mode */

              if (mode == NUCLIDE_PTR_TOTXS)
                {
                  /* Add all partials to total */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                    add = YES;
                }
              else if (mode == NUCLIDE_PTR_SUM_ABSXS)
                {
                  /* Get ty */

                  ty = (long)RDB[rea + REACTION_TY];

                  /* Add partial absorptions to total absorption */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                    if (ty == 0)
                      add = YES;
                }
              else if (mode == NUCLIDE_PTR_PHOTPRODXS)
                {
                  /* Add partial photon production reactions to total */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_SPECIAL)
                    if (mt > 3000)
                      if (rea != rea0)
                        add = YES;
                }
              else if (mode == NUCLIDE_PTR_NFXS)
                {
                  if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                      && (mt < 18 || (mt > 21 && mt != 38)))
                    add = YES;
                }
              else
                Die(FUNCTION_NAME, "WTF?");

              /* Check if mode was added */

              if (add == NO)
                {
                  /* Pointer to next */

                  rea = NextItem(rea);

                  /* Cycle loop */

                  continue;
                }

              /* Pointer to reaction-wise grid */

              erg = (long)RDB[rea + REACTION_PTR_EGRID];
              CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

              /* Get number of points */

              ne = (long)RDB[erg + ENERGY_GRID_NE];

              /* Pointer to energy grid data */

              pte = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

              /* Get pointer to cross section data */

              ptr = (long)RDB[rea + REACTION_PTR_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Reconstruct cross section */

              InterpolateData(E, xs, sz, &RDB[pte], &RDB[ptr], ne, 0, NULL,
                              NULL, NO);

              /* Loop over points and add to totals */

              for (i = 0; i < sz; i++)
                {
                  /* Add value to sum */

                  WDB[loc0 + i] = RDB[loc0 + i] + xs[i];

                  /* Compare to maximum */

                  if (RDB[loc0 + i] > max)
                    max = RDB[loc0 + i];
                }

              /* Add reaction to total list */

              ptr = NewItem(rea0 + REACTION_PTR_PARTIAL_LIST,
                            RLS_DATA_BLOCK_SIZE);

              /* Allocate memory for reaction counter */

              ptp = AllocPrivateData(1, PRIVA_ARRAY);
              WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

              WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
              WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
              WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
              WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Free temporary arrays */

          Mem(MEM_FREE, E);
          Mem(MEM_FREE, xs);

          /*******************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Common energy grid available ********************************/

          /* Put minimum and maximum energy */

          WDB[rea0 + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
          WDB[rea0 + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];

          /* Reset ures energy boundaries */

          WDB[rea0 + REACTION_URES_EMIN] = INFTY;
          WDB[rea0 + REACTION_URES_EMAX] = -INFTY;

          /* Put pointer to grid */

          WDB[rea0 + REACTION_PTR_EGRID] = (double)erg;

          /* Get number of points */

          sz = (long)RDB[erg + ENERGY_GRID_NE];

          /* Put first point and number of points */

          WDB[rea0 + REACTION_XS_I0] = 0.0;
          WDB[rea0 + REACTION_XS_NE] = (double)sz;

          /* Allocate memory for data */

          loc0 = ReallocMem(DATA_ARRAY, sz);

          /* Put pointer */

          WDB[rea0 + REACTION_PTR_XS] = (double)loc0;

          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Reset add flag */

              add = NO;

              /* Get mt */

              mt = (long)RDB[rea + REACTION_MT];

              /* Check type and mode */

              if (mode == NUCLIDE_PTR_TOTXS)
                {
                  /* Add all partials to total */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                    add = YES;
                }
              else if (mode == NUCLIDE_PTR_SUM_ABSXS)
                {
                  /* Get ty */

                  ty = (long)RDB[rea + REACTION_TY];

                  /* Add partial absorptions to total absorption */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                    if (ty == 0)
                      add = YES;
                }
              else if (mode == NUCLIDE_PTR_PHOTPRODXS)
                {
                  /* Add partial photon production reactions to total */

                  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_SPECIAL)
                    if (mt > 3000)
                      if (rea != rea0)
                        add = YES;
                }
              else if (mode == NUCLIDE_PTR_NFXS)
                {
                  if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
                      && (mt < 18 || (mt > 21 && mt != 38)))
                    add = YES;
                }
              else
                Die(FUNCTION_NAME, "WTF?");

              /* Check if mode was added */

              if (add == NO)
                {
                  /* Pointer to next */

                  rea = NextItem(rea);

                  /* Cycle loop */

                  continue;
                }

              /* Get index to first point and number of points */

              i0 = (long)RDB[rea + REACTION_XS_I0];
              ne = (long)RDB[rea + REACTION_XS_NE];

              /* Check values */

              CheckValue(FUNCTION_NAME, "ne", "", i0 + ne, 2, sz);

              /* Get pointer to data */

              ptr = (long)RDB[rea + REACTION_PTR_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over points and add to sum */

              for (i = 0; i < ne; i++)
                {
                  /* Add value to total */

                  WDB[loc0 + i0 + i] = RDB[loc0 + i0 + i] + RDB[ptr + i];

                  /* Compare to maximum */

                  if (RDB[loc0 + i0 + i] > max)
                    max = RDB[loc0 + i0 + i];
                }

              /* Add reaction to total list */

              ptr = NewItem(rea0 + REACTION_PTR_PARTIAL_LIST,
                            RLS_DATA_BLOCK_SIZE);

              /* Allocate memory for reaction counter */

              ptp = AllocPrivateData(1, PRIVA_ARRAY);
              WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

              WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
              WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
              WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
              WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

              /* Next reaction */

              rea = NextItem(rea);
            }

          /*******************************************************************/
        }

      /* Put maximum cross section */

      if (mode == NUCLIDE_PTR_TOTXS)
        WDB[nuc + NUCLIDE_MAX_TOTXS] = max;

      /* Check if list was created */

      if ((ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST]) > VALID_PTR)
        {
          /* Close list */

          CloseList(ptr);

          /* Sort list by minimum energy */

          SortList(ptr, RLS_DATA_EMIN, SORT_MODE_ASCEND);
        }
    }

  /***************************************************************************/

  /***** Multiply photon heating values with total xs ************************/

  if (type == NUCLIDE_TYPE_PHOTON)
    {
      /* Pointer to total */

      rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Pointer to heating xs */

      rea0 = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_HEATPRODXS];
      CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

      /* Check that energy grids match */

      if ((long)RDB[rea + REACTION_PTR_EGRID] !=
          (long)RDB[rea0 + REACTION_PTR_EGRID])
        Die(FUNCTION_NAME, "Mismatch in energy grid");

      /* Pointer to data */

      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      loc0 = (long)RDB[rea0 + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Check that data starts from zero */

      if ((long)RDB[rea + REACTION_XS_I0] != 0)
        Die(FUNCTION_NAME, "total XS does not start from zero");

      if ((long)RDB[rea0 + REACTION_XS_I0] != 0)
        Die(FUNCTION_NAME, "heat XS does not start from zero");

      /* Number of energy points */

      ne = (long)RDB[rea + REACTION_XS_NE];

      /* Loop over array and multiply */

      for (i = 0; i < ne; i++)
        WDB[loc0 + i] = RDB[loc0 + i]*RDB[ptr + i]*MEV;
    }

  /***************************************************************************/

  /***** Multiply neutron heating values with total xs ***********************/

  /* Check nuclide type (OTF-SAB -moodissa jää nuklideja joiden tyyppi   */
  /* ei ole transport. Ei ihan varmaa että pitäisikö tää koko aliohjelma */
  /* silloin skipata). Ton pointterin nollaaminen on myös vähän viritelmä */
  /* joka pitää tehdä jos vapaalla ytimellä on datassa pelkkää nollaa. */

  if (type != NUCLIDE_TYPE_TRANSPORT)
    return;

  /* Data may have been removed by processnuclides.c, in which case the */
  /* pointer in S(a,b) nuclides is not yet removed. */

  if (sab)
    {
      /* Pointer to free data */

      ptr = (long)RDB[nuc + NUCLIDE_SAB_PTR_FREE];
      CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

      /* Check if pointer is zero and set */

      if (RDB[ptr + NUCLIDE_PTR_HEATPRODXS] < VALID_PTR)
        WDB[nuc + NUCLIDE_PTR_HEATPRODXS] = -1.0;
    }

  if ((long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS] > VALID_PTR)
    {
      /* Pointer to total or non-fission cross section (= total - fission) */

      if (sab)
        {
          /* Get pointer to free nuclide */

          ptr = (long)RDB[nuc + NUCLIDE_SAB_PTR_FREE];
          CheckPointer(FUNCTION_NAME, "(free1)", DATA_ARRAY, ptr);

          /* Check for fission  */

          if (fiss)
            {
              rea = (long)RDB[ptr + NUCLIDE_PTR_NFXS];
              CheckPointer(FUNCTION_NAME, "(sab_fiss_rea)", DATA_ARRAY, rea);
            }
          else
            {
              rea = (long)RDB[ptr + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(sab_rea)", DATA_ARRAY, rea);
            }
        }
      else if (fiss)
        {
          rea = (long)RDB[nuc + NUCLIDE_PTR_NFXS];
          CheckPointer(FUNCTION_NAME, "(fiss_rea)", DATA_ARRAY, rea);
        }
      else
        {
          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
        }

      /* Number of energy points and check first point */

      ne = (long)RDB[rea + REACTION_XS_NE];

      if ((long)RDB[rea + REACTION_XS_I0] != 0)
        Die(FUNCTION_NAME, "Total XS does not start from zero");

      /* Pointer to energy grid */

      ptr = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to energies */

      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      E0 = &RDB[ptr];

      /* Pointer to cross section data */

      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      xs0 = &RDB[ptr];

      /* Pointer to KERMA */

      rea = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS];
      CheckPointer(FUNCTION_NAME, "(kerma)", DATA_ARRAY, rea);

      /* Number of energy points and check first point */

      sz = (long)RDB[rea + REACTION_XS_NE];
      if ((long)RDB[rea + REACTION_XS_I0] != 0)
        Die(FUNCTION_NAME, "heat XS does not start from zero");

      /* Allocate memory for temporary array */

      xs = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

      /* Pointer to energy grid */

      ptr = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to energies */

      pte = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

      /* Reconstruct total */

      n = InterpolateData(&RDB[pte], xs, sz, E0, xs0, ne, 0, NULL, NULL, NO);

      /* Check negative points */

      if (n > 0)
        Warn(FUNCTION_NAME,
             "%ld negative total xs points in data (%s)",
             n, GetText(nuc + NUCLIDE_PTR_NAME));

      /* Pointer to cross section data */

      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over array and multiply */

      for (i = 0; i < sz; i++)
        {
          WDB[ptr + i] = RDB[ptr + i]*xs[i]*MEV;
        }

      /* Free memory */

      Mem(MEM_FREE, xs);
    }

  /***************************************************************************/
}

/*****************************************************************************/
