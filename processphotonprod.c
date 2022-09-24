/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processphotonprod.c                            */
/*                                                                           */
/* Created:       2015/11/07 (JLe)                                           */
/* Last modified: 2010/07/12 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes data for photon production from neutron reactions  */
/*                                                                           */
/* Comments: - Kato noi uresit: fissioiden summa ja mf 13 (mcnp kaatuu jos   */
/*             ures data ja tää)                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPhotonProd:"

/*****************************************************************************/

void ProcessPhotonProd(long nuc)
{
  long ace, ptr, erg, nr, n, i, j, rea, i0, i1, ne, loc0, loc1, loc2, sum;
  long erg0, rea0, NXS[16], JXS[32], NTRP, NES, L0, L1, L2, IE, NE, MFTYPE;
  long MTMULT, LOCA, NR, LAW, NP, pte, ptx, imode;
  const double *XSS;
  double *E0, *xs0, *xs1, *xs, mu, EE[2];

  /* Check mode */

  if ((long)RDB[DATA_PHOTON_PRODUCTION] == NO)
    return;

  /* Check nuclide type */

  if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_TRANSPORT)
    return;

 /* Pointer to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];

  /* Read data to NXS array */

  ptr = (long)ACE[ace + ACE_PTR_NXS];

  for (n = 0; n < 16; n++)
    NXS[n] = (long)ACE[ptr++];

  /* Read data to JXS array */

  ptr = (long)ACE[ace + ACE_PTR_JXS];

  for (n = 0; n < 32; n++)
    JXS[n] = (long)ACE[ptr++];

  /* Get pointer to XSS array */

  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

  /* Check mode : NXS(16) */

  if (NXS[15] == -1)
    {
      /* No photon production */

      return;
    }
  else if (NXS[15] != 0)
    Die(FUNCTION_NAME, "NXS(16) = %ld\n", NXS[15]);

  /* Get number of photon producing reactions : NXS(6) */

  if ((NTRP = NXS[5]) < 1)
    return;

  /* Check strange error in some nuclides (6012.03c in JENDL-3.2) */

  for (nr = 0; nr < NTRP; nr++)
    if ((long)XSS[JXS[17] + nr - 1] == 0)
      {
        /* This fails for some nuclides (6012 in JENDL-3.2) */

        Warn(FUNCTION_NAME, "Error in photon production data");

        /* Exit */

        return;
      }

  /***************************************************************************/

  /***** Make global energy grid *********************************************/

  /* Add nuclide energy grid points */

  E0 = NULL;
  NES = 0;

  E0 = AddPts(E0, &NES, &XSS[JXS[0] - 1], NXS[2]);

  /* Add boundary ponts etc. points in array */

  for (nr = 0; nr < NTRP; nr++)
    {
      /* Pointer to data */

      LOCA = (long)XSS[JXS[13] - 1 + nr];

      /* Check */

      if ((nr == 0) && (LOCA != 1))
        Die(FUNCTION_NAME, "Pointer error");

      /* Get pointer to SIGP-block (Table F-16, page F-31) */

      L0 = LOCA + JXS[14] - 1;

      /* Get MFTYPE */

      MFTYPE = (long)XSS[L0 - 1];

      /* Check */

      if ((MFTYPE == 12) || (MFTYPE == 16))
        {
          /* Get number of interpolation regions */

          NR = (long)XSS[L0 + 1];

          /* Number of energies */

          NE = (long)XSS[L0 + 2*NR + 2];
          CheckValue(FUNCTION_NAME, "NE", "", NE, 2, NXS[2]);

          /* Add point so array */

          E0 = AddPts(E0, &NES, &XSS[L0 + 2*NR + 3], NE);

          /* Add point below minimum and above maximum */

          EE[0] = 0.999999*XSS[L0 + 2*NR + 3];
          EE[1] = 1.000001*XSS[L0 + 2*NR + 3 + NE - 1];

          E0 = AddPts(E0, &NES, EE, 2);

          /* Find reaction */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Compare MT */

              if ((long)RDB[rea + REACTION_MT] == (long)XSS[L0])
                break;

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Check pointer (error tested later) */

          if (rea > VALID_PTR)
            {
              /* Add points below and above minimum */

              EE[0] = 0.999999*RDB[rea + REACTION_EMIN];
              EE[1] = 1.000001*RDB[rea + REACTION_EMIN];

              E0 = AddPts(E0, &NES, EE, 2);

              /* Add points below and above maximum */

              EE[0] = 0.999999*RDB[rea + REACTION_EMAX];
              EE[1] = 1.000001*RDB[rea + REACTION_EMAX];

              E0 = AddPts(E0, &NES, EE, 2);
            }

          /* Check number of interpolation regions */

          if (NR > 0)
            {
              /* Allocate memory for temporary array */

              xs = (double *)Mem(MEM_ALLOC, NR, sizeof(double));

              /* Read points */

              for (n = 0; n < NR; n++)
                {
                  i = (long)XSS[L0 + 2 + n] - 1;
                  xs[n] = XSS[L0 + 2*NR + 3 + i];
                }

              /* Add points below */

              for (n = 0; n < NR; n++)
                xs[n] = 0.999999*xs[n];

              E0 = AddPts(E0, &NES, xs, NR);

              /* Add points above */

              for (n = 0; n < NR; n++)
                xs[n] = 1.000001*xs[n]/0.999999;

              E0 = AddPts(E0, &NES, xs, NR);

              /* Free memory */

              Mem(MEM_FREE, xs);

              /* Check for log interpolation */

              for (n = 0; n < NR; n++)
                {
                  /* Check log-log */

                  if ((long)XSS[L0 + 2 + NR] == 5)
                    {
                      /* Make array of 1000 points */

                      xs = MakeArray(XSS[L0 + 2*NR + 3],
                                     XSS[L0 + 2*NR + 3 + NE - 1], 1000, 2);

                      /* Add to grid */

                      E0 = AddPts(E0, &NES, xs, 1000);

                      /* Free memory */

                      Mem(MEM_FREE, xs);
                    }
                }
            }
        }
    }

  /* Make global energy grid for data */

  erg0 = MakeEnergyGrid(NES, 0, 0, -1, E0, EG_INTERP_MODE_LIN);

  /***************************************************************************/

  /***** Common stuff ********************************************************/

  /* Create structure for total */

  loc0 = NewItem(nuc + NUCLIDE_PTR_PHOTON_PROD, PHOTON_PROD_BLOCK_SIZE);

  /* Put mt */

  WDB[loc0 + PHOTON_PROD_MT] = 1.0;

  /* Put minimum energy to -1 to to ensure that reaction remains on top */
  /* after sort */

  WDB[loc0 + PHOTON_PROD_EMIN] = -1.0;

  /* Create reaction */

  rea0 = NewItem(loc0 + PHOTON_PROD_PTR_PRODXS, REACTION_BLOCK_SIZE);

  /* Put mt */

  WDB[rea0 + REACTION_MT] = 1.0;

  /* Allocate memory for previous value */

  AllocValuePair(rea0 + REACTION_PTR_PREV_XS);

  /* Allocate memory for xs data */

  ptx = ReallocMem(DATA_ARRAY, NES);
  WDB[rea0 + REACTION_PTR_XS] = (double)ptx;

  /* Reset data (just to be sure) */

  for (n = 0; n < NES; n++)
    WDB[ptx + n] = 0.0;

  /* Put size and index */

  WDB[rea0 + REACTION_XS_NE] = (double)NES;
  WDB[rea0 + REACTION_XS_I0] = 0.0;

  /* Put pointer to energy grid */

  WDB[rea0 + REACTION_PTR_EGRID] = (double)erg0;
  WDB[rea0 + REACTION_CACHE_OPTI_IDX] = -1.0;

  /* Loop over reactions */

  for (nr = 0; nr < NTRP; nr++)
    {
      /***********************************************************************/

      /***** Cross sections and general stuff ********************************/

      /* Create structure */

      loc0 = NewItem(nuc + NUCLIDE_PTR_PHOTON_PROD, PHOTON_PROD_BLOCK_SIZE);

      /* Store original mt */

      WDB[loc0 + PHOTON_PROD_MT] = XSS[JXS[12] - 1 + nr];

      /* Create reaction structure for cross section */

      rea0 = NewItem(loc0 + PHOTON_PROD_PTR_PRODXS, REACTION_BLOCK_SIZE);

      /* Put mt */

      WDB[rea0 + REACTION_MT] = XSS[JXS[12] - 1 + nr];

      /* Allocate memory for previous value */

      AllocValuePair(rea0 + REACTION_PTR_PREV_XS);

      /* Pointer to data */

      LOCA = (long)XSS[JXS[13] - 1 + nr];

      /* Check */

      if ((nr == 0) && (LOCA != 1))
        Die(FUNCTION_NAME, "Pointer error");

      /* Get pointer to SIGP-block (Table F-16, page F-31) */

      L0 = LOCA + JXS[14] - 1;

      /* Get MFTYPE */

      MFTYPE = (long)XSS[L0 - 1];

      /* Check */

      if ((MFTYPE == 12) || (MFTYPE == 16))
        {
          /*******************************************************************/

          /***** Multipliers to reaction cross sections **********************/

          /* Allocate memory for temporary data */

          xs0 = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

          /* Get reaction in which multiplier is applied to */

          MTMULT = (long)XSS[L0];

          /* Find reaction */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Compare MT */

              if ((long)RDB[rea + REACTION_MT] == MTMULT)
                break;

              /* Next reaction */

              rea = NextItem(rea);
            }

          /* Check */

          if (rea > VALID_PTR)
            {
              /* Pointer to main cross section data array */

              ptr = (long)RDB[rea + REACTION_PTR_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Pointer to corresponding energy array */

              pte = (long)RDB[rea + REACTION_PTR_EGRID];
              CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

              pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

              /* First point in xs data and number of points */

              i0 = (long)RDB[rea + REACTION_XS_I0];
              ne = (long)RDB[rea + REACTION_XS_NE];

              /* Interpolate data into grid */

              InterpolateData(E0, xs0, NES, &RDB[pte + i0], &RDB[ptr], ne, 0,
                              NULL, NULL, NO);

              /* Link pointer to xs if ures data exists */

              if ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR)
                WDB[loc0 + PHOTON_PROD_PTR_URES_REA] = (double)rea;

              /* Copy minimum energy */

              WDB[loc0 + PHOTON_PROD_EMIN] = RDB[rea + REACTION_EMIN];
            }
          else if (MTMULT == 18)
            {
              /* Calculate total fission from partials, allocate memory */

              xs = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

              /* Reset counter */

              i = 0;

              /* Loop over channels */

              for (n = 0; n < 4; n++)
                {

                  if (n == 0)
                    MTMULT = 19;
                  else if (n == 1)
                    MTMULT = 20;
                  else if (n == 2)
                    MTMULT = 21;
                  else if (n == 3)
                    MTMULT = 38;
                  else
                    Die(FUNCTION_NAME, "Overflow");

                  /* Find reaction */

                  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
                  while (rea > VALID_PTR)
                    {
                      /* Compare MT */

                      if ((long)RDB[rea + REACTION_MT] == MTMULT)
                        break;

                      /* Next reaction */

                      rea = NextItem(rea);
                    }

                  /* Check if found */

                  if (rea < VALID_PTR)
                    continue;

                  /* Update count */

                  i++;

                  /* Reset data array (just to be sure) */

                  for (j = 0; j < NES; j++)
                    xs[j] = 0.0;

                  /* Pointer to main cross section data array */

                  ptr = (long)RDB[rea + REACTION_PTR_XS];
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Pointer to corresponding energy array */

                  pte = (long)RDB[rea + REACTION_PTR_EGRID];
                  CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

                  pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
                  CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);

                  /* First point in xs data and number of points */

                  i0 = (long)RDB[rea + REACTION_XS_I0];
                  ne = (long)RDB[rea + REACTION_XS_NE];

                  /* Interpolate data into grid */

                  InterpolateData(E0, xs, NES, &RDB[pte + i0], &RDB[ptr], ne,
                                  0, NULL, NULL, NO);

                  /* Add to total */

                  for (j = 0; j < NES; j++)
                    xs0[j] = xs0[j] + xs[j];
                }

              /* Check count */

              if (i == 0)
                Die(FUNCTION_NAME, "No fission channels found");

              /* Free temprorary array */

              Mem(MEM_FREE, xs);
            }
          else
            {
              /* Free temporary array */

              Mem(MEM_FREE, xs0);

              /* Remove item from list */

              RemoveItem(loc0);

              /* Print warning */

#ifdef DEBUG

              Warn(FUNCTION_NAME, "Reaction mt %ld has no xs points (removed)",
                   (long)XSS[L0]);

#endif

              /* Cycle loop */

              continue;
            }

          /* Avoid compiler warning */

          imode = -1;

          /* Allocate memory for data */

          xs1 = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

          /* Get number of interpolation regions */

          NR = (long)XSS[L0 + 1];

          /* Number of energies */

          NE = (long)XSS[L0 + 2*NR + 2];
          CheckValue(FUNCTION_NAME, "NE", "", NE, 2, NES);

          /* Check number of interpolation regions */

          if (NR < 2)
            {
              /* Set interpolation for InterpolateData() */

              if (NR == 0)
                imode = 0;
              else if ((long)XSS[L0 + 2 + NR] == 2)
                imode = 0;
              else if ((long)XSS[L0 + 2 + NR] == 1)
                imode = 3;
              else if ((long)XSS[L0 + 2 + NR] == 5)
                imode = 5;
              else
                Die(FUNCTION_NAME, "Invalid interpolation mode (%ld %ld)",
                    NR, (long)XSS[L0 + 2 + NR]);

              /* Check that region index matche number of points */

              if (NR == 1)
                if ((long)XSS[L0 + 2] != NE)
                  Die(FUNCTION_NAME, "Mismatch in indexes");

              /* Check first and last points */

              if (XSS[L0 + 2*NR + 3] < XSS[JXS[0] - 1])
                Die(FUNCTION_NAME, "Lower limit below nuclide minimum");
              else if (XSS[L0 + 2*NR + 3 + NE - 1] >
                       XSS[JXS[0] - 1 + NXS[2] - 1])
                Die(FUNCTION_NAME, "Upper limit above nuclide minimum");

#ifdef DEBUG
              if (imode == 5)
                Warn(FUNCTION_NAME, "Log-log interpolation of xs");
#endif

              /* Interpolate multipliers */

              InterpolateData(E0, xs1, NES, &XSS[L0 + 2*NR + 3],
                              &XSS[L0 + 2*NR + 3 + NE], NE, imode, NULL, NULL, NO);
            }
          else
            {

#ifdef DEBUG

              printf("\n-------------------------\n");
              printf("mt = %ld\n", (long)XSS[JXS[12] - 1 + nr]);
              printf("XS interpolation regions (%ld):\n",
                     (long)XSS[L0 + 2*NR + 2]);
              for (n = 0; n < NR; n++)
                printf("NBT = %E INT = %E\n",
                       XSS[L0 + 2 + n], XSS[L0 + 2 + NR + n]);
              printf("-------------------------\n");

#endif

              /* Allocate memory for data */

              xs = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

              /* Reset index */

              i0 = 0;

              /* Loop over regions */

              for (i = 0; i < NR; i++)
                {
                  /* Get interpolation region index */

                  i1 = (long)XSS[L0 + 2 + i] - 1;

                  /* Set interpolation for InterpolateData() */

                  if ((long)XSS[L0 + 2 + NR + i] == 2)
                    imode = 0;
                  else if ((long)XSS[L0 + 2 + NR + i] == 1)
                    imode = 3;
                  else if ((long)XSS[L0 + 2 + NR] == 5)
                    imode = 5;
                  else
                    Die(FUNCTION_NAME, "Invalid interpolation mode (%ld %ld)",
                        NR,(long)XSS[L0 + 2 + NR]);

                  /* Reset temporary array (just to be sure) */

                  for (n = 0; n < NES; n++)
                    xs[n] = 0.0;

                  /* Array size */

                  ne = i1 - i0 + 1;

#ifdef DEBUG
                  if (imode == 5)
                    Warn(FUNCTION_NAME, "Log-log interpolation of xs");
#endif

                  /* Interpolate multipliers */

                  InterpolateData(E0, xs, NES, &XSS[L0 + 2*NR + 3 + i0],
                          &XSS[L0 + 2*NR + 3 + NE + i0], ne, imode, NULL, NULL, NO);


                  /* Add to final array */

                  for (n = 0; n < NES; n++)
                    if (E0[n] <= XSS[L0 + 2*NR + 3 + i1])
                      xs1[n] = xs1[n] + xs[n];

                  /* Update index */

                  i0 = i1;
                }

              /* Free temporary array */

              Mem(MEM_FREE, xs);
            }

          /* Find first non-zero point */

          for (i0 = 0; i0 < NES; i0++)
            if (xs0[i0]*xs1[i0] > 0.0)
              break;

          /* Check if no non-zero points */

          if (i0 == NES)
            {
              /* Free temporary arrays */

              Mem(MEM_FREE, xs0);
              Mem(MEM_FREE, xs1);

              /* Remove item from list */

              RemoveItem(loc0);

              /* Print warning */

#ifdef DEBUG

              Warn(FUNCTION_NAME, "Reaction %ld has zero xs (removed)",
                   (long)XSS[L0]);

#endif

              /* Cycle loop */

              continue;
            }

          /* Add one point for leading zero */

          if (i0 > 0)
            i0--;

          /* Calculate array size */

          ne = NES - i0;

          /* Allocate memory for data */

          ptr = ReallocMem(DATA_ARRAY, ne);
          WDB[rea0 + REACTION_PTR_XS] = (double)ptr;

          /* Put data and add to total */

          for (n = 0; n < ne; n++)
            {
              WDB[ptr++] = xs1[i0 + n]*xs0[i0 + n];
              WDB[ptx + i0 + n] = RDB[ptx + i0 + n] + xs1[i0 + n]*xs0[i0 + n];
            }

          /* Put size and index */

          WDB[rea0 + REACTION_XS_NE] = (double)ne;
          WDB[rea0 + REACTION_XS_I0] = (double)i0;

          /* Put pointer to energy grid */

          WDB[rea0 + REACTION_PTR_EGRID] = (double)erg0;
          WDB[rea0 + REACTION_CACHE_OPTI_IDX] = -1.0;

          /* Free temporary arrays */

          Mem(MEM_FREE, xs0);
          Mem(MEM_FREE, xs1);

          /*******************************************************************/
        }
      else if (MFTYPE == 13)
        {
          /*******************************************************************/

          /***** Direct cross sections ***************************************/

          /* Index to first point */

          IE = (long)XSS[L0] - 1;
          CheckValue(FUNCTION_NAME, "NE", "", IE, 0, NES);

          /* Number of energies */

          NE = (long)XSS[L0 + 1];
          CheckValue(FUNCTION_NAME, "NE", "", NE, 2, NES);

          /* Allocate memory for data */

          xs1 = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

          /* Interpolate data into nuclide energy grid */

          InterpolateData(E0, xs1, NES, &XSS[JXS[0] - 1 + IE], &XSS[L0 + 2],
                          NE, 0, NULL, NULL, NO);

          /* Find first non-zero point */

          for (i0 = 0; i0 < NES; i0++)
            if (xs1[i0] > 0.0)
              break;

          /* Add one point for leading zero */

          if (i0 > 0)
            i0--;

          /* Calculate array size */

          ne = NES - i0;

          /* Allocate memory for data */

          ptr = ReallocMem(DATA_ARRAY, ne);
          WDB[rea0 + REACTION_PTR_XS] = (double)ptr;

          /* Put data and add to total */

          for (n = 0; n < ne; n++)
            {
              WDB[ptr++] = xs1[i0 + n];
              WDB[ptx + i0 + n] = RDB[ptx + i0 + n] + xs1[i0 + n];
            }

          /* Put size and index */

          WDB[rea0 + REACTION_XS_NE] = (double)ne;
          WDB[rea0 + REACTION_XS_I0] = (double)i0;

          /* Put pointer to energy grid */

          WDB[rea0 + REACTION_PTR_EGRID] = (double)erg0;
          WDB[rea0 + REACTION_CACHE_OPTI_IDX] = -1.0;

          /* Free temporary array */

          Mem(MEM_FREE, xs1);

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "MFTYPE = %ld\n", MFTYPE);

      /***********************************************************************/

      /***** Angular distributions *******************************************/

      /* Pointer */

      LOCA = (long)XSS[JXS[15] + nr - 1];

      /* Check isotropic */

      if (LOCA != 0)
        {
          /* Get pointer to LANDP-block (Table F-17, page F-32) */

          L0 = LOCA + JXS[16] - 1;

          /* Energies */

          if ((NE = (long)XSS[L0 - 1]) < 2)
            Warn(FUNCTION_NAME, "No energy points in angular distribution");
          else
            {
              /* Create structure */

              loc1 = NewItem(loc0 + PHOTON_PROD_PTR_ANG, ANG_BLOCK_SIZE);

              /* Make energy grid */

              erg = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L0], EG_INTERP_MODE_LIN);
              WDB[loc1 + ANG_PTR_EGRID] = (double)erg;

              /* Set distribution type and number of bins */

              WDB[loc1 + ANG_TYPE] = (double)ANG_TYPE_EQUIBIN;
              WDB[loc1 + ANG_BINS] = 32.0;

              /* Allocate memory for data */

              ptr = ReallocMem(DATA_ARRAY, 33*NE);

              /* Set pointer */

              WDB[loc1 + ANG_PTR_D0] = (double)ptr;

              /* Loop over bins */

              for (n = 0; n < NE; n++)
                {
                  /* Pointer */

                  L1 = (long)XSS[L0 + NE + n];

                  /* Check isotropic */

                  if (L1 == 0)
                    {
                      /* Isotropic distribution, calculate values */

                      for (i = 0; i < 33; i++)
                        WDB[ptr++] = (double)i/16.0 - 1.0;
                    }
                  else if (L1 > 0)
                    {
                      /* Anisotropic distribution, get pointer to table */

                      L1 = L1 + JXS[16] - 2;

                      /* Check first and last value */

                      if (XSS[L1] != -1.0)
                        Die(FUNCTION_NAME, "Error in min cosine");
                      else if (XSS[L1 + 32] != 1.0)
                        Die(FUNCTION_NAME, "Error in max cosine");

                      /* Loop over data */

                      for (i = 0; i < 33; i++)
                        {
                          /* Get value */

                          mu = XSS[L1 + i];
                          CheckValue(FUNCTION_NAME, "mu", " (mu)", mu,
                                     -1.0, 1.0);

                          /* Truncate */

                          if (mu > 1.0)
                            mu = 1.0;
                          else if (mu < -1.0)
                            mu = -1.0;

                          /* Set value */

                          WDB[ptr++] = mu;
                        }
                    }
                  else
                    Die(FUNCTION_NAME, "Pointer error");
                }
            }
        }

      /***********************************************************************/

      /***** Energy distributions ********************************************/

      /* Create structure */

      loc1 = NewItem(loc0 + PHOTON_PROD_PTR_ERG, ERG_BLOCK_SIZE);

      /* Pointer to nuclide data */

      WDB[loc1 + ERG_PTR_NUCLIDE] = (double)nuc;

      /* Pointer */

      LOCA = (long)XSS[JXS[17] + nr - 1];

      /* Get pointer to LDLWP-block (Table F-14, page F-20) */

      L0 = LOCA + JXS[18] - 1;

      /* Check LNW */

      if ((long)XSS[L0 - 1] != 0)
        {
          /* Print warning */

          Warn(FUNCTION_NAME,
               "LNW = %ld, skipping photon production data for %s mt %ld",
               (long)XSS[L0 - 1], GetText(nuc + NUCLIDE_PTR_NAME),
               (long)RDB[loc0 + PHOTON_PROD_MT]);

          /* Remove mt */

          RemoveItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Get law */

      LAW = (long)XSS[L0];
      WDB[loc1 + ERG_LAW] = (double)LAW;

      /* Get interpolation regions */

      NR = (long)XSS[L0 + 2];

#ifdef DEBUG

      if (NR > 0)
        {
          printf("--------------------------\n");
          printf("mt = %ld\n", (long)XSS[JXS[12] - 1 + nr]);
          printf("ERG interpolation regions (%ld):\n",
                 (long)XSS[L0 + 2*NR + 3]);

          for (n = 0; n < NR; n++)
            printf("NBT = %E INT = %E\n",
                   XSS[L0 + 3 + n], XSS[L0 + NR + 3 + n]);
          printf("--------------------------\n");
        }

#endif

      /* Check */
      /*
      if (NR > 1)
        Warn(FUNCTION_NAME, "Interpolation data (NR = %ld)", NR);
      */
      /* Get number of energies */

      NE = (long)XSS[L0 + 2*NR + 3];

      /* Check probabilities */

      for (n = 0; n < NE; n++)
        if (XSS[L0 + 2*NR + 4 + NE + n] != 1.0)
          Warn(FUNCTION_NAME, "Law probability %E %ld\n",
              XSS[L0 + 2*NR + 4 + NE + n], NE);

      /* Pointer to data */

      L1 = JXS[18] - 1 + (long)XSS[L0 + 1] - 1;

      /* Check */

      if (LAW == 2)
        {
          /*******************************************************************/

          /***** 2 = Discrete photon energy **********************************/

          /* Check indicator */

          if ((XSS[L1] != 0.0) && (XSS[L1] != 1.0) && (XSS[L1] != 2.0))
            Die(FUNCTION_NAME, "Invalid primary indicator");

          /* Check energy */

          CheckValue(FUNCTION_NAME, "EG", "", XSS[L1 + 1], 1E-6, 100.0);

          /* Allocate memory for data */

          ptr = ReallocMem(DATA_ARRAY, 2);
          WDB[loc1 + ERG_PTR_DATA] = (double)ptr;

          /* Put data */

          WDB[ptr] = XSS[L1];
          WDB[ptr + 1] = XSS[L1 + 1];

          /*******************************************************************/
        }
      else if (LAW == 4)
        {
          /*******************************************************************/

          /*****  4 = Continuous tabular distribution ************************/

          /* Number of interpolation regions */

          NR = (long)XSS[L1];

          /* Check */

          if (NR == 0)
            WDB[loc1 + ERG_INTERP] = 0.0;
          else if (NR == 1)
            WDB[loc1 + ERG_INTERP] = XSS[L1 + 1 + NR];
          else if (NR > 1)
            {
              /* Set type according to first region */

              WDB[loc1 + ERG_INTERP] = XSS[L1 + 1 + NR];

              /* Check if remaining regions differ */

              for (n = 0; n < NR; n++)
                if (RDB[loc1 + ERG_INTERP] != XSS[L1 + 1 + NR + n])
                  {
                    /* Set interpolation to -1 to indicate multiple regions */

                    WDB[loc1 + ERG_INTERP] = -1.0;

                    /* Allocate memory for data */

                    loc2 = ReallocMem(DATA_ARRAY, 2*NR + 2);

                    /* Put number of regions */

                    WDB[loc1 + ERG_NR] = (double)NR + 1;

                    /* Put pointer */

                    WDB[loc1 + ERG_PTR_INTERP] = (double)loc2;

                    /* Put first index */

                    WDB[loc2++] = -1.0;

                    /* Read indexes */

                    for (i = 0; i < NR; i++)
                      WDB[loc2++] = XSS[L1 + 1 + i] - 1.0;

                    /* Read interpolations */

                    for (i = 0; i < NR; i++)
                      WDB[loc2++] = XSS[L1 + 1 + NR + i];

                    /* Put last interpolation */

                    WDB[loc2++] = -1.0;

                    /* Break loop */

                    break;
                  }
            }

          /* Number of energies */

          NE = (long)XSS[L1 + 2*NR + 1];

          /* Create energy grid */

          erg = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L1 + 2*NR + 2],
                               EG_INTERP_MODE_LIN);

          /* Allocate memory for types */

          loc2 = ReallocMem(DATA_ARRAY, NE);

          /* Calculate total data size and put types  */

          sum = 0;
          for(i = 0; i < NE; i++)
            {
              /* Get pointer */

              L2 = JXS[18] + (long)XSS[L1 + 2*NR + 2 + NE + i] - 2;

              /* Muutettu 31.3.2016 (2.1.26) */

              /*
              if ((long)XSS[L2] > 9998)
                Die(FUNCTION_NAME, "XSS[L2] = %E\n", XSS[L2]);
              else
                WDB[loc2 + i] = (double)XSS[L2];
              */

              /* Check (ei muista enää mistä toi raja tuli, mutta siinä */
              /* ei oo mitään järkeä) 2.12.2017 / 2.1.30 / JLe. */
              /*
              if ((long)XSS[L2] > 9998)
                Warn(FUNCTION_NAME, "INTT' = %ld", (long)XSS[L2]);
              */
              /* Put type */

              WDB[loc2 + i] = (double)XSS[L2];

              /* Add to sum */

              sum = sum + 3*(long)XSS[L2 + 1] + 1;
            }

          /* Allocate memory for data */

          ptr = ReallocMem(DATA_ARRAY, NE + sum + 2);

          /* Set pointer */

          WDB[loc1 + ERG_PTR_DATA] = (double)ptr;

          /* Set pointer to incident energy grid */

          WDB[ptr++] = (double)erg;

          /* Set pointer to types */

          WDB[ptr++] = (double)loc2;

          /* Pointer to pointers */

          loc2 = ptr;
          ptr = ptr + NE;

          /* Loop over initial energies */

          for(i = 0; i < NE; i++)
            {
              /* Set pointer to distribution */

              WDB[loc2++] = (double)ptr;

              /* Get pointer to data */

              L2 = JXS[18] + (long)XSS[L1 + 2*NR + 2 + NE + i] - 2;

              /* Get number of energy points */

              NP = (long)XSS[L2 + 1];
              WDB[ptr++] = (double)NP;

              /*
              printf("----------------- %E %E\n", XSS[L2], XSS[L2 + 1]);

              for(j = 0; j < NP; j++)
                printf("%E %E %E\n",
                       XSS[L2 + 2 + j],
                       XSS[L2 + 2 + NP + j],
                       XSS[L2 + 2 + 2*NP + j]);
              */
              /* Read energy values. */

              for(j = 0; j < NP; j++)
                WDB[ptr++] = XSS[L2 + 2 + j];

              /* Read pdf values. */

              for(j = 0; j < NP; j++)
                WDB[ptr++] = XSS[L2 + NP + 2 + j];

              /* Check cdf values */

              for(j = 1; j < NP; j++)
                if (XSS[L2 + 2 + 2*NP + j - 1] > XSS[L2 + 2 + 2*NP + j])
                  Die(FUNCTION_NAME, "CDF values not in ascending order");

              /* Read cdf values. */

              for(j = 0; j < NP; j++)
                WDB[ptr++] = XSS[L2 + 2*NP + 2 + j];

              /* Check last value (first may be > 0 if discrete lines) */

              if (NP > 1)
                if (XSS[L2 + 2*NP + 2 + NP - 1] != 1.0)
                  Warn(FUNCTION_NAME, "CDF(NP - 1) = %E",
                       XSS[L2 + 2*NP + 2 + NP - 1]);
            }

          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr - 1);

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "Unsupported law %ld", LAW);

      /***********************************************************************/
    }

  /* Sort list by ascending energy */

  loc0 = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_PROD];
  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

  /* Check list size */

  if (ListSize(loc0) == 1)
    {
      /* No reactions, remove */

      WDB[nuc + NUCLIDE_PTR_PHOTON_PROD] = -1.0;

      /* Print warning */

#ifdef DEBUG

      Warn(FUNCTION_NAME, "All reactions of %s were removed",
           GetText(nuc + NUCLIDE_PTR_NAME));

#endif

    }
  else
    {
      /* Sort list */

      SortList(loc0, PHOTON_PROD_EMIN, SORT_MODE_ASCEND);

      /* Add counter */

      WDB[DATA_N_PHOTON_PROD_DATA] = RDB[DATA_N_PHOTON_PROD_DATA] + 1.0;
    }

  /***************************************************************************/
}

/*****************************************************************************/
