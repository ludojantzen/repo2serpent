/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processedistributions.c                        */
/*                                                                           */
/* Created:       2010/02/05 (JLe)                                           */
/* Last modified: 2019/11/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes ENDF energy distributions in ACE data              */
/*                                                                           */
/* Comments: - Adopted from energydistribution.c in Serpent 1.1.12           */
/*                                                                           */
/*           - 2.1.22 (15.10.2014) LAW 4/44/61 lisättiin toi energiariippuva */
/*             tyyppivektori                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessEDistributions:"

/*****************************************************************************/

void ProcessEDistributions(long rea, long prec)
{
  long mt, nuc, ace, ptr, n, i, j, k, nr, loc0, l0, l1, l2, erg, sum;
  long NXS[16], JXS[32], LED, LDIS, L0, L1, L2, L3, NR, NE, NP, NP2, LAW;
  long NMU, NPEP, NES, L, TY, JED, KY;
  const double *XSS;
  double E, mu, xs;

  /***************************************************************************/

  /***** Get pointers etc. ***************************************************/

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  /* Check mt */

  if ((mt != 5) && (mt != 11) && (mt < 15))
    return;
  else if ((mt > 91) && (mt < 875))
    return;
  else if ((mt > 891) && (mt != 1002) && (mt != 1004))
    return;

  /* Pointer to nuclide data */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

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

  /* Reset pointers to avoid compiler warning */

  LDIS = -1;
  LED = -1;
  nr = -1;
  loc0 = -1;

  /* Check data type */

  if (prec > VALID_PTR)
    {
      /* Delayed neutron energy distribution data */

      LED = JXS[25];
      LDIS = JXS[26];

      /* Get precursor group index */

      nr = (long)RDB[prec + PREC_IDX] - 1;

      /* Check precursor group */

      if (nr > NXS[7] - 1)
        Die(FUNCTION_NAME, "nr > NXS[7] - 1");

      /* Check pointers */

      if ((LED <= 0) || (LDIS <= 0))
        Die(FUNCTION_NAME,
            "Delayed neutron energy distribution not given");

      /* Set distrubution pointer */

      loc0 = prec + PREC_PTR_ERG;
    }
  else if ((mt != 1002) && (mt != 1004))
    {
      /* Energy distribution data */

      LED = JXS[9];
      LDIS = JXS[10];

      /* Get reaction number */

      nr = (long)RDB[rea + REACTION_NR];

      /* Check reaction number (NXS[4] = number of reactions having */
      /* secondary neutrons excluding elastic. */

      if (nr > NXS[4] - 1)
        Die(FUNCTION_NAME, "nr > NXS[4] - 1");

      /* Check pointers */

      if ((LED <= 0) || (LDIS <= 0))
        Die(FUNCTION_NAME, "Energy distribution not given (mt %ld)", mt);

      /* Set distrubution pointer */

      loc0 = rea + REACTION_PTR_ERG;
    }

  /***************************************************************************/

  /***** S(a,b) data *********************************************************/

  if ((mt == 1002) || (mt == 1004))
    {
      /* Allocate memory for block */

      erg = NewItem(rea + REACTION_PTR_ERG, ERG_BLOCK_SIZE);

      /* Pointer to nuclide data */

      WDB[erg + ERG_PTR_NUCLIDE] = (double)nuc;

      /* Check mt */

      if (mt == 1004)
        {
          /*******************************************************************/

          /***** Inelastic mode **********************************************/

          /* Set law */

          WDB[erg + ERG_LAW] = 1004.0;

          /* Pointer data (Table F-23, page F-36) */

          L0 = JXS[0] - 1;

          /* Get number of original energy points */

          NES = (long)XSS[L0];

          /* Pointer to energy array */

          L = JXS[0];

          /* Create energy grid */

          ptr = MakeEnergyGrid(NES, 0, 0, -1, &XSS[L], EG_INTERP_MODE_LIN);

          /* Check type */

          if (NXS[6] < 2)
            {
              /***************************************************************/

              /***** Discrete distribution ***********************************/

              /* Allocate memory for data */

              l0 = ReallocMem(DATA_ARRAY, NES*NXS[3]*(NXS[2] + 2) + 5 + NES);

              /* Set pointer */

              WDB[erg + ERG_PTR_DATA] = (double)l0;

              /* Write mode */

              WDB[l0++] = (double)NXS[1];

              /* Set pointer to incident energy grid */

              WDB[l0++] = (double)ptr;

              /* Write sizes */

              WDB[l0++] = (double)NXS[3];
              WDB[l0++] = (double)(NXS[2] + 1);

              /* Write type */

              WDB[l0++] = (double)NXS[6];

              /* Pointer to data */

              L = JXS[2] - 1;

              /* Pointer to distributions */

              l1 = l0 + NES;

              /* Loop over incoming energies */

              for (n = 0; n < NES; n++)
                {
                  /* Set pointer */

                  WDB[l0++] = (double)l1;

                  /* Loop over energy-cosine bins */

                  for (i = 0; i < NXS[3]; i++)
                    {
                      /* Get energy value */

                      E = XSS[L];

                      /* Check value */

                      CheckValue(FUNCTION_NAME, "E1", " (sab)", E, 1E-13, 1E-3);

                      /* Set value */

                      WDB[l1++] = E;

                      /* Compare to minimum and maximum emission energies */

                      if (E < RDB[rea + REACTION_SAB_MIN_EM_E])
                        WDB[rea + REACTION_SAB_MIN_EM_E] = E;
                      if (E > RDB[rea + REACTION_SAB_MAX_EM_E])
                        WDB[rea + REACTION_SAB_MAX_EM_E] = E;

                      /* Loop over cosines */

                      for (j = 0; j < NXS[2] + 1; j++)
                        {
                          /* Get value */

                          mu = XSS[L + 1 + j];

                          /* Check limits */

                          if (mu > 1.0)
                            {
                              Warn(FUNCTION_NAME, "mu = %E (%s %ld)", mu,
                                   GetText(nuc + NUCLIDE_PTR_NAME), mt);
                              mu = 1.0;
                            }
                          else if (mu < -1.0)
                            {
                              Warn(FUNCTION_NAME, "mu = %E (%s %ld)", mu,
                                   GetText(nuc + NUCLIDE_PTR_NAME), mt);
                              mu = -1.0;
                            }

                          /* Set value */

                          WDB[l1++] = mu;
                        }

                      /* Update pointer */

                      L = L + NXS[2] + 2;
                    }
                }

              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "S(a,b) inl", DATA_ARRAY,  l1 - 1);

              /***************************************************************/
            }
          else
            {
              /***************************************************************/

              /***** Continuous distribution *********************************/

              /* Get pointer to data */

              L = JXS[2] - 1;

              /* Allocate memory for data */

              l0 = ReallocMem(DATA_ARRAY, NES + 5);

              /* Set pointer */

              WDB[erg + ERG_PTR_DATA] = (double)l0;

              /* Write mode */

              WDB[l0++] = (double)NXS[1];

              /* Set pointer to incident energy grid */

              WDB[l0++] = (double)ptr;

              /* Write number of sizes */

              WDB[l0++] = (double)NXS[3];
              WDB[l0++] = (double)(NXS[2] - 1);

              /* Write type */

              WDB[l0++] = (double)NXS[6];

              /* Avoid compiler warning */

              l2 = -1;

              /* Loop over incoming energies */

              for (n = 0; n < NES; n++)
                {
                  /* Get pointer to distribution */

                  L1 = (long)XSS[L + n] - 1;

                  /* Number of outgoing energies */

                  NE = (long)XSS[L + NES + n];

                  /* Allocate memory for data */

                  l1 = ReallocMem(DATA_ARRAY, 4*NE + 1);

                  /* Set pointer */

                  WDB[l0++] = (double)l1;

                  /* Put number of energy points */

                  WDB[l1++] = (double)NE;

                  /* Loop over incident energies */

                  for (i = 0; i < NE; i++)
                    {
                      /* Allocate memory for cosine distributions */

                      l2 = ReallocMem(DATA_ARRAY, NXS[2] - 1);

                      /* Write emission energy */

                      WDB[l1] = XSS[L1 + 1];

                      /* Write pdf */

                      WDB[l1 + NE] = XSS[L1 + 2];

                      /* Write cdf */

                      WDB[l1 + 2*NE] = XSS[L1 + 3];

                      /* Write pointer to cosine distributions */

                      WDB[l1 + 3*NE] = (double)l2;

                      /* Update pointer */

                      l1++;

                      /* Loop over cosines */

                      for (j = 0; j < NXS[2] - 1; j++)
                        {
                          /* Get value */

                          mu = XSS[L1 + 4 + j];

                          /* Write */

                          WDB[l2++] = mu;
                        }

                      /* Update pointer */

                      L1 = L1 + 4 + NXS[2] - 2;
                    }
                }

              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "S(a,b) cont", DATA_ARRAY,  l2 - 1);

              /* Pointer to data */

              l0 = (long)RDB[erg + ERG_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

              /* Reset sum */

              sum = 0;

              /* Loop over incident energies */

              for (n = 0; n < NES; n++)
                {
                  l1 = (long)RDB[l0 + 5 + n];
                  CheckPointer(FUNCTION_NAME, "(l1)", DATA_ARRAY, l1);

                  /* Get number of emission energies */

                  NE = (long)RDB[l1++];
                  CheckValue(FUNCTION_NAME, "NE", "", NE, 2, 5000);

                  /* Loop over values and check */

                  for (i = 0; i < NE; i++)
                    {
                      /* Check emission energy */

                      CheckValue(FUNCTION_NAME, "E", "", RDB[l1 + i],
                                 1E-15, 1E-5);

                      /* Check pdf */

                      CheckValue(FUNCTION_NAME, "pdf", "", RDB[l1 + NE + i],
                                 -1E-12, INFTY);

                      /* Check cdf */

                      CheckValue(FUNCTION_NAME, "cdf", "", RDB[l1 + 2*NE + i],
                                 0.0, 1.0);

                      /* Pointer to cosine distribution */

                      l2 = (long)RDB[l1 + 3*NE + i];
                      CheckPointer(FUNCTION_NAME, "(l2)", DATA_ARRAY, l2);

                      /* Loop over distribution */

                      for (j = 0; j < (long)RDB[l0 + 3]; j++)
                        {
                          /* Check cosine */

                         CheckValue(FUNCTION_NAME, "mu", "", RDB[l2 + j],
                                    -12.0, 12.0);

                         /* Add to count */

                         if ((RDB[l2 + j] < -1.0) || (RDB[l2 + j] > 1.0))
                           sum++;
                        }
                    }
#ifdef DEBUG
                  /* Check first and last value of CDF */

                  if (RDB[l1 + 2*NE] > 1E-2)
                    Warn(FUNCTION_NAME, "First value of CDF = %E\n",
                         RDB[l1 + 2*NE]);

                  if (RDB[l1 + 2*NE + NE - 1] < 1.0)
                    Warn(FUNCTION_NAME, "Last value of CDF = %E\n",
                         RDB[l1 + 2*NE + NE - 1]);
#endif
                }

#ifdef DEBUG
                  /* Check number of invalid cosine values */

                  if (sum > 0)
                    Warn(FUNCTION_NAME, "%ld invalid cosine values", sum);
#endif

              /***************************************************************/
            }

          /*******************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Elastic mode ***********************************************/

          /* Set law */

          WDB[erg + ERG_LAW] = 1002.0;

          /* Pointer data (Table F-23, page F-36) */

          L0 = JXS[3] - 1;

          /* Get number of original energy points */

          NES = (long)XSS[L0];

          /* Pointer to energy array */

          L = JXS[3];

          /* Create energy grid */

          ptr = MakeEnergyGrid(NES, 0, 0, -1, &XSS[L], EG_INTERP_MODE_LIN);

          /* Allocate memory for data */

          if (NXS[4] != 4)
            l0 = ReallocMem(DATA_ARRAY, NES*(NXS[5] + 1) + 3);
          else
            l0 = ReallocMem(DATA_ARRAY, NES + 3);

          /* Set pointer */

          WDB[erg + ERG_PTR_DATA] = (double)l0;

          /* Write mode (use minus sign to indicate continuous */
          /* cosine distribution) */

          if ((NXS[6] == 2) && (NXS[4] != 4))
            WDB[l0++] = (double)(-NXS[4]);
          else
            WDB[l0++] = (double)NXS[4];

          /* Set pointer to incident energy grid */

          WDB[l0++] = (double)ptr;

          /* Write size */

          WDB[l0++] = (double)labs(NXS[5]) + 1.0;

          /* Pointer to data */

          L = JXS[5] - 1;

          /* Check type */

          if (NXS[4] != 4)
            {
              /* Loop over incoming energies and cosine bins */

              for (n = 0; n < NES; n++)
                for(j = 0; j < NXS[5] + 1; j++)
                  {
                    /* Get value */

                    mu = XSS[L++];

                      /* Check limits */

                      if (mu > 1.0)
                        {
                          Warn(FUNCTION_NAME, "mu = %E (%s %ld)", mu,
                               GetText(nuc + NUCLIDE_PTR_NAME), mt);
                          mu = 1.0;
                        }
                      else if (mu < -1.0)
                        {
                          Warn(FUNCTION_NAME, "mu = %E (%s %ld)", mu,
                               GetText(nuc + NUCLIDE_PTR_NAME), mt);
                          mu = -1.0;
                        }

                    /* Set value */

                    WDB[l0++] = mu;
                  }
            }
          else
            {
              /* Exact treatment. Check number of energies */

              if (NXS[5] > -1)
                Die(FUNCTION_NAME, "ERROR");

              /* Write cross section values */

              L = JXS[3] + NES;

              for (n = 0; n < NES; n++)
                {
                  /* Get value */

                  xs = XSS[L++];

                  /* Check value */

                  CheckValue(FUNCTION_NAME, "D/E", " (sab)", xs, 1E-10, 10.0);

                  /* Set value */

                  WDB[l0++] = xs;
                }

              /* Check pointer */

              CheckPointer(FUNCTION_NAME, "S(a,b) ela", DATA_ARRAY, l0 - 1);
            }

          /*******************************************************************/
        }

      /* Exit subroutine */

      return;
    }

  /***************************************************************************/

  /***** ENDF reaction laws **************************************************/

  /* Get location of energy distribution data (Table F-13, Page F-19). */

  L0 = (long)XSS[LED + nr - 1];

  /* Loop over laws */

  while (L0 > 0)
    {
      /* Allocate memory for block */

      erg = NewItem(loc0, ERG_BLOCK_SIZE);

      /* Pointer to nuclide data */

      WDB[erg + ERG_PTR_NUCLIDE] = (double)nuc;

      /* Adjust pointer */

      L0 = L0 - 1 + LDIS;

      /* Get number of interpolation regions. */

      NR = (long)XSS[L0 + 2];

      /* Check interpolation data */

      if ((NR > 1) ||
          ((NR == 1) && (XSS[L0 + 3 + 2*NR + (long)XSS[L0 + 3]] < 10.0)))
        Die(FUNCTION_NAME, "Unable to process E interpolation (%s mt %ld)",
            GetText(nuc + NUCLIDE_PTR_NAME), mt);

      /* Get number of energy-probability pairs */

      NE = (long)XSS[L0 + 3 + 2*NR];

      /* Read energy-probability pairs */

      if (NE > 0)
        {
          /* Create energy grid */

          ptr = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L0 + 4 + 2*NR],
                               EG_INTERP_MODE_LIN);

          WDB[erg + ERG_PTR_EGRID] = (double)ptr;

          /* Allocate memory for probabilities */

          ptr = ReallocMem(DATA_ARRAY, NE);

          /* Set pointer */

          WDB[erg + ERG_PTR_PROB] = (double)ptr;

          /* Read P */

          for (n = 0; n < NE; n++)
            WDB[ptr++] = XSS[L0 + 4 + 2*NR + NE + n];
        }

      /* Get law type */

      LAW = (long)XSS[L0];
      WDB[erg + ERG_LAW] = (double)LAW;

      /* Pointer to data in XSS */

      L1 = LDIS + (long)XSS[L0 + 1] - 2;

      /* Reset interpolation type */

      WDB[erg + ERG_INTERP] = 0.0;

      /* Read number of interpolation schemes (näitä ei käytetä) */

      if ((LAW == 3) || (LAW == 66))
        NR = 0;
      else if ((NR = (long)XSS[L1]) > 0)
        {
          /* Set type according to first region */

          WDB[erg + ERG_INTERP] = XSS[L1 + 1 + NR];

          /* Check remainding regions */

          for (n = 0; n < NR; n++)
            if (RDB[erg + ERG_INTERP] != XSS[L1 + 1 + NR + n])
              {
                /* Multiple interpolation regions are rare, and the routines */
                /* are tested only with laws 7/9/11. */

                if ((LAW != 7) && (LAW != 9) && (LAW != 11))
                  Warn(FUNCTION_NAME,
                       "Multiple interpolations: (law %ld mt %ld)",
                       (long)RDB[erg + ERG_INTERP],
                       (long)XSS[L1 + 1 + NR + n]);
                else
                  {
                    /* Set interpolation to -1 to indicate multiple regions */

                    WDB[erg + ERG_INTERP] = -1.0;

                    /* Allocate memory for data */

                    l0 = ReallocMem(DATA_ARRAY, 2*NR + 2);

                    /* Put number of regions */

                    WDB[erg + ERG_NR] = (double)NR + 1;

                    /* Put pointer */

                    WDB[erg + ERG_PTR_INTERP] = (double)l0;

                    /* Put first index */

                    WDB[l0++] = -1.0;

                    /* Read indexes */

                    for (i = 0; i < NR; i++)
                      WDB[l0++] = XSS[L1 + 1 + i] - 1.0;

                    /* Read interpolations */

                    for (i = 0; i < NR; i++)
                      WDB[l0++] = XSS[L1 + 1 + NR + i];

                    /* Put last interpolation */

                    WDB[l0++] = -1.0;
                  }

                /* Break loop */

                break;
              }
        }

      /* Reset data pointer */

      WDB[erg + ERG_PTR_DATA] = NULLPTR;

      /***********************************************************************/

      /***** Read law data ***************************************************/

      switch(LAW)
        {
        default:
          {
            /*****************************************************************/

            /****** Undefined ENDF law ***************************************/

            /* Print error */

            Die(FUNCTION_NAME, "Unsupported law %ld (%s mt %ld)", LAW,
                GetText(nuc + NUCLIDE_PTR_NAME), mt);

            /* Break case (voi olla -nofatal -moodissa) */

            break;

            /*****************************************************************/
          }
        case 1:
          {
            /*****************************************************************/

            /***** Equi-probable energy bins *********************************/

            /* Get number of incident energy points */

            NE = (long)XSS[L1 + 2*NR + 1];

            /* Create energy grid */

            ptr = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L1 + 2*NR + 2],
                                 EG_INTERP_MODE_LIN);

            /* Get number of out-going energy points */

            NP = (long)XSS[L1 + 2*NR + NE + 2];

            /* Allocate memory for data */

            l0 = ReallocMem(DATA_ARRAY, NE*NP + 2);

            /* Set pointer */

            WDB[erg + ERG_PTR_DATA] = (double)l0;

            /* Set pointer to incident energy grid */

            WDB[l0++] = (double)ptr;

            /* Set number of data points */

            WDB[l0++] = (double)(NP - 1);

            /* Read data points */

            for (n = 0; n < NE*NP; n++)
              WDB[l0++] = XSS[L1 + 2*NR + NE + 3 + n];

            /* Check pointer */

            CheckPointer(FUNCTION_NAME, "l0", DATA_ARRAY, l0 - 1);

            /* Perform a few consistency checks */

            l0 = (long)RDB[erg + ERG_PTR_DATA];
            ptr = (long)RDB[l0++];
            NE = (long)RDB[ptr + ENERGY_GRID_NE];

            /* Get number of data points */

            NP = (long)RDB[l0++] + 1;

            /* Check emission energies */

            for (n = 0; n < NE*NP; n++)
              {
                CheckValue(FUNCTION_NAME, "E1", " (law 1)", RDB[l0], 0.0,
                           E_CHECK_MAX);
                l0++;
              }

            /* Exit */

            break;

            /*****************************************************************/
          }
        case 3:
        case 66:
          {
            /*****************************************************************/

            /*****  3 = Level scattering *************************************/
            /***** 66 = N-body phase space distribution **********************/

            /* NOTE: NR is not defined in these laws */

            /* Allocate memory for data */

            l0 = ReallocMem(DATA_ARRAY, 2);

            /* Set pointer */

            WDB[erg + ERG_PTR_DATA] = (double)l0;

            /* Read data points */

            WDB[l0++] = XSS[L1];
            WDB[l0++] = XSS[L1 + 1];

            /* Check values (XSS[L1] may be negative for isomeric states) */

            CheckValue(FUNCTION_NAME, "XSS[L1]", " (law 3/66)", XSS[L1],
                      -10.0, E_CHECK_MAX);

            if (LAW == 3)
              {
                CheckValue(FUNCTION_NAME, "XSS[L1 + 1]", " (law 3)",
                           XSS[L1 + 1], 0.4, 0.994);
              }
            else
              {
                CheckValue(FUNCTION_NAME, "XSS[L1 + 1]", " (law 66)",
                           XSS[L1 + 1], 1.0, 4.0);
              }

            /* Exit */

            break;

            /*****************************************************************/
          }
        case 4:
        case 44:
        case 61:
          {
            /*****************************************************************/

            /*****  4 = Continuous tabular distribution **********************/
            /***** 44 = Kalbach-87 formalism *********************************/
            /***** 61 = like 44 but tabular angular distr. *******************/

            /* Get number of energy points */

            NE = (long)XSS[L1 + 2*NR + 1];

            /* Create energy grid */

            ptr = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L1 + 2*NR + 2],
                                 EG_INTERP_MODE_LIN);

            /* Allocate memory for types */

            l1 = ReallocMem(DATA_ARRAY, NE);

            /* Calculate total data size and put types  */

            sum = 0;
            for(i = 0; i < NE; i++)
              {
                /* Get pointer */

                L2 = LDIS + (long)XSS[L1 + 2*NR + 2 + NE + i] - 2;

                /* Put type */

                WDB[l1 + i] = (double)XSS[L2];

                if (LAW == 4)
                  sum = sum + 3*(long)XSS[L2 + 1] + 1;
                else if (LAW == 44)
                  sum = sum + 5*(long)XSS[L2 + 1] + 1;
                else
                  sum = sum + 4*(long)XSS[L2 + 1] + 1;
              }

            /* Allocate memory for data */

            l0 = ReallocMem(DATA_ARRAY, NE + sum + 2);

            /* Set pointer */

            WDB[erg + ERG_PTR_DATA] = (double)l0;

            /* Set pointer to incident energy grid */

            WDB[l0++] = (double)ptr;

            /* Set pointer to types */

            WDB[l0++] = (double)l1;

            /* Pointer to pointers */

            l2 = l0;
            l0 = l0 + NE;

            /* Loop over initial energies */

            for(i = 0; i < NE; i++)
              {
                /* Set pointer to distribution */

                WDB[l2++] = (double)l0;

                /* Get pointer to data */

                L2 = LDIS + (long)XSS[L1 + 2*NR + 2 + NE + i] - 2;

                /* Get number of energy points */

                NP = (long)XSS[L2 + 1];
                WDB[l0++] = (double)NP;

                /* Read energy values. */

                for(j = 0; j < NP; j++)
                  WDB[l0++] = XSS[L2 + 2 + j];

                /* Read pdf values. */

                for(j = 0; j < NP; j++)
                  WDB[l0++] = XSS[L2 + NP + 2 + j];

                /* Read cdf values. */

                for(j = 0; j < NP; j++)
                  WDB[l0++] = XSS[L2 + 2*NP + 2 + j];

                /* Angular distributions */

                if (LAW == 44)
                  {
                    /* Read precombound fractions. */

                    for(j = 0; j < NP; j++)
                      WDB[l0++] = XSS[L2 + 3*NP + 2 + j];

                    /* Read angular distribution slope values. */

                    for(j = 0; j < NP; j++)
                      WDB[l0++] = XSS[L2 + 4*NP + 2 + j];
                  }
                else if (LAW == 61)
                  {
                    /* Read angular distributions */

                    for(j = 0; j < NP; j++)
                      {
                        if ((long)XSS[L2 + 3*NP + 2 + j] > 0)
                          {
                            /* Tabular angular distribution. Get pointer */

                            L3 = LDIS + (long)XSS[L2 + 3*NP + 2 + j] - 3;

                            /* Allocate memory for data */

                            l1 =
                              ReallocMem(DATA_ARRAY,3*(long)XSS[L3 + 2] + 2);

                            /* Set pointer */

                            WDB[l0++] = (double)l1;

                            /* Set type */

                            WDB[l1++] = XSS[L3 + 1];

                            /* Set number of energy points */

                            NP2 = (long)XSS[L3 + 2];
                            WDB[l1++] = (double)NP2;

                            /* Read cosine values */

                            for (k = 0; k < NP2; k++)
                              WDB[l1++] = XSS[L3 + 3 + k];

                            /* Read PDF values */

                            for (k = 0; k < NP2; k++)
                              WDB[l1++] = XSS[L3 + 3 + NP2 + k];

                            /* Read CDF values */

                            for (k = 0; k < NP2; k++)
                              WDB[l1++] = XSS[L3 + 3 + 2*NP2 + k];
                          }
                        else
                          {
                            /* Isotropic distribution. Set pointer to null */

                            WDB[l0++] = NULLPTR;
                          }
                      }
                  }
              }

            /* Check pointer */

            if (LAW == 61)
              {
                CheckPointer(FUNCTION_NAME, "l1", DATA_ARRAY, l1 - 1);
              }
            else
              {
                CheckPointer(FUNCTION_NAME, "l0", DATA_ARRAY, l0 - 1);
              }

            /* Perform a few consistency checks */

            l0 = (long)RDB[erg + ERG_PTR_DATA];
            ptr = (long)RDB[l0++];
            NE = (long)RDB[ptr + ENERGY_GRID_NE];

            /* Skip type */

            l0++;

            /* Loop over initial energies */

            for (i = 0; i < NE; i++)
              {
                /* Pointer to data */

                l1 = (long)RDB[l0 + i];

                /* Get number of data points */

                NP = (long)RDB[l1++];

                /* Check emission energies */

                for (j = 0; j < NP; j++)
                  {
                    CheckValue(FUNCTION_NAME, "E1", " (law 4/44)", RDB[l1 + j],
                               -5E-3, E_CHECK_MAX);
                  }

                /* Force limiting values */

                for (j = 0; j < NP; j++)
                  {
                    if (RDB[l1 + j] < 0.0)
                      WDB[l1 + j] = 0.0;
                    else if (RDB[l1 + j] > E_CHECK_MAX)
                      WDB[l1 + j] = E_CHECK_MAX;
                  }

                /* NOTE: PDF distribution is not checked because some */
                /* isotopes contain rubbish values etc. */

                /* Check CDF distribution (some files have cumulative */
                /* probabilities exciiding 1.0) */

                /* NOTE: Muutamalla metastabiililla nuklidilla on */
                /* negatiivinen arvo tossa JENDL-4.0 -kirjastossa */
                /* (24.11.2014 / 2.1.22) */

                for (j = 0; j < NP; j++)
                  {
                    CheckValue(FUNCTION_NAME, "CDF", " (law 4/44)",
                               RDB[l1 + 2*NP + j], -2E+4, 1.003);
                  }

                /* Check first and last value of CDF */

                CheckValue(FUNCTION_NAME, "CDF(0)", " (law 4/55)",
                           RDB[l1 + 2*NP], 0.0, 1E-6);

                CheckValue(FUNCTION_NAME, "CDF(NP - 1)", " (law 4/55)",
                           RDB[l1 + 3*NP - 1], 1.0 - 1E-6, 1.0 + 1E-6);

                /* Force last and first values */

                WDB[l1 + 2*NP] = 0.0;
                WDB[l1 + 3*NP - 1] = 1.0;

                /* Force limiting values */

                for (j = 0; j < NP; j++)
                  {
                    if (RDB[l1 + 2*NP + j] < 0.0)
                      WDB[l1 + 2*NP + j] = 0.0;
                    else if (RDB[l1 + 2*NP + j] > 1.0)
                      WDB[l1 + 2*NP + j] = 1.0;
                  }

                if (LAW == 61)
                  {
                    /* Pointer to data */

                    l1 = (long)RDB[l1 + 3*NP];

                    /* Skip type */

                    l1++;

                    /* Read number of points */

                    NP2 =(long)RDB[l1++];

                    /* Check cosine values */

                    for (k = 0; k < NP2; k++)
                      {
                        CheckValue(FUNCTION_NAME, "cos", " (law 61)",
                                   RDB[l1 + k], -1.0, 1.0);
                      }

                    /* Check first and last cosine value */

                    if ((RDB[l1] != -1.0) || (RDB[l1 + NP2 - 1] != 1.0))
                      Die(FUNCTION_NAME,"Error in cosines (mt %ld, law 61)",
                          mt);

                    /* Check CDF values */

                    for (k = 0; k < NP2; k++)
                      {
                        CheckValue(FUNCTION_NAME, "CDF", " (law 61)",
                                   RDB[l1 + 2*NP2 + k], 0.0, 1.0);
                      }

                    /* Check first and last value of CDF */

                    CheckValue(FUNCTION_NAME, "CDF(0)", " (law 4/55)",
                               RDB[l1 + 2*NP2], 0.0, 1E-6);

                    CheckValue(FUNCTION_NAME, "CDF(NP - 1)", " (law 4/55)",
                               RDB[l1 + 3*NP2 - 1], 1.0 - 1E-6, 1.0 + 1E-6);

                    /* Force last and first values */

                    WDB[l1 + 2*NP2] = 0.0;
                    WDB[l1 + 3*NP2 - 1] = 1.0;
                  }
              }

            /* Exit */

            break;

            /*****************************************************************/
          }
        case 7:
        case 9:
        case 11:
          {
            /*****************************************************************/

            /*****  7 = Maxwell fission spectrum *****************************/
            /*****  9 = evaporation spectrum *********************************/
            /***** 11 = Watt fission spectrum ********************************/

            /* Get number of energy points */

            NE = (long)XSS[L1 + 2*NR + 1];

            /* Create energy grid */

            ptr = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L1 + 2*NR + 2],
                                 EG_INTERP_MODE_LIN);

            /* Allocate memory for energy and data */

            if (LAW == 7)
              l0 = ReallocMem(DATA_ARRAY, NE + 2);
            else if (LAW == 9)
              l0 = ReallocMem(DATA_ARRAY, NE + 2);
            else
              l0 = ReallocMem(DATA_ARRAY, 2*NE + 1);

            /* Set pointer */

            WDB[erg + ERG_PTR_DATA] = (double)l0;

            /* Set pointer to incident energy grid */

            WDB[l0++] = (double)ptr;

            /* Restriction energy */

            if ((LAW == 9) || (LAW == 7))
              WDB[l0++] = XSS[L1 + 2*NR + 2 + 2*NE];

            /* Read data points */

            for (n = 0; n < NE; n++)
              WDB[l0++] = XSS[L1 + 2*NR + 2 + NE + n];

            /* Put energy boundaries to interpolations */

            if ((l1 = (long)RDB[erg + ERG_PTR_INTERP]) > 0)
              {
                /* Get number of regions */

                nr = (long)RDB[erg + ERG_NR];

                /* Loop over regions */

                for (n = 0; n < nr; n++)
                  {
                    /* Put first energy */

                    if (n == 0)
                      WDB[l1] = 1E-11;
                    else
                      {
                        /* Get index */

                        i = (long)RDB[l1 + n];

                        /* Check */

                        if ((i < 0) || (i > NE - 1))
                          Die(FUNCTION_NAME,
                              "Invalid interpolation energy index");

                        /* Put energy */

                        WDB[l1 + n] = XSS[L1 + 2*NR + 2 + i];
                      }
                  }
              }

            /* Check Watt spectrum */

            if (LAW == 11)
              {
                /* Update pointer for param b-values */

                L1 = L1 + 2 + 2*NR + 2*NE;

                /* Check that there are no interpolation regions */

                if ((NR != 0) || ((long)XSS[L1] != 0))
                  Die(FUNCTION_NAME, "Interpolation in Watt spectrum (mt %ld)",
                      mt);

                /* Read data points */

                for (n = 0; n < NE; n++)
                  WDB[l0++] = XSS[L1 + 2*NR + 2 + NE + n];
              }

            /* Check pointer */

            CheckPointer(FUNCTION_NAME, "l0", DATA_ARRAY, l0 - 1);

            /* Perform a few consistency checks */

            l0 = (long)RDB[erg + ERG_PTR_DATA];
            ptr = (long)RDB[l0++];
            NE = (long)RDB[ptr + ENERGY_GRID_NE];

            /* Check restriction energy */

            if ((LAW == 9) || (LAW == 7))
              {
                CheckValue(FUNCTION_NAME, "U", " (law 9)", RDB[l0], -500.0,
                           500.0);
                l0++;
              }

            /* Check data points */

            for (n = 0; n < NE; n++)
              {
                CheckValue(FUNCTION_NAME, "D1", " (law 7/9/11)", RDB[l0],
                           1E-6, 10.0);
                l0++;
              }

            /* Check data points for Watt spectrum */

            if (LAW == 11)
              for (n = 0; n < NE; n++)
                {
                  CheckValue(FUNCTION_NAME, "D2", " (law 11)", RDB[l0],
                             1E-6, 5.0);
                  l0++;
                }

            /* Exit */

            break;

            /*****************************************************************/
          }
        case 67:
          {
            /*****************************************************************/

            /***** Laboratory Angle-Energy Law *******************************/

            Die(FUNCTION_NAME, "%s: ENDF Law 67 not supported",
                GetText(nuc + NUCLIDE_PTR_NAME));

            /* Get number of energy points */

            NE = (long)XSS[L1 + 2*NR + 1];

            /* Create energy grid */

            ptr = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L1 + 2*NR + 2],
                                 EG_INTERP_MODE_LIN);

            /* Allocate memory for data */

            l0 = ReallocMem(DATA_ARRAY, NE + 1);

            /* Set pointer */

            WDB[erg + ERG_PTR_DATA] = (double)l0;

            /* Set pointer to incident energy grid */

            WDB[l0++] = (double)ptr;

            /* Loop over distributions */

            for (i = 0; i < NE; i++)
              {
                /* Pointer to data */

                L2 = LDIS + (long)XSS[L1 + 2*NR + 2 + NE + i] - 2;

                /* Number of secondary cosines */

                NMU = (long)XSS[L2 + 1];

                /* Allocate memory for data */

                l1 = ReallocMem(DATA_ARRAY, 2*NMU + 2);

                /* Set pointer */

                WDB[l0++] = (double)l1;

                /* Interpolation scheme for secondary cosines */

                CheckValue(FUNCTION_NAME, "INTMU", "", XSS[L2], 1.0, 2.0);
                WDB[l1++] = XSS[L2];

                /* Store number of cosines */

                WDB[l1++] = (double)NMU;

                /* Secondary cosines */

                for (j = 0; j < NMU; j++)
                  {
                    CheckValue(FUNCTION_NAME, "cos", "", XSS[L2 + 2 + j],
                               -1.0, 1.0);
                    WDB[l1++] = XSS[L2 + 2 + j];
                  }

                /* Loop over data */

                for (j = 0; j < NMU; j++)
                  {
                    /* Pointer to data */

                    L3 = LDIS + (long)XSS[L2 + 2 + NMU + j] - 2;

                    /* Number of secondary energies */

                    NPEP = (long)XSS[L3 + 1];

                    /* Allocate memory for data */

                    l2 = ReallocMem(DATA_ARRAY, 3*NPEP + 2);

                    /* Set pointer */

                    WDB[l1++] = (double)l2;

                    /* Interpolation between secondary energies */

                    CheckValue(FUNCTION_NAME, "INTEP", "", XSS[L3], 1.0, 2.0);
                    WDB[l2++] = XSS[L3];

                    /* Store number of energies */

                    WDB[l2++] = (double)NPEP;

                    /* Energy distribution */

                    for (k = 0; k < NPEP; k++)
                      {
                        CheckValue(FUNCTION_NAME, "E", "", XSS[L3 + 2 + k],
                                   0.0, 1000.0);
                        WDB[l2++] = XSS[L3 + 2 + k];
                      }

                    /* PDF */

                    for (k = 0; k < NPEP; k++)
                      {
                        CheckValue(FUNCTION_NAME, "pdf", "",
                                   XSS[L3 + 2 + NPEP + k], 0.0, INFTY);
                        WDB[l2++] = XSS[L3 + 2 + NPEP + k];
                      }

                    /* CDF */

                    for (k = 0; k < NPEP; k++)
                      {
                        CheckValue(FUNCTION_NAME, "cdf", "",
                                   XSS[L3 + 2 + 2*NPEP + k], 0.0, 1.0);
                        WDB[l2++] = XSS[L3 + 2 + 2*NPEP + k];
                      }

                    /* Check first and last values */

                    if (XSS[L3 + 2 + 2*NPEP] != 0.0)
                      Die (FUNCTION_NAME, "cdf(1) = %E\n",
                           XSS[L3 + 2 + 2*NPEP]);

                    if (XSS[L3 + 2 + 2*NPEP + NPEP - 1] != 1.0)
                      Die (FUNCTION_NAME, "cdf(NPEP) = %E\n",
                           XSS[L3 + 2 + 2*NPEP + NPEP - 1]);
                  }
              }

            /* Exit */

            break;

            /*****************************************************************/
          }
        }

      /***********************************************************************/

      /* Pointer to next law */

      L0 = (long)XSS[L0 - 1];
    }

  /***************************************************************************/

  /***** Energy-dependent neutron yields *************************************/

  /* Get TY */

  if ((TY = (long)fabs(RDB[rea + REACTION_TY])) < 101)
    return;

  /* Get pointer to DWL array */

  if ((JED = JXS[10]) < 0)
    Die(FUNCTION_NAME, "JED < 0");

  /* Get pointer to data */

  if ((KY = JED + TY - 101 - 1) < 0)
    Die(FUNCTION_NAME, "KY < 0");

  /* Get number of interpolation regions (pp. F-30 in MCNP 4C Appendix F) */

  if ((NR = (long)XSS[KY]) > 0)
    Die(FUNCTION_NAME, "Interpolation regions");

  /* Number of energies */

  NE = (long)XSS[KY + 1];
  CheckValue(FUNCTION_NAME, "NE", "", NE, 2, 1000);

  /* Allocate memory for data */

  l0 = ReallocMem(DATA_ARRAY, NE + 2);
  WDB[rea + REACTION_PTR_MULT] = (double)l0;

  /* Put number of energies */

  WDB[l0++] = (double)NE;

  /* Make energy grid */

  ptr = MakeEnergyGrid(NE, 0, 0, -1, &XSS[KY + 2], EG_INTERP_MODE_LIN);
  WDB[l0++] = (double)ptr;

  /* Reset dubious value count */

  i = 0;

  /* Read values */

  for (n = 0; n < NE; n++)
    if ((WDB[l0++] = XSS[KY + 2 + NE + n]) > 25.0)
      {
        /* There are some errors in the TEND-2014 data file, leading to */
        /* very high multiplicities. The values are truncated to 25 to  */
        /* avoid problems. */

        WDB[l0 - 1] = 25.0;

        /* Update dubious value counter */

        i++;
      }

  /* Print warning message */

  if (i > 0)
    Warn(FUNCTION_NAME,
         "Nuclide %s mt %ld has suspiciously high neutron mutiplicity",
         GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT]);

  /***************************************************************************/
}

/*****************************************************************************/
