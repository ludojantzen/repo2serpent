/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processmudistributions.c                       */
/*                                                                           */
/* Created:       2010/01/19 (JLe)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes ENDF angular distributions in ACE data             */
/*                                                                           */
/* Comments: - Adopted from angulardistribution.c in Serpent 1.1.12          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessMuDistributions:"

/*****************************************************************************/

void ProcessMuDistributions(long rea)
{
  long L, L0, INTT, NE, NP, mt, n, i, j, l0, l1, type, erg, ang, ace, JXS[32];
  long ptr, sens, nuc, count, tblr;
  double mu, P, *XSS;
  double LMS[7] = {0};

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to distribution data */

  if ((L = (long)RDB[rea + REACTION_PTR_ANG]) < 0)
    return;

  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  /* Check mt (Some mt's > 100 have the pointer set, as if there */
  /* is angular distribution data, but the values make no sense. */
  /* The same pointers are read in Serpent 1 as well.) */

  if ((mt != 2) && ((mt < 16) || (mt > 91)))
    {
      /* Reset pointer */

      WDB[rea + REACTION_PTR_ANG] = NULLPTR;

      /* Exit subroutine */

      return;
    }

  /* Get pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Get pointer to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];

  /* Read data to JXS array */

  ptr = (long)ACE[ace + ACE_PTR_JXS];

  for (n = 0; n < 32; n++)
    JXS[n] = (long)ACE[ptr++];

  /* Get pointer to XSS array */

  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

  /* Get and check the number of energy points (Table F-12, page F-18). */

  if ((NE = (long)XSS[L - 1]) < 1)
    Die(FUNCTION_NAME, "mt %ld has no energy points in angular distribution",
        mt);

  /* Find the type of first anisotropic distribution */

  type = 0;

  n = 0;
  while ((n < NE) && ((type = (long)XSS[L + NE + n]) == 0))
    n++;

  /* Check if reaction has only isotropic distributions */

  if (type == 0)
    return;

  /* Check number of energies */

  if (NE < 2)
    Die(FUNCTION_NAME, "NE = %ld (< 2)", NE);

  /* Check energies (NOTE: Some data (Ag-110m in JENDL-3.2) contains */
  /* negative values, which are converted to 1E-11 in Serpent 1). */

  for (n = 0; n < NE; n++)
    if ((XSS[L + n] < 1E-12) || (XSS[L + n]) > 1000.0)
      Warn(FUNCTION_NAME, "Invalid energy value %E\n", XSS[L + n]);

  /* Check order and calculate number of coincident points */

  count = 0;
  for (n = 0; n < NE - 1; n++)
    {
      if (XSS[L + n] > XSS[L + n + 1])
        Die(FUNCTION_NAME, "Energy values are not in ascedning order: %E %E",
            XSS[L + n], XSS[L + n + 1]);
      else if (XSS[L + n] == XSS[L + n + 1])
        {
          /* Add to count */

          count++;

          /* Adjust energy (tää on kauhea viritelmä joka yrittää ottaa */
          /* huomioon sen että identtisiä pisteitä voi olla useampia.  */
          /* Välttämättä toi identtisyys ei edes ole ongelma, mutta se */
          /* on korjattu jossain vaiheessa niin oletetaan että sillä   */
          /* on joku vaikutus johonkin -- JLE / 2.1.22 / 22.10.2014).  */

          XSS[L + n] = XSS[L + n] - ((double)(NE - n))*1E-11/((double)NE);

          if (XSS[L + n] < 0.0)
            XSS[L + n] = 0.0;
        }
    }

  /* Check again */

  for (n = 0; n < NE - 1; n++)
    if (XSS[L + n] >= XSS[L + n + 1])
      Warn(FUNCTION_NAME, "Energy values are not in ascedning order: %E %E",
           XSS[L + n], XSS[L + n + 1]);

  /* Check number of coincident points */

#ifdef DEBUG2

  if (count > 0)
    Warn(FUNCTION_NAME, "%ld coincident energy points in mubar data (%s mt %ld)", count, GetText(nuc + NUCLIDE_PTR_NAME), mt);

#endif

  /* Make energy grid */

  erg = MakeEnergyGrid(NE, 0, 0, -1, &XSS[L], EG_INTERP_MODE_LIN);

  /* Allocate memory for distribution */

  WDB[rea + REACTION_PTR_ANG] = 0.0;
  ang = NewItem(rea + REACTION_PTR_ANG, ANG_BLOCK_SIZE);

  /* Put pointer to energy distribution */

  WDB[ang + ANG_PTR_EGRID] = (double)erg;

  /* Check type */

  if (type > 0)
    {
      /***********************************************************************/

      /***** 32 equi-probable cosine bins ************************************/

      /* This won't work with the scattering moment sensitivities currently */
      /* TODO: Might be possible? */

      if (((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
          ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM))
        Die (FUNCTION_NAME, "Legendre sensitivity error\n"
             "32 equi-probable cosine bin distribution");

     /* Set type */

      WDB[ang + ANG_TYPE] = (double)ANG_TYPE_EQUIBIN;

      /* Allocate memory for data */

      l0 = ReallocMem(DATA_ARRAY, 33*NE);

      /* Set pointer */

      WDB[ang + ANG_PTR_D0] = (double)l0;

      /* Set number of bins */

      WDB[ang + ANG_BINS] = 32.0;

      /* Loop over energies */

      for (n = 0; n < NE; n++)
        {
          /* Get pointer to first bin array */

          L0 = (long)XSS[L + NE + n];

          if (L0 == 0)
            {
              /* Isotropic distribution, calculate values */

              for (i = 0; i < 33; i++)
                WDB[l0++] = (double)i/16.0 - 1.0;
            }
          else if (L0 > 0)
            {
              /* Anisotropic distribution, get values from bins */

              for (i = 0; i < 33; i++)
                {
                  /* Get value */

                  mu = XSS[JXS[8] + L0 - 2 + i];
                  CheckValue(FUNCTION_NAME, "mu", " (mu)", mu, -1.001, 1.001);

                  /* Truncate */

                  if (mu > 1.0)
                    mu = 1.0;
                  else if (mu < -1.0)
                    mu = -1.0;

                  /* Set value */

                  WDB[l0++] = mu;
                }

              /* Check first and last value */

              if (fabs(RDB[l0 - 33] + 1.0) > 5.5E-3)
                Warn(FUNCTION_NAME,
                     "Invalid boundary value %E in distribution (mt = %ld)",
                     RDB[l0 - 33], mt);
              else if (fabs(RDB[l0 - 1] - 1.0) > 4.1E-2)
                Warn(FUNCTION_NAME,
                     "Invalid boundary value %E in distribution (mt = %ld)",
                     RDB[l0 - 1], mt);
              else
                {
                  /* Make sure values are ok. */

                  WDB[l0 - 33] = -1.0;
                  WDB[l0 - 1] = 1.0;

                }
            }
          else
            Die(FUNCTION_NAME, "Mixed angular distributions (mt %ld)", mt);
        }

      /* Check that memory and pointers match */

      CheckPointer(FUNCTION_NAME, "angular mu", DATA_ARRAY, l0 - 1);

      /***********************************************************************/
    }
  else if (type < 0)
    {
      /***** Tabular angular distribution ************************************/

     /* Set type */

      WDB[ang + ANG_TYPE] = (double)ANG_TYPE_TABULAR;

      /* Allocate memory for data pointers */

      l0 = ReallocMem(DATA_ARRAY, NE);

      /* Set pointer */

      WDB[ang + ANG_PTR_D0] = (double)l0;

     /* Get pointer to data */

      L0 = L + 2*NE + 2;

      /* Set data type */

      INTT = (long)XSS[L0 - 2];
      WDB[ang + ANG_INTT] = (double)INTT;

      /* Loop over energies */

      for (i = 0; i < NE; i++)
        {
          /* Check type and get number of points */

          if ((long)XSS[L0 - 2] != INTT)
            Die(FUNCTION_NAME, "mixed mu interpolations (mt %ld)", mt);

          /* Get number of values */

          NP = (long)XSS[L0 - 1];

          /* Allocate memory for data. For Sens, add 7 doubles at the end */
          /* for precalculated Legendre moments */

          if ((long)RDB[DATA_SENS_MODE] == SENS_MODE_NONE)
            l1 = ReallocMem(DATA_ARRAY, 3*NP + 1);
          else
            l1 = ReallocMem(DATA_ARRAY, 3*NP + 1 + 7);

          /* Set pointer */

          WDB[l0++] = (double)l1;

          /* Set number of points */

          WDB[l1++] = (double)NP;

          /* Store pointer to current location of l1 */

          tblr = l1;

          /* Loop over cosine values */

          for (j = 0; j < NP; j++)
            {
              /* Get value */

              mu = XSS[L0 + j];
              CheckValue(FUNCTION_NAME, "mu", " (mu)", mu, -1.001, 1.001);

              /* Truncate */

              if (mu > 1.0)
                mu = 1.0;
              else if (mu < -1.0)
                mu = -1.0;

              /* Set value */

              WDB[l1++] = mu;
            }

          /* Check first and last cosine values */

          /* NOTE: U-237 in kth endfb-6 has a large difference in last value */

          if ((RDB[l1 - NP] != -1.0) || (fabs(RDB[l1 - 1] - 1.0) > 5E-2))
            Die(FUNCTION_NAME, "Error in angular distribution (mt %ld)",  mt);
          else
            {
              /* Adjust these values just to be sure */

              WDB[l1 - NP] = -1.0;
              WDB[l1 - 1] = 1.0;
            }

          /* Loop over PDF  values (no check, values may not be normalised). */

          for (j = 0; j < NP; j++)
            WDB[l1++] = XSS[L0 + NP + j];

          /* Loop over CDF  values */

          for (j = 0; j < NP; j++)
            {
              /* Get value */

              P = XSS[L0 + 2*NP + j];
              CheckValue(FUNCTION_NAME, "CDF", "", P, 0.0, 1.0);

              /* Set value */

              WDB[l1++] = P;
            }

          /* Check first and last cumulative value (NOTE: NJOY-2012 antaa  */
          /* JEFF-3.2 H1:lle suurimmaksi arvoksi 0.9999999 ajettaessa HP:n */
          /* läppärillä, mutta samalla inputilla tasan 1.0 espvm5m100:ssa) */

          if ((RDB[l1 - NP] != 0.0) || (RDB[l1 - 1] != 1.0))
            Die(FUNCTION_NAME, "Error in angular distribution (mt %ld)", mt);

          /* Check that memory and pointers match */

          CheckPointer(FUNCTION_NAME, "angular mu", DATA_ARRAY, l1 - 1);

          /* Update pointer */

          L0 = L0 + 3*NP + 2;

          /********* Sens ******************/

          if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
            {
              /* Precalculate Legendre moments for this energy point */

              CalculateLegendreMoments(NP, tblr, LMS);

              /* Store Legendre moments for this energy point */

              for (j = 0; j < 7; j++)
                WDB[l1 + j] = LMS[j];
            }

          /********************************/
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid angular distribution type %ld (mt %ld)",
        type, mt);
}

/*****************************************************************************/
