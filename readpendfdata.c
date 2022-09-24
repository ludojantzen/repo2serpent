/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readpendfdata.c                                */
/*                                                                           */
/* Created:       2018/03/22 (RTu)                                           */
/* Last modified: 2019/11/08 (RTu)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads and processes additional PENDF data                    */
/*                                                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadPendfData:"

/*****************************************************************************/

void ReadPendfData(long nuc, FILE *fp, long pos, long mt)
{
  long ace, ptr, rea, n, NES, LFC, NFC, IFC, NR, erg, NPLY, np, i;
  long j, NBIN, NUNR, skip, N, M, IFF, LSSF;
  long data_ptr, comp_ptr, type_ptr, tab_ptr, ptr1;
  double *XSS, *NXS, *JXS, *xs0, *E0, *tot_xs, *kerma, *E;
  double f, val;
  char dummy[MAX_STR], *eof;

  /* Set file position to the end of ace data */

  if (fseek(fp, pos, SEEK_SET) != 0)
    Die(FUNCTION_NAME, "fseek error\n");

  /* Reset dummy */

  dummy[71] = '\0';

  /* Loop until given MT */

  do
    {
      eof = fgets(dummy, 82, fp);
      dummy[75] = '\0';
    }
  while ((atoi(&dummy[71]) != mt) && (eof != NULL));

  if (mt == 1458)
    {
      if (eof == NULL)
        {
          Warn(FUNCTION_NAME, "No MT458 data for nuclide %s and edep > 0",
               GetText(nuc + NUCLIDE_PTR_NAME));
        }
      else
        {
          /* Create new structure for MT 458  data */

          data_ptr = ReallocMem(DATA_ARRAY, FISSE_DATA_BLOCK_SIZE);
          WDB[nuc + NUCLIDE_PTR_FISSE_DATA] = (double)data_ptr;

          /* Allocate memory for component array */

          comp_ptr = ReallocMem(DATA_ARRAY, FISSE_COMP_BLOCK_SIZE);
          WDB[data_ptr + FISSE_DATA_COMP] = (double)comp_ptr;

          /* Read LFC and NFC flags */

          LFC = ENDFColI(4, dummy);
          CheckValue(FUNCTION_NAME, "LFC", "", LFC, 0, 1);
          WDB[data_ptr + FISSE_DATA_LFC] = (double)LFC;

          NFC = ENDFColI(6, dummy);
          CheckValue(FUNCTION_NAME, "NFC", "", NFC, 0, FISSE_COMP_BLOCK_SIZE);

          /* Next line */

          ENDFNewLine(FUNCTION_NAME, dummy, fp);

          /* Read NPLY */

          NPLY = ENDFColI(4, dummy);
          CheckValue(FUNCTION_NAME, "NPLY", "", NPLY, 0, 10);
          WDB[data_ptr + FISSE_DATA_NPLY] = (double)NPLY;

          /* Check NPLY */

          if (NPLY == 0)
            {
              /* Check N1 and N2 */

              if ((ENDFColI(5, dummy)) != 18)
                Die(FUNCTION_NAME, "N1 != 18");
              else if ((ENDFColI(6, dummy)) != 9)
                Die(FUNCTION_NAME, "N2 != 9");

              /* Read components of fission energy */

              ENDFNewLine(FUNCTION_NAME, dummy, fp);

              /* EFR */

              val = 1E-6*ENDFColF(1, dummy);
              CheckValue(FUNCTION_NAME, "EFR", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_EFR] = val;

              /* ENP */

              val = 1E-6*ENDFColF(3, dummy);
              CheckValue(FUNCTION_NAME, "ENP", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_ENP] = val;

              /* END */

              val = 1E-6*ENDFColF(5, dummy);
              CheckValue(FUNCTION_NAME, "END", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_END] = val;

              ENDFNewLine(FUNCTION_NAME, dummy, fp);

              /* EGP */

              val = 1E-6*ENDFColF(1, dummy);
              CheckValue(FUNCTION_NAME, "EGP", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_EGP] = val;

              /* EGD */

              val = 1E-6*ENDFColF(3, dummy);
              CheckValue(FUNCTION_NAME, "EGD", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_EGD] = val;

              /* EB */

              val = 1E-6*ENDFColF(5, dummy);
              CheckValue(FUNCTION_NAME, "EB", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_EB] = val;

              ENDFNewLine(FUNCTION_NAME, dummy, fp);

              /* ENU */

              val = 1E-6*ENDFColF(1, dummy);
              CheckValue(FUNCTION_NAME, "ENU", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_ENU] = val;

              /* ER */

              val = 1E-6*ENDFColF(3, dummy);
              CheckValue(FUNCTION_NAME, "ER", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_ER] = val;

              /* ET */

              val = 1E-6*ENDFColF(5, dummy);
              CheckValue(FUNCTION_NAME, "ET", "", val, 0.0, 300);
              WDB[comp_ptr + FISSE_COMP_ET] = val;

              /* Tabular data */

              if (LFC == 1)
                {
                  if (NFC == 0)
                    Die(FUNCTION_NAME, "LFC = 1 but no TAB1 records available\n");

                  /* Allocate memory for component type array */

                  type_ptr = ReallocMem(DATA_ARRAY, FISSE_COMP_BLOCK_SIZE);
                  WDB[data_ptr + FISSE_DATA_COMP_TYPES] = (double)type_ptr;

                  /* Initialize all component types to Sher and Beck */

                  for (n = 0; n < FISSE_COMP_BLOCK_SIZE; n++)
                    WDB[type_ptr + n] = (double)FISSE_TYPE_SHER_BECK;

                  /* Read tabular component data */

                  for (n = 0; n < NFC; n++)
                    {
                      /* Next line */

                      ENDFNewLine(FUNCTION_NAME, dummy, fp);

                      /* Get IFC, NR and NP */

                      IFC = ENDFColI(4, dummy);
                      NR = ENDFColI(5, dummy);
                      np = ENDFColI(6, dummy);

                      /* Check IFC */

                      CheckValue(FUNCTION_NAME, "IFC", "", IFC, 1, FISSE_COMP_BLOCK_SIZE);

                      /* Change component type to tabular */

                      WDB[type_ptr + IFC - 1] = (double)FISSE_TYPE_TAB;

                      /* Allocate memory for tab block */

                      tab_ptr = ReallocMem(DATA_ARRAY, FISSE_TAB_BLOCK_SIZE);
                      WDB[comp_ptr + IFC - 1] = (double)tab_ptr;

                      /* Store NP */

                      CheckValue(FUNCTION_NAME, "NP", "", np, 2, INFTY);
                      WDB[tab_ptr + FISSE_TAB_NP] = (double)np;

                      /* Store NR */

                      CheckValue(FUNCTION_NAME, "NR", "", NR, 1, 100);
                      WDB[tab_ptr + FISSE_TAB_NR] = (double)NR;

                      /* Read interpolation ranges */

                      if (NR == 1)
                        {
                          ENDFNewLine(FUNCTION_NAME, dummy, fp);
                          WDB[tab_ptr + FISSE_TAB_NBT] = (double)ENDFColI(1, dummy);
                          WDB[tab_ptr + FISSE_TAB_INT] = (double)ENDFColI(2, dummy);
                        }
                      else
                        {
                          ptr = ReallocMem(DATA_ARRAY, NR);
                          WDB[tab_ptr + FISSE_TAB_NBT] = (double)ptr;

                          ptr1 = ReallocMem(DATA_ARRAY, NR);
                          WDB[tab_ptr + FISSE_TAB_INT] = (double)ptr1;

                          ENDFNewLine(FUNCTION_NAME, dummy, fp);

                          i = 1;
                          for (j = 0; j < NR; j++)
                            {
                              WDB[ptr + j]  = (double)ENDFColI(i++, dummy);
                              WDB[ptr1 + j] = (double)ENDFColI(i++, dummy);

                              if ((i > 6) && (j < NR - 1))
                                {
                                  ENDFNewLine(FUNCTION_NAME, dummy, fp);
                                  i = 1;
                                }
                            }
                        }

                      /* Incident neutron energy/energy deposition pairs */

                      erg = ReallocMem(DATA_ARRAY, np);
                      WDB[tab_ptr + FISSE_TAB_ERG] = (double)erg;

                      ptr = ReallocMem(DATA_ARRAY, np);
                      WDB[tab_ptr + FISSE_TAB_EDEP] = (double)ptr;

                      /* Read energy-cross section pairs */

                      ENDFNewLine(FUNCTION_NAME, dummy, fp);

                      i = 1;
                      for (j = 0; j < np; j++)
                        {
                          WDB[erg + j] = 1E-6*ENDFColF(i++, dummy);
                          WDB[ptr + j] = 1E-6*ENDFColF(i++, dummy);

                          if ((i > 6) && (j < np - 1))
                            {
                              ENDFNewLine(FUNCTION_NAME, dummy, fp);
                              i = 1;
                            }
                        }
                    }
                }
            }
          else
            {
              /* Polynomial data, allocate memory for coefficient vectors */

              for (i = 0; i < FISSE_COMP_BLOCK_SIZE; i++)
                WDB[comp_ptr + i] = (double)ReallocMem(DATA_ARRAY, NPLY + 1);

              /* Loop over data */

              for (n = 0; n < NPLY + 1; n++)
                {
                  /* Read data */

                  ENDFNewLine(FUNCTION_NAME, dummy, fp);

                  /* EFR */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(1, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_EFR];
                  WDB[ptr + n] = val;

                  /* ENP */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(3, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_ENP];
                  WDB[ptr + n] = val;

                  /* END */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(5, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_END];
                  WDB[ptr + n] = val;

                  ENDFNewLine(FUNCTION_NAME, dummy, fp);

                  /* EGP */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(1, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_EGP];
                  WDB[ptr + n] = val;

                  /* EGD */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(3, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_EGD];
                  WDB[ptr + n] = val;

                  /* EB */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(5, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_EB];
                  WDB[ptr + n] = val;

                  ENDFNewLine(FUNCTION_NAME, dummy, fp);

                  /* ENU */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(1, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_ENU];
                  WDB[ptr + n] = val;

                  /* ER */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(3, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_ER];
                  WDB[ptr + n] = val;

                  /* ET */

                  val = 1E-6*pow(1E6, (double)n)*ENDFColF(5, dummy);
                  ptr = (long)RDB[comp_ptr + FISSE_COMP_ET];
                  WDB[ptr + n] = val;
                }
            }
        }
    }
  else
    {
      /* Get pointers to ACE data */

      ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
      CheckPointer(FUNCTION_NAME, "ace", ACE_ARRAY, ace);

      XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];
      NXS = &ACE[(long)ACE[ace + ACE_PTR_NXS]];
      JXS = &ACE[(long)ACE[ace + ACE_PTR_JXS]];

      if (mt == 3301 || mt == 3318 || mt == 3319)
        {
          /* Check EOF */

          if (eof == NULL)
            Die(FUNCTION_NAME, "No KERMA data for nuclide %s mt = %d",
                GetText(nuc + NUCLIDE_PTR_NAME), mt);

          /* Next line */

          ENDFNewLine(FUNCTION_NAME, dummy, fp);

          /* Get number of energy-cross section pairs */

          np = ENDFColI(6, dummy);

          /* Allocate memory for temporary energy and cross section arrays */

          xs0 =  (double *)Mem(MEM_ALLOC, np, sizeof(double));
          E0  =  (double *)Mem(MEM_ALLOC, np, sizeof(double));

          /* Read energy-cross section pairs */

          ENDFNewLine(FUNCTION_NAME, dummy, fp);
          ENDFNewLine(FUNCTION_NAME, dummy, fp);

          i = 1;
          for (j = 0; j < np; j++)
            {
              E0[j]  = ENDFColF(i++, dummy)/1E6;
              xs0[j] = ENDFColF(i++, dummy)/1E6;

              if ((i > 6) && (j < np - 1))
                {
                  ENDFNewLine(FUNCTION_NAME, dummy, fp);
                  i = 1;
                }
            }

          /* The PENDF kerma data should be on the same energy grid as the */
          /* ACE data so the following linear interpolation is just a precaution. */
          /* For some reason the KERMA data for U-235 processed from ENDF/B-VII.1 */
          /* contains one more energy point than the ACE data at E ~ 1.09 MeV */

          /* Number of energy points in ace data */

          NES = (long)NXS[2];

          /* Energy grid in ace data */

          E = &XSS[(long)JXS[0] - 1];

          /* Check that boundary values are equal. They might not be exactly equal since */
          /* the values are converted from different strings with different functions. */
          /* Make them equal so InterpolateData works properly at boundaries. */

          if (fabs(E0[0] - E[0]) < 1E-9*E[0])
            E0[0] = E[0];

          if (fabs(E0[np - 1] - E[NES - 1]) < 1E-9*E[NES - 1])
            E0[np - 1] = E[NES - 1];

          /* Allocate memory for temporary array */

          kerma = (double *)Mem(MEM_ALLOC, NES, sizeof(double));

          /* Linear interpolation */

          n = InterpolateData(E, kerma, NES, E0, xs0,
                              np, 0, NULL, NULL, YES);

          if (n > 0)
            Warn(FUNCTION_NAME,
                 "%ld negative xs points in KERMA data (%s)",
                 n, GetText(nuc + NUCLIDE_PTR_NAME));

          if (mt == 3301)
            {
              /* Total cross section array */

              tot_xs = &XSS[(long)JXS[0] - 1 + NES];

              /* KERMA reaction pointer */

              rea = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS];

              /* Pointer to XS data */

              ptr = (long)RDB[rea + REACTION_PTR_XS];

              /* Convert non-local KERMA to heating factors and */
              /* replace local KERMA in the XSS array */

              n = 0;

              for (i = 0; i < NES; i++)
                {
                  if (fabs(tot_xs[i]) < ZERO)
                    {
                      n++;
                      XSS[ptr + i] = 0.0;
                    }
                  else
                    XSS[ptr + i] = kerma[i]/tot_xs[i];
                }

              if (n > 0)
                Warn(FUNCTION_NAME,
                     "%ld zero total cross sections in data (%s)",
                     n, GetText(nuc + NUCLIDE_PTR_NAME));
            }
          else
            {
              NFKERMA(nuc, kerma);
            }

          /* Free temporary arrays */

          Mem(MEM_FREE, xs0);
          Mem(MEM_FREE, E0);
          Mem(MEM_FREE, kerma);
        }
      else if (mt == 2153)
        {
          /* Pointer to UNR block */

          if ((ptr = (long)JXS[22] - 1) < 0)
            Die(FUNCTION_NAME, "Ures flag set but ACE block not found");

          /* Number of incident energies */

          N = (long)XSS[ptr];

          /* Number of probabilities */

          M = (long)XSS[ptr + 1];

          /* Factors flag */

          IFF = (long)XSS[ptr + 5];

          /* Check EOF */

          if (eof == NULL)
            Die(FUNCTION_NAME, "No non-local KERMA URES data for nuclide %s",
                GetText(nuc + NUCLIDE_PTR_NAME));

          /* Get NBIN, LSSF and NUNR */

          NBIN = ENDFColI(6, dummy);

          if (M != NBIN)
            Die(FUNCTION_NAME, "Mismatch in the number of probabilities "
                "in KERMA URES data for nuclide %s", GetText(nuc + NUCLIDE_PTR_NAME));

          ENDFNewLine(FUNCTION_NAME, dummy, fp);

          LSSF = ENDFColI(3, dummy);

          if (IFF != LSSF)
            Die(FUNCTION_NAME, "Mismatch in the factors flag "
                "in KERMA URES data for nuclide %s", GetText(nuc + NUCLIDE_PTR_NAME));

          NUNR = ENDFColI(6, dummy);

          if (N != NUNR)
            Die(FUNCTION_NAME, "Mismatch in the number of energies "
                "in KERMA URES data for nuclide %s", GetText(nuc + NUCLIDE_PTR_NAME));

          /* Number of data values per energy */

          np = 1 + 6*NBIN;

          /* Number of skipped data values per energy */

          skip = 1 + 5*NBIN;

          /* Conversion to MeV/reaction if the read values are */
          /* heating values in eV/reaction */

          if (LSSF == 0)
            f = 1E-6;
          else
            f = 1.0;

          /* Update pointer to the start of probability tables */

          ptr = ptr + 6 + N;

          /* Allocate memory for temporary array */

          xs0 = (double *)Mem(MEM_ALLOC, NBIN, sizeof(double));

          /* Read KERMA URES cross sections or factors */

          i = 1;

          ENDFNewLine(FUNCTION_NAME, dummy, fp);

          for (n = 0; n < NUNR; n++)
            {
              for (j = 0; j < np; j++)
                {
                  /* Skip everything else than KERMA data */

                  if (j > (skip - 1))
                    {
                      xs0[j - skip] = f*ENDFColF(i++, dummy);
                    }
                  else
                    i++;

                  if ((i > 6) && !(n == (NUNR - 1) && j == (np - 1)))
                    {
                      ENDFNewLine(FUNCTION_NAME, dummy, fp);
                      i = 1;
                    }
                }
              /* Replace data in the XSS array */

              memcpy(&XSS[ptr + n*6*M + 5*M], xs0, NBIN*sizeof(double));
            }

          /* Free temporary array */

          Mem(MEM_FREE, xs0);
        }
      else
        Die(FUNCTION_NAME, "Unknown mt (%d)", mt);
    }
}

/*****************************************************************************/
