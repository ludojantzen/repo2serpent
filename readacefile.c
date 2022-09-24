/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readacefile.c                                  */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2020/05/14 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Reads data from ACE format cross section library file        */
/*              into ACE data block                                          */
/*                                                                           */
/* Comments: - The new ACE format is described in:                           */
/*                                                                           */
/*             J. Conlin, et al. "Updating the Format of ACE Data Tables."   */
/*             Trans. Am. Nucl. Soc. 107 (2012) 631-633.                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadACEFile:"

/*****************************************************************************/

void ReadACEFile(long nuc)
{
  long ace, ptr, rea, n, sz, NXS[16], JXS[32], NES, L0, L, NTR, nr, mt, nc, I0;
  long pos, fiss, edep, ures, nfix, IZ[16];
  double *XSS, awr, Emin, Emax, T;
  char HZ1[MAX_STR], HZ2[MAX_STR], dummy[MAX_STR], name[MAX_STR];
  char file[MAX_STR], date[MAX_STR];
  FILE *fp;

  /* Check nuclide type */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
    if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TRANSMU_DATA))
      Die(FUNCTION_NAME, "Decay data");

  /* Get pointer to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
  CheckPointer(FUNCTION_NAME, "ace", ACE_ARRAY, ace);

  /* Put name */

  WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_NAME];
  strcpy(name, GetText(DATA_DUMMY));

  /* Get file name */

  WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_FILE];
  strcpy(file, GetText(DATA_DUMMY));

  /* Test format */

  TestDOSFile(GetText(DATA_DUMMY));

  /* Open file for writing */

  fp = OpenDataFile(DATA_DUMMY, "ACE data file");

  /***************************************************************************/

  /***** Read data ***********************************************************/

  /* Reset HZ */

  *HZ1 = '\0';
  *HZ2 = '\0';

  /* Read ZAID and data */

  while (fscanf(fp, "%s", HZ1) != EOF)
    {
      /* Check for new format (assuming here that the character string  */
      /* is '2.0.0' -- this is something that may need to be checked in */
      /* the future). */

      if (!strcmp(HZ1, "2.0.0") || !strcmp(HZ1, "2.0.1"))
        {
          /* New format, read reminder of line */

          if (fgets(dummy, 81, fp) == NULL)
            Die(FUNCTION_NAME, "fgets error");

          /* Get new HZ */

          sscanf(dummy, "%s", HZ2);

          /* Next line */

          if (fgets(dummy, 81, fp) == NULL)
            Die(FUNCTION_NAME, "fgets error");

          /* Reset number of comment lines */

          nc = 0;

          /* Get atomic weight ratio, temperature, date and number of */
          /* comment lines (new ACE format). Only the last entry is   */
          /* used here. */

          sscanf(dummy, "%lf %lf %s %ld", &awr, &T, date, &nc);
          CheckValue(FUNCTION_NAME, "nc", "", nc, 2, 10000);

          /* Skip comment lines */

          for (n = 0; n < nc - 2; n++)
            if (fgets(dummy, 82, fp) == NULL)
              Die(FUNCTION_NAME, "fgets error");

          /* The next line is the first line of 'old-style' ACE, get ZAID */

          if (fscanf(fp, "%s", HZ1) == EOF)
            Die(FUNCTION_NAME, "fscanf error");
        }

      /* Get atomic weight ratio. Temperature is the next entry after */
      /* AWR, but it is not used. The value is taken from the directory */
      /* file. */

      if (fgets(dummy, 81, fp) == NULL)
        Die(FUNCTION_NAME, "fgets error");

      sscanf(dummy, "%lf", &awr);

      /* Skip comment line */

      if (fgets(dummy, 81, fp) == NULL)
        Die(FUNCTION_NAME, "fgets error");

      /* Preserve decay awr */

      if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
        {
          /* Put atomic weight ratio */

          WDB[nuc + NUCLIDE_AWR] = awr;

          /* Use value read from directory file for atomic weight */

          WDB[nuc + NUCLIDE_AW] = ACE[ace + ACE_AW];
        }

      /* Read 16 IZ AW pairs (used with S(a,b) data) */

      for (n = 0; n < 16; n++)
        if (fscanf(fp, "%ld %s", &IZ[n], dummy) == EOF)
          Die(FUNCTION_NAME, "fscanf error");

      /* Read the NXS array */

      for (n = 0; n < 16; n++)
        if (fscanf(fp, "%ld", &NXS[n]) == EOF)
                Die(FUNCTION_NAME, "Error in NXS array (%s)",
              GetText(nuc + NUCLIDE_PTR_NAME));

      /* Get data size (NXS[0]) */

      if ((sz = NXS[0]) < 10)
        Die(FUNCTION_NAME, "Error in ACE file");

      /* Read the JXS array */

      for (n = 0; n < 32; n++)
        if (fscanf(fp, "%ld", &JXS[n]) == EOF)
                Die(FUNCTION_NAME, "Error in JXS array (%s)",
              GetText(nuc + NUCLIDE_PTR_NAME));

      /* Compare ZAIDs */

      if ((!strcmp(HZ1, name)) || (!strcmp(HZ2, name)))
        {
          /* Allocate memory for NXS array */

          ptr = ReallocMem(ACE_ARRAY, 16);
          ACE[ace + ACE_PTR_NXS] = (double)ptr;

          /* Copy data */

          for (n = 0; n < 16; n++)
            ACE[ptr++] = (double)NXS[n];

          /* Allocate memory for JXS array */

          ptr = ReallocMem(ACE_ARRAY, 32);
          ACE[ace + ACE_PTR_JXS] = (double)ptr;

          /* Copy data */

          for (n = 0; n < 32; n++)
            ACE[ptr++] = (double)JXS[n];

          /* Allocate memory for XSS array */

          ptr = ReallocMem(ACE_ARRAY, sz);

          /* Set pointer */

          ACE[ace + ACE_PTR_XSS] = (double)ptr;

          /* Read data */

          for (n = 0; n < sz; n++)
            if (fscanf(fp, "%lf", &ACE[ptr + n]) == EOF)
              {
                /* Print warning */

                Warn(FUNCTION_NAME, "Error in XSS array (%s)",
                     GetText(nuc + NUCLIDE_PTR_NAME));

                /* Break */

                break;
              }

          /* Break loop */

          break;
        }


      /* Count number of full lines */

      n = (long)((double)sz/4);
      if (!(sz % 4))
        n = n - 1;

      /* Move to last row */

      fseek(fp, n*81 + 1, SEEK_CUR);

      /* NOTE: Tämä oli ennen 20.11.2019 (2.1.32). Toimii muuten, mutta */
      /* feilaa MCNP5:n sab2002 kirjaston jollain nuklideilla */

      /*
      fseek(fp, 81*sz/4, SEEK_CUR);
      */

      /* Read last row to dummy variable */

      if (fgets(dummy, 81, fp) == NULL)
        Die(FUNCTION_NAME, "fgets error");
    }

  /* Check that data was found */

  if ((strcmp(HZ1, name)) && (strcmp(HZ2, name)))
    Die(FUNCTION_NAME, "Unable to find isotope %s in file %s", name, file);

  /* Pointer to XSS array */

  ptr = (long)ACE[ace + ACE_PTR_XSS];
  XSS = &ACE[ptr];

  /* Set ures energy boundaries */

  if ((L = JXS[22] - 1) > 0)
    {
      /* Get number of energy points */

      NES = (long)XSS[L];

      /* Set minimum and maximum energies */

      WDB[nuc + NUCLIDE_URES_EMIN] = XSS[L + 6];
      WDB[nuc + NUCLIDE_URES_EMAX] = XSS[L + 6 + NES - 1];
    }
  else
    {
      /* Reset boundaries */

      WDB[nuc + NUCLIDE_URES_EMIN] = INFTY;
      WDB[nuc + NUCLIDE_URES_EMAX] = -INFTY;
    }

  /**************************************************************************/

  /***** Add reaction channels for transport data ***************************/

  /* Check type */

  if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) ||
      ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DBRC))
    {
      /**********************************************************************/

      /***** Transport nuclide **********************************************/

      /* Number of energy points */

      NES = NXS[2];

      /* Put pointer to nuclide energy grid and number of points */

      WDB[nuc + NUCLIDE_PTR_EGRID] = (double)(JXS[0] - 1);
      WDB[nuc + NUCLIDE_EGRID_NE] = (double)NES;

      /* Get number of reactions (minus elastic scattering). Include */
      /* only elastic scattering for DBRC nuclides. */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DBRC)
        NTR = 0;
      else
        NTR = NXS[3];

      /* Compare minimum and maximum value to XS limits */

      if (XSS[JXS[0] - 1] < RDB[DATA_NEUTRON_XS_EMIN])
        WDB[DATA_NEUTRON_XS_EMIN] = XSS[JXS[0] - 1];

      if (XSS[JXS[0] + NES - 2] > RDB[DATA_NEUTRON_XS_EMAX])
        WDB[DATA_NEUTRON_XS_EMAX] = XSS[JXS[0] + NES - 2];

      /* Number of delayed neutron precursor groups */

      WDB[nuc + NUCLIDE_ACE_PREC_GROUPS] = (double)NXS[7];

      /* Add elastic scattering */

      rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

      /* Put nuclide pointer */

      WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

      /* Reaction index (used for energy distributions) */

      WDB[rea + REACTION_NR] = -1.0;

      /* Put type and MT */

      WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
      WDB[rea + REACTION_MT] = 2.0;

      /* Put awr */

      WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

      /* Set minimum and maximum energy */

      WDB[rea + REACTION_EMIN] = XSS[JXS[0] - 1];
      WDB[rea + REACTION_EMAX] = XSS[JXS[0] + NES - 2];

      WDB[nuc + NUCLIDE_EMIN] = XSS[JXS[0] - 1];
      WDB[nuc + NUCLIDE_EMAX] = XSS[JXS[0] + NES - 2];

      /* Store number of energy points, pointer to grid and */
      /* index to first point */

      WDB[rea + REACTION_PTR_EGRID] = (double)(JXS[0] - 1);
      WDB[rea + REACTION_XS_NE] = (double)NES;
      WDB[rea + REACTION_XS_I0] = 0.0;

      /* Store pointer to XS data */

      WDB[rea + REACTION_PTR_XS] = (double)(JXS[0] - 1 + 3*NES);

      /* Set ures energy boundaries */

      WDB[rea + REACTION_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];
      WDB[rea + REACTION_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMIN];

      /* Set Q-value */

      WDB[rea + REACTION_Q] = 0.0;

      /* Multiplication and frame of reference */

      WDB[rea + REACTION_TY] = -1.0;
      WDB[rea + REACTION_WGT_F] = 1.0;

      /* Set branching fraction to 1.0 */

      WDB[rea + REACTION_BR] = 1.0;

      /* Set interpolation mode */

      WDB[rea + REACTION_ITP] = 0.0;

      /* Pointer to angular distribution (Tables F-11 and F-12) */

      if ((L = (long)XSS[JXS[7] - 1]) > 0)
        WDB[rea + REACTION_PTR_ANG] = (double)(L - 1 + JXS[8]);
      else
        WDB[rea + REACTION_PTR_ANG] = NULLPTR;

      /* Include heat production */

      if ((long)RDB[DATA_INCLUDE_HEAT_PROD_XS] == YES)
        {
          rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);
          WDB[nuc + NUCLIDE_PTR_HEATPRODXS] = (double)rea;

          /* Put nuclide pointer */

          WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

          /* Put type and MT */

          WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;
          WDB[rea + REACTION_MT] = 301.0;

          /* Put awr */

          WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

          /* Set minimum and maximum energy */

          WDB[rea + REACTION_EMIN] = XSS[JXS[0] - 1];
          WDB[rea + REACTION_EMAX] = XSS[JXS[0] + NES - 2];

          /* Store number of energy points, pointer to grid and */
          /* index to first point */

          WDB[rea + REACTION_PTR_EGRID] = (double)(JXS[0] - 1);
          WDB[rea + REACTION_XS_NE] = (double)NES;
          WDB[rea + REACTION_XS_I0] = 0.0;

          /* Store pointer to XS data */

          WDB[rea + REACTION_PTR_XS] = (double)(JXS[0] - 1 + 4*NES);

          /* Set ures energy boundaries (mites nää?) */

          WDB[rea + REACTION_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];
          WDB[rea + REACTION_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMIN];

          /* Multiplication (not used but must be set to avoid error) */

          WDB[rea + REACTION_WGT_F] = 1.0;
        }

      /* Include photon production from total block */

      if (((long)RDB[DATA_INCLUDE_PHOT_PROD_XS] == YES) && (JXS[11] > 0))
          {
            /* Add photon production */

            rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);
            WDB[nuc + NUCLIDE_PTR_PHOTPRODXS] = (double)rea;

            /* Put nuclide pointer */

            WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

            /* Put type and MT */

            WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;
            WDB[rea + REACTION_MT] = 202.0;

            /* Put awr */

            WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

            /* Set minimum and maximum energy */

            WDB[rea + REACTION_EMIN] = XSS[JXS[0] - 1];
            WDB[rea + REACTION_EMAX] = XSS[JXS[0] + NES - 2];

            /* Store number of energy points, pointer to grid and */
            /* index to first point */

            WDB[rea + REACTION_PTR_EGRID] = (double)(JXS[0] - 1);
            WDB[rea + REACTION_XS_NE] = (double)NES;
            WDB[rea + REACTION_XS_I0] = 0.0;

            /* Store pointer to XS data */

            WDB[rea + REACTION_PTR_XS] = (double)(JXS[11] - 1);

            /* Set ures energy boundaries (mites nää?) */

            WDB[rea + REACTION_URES_EMIN] = INFTY;
            WDB[rea + REACTION_URES_EMAX] = -INFTY;

            /* Multiplication (not used but must be set to avoid error) */

            WDB[rea + REACTION_WGT_F] = 1.0;
          }

      /* Loop over reaction channels in SIG block*/

      for (nr = 0; nr < NTR; nr++)
        {
          /* Get pointer to SIG-block (Table F-10, page F-17) */

          L = (long)XSS[JXS[5] - 1 + nr] + JXS[6] - 1;

          /* Get number of energy points */

          NES = (long)XSS[L];

          /* Pointer to energy array */

          L0 = JXS[0] - 1 + NXS[2] - NES;

          /* Tässä oli 22.8.2011 asti tarkistus > 2, mutta se jättää */
          /* mt 37:n 94244 / endfb68 pois */

          if (NES > 0)
            {
              /* Allocate memory */

              rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

              /* Put nuclide pointer */

              WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

              /* Reaction index (used for energy distributions) */

              WDB[rea + REACTION_NR] = (double)nr;

              /* Get mt */

              mt = (long)XSS[JXS[2] + nr - 1];
              WDB[rea + REACTION_MT] = (double)mt;

              /* Check fissile */

              if (((mt > 17) && (mt < 22)) || (mt == 38))
                {
                  /* Set flag */

                  SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_FISSILE);

                  /* Set default energy boundaries */

                  WDB[rea + REACTION_FISSY_IE0] = -INFTY;
                  WDB[rea + REACTION_FISSY_IE1] = 1E+6;
                  WDB[rea + REACTION_FISSY_IE2] = 1E+9;
                }

              /* Store minimum and maximum energy */

              WDB[rea + REACTION_EMIN] = XSS[L0];
              WDB[rea + REACTION_EMAX] = XSS[L0 + NES - 1];

              /* Store number of energy points, pointer to grid and */
              /* index to first point */

              WDB[rea + REACTION_PTR_EGRID] = (double)(JXS[0] - 1);
              WDB[rea + REACTION_XS_NE] = (double)NES;
              WDB[rea + REACTION_XS_I0] = (double)(NXS[2] - NES);

              /* Store pointer to XS data */

              WDB[rea + REACTION_PTR_XS] = (double)(L + 1);

              /* Set ures energy boundaries */

              if ((mt == 102) || (mt == 18) || (mt == 19))
                {
                  /* Set ures energy boundaries */

                  WDB[rea + REACTION_URES_EMIN] =
                    RDB[nuc + NUCLIDE_URES_EMIN];
                  WDB[rea + REACTION_URES_EMAX] =
                    RDB[nuc + NUCLIDE_URES_EMIN];
                }
              else
                {
                  /* Reset boundaries */

                  WDB[rea + REACTION_URES_EMIN] = INFTY;
                  WDB[rea + REACTION_URES_EMAX] = -INFTY;
                }

              /* Put awr */

              WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

              /* Store Q-value */

              WDB[rea + REACTION_Q] = XSS[JXS[3] - 1 + nr];

              /* Multiplication and frame of reference */

              WDB[rea + REACTION_TY] = XSS[JXS[4] - 1 + nr];

              /* Check known error in ACE files */

              if ((XSS[JXS[4] - 1 + nr] == 0.0) &&
                  (((mt > 17) && (mt < 22)) || (mt == 38)))
                {
                  /* Print warning */

#ifdef DEBUG

                  Warn(FUNCTION_NAME,
                       "Conflicting reaction type for fission (%s mt %ld)",
                       GetText(nuc + NUCLIDE_PTR_NAME), mt);

#endif
                  /* Set ty for fission */

                  WDB[rea + REACTION_TY] = 19.0;
                }

              /* Put weight multiplicator */

              if((RDB[rea + REACTION_TY] == 0.0) ||
                 (RDB[rea + REACTION_TY] == 19.0) ||
                 (fabs(RDB[rea + REACTION_TY]) > 100.0))
                WDB[rea + REACTION_WGT_F] = 1.0;
              else if (fabs(RDB[rea + REACTION_TY]) < 7.0)
                WDB[rea + REACTION_WGT_F] = fabs(RDB[rea + REACTION_TY]);
              else
                Die(FUNCTION_NAME, "Invalid TYR value: %ld\n",
                    (long)RDB[rea + REACTION_TY]);

              /* Override fission if switched off (18.7.2013 / 2.1.15) */

              if (((long)RDB[DATA_NPHYS_SAMPLE_FISS] == NO) &&
                  ((long)RDB[rea + REACTION_TY] == 19.0))
                WDB[rea + REACTION_TY] = 0.0;

              /* Set interpolation mode */

              WDB[rea + REACTION_ITP] = 0.0;

              /* Set branching fraction to 1.0 */

              WDB[rea + REACTION_BR] = 1.0;

              /* Pointer to angular distribution (NOTE: L is re-used) */

              if ((L = (long)XSS[JXS[7] + nr]) > 0)
                WDB[rea + REACTION_PTR_ANG] = (double)(L - 1 + JXS[8]);
              else
                WDB[rea + REACTION_PTR_ANG] = NULLPTR;

              /* Check type */

              if (((mt > 10) && (mt < 100)) || ((mt > 101) && (mt < 200)) ||
                  ((mt > 599) && (mt < 851)) || ((mt > 874) && (mt < 892)))
                {
                  /* Partial reaction */

                  WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
                }
              else if (mt == 4)
                {
                  /* Total inelastic, needed for forming some metastable  */
                  /* states in burnup calculation. Type set to special in */
                  /* addnuclide.c (JLe / 9.8.2016 / 2.1.26) */

                  WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
                }
              else if (mt == 5)
                {
                  /* Combination of multiple inelastic channels */
                  /* (this used to be a problem) */

                  WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
                }
              else
                {
                  /* Special */

                  WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;
                }
            }
        }

      /* Get number of photon production reactions */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DBRC)
        NTR = 0;
      else
        NTR = NXS[5];

      /* Loop over reaction channels in SIGP block*/

      for (nr = 0; nr < NTR; nr++)
        {
          /* Get pointer to SIGP-block (Table F-10, page F-17) */

          L = (long)XSS[JXS[13] - 1 + nr] + JXS[14] - 1;

          /* Tää on ihan vaiheessa */

          I0 = (long)XSS[L - 1];
          NES = (long)XSS[L];
        }

      /***********************************************************************/

      /***** URES data *******************************************************/

      /* Check if ures probability table data is available (ures data */
      /* is now read for DBRC nuclides as well). */

      if ((L = JXS[22] - 1) > 0)
        {
          /* Add to ures counter */

          WDB[DATA_URES_AVAIL] = RDB[DATA_URES_AVAIL] + 1.0;

          /* Set available flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_AVAIL);

          /* Check ures option */

          if ((long)RDB[DATA_USE_URES] == NO)
            {
              /* Pointer to list */

              if ((ptr = (long)RDB[DATA_URES_PTR_USE_LIST]) < VALID_PTR)
                L = -1;
              else
                {
                  /* Loop over list and compare */

                  while ((long)RDB[ptr] > 0)
                    {
                      if (!strcmp(GetText(ptr),
                                  GetText(nuc + NUCLIDE_PTR_NAME)))
                        {
                          L = -1;
                          break;
                        }

                      ptr++;
                    }
                }
            }
          else if ((ptr = (long)RDB[DATA_URES_PTR_USE_LIST]) > VALID_PTR)
            {
              /* Loop over list and compare */

              while ((long)RDB[ptr] > 0)
                {
                  if (!strcmp(GetText(ptr), GetText(nuc + NUCLIDE_PTR_NAME)))
                    break;

                  ptr++;
                }

              if ((long)RDB[ptr] < 1)
                L = -1;
            }
        }

      /* Check if data is available and used */

      if (L > 0)
        {
          /* Set flag */

          SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_USED);

          /* Add counter */

          WDB[DATA_URES_USED] = RDB[DATA_URES_USED] + 1.0;
        }

      /**********************************************************************/
    }
  else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
    {
      /**********************************************************************/

      /***** Dosimetry data *************************************************/

      /* Get number of reactions */

      NTR = NXS[3];

      /* Reset nuclide-wise minimum and maximum energy */

      WDB[nuc + NUCLIDE_EMIN] = INFTY;
      WDB[nuc + NUCLIDE_EMAX] = -INFTY;

      /* Loop over reaction channels */

      for (nr = 0; nr < NTR; nr++)
        {
          /* Get pointer to SIGD-block (Table F-22, page F-35) */

          L0 = (long)XSS[JXS[5] - 1 + nr] + JXS[6] - 1;

          /* Get reaction MT (Table F-6, page F-15) */

          mt = (long)XSS[JXS[2] - 1 + nr];

          /* Check number of interpolation regions */

          if ((long)XSS[L0 - 1] > 0)
            Die(FUNCTION_NAME, "Non-linear interpolation (%s mt %ld)",
                GetText(nuc + NUCLIDE_PTR_NAME), mt);

          /* Get number original energy points */

          NES = (long)XSS[L0];
          CheckValue(FUNCTION_NAME, "NES", " (dosimetry)", NES, 0, INFTY);

          /* Check and modify coincident energy points */

          nfix = 0;

          /* Loop over grid */

          for (n = 0; n < NES - 1; n++)
            if (XSS[L0 + n + 1] == XSS[L0 + n + 2])
              {
                /* Is there another, non-coincident point at higher energy? */

                if ((n < NES - 2) && (XSS[L0 + n + 1] != XSS[L0 + n + 3]))
                  XSS[L0 + n + 2] = XSS[L0 + n + 2] +
                    0.000000001*(XSS[L0 + n + 3] - XSS[L0 + n + 2]);

                /* How about at lower energy? */

                else if (n != 0)
                  XSS[L0 + n + 1] = XSS[L0 + n + 1] -
                    0.000000001*(XSS[L0 + n + 1] - XSS[L0 + n]);
                else
                  Die(FUNCTION_NAME, "Coincident energies in two-point grid");

                /* Update count */

                nfix++;
              }

          if (nfix > 0)
            Note(0, "%ld coincident energy points adjusted (%s mt %ld)", nfix,
                 GetText(nuc + NUCLIDE_PTR_NAME), mt);

          /* Allocate memory for data */

          rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

          /* Reaction index (used for energy distributions) */

          WDB[rea + REACTION_NR] = (double)nr;

          /* Set mt */

          WDB[rea + REACTION_MT] = (double)mt;

          /* Put awr */

          WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

          /* Store minimum and maximum energy */

          WDB[rea + REACTION_EMIN] = XSS[L0 + 1];
          WDB[rea + REACTION_EMAX] = XSS[L0 + NES];

          /* Compare to nuclide-wise values */

          if (RDB[rea + REACTION_EMIN] < RDB[nuc + NUCLIDE_EMIN])
            WDB[nuc + NUCLIDE_EMIN] = RDB[rea + REACTION_EMIN];

          if (RDB[rea + REACTION_EMAX] > RDB[nuc + NUCLIDE_EMAX])
            WDB[nuc + NUCLIDE_EMAX] = RDB[rea + REACTION_EMAX];

          /* Store number of energy points, pointer to grid and */
          /* index to first point */

          WDB[rea + REACTION_PTR_EGRID] = (double)(L0 + 1);
          WDB[rea + REACTION_XS_NE] = (double)NES;

          /* Store pointer to XS data */

          WDB[rea + REACTION_PTR_XS] = (double)(L0 + NES + 1);

          /* Set reaction type */

          WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;
        }

      /**********************************************************************/
    }

  else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB)
    {
      /**********************************************************************/

      /***** S(a,b) data ****************************************************/

      /* Check ZA is found in IZ list */

      for (n = 0; n < 16; n++)
        if (IZ[n] == (long)RDB[nuc + NUCLIDE_ZA])
          break;

      /* Check if not found */

      if (n == 16)
        {
          /* Print error message */

          fprintf(outp, "S(a,b) library \"%s\" can be associated with", HZ1);
          fprintf(outp, " the following nuclides:\n\n");

          n = 0;
          while (IZ[n] > 0)
            fprintf(outp, "%ld\n", IZ[n++]);

          Error(0, "Invald ZA %ld in material card",
                (long)RDB[nuc + NUCLIDE_ZA]);
        }

      /* Set S(a,b) flag */

      SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_SAB_DATA);

      /* Number of reaction modes */

      if (JXS[3] - 1 > 0)
        NTR = 2;
      else
        NTR = 1;

      /* Reset nuclide energy boundaries */

      WDB[nuc + NUCLIDE_EMIN] = INFTY;
      WDB[nuc + NUCLIDE_EMAX] = -INFTY;

      /* Avoid compiler warning */

      Emax = -1.0;

      /* Loop over reaction channels */

      for (nr = 0; nr < NTR; nr++)
        {
          /* Pointer to (in)elastic data (Table F-23, page F-36) */

          if (nr == 0)
            L0 = JXS[0] - 1;
          else
            L0 = JXS[3] - 1;

          /* Get number of energy points */

          NES = (long)XSS[L0];

          /* Allocate memory for reaction data */

          rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

          /* Put awr */

          WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

          /* Set mt */

          if (nr == 0)
            WDB[rea + REACTION_MT] = 1004.0;
          else
            WDB[rea + REACTION_MT] = 1002.0;

          /* Set interpolation mode */

          if ((nr == 1) && (NXS[4] == 4))
            WDB[rea + REACTION_ITP] = 4.0;
          else
            WDB[rea + REACTION_ITP] = 0.0;

          /* Store minimum and maximum energy */

          WDB[rea + REACTION_EMIN] = XSS[L0 + 1];
          WDB[rea + REACTION_EMAX] = XSS[L0 + NES];

          /* Reset minimum and maximum emission energy */

          WDB[rea + REACTION_SAB_MIN_EM_E] = INFTY;
          WDB[rea + REACTION_SAB_MAX_EM_E] = -INFTY;

          /* Maximum S(a,b) energy (needed to get the extra point in */
          /* elastic channel */

          if (nr == 0)
            Emax = XSS[L0 + NES];
          else if (Emax < XSS[L0 + NES])
            Emax = XSS[L0 + NES];

          /* Store */

          WDB[rea + REACTION_SAB_EMAX] = Emax;

          /* Compare to nuclide minimum */

          if (RDB[rea + REACTION_EMIN] < RDB[nuc + NUCLIDE_EMIN])
            WDB[nuc + NUCLIDE_EMIN] = RDB[rea + REACTION_EMIN];

          /* Compare to nuclide maximum */

          if (RDB[rea + REACTION_EMAX] > RDB[nuc + NUCLIDE_EMAX])
            WDB[nuc + NUCLIDE_EMAX] = RDB[rea + REACTION_EMAX];

          /* Store number of energy points, pointer to grid and */
          /* index to first point */

          WDB[rea + REACTION_PTR_EGRID] = (double)(L0 + 1);
          WDB[rea + REACTION_XS_NE] = (double)NES;
          WDB[rea + REACTION_XS_I0] = 0.0;

          /* Store pointer to XS data */

          WDB[rea + REACTION_PTR_XS] = (double)(L0 + 1 + NES);

          /* Reset ures energy boundaries */

          WDB[rea + REACTION_URES_EMIN] = INFTY;
          WDB[rea + REACTION_URES_EMAX] = -INFTY;

          /* Store Q-value */

          WDB[rea + REACTION_Q] = 0.0;

          /* Multiplication and frame of reference */

          WDB[rea + REACTION_TY] = 1.0;
          WDB[rea + REACTION_WGT_F] = 1.0;

          /* Set branching fraction to 1.0 */

          WDB[rea + REACTION_BR] = 1.0;

          /* Set type */

          WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
        }

      /**********************************************************************/
    }
  else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
    {
      /**********************************************************************/

      /***** Photon interaction data ****************************************/

      /* Check Z */

      if ((long)RDB[nuc + NUCLIDE_Z] > 99)
        Error(0, "Photon physics models not available for %s (%s)",
              ZAItoIso((long)RDB[nuc + NUCLIDE_Z], 1),
              GetText(nuc + NUCLIDE_PTR_NAME));

      /* Number of energy points */

      NES = NXS[2];

      /* Put pointer to nuclide energy grid and nuber of points */

      WDB[nuc + NUCLIDE_PTR_EGRID] = (double)(JXS[0] - 1);
      WDB[nuc + NUCLIDE_EGRID_NE] = (double)NES;

      /* Compare minimum and maximum value to XS limits */

      if (exp(XSS[JXS[0] - 1]) < RDB[DATA_PHOTON_XS_EMIN])
        WDB[DATA_PHOTON_XS_EMIN] = exp(XSS[JXS[0] - 1]);

      if (exp(XSS[JXS[0] + NES - 2]) > RDB[DATA_PHOTON_XS_EMAX])
        WDB[DATA_PHOTON_XS_EMAX] = exp(XSS[JXS[0] + NES - 2]);

      /* Loop over 4 reaction modes (incoherent, coherent, photoelectric */
      /* and pair production) and average heating numbers. */

      for (nr = 0; nr < 5; nr++)
        {
          /* Add reaction */

          rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

          /* Put nuclide pointer */

          WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

          /* Reaction index (used for energy distributions) */

          WDB[rea + REACTION_NR] = -1.0;

          /* Put type */

          if (nr < 4)
            WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
          else
            WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

          /* Put mt */

          if (nr == 0)
            WDB[rea + REACTION_MT] = 504.0;
          else if (nr == 1)
            WDB[rea + REACTION_MT] = 502.0;
          else if (nr == 2)
            WDB[rea + REACTION_MT] = 522.0;
          else if (nr == 3)
            WDB[rea + REACTION_MT] = 516.0;
          else if (nr == 4)
            WDB[rea + REACTION_MT] = 301.0;

          /* Set minimum and maximum energy */

          WDB[rea + REACTION_EMIN] = exp(XSS[JXS[0] - 1]);
          WDB[rea + REACTION_EMAX] = exp(XSS[JXS[0] + NES - 2]);

          WDB[nuc + NUCLIDE_EMIN] = exp(XSS[JXS[0] - 1]);
          WDB[nuc + NUCLIDE_EMAX] = exp(XSS[JXS[0] + NES - 2]);

          /* Store number of energy points, pointer to grid and */
          /* index to first point */

          WDB[rea + REACTION_PTR_EGRID] = (double)(JXS[0] - 1);
          WDB[rea + REACTION_XS_NE] = (double)NES;
          WDB[rea + REACTION_XS_I0] = 0.0;

          /* Store pointer to XS data */

          if (nr < 4)
            WDB[rea + REACTION_PTR_XS] = (double)(JXS[0] - 1 + (nr + 1)*NES);
          else
            WDB[rea + REACTION_PTR_XS] = (double)(JXS[4] - 1);

          /* Multiplication */

          WDB[rea + REACTION_WGT_F] = 1.0;
        }

      /**********************************************************************/
    }
  else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
    {
      /**********************************************************************/

      /***** Transmutation data *********************************************/

      /* Check ace type */

      if ((long)ACE[ace + ACE_TYPE] != NUCLIDE_TYPE_TRANSMUXS)
        Die(FUNCTION_NAME, "Invalid ace type");

      /* Number of energy points */

      NES = NXS[2];

      /* Put pointer to nuclide energy grid and nuber of points */

      WDB[nuc + NUCLIDE_PTR_EGRID] = (double)(JXS[0] - 1);
      WDB[nuc + NUCLIDE_EGRID_NE] = (double)NES;

      /* Set minimum and maximum energy */

      WDB[nuc + NUCLIDE_EMIN] = XSS[JXS[0] - 1];
      WDB[nuc + NUCLIDE_EMAX] = XSS[JXS[0] + NES - 2];

      /* Get number of reactions (minus elastic scattering). Include */
      /* only elastic scattering for DBRC nuclides. */

      NTR = NXS[3];

      /* Loop over reaction channels */

      for (nr = 0; nr < NTR; nr++)
        {
          /* Get pointer to SIG-block (Table F-10, page F-17) */

          L = (long)XSS[JXS[5] - 1 + nr] + JXS[6] - 1;

          /* Get number of energy points */

          NES = (long)XSS[L];

          /* Pointer to energy array (tää vaikuttaa epäilyttävältä) */

          L0 = JXS[0] - 1 + NXS[2] - NES;

          /* Get mt */

          mt = (long)XSS[JXS[2] + nr - 1];

          /* Check number of energy points and mt (ei fissiota nyt) */

          if ((NES > 0) && ((mt > 101) && (mt < 200)))
            {
              /* Allocate memory */

              rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

              /* Put nuclide pointer */

              WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

              /* Put mt */

              WDB[rea + REACTION_MT] = (double)mt;

              /* Store minimum and maximum energy */

              WDB[rea + REACTION_EMIN] = XSS[L0];
              WDB[rea + REACTION_EMAX] = XSS[L0 + NES - 1];

              /* Store number of energy points, pointer to grid and */
              /* index to first point */

              WDB[rea + REACTION_PTR_EGRID] = (double)(JXS[0] - 1);
              WDB[rea + REACTION_XS_NE] = (double)NES;
              WDB[rea + REACTION_XS_I0] = (double)(NXS[2] - NES);

              /* Store pointer to XS data */

              WDB[rea + REACTION_PTR_XS] = (double)(L + 1);

              /* Put awr */

              WDB[rea + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

              /* Store Q-value */

              WDB[rea + REACTION_Q] = XSS[JXS[3] - 1 + nr];

              /* Multiplication and frame of reference */

              WDB[rea + REACTION_TY] = XSS[JXS[4] - 1 + nr];

              /* Check known error in ACE files (jos fissio joskus lisätään) */

              if ((XSS[JXS[4] - 1 + nr] == 0.0) &&
                  (((mt > 17) && (mt < 22)) || (mt == 38)))
                {
                  /* Print warning */

#ifdef DEBUG

                  Warn(FUNCTION_NAME,
                       "Conflicting reaction type for fission (%s mt %ld)",
                       GetText(nuc + NUCLIDE_PTR_NAME), mt);

#endif
                  /* Set ty for fission */

                  WDB[rea + REACTION_TY] = 19.0;
                }

              /* Set branching fraction to 1.0 */

              WDB[rea + REACTION_BR] = 1.0;

              /* Set type */

              WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;
            }
        }

      /**********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid nuclide type (%s)", name);

  /**************************************************************************/

  /***** Add total inelastic ************************************************/

  /* Reset minimum and maximum energy */

  Emin = INFTY;
  Emax = -INFTY;

  /* Count number of inelastic modes */

  n = 0;

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Check mt */

      if ((mt > 50) && (mt < 92))
        {
          /* Update count */

          n++;

          /* Compare energies */

          if (RDB[rea + REACTION_EMIN] < Emin)
            Emin = RDB[rea + REACTION_EMIN];
          if (RDB[rea + REACTION_EMAX] > Emax)
            Emax = RDB[rea + REACTION_EMAX];
        }

      else if ((long)RDB[rea + REACTION_MT] == 4)
        {
          /* Total inelastic already exists, reset count */

          n = 0;

          /* Break loop */

          break;
        }

      /* Next */

      rea = NextItem(rea);
    }

  /* Check count */

  if (n > 0)
    {
      /* Allocate memory */

      rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

      /* Put nuclide pointer */

      WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

      /* Put mt */

      WDB[rea + REACTION_MT] = 4.0;

      /* Put type */

      WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_PARTIAL;

      /* Put minimum and maximum energy */

      WDB[rea + REACTION_EMIN] = Emin;
      WDB[rea + REACTION_EMAX] = Emax;

      /* Set branching fraction to 1.0 */

      WDB[rea + REACTION_BR] = 1.0;

      /* Make sure that cross sections are not processed. */

      WDB[rea + REACTION_PTR_EGRID] = -INFTY;
    }

  /**************************************************************************/

  /***** Remove redundant reaction modes ************************************/

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Redundant (n,p) */

      if ((long)RDB[rea + REACTION_MT] == 103)
        {
          /* Loop over reactions */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (ptr > VALID_PTR)
            {
              /* Check mt and set type to special */

              if (((long)RDB[ptr + REACTION_MT] > 599) &&
                  ((long)RDB[ptr + REACTION_MT] < 650))
                WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

              /* Next reaction */

              ptr = NextItem(ptr);
            }
        }

      /* Redundant (n,d) */

      if ((long)RDB[rea + REACTION_MT] == 104)
        {
          /* Loop over reactions */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (ptr > VALID_PTR)
            {
              /* Check mt and set type to special */

              if (((long)RDB[ptr + REACTION_MT] > 649) &&
                  ((long)RDB[ptr + REACTION_MT] < 700))
                WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

              /* Next reaction */

              ptr = NextItem(ptr);
            }
        }

      /* Redundant (n,t) */

      if ((long)RDB[rea + REACTION_MT] == 105)
        {
          /* Loop over reactions */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (ptr > VALID_PTR)
            {
              /* Check mt and set type to special */

              if (((long)RDB[ptr + REACTION_MT] > 699) &&
                  ((long)RDB[ptr + REACTION_MT] < 750))
                WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

              /* Next reaction */

              ptr = NextItem(ptr);
            }
        }

      /* Redundant (n,He-3) */

      if ((long)RDB[rea + REACTION_MT] == 106)
        {
          /* Loop over reactions */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (ptr > VALID_PTR)
            {
              /* Check mt and set type to special */

              if (((long)RDB[ptr + REACTION_MT] > 749) &&
                  ((long)RDB[ptr + REACTION_MT] < 800))
                WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

              /* Next reaction */

              ptr = NextItem(ptr);
            }
        }

      /* Redundant (n,a) */

      if ((long)RDB[rea + REACTION_MT] == 107)
        {
          /* Loop over reactions */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (ptr > VALID_PTR)
            {
              /* Check mt and set type to special */

              if (((long)RDB[ptr + REACTION_MT] > 799) &&
                  ((long)RDB[ptr + REACTION_MT] < 850))
                WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_SPECIAL;

              /* Next reaction */

              ptr = NextItem(ptr);
            }
        }

      /* Next reaction */

      rea = NextItem(rea);
    }

  /**************************************************************************/

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
    {
      /* Store position */

      if ((pos = ftell(fp)) == -1)
        Die(FUNCTION_NAME, "ftell error\n");

      /* Set fissile flag */

      fiss = (long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE;

      /* Set ures flag */

      ures = (long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED;

      /* Get edep mode */

      edep = (long)RDB[DATA_EDEP_MODE];

      /* Fission energy deposition data MF1 MT458 */

      if (fiss && (edep > EDEP_MODE_CONSTANT))
        ReadPendfData(nuc, fp, pos, 1458);

      /* Non-local KERMA MF3 MT301 and URES data MF2 MT153 */

      if (edep == EDEP_MODE_NEUTRON_PHOTON)
        {
          ReadPendfData(nuc, fp, pos, 3301);

          if (ures)
            ReadPendfData(nuc, fp, pos, 2153);
        }

      /* Fission KERMA MF3 MT318 */

      if (fiss && (edep > EDEP_MODE_MT458))
        {
          if (edep == EDEP_MODE_LOCAL_PHOTON)
            ReadPendfData(nuc, fp, pos, 3319);
          else
            ReadPendfData(nuc, fp, pos, 3318);
        }
    }
  /* Close file  */

  fclose(fp);
}

/*****************************************************************************/
