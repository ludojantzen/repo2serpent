/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processrelaxation.c                            */
/*                                                                           */
/* Created:       2014/11/23 (TKa)                                           */
/* Last modified: 2017/04/24 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Processes atomic relaxation data                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessRelaxation:"

/*****************************************************************************/

void ProcessRelaxation() {

  long i, j, elementfound, Nss, subi, ntr, eln, maxdesignator, ptrd2i, ptrntr,
      ptrebi, ptreln, ptretr, ptrletr, ptridxrad, ptrsubi, ptrsubj, ptrsubk,
      ptrsubj2, ptrsubk2, ptretr2, ptrletr2, ptridxrad2, ptrelncdf, nuc, loc0,
      ptralias11, ptralias21, ptrcutoff1, ptralias12, ptralias22, ptrcutoff2,
      nradtr, nnonradtr, nradtrEmin;
  double ebi, Etrtmp, ftrtmp;
  double *ftrarr, *Erad;
  const int nbuf = 500;    /* NOTE: same as linebuf length */
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100];
  FILE *fp;

  /* Number of transitions */
  nradtr = nnonradtr = nradtrEmin = 0;

  /* Array for all the radiative transitions above DATA_PHOTON_EMIN */
  Erad = (double *)Mem(MEM_ALLOC, PHOTON_NRADTR_MAX, sizeof(double));

  fprintf(outp, "Processing atomic relaxation data...\n");

  /* Set file name */
  sprintf(fname, "%s%s", GetText(DATA_PHOTON_DATA_DIR),
          GetText(DATA_PHOTON_RELAX_FNAME));

  /* Open file */
  if (!(fp = fopen(fname, "r")))
    Die(FUNCTION_NAME, "Unable to open file for reading");

  /* Loop over nuclides */
  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR) {

    if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) {

      /* Rewind the file because the data might not be ordered by Z */
      rewind(fp);

      /* Set Z string */
      sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);

      elementfound = 0;

      while (fgets(linebuf, nbuf, fp)) {

        if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
           (strcmp(tmpstr, strZ) == 0)) {

          elementfound = 1;

          /* Allocate memory for atomic relaxation block */
          loc0 = NewItem(nuc + NUCLIDE_PTR_RELAX, RELAX_BLOCK_SIZE);

          if (!fgets(linebuf, nbuf, fp))
            Die(FUNCTION_NAME, "Can't read relaxation data");

          /* Read and store number of subshells */
          if (sscanf(linebuf, "NSS %ld\n", &Nss) == 0)
            Die(FUNCTION_NAME, "Can't read relaxation data NSS");
          CheckValue(FUNCTION_NAME, "Nss", "", (double)Nss, 0.0, 30.0);
          WDB[loc0 + RELAX_NSS] = (double)Nss;

          /* Create a map from subshell designator to subshell index */
          /* (Designator from ENDF, index 0,1,2 etc) */
          maxdesignator = 599 - 534 + 1;
          ptrd2i = ReallocMem(DATA_ARRAY, maxdesignator + 1);
          WDB[loc0 + RELAX_D2IMAP] = (double)ptrd2i;

          /* Initialize map */
          for (i = 0; i < maxdesignator + 1; i++)
            WDB[ptrd2i + i] = -1.0;

          /* Read and init subshell data */
          if (Nss > 0) {

            /* Subshell data array size */
            ptrntr = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_NTR] = (double)ptrntr;

            /* Binding energy array */
            ptrebi = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_EBI] = (double)ptrebi;

            /* Number of electrons in subshell */
            ptreln = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_ELN] = (double)ptreln;

            /* Transition energy */
            ptretr = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_ETR] = (double)ptretr;

            /* Log-transition energy */
            ptrletr = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_LETR] = (double)ptrletr;

            /* Indexes for radiative transitions */
            ptridxrad = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_IDX_RAD] = (double)ptridxrad;

            /* Subshell designator array */
            ptrsubi = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_SUBI] = (double)ptrsubi;

            /* Secondary subshell designator array */
            ptrsubj = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_SUBJ] = (double)ptrsubj;

            /* Tertiary subshell designator array */
            ptrsubk = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_SUBK] = (double)ptrsubk;

            /* Electron number CDF */
            ptrelncdf = ReallocMem(DATA_ARRAY, Nss+1);
            WDB[loc0 + RELAX_ELNCDF] = (double)ptrelncdf;
            WDB[ptrelncdf] = 0.0;

            /* Alias and cutoff arrays for Walker's alias method */
            ptralias11 = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_ALIAS1] = (double)ptralias11;

            ptralias21 = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_ALIAS2] = (double)ptralias21;

            ptrcutoff1 = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_CUTOFF] = (double)ptrcutoff1;


            /* Loop over shells */
            for (i = 0; i < Nss; i++) {

              /* Read and store SUBI*/
              if (fgets(linebuf, nbuf, fp)) {
                if (sscanf(linebuf, "SUBI %ld\n", &subi) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data SUBI");

                /* Check subi*/
                CheckValue(FUNCTION_NAME, "subi", "", (double)subi, 1.0, 
                           (double)maxdesignator);

                WDB[ptrsubi++] = (double)subi;

                /* Map SUBI to index */
                WDB[ptrd2i + subi] = (double)i;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Read number of transitions */
              if (fgets(linebuf, nbuf, fp)) {
                if (sscanf(linebuf, "NTR %ld\n", &ntr) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data NTR");

                CheckValue(FUNCTION_NAME, "ntr", "", (double)ntr, 0, 1.0e4);
                WDB[ptrntr++] = (double)ntr;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Read electron binding energy */
              if (fgets(linebuf, nbuf, fp)) {
                if (sscanf(linebuf, "EBI %lf\n", &ebi) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data EBI");
                ebi /= 1.0e6; /* ev to MeV */
                CheckValue(FUNCTION_NAME, "ebi", "", ebi, 0, 0.5);
                WDB[ptrebi++] = ebi;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Read number of electrons in subshell */
              if (fgets(linebuf, nbuf, fp)) {
                if (sscanf(linebuf, "ELN %ld\n", &eln) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data ELN");
                CheckValue(FUNCTION_NAME, "eln", "", (double)eln, 2.0, 6.0);
                WDB[ptreln++] = (double)eln;
                WDB[ptrelncdf+i+1] = RDB[ptrelncdf+i] + (double)eln;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Secondary subshell designator array */
              ptrsubj2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptrsubj++] = (double)ptrsubj2;

              /* Tertiary subshell designator array */
              ptrsubk2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptrsubk++] = (double)ptrsubk2;

              /* Transition energy */
              ptretr2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptretr++] = (double)ptretr2;

              /* Log-transition energy */
              ptrletr2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptrletr++] = (double)ptrletr2;

              /* Radiative transition index */
              ptridxrad2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptridxrad++] = (double)ptridxrad2;

              /* Transition probability array */
              ftrarr = (double *)Mem(MEM_ALLOC, ntr, sizeof(double));

              /* Read transitions */
              for (j = 0; j < ntr; j++) {

                if (fgets(linebuf, nbuf, fp) != NULL) {
                  if (sscanf(linebuf, "%lf %lf %lf %lf\n", &WDB[ptrsubj2++],
                      &WDB[ptrsubk2++], &Etrtmp, &ftrtmp) != 4)
                    Die(FUNCTION_NAME, "Can't read relaxation data, line:\n %s", linebuf);

                  /* ev to MeV */
                  Etrtmp /= 1.0e6;

                  /* Check and set transition energy */
                  CheckValue(FUNCTION_NAME, "Etr", "", Etrtmp, 0.0, 0.5);
                  WDB[ptretr2++] = Etrtmp;
                  WDB[ptrletr2++] = log(Etrtmp);

                  /* Check and set transition probability */
                  CheckValue(FUNCTION_NAME, "ftr", "", ftrtmp, 0.0, 1.0);
                  ftrarr[j] = ftrtmp;

                  /* Count transitions */
                  if (RDB[ptrsubk2-1] == 0) {

                    /* Set radiative transition index */
                    if (Etrtmp > RDB[DATA_PHOTON_EMIN]) {
                      WDB[ptridxrad2++] = (double)nradtrEmin;
                      Erad[nradtrEmin] = Etrtmp;
                      nradtrEmin++;

                      if (nradtrEmin > PHOTON_NRADTR_MAX)
                        Die(FUNCTION_NAME, "Maximum number of radiative "
                            "transitions exceeded (%ld)", PHOTON_NRADTR_MAX);
                    }
                    else {
                      /* -1 for detecting bugs */
                      WDB[ptridxrad2++] = -1;
                    }

                    /* Counter for all radiative transitions */
                    nradtr++;
                  }
                  else {
                    nnonradtr++;
                  }
                }
                else
                  Die(FUNCTION_NAME, "Can't read relaxation data");
              }

              /* Walker's alias arrays */
              ptralias12 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptralias11++] = (double)ptralias12;

              ptralias22 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptralias21++] = (double)ptralias22;

              ptrcutoff2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptrcutoff1++] = (double)ptrcutoff2;

              /* Initialize Walker's alias method */
              WalkerAliasInit(ftrarr, ntr, &WDB[ptralias12], &WDB[ptralias22],
                              &WDB[ptrcutoff2]);

              /* Free transtion probabilities */
              Mem(MEM_FREE, ftrarr);

            }
          }

          /* Exit loop */
          break;
        }
      }

      if (!elementfound)
        Die(FUNCTION_NAME, "Atomic relaxation data not found for element Z=%s",
            strZ);
    }

    nuc = NextItem(nuc);
  }

  /* Close the file */
  fclose(fp);

  /* Calculate CDFs of shell probabilities of photoelectric effect at
   * fluorescence energies */
  if (nradtrEmin > 0)
    ProcessPhotoelectricFluorescenceCDF(Erad, nradtrEmin);


  /* Free memory */
  Mem(MEM_FREE, Erad);

  fprintf(outp, " - %ld radiative transitions (fluorescence)\n", nradtr);
  fprintf(outp, " - %ld non-radiative transitions (Auger electrons)\n",
          nnonradtr);

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/

}

/*****************************************************************************/
