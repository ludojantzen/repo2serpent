/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : xsplotter.c                                    */
/*                                                                           */
/* Created:       2010/12/15 (JLe)                                           */
/* Last modified: 2018/05/25 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Writes cross section data into matlab-format file            */
/*                                                                           */
/* Comments: - pitäiskö tää linkittää user energy gridiin?                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "XSPlotter:"

/*****************************************************************************/

void XSPlotter()
{
  double Emin, Emax, *E;
  long nuc, mat, rea, ne, n, i, ptr, rea1, ncol, ang, erg, l0, Ane;
  long ene, m, np;
  char tmpstr[MAX_STR];
  FILE *fp;
  unsigned long seed;

  /* Check if plot is defined */

  if ((long)RDB[DATA_XSPLOT_NE] < 1)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Init random number sequence */

  seed = ReInitRNG(0);
  SEED[0] = seed;

  /* Set filename */

  sprintf(tmpstr, "%s_xs%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          (long)RDB[DATA_BURN_STEP]);

  fprintf(outp, "Creating XS plot file \"%s\"...\n", tmpstr);

  /* Open file for writing */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Remember collision count */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, 0);

  /***************************************************************************/

  /***** Generate energy grid ************************************************/

  /* Boundaries and grid size */

  if (RDB[DATA_XSPLOT_EMIN] > 0.0)
    Emin = RDB[DATA_XSPLOT_EMIN];
  else
    Emin = RDB[DATA_NEUTRON_EMIN];

  if (RDB[DATA_XSPLOT_EMAX] > 0.0)
    Emax = RDB[DATA_XSPLOT_EMAX];
  else
    Emax = RDB[DATA_NEUTRON_EMAX];

  ne = (long)RDB[DATA_XSPLOT_NE];

  /* Generate array */

  E = MakeArray(Emin, Emax, ne, 2);

  /* Print array */

  fprintf(fp, "E = [\n");

  for (n = 0; n < ne; n++)
    fprintf(fp, "%1.5E\n", E[n]);

  fprintf(fp, "];\n\n");

  /***************************************************************************/

  /***** Print microscopic cross sections ************************************/

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while(nuc > VALID_PTR)
    {
      /* Check type */

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
        {
          /* Set name */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
            tmpstr[strlen(tmpstr) - 4] = '_';

          /* Print reaction mt's */

          fprintf(fp, "%s_mt = [\n", tmpstr);

          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Print mt */

              fprintf(fp, "%4ld %% %s\n", (long)RDB[rea + REACTION_MT],
                      ReactionMT((long)RDB[rea + REACTION_MT], YES));

              /* Next */

              rea = NextItem(rea);
            }

          fprintf(fp, "];\n\n");

          /* Print cross sections */

          fprintf(fp, "%s_xs = [\n", tmpstr);

          for (n = 0; n < ne; n++)
            {
              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Get cross section */

                  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
                    fprintf(fp, "%1.5E ", PhotonMicroXS(rea,  E[n], 0));
                  else
                    fprintf(fp, "%1.5E ", MicroXS(rea, E[n], 0));

                  /* 0K data */

                  if ((rea1 = (long)RDB[rea + REACTION_PTR_0K_DATA])
                      > VALID_PTR)
                    fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));

                  /* Majorant */

                  if ((rea1 = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT])
                      > VALID_PTR)
                    fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));

                  /* Next reaction */

                  rea = NextItem(rea);
                }

              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");
        }
      else if (((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY) ||
               ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                NUCLIDE_FLAG_TRANSMU_DATA))
        {
          /* Set name */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
            tmpstr[strlen(tmpstr) - 4] = '_';

          /* Print reaction mt's */

          fprintf(fp, "%s_mt = [\n", tmpstr);

          if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
            fprintf(fp, "%4ld %% Total\n", (long)501);
          else
            {
              fprintf(fp, "%4ld %% Total\n", (long)1);
              if ((long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS] > VALID_PTR)
                fprintf(fp, "%4ld %% Sum of absorption\n", (long)101);
            }

          if ((long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS] > VALID_PTR)
            fprintf(fp, "%4ld %% Heat production\n", (long)301);

          if ((long)RDB[nuc + NUCLIDE_PTR_PHOTPRODXS] > VALID_PTR)
            fprintf(fp, "%4ld %% Photon production\n", (long)202);

          /* Pointer total xs */

          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          if ((long)RDB[rea + REACTION_PTR_0K_DATA] > VALID_PTR)
            fprintf(fp, "%4ld %% %s (0K)\n",
                    (long)RDB[rea + REACTION_MT],
                    ReactionMT((long)RDB[rea + REACTION_MT], YES));

          if ((long)RDB[rea + REACTION_PTR_TMP_MAJORANT] > VALID_PTR)
            fprintf(fp, "%4ld %% %s (temperature Majorant)\n",
                    (long)RDB[rea + REACTION_MT],
                    ReactionMT((long)RDB[rea + REACTION_MT], YES));

          /* Pointer to partial list */

          ptr = (long)RDB[rea + REACTION_PTR_PARTIAL_LIST];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over partials */

          i = 0;
          while ((rea = ListPtr(ptr, i++)) > VALID_PTR)
            {
              /* Pointer to reaction data */

              rea = (long)RDB[rea + RLS_DATA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Print mt */

              fprintf(fp, "%4ld %% %s\n", (long)RDB[rea + REACTION_MT],
                      ReactionMT((long)RDB[rea + REACTION_MT], YES));

              if ((long)RDB[rea + REACTION_PTR_0K_DATA] > VALID_PTR)
                fprintf(fp, "%4ld %% %s (0K)\n",
                        (long)RDB[rea + REACTION_MT],
                        ReactionMT((long)RDB[rea + REACTION_MT], YES));

              if ((long)RDB[rea + REACTION_PTR_TMP_MAJORANT] > VALID_PTR)
                fprintf(fp, "%4ld %% %s (temperature Majorant)\n",
                        (long)RDB[rea + REACTION_MT],
                        ReactionMT((long)RDB[rea + REACTION_MT], YES));
            }

          fprintf(fp, "];\n\n");

          /* Print cross sections */

          fprintf(fp, "%s_xs = [\n", tmpstr);

          for (n = 0; n < ne; n++)
            {
              /* Pointer to total xs */

              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Get cross section */

              if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
                fprintf(fp, "%1.5E ", PhotonMicroXS(rea, E[n], 0));
              else
                fprintf(fp, "%1.5E ", MicroXS(rea, E[n], 0));

              /* Pointer to total absorption xs */

              if ((rea1 = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS]) > VALID_PTR)
                fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));

              /* Heat production */

              if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
                {
                  ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_HEATPRODXS];
                  fprintf(fp, "%1.5E ", PhotonMicroXS(ptr, E[n], 0));
                }
              else if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_HEATPRODXS]) >
                       VALID_PTR)
                fprintf(fp, "%1.5E ", MicroXS(ptr, E[n], 0));

              /* Photon production */

              if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
                if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTPRODXS]) > VALID_PTR)
                  fprintf(fp, "%1.5E ", MicroXS(ptr, E[n], 0));

              /* 0K data */

              if ((rea1 = (long)RDB[rea + REACTION_PTR_0K_DATA])
                  > VALID_PTR)
                fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));

              /* Majorant */

              if ((rea1 = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT])
                  > VALID_PTR)
                fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));

              /* Pointer to partial list */

              ptr = (long)RDB[rea + REACTION_PTR_PARTIAL_LIST];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Loop over partials */

              i = 0;
              while ((rea = ListPtr(ptr, i++)) > VALID_PTR)
                {
                  /* Pointer to reaction data */

                  rea = (long)RDB[rea + RLS_DATA_PTR_REA];
                  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

                  /* Get cross section */

                  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
                    fprintf(fp, "%1.5E ", PhotonMicroXS(rea,  E[n], 0));
                  else
                    fprintf(fp, "%1.5E ", MicroXS(rea, E[n], 0));

                  /* 0K data */

                  if ((rea1 = (long)RDB[rea + REACTION_PTR_0K_DATA])
                      > VALID_PTR)
                    fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));

                  /* Majorant */

                  if ((rea1 = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT])
                      > VALID_PTR)
                    fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));
                }

              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");

          /* Print multi-group cross sections */

          if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
            {
              fprintf(fp, "%s_mg_xs = [\n", tmpstr);

              for (n = 0; n < ne; n++)
                {
                  /* Pointer total xs */

                  rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
                  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

                  /* Get cross section */

                  fprintf(fp, "%1.5E ", MGXS(rea, E[n], -1));

                  fprintf(fp, "\n");
                }

              fprintf(fp, "];\n\n");
            }
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Print branch fractions **********************************************/

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while(nuc > VALID_PTR)
    {
      /* Check type */

      if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
          ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BRA_DATA))
        {
          /* Set name */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
            tmpstr[strlen(tmpstr) - 4] = '_';

          /* Print reaction mt's */

          fprintf(fp, "%s_bra_mt = [\n", tmpstr);

          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check pointer to branch data and Print mt */

              if (((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR) &&
                  ((long)RDB[rea + REACTION_RFS] == 1))
                fprintf(fp, "%4ld %% %s\n", (long)RDB[rea + REACTION_MT],
                        ReactionMT((long)RDB[rea + REACTION_MT], YES));

              /* Next */

              rea = NextItem(rea);
            }

          fprintf(fp, "];\n\n");

          /* Print cross sections */

          fprintf(fp, "%s_bra_f = [\n", tmpstr);

          for (n = 0; n < ne; n++)
            {
              /* Loop over reactions */

              rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
              while (rea > VALID_PTR)
                {
                  /* Check pointer to branch data and Print fraction */

                  if (((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR) &&
                      ((long)RDB[rea + REACTION_RFS] == 1))
                    fprintf(fp, "%1.5E ", BranchFrac(rea, 1, E[n], 0));

                  /* Next reaction */

                  rea = NextItem(rea);
                }

              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Print photon yields *************************************************/

  if ((long)RDB[DATA_PHOTON_PRODUCTION] != NO)
    {
      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while(nuc > VALID_PTR)
        {
          /* Check if reactions exist */

          if ((long)RDB[nuc + NUCLIDE_PTR_PHOTON_PROD] < VALID_PTR)
            {
              /* Next nuclide */

              nuc = NextItem(nuc);

              /* Cycle loop */

              continue;
            }

          /* Print mt's */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));
          tmpstr[strlen(tmpstr) - 4] = '_';
          fprintf(fp, "%s_pprod_mt = [\n", tmpstr);

          /* Loop over reactions */

          ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_PROD];
          while (ptr > VALID_PTR)
            {
              fprintf(fp, "%ld\n", (long)RDB[ptr + PHOTON_PROD_MT]);
              ptr = NextItem(ptr);
            }

          fprintf(fp, "];\n\n");

          /* Print cross sections */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));
          tmpstr[strlen(tmpstr) - 4] = '_';
          fprintf(fp, "%s_pprod_xs = [\n", tmpstr);

          for (n = 0; n < ne; n++)
            {
              /* Loop over reactions */

              ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_PROD];
              while (ptr > VALID_PTR)
                {
                  rea1 = (long)RDB[ptr + PHOTON_PROD_PTR_PRODXS];
                  CheckPointer(FUNCTION_NAME, "(rea1)", DATA_ARRAY, rea1);

                  fprintf(fp, "%1.5E ", MicroXS(rea1, E[n], 0));
                  ptr = NextItem(ptr);
                }
              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");

          /* Next */

          nuc = NextItem(nuc);
        }
    }

  /***************************************************************************/

  /***** Print photon line spectra *******************************************/

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while(nuc > VALID_PTR)
    {
      /* Loop over radiations */

      ptr = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
      while (ptr > VALID_PTR)
        {
          /* Check type */

          if ((long)RDB[ptr + NUCLIDE_RAD_TYPE] == PARTICLE_TYPE_GAMMA)
            {
              /* Get pointer to spectrum */

              ptr = (long)RDB[ptr + NUCLIDE_RAD_PTR_SPEC];

              /* NOTE: Tässä pitäisi varmaan olla break (25.5.2018/2.1.31) */

              break;
            }

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Check if line spectra exists */

      if (ptr > VALID_PTR)
        {
          /* Set name */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
            tmpstr[strlen(tmpstr) - 4] = '_';

          /* Print spectra */

          fprintf(fp, "%s_pspec = [\n", tmpstr);

          while (ptr > VALID_PTR)
            {
              /* Check spectrum type */

              if ((long)RDB[ptr + DECAY_SPEC_TYPE] != DECAY_SPEC_LINE)
                {
                  /* Next */

                  ptr = NextItem(ptr);

                  /* Cycle loop */

                  continue;
                }

              /* Print */

              fprintf(fp, "%1.5E %1.5E\n", RDB[ptr + DECAY_SPEC_LINE_E],
                      RDB[ptr + DECAY_SPEC_RI]);

              /* Next line */

              ptr = NextItem(ptr);
            }

          fprintf(fp, "];\n\n");
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Print nubar data ****************************************************/

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while(nuc > VALID_PTR)
    {
      /* Pointer to fission channel */

      if ((rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS]) > VALID_PTR)
        {
          /* Set name */

          sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));

          if ((long)RDB[nuc + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY)
            tmpstr[strlen(tmpstr) - 4] = '_';

          /* Print nubars */

          fprintf(fp, "%s_nu = [\n", tmpstr);

          for (n = 0; n < ne; n++)
            {
              /* Prompt nubar data */

              if ((ptr = (long)RDB[rea + REACTION_PTR_TNUBAR])> VALID_PTR)
                fprintf(fp, "%1.5f ", Nubar(ptr, E[n], 0));

              /* Delayed nubar data */

              if ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR])> VALID_PTR)
                fprintf(fp, "%1.5f ", Nubar(ptr, E[n], 0));

              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Print energy-dependent isomeric branching ratios ********************/

  /* Tää tehtiinkin eri tavalla */

#ifdef mmmmmmmmmmmmmmmmmmmmmm

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check pointer to branching data */

      if ((long)RDB[nuc + NUCLIDE_BRA_TYPE] != BRA_TYPE_ENE)
        {
          /* Next nuclide */

          nuc = NextItem(nuc);

          /* Cycle loop */

          continue;
        }

      /* Set name */

      sprintf(tmpstr, "i%s", GetText(nuc + NUCLIDE_PTR_NAME));
      tmpstr[strlen(tmpstr) - 4] = '_';

      /* Print data */

      fprintf(fp, "%s_bra = [\n", tmpstr);

      for (n = 0; n < ne; n++)
        {
          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check pointer to branching data */

              if ((long)RDB[rea + REACTION_RFS] > 0)
                fprintf(fp, "%1.5E ", MicroXS(rea, E[n], 0));

              /* Next reaction */

              rea = NextItem(rea);
            }

          fprintf(fp, "\n");

        }

      fprintf(fp, "];\n\n");

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

#endif

  /***************************************************************************/

  /***** Print macroscopic cross sections ************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while(mat > VALID_PTR)
    {
      /* Set name */

      sprintf(tmpstr, "m%s", GetText(mat + MATERIAL_PTR_NAME));

      /* Print reaction mt's */

      fprintf(fp, "%s_mt = [\n", tmpstr);

      if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_ELAXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_ABSXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_HEATTXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_PHOTPXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_PROTPXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_DEUTPXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_TRITPXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_HE3PXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_HE4PXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_NSF]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_INLPXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);
      if ((rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS]) > VALID_PTR)
        fprintf(fp, "%4ld\n", (long)RDB[rea + REACTION_MT]);

      fprintf(fp, "];\n\n");

      /* Print cross sections */

      fprintf(fp, "%s_xs = [\n", tmpstr);

      for (n = 0; n < ne; n++)
        {
          /* Add to counter */

          ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
          AddPrivateData(ptr, 1.0, 0);

          if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_ELAXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_ABSXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_HEATTXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_PHOTPXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_PROTPXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_DEUTPXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_TRITPXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_HE3PXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_HE4PXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_NSF]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_INLPXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", PhotonMacroXS(rea, E[n], 0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", PhotonMacroXS(rea, E[n], 0));
          if ((rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS]) > VALID_PTR)
            fprintf(fp, "%1.5E ", MacroXS(rea, E[n],0));

          fprintf(fp, "\n");
        }

      fprintf(fp, "];\n\n");

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** Print multi-group macroscopic xs ************************************/

  /* Check mode */

  if ((long)RDB[DATA_OPTI_MG_MODE] == YES)
    {
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while(mat > VALID_PTR)
        {
          /* Set name */

          sprintf(tmpstr, "m%s", GetText(mat + MATERIAL_PTR_NAME));

          /* Print cross sections */

          fprintf(fp, "%s_mg_xs = [\n", tmpstr);

          for (n = 0; n < ne; n++)
            {
              if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTXS]) > VALID_PTR)
                {
                  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
                    fprintf(fp, "%1.5E ", MGXS(rea, E[n], -1));
                }
              if ((rea = (long)RDB[mat + MATERIAL_PTR_ELAXS]) > VALID_PTR)
                {
                  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
                    fprintf(fp, "%1.5E ", MGXS(rea, E[n], -1));
                }
              if ((rea = (long)RDB[mat + MATERIAL_PTR_ABSXS]) > VALID_PTR)
                {
                  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
                    fprintf(fp, "%1.5E ", MGXS(rea, E[n], -1));
                }
              if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
                {
                  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
                    fprintf(fp, "%1.5E ", MGXS(rea, E[n], -1));
                }
              if ((rea = (long)RDB[mat + MATERIAL_PTR_INLPXS]) > VALID_PTR)
                {
                  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
                    fprintf(fp, "%1.5E ", MGXS(rea, E[n], -1));
                }

              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");

          /* Next material */

          mat = NextItem(mat);
        }
    }

  /***************************************************************************/

  /***** Print Legendre moments for angular distributions ********************/

  if ((long)RDB[DATA_SENS_MODE] != SENS_MODE_NONE)
    {
      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      while(nuc > VALID_PTR)
        {
          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Get pointer to distribution data or skip reaction */

              if ((ang = (long)RDB[rea + REACTION_PTR_ANG]) < VALID_PTR)
                {
                  rea = NextItem(rea);

                  continue;
                }

              /* Legendre moments are currently calculated only for */
              /* tabular type. */

              if ((long)RDB[ang + ANG_TYPE] != ANG_TYPE_TABULAR)
                {
                  rea = NextItem(rea);

                  continue;
                }

              /* Get reaction energy grid */

              erg = (long)RDB[ang + ANG_PTR_EGRID];
              CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

              /* Pointer to data */

              l0 = (long)RDB[ang + ANG_PTR_D0];
              CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

              /* Get number of energy points */

              Ane = (long)RDB[erg + ENERGY_GRID_NE];
              CheckValue(FUNCTION_NAME, "ne", "", ne, 2, 100000);

              /* Get energy grid data */

              ene = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(ene)", DATA_ARRAY, ene);

              /* Set name */

              sprintf(tmpstr, "i%s_MT_%ld", GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT]);

              for (n = 0; n < strlen(tmpstr); n++)
                if (tmpstr[n] == '.')
                  tmpstr[n] = '_';

              /* Print reaction mt's */

              fprintf(fp, "%s_leg = [\n", tmpstr);

              /* Print header line */

              fprintf(fp, "%% %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n",
                      "Energy",
                      "Leg.Mom.P1",
                      "Leg.Mom.P2",
                      "Leg.Mom.P3",
                      "Leg.Mom.P4",
                      "Leg.Mom.P5",
                      "Leg.Mom.P6",
                      "Leg.Mom.P7");

              for (n = 0; n < Ane; n++)
                {
                  /* Print energy */

                  fprintf(fp, " % 1.5E ", RDB[ene + n]);

                  /* Get pointer to data for current energy point */

                  ptr = (long)RDB[l0 + n];

                  /* Get number of cosine points for this point */

                  np = (long)RDB[ptr];

                  /* Print legendre moments */

                  for (m = 0; m < 7; m++)
                    fprintf(fp, "% 1.5E ", RDB[ptr + 3*np + 1 + m]);

                  fprintf(fp, "\n");

                }

              fprintf(fp, "];\n\n");

              /* Next reaction */

              rea = NextItem(rea);
            }

          nuc = NextItem(nuc);
        }
    }

  /***************************************************************************/

  /***** Print majorant ******************************************************/

  /* TODO: muuta tohon ne simulation moodit instead */

  if ((rea = (long)RDB[DATA_PTR_MAJORANT]) > VALID_PTR)
    {
      fprintf(fp, "majorant_xs = [\n");

      for (n = 0; n < ne; n++)
        fprintf(fp, "%1.5E\n", DTMajorant(PARTICLE_TYPE_NEUTRON, E[n],0));

      fprintf(fp, "];\n\n");
    }

  if ((rea = (long)RDB[DATA_PTR_PHOTON_MAJORANT]) > VALID_PTR)
    {
      fprintf(fp, "photon_majorant_xs = [\n");

      for (n = 0; n < ne; n++)
        fprintf(fp, "%1.5E\n", DTMajorant(PARTICLE_TYPE_GAMMA, E[n],0));

      fprintf(fp, "];\n\n");
    }

  /***************************************************************************/

  /* Free memory */

  Mem(MEM_FREE, E);

  /* Close file */

  fclose(fp);

  /* Restore collision count */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  PutPrivateData(ptr, ncol, 0);

  /* Exit */

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
