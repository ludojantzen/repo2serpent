/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsensitivities.c                         */
/*                                                                           */
/* Created:       2017/05/04 (VVa)                                           */
/* Last modified: 2018/09/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Does all sensitivity calculation processing.                 */
/*                                                                           */
/* Comments: -Memory for statistics is allocated in processsensstats.c       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSensitivities:"

/*****************************************************************************/

void ProcessSensitivities()
{
  long ptr, zai, nuc, lst, rls, rea, *zaiarr, nzai, i, j;
  long nmt, *mtarr, mt, n, maxi, maxmt, det, tbins, rbins;
  long nrea, nene, nmu, np, sens, loc0, loc1, loc2, ptr1, ptr2, type, idx;
  long ene, mat, nmat, include, includesum, msh, bin, arr, resp;
  double E, val, tmpdbl, P;
  char tmpstr[MAX_STR], fname[MAX_STR], name[MAX_STR];
  FILE *fp;
  long reaindices[] = {
    ELA_SCATT_IDX,
    SAB_SCATT_IDX,
    INL_SCATT_IDX,
    CAPT_IDX,
    FISS_IDX,
    NXN_SCATT_IDX,
    TEMPERATURE_IDX,
    SCATT_TEMP_IDX,
    NUBAR_TOT_IDX,
    NUBAR_PRO_IDX,
    NUBAR_DEL_IDX,
    CHI_TOT_IDX,
    CHI_PRO_IDX,
    CHI_DEL_IDX,
    ELA_P1_IDX,
    ELA_P2_IDX,
    ELA_P3_IDX,
    ELA_P4_IDX,
    ELA_P5_IDX,
    ELA_P6_IDX,
    ELA_P7_IDX,
    0
  };
  char *reanames[] = {
    "ela",
    "sab",
    "inl",
    "capt",
    "fiss",
    "nxn",
    "temperature",
    "nuTot",
    "nuPro",
    "nuDel",
    "chiTot",
    "chiPro",
    "chiDel",
    "elaP1",
    "elaP2",
    "elaP3",
    "elaP4",
    "elaP5",
    "elaP6",
    "elaP7",
    "\0"
  };
  char *sumreanames[] = {
    "ela",
    "sab",
    "inl",
    "capt",
    "fiss",
    "nxn",
    "\0"
  };


  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  fprintf(outp, "Processing sensitivity calculation data...\n");

  /****************************************************************/
  /* Link the energy group structure used for the Sens calculation */
  /****************************************************************/

  /* Check pointer to text */

  if ((long)RDB[sens + SENS_PTR_EGRID] < VALID_PTR)
    Die(FUNCTION_NAME, "Sensitivity calculation energy grid not set, use \"set sensopt egrid <name\"");

  /* Get pointer to energy grid */

  ene = (long)RDB[DATA_PTR_ENE0];

  /* Loop over energy grid definition */

  while (ene > VALID_PTR)
    {

      /* Compare energy grid names */

      if (CompareStr(sens + SENS_PTR_EGRID, ene + ENE_PTR_NAME))
        break;

      /* Next energy grid definition */

      ene = NextItem(ene);
    }

  /* Check if found */

  if (ene < VALID_PTR)
    {
      /* User defined energy grid not found */

      Die(FUNCTION_NAME, "Energy grid \"%s\" used for sensitivity calculation not defined",
          GetText(sens + SENS_PTR_EGRID));
    }
  else
    {
      /* Energy grid was found, link grid */

      WDB[sens + SENS_PTR_EGRID] = RDB[ene + ENE_PTR_GRID];

      /* Put number of groups */

      WDB[sens + SENS_N_ENE] = RDB[ene + ENE_NB];
    }


  if ((long)RDB[sens + SENS_PTR_MESH] > VALID_PTR)
    {
      /* Get pointer to energy grid */

      msh = (long)RDB[DATA_PTR_MESH0];

      /* Loop over energy grid definition */

      while (msh > VALID_PTR)
        {

          /* Compare energy grid names */

          if (CompareStr(sens + SENS_PTR_MESH, msh + MESH_PTR_NAME))
            break;

          /* Next energy grid definition */

          msh = NextItem(msh);
        }

      /* Check if found */

      if (msh < VALID_PTR)
        {
          /* User defined energy grid not found */

          Die(FUNCTION_NAME, "Datamesh \"%s\" used for sensitivity calculation not defined",
              GetText(sens + SENS_PTR_MESH));
        }
      else
        {
          /* Energy grid was found, link grid */

          WDB[sens + SENS_PTR_MESH] = msh;

          /* Check number of groups */

          if ((long)RDB[sens + SENS_N_ENE] !=  (long)RDB[msh + MESH_N0]*(long)RDB[msh + MESH_N1]*(long)RDB[msh + MESH_N2])
            Die(FUNCTION_NAME, "Must have same number of bins in ene-mesh and data-mesh");
        }
    }

  /****************************************************************************/
  /**************** Materials *************************************************/
  /****************************************************************************/

  /******************/
  /* Link materials */
  /******************/

  if ((ptr = (long)RDB[sens + SENS_PTR_MAT_ARR]) > VALID_PTR)
    {
      /* Loop over material list */

      while ((long)RDB[ptr] > VALID_PTR)
        {

          /* Find this material from global materials */

          mat = (long)RDB[DATA_PTR_M0];
          if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, GetText(ptr)))
              < VALID_PTR)
            Die(FUNCTION_NAME, "Unable to find material \"%s\" linked to "
                "sensitivity calculations."
                , GetText(ptr));

          /* Swap material name pointer for material pointer */

          WDB[ptr] = (double)mat;

          /* Next material */

          ptr++;
        }
    }
  else if (ptr == -1)
    {
      /* Include all materials (parents) */

      /* Loop over materials to count them */

      mat = (long)RDB[DATA_PTR_M0];
      nmat = 0;
      while (mat > VALID_PTR)
        {

          if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR)
            {
              /* This is not a child material  */
              /* Increment number of materials */

              nmat++;
            }

          /* Next nuclide */

          mat = NextItem(mat);
        }

      /* Store number of linked materials */

      WDB[sens + SENS_N_MAT] = (double)nmat;

      /* Allocate memory for material list */

      ptr = ReallocMem(DATA_ARRAY, nmat + 1);

      /* Store pointer to material list */

      WDB[sens + SENS_PTR_MAT_ARR] = ptr;

      /* Loop over materials to store them to material list */

      mat = (long)RDB[DATA_PTR_M0];

      while (mat > VALID_PTR)
        {

          if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] < VALID_PTR)
            {
              /* This is not a child material  */
              /* Store material */

              WDB[ptr++] = (double)mat;
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Put nullpointer */

      WDB[ptr] = NULLPTR;

    }

  /********************************************/
  /* Set up the indices for the sum materials */
  /********************************************/

  /* Allocate memory for index mapping */

  ptr = ReallocMem(DATA_ARRAY, N_SENS_MAT);

  /* Store pointer */

  WDB[sens + SENS_PTR_MAT_INDICES] = (double)ptr;

  /* Increment total number of materials for "lost" */

  WDB[sens + SENS_N_MAT] += (double)1.0;

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
    {
      /* Incerement total number of materials for "total" */

      WDB[sens + SENS_N_MAT] += (double)1.0;

      /* Increment index for "lost" material due to total */

      WDB[ptr + SENS_NON_MAT_IDX] += 1.0;

      /* Increment index for "sum" material due to total */

      WDB[ptr + SENS_SUM_MAT_IDX] += 1.0;
    }

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_SUM)
    {
      /* Incerement total number of materials for "sum" */

      WDB[sens + SENS_N_MAT] += (double)1.0;

      /* Increment index for "lost" material due to sum */

      WDB[ptr + SENS_NON_MAT_IDX] += 1.0;
    }

  /* Check that at least one material is included in addition to "lost" */

  if ((long)RDB[sens + SENS_N_MAT] <= 1)
    {
      Warn(FUNCTION_NAME, "No materials linked for sensitivity calculation, including total.");
      SetOption(sens + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_MAT_TOT);
      WDB[sens + SENS_N_MAT] += (double)1.0;

      /* Increment index for "lost" material due to total */

      WDB[ptr + SENS_NON_MAT_IDX] += 1.0;
    }

  /****************************************************************************/
  /**************** ZAIs ******************************************************/
  /****************************************************************************/

  /***************/
  /* Check ZAIs  */
  /***************/

  if ((ptr = (long)RDB[sens + SENS_PTR_ZAI_ARR]) > VALID_PTR)
    {
      /* Loop over ZAI list */

      while ((zai = (long)RDB[ptr]) >= 10010)
        {
          /* Loop over nuclides */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {

              /* Check if found */

              if ((long)RDB[nuc + NUCLIDE_ZAI] == zai)
                break;

              /* Not found, check next material */

              nuc = NextItem(nuc);
            }

          /* Check if found */

          if (nuc < VALID_PTR)
            Die(FUNCTION_NAME, "Could not find ZAI %ld linked to sensitivity calculation from "
                "nuclides present in calculation.", zai);

          /* Next ZAI */
          ptr++;
        }
    }
  else if (ptr == -1)
    {
      /* Include all ZAIs */

      /* Create an array for unique ZAIs */

      zaiarr = (long *)Mem(MEM_ALLOC, 1, sizeof(long));

      /* Reset number of unique ZAIs */

      nzai = 0;

      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];

      while (nuc > VALID_PTR)
        {
          /* Get zai */

          zai = (long)RDB[nuc + NUCLIDE_ZAI];

          /* Check if the material name has been stored */

          for (i = 0; i < nzai; i++)
            {
              if (zaiarr[i] == zai)
                break;
            }

          if (i >= nzai)
            {

              /* New zai */

              /* Reallocate zai array */

              zaiarr = (long*)Mem(MEM_REALLOC, zaiarr,  (nzai + 1)*sizeof(long));

              /* Store zai */

              zaiarr[nzai] = zai;

              /* Increment number of found ZAIs */

              nzai++;
            }

          /* Next nuclide */

          nuc = NextItem(nuc);
        }

      /* Store number of linked ZAIs */

      WDB[sens + SENS_N_ZAI] = (double)nzai;

      /* Allocate memory for ZAI list */

      ptr = ReallocMem(DATA_ARRAY, nzai + 1);

      /* Store ZAI list */
      fprintf(outp, "Storing %ld for ZAI array\n", ptr);
      WDB[sens + SENS_PTR_ZAI_ARR] = ptr;

      /* Store ZAIs */

      for (i = 0; i < nzai; i++)
        WDB[ptr++] = zaiarr[i];

      /* Put nullpointer */

      WDB[ptr] = NULLPTR;

      /* Free temporary array */

      Mem(MEM_FREE, zaiarr);
    }

  /***************************************/
  /* Set up the indices for the sum zais */
  /***************************************/

  /* Allocate memory for index mapping */

  ptr = ReallocMem(DATA_ARRAY, N_SENS_ZAI);

  /* Store pointer */

  WDB[sens + SENS_PTR_ZAI_INDICES] = (double)ptr;

  /* Increment total number of zais for "lost" */

  WDB[sens + SENS_N_ZAI] += (double)1.0;

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_TOT)
    {
      /* Incerement total number of zais for "total" */

      WDB[sens + SENS_N_ZAI] += (double)1.0;

      /* Increment index for "lost" zai due to total */

      WDB[ptr + SENS_NON_ZAI_IDX] += 1.0;

      /* Increment index for "sum" zai due to total */

      WDB[ptr + SENS_SUM_ZAI_IDX] += 1.0;
    }

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_SUM)
    {
      /* Incerement total number of zais for "sum" */

      WDB[sens + SENS_N_ZAI] += (double)1.0;

      /* Increment index for "lost" zai due to sum */

      WDB[ptr + SENS_NON_ZAI_IDX] += 1.0;
    }

  /* Check that at least one ZAI is included in addition to "lost" */

  if ((long)RDB[sens + SENS_N_ZAI] <= 1)
    {
      Warn(FUNCTION_NAME, "No zais linked for sensitivity calculation, including total");
      SetOption(sens + SENS_SCORE_FLAGS, SENS_SCORE_FLAG_ZAI_TOT);
      WDB[sens + SENS_N_ZAI] += (double)1.0;

      /* Increment index for "lost" zai due to total */

      WDB[ptr + SENS_NON_ZAI_IDX] += 1.0;
    }

  /****************************************************************************/
  /**************** Perturbations *********************************************/
  /****************************************************************************/

  /**************************/
  /* Standard perturbations */
  /**************************/

  /* Allocate memory for index mapping */

  ptr = ReallocMem(DATA_ARRAY, N_SENS_PERT);

  /* Store pointer */

  WDB[sens + SENS_PTR_PERT_INDICES] = (double)ptr;

  /* Reset current index */

  nrea = 0;

  /* Put total index */

  WDB[ptr + TOT_REA_IDX  ] = nrea++;

  /* Put indices for normal XS perturbations */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS)
    {
      /* Check if a list of sum reaction modes was given */

      if ((ptr1 = (long)RDB[sens + SENS_PTR_REA_LIST])  > VALID_PTR)
        {
          /* Was given, include only some of the sum reactions */

          while ((long)RDB[ptr1] > VALID_PTR)
            {
              /* Get name of sum reaction */

              sprintf(name, "%s", GetText(ptr1));

              /* Compare reaction name to allowed reaction names */

              n = 0;
              while (sumreanames[n][0] != '\0')
                if (!strcasecmp(sumreanames[n++], name))
                  {
                    n--;
                    break;
                  }

              /* Check if reaction was not found */

              if (sumreanames[n][0] == '\0')
                {
                  /* Print list of allowed reaction names */

                  sprintf(tmpstr, " Allowed sum reaction mode names:");

                  n = 0;

                  while (sumreanames[n][0] != '\0')
                    sprintf(tmpstr + strlen(tmpstr), " %s", sumreanames[n++]);

                  /* Die with input error */

                  Die(FUNCTION_NAME, "Unknown sum reaction mode name \"%s\" for \"sens pert xs realist\". %s", name, tmpstr);
                }

              /* Get reaction index */

              idx = reaindices[n];

              /* Check if reaction was doubly defined */

              if (RDB[ptr + idx] > 0)
                {
                  /* Was previously set on, give warning */

                  Warn(FUNCTION_NAME, "Sum reaction mode \"%s\" included multiple "
                       "times in \"sens pert xs realist\".", name);
                }
              else
                {
                  /* Was not previously set on, set it on now  */

                  WDB[ptr + idx] = nrea++;

                  /* Print output */

                  fprintf(outp, "Including perturbation of sum reaction mode \"%s\".\n", name);
                }

              /* Next reaction name */

              ptr1++;
            }
        }
      else
        {
          /* Sum reaction list was not given, include all of them */

          n = 0;

          while (sumreanames[n][0] != '\0')
            {
              /* Print output */

              fprintf(outp, "Including perturbation of sum reaction mode \"%s\".\n", sumreanames[n]);

              /* Get index of sum reaction n */

              idx = reaindices[n];

              /* Set perturbation on */

              WDB[ptr + idx] = nrea++;

              /* Increment n */

              n++;
            }
        }
    }
  else if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
    {

      /********************************************/
      /* Process perturbed partial reaction modes */
      /********************************************/

      /* Figure out all unique MT's */

      /* Create an array for unique ZAIs */

      mtarr = (long *)Mem(MEM_ALLOC, 1, sizeof(long));

      /* Reset number of unique MTs */

      nmt = 0;

      /* Reset maximum MT */

      maxmt = 0;

      /* Extend maximum mt number due to chi or nubar sensitivity */

      if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
        maxmt = 1018;
      else if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
        maxmt = 456;

      /* Loop over nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];

      while (nuc > VALID_PTR)
        {
          /* Get reaction list */

          lst = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
          CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

          /* Loop over partial reactions */

          j = 0;
          while ((rls = ListPtr(lst, j++)) > VALID_PTR)
            {
              /* Pointer to reaction data */

              rea = (long)RDB[rls + RLS_DATA_PTR_REA];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              mt = (long)RDB[rea + REACTION_MT];

              /* Check if the MT has been stored */

              for (i = 0; i < nmt; i++)
                {
                  if (mtarr[i] == mt)
                    break;
                }

              if (i >= nmt)
                {

                  /* New mt */

                  /* Reallocate mt array */

                  mtarr = (long*)Mem(MEM_REALLOC, mtarr,  (nmt + 1)*sizeof(long));

                  /* Store mt */

                  mtarr[nmt] = mt;

                  /* Check maximum */

                  if (mt > maxmt)
                    maxmt = mt;

                  /* Increment number of found mts */

                  nmt++;
                }

              /* Next reaction */
            }

          /* Next nuclide */

          nuc = NextItem(nuc);
        }

      /* Store maximum MT */

      WDB[sens + SENS_MAX_MT] = (double)maxmt;

      /* Allocate memory for MT <-> index list */

      ptr = ReallocMem(DATA_ARRAY, maxmt + 1);

      /* Store pointer to list */

      WDB[sens + SENS_PTR_MT_INDICES] = (double)ptr;

      /* Loop over all possible MT's up to maxmt to link */

      for (mt = 0; mt < maxmt + 1; mt++)
        {
          /* Loop over mtarr or pre-defined list to figure out whether to link this mt */

          include = 0;
          includesum = 0;

          if  ((ptr1 = (long)RDB[sens + SENS_PTR_MT_LIST]) > VALID_PTR)
            {
              /* User given list exists, loop over it to see if mt is on that list */

              while ((long)RDB[ptr1] != NULLPTR)
                {
                  /* Break loop if found */

                  if (mt == (long)RDB[ptr1])
                    break;

                  /* Get next mt from user given list */

                  ptr1++;
                }

              /* Check if mt was present on list */

              if ((long)RDB[ptr1] != NULLPTR)
                {
                  /* mt was on user given list */

                  /* Warn if mt is not actually present in any nuclide */

                  for (i = 0; i < nmt; i++)
                    {
                      if (mtarr[i] == mt)
                        break;
                    }

                  /* Check if mt was present in some nuclide */

                  if (i < nmt)
                    include = 1;
                  else
                    {
                      /* Check if mt is one for a sum reaction mode */

                      if ( (mt == 4) || (mt == 18))
                        includesum = 1;
                      else
                        Warn(FUNCTION_NAME, "MT %ld set to be perturbed but not "
                             "actually present for any nuclide.", mt);

                    }
                }
            }
          else
            {
              /* User given list does not exist, link all mt's present */

              for (i = 0; i < nmt; i++)
                {
                  if (mtarr[i] == mt)
                    break;
                }

              /* Check if mt was present in some nuclide */

              if (i < nmt)
                include = 1;
            }

          /* Switch perturbation for mt on or off based on whether to include */

          if (include || includesum)
            {
              /* Print output */

              if (includesum)
                fprintf(outp, "Including perturbation of sum reaction mode mt %ld.\n", mt);
              else
                fprintf(outp, "Including perturbation of reaction mode mt %ld.\n", mt);

              /* Was found */

              WDB[ptr + mt] = (double)nrea++;
            }
          else
            {
              /* Was not found */

              WDB[ptr + mt] = (double)0;
            }
        }
    }

  /* Store largest index for XS perturbation */

  WDB[sens + SENS_MAX_XS_INDEX] = (double)(nrea-1);

  /* Get pointer for basic perturbation indices */

  ptr = (long)RDB[sens + SENS_PTR_PERT_INDICES];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Put index for temperature perturbation */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE)
    {
      WDB[ptr + TEMPERATURE_IDX] = nrea++;
      WDB[ptr + SCATT_TEMP_IDX] = nrea++;
    }

  /* Put indices for nubar perturbations */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
    {
      WDB[ptr + NUBAR_TOT_IDX] = nrea++;
      WDB[ptr + NUBAR_PRO_IDX] = nrea++;
      WDB[ptr + NUBAR_DEL_IDX] = nrea++;

      /* Put indices also to the mt index array */

      if ((ptr1 = (long)RDB[sens + SENS_PTR_MT_INDICES]) > VALID_PTR)
        {
          /* Check maximum MT number in index array */

          if ((long)RDB[sens + SENS_MAX_MT] < 456)
            Die(FUNCTION_NAME, "No space for nubar indices in MT index "
                "array. This should have been allocated...");

          /* Store reaction indices */

          WDB[ptr1 + 452] = RDB[ptr + NUBAR_TOT_IDX];
          WDB[ptr1 + 455] = RDB[ptr + NUBAR_DEL_IDX];
          WDB[ptr1 + 456] = RDB[ptr + NUBAR_PRO_IDX];
        }
    }

  /* Put indices for fission spectrum perturbations */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
    {
      WDB[ptr + CHI_TOT_IDX] = nrea++;
      WDB[ptr + CHI_PRO_IDX] = nrea++;
      WDB[ptr + CHI_DEL_IDX] = nrea++;

      /* Put indices also to the mt index array */

      if ((ptr1 = (long)RDB[sens + SENS_PTR_MT_INDICES]) > VALID_PTR)
        {
          /* Check maximum MT number in index array */

          if ((long)RDB[sens + SENS_MAX_MT] < 1018)
            Die(FUNCTION_NAME, "No space for chi indices in MT index "
                "array. This should have been allocated...");

          /* Not an official MT number, but used in e.g. SCALE */
          /* covariance data */

          WDB[ptr1 + 1018] = RDB[ptr + CHI_TOT_IDX];
        }
    }

  /* Put indices for elastic mu perturbation */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_ELA_MU)
    {
      WDB[ptr + ELA_MU_IDX] = nrea++;
    }

  /* Put indices for inelastic mu perturbation */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_INL_MU)
    {
      WDB[ptr + INL_MU_IDX] = nrea++;
    }

  /* Put indices for elastic legendre moment perturbations */

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM)
    {
      maxi = (long)RDB[sens + SENS_MAX_SCATT_MOM];

      for (i = 0; i < maxi; i++)
        WDB[ptr + ELA_P1_IDX + i] = nrea++;
    }

  /*************************************/
  /* User defined custom perturbations */
  /*************************************/

  /* Get pointer to first custom perturbation */

  loc0 = (long)RDB[sens + SENS_PTR_PERT0];

  /* Loop over custom perturbations */

  while (loc0 > VALID_PTR)
    {
      /* Print perturbation name */

      fprintf(outp, "Custom perturbation %s:\n", GetText(loc0 + SENS_PERT_PTR_NAME));

      /* Check perturbation energy data */

      if ((long)WDB[loc0 + SENS_PERT_PTR_EGRID] > VALID_PTR)
        {
          /* Copy energy grid data file name */

          snprintf(fname, MAX_STR, "%s", GetText(loc0 + SENS_PERT_PTR_EGRID));

          /* Read energy grid data from file */

          fp = fopen(fname, "r");

          /* Check if file open failed */

          if (fp == NULL)
            Error(loc0, "Could not open energy data file \"%s\"\n", fname);

          /* Get interpolation string from file */

          np = fscanf(fp, "%s", tmpstr);

          if (np != 1)
            Error(loc0, "Could not read interpolation string from file \"%s\"",
                  fname);

          /* Check interpolation string */

          if (!strcasecmp(tmpstr, "lin-lin"))
            WDB[loc0 + SENS_PERT_EGRID_INTERP] = (double)ENDF_INTERP_TYPE_LIN_LIN;
          else if (!strcasecmp(tmpstr, "lin-log"))
            WDB[loc0 + SENS_PERT_EGRID_INTERP] = (double)ENDF_INTERP_TYPE_LIN_LOG;
          else if (!strcasecmp(tmpstr, "log-lin"))
            WDB[loc0 + SENS_PERT_EGRID_INTERP] = (double)ENDF_INTERP_TYPE_LOG_LIN;
          else if (!strcasecmp(tmpstr, "log-log"))
            WDB[loc0 + SENS_PERT_EGRID_INTERP] = (double)ENDF_INTERP_TYPE_LOG_LOG;
          else
            Error(loc0, "Invalid interpolation mode \"%s\" read from file \"%s\"",
                  tmpstr, fname);

          /* Get number of value pairs */

          np = fscanf(fp, "%ld", &n);
          CheckValue(FUNCTION_NAME, "n", "Number of value pairs", n, 2, INFTY);

          if (np != 1)
            Error(loc0, "Could not read number of value pairs from file \"%s\"",
                  fname);

          /* Create energy grid structure */
          WDB[loc0 + SENS_PERT_PTR_EGRID] = 0.0;
          loc1 = NewItem(loc0 + SENS_PERT_PTR_EGRID, ENE_BLOCK_SIZE);

          /* Put energy grid name */

          sprintf(tmpstr, "Custom perturbation %s energy grid.", GetText(loc0 + SENS_PERT_PTR_NAME));
          WDB[loc1 + ENE_PTR_NAME] = PutText(tmpstr);

          /* Put number of bins for energy grid */

          WDB[loc1 + ENE_NB] = n - 1;

          /* Allocate memory for grid */

          ptr1 = ReallocMem(DATA_ARRAY, n);

          /* Store pointer */

          WDB[loc1 + ENE_PTR_GRID] = ptr1;

          /* Allocate memory for values */

          ptr2 = ReallocMem(DATA_ARRAY, n);

          /* Store pointer */

          WDB[loc0 + SENS_PERT_PTR_EGRID_VAL] = ptr2;

          /* Loop over grid points */

          tmpdbl = -1.0;

          for (i = 0; i < n; i++)
            {
              /* Get value pair */

              np = fscanf(fp, "%lf %lf", &E, &val);

              if (np != 2)
                Error(loc0, "Could not read value pair %ld from file \"%s\"",
                      i+1, fname);

              /* Check that energy is larger than previous */

              if (E < tmpdbl)
                Error(loc0, "Energies should be in ascending order in file "
                      "\"%s\". Got %E after %E.", fname, E, tmpdbl);

              /* Store values */

              WDB[ptr1 + i] = E;
              WDB[ptr2 + i] = val;

              /* Store previous energy */

              tmpdbl = E;
            }

          /* Get energy grid interpolation */

          type = (long)RDB[loc0 + SENS_PERT_EGRID_INTERP];

          /* Create energy grid */

          if ((type == ENDF_INTERP_TYPE_LOG_LIN) ||
              (type == ENDF_INTERP_TYPE_LOG_LIN))
            ptr = MakeEnergyGrid(n, 0, 0, -1, &RDB[ptr1], EG_INTERP_MODE_LOG);
          else
            ptr = MakeEnergyGrid(n, 0, 0, -1, &RDB[ptr1], EG_INTERP_MODE_LIN);

          /* Put pointer */

          WDB[loc1 + ENE_PTR_GRID] = (double)ptr;

          /* Put minimum and maximum energies */

          WDB[ene + ENE_EMIN] = RDB[ptr + ENERGY_GRID_EMIN];
          WDB[ene + ENE_EMAX] = RDB[ptr + ENERGY_GRID_EMAX];

          /* Close file */

          fclose(fp);
        }

      /* Check that linked ZAIs are found */

      if ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_ZAI_LIST]) > VALID_PTR)
        {
          /* Loop over ZAI list */

          fprintf(outp, "ZAIs linked to perturbation:");

          while ((zai = (long)RDB[ptr]) >= 10010)
            {
              /* Loop over nuclides */

              nuc = (long)RDB[DATA_PTR_NUC0];
              while (nuc > VALID_PTR)
                {

                  /* Check if found */

                  if ((long)RDB[nuc + NUCLIDE_ZAI] == zai)
                    break;

                  /* Not found, check next material */

                  nuc = NextItem(nuc);
                }

              /* Check if found */

              if (nuc < VALID_PTR)
                Error(loc0, "Could not find ZAI %ld linked to custom perturbation from "
                    "nuclides present in calculation.", zai);

              /* Print ZAI */

              fprintf(outp, " %ld", zai);

              /* Next ZAI */

              ptr++;
            }

          fprintf(outp, "\n");
        }

      /* Swap material names to material pointers */

      if ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_MAT_LIST]) > VALID_PTR)
        {
          /* Loop over list */

          fprintf(outp, "Materials linked to perturbation:");

          /* Loop over material list */

          while ((long)RDB[ptr] > VALID_PTR)
            {
              /* Find this material from global materials */

              mat = (long)RDB[DATA_PTR_M0];
              if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, GetText(ptr)))
                  < VALID_PTR)
                Error(loc0, "Could not find material \"%s\" linked to "
                      "custom perturbation from "
                      "materials present in calculation.", GetText(ptr));

              /* Print material name */

              fprintf(outp, " %s", GetText(ptr));

              /* Swap material name pointer for material pointer */

              WDB[ptr] = (double)mat;

              /* Next material */

              ptr++;
            }

          fprintf(outp, "\n");
        }

      /* Swap reaction names to reaction indices */

      if ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_REA_LIST]) > VALID_PTR)
        {
          fprintf(outp, "Reactions linked to perturbation:");

          /* Loop over list */

          while ((long)RDB[ptr] > VALID_PTR)
            {
              /* Compare reaction name to allowed reaction names */

              n = 0;
              while (reanames[n][0] != '\0')
                if (!strcasecmp(reanames[n++], GetText(ptr)))
                  {
                    n--;
                    break;
                  }


              /* Check if reaction was not found */

              if (reanames[n][0] == '\0')
                {
                  /* Print list of allowed reaction names */

                  sprintf(tmpstr, " Allowed reaction mode names:");

                  n = 0;

                  while (reanames[n][0] != '\0')
                    sprintf(tmpstr + strlen(tmpstr), " %s", reanames[n++]);

                  fprintf(outp, ".\n");

                  /* Die with input error */

                  Error(loc0, "Unknown reaction name \"%s\". %s", GetText(ptr), tmpstr);
                }

              /* Print reaction name */

              fprintf(outp, " %s", GetText(ptr));

              /* Swap raction name for reaction index */

              WDB[ptr] = (double)reaindices[n];

              /* Next item */

              ptr++;
            }

          fprintf(outp, "\n");
        }

      /* Put index of current perturbation */

      WDB[loc0 + SENS_PERT_INDEX] = nrea++;

      fprintf(outp, "\n");

      /* Next custom perturbation */

      loc0 = NextItem(loc0);
    }

  /* Store total number of perturbations */

  WDB[sens + SENS_N_PERT] = (double)(nrea - 1);

  /*******************************************/
  /* Link void fraction sensitivity material */
  /*******************************************/

  if ((long)RDB[sens + SENS_PTR_VOID_MAT] > VALID_PTR)
    {
      /* Find this material from global materials */

      mat = (long)RDB[DATA_PTR_M0];
      if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, GetText(sens + SENS_PTR_VOID_MAT)))
          < VALID_PTR)
        Die(FUNCTION_NAME, "Unable to find material \"%s\" linked to "
            "void coefficient sensitivity calculation."
            , GetText(sens + SENS_PTR_VOID_MAT));

      /* Swap material name pointer for material pointer */

      WDB[sens + SENS_PTR_VOID_MAT] = (double)mat;
    }

  /****************************************************/
  /* Link covariance data for uncertainty propagation */
  /****************************************************/

#ifdef XYZ
  if ((loc0 = (long)RDB[DATA_PTR_COVMTX0]) > VALID_PTR)
    {
      /* Link these directly to reactions */

      /* Get ZAIs and mts */

      ZAI1 = (long)RDB[loc0 + COVMTX_ZAI1];
      ZAI2 = (long)RDB[loc0 + COVMTX_ZAI2];

      mt1 = (long)RDB[loc0 + COVMTX_MT1];
      mt2 = (long)RDB[loc0 + COVMTX_MT2];

      /* Loop over nuclides to link all correct reactions */

      nuc = (long)RDB[DATA_PTR_NUC0];

      while (nuc > VALID_PTR)
        {
          /* Get ZAI from nuclide */

          ZAI = (long)RDB[nuc + NUCLIDE_ZAI];

          /* Check if ZAI matches one of the ones in covariance matrix */

          if ((ZAI != ZAI1) && (ZAI != ZAI2))
            {
              /* Get next nuclide */

              nuc = NextItem(nuc);
              continue;
            }

          /* Nuclide is a part of the covariance matrix */


        }

      loc0 = NextItem(loc0);
    }
#endif

  /****************************************************************************/
  /******** Detector ratios ***************************************************/
  /****************************************************************************/

  nmat = (long)RDB[sens + SENS_N_MAT];
  nzai = (long)RDB[sens + SENS_N_ZAI];
  nrea = (long)RDB[sens + SENS_N_PERT] + 1;
  nene = (long)RDB[sens + SENS_N_ENE] + 1;
  nmu = (long)RDB[sens + SENS_N_MU];

  maxi = (1 + (nmat*nzai*nrea*nene*nmu));

  /* Store maximum possible label */

  WDB[sens + SENS_MAX_LABEL] = (double)maxi;

  /* Get last generation to be scored */

  np = (long)RDB[DATA_SENS_LAST_GEN];

  /* Get first response */

  loc0 = (long)RDB[sens + SENS_PTR_RESP0];

  loc1 = -1;

  while (loc0 > VALID_PTR)
    {
      /* Skip non-ratio responses here */

      if ((long)RDB[loc0 + SENS_RESP_TYPE] != SENS_RESP_TYPE_RATIO)
        {
          loc0 = NextItem(loc0);
          continue;
        }

      /* Skip responses added in this loop */

      if ((long)RDB[loc0 + SENS_RESP_DET_BIN_IDX] >= 0)
        {
          loc0 = NextItem(loc0);
          continue;
        }

      /* Process and link both detectors */

      for (i = 0; i < 2; i++)
        {
          switch (i)
            {
            case 0:
              loc1 = loc0 + SENS_RESP_PTR_DET1;
              loc2 = loc0 + SENS_RESP_DET1_N_BIN;
              P = RDB[loc0 + SENS_RESP_DET1_P];
              break;
            case 1:
              loc1 = loc0 + SENS_RESP_PTR_DET2;
              loc2 = loc0 + SENS_RESP_DET2_N_BIN;
              P = RDB[loc0 + SENS_RESP_DET2_P];
              break;
            default:
              Die(FUNCTION_NAME, "WTF!");
            }

          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Find detector from memory */

          det = (long)RDB[DATA_PTR_DET0];

          if ((det = SeekListStr(det, DET_PTR_NAME,
                                 GetText(loc1)))
              < VALID_PTR)
            Die(FUNCTION_NAME, "Unable to find detector \"%s\" linked to "
                "sensitivity calculations."
                , GetText(loc1));

          /* Check that detector has exactly one bin */

          tbins = (long)RDB[det + DET_N_TOT_BINS];
          rbins = (long)RDB[det + DET_N_RBINS];

          if (rbins != 1)
            Error(loc0, "Detector %s linked to sensitivity calculation has more than one response",
                  GetText(det + DET_PTR_NAME));

          /* Store number of bins in detector */

          WDB[loc2] = (double)tbins;

          /* Link detector to sensitivity response */

          WDB[loc1] = (double)det;

        }

      /* Check that the detector have an equal number of bins */

      if ((long)RDB[loc0 + SENS_RESP_DET1_N_BIN] !=
          (long)RDB[loc0 + SENS_RESP_DET2_N_BIN])
        Error(loc0, "Detectors linked to sensitivity calculation have a different"
              "number of bins (%ld vs %ld)",
              (long)RDB[loc0 + SENS_RESP_DET1_N_BIN],
              (long)RDB[loc0 + SENS_RESP_DET2_N_BIN]);

      tbins = (long)RDB[loc0 + SENS_RESP_DET1_N_BIN];

      /* Create new responses for the different bins */

      for (bin = 0; bin < tbins; bin++)
        {
          /* Create a new response for this bin */

          resp = NewItem(sens + SENS_PTR_RESP0, SENS_RESP_BLOCK_SIZE);

          /* Copy basic information from the parent response including  */
          /* PARAM_N_COMMON */

          memcpy(&WDB[resp + LIST_DATA_SIZE], &RDB[loc0 + LIST_DATA_SIZE],
                 (SENS_RESP_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

          /* Store the detector bin number for this response */

          WDB[resp + SENS_RESP_DET_BIN_IDX] = bin;

          /* Rename this response */

          sprintf(tmpstr, "%s_BIN_%ld", GetText(resp + SENS_RESP_PTR_NAME), bin);

          WDB[resp + SENS_RESP_PTR_NAME] = PutText(tmpstr);

          /* Allocate memory for the buffers and statistics */

          for (i = 0; i < 2; i++)
            {
              switch (i)
                {
                case 0:
                  det = (long)RDB[resp + SENS_RESP_PTR_DET1];
                  P = RDB[resp + SENS_RESP_DET1_P];
                  break;
                case 1:
                  det = (long)RDB[resp + SENS_RESP_PTR_DET2];
                  P = RDB[resp + SENS_RESP_DET2_P];
                  break;
                default:
                  Die(FUNCTION_NAME, "WTF!");
                }

              CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

              /* Allocate statistic array if needed */

              if ((arr = (long)RDB[det + DET_PTR_SENS_STAT_ARRAY]) < VALID_PTR)
                {
                  /* Allocate memory for the sensitivity statistic array */

                  arr = ReallocMem(DATA_ARRAY, tbins);
                  WDB[det + DET_PTR_SENS_STAT_ARRAY] = (double)arr;
                }

              /* Allocate memory for the sensitivity statistic if needed */

              if ((long)RDB[arr + bin] < VALID_PTR)
                {
                  /* Allocate memory for the sensitivity statistic array */

                  sprintf(tmpstr, "ADJ_PERT_DET_%s_SENS", GetText(resp + SENS_RESP_PTR_NAME));

                  loc1 = NewStat(tmpstr, 2, np, maxi);
                  WDB[arr + bin] = (double)loc1;
                }

              /* Allocate direct statistic for responses (currently only for 1 bin) */

              if ((loc1 = (long)RDB[det + DET_PTR_RBINS]) > VALID_PTR)
                {
                  /* Allocate statistic array if needed */

                  if ((arr = (long)RDB[loc1 + DET_RBIN_PTR_SENS_DIRECT_STAT_ARRAY]) < VALID_PTR)
                    {
                      /* Allocate memory for the sensitivity statistic array */

                      arr = ReallocMem(DATA_ARRAY, tbins);
                      WDB[loc1 + DET_RBIN_PTR_SENS_DIRECT_STAT_ARRAY] = (double)arr;
                    }

                  /* This may have been allocated as a part of previous reaction rate ratio */

                  if ((long)RDB[arr + bin] < VALID_PTR)
                    {
                      sprintf(tmpstr, "ADJ_PERT_DET_%s_SENS_DIRECT", GetText(resp + SENS_RESP_PTR_NAME));

                      ptr = NewStat(tmpstr, 1, maxi);
                      WDB[arr + bin] = (double)ptr;
                    }

                  if (NextItem(loc1) > VALID_PTR)
                    Error(resp, "Multiple response bins in detector %s in reaction rate ratio sensitivity",
                          GetText(det + DET_PTR_NAME));
                }

              /* Write skip fraction for */

              WDB[det + DET_SKIP_FRAC] = 1.0 - P;
            }

          /* Allocate sensitivity statistic  */

          sprintf(tmpstr, "ADJ_PERT_%s_SENS", GetText(resp + SENS_RESP_PTR_NAME));

          ptr = NewStat(tmpstr, 2, np, maxi);

          /* Store pointer to statistic */

          WDB[resp + SENS_RESP_PTR_STAT] = (double)ptr;
        }

      /* Flag this response to have partials so that it won't be collected itself */

      WDB[loc0 + SENS_RESP_HAS_PARTIALS] = (double)YES;

      /* Get next response */

      loc0 = NextItem(loc0);
    }

  fprintf(outp, "OK.\n\n");
}
