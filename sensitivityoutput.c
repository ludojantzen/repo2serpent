/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sensitivityoutput.c                            */
/*                                                                           */
/* Created:       2017/04/06 (VVa)                                           */
/* Last modified: 2019/01/18 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Writes sensitivity calculation output to file.               */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SensitivityOutput:"

/*****************************************************************************/

void SensitivityOutput()
{
  long nmat, nzai, nrea, nene, nmu, ngen, nbin, ptr, nmg, n, i, j, maxi;
  long sens, res, mat, zai, gen, maxmt, loc0, resp, nblock, block, idx, N;
  double val;
  char outfile[MAX_STR], tmpstr[MAX_STR];
  FILE *fp;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Get pointer to sensitivity block or return */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Check corrector step */

  if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
      (RDB[DATA_BURN_SIE] == (double)NO))
    return;

  /* Check if in active cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Reduce private results array */

  ReducePrivateRes();

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

  /* Check branch index */

  if ((long)RDB[DATA_COEF_CALC_IDX] < 0)
    {
      if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
        sprintf(outfile, "%s_sens%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                (long)RDB[DATA_DYN_TB]);
      else
        sprintf(outfile, "%s_sens%ld.m", GetText(DATA_PTR_INPUT_FNAME),
                (long)RDB[DATA_BURN_STEP]);
    }
  else
    sprintf(outfile, "%s_sens%ldb%ld.m", GetText(DATA_PTR_INPUT_FNAME),
            (long)RDB[DATA_COEF_CALC_BU_IDX],
            (long)RDB[DATA_COEF_CALC_IDX]);

  /* Open file for writing */

  fp = fopen(outfile, "w");

  /******************************/
  /**** Move on to printing *****/
  /******************************/

  /* Get sizes of sub bins */

  nmat = (long)RDB[sens + SENS_N_MAT];
  nzai = (long)RDB[sens + SENS_N_ZAI];
  nrea = (long)RDB[sens + SENS_N_PERT] + 1;
  nene = (long)RDB[sens + SENS_N_ENE] + 1;
  nmu = (long)RDB[sens + SENS_N_MU];
  nblock = (long)RDB[sens + SENS_N_COV_BLOCK];

  fprintf(fp, "%% Number of different bins in sensitivity calculation:\n\n");

  fprintf(fp, "SENS_N_MAT = %ld;\n", nmat - 1);
  fprintf(fp, "SENS_N_ZAI = %ld;\n", nzai - 1);
  fprintf(fp, "SENS_N_PERT = %ld;\n", nrea);
  fprintf(fp, "SENS_N_ENE = %ld;\n", nene - 1);
  fprintf(fp, "SENS_N_MU = %ld;\n", nmu);

  fprintf(fp, "\n");

  /* Calculate maximum bin index */

  nbin = (1 + (nmat*nzai*nrea*nene*nmu));

  /* Get maximum generation index */

  ngen = (long)RDB[DATA_SENS_LAST_GEN];

  /***********************/
  /* Print material list */
  /***********************/

  fprintf(fp, "%% Materials included in sensitivity calculation:\n\n");

  fprintf(fp, "SENS_MAT_LIST = [\n");

  /* Include total-material if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
    fprintf(fp, "\'%-30s\'\n", "total");

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_SUM)
    fprintf(fp, "\'%-30s\'\n", "sum");

  if ((ptr = (long)RDB[sens + SENS_PTR_MAT_ARR]) > VALID_PTR)
    {
      /* Loop over the material list to print material names */

      i = 0;

      while ((mat = (long)RDB[ptr + i]) > VALID_PTR)
        {

          fprintf(fp, "\'%-30s\'\n", GetText(mat + MATERIAL_PTR_NAME));

          /* Increment index */

          i++;
        }

    }

  fprintf(fp, "];\n\n");

  fprintf(fp, "%% Indices for different materials:\n\n");

  /* Reset index */

  n = 1;

  /* Include total-material if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
      fprintf(fp, "iSENS_MAT_%-16s = %ld;\n", "TOT", n++);

  /* Include sum material if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_SUM)
    fprintf(fp, "iSENS_MAT_%-16s = %ld;\n", "SUM", n++);

  if ((ptr = (long)RDB[sens + SENS_PTR_MAT_ARR]) > VALID_PTR)
    {
      /* Loop over the material list to print material names */

      i = 0;

      while ((mat = (long)RDB[ptr + i]) > VALID_PTR)
        {

          fprintf(fp, "iSENS_MAT_%-16s = %ld;\n",
                  GetText(mat + MATERIAL_PTR_NAME), n++);

          /* Increment index */

          i++;
        }

    }

  fprintf(fp, "\n");

  /**********************/
  /* Print nuclide list */
  /**********************/

  fprintf(fp, "%% Nuclides included in sensitivity calculation:\n\n");

  fprintf(fp, "SENS_ZAI_LIST = [\n");

  /* Include total-ZAI if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_TOT)
    fprintf(fp, "%-8d %% total\n", 0);

  /* Include sum-ZAI if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_SUM)
    fprintf(fp, "%-8d %% sum\n", 1);

  if ((ptr = (long)RDB[sens + SENS_PTR_ZAI_ARR]) > VALID_PTR)
    {
      /* Loop over the material list to print material names */

      i = 0;

      while ((zai = (long)RDB[ptr + i]) >= 10010)
        {

          fprintf(fp, "%-8ld\n", zai);

          /* Increment index */

          i++;
        }
    }

  fprintf(fp, "];\n\n");

  fprintf(fp, "%% Indices for different ZAIs:\n\n");

  n = 1;

  /* Include total-ZAI if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_TOT)
    fprintf(fp, "iSENS_ZAI_%-8s = %ld;\n", "TOT", n++);

  /* Include sum-ZAI if requested */

  if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_SUM)
    fprintf(fp, "iSENS_ZAI_%-8s = %ld;\n", "SUM", n++);

  if ((ptr = (long)RDB[sens + SENS_PTR_ZAI_ARR]) > VALID_PTR)
    {
      /* Loop over the material list to print material names */

      i = 0;

      while ((zai = (long)RDB[ptr + i]) >= 10010)
        {

          fprintf(fp, "iSENS_ZAI_%-8ld = %ld;\n", zai, n++);

          /* Increment index */

          i++;
        }
    }

  fprintf(fp, "\n");

  /***********************/
  /* Print reaction list */
  /***********************/

  fprintf(fp, "%% Reactions included in sensitivity calculation:\n\n");

  fprintf(fp, "SENS_PERT_LIST = [\n");
  fprintf(fp, "\'%-30s\'\n", "total xs");

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS)
    {
      /* Get pointer to reaction indices */

      ptr = (long)RDB[sens + SENS_PTR_PERT_INDICES];

      /* Check if reaction was set on for each of the reactions */

      if (RDB[ptr + ELA_SCATT_IDX] > 0)
        fprintf(fp, "\'%-30s\'\n", "ela scatt xs");
      if (RDB[ptr + SAB_SCATT_IDX] > 0)
        fprintf(fp, "\'%-30s\'\n", "sab scatt xs");
      if (RDB[ptr + INL_SCATT_IDX] > 0)
        fprintf(fp, "\'%-30s\'\n", "inl scatt xs");
      if (RDB[ptr + CAPT_IDX] > 0)
        fprintf(fp, "\'%-30s\'\n", "capture xs");
      if (RDB[ptr + FISS_IDX] > 0)
        fprintf(fp, "\'%-30s\'\n", "fission xs");
      if (RDB[ptr + NXN_SCATT_IDX] > 0)
        fprintf(fp, "\'%-30s\'\n", "nxn xs");
    }
  else if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
    {
      /* Get maximum perturbed MT */

      maxmt = (long)RDB[sens + SENS_MAX_MT];

      /* Get pointer to indices */

      ptr = (long)RDB[sens + SENS_PTR_MT_INDICES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Do not include nubar or chi to MTs here */

      if (maxmt >= 452)
        maxmt = 451;

      /* Loop over all possible MT numbers */

      for (i = 0; i < maxmt + 1; i++)
        {
          /* Print out if included in perturbations */

          if ((long)RDB[ptr + i] > 0)
            {
              sprintf(tmpstr, "mt %ld xs", i);
              fprintf(fp, "\'%-30s\'\n", tmpstr);
            }
        }
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE)
    {
      fprintf(fp, "\'%-30s\'\n", "temperature");
      fprintf(fp, "\'%-30s\'\n", "free gas temperature");
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
    {
      fprintf(fp, "\'%-30s\'\n", "nubar total");
      fprintf(fp, "\'%-30s\'\n", "nubar prompt");
      fprintf(fp, "\'%-30s\'\n", "nubar delayed");
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
    {
      fprintf(fp, "\'%-30s\'\n", "chi total");
      fprintf(fp, "\'%-30s\'\n", "chi prompt");
      fprintf(fp, "\'%-30s\'\n", "chi delayed");
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_ELA_MU)
    {
      fprintf(fp, "\'%-30s\'\n", "elastic scattering mu (CM)");
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_INL_MU)
    {
      fprintf(fp, "\'%-30s\'\n", "inelastic scattering cosine (CM)");
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM)
    {
      /* Get maximum scattering moment */

      maxi = (long)RDB[sens + SENS_MAX_SCATT_MOM];

      for (i = 0; i < maxi; i++)
            {
              sprintf(tmpstr, "%s %ld", "ela leg mom", i+1);
              fprintf(fp, "\'%-30s\'\n", tmpstr);
            }
    }

  /* Custom perturbations */

  /* Get pointer to first custom perturbation */

  loc0 = (long)RDB[sens + SENS_PTR_PERT0];

  /* Loop over custom perturbations */

  while (loc0 > VALID_PTR)
    {
      fprintf(fp, "\'%-30s\'\n", GetText(loc0 + SENS_PERT_PTR_NAME));

      /* Next item */

      loc0 = NextItem(loc0);
    }

  fprintf(fp, "];\n\n");

  /******************************/
  /* Print perturbation indices */
  /******************************/

  fprintf(fp, "%% Indices for different perturbations:\n\n");

  ptr =   (long)RDB[sens + SENS_PTR_PERT_INDICES];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_TOT_XS",         1 + (long)RDB[ptr + TOT_REA_IDX]);

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS)
    {
      /* Check if reaction was set on for each of the reactions */

      if (RDB[ptr + ELA_SCATT_IDX] > 0)
        fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_ELA_XS",    1 + (long)RDB[ptr + ELA_SCATT_IDX]);
      if (RDB[ptr + SAB_SCATT_IDX] > 0)
        fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_SAB_XS",    1 + (long)RDB[ptr + SAB_SCATT_IDX]);
      if (RDB[ptr + INL_SCATT_IDX] > 0)
        fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_INL_XS",    1 + (long)RDB[ptr + INL_SCATT_IDX]);
      if (RDB[ptr + CAPT_IDX] > 0)
        fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_CAPT_XS",         1 + (long)RDB[ptr + CAPT_IDX]);
      if (RDB[ptr + FISS_IDX] > 0)
        fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_FISS_XS",         1 + (long)RDB[ptr + FISS_IDX]);
      if (RDB[ptr + NXN_SCATT_IDX] > 0)
        fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_NXN_XS",          1 + (long)RDB[ptr + NXN_SCATT_IDX]);
    }
  else if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
    {
      /* Get maximum perturbed MT */

      maxmt = (long)RDB[sens + SENS_MAX_MT];

      /* Get pointer to indices */

      ptr = (long)RDB[sens + SENS_PTR_MT_INDICES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over all possible MT numbers */

      for (i = 0; i < maxmt + 1; i++)
        {
          /* Print out if included in perturbations */

          if ((long)RDB[ptr + i] > 0)
            {
              sprintf(tmpstr, "SENS_PERT_MT_%ld_XS", i);
              fprintf(fp, "i%-22s = %ld;\n", tmpstr,          1 + (long)RDB[ptr + i]);
            }
        }
    }

  ptr =  (long)RDB[sens + SENS_PTR_PERT_INDICES];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE)
    {
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_TEMPERATURE",    1 + (long)RDB[ptr + TEMPERATURE_IDX]);
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_TEMPERATURE",    1 + (long)RDB[ptr + SCATT_TEMP_IDX]);
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR)
    {
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_NUBAR_TOT",    1 + (long)RDB[ptr + NUBAR_TOT_IDX]);
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_NUBAR_PROMPT", 1 + (long)RDB[ptr + NUBAR_PRO_IDX]);
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_NUBAR_DEL",    1 + (long)RDB[ptr + NUBAR_DEL_IDX]);
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_CHI)
    {
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_CHI_TOT",      1 + (long)RDB[ptr + CHI_TOT_IDX]);
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_CHI_PROMPT",   1 + (long)RDB[ptr + CHI_PRO_IDX]);
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_CHI_DEL",      1 + (long)RDB[ptr + CHI_DEL_IDX]);
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_ELA_MU)
    {
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_ELA_MU",      1 + (long)RDB[ptr + ELA_MU_IDX]);
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_INL_MU)
    {
      fprintf(fp, "i%-22s = %ld;\n", "SENS_PERT_INL_MU",      1 + (long)RDB[ptr + INL_MU_IDX]);
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM)
    {
      /* Get maximum scattering moment */

      maxi = (long)RDB[sens + SENS_MAX_SCATT_MOM];

      for (i = 0; i < maxi; i++)
        {
          sprintf(tmpstr, "%s%ld", "SENS_PERT_ELA_P", i+1);
          fprintf(fp, "i%-22s = %ld;\n", tmpstr,       1 + (long)RDB[ptr + ELA_P1_IDX + i]);
        }
    }

  /* Custom perturbations */

  /* Get pointer to first custom perturbation */

  loc0 = (long)RDB[sens + SENS_PTR_PERT0];

  /* Loop over custom perturbations */

  while (loc0 > VALID_PTR)
    {
      fprintf(fp, "iSENS_PERT_%-12s = %ld;\n", GetText(loc0 + SENS_PERT_PTR_NAME),
              1 + (long)RDB[loc0 + SENS_PERT_INDEX]);

      /* Next item */

      loc0 = NextItem(loc0);
    }

  fprintf(fp, "\n");

  /**************************************/
  /* Print covariance block information */
  /**************************************/

  if (nblock > 0)
    {
      fprintf(fp, "%% Information on covariance blocks\n\n");

      fprintf(fp, "COV_N_BLOCK = %ld;\n\n", nblock);

      /* Get pointer to first block */

      block = (long)RDB[DATA_PTR_COVBLOCK0];
      CheckPointer(FUNCTION_NAME, "(block)", DATA_ARRAY, block);

      idx = 1;
      while (block > VALID_PTR)
        {
          /* Get block order (number of ZAIMTs) */

          N = (long)RDB[block + COVBLOCK_ORDER];

          /* Pointer to the ZAIMT list for the block */

          ptr = (long)RDB[block + COVBLOCK_PTR_ZAIMT_ARRAY];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Print the ZAIMT list for the block */

          fprintf(fp, "COV_BLOCK_%ld_ZAIMTS = [", idx);

          for (i = 0; i < N; i++)
            fprintf(fp, " %ld ", (long)RDB[ptr + i]);

          fprintf(fp, "];\n\n");

          /* Pointer to the ZAIMT list for the block */

          ptr = (long)RDB[block + COVBLOCK_PTR_COVMTX_ARRAY];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Print the found data matrix */

          fprintf(fp, "COV_BLOCK_%ld_CONTAINS_DATA = [\n", idx);

          for (i = 0; i < N; i++)
            {
              for (j = 0; j < N; j++)
                {
                  if ((long)RDB[ptr + i*N + j] > VALID_PTR)
                    fprintf(fp, " 1 ");
                  else
                    fprintf(fp, " 0 ");
                }
              fprintf(fp, "\n");
            }

          fprintf(fp, "];\n\n");

          /* Next block */

          block = NextItem(block);
          idx++;
        }
    }

  /*********************/
  /* Print energy grid */
  /*********************/

  /* Pointer to group structure */

  ptr = (long)RDB[sens + SENS_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of group boundaries */

  nmg = (long)RDB[sens + SENS_N_ENE] + 1;

  /* Pointer to data */

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  fprintf(fp, "%% Sensitivity calculation energy group boundaries:\n\n");

  fprintf(fp, "SENS_E = [ ");

  for (n = 0; n < nmg; n++)
    fprintf(fp, "%12.5E ", RDB[ptr + n]);

  fprintf(fp, "];\n");

  fprintf(fp, "\n");

  /*************************************/
  /* Print lethargy widths of each bin */
  /*************************************/

  /* Pointer to group structure */

  ptr = (long)RDB[sens + SENS_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of group boundaries */

  nmg = (long)RDB[sens + SENS_N_ENE] + 1;

  /* Pointer to data */

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  fprintf(fp, "%% Sensitivity calculation energy group lethargy widths:\n\n");

  fprintf(fp, "SENS_LETHARGY_WIDTHS = [ ");

  for (n = 0; n < nmg - 1; n++)
    {
      val = log(RDB[ptr + nmg - 1]/RDB[ptr + n]) - log(RDB[ptr + nmg - 1]/RDB[ptr + n + 1]);
      fprintf(fp, "%12.5E ", val);
    }

  fprintf(fp, "];\n");

  fprintf(fp, "\n");

  /***********************************************************/
  /* Print results with maximum number of latent generations */
  /***********************************************************/

  gen = ngen - 1;

  fprintf(fp, "%% Sensitivities with %ld latent generations:\n\n", gen);

  resp = (long)RDB[sens + SENS_PTR_RESP0];

  while (resp > VALID_PTR)
    {
      /* Do not output responses that have partials */

      if ((long)RDB[resp + SENS_RESP_HAS_PARTIALS] == YES)
        {
          /* Process next response */

          resp = NextItem(resp);
          continue;
        }

      if ((long)RDB[resp + SENS_RESP_TYPE] == SENS_RESP_TYPE_KEFF)
        {
          fprintf(fp, "%% Effective multiplication factor:\n\n");
        }
      else if ((long)RDB[resp + SENS_RESP_TYPE] == SENS_RESP_TYPE_BEFF)
        {
          fprintf(fp, "%% Effective delayed neutron fraction:\n\n");
        }
      else if ((long)RDB[resp + SENS_RESP_TYPE] == SENS_RESP_TYPE_LEFF)
        {
          fprintf(fp, "%% Effective prompt neutron generation time:\n\n");
        }
      else if ((long)RDB[resp + SENS_RESP_TYPE] == SENS_RESP_TYPE_VOID)
        {
          mat = (long)RDB[sens + SENS_PTR_VOID_MAT];

          fprintf(fp, "%% Void coefficient of material %s:\n\n", GetText(mat + MATERIAL_PTR_NAME));
        }
      else if ((long)RDB[resp + SENS_RESP_TYPE] == SENS_RESP_TYPE_RATIO)
        {
          fprintf(fp, "%% Reaction rate ratio:\n\n");
        }

      res = (long)RDB[resp + SENS_RESP_PTR_STAT];
      CheckPointer(FUNCTION_NAME, "(res)", DATA_ARRAY, res);

      WriteSensResult(fp, NULL, res, nmat, nzai, nrea, nene, nmu, gen);

      fprintf(fp, "\n");

      if ((ptr = (long)RDB[resp + SENS_RESP_PTR_UNC_STAT]) > VALID_PTR)
        {
          /* Write uncertainty results */

          fprintf(fp, "%-26s = [", GetText(ptr + SCORE_PTR_NAME));

          for (n = 0; n < nblock+1; n++)
            fprintf(fp, " %12.5E %7.1E ", Mean(ptr, ngen-1, n),
                    RelErr(ptr, ngen-1, n));

          fprintf(fp, "];\n\n");
        }

      if ((ptr = (long)RDB[resp + SENS_RESP_PTR_VAR_STAT]) > VALID_PTR)
        {
          /* Write uncertainty results */

          fprintf(fp, "%-26s = [", GetText(ptr + SCORE_PTR_NAME));

          for (n = 0; n < nblock+1; n++)
            fprintf(fp, " %12.5E %7.1E ", Mean(ptr, ngen-1, n),
                    RelErr(ptr, ngen-1, n));

          fprintf(fp, "];\n\n");
        }

      resp = NextItem(resp);
    }

  /*************************************************************/
  /* Print results with a smaller number of latent generations */
  /*************************************************************/

  if ((long)RDB[sens + SENS_RESP_FLAGS] & SENS_SCORE_FLAG_HIS)
    {

      /* Loop over the lower generations */

      for (gen  = 0; gen < ngen - 1; gen++)
        {

          fprintf(fp, "%% Sensitivities with %ld latent generations:\n\n", gen);

          if ((long)RDB[sens + SENS_RESP_FLAGS] & SENS_RESP_FLAG_KEFF)
            {
              fprintf(fp, "%% Effective multiplication factor:\n\n");

              res = (long)RDB[RES_ADJ_PERT_KEFF_SENS];
              CheckPointer(FUNCTION_NAME, "(res)", DATA_ARRAY, res);

              sprintf(tmpstr, "%s_%ld_GEN", GetText(res + SCORE_PTR_NAME), gen);

              WriteSensResult(fp, tmpstr, res, nmat, nzai, nrea, nene, nmu, gen);
            }

          if ((long)RDB[sens + SENS_RESP_FLAGS] & SENS_RESP_FLAG_BEFF)
            {
              fprintf(fp, "%% Effective delayed neutron fraction:\n\n");

              res = (long)RDB[RES_ADJ_PERT_BEFF_SENS];
              CheckPointer(FUNCTION_NAME, "(res)", DATA_ARRAY, res);

              sprintf(tmpstr, "%s_%ld_GEN", GetText(res + SCORE_PTR_NAME), gen);

              WriteSensResult(fp, tmpstr, res, nmat, nzai, nrea, nene, nmu, gen);
            }

          if ((long)RDB[sens + SENS_RESP_FLAGS] & SENS_RESP_FLAG_LEFF)
            {
              fprintf(fp, "%% Effective prompt neutron generation time:\n\n");

              res = (long)RDB[RES_ADJ_PERT_LEFF_SENS];
              CheckPointer(FUNCTION_NAME, "(res)", DATA_ARRAY, res);

              sprintf(tmpstr, "%s_%ld_GEN", GetText(res + SCORE_PTR_NAME), gen);

              WriteSensResult(fp, tmpstr, res, nmat, nzai, nrea, nene, nmu, gen);
            }

          if ((long)RDB[sens + SENS_RESP_FLAGS] & SENS_RESP_FLAG_VOID)
            {
              mat = (long)RDB[sens + SENS_PTR_VOID_MAT];

              fprintf(fp, "%% Void coefficient of material %s:\n\n", GetText(mat + MATERIAL_PTR_NAME));

              res = (long)RDB[RES_ADJ_PERT_VOID_SENS];
              CheckPointer(FUNCTION_NAME, "(res)", DATA_ARRAY, res);

              sprintf(tmpstr, "%s_%ld_GEN", GetText(res + SCORE_PTR_NAME), gen);

              WriteSensResult(fp, tmpstr, res, nmat, nzai, nrea, nene, nmu, gen);
            }

          resp = (long)RDB[sens + SENS_PTR_RESP0];

          while (resp > VALID_PTR)
            {
              /* Skip non-ratio responses here */

              if ((long)RDB[resp + SENS_RESP_TYPE] != SENS_RESP_TYPE_RATIO)
                {
                  resp = NextItem(resp);
                  continue;
                }

              /* Do not output responses that have partials */

              if ((long)RDB[resp + SENS_RESP_HAS_PARTIALS] == YES)
                {
                  /* Process next response */

                  resp = NextItem(resp);
                  continue;
                }

              fprintf(fp, "%% Reaction rate ratio %s:\n\n", GetText(resp + SENS_RESP_PTR_NAME));

              res = (long)RDB[resp + SENS_RESP_PTR_STAT];
              CheckPointer(FUNCTION_NAME, "(res)", DATA_ARRAY, res);

              sprintf(tmpstr, "%s_%ld_GEN", GetText(res + SCORE_PTR_NAME), gen);

              WriteSensResult(fp, tmpstr, res, nmat, nzai, nrea, nene, nmu, gen);

              resp = NextItem(resp);
            }

          fprintf(fp, "\n");
        }
    }

  fclose(fp);
}
