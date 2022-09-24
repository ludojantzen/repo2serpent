/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : makeburnmatrix.c                               */
/*                                                                           */
/* Created:       2011/01/25 (JLe)                                           */
/* Last modified: 2018/11/08 (Ari)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Generates burnup matrix for material                         */
/*                                                                           */
/* Comments: - Printataan data fileen testausmielessä                        */
/*                                                                           */
/*           - Flux is truncated to 6 decimals to avoid round-off errors in  */
/*             reproducible MPI mode.                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MakeBurnMatrix:"

/* Väliaikainen funktio palamaluennon kuvien plottausta varten */

void Funktio(long, long, double, long);

/*****************************************************************************/

struct ccsMatrix *MakeBurnMatrix(long mat, long id)
{
  long n, lst, iso, nuc, tgt, ptr, rea, yld, i, j, sz, type, nnz, nsz;
  long *row, *col;
  double br, lambda, rr, *vec, f, flx;
  complex *val;
  struct ccsMatrix *A;

  /* Check burnup mode and burn flag */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) ||
      (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)))
    return NULL;

  /***************************************************************************/

  /***** Set indexes *********************************************************/

  /* Loop over composition list and set indexes */

  i = 0;

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Store index */

      StoreValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX, (double)mat,
                     (double)i, id);

      /* Update index */

      i++;

      /* Next nuclide */

      iso = NextItem(iso);
    }

  /* Avoid compiler warning */

  rr = -1.0;
  flx = -1.0;

  /* Composition size */

  sz = i;

  /* Allocate memory for reaction rate vector */

  vec = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

  /* Get number of non-zero elements */

  nsz = BurnMatrixSize(mat);

  /* Allocate memory for matrix */

  A = ccsMatrixNew(sz, sz, nsz);

  /* Ota pointterit talteen */

  val = A->values;
  row = A->rowind;
  col = A->colptr;

  /* todellinen nnz lasketaan matriisia muodostettaessa */

  nnz = 0;

  /* Asetetaan 1. sarakepointteri osoittamaan alkuun */

  col[0] = 0;

  /***************************************************************************/

  /***** Create matrix *******************************************************/

  /* Get flux */

  flx = Truncate(RDB[mat + MATERIAL_BURN_FLUX_SSA], 6);

  /* Check value */

  CheckValue(FUNCTION_NAME, "flx", "", flx, 0.0, INFTY);

  /* Loop over composition */

  i = 0;
  lst = (long)RDB[mat + MATERIAL_PTR_COMP];

  while ((iso = ListPtr(lst, i)) > VALID_PTR)
    {
      /* Reset reaction rate vector */

      memset(vec, 0.0, sz*sizeof(double));

      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Get decay constant */

      lambda = RDB[nuc + NUCLIDE_LAMBDA];

      /***********************************************************************/

      /***** Decay and transmutation reactions in list ***********************/

      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Get reaction type */

          type = (long)RDB[rea + REACTION_TYPE];

          /* Get branching ratio */

          br = RDB[rea + REACTION_BR];

          /* Tämä printtaa noiden reaktioiden makroskooppiset vaikutusalat */
          /*
          if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
            if ((RDB[iso + COMPOSITION_ADENS]*
                 TestValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, id)
                 > 0.0) && ((long)RDB[rea + REACTION_MT] < 900) &&
                (type == REACTION_TYPE_PARTIAL))
              printf("%ld %ld %E\n", (long)RDB[nuc + NUCLIDE_ZAI],
                     (long)RDB[rea + REACTION_MT],
                     RDB[iso + COMPOSITION_ADENS]*
                     TestValuePair(rea + REACTION_PTR_TRANSMUXS,
                                   (double)mat, id));
          */

          /* Check reaction type */

          if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
            {
              /***************************************************************/

              /***** Fission *************************************************/

              /* Check type and get reaction rate */

              if (type == REACTION_TYPE_DECAY)
                {
                  /* Spontaneous fission */

                  rr = lambda;
                }
              else
                {
                  /* Reset rate */

                  rr = 0.0;

                  /* Neutron-induced fission (pointteri on PRIVA-blokkiin) */

                  if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                    rr = TestValuePair(rea + REACTION_PTR_TRANSMUXS,
                                       (double)mat, id);
                  else
                    Die(FUNCTION_NAME, "Pointer error 1 (%s %ld)",
                        GetText(nuc + NUCLIDE_PTR_NAME),
                        (long)RDB[rea + REACTION_MT]);

                  /* Check and convert units */

                  if (rr > 0.0)
                    rr = flx*rr*BARN;
                  else
                    rr = 0.0;

                  /* Branching ratio should be one */

                  if (br != 1.0)
                    Die(FUNCTION_NAME, "br = %E", br);
                }

              /* Check rate (may be zero if fission is cut off by user) */

              if (rr == 0.0)
                {
                  /* Next reaction */

                  rea = NextItem(rea);

                  /* Cycle loop */

                  continue;
                }

              /* Add total rate to diagonal diagonal */

              vec[i] = vec[i] - br*rr;

              /* Get pointer to distribution */

              yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];
              CheckPointer(FUNCTION_NAME, "yld", DATA_ARRAY, yld);

              /* Loop over distribution */

              n = 0;
              while ((ptr = ListPtr(yld, n++)) > 0)
                {
                  /* Get pointer to target */

                  tgt = (long)RDB[ptr + FY_PTR_TGT];
                  CheckPointer(FUNCTION_NAME, "(tgt1)", DATA_ARRAY, tgt);

                  /* Get target nuclide index */

                  j = (long)TestValuePair(tgt + NUCLIDE_PTR_MATRIX_IDX,
                                          (double)mat, id);

                  /* Check */

                  if (j < 0)
                    Die(FUNCTION_NAME, "Nuclide %s not found in composition",
                        GetText(tgt + NUCLIDE_PTR_NAME));

                  /* Get yield */

                  f = RDB[ptr + FY_INDEPENDENT_FRAC];

                  /* Check value */

                  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

                  /* Add to rate */

                  vec[j] = vec[j] + br*rr*f;
                }

              /***************************************************************/
            }
          else if ((tgt = (long)RDB[rea + REACTION_PTR_TGT]) > VALID_PTR)
            {
              /***************************************************************/

              /***** Decay or transmutation **********************************/

              /* Check type */

              if ((type == REACTION_TYPE_PARTIAL) ||
                  (type == REACTION_TYPE_TRA_BRANCH))
                {
                  /* Transmutation reaction (pointteri on PRIVA-blokkiin */

                  if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                    rr = TestValuePair(rea + REACTION_PTR_TRANSMUXS,
                                       (double)mat, id);
                  else
                    Die(FUNCTION_NAME, "Pointer error 2 (%s %ld)",
                        GetText(nuc + NUCLIDE_PTR_NAME),
                        (long)RDB[rea + REACTION_MT]);

                  /* Check and convert units */

                  if (rr > 0.0)
                    rr = flx*rr*BARN;
                  else
                    rr = 0.0;
                }
              else if ((type == REACTION_TYPE_DECAY) ||
                       (type == REACTION_TYPE_DEC_BRANCH))
                {
                  /* Decay reaction */

                  rr = lambda;
                }
              else
                Die(FUNCTION_NAME, "Invalid type");

              /* Get target nuclide index */

              j = (long)TestValuePair(tgt + NUCLIDE_PTR_MATRIX_IDX,
                                      (double)mat, id);

              /* Check */

              if (j < 0)
                Die(FUNCTION_NAME, "Nuclide %s not found in composition",
                    GetText(tgt + NUCLIDE_PTR_NAME));

              /* Add to vector */

              vec[j] = vec[j] + br*rr;

              /* Add to diagonal */

              if ((type == REACTION_TYPE_PARTIAL) ||
                  ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR))
                vec[i] = vec[i] - rr;
              else if (type == REACTION_TYPE_DECAY)
                vec[i] = vec[i] - br*rr;

              /***************************************************************/
            }

          /* Next reaction */

          rea = NextItem(rea);
        }

      /***********************************************************************/

      /***** Total fission ***************************************************/

      /* Check if fission yield data is not defined */

      if ((long)RDB[nuc + NUCLIDE_PTR_NFY_DATA] < VALID_PTR)
        {
          /* Check pointer to total fission */

          if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
            {
              /* Get reaction rate (pointteri on PRIVA-blokkiin) */

              if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
                rr = TestValuePair(rea + REACTION_PTR_TRANSMUXS,
                                   (double)mat, id);
              else
                Die(FUNCTION_NAME, "Pointer error 3 (%s %ld)",
                    GetText(nuc + NUCLIDE_PTR_NAME),
                    (long)RDB[rea + REACTION_MT]);

              /* Check and convert units */

              if (rr > 0.0)
                rr = flx*rr*BARN;
              else
                rr = 0.0;

              /* Pointer to target */

              tgt = (long)RDB[rea + REACTION_PTR_TGT];
              CheckPointer(FUNCTION_NAME, "(tgt)", DATA_ARRAY, tgt);

              /* Check pointer */

              if (tgt != (long)RDB[DATA_PTR_NUCLIDE_LOST])
                Die(FUNCTION_NAME, "Pointer error");

              /* Get target nuclide index */

              j = (long)TestValuePair(tgt + NUCLIDE_PTR_MATRIX_IDX,
                                      (double)mat, id);

              /* Check */

              if (j < 0)
                Die(FUNCTION_NAME, "Nuclide %s not found in composition",
                    GetText(tgt + NUCLIDE_PTR_NAME));

              /* Add to vector */

              vec[j] = vec[j] + rr;

              /* Add to diagonal */

              vec[i] = vec[i] - rr;
            }
        }

      /***********************************************************************/

      /***** Store values in matrix ******************************************/

      /* Loop over reaction rate vector */

      for (j = 0; j < sz; j++)
        if (vec[j] != 0.0)
          {
            /* Check value */

            CheckValue(FUNCTION_NAME, "vec[j]", "", vec[j], -INFTY, INFTY);

            if ((i != j) && (vec[j] < 0.0))
              Die(FUNCTION_NAME, "Negative off-diagonal value A(%ld,%ld) = %E",
                  i, j, vec[j]);
            else if ((i == j) && (vec[j] > 0.0))
              Die(FUNCTION_NAME, "Positive diagonal value A(%ld,%ld) = %E",
                  i, j, vec[j]);

            /* Put value in matrix */

            val[nnz].re = vec[j];
            val[nnz].im = 0.0;

            /* Check value */

            CheckValue(FUNCTION_NAME, "val[nnz].re", "", val[nnz].re,
                       -INFTY, INFTY);

            /* Put row index */

            row[nnz] = j;

            /* Add counter */

            nnz++;

            /* Check value */

            if (nnz > nsz)
              Die(FUNCTION_NAME, "nnz > nsz");
          }

      /* seuraava sarake alkaa tästä indeksistä */

      col[i+1] = nnz;

      /* Next nuclide */

      i++;

      /***********************************************************************/
    }

  /***************************************************************************/

  /* Free temporary array */

  Mem(MEM_FREE, vec);

  /* Set matrix size */

  A->nnz = nnz;

  /* Palauta matriisi A */

  return A;
}

/*****************************************************************************/

void Funktio(long rea, long mat, double flx, long id)
{
  long mt, nuc, ZAI, step, i;
  char tmpstr[MAX_STR];
  FILE *fp;

  mt = (long)RDB[rea + REACTION_MT];
  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  ZAI = (long)RDB[nuc + NUCLIDE_ZAI];
  step = (long)RDB[DATA_BURN_STEP] + 1;

  sprintf(tmpstr, "%s_dataa.m", GetText(DATA_PTR_INPUT_FNAME));

  if (strcmp(GetText(mat + MATERIAL_PTR_NAME), "UO2Gdp2r10"))
    return;

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    i = 2;
  else
    i = 1;

  if (((ZAI == 541350) || (ZAI == 942390) || (ZAI == 641550)
       || (ZAI == 922380) || (ZAI == 531310)) && (mt == 102))
    {
      fp = fopen(tmpstr, "a");

      fprintf(fp, "m%sflx(%ld,%ld) = %1.5E;\n",
              GetText(mat + MATERIAL_PTR_NAME),
              step, i, flx);
      fprintf(fp, "m%s%ld_mt%ld(%ld,%ld) = %1.5E;\n",
              GetText(mat + MATERIAL_PTR_NAME), ZAI, mt, step, i,
              TestValuePair(rea + REACTION_PTR_TRANSMUXS, (double)mat, id));

      fclose(fp);
    }
}

/*****************************************************************************/
