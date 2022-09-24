/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : makeburnmatrix.c                               */
/*                                                                           */
/* Created:       2015/03/28 (JLe)                                           */
/* Last modified: 2017/04/03 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Generates burnup matrix for material                         */
/*                                                                           */
/* Comments: - Separate subroutine for testing continuous reprocessing       */
/*             for MSR's                                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MakeBurnMatrixMSR:"

/*****************************************************************************/

struct ccsMatrix *MakeBurnMatrixMSR(long mat0, long dep, double **R, long id)
{
  long n, lst, iso, nuc, tgt, ptr, rea, yld, i, j, sz, type, nnz, nsz, im;
  long *row, *col, nm, mat, idx1, idx2, ii, mat1, loc0, iso1, nuc1, im1;
  long mat2, im2, nuc2, iso2, ZAI, loc1, rep;
  double br, lambda, rr, *vec, f, flx, tot;
  complex *val; 
  struct ccsMatrix *A; 

  /* Check burnup mode and burn flag */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) ||
      (!((long)RDB[mat0 + MATERIAL_OPTIONS] & OPT_BURN_MAT)))
    return NULL;

  /* Get pointer to reprocessor */

  rep = (long)RDB[dep + DEP_HIS_PTR_REPROC];
  CheckPointer(FUNCTION_NAME, "(rep)", DATA_ARRAY, rep);
  
  /***************************************************************************/
  
  /***** Set indexes *********************************************************/

  /* Check that material is first in chain */

  if ((long)RDB[mat0 + MATERIAL_FLOW_IDX] != 1)
    Die(FUNCTION_NAME, "Error in index");
 
  /* Get number of materials in chain */

  nm = (long)RDB[mat0 + MATERIAL_FLOW_N];
  CheckValue(FUNCTION_NAME, "nm", "", nm, 1, 20);

  /* Loop over composition list and set indexes */

  i = 0;

  iso = (long)RDB[mat0 + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      
      /* Store index */

      StoreValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX, (double)mat0, 
                     (double)i, id);

      /* Update index */

      i = i + nm;

      /* Next nuclide */

      iso = NextItem(iso);
    }

  /* Avoid compiler warning */

  rr = -1.0;
  flx = -1.0;

  /* Composition size */

  sz = i;

  /* Allocate memory for reaction rate vector */

  vec = (double *)Mem(MEM_ALLOC, sz + 1, sizeof(double));
  
  /* Get number of non-zero elements (NOTE: tää!!!!) */

  nsz = BurnMatrixSize(mat0) + sz;
  nsz = nsz*nm;

  /* Allocate memory for matrix */

  A = ccsMatrixNew(sz + 1, sz + 1, nsz);

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

  /* Reset index */

  ii = 0;

  /* Loop over columns (rows?) */

  for (i = 0; i < sz; i++)
    {
      /* Reset reaction rate vector */

      memset(vec, 0.0, (sz + 1)*sizeof(double));

      /* Find corresponding material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check index and that material belongs to same chain */

          if ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] == mat0)
            if ((i % nm) == (long)RDB[mat + MATERIAL_FLOW_IDX] - 1)
              break;
                 
          /* Pointer to next */

          mat = NextItem(mat);
        }

      /* Check pointer */
      
      if (mat < VALID_PTR)
        Die(FUNCTION_NAME, "Pointer error");

      /* Get material index */

      im = (long)RDB[mat + MATERIAL_FLOW_IDX] - 1;
      CheckValue(FUNCTION_NAME, "im", "", im, 0, nm - 1);

      /* Get flux */

      flx = Truncate(RDB[mat + MATERIAL_BURN_FLUX_SSA], 6);
      CheckValue(FUNCTION_NAME, "flx", "", flx, 0.0, INFTY);

      /* Get pointer to composition list */
      
      lst = (long)RDB[mat + MATERIAL_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

      /* Nuclide index */

      idx1 = (long)(((double)i)/((double)nm));

      /* Get pointer to composition */

      iso = ListPtr(lst, idx1);
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

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
                  /* Check material pointer */

                  if ((long)R[0][ii] != mat)
                    Die(FUNCTION_NAME, "Mismatch in material pointer");

                  /* Check nuclide ZAI */

                  if ((long)RDB[(long)R[1][ii] + NUCLIDE_ZAI] != 
                      (long)RDB[nuc + NUCLIDE_ZAI])
                    Die(FUNCTION_NAME, "Mismatch in target nuclide");

                  /* Check reaction pointer */

                  if (((long)R[1][ii] == nuc) && ((long)R[2][ii] != rea))
                    Die(FUNCTION_NAME, "Mismatch in reaction pointer");

                  /* Get reaction rate */

                  rr = R[3][ii++];

                  /* Check and convert units */
                  
                  if (rr > 0.0)
                    rr = flx*rr*BARN;
                  else
                    rr = 0.0;

                  /* Branching ratio should be one */

                  if (br != 1.0)
                    Die(FUNCTION_NAME, "br = %E", br);
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
                  CheckPointer(FUNCTION_NAME, "(tgt)", DATA_ARRAY, tgt);

                  /* Get target nuclide index */

                  idx2 = (long)TestValuePair(tgt + NUCLIDE_PTR_MATRIX_IDX, 
                                             (double)mat0, id);
                  
                  /* Convert */

                  j = idx2 + im;

                  /* Check */
                  
                  if (j < 0)
                    Die(FUNCTION_NAME, "Nuclide %s not found in composition",
                        GetText(tgt + NUCLIDE_PTR_NAME));
                  else if (j > sz - 1)
                    Die(FUNCTION_NAME, "Indexing error");

                  /* Get yield */
                  
                  f = RDB[ptr + FY_INDEPENDENT_FRAC];
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
                  /* Check material pointer */

                  if ((long)R[0][ii] != mat)
                    Die(FUNCTION_NAME, "Mismatch in material pointer");

                  /* Check nuclide ZAI */

                  if ((long)RDB[(long)R[1][ii] + NUCLIDE_ZAI] != 
                      (long)RDB[nuc + NUCLIDE_ZAI])
                    Die(FUNCTION_NAME, "Mismatch in target nuclide");

                  /* Check reaction pointer */

                  if (((long)R[1][ii] == nuc) && ((long)R[2][ii] != rea))
                    Die(FUNCTION_NAME, "Mismatch in reaction pointer");

                  /* Get reaction rate */

                  rr = R[3][ii++];

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
              
              idx2 = (long)TestValuePair(tgt + NUCLIDE_PTR_MATRIX_IDX, 
                                         (double)mat0, id);
              
              /* Convert */

              j = idx2 + im;

              /* Check */

              if (j < 0)
                Die(FUNCTION_NAME, "Nuclide %s not found in composition",
                    GetText(tgt + NUCLIDE_PTR_NAME));
              else if (j > sz - 1)
                Die(FUNCTION_NAME, "Indexing error");

              /* Add to vector */

              vec[j] = vec[j] + br*rr;

              /* Add to diagonal */

              if (type == REACTION_TYPE_PARTIAL)
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
              /* Check material pointer */

              if ((long)R[0][ii] != mat)
                Die(FUNCTION_NAME, "Mismatch in material pointer");

              /* Check nuclide ZAI */
              
              if ((long)RDB[(long)R[1][ii] + NUCLIDE_ZAI] != 
                  (long)RDB[nuc + NUCLIDE_ZAI])
                Die(FUNCTION_NAME, "Mismatch in target nuclide");
              
              /* Check reaction pointer */
              
              if (((long)R[1][ii] == nuc) && ((long)R[2][ii] != rea))
                Die(FUNCTION_NAME, "Mismatch in reaction pointer");

              /* Get reaction rate */

              rr = R[3][ii++];
              
              /* Check and convert units */
              
              if (rr > 0.0)
                rr = flx*rr*BARN;
              else
                rr = 0.0;
              
              /* Pointer to target */
              
              tgt = (long)RDB[rea + REACTION_PTR_TGT];
              
              /* Check pointer */

              if (tgt != (long)RDB[DATA_PTR_NUCLIDE_LOST])
                Die(FUNCTION_NAME, "Pointer error");
              
              /* Get target nuclide index */
              
              idx2 = (long)TestValuePair(tgt + NUCLIDE_PTR_MATRIX_IDX, 
                                         (double)mat0, id);

              /* Convert */

              j = idx2 + im;

              /* Check */
              
              if (j < 0)
                Die(FUNCTION_NAME, "Nuclide %s not found in composition",
                    GetText(tgt + NUCLIDE_PTR_NAME));
              else if (j > sz - 1)
                Die(FUNCTION_NAME, "Indexing error");

              /* Add to vector */

              vec[j] = vec[j] + rr;
              
              /* Add to diagonal */
              
              vec[i] = vec[i] - rr;
            }
        }

      /***********************************************************************/

      /***** Reprocessing mass flows *****************************************/
    
      /* Loop over outflow */

      loc0 = (long)RDB[mat + MATERIAL_PTR_OUTFLOW];
      while (loc0 > VALID_PTR)
        {
          /* Check pointer */

          if ((long)RDB[loc0 + REPROC_CON_PTR_MAT1] != mat)
            Die(FUNCTION_NAME, "Pointer error");

          /* Check mode */

          if ((long)RDB[loc0 + REPROC_CON_MODE] > 1)
            {
              /* Pointer to next */
              
              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Check that pointer to reprocessor matches that given for */
          /* depletion interval */

          if (rep != (long)RDB[loc0 + REPROC_CON_PTR_REP])
            {
              /* Pointer to next */
              
              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Loop over list */

          ptr = (long)RDB[loc0 + REPROC_CON_PTR_MFLOW];
          while (ptr > VALID_PTR)
            {
              /* Get removal rate */
              
              rr = RDB[ptr + MFLOW_LIST_RATE];
              CheckValue(FUNCTION_NAME, "rr", "", rr, 0.0, INFTY);

              /* Pointer to composition */

              iso = (long)RDB[ptr + MFLOW_LIST_PTR_ISO0];
              CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
              
              /* Check nuclide pointer */
              
              if ((long)RDB[iso + COMPOSITION_PTR_NUCLIDE] == nuc)
                {
                  /* Get target material and nuclide pointers */

                  mat1 = (long)RDB[loc0 + REPROC_CON_PTR_MAT2];
                  CheckPointer(FUNCTION_NAME, "(mat1)", DATA_ARRAY, mat1);

                  iso1 = (long)RDB[ptr + MFLOW_LIST_PTR_ISO1];
                  CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);

                  nuc1 = (long)RDB[iso1 + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY, nuc1);

                  /* Check that ZAI's match */

                  if ((long)RDB[nuc1 + NUCLIDE_ZAI] != 
                      (long)RDB[nuc + NUCLIDE_ZAI])
                    Die(FUNCTION_NAME, "Mismatch in ZAI");
                    
                  /* Get target material index */
                  
                  im1 = (long)RDB[mat1 + MATERIAL_FLOW_IDX] - 1;
                  CheckValue(FUNCTION_NAME, "im1", "", im1, 0, nm - 1);

                  /* Get target nuclide index */
                  
                  idx2 = (long)TestValuePair(nuc1 + NUCLIDE_PTR_MATRIX_IDX, 
                                             (double)mat0, id);
                  
                  /* Convert */
                  
                  j = idx2 + im1;

                  /* Check indexes */

                  if (i == j)
                    Die(FUNCTION_NAME, "Shouldn't happen");

                  /* Add to vector */

                  vec[j] = vec[j] + rr;
              
                  /* Add to diagonal */
                  
                  if ((long)RDB[loc0 + REPROC_CON_MODE] == 1)
                    vec[i] = vec[i] - rr;                  
                } 
              
              /* Next in list */

              ptr = NextItem(ptr);
            }
          
          /* Next */

          loc0 = NextItem(loc0);
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
            
            if (nnz > nsz - 1)
              Die(FUNCTION_NAME, "Index error (1)");
            else
              row[nnz] = j; 

            /* Add counter */

            nnz++;
          
            /* Check value */

            if (nnz > nsz)
              Die(FUNCTION_NAME, "nnz > nsz");
          }

      /* seuraava sarake alkaa tästä indeksistä */
      
      if (i + 1 > sz)
        Die(FUNCTION_NAME, "Index error (2)");
      else
        col[i + 1] = nnz; 

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Constant flows ******************************************************/

  /* Reset reaction rate vector */

  memset(vec, 0.0, (sz + 1)*sizeof(double));

  /* Loop over materials */

  mat1 = (long)RDB[DATA_PTR_M0];
  while (mat1 > VALID_PTR)
    {
      /* Check index and that material belongs to same chain */
      
      if ((long)RDB[mat1 + MATERIAL_FLOW_PTR_FIRST] != mat0)
        {
          /* Pointer to next */

          mat1 = NextItem(mat1);

          /* Cycle loop */

          continue;
        }

      /* Loop over outflow */
      
      loc0 = (long)RDB[mat1 + MATERIAL_PTR_OUTFLOW];
      while (loc0 > VALID_PTR)
        {
          /* Check mode */

          if ((long)RDB[loc0 + REPROC_CON_MODE] != 2)
            {
              /* Pointer to next */
              
              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Check that pointer to reprocessor matches that given for */
          /* depletion interval */

          if (rep != (long)RDB[loc0 + REPROC_CON_PTR_REP])
            {
              /* Pointer to next */
              
              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Reset totals */

          loc1 = (long)RDB[loc0 + REPROC_CON_PTR_MFLOW];
          while (loc1 > VALID_PTR)
            {
              /* Pointer to original */

              ptr = (long)RDB[loc1 + MFLOW_LIST_PTR_ORIG];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              
              /* Reset value */

              WDB[ptr + MFLOW_LIST_TOT] = 0.0;

              /* Pointer to next */

              loc1 = NextItem(loc1);
            }

          tot = 0.0;

          /* Calculate totals */

          loc1 = (long)RDB[loc0 + REPROC_CON_PTR_MFLOW];
          while (loc1 > VALID_PTR)
            {
              /* Pointer to nuclide */

              iso1 = (long)RDB[loc1 + MFLOW_LIST_PTR_ISO0];
              CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);

              /* Pointer to original */

              ptr = (long)RDB[loc1 + MFLOW_LIST_PTR_ORIG];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              /* Add to total */

              WDB[ptr + MFLOW_LIST_TOT] = RDB[ptr + MFLOW_LIST_TOT] +
                RDB[iso1 + COMPOSITION_ADENS];
              
              tot = tot + RDB[iso1 + COMPOSITION_ADENS];

              /* Pointer to next */

              loc1 = NextItem(loc1);
            }

          /* Check if total is zero */

          if (tot == 0.0)
            {
              /* Pointer to next */
              
              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }

          /* Loop over list */

          loc1 = (long)RDB[loc0 + REPROC_CON_PTR_MFLOW];
          while (loc1 > VALID_PTR)
            {
              /* Get ZAI and removal rate */
              
              ZAI = (long)RDB[loc1 + MFLOW_LIST_ZAI];
              rr = RDB[loc1 + MFLOW_LIST_RATE];
              CheckValue(FUNCTION_NAME, "rr", "", rr, 0.0, INFTY);

              /* Check source pointer */

              if ((long)RDB[loc0 + REPROC_CON_PTR_MAT1] != mat1)
                Die(FUNCTION_NAME, "Pointer error");

              /* Get source nuclide pointer */

              iso1 = (long)RDB[loc1 + MFLOW_LIST_PTR_ISO0];
              CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);
              
              nuc1 = (long)RDB[iso1 + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY, nuc1);

              /* Get target material and nuclide pointers */

              mat2 = (long)RDB[loc0 + REPROC_CON_PTR_MAT2];
              CheckPointer(FUNCTION_NAME, "(mat1)", DATA_ARRAY, mat2);
              
              iso2 = (long)RDB[loc1 + MFLOW_LIST_PTR_ISO1];
              CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso2);

              nuc2 = (long)RDB[iso2 + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY, nuc2);

              /* Check that ZAI's match */

              if ((long)RDB[nuc1 + NUCLIDE_ZAI] != 
                  (long)RDB[nuc2 + NUCLIDE_ZAI])
                Die(FUNCTION_NAME, "Mismatch in ZAI");

              /* Get material indexes */
          
              im1 = (long)RDB[mat1 + MATERIAL_FLOW_IDX] - 1;
              CheckValue(FUNCTION_NAME, "im1", "", im1, 0, nm - 1);
              
              im2 = (long)RDB[mat2 + MATERIAL_FLOW_IDX] - 1;
              CheckValue(FUNCTION_NAME, "im2", "", im2, 0, nm - 1);
              
              /* Get nuclide indexes */
              
              idx1 = (long)TestValuePair(nuc1 + NUCLIDE_PTR_MATRIX_IDX, 
                                         (double)mat0, id);
              idx2 = (long)TestValuePair(nuc2 + NUCLIDE_PTR_MATRIX_IDX, 
                                         (double)mat0, id);
              
              /* Convert */
                      
              i = idx1 + im1;
              j = idx2 + im2;
                      
              /* Check indexes */
                      
              if (i == j)
                Die(FUNCTION_NAME, "Shouldn't happen");

              /* Fix by GRi 2017/03/31 */

              rr = rr*RDB[mat1 + MATERIAL_VOLUME]*RDB[mat1 + MATERIAL_ADENS];

              /* Pointer to original */

              ptr = (long)RDB[loc1 + MFLOW_LIST_PTR_ORIG];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
                      
              /* Add to vector */

              if (ZAI == -1)
                {
                  vec[j] = vec[j] + rr*RDB[iso1 + COMPOSITION_ADENS]/tot;
                  vec[i] = vec[i] - rr*RDB[iso1 + COMPOSITION_ADENS]/tot;
                }
              else if (RDB[ptr + MFLOW_LIST_TOT] > 0.0)
                {
                  vec[j] = vec[j] + rr*RDB[iso1 + COMPOSITION_ADENS]/
                    RDB[ptr + MFLOW_LIST_TOT];
                  vec[i] = vec[i] - rr*RDB[iso1 + COMPOSITION_ADENS]/
                    RDB[ptr + MFLOW_LIST_TOT];
                }

              /* Pointer to next */

              loc1 = NextItem(loc1);
            }
          
          /* Next */
          
          loc0 = NextItem(loc0);
        }

      /* Pointer to next material */

      mat1 = NextItem(mat1);
    }

  /* Read vector into matrix */
  
  for (j = 0; j < sz; j++)
    if (vec[j] != 0.0)
      {
        /* Check value */

        CheckValue(FUNCTION_NAME, "vec[j]", "", vec[j], -INFTY, INFTY);
        
        /* Put value in matrix */

        val[nnz].re = vec[j];
        val[nnz].im = 0.0;

        /* Put row index */
            
        if (nnz > nsz - 1)
          Die(FUNCTION_NAME, "Index error (1)");
        else
          row[nnz] = j; 

        /* Add counter */
        
        nnz++;
        
        /* Check value */
        
        if (nnz > nsz)
          Die(FUNCTION_NAME, "nnz > nsz");
      }

  /* Put column index */

  col[sz + 1] = nnz; 
  
  /***************************************************************************/

  /* Free temporary array */

  Mem(MEM_FREE, vec);

  /* Set matrix size */
  
  A->nnz = nnz;

  /* Palauta matriisi A */

  return A;
}

/*****************************************************************************/
