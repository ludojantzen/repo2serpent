/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printdepmatrix.c                               */
/*                                                                           */
/* Created:       2011/06/01 (JLe)                                           */
/* Last modified: 2019/11/08 (ARi)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Prints depletion matrix for debugging                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintDepMatrix:"

/*****************************************************************************/

void PrintDepMatrix(long mat, struct ccsMatrix *A, double t, double *N0, 
                    double *N1, long id)
{
  long iso, nuc, i, j;
  char tmpstr[MAX_STR];
  FILE *fp;
 
  /* Check print flag */

  if ((long)RDB[DATA_BURN_PRINT_DEPMTX] == NO)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* File name */
          
  sprintf(tmpstr,"%s_depmtx_%s%ld.m", GetText(DATA_PTR_INPUT_FNAME),
          GetText(mat + MATERIAL_PTR_NAME), (long)RDB[DATA_BURN_STEP]);

  /* Open file for writing */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");
  
  /* Check that matrix is square */

  if (A->n != A->m)
    Die(FUNCTION_NAME, "Matrix is not square");

  /* Print time step */

  /*  fprintf(fp, "t = %E;\n", t);*/
  fprintf(fp, "t = %22.15E;\n", t);

  /* Print size */

  /*
  fprintf(fp, "%ld %ld\n", A->n, A->m);
  */
  
  /* Print flux */

  fprintf(fp, "flx = %22.15E;\n", Truncate(RDB[mat + MATERIAL_BURN_FLUX_SSA], 6));

  /* Loop over composition and print */

  fprintf(fp, "N0 = zeros(%ld, 1);\n", A->n);

  iso = (long)RDB[mat + MATERIAL_PTR_COMP]; 
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Get nuclide index */
      
      i = (long)TestValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX, (double)mat, id);

      /* Check and print */

      if ((i < 0) || (i > A->n - 1))
        Die(FUNCTION_NAME, "Something wrong here...");
      else 
        fprintf(fp, "N0(%3ld, 1) = %22.15E; %% %s\n", i + 1, 
                /*        fprintf(fp, "N0(%3ld, 1) =  %1.10E; %% %s\n", i + 1, */
                RDB[iso + COMPOSITION_ADENS],
                GetText(nuc + NUCLIDE_PTR_NAME));

      /* Compare to initial composition */

      if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
        if (RDB[iso + COMPOSITION_ADENS] != N0[i])
          Die(FUNCTION_NAME, "Error in N0");

      /* Next isotope */

      iso = NextItem(iso);
    }

  fprintf(fp, "ZAI = zeros(%ld, 1);\n", A->n);

  iso = (long)RDB[mat + MATERIAL_PTR_COMP]; 
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Get nuclide index */
      
      i = (long)TestValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX, (double)mat, id);

      /* Check and print */

      if ((i < 0) || (i > A->n - 1))
        Die(FUNCTION_NAME, "Something wrong here...");
      else 
        fprintf(fp, "ZAI(%3ld) =  %ld;\n", i + 1,
                (long)RDB[nuc + NUCLIDE_ZAI]);

      /* Compare to initial composition */
 
      if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
        if (RDB[iso + COMPOSITION_ADENS] != N0[i])
          Die(FUNCTION_NAME, "Error in N0");

      /* Next isotope */

      iso = NextItem(iso);
    }

  /* Print matrix size and number of non-zero elements */
  /*
  fprintf(fp, "%ld %ld\n", A->n, A->nnz);
  */
  /* Print matrix */

  fprintf(fp, "A = zeros(%ld, %ld);\n", A->n, A->m);

  for (i = 0; i < A->m; i++)
    for (j = A->colptr[i]; j < A->colptr[i + 1]; j++)
      /*      fprintf(fp, "A(%4ld, %4ld) = %17.10E;\n", A->rowind[j] + 1, i + 1, */
      fprintf(fp, "A(%4ld, %4ld) = %22.15E;\n", A->rowind[j] + 1, i + 1, 
              A->values[j].re);

  /* Print final composition */

  fprintf(fp, "N1 = zeros(%ld, 1);\n", A->n);

  for (i = 0; i < A->n; i++)
    /*    fprintf(fp, "N1(%3ld, 1) = %1.10E;\n", i + 1, N1[i]);*/
    fprintf(fp, "N1(%3ld, 1) = %22.15E;\n", i + 1, N1[i]);
          
  /* Close file */
          
  fclose(fp);
}

/*****************************************************************************/
