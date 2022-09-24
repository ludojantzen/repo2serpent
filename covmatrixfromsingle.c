/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : covmatrixfromsingle.c                          */
/*                                                                           */
/* Created:       2018/07/04 (VVa)                                           */
/* Last modified: 2018/10/19 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates a covariance matrix based on single ZAI-mt pair pair */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CovMatrixFromSingle:"

/*****************************************************************************/

void CovMatrixFromSingle(long *ZAImtArray, long **COVxy, long i, long j,
                         long matrix, long blockidx)
{
  long cmtx, zaimt1, zaimt2, zai1, zai2, mt1, mt2, ng, x, y, ptr, datapos;
  long rowpos, cursz, transpose = 0, ii, jj;
  double val;
#ifdef XYZAAA
  char tmpstr[MAX_STR];
  FILE *fout;
#endif

  /* Get pointer to the first covariance matrix or return */

  if ((cmtx = (long)RDB[DATA_PTR_COVMTX0]) < VALID_PTR)
    return;

  /* Get number of energy groups */

  ng = (long)RDB[cmtx + COVMTX_NG];

  /* Loop over the block */

  /* Get the zaimt pair for this covariance */
  /* the column (j) corresponds to the x in the COV(x,y) */
  /* the row    (i) corresponds to the y in the COV(x,y) */

  zaimt1 = ZAImtArray[j];
  zaimt2 = ZAImtArray[i];

  /* Get the covariance matrix that constitutes the current        */
  /* element in the block we are processing, the block is not full */
  /* so we'll simply skip covariances for which no data exists     */


  if (i == j)
    cursz = 1;
  else
    cursz = 2;

  for (ii = 0; ii < cursz; ii++)
    {
      for (jj = 0; jj < cursz; jj++)
        {
          /* Skip diagonal for off-diagonal blocks */

          if (cursz > 1)
            if (ii == jj)
              continue;

          /* Check that the correct covariance matrix was found */

          if (ii > jj)
            {
              cmtx = COVxy[j][i];
              zaimt1 = ZAImtArray[i];
              zaimt2 = ZAImtArray[j];
            }
          else
            {
              cmtx = COVxy[i][j];
              zaimt1 = ZAImtArray[j];
              zaimt2 = ZAImtArray[i];
            }

          if (cmtx < VALID_PTR)
            Die(FUNCTION_NAME, "Covariance between %ld %ld not found from full block",
                zaimt1, zaimt2);

          /* Get first ZAI and mt*/

          zai1 = (long)RDB[cmtx + COVMTX_ZAI1];
          mt1 = (long)RDB[cmtx + COVMTX_MT1];

          /* Get second ZAI and mt*/

          zai2 = (long)RDB[cmtx + COVMTX_ZAI2];
          mt2 = (long)RDB[cmtx + COVMTX_MT2];

          /* Check if the pair matches without inversion (transpose) */
          /* with it or not at all. */

          if ((zai1*10000+mt1 == zaimt1) && (zai2*10000+mt2 == zaimt2))
            transpose = 0;
          else if ((zai1*10000+mt1 == zaimt2) && (zai2*10000+mt2 == zaimt1))
            transpose = 1;
          else
            Die(FUNCTION_NAME, "Covariance matrix contained (%ld, %ld) but was "
                "looking for (%ld, %ld).", zai1*10000+mt1, zai2*10000+mt2,
                zaimt1, zaimt2);

          /**************************************************************************/
          /* Copy data from the covariance matrix to the big matrix we are building */
          /**************************************************************************/

          /*******************************************************************/
          /* As we are copying the data, we will at the same time generate   */
          /* the closest symmetric matrix: Self-covariance (COV(x,x)) should */
          /* be symmetric but typically isn't quite (due to finite precision */
          /* used in data generation calculations?).                         */
          /*******************************************************************/

          /* Get pointer to the data of the covariance matrix */

          ptr = (long)RDB[cmtx + COVMTX_PTR_DATA];

          /* Calculate beginning of this matrix in the bigger matrix */
          /* Column jj, row ii */

          /* We have written ii rows corresponding to ii*cursz */
          /* covariance matrices that each have a size ng*ng */

          /* The position on current line is jj*ng (first line for this matrix) */

          datapos = ii*cursz*ng*ng + jj*ng;

          for (y = 0; y < ng; y++)
            {
              /* Get beginning position of this row in the larger matrix */

              rowpos = datapos + y*cursz*ng;

              for (x = 0; x < ng; x++)
                {

                  if (transpose == 0)
                    WDB[matrix + rowpos + x] = RDB[ptr + ng*y + x];
                  else
                    WDB[matrix + rowpos + x] = RDB[ptr + ng*x + y];

                  /* If we are working on the diagonal, we'll correct the */
                  /* asymmetry here as we have the indices at hand */

                  if (i == j)
                    {
                      /* Get current value */

                      val = RDB[matrix + rowpos + x];

                      /* Take average of matrix and transpose */

                      if (transpose == 0)
                        WDB[matrix + rowpos + x] = 0.5*(val + RDB[ptr + ng*x + y]);
                      else
                        WDB[matrix + rowpos + x] = 0.5*(val + RDB[ptr + ng*y + x]);
                    }
                }
            }
        }
    }
#ifdef XYZAAA
  /* Print out the data to a file */

  sprintf(tmpstr, "COVblck_%ld.txt", blockidx+1);
  fout = fopen(tmpstr, "w");

  /* Loop over the block (rows first to operate in column major order) */

  /* Calculate beginning of this matrix in the bigger matrix */

  for (x = 0; x < ng*cursz*ng*cursz; x++)
    fprintf(fout, "%E ", RDB[matrix + x]);

  fclose(fout);
#endif
}
