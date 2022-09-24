/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : covmatrixfromblock.c                           */
/*                                                                           */
/* Created:       2018/06/18 (VVa)                                           */
/* Last modified: 2018/10/19 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates a covariance matrix based on information of the      */
/*              connected block. The covariance matrix created is a simple   */
/*              data array that contains the covariance data.                */
/*                                                                           */
/* Comments: -  Output matrix is still in a row major order                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CovMatrixFromBlock:"

/*****************************************************************************/

void CovMatrixFromBlock(long *ZAImtArray, long **COVxy, long pos,
                        long cursz, long matrix, long blockidx)
{
  long cmtx, zaimt1, zaimt2, zai1, zai2, mt1, mt2, i, j, ng, x, y, ptr, datapos;
  long rowpos, transpose = 0;
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

  for (i = pos; i < pos+cursz; i++)
    {
      for (j = pos; j < pos+cursz; j++)
        {
          /* Get the zaimt pair for this covariance */
          /* the column (j) corresponds to the x in the COV(x,y) */
          /* the row    (i) corresponds to the y in the COV(x,y) */

          zaimt1 = ZAImtArray[j];
          zaimt2 = ZAImtArray[i];

          /* Get the covariance matrix that constitutes the current        */
          /* element in the block we are processing, the block is not full */
          /* so we'll simply skip covariances for which no data exists     */

          if ((cmtx = COVxy[i][j]) < VALID_PTR)
            continue;

          /* Check that the correct covariance matrix was found */

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
          /* Column j-pos, row i-pos*/

          /* We have written (i-pos) rows corresponding to (i-pos)*cursz */
          /* covariance matrices that each have a size ng*ng */

          /* The position on current line is (j-pos)*ng (first line for this matrix) */

          datapos = (i-pos)*cursz*ng*ng + (j-pos)*ng;

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
