/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixisort.c                               */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2011/10/28 (MPu)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Sorttaa matriisin                                            */
/*                                                                           */
/* Comments: - Perustuu insertion sort -algoritmiin                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixIsort:"

/*****************************************************************************/

void ccsMatrixIsort(struct ccsMatrix *cmat)
{
  long l, i, j, itmp, n, nnz, *col, *row;
  complex ctmp, *val;

  n   = cmat->n;
  nnz = cmat->nnz;

  /* Funktio kirjoitettu alunperin crs-matriisille */
  
  col     = cmat->rowind; /* crs -> ccs */
  row     = cmat->colptr; /* crs -> ccs */
  val     = cmat->values;

  for (i = 0; i < n; i++)
    {
      /* järjestä i. rivi */
      
      for (l = row[i] + 1; l < row[i + 1] ;l++)
	{
	  if (col[l] < 0)
	    break; 
      
	  itmp    = col[l]; 
	  ctmp.re = val[l].re;
	  ctmp.im = val[l].im;
	  
	  j = l; 
    
	  while((col[j - 1] > itmp) && (j > row[i]))
	    {
	      if (j >= nnz)
		Die(FUNCTION_NAME, "j >= nnz");

	      if (i >= n)
		Die(FUNCTION_NAME, "i >= n");

	      col[j]    = col[j-1]; 
	      val[j].re = val[j-1].re; 
	      val[j].im = val[j-1].im; 
	      
	      j--;
	      
	      if (j < 1)
		break;
	    }     

	  col[j]    = itmp; 
	  val[j].re = ctmp.re;
	  val[j].im = ctmp.im; 
	}
    } 
}

/*****************************************************************************/
