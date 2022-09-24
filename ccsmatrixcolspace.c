/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixcolspace.c                            */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2011/12/16 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Varaa lisää tilaa sarakkeille                                */
/*                                                                           */
/* Comments: - Kutsutaan SymbolicLU:sta                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixColspace:"

/*****************************************************************************/

void ccsMatrixColSpace(struct ccsMatrix *cmat, long i, long rs)
{
  long j, nnzi, n, nnz;
  long *tempc;
  complex *tempv;

  /* JLE: jos matriisin koko on pieni, toi rs voi mennä nollaksi, mikä */
  /* aiheuttaa segmentation faultin myöhemmin. Korjataan se tässä      */
  /* ykköseksi (en tiedä toimiiko kaikissa mahdollisissa tilanteissa)  */
  /* 22.5.2011. */

  if (rs == 0)
    rs = 1;

  n   = cmat->n;
  nnz = cmat->nnz;

  /* lisää tilaa rs uudelle fill-inille */

  nnz  = nnz + rs; 

  /* nnzi = nollasta eroavien lkm sarakkeessa i */

  nnzi = nnz - cmat->colptr[i+1];   

  cmat->values = (complex *)Mem(MEM_REALLOC, cmat->values, nnz*sizeof(complex));
  cmat->rowind = (long *)Mem(MEM_REALLOC, cmat->rowind, nnz*sizeof(long));

  /* viimeiset nnzi arvoa */

  tempv = (complex *)Mem(MEM_ALLOC, nnzi, sizeof(complex));

  memcpy(tempv + rs, cmat->values + cmat->colptr[i + 1], 
	 (nnzi - rs)*sizeof(complex));

  memcpy(cmat->values + cmat->colptr[i + 1], tempv, nnzi*sizeof(complex));

  Mem(MEM_FREE, tempv);
  
  tempc = (long *)Mem(MEM_ALLOC, nnzi, sizeof(long));

  memset(tempc, -1, rs*sizeof(long));
  
  memcpy(tempc + rs, cmat->rowind + cmat->colptr[i + 1], 
	 (nnzi - rs)*sizeof(long));

  memcpy(cmat->rowind + cmat->colptr[i + 1], tempc, nnzi*sizeof(long));

  Mem(MEM_FREE, tempc); 

  /* aseta sarakkeiden osoittimet */

  for (j=i+1 ; j <= n ; j++){
    cmat->colptr[j] = cmat->colptr[j] + rs;
  }  

  cmat->nnz = nnz;
}

/*****************************************************************************/
