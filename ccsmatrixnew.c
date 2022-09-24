/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixnew.c                                 */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2011/12/16 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Varaa tilaa CCS-muotoiselle matriisille                      */
/*                                                                           */
/* Comments: - colind, rowptr ja next varataan findrowindexes.c:ssä          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixNew:"

/*****************************************************************************/

struct ccsMatrix *ccsMatrixNew(long m, long n, long nnz)
{
  struct ccsMatrix *ptr;

  /* Varataan muisti matriisille */

  ptr = (struct ccsMatrix *)Mem(MEM_ALLOC, 1, sizeof(struct ccsMatrix));
  
  /* Varataan muisti datalle */
  
  ptr->values = (complex *)Mem(MEM_ALLOC, nnz, sizeof(complex)); 
  ptr->rowind = (long *)Mem(MEM_ALLOC, nnz, sizeof(long));
  ptr->colptr = (long *)Mem(MEM_ALLOC, n + 1, sizeof(long));

  /* Asetetaan koko */
  
  ptr->m   = m;
  ptr->n   = n;
  ptr->nnz = nnz;

  /* Asetetaan loput pointterit nulliin */

  ptr->rowptr = NULL;
  ptr->next   = NULL;
  ptr->colind = NULL;

  /* Palautetaan pointteri matriisiin */

  return ptr;
}

/*****************************************************************************/

