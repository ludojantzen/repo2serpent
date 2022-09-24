/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixfree.c                                */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Vapauttaa matriisille varatun muistin                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixFree:"

/*****************************************************************************/

void ccsMatrixFree(struct ccsMatrix *cmat)
{
  if (cmat == NULL)
    {
      Warn(FUNCTION_NAME, "cmat = NULL");
      return;
    }

  if (cmat->colptr == NULL)
    Warn(FUNCTION_NAME, "cmat->colptr = NULL");
  else
    Mem(MEM_FREE, cmat->colptr);
  
  if (cmat->rowind == NULL)
    Warn(FUNCTION_NAME, "cmat->rowind = NULL");
  else
    Mem(MEM_FREE, cmat->rowind);

  if (cmat->colind != NULL)
    Mem(MEM_FREE, cmat->colind);

  if (cmat->rowptr != NULL)
    Mem(MEM_FREE, cmat->rowptr);

  if (cmat->next != NULL)
    Mem(MEM_FREE, cmat->next);

  if (cmat->values == NULL)
    Warn(FUNCTION_NAME, "cmat->values = NULL");
  else
    Mem(MEM_FREE, cmat->values);

  Mem(MEM_FREE, cmat);
  cmat = NULL;
}

/*****************************************************************************/
