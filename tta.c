/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tta.c                                          */
/*                                                                           */
/* Created:       2011/06/10 (AIs)                                           */
/* Last modified: 2012/05/31 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: solves depletion steps by the (variation) TTA method.        */
/*                                                                           */
/* Comments: This function is intented for development purposes only.        */
/*           Use CRAM (MatrixExponential) instead.                           */
/*                                                                           */
/*           Tää on paljon hitaampi kuin 1.1.** versioissa, johtune A:sta    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TTA:"

#define MAX_TTA_CHAIN 500  /* ketjun maksimipituus */

/*****************************************************************************/

double *TTA(struct ccsMatrix *A, double *N0, double t)
{
  long iso,j;
  double adens, tot_adens, lambda;
  double l[MAX_TTA_CHAIN];
  double Plimit;
  double *N;

  /*varaa muisti tuloksille*/

  N = (double *)Mem(MEM_ALLOC, A->n, sizeof(double));

  /*
  N = WorkArray(DATA_PTR_WORK_COMP2, PRIVA_ARRAY, A->n, OMP_THREAD_NUM);
  */

  /*kokonaisatomitiheys*/
  tot_adens=0;
  for (iso=0; iso< A->n; iso++)
    {
      tot_adens=tot_adens+N0[iso];
      N[iso]=0.0;
    }

  /* Silmukka nuklidikoostumuksen yli */

  for (iso=0; iso< A->n; iso++)
    {
      /* Atomitiheys aika-askeleen alussa */
      
      adens = N0[iso];

      /* transmutaatiokerroin */

      lambda=0.0;
      for (j=A->colptr[iso] ; j < A->colptr[iso+1] ; j++)
        if(A->rowind[j]==iso)
          {
            lambda=-A->values[j].re;
            break;
          }

      if (adens > 0.0)
        {
          
          /* Ketjun katkaisuraja */
          
          Plimit=RDB[DATA_DEP_TTA_CUTOFF]*tot_adens/adens;
          
          /* Poltetaan... */

          if (Plimit < 1)
            TTALoop(iso, t, 0, l, 1.0, 1.0, adens, Plimit, N, A,iso);
          else
            N[iso]=N[iso] + adens*exp(-lambda*t);
          
        }
    }
  
  return N;
}

/*****************************************************************************/
