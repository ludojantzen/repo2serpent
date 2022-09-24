/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : walkeralias.c                                  */
/*                                                                           */
/* Created:       2015/09/18 (TKa)                                           */
/* Last modified: 2017/02/23 (TKa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Walker's alias method for sampling a discrete distribution   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WalkerAlias:"

/*****************************************************************************/

void WalkerAliasInit(const double *p0, long N, double *alias1, double *alias2,
                       double *cutoff) {
  /* Initialize Walker's alias method - creates the alias and cutoff arays */

  long i, j, idxmin, idxmax;
  long *idxAvailable;
  double sum, pN, pmin, pmax;
  double *p;

  /* Check size */
  if (N <= 0)
    Die(FUNCTION_NAME, "non-positive array size");

  /* Allocate memory */
  p = (double *)Mem(MEM_ALLOC, N, sizeof(double));
  idxAvailable = (long *)Mem(MEM_ALLOC, N, sizeof(long));

  /* Normalize distribution */
  sum = 0.0;
  for (i = 0; i < N; i++) {
    if (p0[i] < 0.0)
      Die(FUNCTION_NAME, "negative probability");
    sum += p0[i];
  }
  for (i = 0; i < N; i++)
    p[i] = p0[i]/sum;

  /* Initialize idx array */
  for (i = 0; i < N; i++)
    idxAvailable[i] = 1;

  /* Probability per element */
  pN = 1.0/(double)N;

  /* Create alias and cutoff tables */
  for (i = 0; i < N; i++) {

    /* Find minimum and maximum */
    pmin = INFTY;
    pmax = -INFTY;
    idxmin = -1;
    idxmax = -1;

    for (j = 0; j < N; j++) {
      if (idxAvailable[j]) {
        if (p[j] < pmin) {
          pmin = p[j];
          idxmin = j;
        }
        if (p[j] > pmax) {
          pmax = p[j];
          idxmax = j;
        }
      }
    }

    /* Check minimum and maximum */
    if (idxmin == -1)
      Die(FUNCTION_NAME, "minimum not found");
    if (idxmax == -1)
      Die(FUNCTION_NAME, "maximum not found");

    /* Set cutoff and alias */
    cutoff[i] = pmin/pN;
    alias1[i] = (double)idxmin;
    alias2[i] = (double)idxmax;

    /* Update probability and available indeces */
    p[idxmax] = p[idxmax] - (pN - p[idxmin]);
    idxAvailable[idxmin] = 0;
  }


  /* Check output */
  for (i = 0; i < N; i++) {
    /* Check the sign of probabilities */
    if (p[i] < 0.0)
      Die(FUNCTION_NAME, "negative probability in p");

    /* Check alias arrays */
    if ((alias1[i] < 0) || (alias1[i] > N-1))
      Die(FUNCTION_NAME, "wrong index in alias1");
    if ((alias2[i] < 0) || (alias2[i] > N-1))
      Die(FUNCTION_NAME, "wrong index in alias2");

    /* Check cutoff */
    if ((cutoff[i] < 0) || (cutoff[i] - 1.0 > 1e-10))
      Die(FUNCTION_NAME, "cutoff out of limits %E %ld %ld", cutoff[i]-1, N, i);
  }


  /* Free memory */
  Mem(MEM_FREE, p);
  Mem(MEM_FREE, idxAvailable);
}

/*****************************************************************************/


/*****************************************************************************/

long WalkerAliasSample(const double *alias1, const double *alias2,
                       const double *cutoff, long N, long id) {
  /* Samples from a discrete distribution using Walker's alias method */

  long k;
  double R;

  R = RandF(id)*(double)N;
  k = (long)R;

  if (R - (double)k < cutoff[k])
    return (long)alias1[k];
  else
    return (long)alias2[k];
}

/*****************************************************************************/


/*****************************************************************************/


