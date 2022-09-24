/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : SampleDelnu.c                                  */
/*                                                                           */
/* Created:       2015/05/05 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Samples delayed neutrons from precursor concentrations at    */
/*              the beginning of an interval and adds them to que or source  */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleDelnu:"

/*****************************************************************************/

void SampleDelnu()
{
  long id, loc0, toemit, nn, idx;

  /***************************************************************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  if (RDB[loc0 + PRECDET_N_EMIT] == 0)
    return;

#ifdef DNPRINT
  fprintf(outp, "sampledelnu.c -->\n");
#endif

  /* Get number of neutrons to emit */
  /* Number of live fission chains has to be taken in account */

  toemit = (long)RDB[loc0 + PRECDET_N_EMIT];

#ifdef DNPRINT
  fprintf(outp, "Weight of neutrons to emit %f\n",
          RDB[loc0 + PRECDET_W_EMIT]/RDB[loc0 + PRECDET_N_EMIT]);
  fprintf(outp, "Sampling %ld delayed neutrons for the interval\n", toemit);
#endif

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    {
      /* Open buffer for writing */

      WDB[DATA_BUF_REDUCED] = (double)NO;

      /* Start parallel timer */

      StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn)
#endif
      {
        /* Get Open MP thread id */

        id = OMP_THREAD_NUM;

#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif
        /* Loop over source neutrons */

        for (nn = 0; nn < toemit; nn++)
          {

            /* Calculate particle index */
            /* TODO: These have to be calculated */

            idx = (long)RDB[DATA_NHIST_TOT];
            idx = idx + (long)(RDB[DATA_SRC_POP]) + nn;

            /* Sample source point */

            SamplePointDelnu(id, nn, idx);

          }
      }

      /* Stop parallel timer */

      StopTimer(TIMER_OMP_PARA);
    }
  else
    {

      /* Start parallel timer */

      StartTimer(TIMER_OMP_PARA);


#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn)
#endif
      {
        /* Get Open MP thread id */

        id = OMP_THREAD_NUM;

#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif
        /* Loop over source neutrons */

        for (nn = 0; nn < toemit; nn++)
          {

            /* Calculate particle index */
            /* TODO: These have to be calculated */

            idx = (long)RDB[DATA_NHIST_TOT];
            idx = idx + (long)(RDB[DATA_SRC_POP]) + nn;

            /* Sample source point */

            SampleMeshDelnu(id, nn, idx);

          }
      }

      /* Stop parallel timer */

      StopTimer(TIMER_OMP_PARA);

    }

#ifdef DNPRINT
  fprintf(outp, "<-- sampledelnu.c\n\n");
#endif
  /***************************************************************************/
}

/*****************************************************************************/
