/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizeprecdet.c                             */
/*                                                                           */
/* Created:       2015/11/05 (VVa)                                           */
/* Last modified: 2018/09/28 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets simulation normalization and scales weights of live     */
/*              neutrons and precursors at the beginning of dynamic          */
/*              simulation                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizePrecDet:"
#define DNPRINT
/*****************************************************************************/

void NormalizePrecDet()
{
  long loc0, ptr, part;
  long stp, j, k;
  long ng, tb, pop, toemit, tolive;
  long n0, n1, n2;
  double val, lambda, t0, t1, dt, Wemit, Wlive, Wall;
  double Na;
  double wgt, Wemit_list;


  /* Get pointer to precursor detector or return */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* We don't need to do this for criticality source mode */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    return;

#ifdef DNPRINT
  fprintf(outp, "normalizeprecdet.c -->\n");
#endif

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  /* Calculate length of time interval */

  dt = t1 - t0;

  /* Get current time bin */

  tb = (long)RDB[DATA_DYN_TB];

  /************************************************************/
  /*   Calculate number of physical delayed neutrons to emit  */
  /*   based on the stable populations in the mesh detector   */
  /************************************************************/
  /* Read group binning from precursor detector */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Get pointer to precursor spatial mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of mesh bins */

  n0 = (long)RDB[ptr + MESH_N0];
  n1 = (long)RDB[ptr + MESH_N1];
  n2 = (long)RDB[ptr + MESH_N2];

  /* Get pointer to decay constant array */

  ptr = (long)RDB[loc0 + PRECDET_PTR_LAM_ARRAY];

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];

  /* Reset total weight to emit */

  Wall = 0.0;
  Wemit = 0.0;
  Na = 0;

  /* Close buffer for reading */

  WDB[DATA_BUF_REDUCED] = (double)YES;

  /* Loop over spatial and group bins to calculate */
  /* values to emit and tally during interval      */

  for (j = 0; j < ng; j++)
    for (k = 0; k < n0*n1*n2; k++)
      {

        /* Get value of the buffer */

        val = BufVal(stp, tb, j, k);

        /* Get decay constant for this group */

        lambda = RDB[ptr + j];

        /* Calculate amount to emit on next interval */

        Wemit = val*(1-exp(-lambda*dt));

        /* Add weight to emit to total counter */

        Wall = Wall + Wemit;
      }

  /* Printing for debugging */

#ifdef DNPRINT
  fprintf(outp, "Should emit %E delayed neutrons in total\n", Wall);
#endif

  /* Loop over precursor list to calculate total weight to be emitted from them */
  /* Maybe put this if-block in a separate subroutine? */

  /* Reset weight to emit from list */

  Wemit_list = 0.0;

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    {
      /* Get pointer to precursor list */
      /* (first item is a dummy)       */

      part = (long)RDB[DATA_PART_PTR_PSOURCE];
      CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

      /* First precursor after dummy */

      part = NextItem(part);

      while (part > VALID_PTR)
        {

          /* Get weight of precursor */

          wgt = RDB[part + PARTICLE_WGT];

          /* Get decay constant */

          lambda = RDB[part + PARTICLE_DN_LAMBDA];

          /* Calculate weight to emit on next interval */

          wgt = wgt*(1-exp(-lambda*dt));

          /* Add to weight to emit from list */

          Wemit_list += wgt;

          /* Store weight to emit on first interval */

          WDB[part + PARTICLE_U] = wgt;

          /* Next precursor */

          part = NextItem(part);
        }

#ifdef DNPRINT
      fprintf(outp, "Would emit %E delayed neutrons from list\n", Wemit_list);
      fprintf(outp, "List weights have to be scaled by %E\n", Wall/Wemit_list);
#endif

      /* Loop over precursor list to scale weights */

      /* Get pointer to precursor list */
      /* (first item is a dummy)       */

      part = (long)RDB[DATA_PART_PTR_PSOURCE];
      CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

      /* Reset number of neutrons calculated for activity */

      Na = 0;

      /* First precursor after dummy */

      part = NextItem(part);

      while (part > VALID_PTR)
        {

          /* Get weight of precursor */

          wgt = RDB[part + PARTICLE_WGT];

          /* Scale weight */

          WDB[part + PARTICLE_WGT] = wgt*Wall/Wemit_list;

          /* Add to activity count */

          Na++;

          /* Get weight to emit on next interval */

          wgt = RDB[part + PARTICLE_U];

          /* Scale weight to emit on next interval */
          /* This is actually proportion of weight emitted by this precursor  */
          /* of the weight emitted by all precursors on the next time interval */

          WDB[part + PARTICLE_U] = wgt/Wemit_list;

          if (wgt/Wemit_list > 1.0)
            Die(FUNCTION_NAME, "WTF?!");

          /* Next precursor */

          part = NextItem(part);
        }

      Wemit_list = Wall;

    }

  /* Calculate number of live neutrons and number to emit for first interval */
  /* For further intervals this is done in resizedynsrc.c */

  /* Get physical number of live neutrons */

  Wlive = RDB[loc0 + PRECDET_W_LIVE];

  /* Get total initial population */
  /* MPI mode is taken in account in transportcycle.c */

  pop = (long)RDB[DATA_SRC_POP];

  /* Calculate number of delayed neutrons to emit */

  toemit = (long)round(Wall/(Wlive + Wall)*(double)pop);

  /* Calculate number of live neutrons to get from srcfile */

  tolive = (long)round(Wlive/(Wlive + Wall)*(double)pop);

  /* NB: Since we are rounding the numbers there, our live neutrons and */
  /* delayed neutrons must have slightly different weights to keep      */
  /* consistency. */

#ifdef DNPRINT
  fprintf(outp, "To emit %E (%f %%): %ld\n", Wall, Wall/(Wlive + Wall)*100, toemit);
  fprintf(outp, "To live %E (%f %%): %ld\n", Wlive, Wlive/(Wlive + Wall)*100, tolive);
#endif

  /* Store numbers */

  WDB[loc0 + PRECDET_N_EMIT] = (double)toemit;
  WDB[loc0 + PRECDET_N_LIVE] = (double)tolive;

  /* Live neutrons will have weight 1.0 initially, calculate   */
  /* normalization based on that. Weight of delayed neutrons   */
  /* will be calculated in SampleDelNu in a way that conserves */
  /* the physical number */

  if (tolive > 0)
    {

      /* Set normalization */

      WDB[DATA_NORM_COEF_N] = Wlive / (double)tolive;

    }
  else
    {

      /* If there are no live neutrons initially we can set */
      /* normalization based on delayed neutrons */

      WDB[DATA_NORM_COEF_N] = Wall / (double)toemit;

    }

  /* Store live weight (PRECDET_W_LIVE will no longer be physical */
  /* neutrons */

  WDB[loc0 + PRECDET_W_LIVE] = Wlive/RDB[DATA_NORM_COEF_N];

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    {
      /* Use weights calculated based on precursor list */

      /* Store weight to emit */

      WDB[loc0 + PRECDET_W_EMIT] = Wemit_list/RDB[DATA_NORM_COEF_N];

      /* Store average emission */
      /* Used as a threshold for Russian Roulette for */
      /* saving new precursors in precdet.c */

      WDB[loc0 + PRECDET_AVE_EMIT] = Wemit_list/Na;

#ifdef DNPRINT
      fprintf(outp,
              "Normalization: 1.0 weight corresponds to %E neutrons\n",
              RDB[DATA_NORM_COEF_N]);
      fprintf(outp, "Weight to emit: %f (%E neutrons)\n",
              Wemit_list/RDB[DATA_NORM_COEF_N], Wall);
      fprintf(outp, "Average neutrons emitted was %E\n",
              Wemit_list/Na);
#endif
    }
  else
    {
      /* Use weight calculated from spatial mesh */

      WDB[loc0 + PRECDET_W_EMIT] = Wall/RDB[DATA_NORM_COEF_N];

#ifdef DNPRINT
      fprintf(outp, "Normalization: 1.0 weight corresponds to %E neutrons\n",
        RDB[DATA_NORM_COEF_N]);
      fprintf(outp, "Weight to emit: %f (%E neutrons)\n",
        Wall/RDB[DATA_NORM_COEF_N], Wall);
#endif
    }

  /* Print some things */

#ifdef DNPRINT
  fprintf(outp, "Weight of live neutrons    %f\n", 1.0);
  fprintf(outp, "<-- normalizeprecdet.c\n\n");
#endif

}
