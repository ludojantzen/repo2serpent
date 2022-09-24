/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : decaypointprecdet.c                            */
/*                                                                           */
/* Created:       2015/05/05 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Handles decay of pointwise precursors over current time      */
/*              interval                                                     */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DecayPointPrecdet:"

/*****************************************************************************/

void DecayPointPrecDet()
{
  long loc0, part;
  double dt, lambda, wgt;
  double t0, t1;
  

  /***************************************************************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* Check precursor tracking mode */

  if (RDB[DATA_PRECURSOR_TRANSPORT_MODE] != PREC_MODE_POINT)
    return;

#ifdef DNPRINT
  fprintf(outp, "decaypointprecdet.c -->\n");  
#endif

  /* Now that we have done all the sampling, we should decay the weights for the EOI */
  /* Loop over precursor list to decay the weights over this interval */

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  /* Calculate time interval length */

  dt = t1 - t0;

  /* Get pointer to first precursor (dummy) */

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

      /* Calculate weight to survive next interval */

      wgt = wgt*exp(-lambda*dt);

      /* Store new weight */

      WDB[part + PARTICLE_WGT] = wgt;

      /* Set current time of precursor */

      WDB[part + PARTICLE_T] = t1;

      /* Next precursor */

      part = NextItem(part);
    }    

#ifdef DNPRINT
  fprintf(outp, "<-- decaypointprecdet.c\n\n");  
#endif

  /***************************************************************************/
}

/*****************************************************************************/
