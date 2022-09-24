/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectprecdet.c                               */
/*                                                                           */
/* Created:       2015/05/11 (VVa)                                           */
/* Last modified: 2017/02/01 (VVa)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Collects precursor detectors after cycle or batch            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectPrecDet:"

/*****************************************************************************/

void CollectPrecDet()
{
  long loc0, ptr, stp, i, j, k, nt, ng, nm;
  double val, norm;

  /***************************************************************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Get number of time bins */

  nt = (long)RDB[loc0 + PRECDET_NT];

  /* Get number of group bins */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Get pointer to mesh */

  ptr = (long)RDB[loc0 + PRECDET_PTR_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Calculate number of mesh bins */

  nm = 1;
  nm = nm*(long)RDB[ptr + MESH_N0];
  nm = nm*(long)RDB[ptr + MESH_N1];
  nm = nm*(long)RDB[ptr + MESH_N2];

  /* Get pointer to statistics */

  stp = (long)RDB[loc0 + PRECDET_PTR_STAT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  /* Get normalization */

  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {

      norm = NormCoef(PARTICLE_TYPE_NEUTRON);
      CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);

    }
  else
    {
      /* In dynamic mode the values are */
      /* multiplied by normalization already in precdet.c */

      norm = 1.0;
    }

  /* Loop over bins to add stats */

  for (i = 0; i < nt; i++)
    for (j = 0; j < ng; j++)
      for (k = 0; k < nm; k++)
        {
          /* Get value from buffer */

          val = BufVal(stp, i, j, k);

          /* Add a statistic */

          AddStat(val*norm, stp, i, j, k);
        }

}
