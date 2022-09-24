/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getsenseventtype.c                             */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Figures out the event type (index) based on the sampled      */
/*              reaction.                                                    */
/*                                                                           */
/* Comments: - The structure of this subroutine depends on collision.c.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetSensEventType"

/*****************************************************************************/

void GetSensEventType(long mat, long rea, double E, long *reaIdx, long *reaType,
                     long *isfission)
{
  long ptr, sens, mt;

  /* Get pointer to sensitivity block or return */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  /* Set fission flag */

  if (((mt > 17) && (mt < 22)) || (mt == 38))
    *isfission = 1;
  else
    *isfission = 0;

  /* Reset reaction index */

  *reaIdx = 0;
  *reaType = 0;

  /* Check reaction type */

  if (mt == 2)
    {
      /* Check sampling */

      if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
        return;

      *reaType = (long)ELA_SCATT_IDX;

    }
  else if ((mt == 1002) || (mt == 1004))
    {
      /* Check sampling */

      if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
        return;

      *reaType = (long)SAB_SCATT_IDX;

    }
  else if ((mt == 2002) || (mt == 2004))
    {
      /* S(a,b) scattering with on-the-fly interpolation */

      /* Check sampling */

      if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
        return;

      Die(FUNCTION_NAME, "On-the-fly S(a,b) not tested with sensitivity calculations");

      *reaType = (long)SAB_SCATT_IDX;
    }
  else if (RDB[rea + REACTION_TY] == 0.0)
    {
      /* Capture */

      if ((long)RDB[DATA_NPHYS_SAMPLE_CAPT] == NO)
        return;

      *reaType = (long)CAPT_IDX;

    }
  else if (fabs(RDB[rea + REACTION_TY]) > 100.0)
    {
      /* Complex reaction */

      Warn(FUNCTION_NAME, "Complex reactions not tested with Sens");

      return;
    }
  else if (((mt > 17) && (mt < 22)) || (mt == 38))
    {
      *reaType = (long)FISS_IDX;
    }
  else if ((mt > 50) && (mt < 91))
    {
      /* Inelastic level scattering */

      if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
        return;

      *reaType = (long)INL_SCATT_IDX;

    }
  else if (RDB[rea + REACTION_WGT_F] > 1.0)
    {
      /* NxN reaction */

      if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
        return;

      *reaType = (long)NXN_SCATT_IDX;
    }
  else if (mt < 100)
    {
      /* Continuum single neutron inelastic scattering */

      if ((long)RDB[DATA_NPHYS_SAMPLE_SCATT] == NO)
        return;

      *reaType = (long)INL_SCATT_IDX;
    }
  else
    {
      /* Unknown reaction mode */

      return;
    }

  if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS)
    {
      /* Get pointer to reaction type index list */

      ptr = (long)RDB[sens + SENS_PTR_PERT_INDICES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      *reaIdx = (long)RDB[ptr + *reaType];
    }
  else if ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT)
    {
      /* Check if perturbations are by partial reactions */

      ptr = (long)RDB[sens + SENS_PTR_MT_INDICES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Read index from pre-built array */

      *reaIdx = (long)RDB[ptr + mt];

    }
}
