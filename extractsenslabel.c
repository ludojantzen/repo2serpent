/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : extractsenslabel.c                             */
/*                                                                           */
/* Created:       2017/03/21 (VVa)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Extracts sensitivity event sub indices based on compressed   */
/*              label.                                                       */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ExtractSensLabel:"

/*****************************************************************************/

void ExtractSensLabel(long label, long *imat,
                  long *izai, long *irea, long *iene, double *hitmiss)
{
  long loc0, nzai, nrea, nene, nmu, nxyz, imu, ixyz;

  /* Get pointer to sensitivity block or return */

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Get sizes of sub bins */

  nzai = (long)RDB[loc0 + SENS_N_ZAI];
  nrea = (long)RDB[loc0 + SENS_N_PERT] + 1;
  nene = (long)RDB[loc0 + SENS_N_ENE] + 1;
  nmu  = (long)RDB[loc0 + SENS_N_MU];
  nxyz = (long)RDB[loc0 + SENS_N_SPAT];

  /* extract label sub indices based on the collapsed index */

  if (label > 0)
    {
      *hitmiss = 1.0;
    }
  else
    {
      *hitmiss = -1.0;
      label = -label;
    }

  label = label - 1;

  /* Extract material index */

  *imat = label / (nzai*nrea*nene*nxyz*nmu);
  label = label % (nzai*nrea*nene*nxyz*nmu);

  /* Extract nuclide index */

  *izai = label / (nrea*nene*nxyz*nmu);
  label = label % (nrea*nene*nxyz*nmu);

  /* Extract reaction index */

  *irea = label / (nene*nxyz*nmu);
  label = label % (nene*nxyz*nmu);

  /* Extract energy index */

  *iene = label / (nxyz*nmu);
  label = label % (nxyz*nmu);

  /* Extract spatial index */

  ixyz = label / (nmu);
  label = label % (nmu);

  /* The rest of the label is the angular index */

  imu = label;

  /* Check the sub indices */

  CheckValue(FUNCTION_NAME, "imat", "", *imat, 0, (long)RDB[loc0 + SENS_N_MAT]);
  CheckValue(FUNCTION_NAME, "izai", "", *izai, 0, nzai);
  CheckValue(FUNCTION_NAME, "irea", "", *irea, 0, nrea);
  CheckValue(FUNCTION_NAME, "iene", "", *iene, 0, nene);
  CheckValue(FUNCTION_NAME, "ixyz",  "", ixyz, 0, nxyz);
  CheckValue(FUNCTION_NAME, "imu",  "", imu, 0, nmu);
}
