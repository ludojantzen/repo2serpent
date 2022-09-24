/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : targetvelocity.c                               */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2018/05/25 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Samples velocity of target nuclide with or withoutd DBRC     */
/*                                                                           */
/* Comments: - From Serpent 1.1.14 (15.11.2010)                              */
/*                                                                           */
/*           - DBRC methodology based on reference:                          */
/*                                                                           */
/*             B. Becker, R. Dagan and G. Lohnert. "Proof and implementation */
/*             of the stochastic formula for ideal gas, energy dependent     */
/*             scattering kernel." Ann. Nucl. Energy, 36 (2009) 470-474.     */
/*                                                                           */
/*           - Routine uses 400*kT as the threshold value for compatibility  */
/*             with MCNP. Dagan defines the limit to 210 eV, which has a     */
/*             noticeable impact on the results.                             */
/*                                                                           */
/*           - 0K datan prosessointia ja pointtereita muutettu 7.2.2013 /    */
/*             2.1.13 (JLe)                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TargetVelocity:"

#define FREEGAS_THRESHOLD  400.0
#define MAX_RESAMPLE 100000000

/*****************************************************************************/

void TargetVelocity(long rea, double E0, double *Vx, double *Vy, double *Vz,
                    double u, double v, double w, double kT, long id)
{
  double awr, V, ar, ycn, r1, z2, rnd1, rnd2, s, z, c, x2, E, max, xs;
  long ptr, nuc, ncol, n;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  ncol = (long)GetPrivateData(ptr, id);

  /* Check if velocity is already sampled for this collision in dopmicroxs */

  if ((z2 = TestValuePair(nuc + NUCLIDE_PREV_COL_Z2, (double)ncol, id)) > 0.0)
    {
      /* Check temperature */

      if (RDB[nuc + NUCLIDE_XS_TEMP] > 0.0)
        Die(FUNCTION_NAME, "Velocity stored for non-zero temperature");

      /* Get cosine */

      c = TestValuePair(nuc + NUCLIDE_PREV_COL_COS, (double)ncol, id);

      /* Sanity check for mu and direction vectors (for NAN's etc.) */

      CheckValue(FUNCTION_NAME, "mu", "", c, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "u", "", u, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "v", "", v, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "w", "", w, -1.01, 1.01);

      /* Rotate direction cosines */

      AziRot(c, &u, &v, &w, id);

      /* Calculate velocity components (z2 on jaettu oikealla ar:llä */
      /* dopmicroxs.c:ssä) */

      V = sqrt(z2);

      *Vx =  u*V;
      *Vy =  v*V;
      *Vz =  w*V;

      /* Check (tohon jotkut järkevämmät rajat) */

      CheckValue(FUNCTION_NAME, "Vx", "", *Vx, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "Vy", "", *Vy, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "Vz", "", *Vz, -INFTY, INFTY);

      /* Exit subroutine */

      return;
    }

  /* Check if velocity is already sampled for this collision in dopmicroxs */

  if ((z2 = TestValuePair(nuc + NUCLIDE_PREV_COL_TV_Z2, (double)ncol, id))
      > 0.0)
    {
      /* Get cosine */

      c = TestValuePair(nuc + NUCLIDE_PREV_COL_TV_COS, (double)ncol, id);

      /* Sanity check for mu and direction vectors (for NAN's etc.) */

      CheckValue(FUNCTION_NAME, "mu", "", c, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "u", "", u, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "v", "", v, -1.01, 1.01);
      CheckValue(FUNCTION_NAME, "w", "", w, -1.01, 1.01);

      /* Rotate direction cosines */

      AziRot(c, &u, &v, &w, id);

      /* Calculate velocity components (z2 on jaettu oikealla ar:llä */
      /* dopmicroxs.c:ssä) */

      V = sqrt(z2);

      *Vx =  u*V;
      *Vy =  v*V;
      *Vz =  w*V;

      /* Check (tohon jotkut järkevämmät rajat) */

      CheckValue(FUNCTION_NAME, "Vx", "", *Vx, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "Vy", "", *Vy, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "Vz", "", *Vz, -INFTY, INFTY);

      /* Exit subroutine */

      return;
    }

  /* Get atomic weight ratio */

  awr = RDB[nuc + NUCLIDE_AWR];

  /* Reset DBRC maximum cross section */

  max = 0.0;

  /* Check if dbrc is used and get maximum cross section */

  if ((long)RDB[DATA_USE_DBRC] == YES)
    if ((ptr = (long)RDB[rea + REACTION_PTR_0K_DATA]) > VALID_PTR)
      if ((E0 > RDB[DATA_DBRC_EMIN]) && (E0 < RDB[DATA_DBRC_EMAX]))
        {
          /* Pointer to majorant */

          ptr = (long)RDB[ptr + REACTION_PTR_TMP_MAJORANT];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Get cross section */

          max = MicroMajorantXS(ptr, E0, id);
        }

  /* Compare to free-gas threshold criteria (MCNP limits) */

  if (((max == 0.0) && (E0 > FREEGAS_THRESHOLD*kT) && (awr > 1.0)) ||
      (kT == 0.0))

    /*

  if (((max == 0.0) && (E0 > 210E-6) && (awr > 1.0)) ||
      (kT == 0.0))
    */

    {
      /* Target velocity insignificant, set components to zero */

      *Vx = 0.0;
      *Vy = 0.0;
      *Vz = 0.0;

      StoreValuePair(nuc + NUCLIDE_PREV_COL_TV_ET, (double)ncol, 0.0, id);

      /* Exit subroutine */

      return;
    }

  /* Additional variables */

  ar = awr/kT;
  ycn = sqrt(E0*ar);

  /* Check */

  CheckValue(FUNCTION_NAME, "ar", "", ar, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "ycn", "", ycn, ZERO, INFTY);

  /* Rejection sampling loop */

  for (n = 0; n < MAX_RESAMPLE; n++)
    {
      /***********************************************************************/

      /***** Rejection by target velocity ************************************/

      /* Algorithm copied from MCNP5 tgtvel subroutine. Samples target  */
      /* energy z2/ar and cosine c between target and neutron velocity. */

      do
        {
          if (RandF(id)*(ycn + 1.12837917) > ycn)
            {
              r1 = RandF(id);
              z2 = -log(r1*RandF(id));
            }
          else
            {
              do
                {
                  rnd1 = RandF(id);
                  rnd2 = RandF(id);

                  r1 = rnd1*rnd1;
                  s = r1 + rnd2*rnd2;
                }
              while (s > 1.0);

              z2 = -r1*log(s)/s - log(RandF(id));
            }

          z = sqrt(z2);
          c = 2.0*RandF(id) - 1.0;

          x2 = ycn*ycn + z2 - 2*ycn*z*c;

          rnd1 = RandF(id)*(ycn + z);
        }
      while (rnd1*rnd1 > x2);

      /* Break loop if no DBRC */

      if (max == 0.0)
        break;

      /***********************************************************************/

      /***** Doppler-broadening rejection correction *************************/

      /* Energy relative to target nuclide */

      E = x2/ar;

      /* Get pointer to zero kelvin data */

      ptr = (long)RDB[rea + REACTION_PTR_0K_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get zero kelvin cross section */

      xs = MicroXS(ptr, E, id);

      /* Add to total count */

      ptr = (long)RDB[DATA_PTR_DBRC_COUNT];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      AddPrivateData(ptr, 1.0, id);

      /* Check value and add error count */

      if (xs > max)
        {
          ptr = (long)RDB[DATA_PTR_DBRC_EXCEED_COUNT];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          AddPrivateData(ptr, 1.0, id);
        }

      /* Rejection sampling */

      if (RandF(id) < xs/max)
        break;

      /***********************************************************************/
    }

  /* Check for infinite loop */

  if (n == MAX_RESAMPLE)
    Warn(FUNCTION_NAME, "Infinite loop: nuc = %s, E = %E, max = %E\n",
         GetText(nuc + NUCLIDE_PTR_NAME), E0, max);

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", c, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", w, -1.01, 1.01);

  StoreValuePair(nuc + NUCLIDE_PREV_COL_TV_ET, (double)ncol, z2/ar*awr, id);
  StoreValuePair(nuc + NUCLIDE_PREV_COL_TV_T, (double)ncol,  kT/KELVIN, id);
  StoreValuePair(nuc + NUCLIDE_PREV_COL_TV_Z2, (double)ncol, z2/ar, id);
  StoreValuePair(nuc + NUCLIDE_PREV_COL_TV_COS, (double)ncol, c, id);

  /* Rotate direction cosines */

  AziRot(c, &u, &v, &w, id);

  /* Calculate velocity components */

  V = sqrt(z2/ar);

  *Vx =  u*V;
  *Vy =  v*V;
  *Vz =  w*V;
}

/*****************************************************************************/
