/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplemu.c                                     */
/*                                                                           */
/* Created:       2010/01/19 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Samples scattering cosine from ENDF distributions            */
/*                                                                           */
/* Comments: - Adopted from angulardistribution.c in Serpent 1.1.12          */
/*           - Serpent 1.1.12 has sections of code that are commented out    */
/*             and changed, but the routine should be OK.                    */
/*                                                                           */
/*           - Dubious value check may terminate calculation (some limit     */
/*             hit and mu = 0.0?)                                            */
/*                                                                           */
/*           - Added new mode (16/01/30, 2.1.25) for sampling direction of   */
/*             emitted photon in coupled neutro/gamma simulation. Pointer    */
/*             to angular distribution is then given directly, or null       */
/*             if isotropic. Reaction pointer must be null.                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleMu:"

/*****************************************************************************/

double SampleMu(long rea, long ang, double E, double LMS[7],
                    double *pdf, long id)
{
  double mu, r, rnd, d1, d2, p1, p2, c1, c2, a, b, Emin, Emax;
  long sens, type, i, j, k, ne, np, nb, erg, l0;
  long ptrlow, ptrhigh, nplow, nphigh;
  double pdflow, pdfhigh;

  /* Check for mode (reaction pointer is given for neutron reactions */
  /* and null for photon production). */

  if (rea > VALID_PTR)
    {

#ifdef DEBUG

      /* Pointer must be null */

      if (ang > VALID_PTR)
        Die(FUNCTION_NAME, "ang = %ld", ang);

#endif

      /* Get pointer to angular distribution. Use isotropic distribution */
      /* for hydrogen to speed up calculation */

      if ((RDB[rea + REACTION_AWR] < 1.0) && (E < 1E-3))
        ang = -1;
      else
        ang = (long)RDB[rea + REACTION_PTR_ANG];
    }

  /* Check angular distribution */

  if (ang < VALID_PTR)
    {
      /* No distribution. Sample isotropic */

      mu = 1.0 - 2.0*RandF(id);

      /* Adjust values if > +- 1.0 */

      if (mu > 1.0)
        mu = 1.0;
      else if (mu < -1.0)
        mu = -1.0;

      /* Return mu */

      return mu;
    }

  /* Get pointer to energy grid */

  erg = (long)RDB[ang + ANG_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Get number of energy points */

  ne = (long)RDB[erg + ENERGY_GRID_NE];
  CheckValue(FUNCTION_NAME, "ne", "", ne, 2, 1000000);

  /* Get minimum and maximum energy */

  Emin = RDB[erg + ENERGY_GRID_EMIN];
  Emax = RDB[erg + ENERGY_GRID_EMAX];

  /* Get interpolation factor */

  if ((r = GridFactor(erg, E, id)) < 0.0)
    {
      /* Avoid compiler warning */

      i = 0;
      r = 0.0;

      /* Check if energy is above or below limits */

      if (E > Emax)
        {
          i = ne - 2;
          r = 1.0;
        }
      else if (E < Emin)
        {
          i = 0;
          r = 0.0;
        }
      else
        Die(FUNCTION_NAME, "wtf?!? (%E %E %E)", Emin, E, Emax);
    }
  else
    {
      /* Check incident energy */

      CheckValue(FUNCTION_NAME, "E", "", E, Emin, Emax);

      /* Get bin index and adjust factor */

      i = (long)r;
      r = r - (double)i;
    }

  /* Check values */

  CheckValue(FUNCTION_NAME, "i", "", i, 0, ne - 1);
  CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);

  /* Avoid compiler warnings */

  mu = 0.0;

  /* Check type */

  if ((type = (long)RDB[ang + ANG_TYPE]) == ANG_TYPE_EQUIBIN)
    {
      /***********************************************************************/

      /***** 32 equi-probable cosine bin distribution ************************/

      if (((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
          ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM))
        Die (FUNCTION_NAME, "Legendre sensitivity error\n"
             "32 equi-probable cosine bin distribution mt = %ld, E = %E, nucZAI = %ld ",
             (long)RDB[rea + REACTION_MT], E,
             (long)RDB[((long)RDB[rea + REACTION_PTR_NUCLIDE] + NUCLIDE_ZAI)]);

      /* Get number of bins */

      nb = (long)RDB[ang + ANG_BINS] + 1;

      /* Get pointer to data array */

      l0 = (long)RDB[ang + ANG_PTR_D0];
      CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

      /* Sample cosine interval */

      if (RandF(id) > r)
        l0 = l0 + i*nb;
      else
        l0 = l0 + (i + 1)*nb;

     /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

      /* Sample bin */

      rnd = RandF(id);

      k = (long)(((double)nb - 1.0)*rnd);

      /* Get data */

      d1 = RDB[l0 + k];
      d2 = RDB[l0 + k + 1];

      /* Check values */

      CheckValue(FUNCTION_NAME, "d1", "", d1, -1.0, 1.0);
      CheckValue(FUNCTION_NAME, "d2", "", d2, -1.0, 1.0);

      /* Get value */

      mu = d1 + (rnd*((double)nb - 1.0) - (double)k)*(d2 - d1);

      /***********************************************************************/
    }
  else if (type == ANG_TYPE_TABULAR)
    {
      /***********************************************************************/

      /***** Tabular angular distribution ************************************/

      /* Pointer to data */

      l0 = (long)RDB[ang + ANG_PTR_D0];
      CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

      if (((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
          ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM) &&
          (LMS != NULL))

        {
          /* pointers to the "np" value, i.e. before first cosine value (d1) */

          ptrlow = (long)RDB[l0 + i];
          ptrhigh = (long)RDB[l0 + i + 1];

          /* Number of cosine points in the low and high energy grids */

          nplow = (long)RDB[ptrlow];
          nphigh = (long)RDB[ptrhigh];

          /* Interpolate f^x_n (E) between low and high energy grids  */
          /* This assumes that the data in the grid is pre-calculated */
          /* Legendre moment factors (stored in processmudistributions.c) */

          LMS[0] = (1-r)*RDB[ptrlow + 3*nplow +1] + r*RDB[ptrhigh + 3*nphigh +1];
          LMS[1] = (1-r)*RDB[ptrlow + 3*nplow +2] + r*RDB[ptrhigh + 3*nphigh +2];
          LMS[2] = (1-r)*RDB[ptrlow + 3*nplow +3] + r*RDB[ptrhigh + 3*nphigh +3];
          LMS[3] = (1-r)*RDB[ptrlow + 3*nplow +4] + r*RDB[ptrhigh + 3*nphigh +4];
          LMS[4] = (1-r)*RDB[ptrlow + 3*nplow +5] + r*RDB[ptrhigh + 3*nphigh +5];
          LMS[5] = (1-r)*RDB[ptrlow + 3*nplow +6] + r*RDB[ptrhigh + 3*nphigh +6];
          LMS[6] = (1-r)*RDB[ptrlow + 3*nplow +7] + r*RDB[ptrhigh + 3*nphigh +7];
        }
      else
        {
          ptrlow = -1;
          ptrhigh = -1;
        }

      /* Sample distribution */

      if (RandF(id) > r)
        l0 = (long)RDB[l0 + i];
      else
        l0 = (long)RDB[l0 + i + 1];

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

      /* Get number of points */

      np = (long)RDB[l0++];

      /* Check value */

      CheckValue(FUNCTION_NAME, "np", "", np, 2, 1000000);

      /* Re-sampling loop (give a second chance) */

      j = 0;

      do
        {
          /* Sample secondary bin */

          rnd = RandF(id);

          k = 0;
          while (rnd > RDB[l0 + 2*np + k])
            k++;

          k--;

          /* Check values */

          if ((k < 0) || (k > np - 2) || (rnd < RDB[l0 + 2*np + k]) ||
              (rnd > RDB[l0 + 2*np + k + 1]))
            Die(FUNCTION_NAME, "Unable to sample secondary cosine bin");

          /* Get values */

          d1 = RDB[l0 + k];
          d2 = RDB[l0 + k + 1];
          p1 = RDB[l0 + np + k];
          p2 = RDB[l0 + np + k + 1];
          c1 = RDB[l0 + 2*np + k];
          c2 = RDB[l0 + 2*np + k + 1];

          /* Check values */

          CheckValue(FUNCTION_NAME, "d1", "", d1, -1.0, 1.0);
          CheckValue(FUNCTION_NAME, "d2", "", d2, -1.0, 1.0);
          CheckValue(FUNCTION_NAME, "p1", "", p1, 0.0, 380.0);
          CheckValue(FUNCTION_NAME, "p2", "", p2, 0.0, 380.0);
          CheckValue(FUNCTION_NAME, "c1", "", c1, 0.0, 1.0);
          CheckValue(FUNCTION_NAME, "c2", "", c2, 0.0, 1.0);

          /* Check distribution type */

          if ((long)RDB[ang + ANG_INTT] == 1)
            {
              /* Check Legendre sensitivity calculation mode */

              if (((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
                  ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM))
                Die(FUNCTION_NAME, "Legendre sensitivity error\n"
                    "Histogram distributions do not currently work with "
                    "Legendre sensitivity calculations");

              /* Histogram */

              mu = d1 + (rnd - c1)/p1;
            }
          else if ((long)RDB[ang + ANG_INTT] == 2)
            {
              /* Linear-linear */

              a = (p2 - p1)/(d2 - d1);
              b = p1*p1 + 2.0*a*(rnd - c1);

              /* Check values */

              CheckValue(FUNCTION_NAME, "a", "", a, -INFTY, INFTY);
              CheckValue(FUNCTION_NAME, "b", "", b, -1E-3, 1.5E+5);

                if (a == 0.0)
                mu = d1  + (rnd - c1)/p1;
              else if (b > 0.0)
                mu = d1 + (sqrt(b) - p1)/a;
              else
                mu = d1 - p1/a;
            }
          else
            Die(FUNCTION_NAME, "Invalid interpolation type");

          /* Increase counter */

          j++;
        }
      while ((fabs(mu) > 1.0) && (j < 5));

      if (j == 5)
        Die(FUNCTION_NAME, "j = 5");

      if (((sens = (long)RDB[DATA_PTR_SENS0]) > VALID_PTR) &&
          ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_SCATT_MOM) &&
          (pdf != NULL))
        {
          /* Interpolate f^x (mu|E) between low and high energy grids */

          /* f^x (mu|Elow) */

          pdflow = InterpolateMuPDF(ptrlow, mu);

          /* f^x (mu|Ehigh) */

          pdfhigh = InterpolateMuPDF(ptrhigh, mu);

          /* Interpolate */

          *pdf = (1-r)*pdflow + r*pdfhigh;
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid distribtion type %ld", type);

  /* Check value */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.0, 1.0);

#ifdef DEBUG

  /* Check dubious values */

  if ((mu == -1.0) || (mu == 0.0) || (mu == 1.0))
    {
      if (rea > VALID_PTR)
        Warn (FUNCTION_NAME, "Dubious value mu = %E (mt = %ld, E = %E)",
              mu, (long)RDB[rea + REACTION_MT], E);
      else
        Warn (FUNCTION_NAME, "Dubious value mu = %E (photon, E = %E)",
              mu, E);
    }

#endif

  /* Adjust values if > +- 1.0 */

  if (mu > 1.0)
    mu = 1.0;
  else if (mu < -1.0)
    mu = -1.0;

  /* Exit */

  return mu;
}

/*****************************************************************************/
