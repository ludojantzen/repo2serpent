/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sabscattering.c                                */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2019/11/28 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Handles S(a,b) scattering laws                               */
/*                                                                           */
/* Comments: - Testattu MCNP6:n diskreeteillä kirjastoilla 2.10.2016/2.1.28. */
/*                                                                           */
/*           - Continuous out-going energy distribution described in OpenMC  */
/*             manual sections 5.11.4.3 and 5.8.1.2.2. url:                  */
/*                                                                           */
/*             docs.openmc.org/en/stable/methods/neutron_physics.html        */
/*                                                                           */
/*           - NOTE: Jos tähän tehdään muutoksia niin samat pitää tehdä myös */
/*                   otfsabscattering.c:hen.                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SabScattering:"

/*****************************************************************************/

void SabScattering(long rea, double *E, double *u, double *v, double *w,
                   long id)
{
  double mu, E0, Ei0, r, d1, d2, a, b, rnd, p1, p2, c1, c2, f;
  long erg, law, l0, l1, l2, nc2, ctype, ne, ne2, n, i, j, k, l, cyc;
  long ptr;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to energy distribution */

  erg = (long)RDB[rea + REACTION_PTR_ERG];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Initial energy */

  E0 = *E;

  /* Avoid compiler warning */

  mu = 1.0;

  /* Get pointer to data */

  l0 = (long)RDB[erg + ERG_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

  /* Get scattering mode (elastic or inelastic) */

  law = (long)RDB[erg + ERG_LAW];

  /* Get cosine distribution type */

  ctype = (long)RDB[l0++] + 1;

  /* Get pointer to incident energy grid (recycle pointer) */

  erg = (long)RDB[l0++];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Check for distribution type */

  if ((law == 1004) && ((long)RDB[l0 + 2] == 2))
    {
      /***********************************************************************/

      /***** Continous inelastic *********************************************/

      /* Get pointer to grid data */

      ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get interpolation factor */

      if ((r = GridFactor(erg, E0, id)) < 0.0)
        {
          /* Get number of energies */

          ne = (long)RDB[erg + ENERGY_GRID_NE];
          i = -1;

          /* Check energy */

          if (E0 <= RDB[erg + ENERGY_GRID_EMIN])
            {
              i = 0;
              r = 0.0;
            }
          else if (E0 >= RDB[erg + ENERGY_GRID_EMAX])
            {
              i = ne - 2;
              r = 1.0;
            }
        }
      else
        {
          /* Separate integer and decimal parts */

          i = (long)r;
          r = r - (double)i;
        }

      /* Check interpolation factor and energy */

      CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);
      CheckValue(FUNCTION_NAME, "E0", "", E0, 0.99999*RDB[ptr + i],
                 1.00001*RDB[ptr + i  +1]);

      /* Skip number of secondary energy bins */

      l0++;

      /* Get number of secondary cosine bins */

      nc2 = (long)RDB[l0++];

      /* Skip type */

      l0++;

      /* Pick closest distribution */

      if (r < 0.5)
        {
          Ei0 = RDB[ptr + i];
          l = (long)RDB[l0 + i];
        }
      else
        {
          Ei0 = RDB[ptr + i + 1];
          l = (long)RDB[l0 + i + 1];
        }

      /* Check values */

      CheckValue(FUNCTION_NAME, "Ei0", "", Ei0, ZERO, INFTY);
      CheckPointer(FUNCTION_NAME, "(l)", DATA_ARRAY, l);

      /* Number of energies */

      ne = (long)RDB[l++];
      CheckValue(FUNCTION_NAME, "ne", "", ne, 2, 2000);

      /* Avoid compiler warning */

      d1 = -1;
      d2 = -1;

      /* Re-sampling loop */

      for (cyc = 0; cyc < 100; cyc++)
        {
          /* Sample secondary bins*/

          rnd = RandF(id);

          j = 0;
          while (rnd > RDB[l + 2*ne + j])
            j++;

          j--;

          /* Check index */

          if ((j < 0) || (j > ne - 2))
            continue;

          /* Get values */

          d1 = RDB[l + j];
          d2 = RDB[l + j + 1];
          p1 = RDB[l + ne + j];
          p2 = RDB[l + ne + j + 1];
          c1 = RDB[l + 2*ne + j];
          c2 = RDB[l + 2*ne + j + 1];

          /* Check values */

          CheckValue(FUNCTION_NAME, "d1", "", d1, 1E-15, 1E-5);
          CheckValue(FUNCTION_NAME, "d2", "", d2, 1E-15, 1E-5);
          CheckValue(FUNCTION_NAME, "p1", "", p1, -1E-15, INFTY);
          CheckValue(FUNCTION_NAME, "p2", "", p2, -1E-15, INFTY);
          CheckValue(FUNCTION_NAME, "c1", "", c1, 0.0, rnd);
          CheckValue(FUNCTION_NAME, "c2", "", c2, rnd, 1.0);

          /* Linear-linear interpolation */

          a = (p2 - p1)/(d2 - d1);
          b = p1*p1 + 2.0*a*(rnd - c1);

          /* Check values */

          CheckValue(FUNCTION_NAME, "a", "", a, -INFTY, INFTY);
          CheckValue(FUNCTION_NAME, "b", "", b, -INFTY, INFTY);

          if (a == 0.0)
            *E = d1  + (rnd - c1)/p1;
          else if (b > 0.0)
            *E = d1 + (sqrt(b) - p1)/a;
          else
            *E = d1 - p1/a;

          /* Interpolation between bins */

          if (*E < 0.5*Ei0)
            *E = *E*(2.0*E0/Ei0 - 1.0);
          else
            *E = *E + E0 - Ei0;

          /* Calculate interpolation factor mu sampling */

          f = (rnd - c1)/(c2 - c1);
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);

          /* Pointers to distribution */

          l1 = (long)RDB[l + 3*ne + j];
          CheckPointer(FUNCTION_NAME, "(l1)", DATA_ARRAY, l1);

          l2 = (long)RDB[l + 3*ne + j + 1];
          CheckPointer(FUNCTION_NAME, "(l2)", DATA_ARRAY, l2);

          /* Sample bin */

          k = (long)(((double)nc2)*RandF(id));

          /* k value */

          mu = RDB[l1 + k] + f*(RDB[l2 + k] - RDB[l1 + k]);

          /* k - 1 value */

          if (k == 0)
            d1 = -1.0 - (mu + 1.0);
          else
            d1 = RDB[l1 + k - 1] + f*(RDB[l2 + k - 1] - RDB[l1 + k - 1]);

          /* k + 1 value */

          if (k == nc2 - 1)
            d2 = 1.0 + (1.0 - mu);
          else
            d2 = RDB[l1 + k + 1] + f*(RDB[l2 + k + 1] - RDB[l1 + k + 1]);

          /* Smear */

          if ((mu - d1) < (d2 - mu))
            mu = mu + (mu - d1)*(RandF(id) - 0.5);
          else
            mu = mu + (d2 - mu)*(RandF(id) - 0.5);

          /* Check */

          if ((mu >= -1.0) && (mu <= 1.0))
            break;
        }

      /* Check */

      if (cyc == 100)
        Die(FUNCTION_NAME, "Unable to sample energy or cosine");
      else if (cyc > 5)
        Warn(FUNCTION_NAME, "Energy and angle re-sampled %ld times", cyc);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Discrete distributions ******************************************/

      /* Sample energy and/or scattering cosine */

      if (ctype == 5)
        {
          /*******************************************************************/

          /***** Exact treatment of elastic scattering ***********************/

          /* Check scattering mode */

          if (law != 1002)
            Die(FUNCTION_NAME, "Invalid scattering mode");

          /* Get index to energy grid */

          if ((i = GridSearch(erg, E0)) < 0)
            {
              /* Energy below or above grid, get grid size */

              ne = (long)RDB[erg + ENERGY_GRID_NE];

              /* Avoid compiler warning */

              i = -1;

              /* Check incident energy */

              if (E0 < RDB[erg + ENERGY_GRID_EMIN])
                i = 0;
              else if (E0 > RDB[erg + ENERGY_GRID_EMAX])
                i = ne - 2;
              else
                Die(FUNCTION_NAME, "law 1002: i = %E", i);
            }

          /* Pointers to data */

          l1 = l0;
          l2 = l0 + i + 1;

          /* Sample */

          d1 = RandF(id)*RDB[l2];

          /* Search */

          while (l2 != l1 + 1)
            {
              n = (long)((double)(l1 + l2)/2.0);

              if (d1 < RDB[n])
                l2 = n;
              else
                l1 = n;
            }

          /* Pointer to incident energy grid data (recycle pointer) */

          erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

          /* Pointer to value */

          erg = erg + l2 - l0 - 1;

          /* Calculate mu */

          mu = 1.0 - 2.0*RDB[erg]/E0;

          /*******************************************************************/
        }
      else
        {
          /*******************************************************************/

          /***** Other elastic or inelastic modes ****************************/

          /* Get interpolation factor */

          if ((r = GridFactor(erg, E0, id)) < 0.0)
            {
              /* Get number of energies */

              ne = (long)RDB[erg + ENERGY_GRID_NE];
              i = -1;

              /* Check energy */

              if (E0 <= RDB[erg + ENERGY_GRID_EMIN])
                {
                  i = 0;
                  r = 0.0;
                }
              else if (E0 >= RDB[erg + ENERGY_GRID_EMAX])
                {
                  i = ne - 2;
                  r = 1.0;
                }
            }
          else
            {
              /* Separate integer and decimal parts */

              i = (long)r;
              r = r - (double)i;
            }

          /* Check interpolation factor */

          CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);

          /* Avoid compiler warning */

          nc2 = -1;
          l1 = -1;
          l2 = -1;

          /* Check scattering mode */

          if (law == 1004)
            {
              /* Inelastic mode, get number of secondary energies and cosines */

              ne2 = (long)RDB[l0++];
              nc2 = (long)RDB[l0++];

              /* Get secondary energy type */

              n = (long)RDB[l0++];
              CheckValue(FUNCTION_NAME, "n", "", n, 0, 1);

              /* Get index to energy bin */

              if (n == 0)
                k = (long)(RandF(id)*(double)ne2);
              else if ((a = RandF(id)*((double)ne2 - 3.0)) > 1.0)
                k = (long)a + 1;
              else if (a > 0.6)
                k = ne2 - 2;
              else if (a > 0.5)
                k = ne2 - 1;
              else if (a > 0.1)
                k = 1;
              else
                k = 0;

              /* Check value */

              CheckValue(FUNCTION_NAME, "k", " (law 1004)", k, 0, ne2 - 1);

              /* Pointers to distribution */

              l1 = (long)RDB[l0 + i] + k*(nc2 + 1);
              l2 = (long)RDB[l0 + i + 1] + k*(nc2 + 1);

              /* Get energy values */

              d1 = RDB[l1];
              d2 = RDB[l2];

              /* Interpolate energy value */

              *E = d1 + r*(d2 - d1);

              /* Update pointers (cosines are listed after each energy) */

              l1++;
              l2++;
            }
          else if (law == 1002)
            {
              /* NOTE: tää on kohtuullisen harvinainen sirontalaki, johon */
              /* törmää esim. h/zr ja zr/h -kirjastoissa. Ei ole kunnolla */
              /* testattu. */

              /* Elastic, get number of cosines */

              nc2 = (long)RDB[l0++];

              /* Get pointers */

              l1 = l0 + i*nc2;
              l2 = l0 + (i + 1)*nc2;
            }
          else
            Die(FUNCTION_NAME, "Invalid scattering mode");

          /* Check cosine distribution type */

          if (ctype == -2)
            {
              /* Re-sampling loop */

              for (cyc = 0; cyc < 100; cyc++)
                {
                  /* Sample bin */

                  k = (long)(((double)nc2)*RandF(id));

                  /* k value */

                  mu = RDB[l1 + k] + r*(RDB[l2 + k] - RDB[l1 + k]);

                  /* k - 1 value */

                  if (k == 0)
                    d1 = -1.0 - (mu + 1.0);
                  else
                    d1 = RDB[l1 + k - 1] + r*(RDB[l2 + k - 1] -
                                              RDB[l1 + k - 1]);

                  /* k + 1 value */

                  if (k == nc2 - 1)
                    d2 = 1.0 + (1.0 - mu);
                  else
                    d2 = RDB[l1 + k + 1] + r*(RDB[l2 + k + 1] -
                                              RDB[l1 + k + 1]);

                  /* Smear */

                  if ((mu - d1) < (d2 - mu))
                    mu = mu + (mu - d1)*(RandF(id) - 0.5);
                  else
                    mu = mu + (d2 - mu)*(RandF(id) - 0.5);

                  /* Check */

                  if ((mu >= -1.0) && (mu <= 1.0))
                    break;
                }

              /* Check */

              if (cyc == 100)
                Die(FUNCTION_NAME, "Unable to sample cosine");
              else if (cyc > 5)
                Warn(FUNCTION_NAME, "Cosine re-sampled %ld times", cyc);
            }
          else if (ctype == 4)
            {
              /* Sample from discrete cosines */

              l = (long)(((double)nc2)*RandF(id));

              /* Get cosine values */

              d1 = RDB[l1 + l];
              d2 = RDB[l2 + l];

              /* Check values */

              CheckValue(FUNCTION_NAME, "d1 (law 1004)", "", d1, -1.0, 1.0);
              CheckValue(FUNCTION_NAME, "d2 (law 1004)", "", d2, -1.0, 1.0);

              /* Interpolate mu value */

              mu = d1 + r*(d2 - d1);
            }
          else
            Die(FUNCTION_NAME,
                "Invalid or unsupported cosine distribution mode");

          /*******************************************************************/
        }

      /***********************************************************************/
    }

  /* Check energy and cosines */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r", "", *u**u+*v**v+*w**w - 1.0, -1E-5, 1E-5);

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.01, 1.01);

  /* Rotate direction cosines around a random azimuthal angle */

  AziRot(mu, u, v, w, id);
}

/*****************************************************************************/
