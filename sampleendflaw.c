/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sampleendflaw.c                                */
/*                                                                           */
/* Created:       2011/02/11 (JLe)                                           */
/* Last modified: 2019/11/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Processes ENDF energy distributions in ACE data              */
/*                                                                           */
/* Comments: - Adopted from sampleenergy.c in Serpent 1.1.12                 */
/*                                                                           */
/*           - Serpent 1:ss‰ voi olla jotain noita interpolaatioita v‰‰rin   */
/*             p‰in, esim. 4/44/61 (tai sitten ei), ja joku toinen viel‰     */
/*                                                                           */
/*           - Interpolation factorit voi olla eri p‰in                      */
/*                                                                           */
/*           - 2.1.22 (15.10.2014) LAW 4/44/61 lis‰ttiin toi energiariippuva */
/*             tyyppivektori                                                 */
/*                                                                           */
/*           - Ton reaktiopointterin v‰litt‰minen t‰nne tuntuu turhalta      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleENDFLaw:"

/*****************************************************************************/

void SampleENDFLaw(long rea, long erg, double Ein, double *Eout, double *mu,
                   long id)
{
  long law, l0, l1, l2, ld1, ld2, ne, nb, nr, type, nd, i, j, k, l, m, n;
  long l4, nuc, np, mt, K1, K2, ptr, np1, np2;
  double rnd1, rnd2, rnd3, rnd4, d0, d1, d2, U, kT, a, b, g, d, c;
  double r, El1, Elk, p1, p2, c1, c2, R, R1, R2, A, A1, A2, T, E0, EE1, EEk;
  double p, x, y, awr, Q ,rd;

  /***************************************************************************/

  /***** Get pointers etc. ***************************************************/

  /* Check incident energy value */

  CheckValue(FUNCTION_NAME, "Ein", "", Ein, ZERO, INFTY);

  /* Get pointer to energy distribution if not given */

  if (erg < VALID_PTR)
    {
      /* Check reaction channel pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get pointer */

      erg = (long)RDB[rea + REACTION_PTR_ERG];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
    }

  /* Reaction pointer is not given for photon emission reactions */

  if (rea > VALID_PTR)
    {
      /* Check boundaries and pointer to energy distribution */

      if ((Ein < RDB[rea + REACTION_EMIN]) || (Ein > RDB[rea + REACTION_EMAX])
          || (erg < VALID_PTR))
        {
          /* No distribution given. No change in energy and cosine. */

          if (Eout != NULL)
            *Eout = Ein;

          if (mu != NULL)
            *mu = 0.0;

          /* Exit */

          return;
        }

      /* Get reaction mt */

      mt = (long)RDB[rea + REACTION_MT];
    }
  else
    mt = -1;

  /* Get pointer to nuclide data */

  nuc = (long)RDB[erg + ERG_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /***************************************************************************/

  /***** Sample distribution law if multiple are laws given ******************/

  /* Get pointer to energy grid */

  if (NextItem(erg) > 0)
    {
      /* Sample random variable */

      rnd1 = RandF(id);

      do
        {
          /* Get pointer to energy grid */

          ptr = (long)RDB[erg + ERG_PTR_EGRID];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Find interval */

          if ((n = GridSearch(ptr, Ein)) < 0)
            {
              /* Get number of energies */

              ne = (long)RDB[ptr + ENERGY_GRID_NE];

              /* Check energy */

              if (Ein < RDB[ptr + ENERGY_GRID_EMIN])
                n = 0;
              else if (Ein > RDB[ptr + ENERGY_GRID_EMAX])
                n = ne - 1;
              else
                Die(FUNCTION_NAME, "wtf?");
            }

          /* Pointer to probability data */

          ptr = (long)RDB[erg + ERG_PTR_PROB];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Compare probability (NOTE: Threshold must be > 0 because some */
          /* nuclides have incorrectly normalised data. */

          /* Tuleeko t‰h‰n interpolaatio? */

          if ((rnd1 = rnd1 - RDB[ptr + n]) < 1E-5)
            break;

          /* Pointer to next law */

          erg = NextItem(erg);

          /* Check pointer */

          if (erg < VALID_PTR)
            Die(FUNCTION_NAME, "Unable to sample law (%s, mt %ld)",
                GetText(nuc + NUCLIDE_PTR_NAME), mt);
        }
      while (1 != 2);
    }

  /***************************************************************************/

  /***** Sample from distributions *******************************************/

  law = (long)RDB[erg + ERG_LAW];

  switch(law)
    {
    default:
      {
        /***** Unsupported law ***********************************************/

        Die(FUNCTION_NAME, "Unsupported ENDF law %ld (%s)",
            law, GetText(nuc + NUCLIDE_PTR_NAME));

        /*********************************************************************/
      }
    case 1:
      {
        /***** Equi-probable bin distribution ********************************/

        /* Get pointer to data */

        l0 = (long)RDB[erg + ERG_PTR_DATA];
        CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

        /* Get pointer to energy grid */

        ptr = (long)RDB[l0++];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get grid factor */

        if ((r = GridFactor(ptr, Ein, id)) < 0)
          {
            /* Avoid compiler warning */

            i = -1;

            /* Compare to minimum */

            if (Ein < RDB[ptr + ENERGY_GRID_EMIN])
              {
                i = 0;
                r = 0.0;
              }
            else
              Die(FUNCTION_NAME, "wtf? (law 1)");
          }
        else
          {
            /* Separate integer and decimal parts */

            i = (long)r;
            r = r - (double)i;
          }

        /* Get number of bins */

        nb = (long)RDB[l0++] + 1;

        /* Check factor */

        CheckValue(FUNCTION_NAME, "r (law 1)", "", r, 0.0, 1.0);

        /* Sample energy interval */

        if (RandF(id) > r)
          l0 = l0 + i*nb;
        else
          l0 = l0 + (i + 1)*nb;

        /* Sample bin */

        rnd1 = RandF(id);
        k = (long)(((double)nb - 1.0)*rnd1);

        /* Get data */

        d0 = RDB[l0 + k];
        d1 = RDB[l0 + k + 1];

        /* Check values */

        CheckValue(FUNCTION_NAME, "d0 (law 1)", "", d0, 0.0, 200.0);
        CheckValue(FUNCTION_NAME, "d1 (law 1)", "", d1, 0.0, 200.0);

        /* Calculate energy */

        *Eout = d0 + (rnd1*((double)nb - 1.0) - (double)k)*(d1 - d0);

        /* Overwrite initial energy to avoid errors in checking. For some  */
        /* reason the data for some isotopes begin at a much higher energy */
        /* than the channel maximum. */

        Ein = d1;

        /* Break case */

        break;

        /*********************************************************************/
      }
    case 2:
      {
        /***** Discrete photon energy ****************************************/

        /* Get pointer to data */

        l0 = (long)RDB[erg + ERG_PTR_DATA];
        CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

        /* Get parameters */

        n = (long)RDB[l0];
        E0 = RDB[l0 + 1];

        /* Check value */

        CheckValue(FUNCTION_NAME, "E0", "", E0, ZERO, INFTY);

        /* Check type */

        if ((n == 0) || (n == 1))
          {
            /* Value is photon energy */

            *Eout = E0;
          }
        else if (n == 2)
          {
            /* Value is binding energy, get awr */

            awr = RDB[nuc + NUCLIDE_AWR];
            CheckValue(FUNCTION_NAME, "awr", "", awr, 0.99, 300.0);

            /* Calculate photon energy */

            *Eout = E0 + awr/(awr + 1.0)*Ein;
          }
        else
          Die(FUNCTION_NAME, "Invalide mode (law 2)");

        /* Break case */

        break;

        /*********************************************************************/
      }
    case 3:
      {
        /***** Level scattering **********************************************/

        /* Get pointer to data */

        l0 = (long)RDB[erg + ERG_PTR_DATA];
        CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

        /* Get parameters */

        d0 = RDB[l0];
        d1 = RDB[l0 + 1];

        /* Check values (isomeric states may have d0 < 0) */

        CheckValue(FUNCTION_NAME, "d0 (law 3)", "", d0, -1.0, 35.0);
        CheckValue(FUNCTION_NAME, "d1 (law 3)", "", d1,  0.3, 1.0);

        if (Ein < d0)
          {
            /* Check if isotope is an isomeric state. */

            if ((long)RDB[nuc + NUCLIDE_I] > 0)
              *Eout = Ein;
            else
              Die(FUNCTION_NAME, "E0 below threshold in level scattering");
          }
        else
          *Eout = d1*(Ein - d0);

        /* Break case */

        break;

        /********************************************************************/
      }
    case 4:
    case 44:
    case 61:
      {
        /*****  4 = Continuous tabular distribution **************************/
        /***** 44 = Kalbach-87 formalism *************************************/
        /***** 61 = like 44 but tabular angular distr. ***********************/

        /* Get pointer to data */

        l0 = (long)RDB[erg + ERG_PTR_DATA];
        CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

        /* Get pointer to energy grid */

        ptr = (long)RDB[l0++];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get number of values */

        ne = (long)RDB[ptr + ENERGY_GRID_NE];
        CheckValue(FUNCTION_NAME, "ne", "", ne, 1, 1000000);

        /* Check limiting values or find bin interval */

        if (Ein <= RDB[ptr + ENERGY_GRID_EMIN])
          {
            i = 0;
            r = 0.0;
          }
        else if (Ein >= RDB[ptr + ENERGY_GRID_EMAX])
          {
            i = ne - 1;
            r = 0.0;
          }
        else
          {
            /* Get grid factor */

            if ((r = GridFactor(ptr, Ein, id)) < 0.0)
              Die(FUNCTION_NAME, "Error in grid search");

            /* Separate integer and decimal parts  */

            i = (long)r;
            r = r - (double)i;
          }

        /* Check interpolation facror */

        CheckValue(FUNCTION_NAME, "r (law 4/44/61)", "", r, 0.0, 1.0);

        /* Get pointer to distribution types */

        l1 = (long)RDB[l0++];
        CheckPointer(FUNCTION_NAME, "(l1)", DATA_ARRAY, l1);

        /* Get type (if > 9, the value also gives the number of discrete */
        /* photon lines) */

        if (RDB[l1 + i] > 9)
          {
            nd = (long)(RDB[l1 + i]/10.0);
            type = (long)RDB[l1 + i] - 10*nd;

            /* Check interpolation type (can be 0 if only discrete lines) */

            CheckValue(FUNCTION_NAME, "type", "", type, 0, 2);

            /* Check second type */

            if ((long)(RDB[l1 + i + 1]/10.0) != nd)
              Die(FUNCTION_NAME, "Mismatch in number of discrete lines");

            /* Check that this really is photon production data */

            if (mt != -1)
              Die(FUNCTION_NAME, "INTT' = %ld, m = %ld\n",
                  (long)RDB[l1 + i], mt);

            /* Pointers to discrete line data */

            ld1 = (long)RDB[l0 + i];
            CheckPointer(FUNCTION_NAME, "(ld1)", DATA_ARRAY, ld1);

            if (i == ne - 1)
              Die(FUNCTION_NAME, "i = ne - 1 (discrete lines)");

            ld2 = (long)RDB[l0 + i + 1];
            CheckPointer(FUNCTION_NAME, "(ld2)", DATA_ARRAY, ld2);

            /* Array sizes */

            np1 = (long)RDB[ld1++];
            np2 = (long)RDB[ld2++];

            CheckValue(FUNCTION_NAME, "np1", "", np1, 1, 1000000);
            CheckValue(FUNCTION_NAME, "np2", "", np2, 1, 1000000);

            /* Remember interpolation factor (in case re-used later) */

            rd = r;
          }
        else
          {
            /* Set interpolation type */

            type = (long)RDB[l1 + i];

            /* Check type */

            CheckValue(FUNCTION_NAME, "type", "", type, 1, 2);

            /* Avoid compiler warning */

            nd = 0;
            ld1 = -1;
            ld2 = -1;
            np1 = -1;
            np2 = -1;
            rd = -1.0;
          }

        /* Check index and interpolation type */

        if ((i < ne - 1) && ((long)RDB[erg + ERG_INTERP] != 1))
          {
            /* Get limiting points for final interpolation */

            l1 = (long)RDB[l0 + i];
            CheckPointer(FUNCTION_NAME, "(l1)", DATA_ARRAY, l1);

            l2 = (long)RDB[l0 + i + 1];
            CheckPointer(FUNCTION_NAME, "(l2)", DATA_ARRAY, l2);

            if ((K1 = (long)RDB[l1++]) < 0)
              Die(FUNCTION_NAME, "K1 < 0");

            if ((K2 = (long)RDB[l2++]) < 0)
              Die(FUNCTION_NAME, "K2 < 0");

            if (nd < K1)
              EE1 = RDB[l1 + nd] + r*(RDB[l2 + nd] - RDB[l1 + nd]);
            else
              EE1 = 0.0;

            /* JLe: T‰h‰n lis‰ttiin nd (22.11.2019 / 2.1.32) */

            if (nd < K2)
              EEk = RDB[l1 + nd + K1 - 1] + r*(RDB[l2 + nd + K2 - 1]
                                               - RDB[l1 + nd + K1 - 1]);
            else
              EEk = 0.0;

            /* Check values */

            CheckValue(FUNCTION_NAME, "EE1 (law 4/44/61)", "", EE1, 0, 500.0);
            CheckValue(FUNCTION_NAME, "EEk (law 4/44/61)", "", EEk, 0, 500.0);

            /* Sample distribution */

            if (RandF(id) > r)
              l = i;
            else
              l = i + 1;

            /* Pointer to distribution */

            l1 = (long)RDB[l0 + l];
            CheckPointer(FUNCTION_NAME, "(l1)", DATA_ARRAY, l1);

            /* Number of energy points */

            ne = (long)RDB[l1++];

            /* Check number of discrete lines */

            CheckValue(FUNCTION_NAME, "nd", "", nd, 0, ne);

            /* Minimum and maximum energies */

            if (nd < ne)
              El1 = RDB[l1 + nd];
            else
              El1 = 0.0;

            Elk = RDB[l1 + ne - 1];

            /* 7.2.2016 / 2.1.25 - T‰‰ seuraava p‰tk‰ on mahdollisesti    */
            /* joku viritelm‰ joka on j‰‰nyt kun s‰mpl‰yksess‰ on ollut   */
            /* jotain muuta vikaa ja se on korjattu n‰in. Erityisesti tuo */
            /* testin toistaminen on v‰h‰n kummallinen juttu. T‰m‰ joka   */
            /* tapauksessa aiheuttaa errorin noiden diskreettien fotoni-  */
            /* viivojen kanssa. J‰‰kˆˆn muuten viel‰ paikalleen. */

            if (nd == 0)
              {
                /* Select the other distribution if values are identical */

                if (El1 == Elk)
                  for (l = i; l < i + 2; l++)
                    {
                      /* Pointer to distribution */

                      l1 = (long)RDB[l0 + l];

                      /* Number of energy points */

                      ne = (long)RDB[l1++];

                      /* Minimum and maximum energies */

                      El1 = RDB[l1];
                      Elk = RDB[l1 + ne - 1];

                      CheckValue(FUNCTION_NAME, "El1", "", El1, 0, 500.0);
                      CheckValue(FUNCTION_NAME, "Elk", "", Elk, 0, 500.0);

                      /* Compare and break */

                      if (El1 != Elk)
                        break;
                    }

                /* Check errors */

                if (El1 == Elk)
                  {
                    for (l = i; l < i + 2; l++)
                      {
                        /* Pointer to distribution */

                        l1 = (long)RDB[l0 + l];

                        /* Number of energy points */

                        ne = (long)RDB[l1++];

                        /* Minimum and maximum energies */

                        El1 = RDB[l1];
                        Elk = RDB[l1 + ne - 1];

                        CheckValue(FUNCTION_NAME, "El1", "", El1, 0, 500.0);
                        CheckValue(FUNCTION_NAME, "Elk", "", Elk, 0, 500.0);

                        /* Compare and break */

                        if (El1 != Elk)
                          break;
                      }

                    Die(FUNCTION_NAME, "El1 = Elk = %1.5E",El1);
                  }
              }

            /* Check energies */

            CheckValue(FUNCTION_NAME, "El1 (law 4/44/61)", "", El1, 0.0, 500.0);
            CheckValue(FUNCTION_NAME, "Elk (law 4/44/61)", "", Elk, 0.0, 500.0);
          }
        else
          {
            /* Pointer to distribution */

            l1 = (long)RDB[l0 + i];
            CheckPointer(FUNCTION_NAME, "(l1)", DATA_ARRAY, l1);

            /* Number of energy points */

            ne = (long)RDB[l1++];

            /* Set interpolation factor to zero */

            r = 0.0;

            /* Avoid compiler warning (this should also cause NAN's in the */
            /* final interpolation if the values are used for some reason) */

            El1 = 0.0;
            Elk = 0.0;
            EE1 = 0.0;
            EEk = 0.0;
          }

        /* Check number of energy points */

        if ((nd < ne) && (ne < 2))
          Die(FUNCTION_NAME, "ne = %ld (< 2) (law 4/44/61)", ne);

        /* Re-sampling loop (NOTE: noi taulukoidut arvot on joillain */
        /* nuklideilla (23000 @ JEFF-3.1.1) sellasia ett‰ ne antaa   */
        /* ihan h‰m‰ri‰ tuloksia (b -1E+6). */

        for (j = 0; j < 100; j++)
          {
            /* Sample secondary bin */

            rnd1 = RandF(id);

            /* Discrete part */

            if (nd > 0)
              {
                /* Avoid some compiler warnings */

                d1 = -1;
                d2 = -1;
                m = -1;

                /* Check if possible */

                if ((rnd1 > RDB[ld1 + 2*np1 + nd - 1]) &&
                    (rnd1 > RDB[ld2 + 2*np2 + nd - 1]))
                  {
                    /* Calculate interpolated probability */

                    R = rd*(RDB[ld2 + 2*np2 + nd - 1] -
                            RDB[ld1 + 2*np1 + nd - 1])
                      + RDB[ld1 + 2*np1 + nd - 1];

                    CheckValue(FUNCTION_NAME, "R", "", R, 0.0, 1.0);

                    /* Skip loop */

                    k = nd;
                  }
                else
                  {
                    /* Loop over lines */

                    for (k = 0; k < nd; k++)
                      {
                        /* Calculate interpolated probability */

                        R = rd*(RDB[ld2 + 2*np2 + k] - RDB[ld1 + 2*np1 + k])
                          + RDB[ld1 + 2*np1 + k];

                        CheckValue(FUNCTION_NAME, "R", "", R, 0.0, 1.0);

                        /* Compare */

                        if (rnd1 < R)
                          {
                            /* Set energy */

                            E0 = RDB[ld1 + k] + rd*(RDB[ld2 + k] -
                                                    RDB[ld1 + k]);

                            /* Negatiivinen arvo voi meinata sidosenergiaa? */
                            /* (kts. LAW 2) */

                            if (E0 < 0.0)
                              Die(FUNCTION_NAME, "Discrete line %E\n", E0);

                            /* Reset interpolation factor */

                            r = 0.0;

                            /* Break loop */

                            break;
                          }
                      }
                  }

                /* Break re-sampling loop if line found */

                if (k < nd)
                  break;

                /* Check if continuous part */

                if (nd == ne)
                  {
                    /* This is to avoid compiler warning */

                    E0 = 0.0;

                    /* Try resampling */

                    if (j < 5)
                      continue;
                    else
                      {
                        /* Print warning */

                        Warn(FUNCTION_NAME,
                             "Unable to sample discrete photon line (%s)",
                             GetText(nuc + NUCLIDE_PTR_NAME));

                        /* Set energy */

                        *Eout = 0.0;

                        /* Exit */

                        return;
                      }
                  }

                /* Adjust random number */

                rnd1 = 1.0 - (1.0 - rnd1)*(1.0 - RDB[l1 + 2*ne + nd])
                  /(1.0 - R);
                CheckValue(FUNCTION_NAME, "rnd1", "", rnd1, 0.0, INFTY);

                /* Pick next point or search */

                if (rnd1 >= 1.0)
                  k = nd;
                else
                  k = SearchArray(&RDB[l1 + 2*ne], rnd1, ne);
              }
            else
              {
                /* Find bin */

                k = SearchArray(&RDB[l1 + 2*ne], rnd1, ne);
              }

            /* Check index */

            if ((k < 0) || (k > ne - 2))
              {
                /* Print warning */

                Warn(FUNCTION_NAME,
                     "Unable to sample secondary energy bin (%s)",
                     GetText(nuc + NUCLIDE_PTR_NAME));

                /* Set energy */

                *Eout = 0.0;

                /* Exit */

                return;
              }
            /*
              Die(FUNCTION_NAME,
                  "Unable to sample secondary energy bin %ld %ld %E %E",
                  k, ne-2, rnd1 - 1, RDB[l1 + 2*ne + k]);
            */

            /* Choose nearest angular distribution for law 61 */

            if (rnd1 - RDB[l1 + 2*ne + k] < RDB[l1 + 2*ne + k + 1] - rnd1)
              m = k;
            else
              m = k + 1;

            /* Sample from first distribution if histogram */

            if (type == 1)
              m = k;

            /* Get values */

            d1 = RDB[l1 + k];
            d2 = RDB[l1 + k + 1];
            p1 = RDB[l1 + ne + k];
            p2 = RDB[l1 + ne + k + 1];
            c1 = RDB[l1 + 2*ne + k];
            c2 = RDB[l1 + 2*ne + k + 1];

            /* NOTE: Noissa p-arvoissa voi olla jotain selv‰sti virheellisi‰ */
            /* lukuja, esim. p1 = -2.14748E+09 26056.60c / dlc200 */

            CheckValue(FUNCTION_NAME, "p1 (law 4/44/61)", "", p1, -INFTY,1E+15);
            CheckValue(FUNCTION_NAME, "p2 (law 4/44/61)", "", p2, -INFTY,1E+15);
            CheckValue(FUNCTION_NAME, "c1 (law 4/44/61)", "", c1, 0.0, 1.0);
            CheckValue(FUNCTION_NAME, "c2 (law 4/44/61)", "", c2, 0.0, 1.0);
            CheckValue(FUNCTION_NAME, "d1 (law 4/44/61)", "", d1, 0.0, 1E+3);
            CheckValue(FUNCTION_NAME, "d2 (law 4/44/61)", "", d2, 0.0, 1E+3);

            /* Avoid compiler warning */

            E0 = 0.0;

            /* Check type */

            if (k < nd)
              {
                /* Discrete photon line */

                Die(FUNCTION_NAME, "Sampled discrete line");
              }
            else if (type == 1)
              {
                /* Histogram */

                if (p1 > 0.0)
                  E0 = d1 + (rnd1 - c1)/p1;
                else if (rnd1 - c1 == 0.0)
                  E0 = d1;
                else
                  {
                    /* Some distributions may have cause errors if */
                    /* rnd1 - c1 is small and p1 zero. */

                    Warn(FUNCTION_NAME, "rnd1 - c1 = %1.5E, p1 = %1.5E",
                         rnd1 - c1, p1);
                    E0 = d1;
                  }
              }
            else if (type == 2)
              {
                /* Linear-linear */

                a = (p2 - p1)/(d2 - d1);
                b = p1*p1 + 2.0*a*(rnd1 - c1);

                if (a == 0.0)
                  E0 = d1 + (rnd1 - c1)/p1;
                else if (b > 0.0)
                  E0 = d1 + (sqrt(b) - p1)/a;
                else if (b > -1E-9)
                  E0 = d1 - p1/a;
                else
                  {
                    /* Re-sample */

                    continue;
                  }
              }
            else
              Die(FUNCTION_NAME, "Invalid interpolation type %ld", type);

            /* Break re-sampling loop */

            break;
          }

        /* Check number of re-samples */

        if (j == 100)
          Die(FUNCTION_NAME, "Unable to sample final energy (mt %ld law %ld)",
              mt, law);

        /* Calculate output value */

        if (r > 0.0)
          {
            /* JLe 10.3.2017 / 2.1.29: T‰m‰ on viel‰ toistaiseksi */
            /* selvitt‰m‰tˆn ongelma (esim. *FENDL-3.0 6013 fotonituotto.) */

            if (Elk == El1)
              {
                /* Print warning */
#ifdef DEBUG
                Warn(FUNCTION_NAME, "Elk = %E, El1 = %E (%s, mt %ld, law %ld)",
                     Elk, El1, GetText(nuc + NUCLIDE_PTR_NAME), mt, law);
#endif
                /* Set value */

                *Eout = E0;
              }
            else
              *Eout = EE1 + (E0 - El1)*(EEk - EE1)/(Elk - El1);
          }
        else
          *Eout = E0;

        /* Check angle (22.9.2009). Some nuclides use law 61 to sample */
        /* fission energy. The variable is not given in Fission(). */

        if (mu == NULL)
          break;

        /* Sample angular distribution if law 44/61 */

        if (law == 44)
          {
            /***** 44 = Kalbach-87 formalism *********************************/

            /* Check discrete lines */

            if (nd > 0)
              Die(FUNCTION_NAME, "law = 44, nd = %ld", nd);

            /* Get parameters */

            R1 = RDB[l1 + 3*ne + k];
            R2 = RDB[l1 + 3*ne + k + 1];
            A1 = RDB[l1 + 4*ne + k];
            A2 = RDB[l1 + 4*ne + k + 1];

            /* Check values */

            CheckValue(FUNCTION_NAME, "R1 (law 44)", "", R1, -1E-3, 1.0005);
            CheckValue(FUNCTION_NAME, "R2 (law 44)", "", R2, -1E-3, 1.0005);
            CheckValue(FUNCTION_NAME, "A1 (law 44)", "", A1, 0.0, 200.0);
            CheckValue(FUNCTION_NAME, "A2 (law 44)", "", A2, 0.0, 200.0);

            /* Check type */

            if (type == 1)
              {
                /* Histogram */

                R = R1;
                A = A1;
              }
            else
              {
                /* Linear-linear */

                R = R1 + (R2 - R1)*(E0 - d1)/(d2 - d1);
                A = A1 + (A2 - A1)*(E0 - d1)/(d2 - d1);
              }

            /* Check R */

            CheckValue(FUNCTION_NAME, "R (law 44)", "", R, -1E-3, 1.0005);

            rnd1 = RandF(id);
            rnd2 = RandF(id);

            if (rnd2 > R)
              {
                T = (2.0*rnd1 - 1.0)*sinh(A);

                /* Calculate mu */

                *mu = log(T + sqrt(T*T + 1.0))/A;
              }
            else
              *mu = log(rnd1*exp(A) + (1.0 - rnd1)*exp(-A))/A;

            /*****************************************************************/
          }
        else if (law == 61)
          {
            /***** 61 = like 44 but tabular angular distr. *******************/

            /* Check discrete lines */

            if (nd > 0)
              Die(FUNCTION_NAME, "law = 61, nd = %ld", nd);

            /* Get pointer to cosine distribution */

            l4 = l1 + 3*ne + m;

            /* Check pointer */

            if (l4 < erg)
              Die(FUNCTION_NAME, "Pointer error");

            /* Get pointer to data */

            if ((l0 = (long)RDB[l4]) < 0)
              {
                /* Isotropic distribution */

                *mu = 1.0 - 2.0*RandF(id);
              }
            else
              {
                /* Get data type */

                type = (long)RDB[l0++];

                /* Get number of tabulated points */

                np = (long)RDB[l0++];

                /* Sample secondary bin */

                rnd1 = RandF(id);

                k = 0;
                while (rnd1 > RDB[l0 + 2*np + k])
                  k++;

                k--;

                /* Check values */

                if ((k < 0) || (k > np - 2) || (rnd1 < RDB[l0 + 2*np + k]) ||
                    (rnd1 > RDB[l0 + 2*np + k + 1]))
                  Die(FUNCTION_NAME, "Unable to sample secondary cosine bin");

                /* Get values */

                d1 = RDB[l0 + k];
                d2 = RDB[l0 + k + 1];
                p1 = RDB[l0 + np + k];
                p2 = RDB[l0 + np + k + 1];
                c1 = RDB[l0 + 2*np + k];
                c2 = RDB[l0 + 2*np + k + 1];

                /* Check values */

                CheckValue(FUNCTION_NAME, "d1 (law 61)", "", d1, -1.0, 1.0);
                CheckValue(FUNCTION_NAME, "d2 (law 61)", "", d2, -1.0, 1.0);
                CheckValue(FUNCTION_NAME, "p1 (law 61)", "", p1, 0.0, 1000.0);
                CheckValue(FUNCTION_NAME, "p2 (law 61)", "", p2, 0.0, 1000.0);
                CheckValue(FUNCTION_NAME, "c1 (law 61)", "", c1, 0.0, 1.0);
                CheckValue(FUNCTION_NAME, "c2 (law 61)", "", c2, 0.0, 1.0);

                /* Check errors and distribution type */

                if ((p1 == 0) && (p2 == 0))
                  *mu = 1.0 - 2.0*RandF(id);
                else if (type == 1)
                  {
                    /* Histogram */

                    *mu = d1 + (rnd1 - c1)/p1;
                  }
                else if (type == 2)
                  {
                    /* Linear-linear */

                    a = (p2 - p1)/(d2 - d1);
                    b = p1*p1 + 2.0*a*(rnd1 - c1);

                    /* Check value */

                    CheckValue(FUNCTION_NAME, "b", "", b, -1E-6, 5E+5);

                    if (a == 0.0)
                      *mu = d1  + (rnd1 - c1)/p1;
                    else if (b > 0)
                      *mu = d1 + (sqrt(b) - p1)/a;
                    else if (b > -1E-6)
                      *mu = d1 - p1/a;
                    else
                      Die(FUNCTION_NAME, "b = %E (2)", b);
                  }
                else
                  Die(FUNCTION_NAME, "Invalid interpolation type %ld", type);
              }

            /*****************************************************************/
          }

        break;

        /*********************************************************************/
      }
    case 7:
    case 9:
    case 11:
      {
        /*****  7 = Maxwell fission spectrum *********************************/
        /*****  9 = evaporation spectrum *************************************/
        /***** 11 = Watt fission spectrum ************************************/

        /* Get pointer to data */

        l0 = (long)RDB[erg + ERG_PTR_DATA];
        CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

        /* Get interpolation type */

        if ((type = (long)RDB[erg + ERG_INTERP]) < 0)
          {
            /* Multiple regions, get number and pointer to data */

            nr = (long)RDB[erg + ERG_NR];
            l1 = (long)RDB[erg + ERG_PTR_INTERP];

            /* Search interval */

            if ((n = SearchArray(&RDB[l1], Ein, nr)) < 0)
              type = 0;
            else
              type = (long)RDB[l1 + nr + n];
          }

        /* Get pointer to energy grid */

        ptr = (long)RDB[l0++];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get interval */

        if ((n = GridSearch(ptr, Ein)) < 0)
          {
            /* Get number of energies */

            ne = (long)RDB[ptr + ENERGY_GRID_NE];
            n = -1;

            /* Check energy */

            if (Ein < RDB[ptr + ENERGY_GRID_EMIN])
              n = 0;
            else if (Ein > RDB[ptr + ENERGY_GRID_EMAX])
              n = ne - 1;
            else
              Die(FUNCTION_NAME, "wtf?");
          }

        /* Read restriction energy if evaporation spectrum */

        if ((law == 9) || (law == 7))
          {
            U = RDB[l0++];

            /* Check Value */

            CheckValue(FUNCTION_NAME, "U (law 9)", "", U, -500.0, 500.0);
          }
        else
          {
            /* This is to avoid compiler warning */

            U = 0;
          }

        /* Check law type */

        if (law == 7)
          {
            /***** Maxwell fission spectrum **********************************/

            /* Set pointers */

            l1 = l0 + n;
            l0 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA] + n;

            /* Interpolate kT */

            kT = ENDFInterp(type, Ein, RDB[l0], RDB[l0 + 1], RDB[l1],
                            RDB[l1 + 1]);

            /* Check value */

            CheckValue(FUNCTION_NAME, "kT (law 7)", "", kT, 0, 10.0);

            /* Check restriction energy */

            if (Ein < U)
              {
                /* print warning */

                Warn(FUNCTION_NAME, "Ein < U (law 7)");

                /* Set energy */

                *Eout = 0.0;

                /* Break loop */

                break;
              }

            /* Sampling loop */

            do
              {
                /* Sample two random numbers */

                do
                  {
                    rnd1 = RandF(id);
                    rnd2 = RandF(id);

                    rnd1 = rnd1*rnd1;
                    rnd2 = rnd2*rnd2;

                    c = rnd1 + rnd2;
                  }
                while (c > 1.0);

                /* Calculate energy */

                *Eout = -kT*((rnd1*log(RandF(id)))/c + log(RandF(id)));
              }
            while (*Eout > Ein - U);

            /* Break case */

            break;

            /*****************************************************************/
          }
        else if (law == 9)
          {
            /***** Evaporation spectrum **************************************/

            /* Set pointers */

            l1 = l0 + n;
            l0 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA] + n;

            /* Interpolate kT */

            kT = ENDFInterp(type, Ein, RDB[l0], RDB[l0 + 1], RDB[l1],
                            RDB[l1 + 1]);

            /* Check value */

            CheckValue(FUNCTION_NAME, "kT (law 9)", "", kT, 0, 10.0);

            /* Check restriction energy */

            if (Ein < U)
              {
                /* print warning */

                Warn(FUNCTION_NAME, "Ein < U (law 9)");

                /* Set energy */

                *Eout = 0.0;

                /* Break loop */

                break;
              }

            /* Check low sampling efficiency */

            if (Ein - U < 0.01*kT)
              {
                *Eout = Ein - U;

                break;
              }

            /* Sample energy (p. 2-44 in MCNP4C manual) */

            do
              *Eout = -kT*log(RandF(id)*RandF(id));
            while (*Eout > Ein - U);

            /* Break case */

            break;

            /*****************************************************************/
          }
        else if (law == 11)
          {
            /***** Watt fission spectrum *************************************/

            /* Set pointers */

            l1 = l0 + n;
            l2 = l1 + (long)RDB[ptr + ENERGY_GRID_NE];
            l0 = (long)RDB[ptr + ENERGY_GRID_PTR_DATA] + n;

            /* Interpolate a and b */

            a = ENDFInterp(type, Ein, RDB[l0], RDB[l0 + 1], RDB[l1],
                           RDB[l1 + 1]);

            b = ENDFInterp(type, Ein, RDB[l0], RDB[l0 + 1], RDB[l2],
                           RDB[l2 + 1]);

            /* Check values */

            CheckValue(FUNCTION_NAME, "a (law 11)", "", a, 0, 2.0);
            CheckValue(FUNCTION_NAME, "b (law 11)", "", b, 0, 10.0);

            /* Calculate a few constants */

            c = 1.0 + a*b/8.0;
            CheckValue(FUNCTION_NAME, "c (law 11)", "", c, 1.0, INFTY);

            g = sqrt(c*c - 1.0) + c;

            /* Sample energy (p. 2-45 in MCNP4C manual) */

            do
              {
                rnd1 = log(RandF(id));
                d = (1.0 - g)*(1.0 - rnd1) - log(RandF(id));
                *Eout = -a*g*rnd1;
              }
            while (d*d > b*(*Eout));

            /* Break case */

            break;

            /*****************************************************************/
          }
      }
    case 66:
      {
        /***** N-body phase space distribution *******************************/

        /* Check reaction pointer */

        CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

        /* Get reaction Q-value */

        Q = RDB[rea + REACTION_Q];

        /* Get pointer to data */

        l0 = (long)RDB[erg + ERG_PTR_DATA];
        CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

        /* Get number of bodies */

        n = (long)RDB[l0];

        /* Get total mass of system (NOTE: Tota pit‰‰ skaalata jollain */
        /* kertoimella jos lakia k‰ytet‰‰n muille kuin neutroneille).  */

        d0 = RDB[l0 + 1];
        CheckValue(FUNCTION_NAME, "d0 (law 66)", "", d0, 2.5, 300.0);

        /* Target awr */

        awr = RDB[nuc + NUCLIDE_AWR];

        /* Calculate maximum energy */

        *Eout = ((d0 - 1.0)/d0)*((awr/(awr + 1.0))*Ein + Q);

        /* Sample random numbers */

        do
          {
            rnd1 = RandF(id);
            rnd1 = rnd1*rnd1;
            rnd2 = RandF(id);
            rnd2 = rnd1 + rnd2*rnd2;
          }
        while (rnd2 > 1.0);

        do
          {
            rnd3 = RandF(id);
            rnd3 = rnd3*rnd3;
            rnd4 = RandF(id);
            rnd4 = rnd3 + rnd4*rnd4;
          }
        while (rnd4 > 1.0);

        /* Avoid compiler warning */

        p = 0.0;

        /* Get value of p */

        if (n == 3)
          p = RandF(id);
        else if (n == 4)
          p = RandF(id)*RandF(id);
        else if (n == 5)
          p = RandF(id)*RandF(id)*RandF(id)*RandF(id);
        else
          Die(FUNCTION_NAME, "Invalid n-value");

        /* Get x and y */

        x = -rnd1*log(rnd2)/rnd2 - log(RandF(id));
        y = -rnd3*log(rnd4)/rnd4 - log(p);

        /* Final emission energy */

        *Eout = (x/(x + y))*(*Eout);

        /* Mu is sampled isotropically */

        *mu = 1.0 - 2.0*RandF(id);

        /* Break case */

        break;

        /********************************************************************/
      }
    }

  /***************************************************************************/

#ifdef DEBUG

  /* Check that sampled energy is reasonable (energy may be zero */
  /* in photon production, limits are checked later) */

  if (Eout != NULL)
    {
      /* This checks for NAN and INF */

      if (!((*Eout > -INFTY) && (*Eout < INFTY)))
        Die(FUNCTION_NAME, "Ein = %E, Eout = %E (%s, mt %ld, law %ld)", Ein,
            *Eout, GetText(nuc + NUCLIDE_PTR_NAME), mt, law);

      /* This checks for reasomable limits */

      if ((!((*Eout > 1E-18) && (*Eout < 1000.0))) && (mt > 0))
        Warn(FUNCTION_NAME, "Ein = %E, Eout = %E (%s, mt %ld, law %ld)", Ein,
             *Eout, GetText(nuc + NUCLIDE_PTR_NAME), mt, law);
    }

  /* Check mu */

  if (mu != NULL)
    {
      /* Check value */

      if (fabs(*mu) - 1.0 > 1E-4)
        Die(FUNCTION_NAME, "mu = %f (|mu| - 1 = %E)", *mu, fabs(*mu) - 1.0);

      /* Check dubious values */

      if ((*mu == -1.0) || (*mu == 0.0) || (*mu == 1.0))
        Die(FUNCTION_NAME, "Dubious value mu = %1.4f sampled (%s, law %ld)",
             *mu, GetText(nuc + NUCLIDE_PTR_NAME), law);
    }

#endif

  /* Adjust mu */

  if (mu != NULL)
    {
      if (*mu > 1.0)
        *mu = 1.0;
      else if (*mu < -1.0)
        *mu = -1.0;
    }
}

/*****************************************************************************/
