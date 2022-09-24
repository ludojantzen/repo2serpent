/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normcoef.c                                     */
/*                                                                           */
/* Created:       2011/05/05 (JLe)                                           */
/* Last modified: 2016/09/29 (JLe)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Calculates normalization coefficients for reaction rates     */
/*                                                                           */
/* Comments: - Coefficient can be zero in burnup calculation                 */
/*                                                                           */
/*           - Spontaaniin fissioon normeeraus menee eri tavalla kuin        */
/*             ykkösessä                                                     */
/*                                                                           */
/*           - Noi kertoimet pitää laittaa negatiivisiksi kaikista muista    */
/*             paitsi siitä mitä käytetään                                   */
/*                                                                           */
/*           - Ton viimeksi luetun arvon voisi poimia muistiin ja hakea      */
/*             sieltä                                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormCoef:"

/*****************************************************************************/

double NormCoef(long type)
{
  long loc0, ptr, i, mat;
  double div, norm, dh, fiss, fissE, capt, leak, flx, src, sf, fmass, nsf;
  double val, heat, cut, dt;

  /* Check batch counter */

  if ((long)RDB[DATA_BATCH_COUNT] != (long)RDB[DATA_BATCH_INTERVAL])
    Die(FUNCTION_NAME, "Mismatch in batch count");

  /* Get time interval */

  if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
    dt = RDB[DATA_TIME_CUT_TMAX] - RDB[DATA_TIME_CUT_TMIN];
  else
    dt = 1.0;

  /* Check multi-particle transport mode (NOTE: tää olettaa nyt että  */
  /* neutronit on primäärilähde jonka mukaan myös normeeraus tehdään) */

  if ((long)RDB[DATA_MULTI_PARTICLE_TRANSPORT] == YES)
    type = PARTICLE_TYPE_NEUTRON;

  /* Check previous */

  if ((type == PARTICLE_TYPE_NEUTRON) && (RDB[DATA_NORM_COEF_N] > 0.0))
    return RDB[DATA_NORM_COEF_N];
  else if ((type == PARTICLE_TYPE_GAMMA) && (RDB[DATA_NORM_COEF_G] > 0.0))
    return RDB[DATA_NORM_COEF_G];

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Check track plotter mode */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    return 1.0;

  /* Avoid compiler warning */

  val = -1.0;
  div = -1.0;
  cut = -1.0;
  norm = -1.0;

  /* Check radioactive decay source */

  if ((long)RDB[DATA_USE_DECAY_SRC] == YES)
    {
      /***********************************************************************/

      /***** Radioactive decay source ****************************************/

      /* Get decay source rate used for normalization */

      if ((val = RDB[DATA_NORM_DECAY_SRC_RATE]) < ZERO)
        Die(FUNCTION_NAME, "Zero source rate for normalization");

      /* Check transport mode (either should be OK) */

      if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
        {
          /* Normalize to source rate */

          ptr = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            div = BufVal(ptr, 0);

          /* Check divisor */

          if (div < ZERO)
            Die(FUNCTION_NAME, "Zero neutron source rate sampled");
        }
      else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
        {
          /* Normalize to source rate */

          ptr = (long)RDB[RES_TOT_PHOTON_SRCRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
            div = BufVal(ptr, 0);

          /* Check divisor */

          if (div < ZERO)
            Die(FUNCTION_NAME, "Zero photon source rate sampled");
        }
      else
        Die(FUNCTION_NAME, "Invalid type");

      /***********************************************************************/
    }

  /* Check type */

  else if (type == PARTICLE_TYPE_NEUTRON)
    {
      /***********************************************************************/

      /***** Calculate reaction rates for normalization **********************/

      /* Check transport mode */

      if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
        return 0.0;

      /* Get pointer to normalization */

      loc0 = (long)RDB[DATA_PTR_NORM];

      /* Get pointer to material (tätä pointteria käytetään nyt tarkistamaan */
      /* että onko normalisointi kiinnitetty johonkin tiettyyn materiaaliin. */
      /* myöhemmin noi arvot luetaan aina sieltä) */

      if (loc0 > VALID_PTR)
        mat = (long)RDB[loc0 + NORM_PTR_MAT];
      else
        mat = -1;

      /* Avoid compiler warning */

      i = -1;
      fmass = -INFTY;
      sf = -1.0;
      dh = -1.0;

      /* Check normalization type */

      if ((long)RDB[DATA_NORM_BURN] == BURN_NORM_ALL)
        {
          fmass = RDB[DATA_INI_FMASS];
          sf = RDB[DATA_TOT_SFRATE];
          dh = RDB[DATA_TOT_DECAY_HEAT];
          i = 0;
        }
      else if ((long)RDB[DATA_NORM_BURN] == BURN_NORM_BURN)
        {
          fmass = RDB[DATA_INI_BURN_FMASS];
          sf = RDB[DATA_BURN_SFRATE];
          dh = RDB[DATA_BURN_DECAY_HEAT];
          i = 1;
        }
      else if ((long)RDB[DATA_NORM_BURN] == BURN_NORM_NOT_BURN)
        {
          fmass = RDB[DATA_INI_FMASS] - RDB[DATA_INI_BURN_FMASS];
          sf = RDB[DATA_TOT_SFRATE] - RDB[DATA_BURN_SFRATE];
          dh = RDB[DATA_TOT_DECAY_HEAT] - RDB[DATA_BURN_DECAY_HEAT];
          i = 2;
        }
      else
        Die(FUNCTION_NAME, "Invalid normalization");

      /* Override mass if normalization fixed to single material */

      if (mat > VALID_PTR)
        fmass = RDB[mat + MATERIAL_INI_FMASS];

      /* Fission rate */

      if (mat > VALID_PTR)
        {
          ptr = (long)RDB[loc0 + NORM_PTR_FISSRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          fiss = BufVal(ptr, 0);
        }
      else
        {
          ptr = (long)RDB[RES_TOT_FISSRATE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          fiss = BufVal(ptr, i);
        }

      /* Fission neutron production rate rate */

      if (mat > VALID_PTR)
        {
          ptr = (long)RDB[loc0 + NORM_PTR_NSF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          nsf = BufVal(ptr, 0);
        }
      else
        {
          ptr = (long)RDB[RES_TOT_NSF];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          nsf = BufVal(ptr, i);
        }

      /* Fission energy production rate rate */

      if (mat > VALID_PTR)
        {
          ptr = (long)RDB[loc0 + NORM_PTR_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          fissE = BufVal(ptr, 0);
        }
      else
        {
          ptr = (long)RDB[RES_TOT_FISSE];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          fissE = BufVal(ptr, i);
        }

      /* Capture rate */

      ptr = (long)RDB[RES_TOT_CAPTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      capt = BufVal(ptr, i);

      /* Leak rate (no binning) */

      ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      leak = BufVal(ptr, 0);

      /* Flux */

      ptr = (long)RDB[RES_TOT_NEUTRON_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      flx = BufVal(ptr, i);

      /* Source rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      src = BufVal(ptr, i);

      /* Cut-off rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      cut = BufVal(ptr, 0);

      /* Reset decay heat if not used */

      if ((long)RDB[DATA_NORM_INCLUDE_DH] == NO)
        dh = 0.0;

      /***********************************************************************/

      /***** Calculate normalization coefficient *****************************/

      /* Check if normalization is defined */

      if (loc0 > VALID_PTR)
        {
          /* Check mode */

          if ((val = RDB[loc0 + NORM_POWER]) >= 0.0)
            {
              /* Normalize to power */

              val = val + dh;
              div = fissE;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero fission power");
            }
          else if ((val = RDB[loc0 + NORM_POWDENS]) >= 0.0)
            {
              /* Normalize to power density */

              val = val*fmass*1E+6 + dh;
              div = fissE;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero fission power");
            }
          else if ((val = RDB[loc0 + NORM_GENRATE]) >= 0.0)
            {
              /* Normalize to neutron generation rate */

              div = nsf;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero neutron production");
            }
          else if (RDB[loc0 + NORM_SFRATE] >= 0.0)
            {
              /* Normalize to spontaneous fission rate */

              val = sf;
              div = src;

              /* Check divisor */

              if (div < ZERO)
                Die(FUNCTION_NAME, "Zero source rate3");
            }
          else if ((val = RDB[loc0 + NORM_FISSRATE]) >= 0.0)
            {
              /* Normalize to fission rate */

              div = fiss;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero fission rate");
            }
          else if ((val = RDB[loc0 + NORM_ABSRATE]) >= 0.0)
            {
              /* Normalize to absorption rate */

              div = fiss + capt;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero absorption rate");
            }
          else if ((val = RDB[loc0 + NORM_LOSSRATE]) >= 0.0)
            {
              /* Normalize to loss rate */

              div = fiss + capt + leak;

              /* Check divisor */

              if (div < ZERO)
                Die(FUNCTION_NAME, "Zero loss rate");
            }
          else if ((val = RDB[loc0 + NORM_FLUX]) >= 0.0)
            {
              /* Normalize to flux */

              div = flx;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero flux");
            }
          else if ((val = RDB[loc0 + NORM_SRCRATE]) >= 0.0)
            {
              /* Normalize to source rate */

              div = src;

              /* Check divisor */

              if (div < ZERO)
                Die(FUNCTION_NAME, "Zero source rate4");
            }
          else
            Die(FUNCTION_NAME, "Error in normalization");
        }
      else
        {
          /* If not defined, normalize to lossrate = 1 */

          val = 1.0;
          div = fiss + capt + leak + cut;
        }

      /* Multiply by time interval */

      val = val*dt;

      /***********************************************************************/
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      /***********************************************************************/

      /***** Calculate reaction rates for normalization **********************/

      /* Check transport mode */

      if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == NO)
        return 0.0;

      /* Leak rate */

      ptr = (long)RDB[RES_TOT_PHOTON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      leak = BufVal(ptr, 0);

      /* Flux */

      ptr = (long)RDB[RES_TOT_PHOTON_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      flx = BufVal(ptr, 0);

      /* Source rate */

      ptr = (long)RDB[RES_TOT_PHOTON_SRCRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      src = BufVal(ptr, 0);

      /* Heating rate */

      ptr = (long)RDB[RES_TOT_PHOTON_HEATRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      heat = BufVal(ptr, 0);

      /* Energy cut-off rate */

      ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      cut = BufVal(ptr, 0);

      /***********************************************************************/

      /***** Calculate normalization coefficient *****************************/

      /* Check if normalization is defined */

      if ((loc0 = (long)RDB[DATA_PTR_NORM]) > VALID_PTR)
        {
          /* Check mode */

          if ((val = RDB[loc0 + NORM_POWER]) >= 0.0)
            {
              /* Normalize to power */

              div = heat;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero gamma heating power");
            }
          else if ((val = RDB[loc0 + NORM_LOSSRATE]) >= 0.0)
            {
              /* Normalize to loss rate */

              div = leak + cut;

              /* Check divisor */

              if (div < ZERO)
                Die(FUNCTION_NAME, "Zero loss rate");
            }
          else if ((val = RDB[loc0 + NORM_FLUX]) >= 0.0)
            {
              /* Normalize to flux */

              div = flx;

              /* Check divisor */

              if (div < ZERO)
                Error(0, "Error in normalization: zero flux");
            }
          else if ((val = RDB[loc0 + NORM_SRCRATE]) >= 0.0)
            {
              /* Normalize to source rate */

              div = src;

              /* Check divisor */

              if (div < ZERO)
                Die(FUNCTION_NAME, "Zero source rate6");
            }
          else
            Die(FUNCTION_NAME, "Error in normalization");
        }
      else
        {
          /* If not defined, normalize to source rate = 1 */

          val = 1.0;
          div = src;
        }

      /* Multiply by time interval */

      val = val*dt;

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid type");

  /* Calculate coefficient */

  if (val == 0.0)
    norm = 0.0;
  else if (div > 0.0)
    norm = val/div;
  else if (cut > 0.0)
    norm = 1.0;
  else
    Die(FUNCTION_NAME, "Division by zero");

  /* Check value */

  if (!((norm >= 0.0) && (norm < INFTY)))
    Die(FUNCTION_NAME, "Error in normalization (norm = %E)", norm);

  /* Store value */

  if (type == PARTICLE_TYPE_NEUTRON)
    WDB[DATA_NORM_COEF_N] = norm;
  else
    WDB[DATA_NORM_COEF_G] = norm;

  /* Return coefficient */

  return norm;
}

/*****************************************************************************/
