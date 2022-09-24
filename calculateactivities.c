/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculateactivities.c                          */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2019/04/12 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates activities, decay heat and spontaneous fission    */
/*              rates for materials                                          */
/*                                                                           */
/* Comments: - Fotonitransportmoodissa missä materiaalikoostumukset on       */
/*             annettu neutroni / decay -datana alkuperäiset listat on       */
/*             korvattu replacephotondata.c:ssä. Alkuperäinen lista luetaan  */
/*             erillisestä paikasta.                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateActivities:"

/*****************************************************************************/

void CalculateActivities()
{
  long mat, iso, nuc, ptr, rad, type, loc0;
  double vol, adens, sum1, sum2, A, H, SF, HT, GT, I, f;

  /* Check rmx self-adapted mode */

  if ((long)RDB[DATA_RMTX_MFP_CALC] == YES)
    return;

  /* Check if decay data file is given */

  if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] < 1)
    return;

  fprintf(outp, "Calculating activities...\n");

  /***************************************************************************/

  /***** Normalization to activity *******************************************/

  /* Reset coefficient */

  f = 1.0;

  /* Check normalization */

  A = -1.0;
  if ((loc0 = (long)RDB[DATA_PTR_NORM]) > VALID_PTR)
    A = RDB[loc0 + NORM_ACTI_A];

  /* Check value */

  if (A > 0)
    {
      /* Reset sum */

      sum1 = 0.0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Check physical flag */

          if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
            {
              /* Pointer to next */

              mat = NextItem(mat);

              /* Cycle loop */

              continue;
            }

          /* Check if normalization is fixed to material */

          if ((long)RDB[loc0 + NORM_PTR_MAT] > VALID_PTR)
            if ((long)RDB[loc0 + NORM_PTR_MAT] != mat)
              {
                /* Pointer to next */

                mat = NextItem(mat);

                /* Cycle loop */

                continue;
              }

          /* Check divisor type */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
            Die(FUNCTION_NAME, "Divided parent material");

          /* Get volume */

          vol = RDB[mat + MATERIAL_VOLUME];

          /* Check (tällä on korvattu se fail-hässäkkä 20.11.2012 / 2.1.10) */

          if ((vol > 0.0) && (vol < 1E+18))
            {
              /* Loop over composition */

              if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
                  < VALID_PTR)
                iso = (long)RDB[mat + MATERIAL_PTR_COMP];

              while (iso > VALID_PTR)
                {
                  /* Pointer to nuclide */

                  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Get atomic density */

                  adens = RDB[iso + COMPOSITION_ADENS];

                  /* Add to activity */

                  if (((long)RDB[loc0 + NORM_ACTI_ZAI] == -1) ||
                      ((long)RDB[nuc + NUCLIDE_ZAI] ==
                       (long)RDB[loc0 + NORM_ACTI_ZAI]))
                    sum1 = sum1 + vol*adens*RDB[nuc + NUCLIDE_LAMBDA]*1E+24;

                  /* Next */

                  iso = NextItem(iso);
                }
            }

          /* Next material */

          mat = NextItem(mat);
        }

      /* Check sum */

      if (sum1 > 0.0)
        f = A/sum1;
      else if ((mat = (long)RDB[loc0 + NORM_PTR_MAT]) > VALID_PTR)
        Error(loc0, "No activity for %ld in material %s",
              (long)RDB[loc0 + NORM_ACTI_ZAI],
              GetText(mat + MATERIAL_PTR_NAME));
      else
        Error(loc0, "No activity for %ld", (long)RDB[loc0 + NORM_ACTI_ZAI]);
    }

  /***************************************************************************/

  /***** Calculate total activities ******************************************/

  /* Check if source material is given */

  if ((mat = (long)RDB[DATA_NORM_PTR_RAD_SRC_MAT]) > VALID_PTR)
    if (RDB[mat + MATERIAL_VOLUME] == 0.0)
      Error(mat, "Volume of source material \"%s\" is not defined",
            GetText(mat + MATERIAL_PTR_NAME));

  /* Reset global data */

  WDB[DATA_TOT_ACTIVITY] = 0.0;
  WDB[DATA_TOT_SFRATE] = 0.0;
  WDB[DATA_TOT_DECAY_HEAT] = 0.0;
  WDB[DATA_ACT_ACTIVITY] = 0.0;
  WDB[DATA_ACT_DECAY_HEAT] = 0.0;
  WDB[DATA_FP_ACTIVITY] = 0.0;
  WDB[DATA_FP_DECAY_HEAT] = 0.0;
  WDB[DATA_TOT_ING_TOX] = 0.0;
  WDB[DATA_TOT_INH_TOX] = 0.0;
  WDB[DATA_ACT_ING_TOX] = 0.0;
  WDB[DATA_ACT_INH_TOX] = 0.0;
  WDB[DATA_FP_ING_TOX] = 0.0;
  WDB[DATA_FP_INH_TOX] = 0.0;

  WDB[DATA_SR90_ACTIVITY] = 0.0;
  WDB[DATA_TE132_ACTIVITY] = 0.0;
  WDB[DATA_I131_ACTIVITY] = 0.0;
  WDB[DATA_I132_ACTIVITY] = 0.0;
  WDB[DATA_CS134_ACTIVITY] = 0.0;
  WDB[DATA_CS137_ACTIVITY] = 0.0;

  WDB[DATA_TOT_PHOTON_DEC_SRC_RATE] = 0.0;
  WDB[DATA_PHOTON_DEC_SRC_MAX_I] = 0.0;
  WDB[DATA_PHOTON_DEC_SRC_VOL] = 0.0;
  WDB[DATA_TOT_NEUTRON_DEC_SRC_RATE] = 0.0;
  WDB[DATA_NEUTRON_DEC_SRC_MAX_I] = 0.0;
  WDB[DATA_NEUTRON_DEC_SRC_VOL] = 0.0;

  WDB[DATA_TOT_ALPHA_DEC_SRC_RATE] = 0.0;
  WDB[DATA_TOT_BETA_DEC_SRC_RATE] = 0.0;

  /* Reset material-wise values */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset values */

      WDB[mat + MATERIAL_ACTIVITY] = 0.0;
      WDB[mat + MATERIAL_SFRATE] = 0.0;
      WDB[mat + MATERIAL_DECAY_HEAT] = 0.0;
      WDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE] = 0.0;
      WDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE] = 0.0;

      /* Pointer to next */

      mat = NextItem(mat);
    }

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check physical flag */

      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
        {
          /* Pointer to next */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Check divisor type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
        Die(FUNCTION_NAME, "Divided parent material");

      /* Get volume */

      vol = RDB[mat + MATERIAL_VOLUME];

      /* Check (tällä on korvattu se fail-hässäkkä 20.11.2012 / 2.1.10) */

      if ((vol > 0.0) && (vol < 1E+18))
        {
          /* Loop over composition */

          if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP])
              < VALID_PTR)
            iso = (long)RDB[mat + MATERIAL_PTR_COMP];

          while (iso > VALID_PTR)
            {
              /* Pointer to nuclide */

              nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
              CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

              /* Get atomic density */

              adens = RDB[iso + COMPOSITION_ADENS];

              /* Calculate activity, decay heat toxicities */

              A = vol*adens*RDB[nuc + NUCLIDE_LAMBDA]*1E+24*f;
              H = A*RDB[nuc + NUCLIDE_DECAY_E]*MEV;
              SF = A*RDB[nuc + NUCLIDE_SF_BR];
              HT = A*RDB[nuc + NUCLIDE_SPEC_INH_TOX];
              GT = A*RDB[nuc + NUCLIDE_SPEC_ING_TOX];

              /* Check values */

              CheckValue(FUNCTION_NAME, "A", "", A, 0.0, 1E+30);
              CheckValue(FUNCTION_NAME, "H", "", H, 0.0, 1E+30);
              CheckValue(FUNCTION_NAME, "SF", "", SF, 0.0, A);

              /* Add to totals */

              WDB[DATA_TOT_ACTIVITY] = RDB[DATA_TOT_ACTIVITY] + A;
              WDB[DATA_TOT_DECAY_HEAT] = RDB[DATA_TOT_DECAY_HEAT] + H;
              WDB[DATA_TOT_SFRATE] = RDB[DATA_TOT_SFRATE] + SF;
              WDB[DATA_TOT_INH_TOX] = RDB[DATA_TOT_INH_TOX] + HT;
              WDB[DATA_TOT_ING_TOX] = RDB[DATA_TOT_ING_TOX] + GT;

              /* Loop over radiations */

              rad = (long)RDB[nuc + NUCLIDE_PTR_RADIATIONS];
              while (rad > VALID_PTR)
                {
                  /* Get intensity */

                  I = A*RDB[rad + NUCLIDE_RAD_SPEC_I];
                  CheckValue(FUNCTION_NAME, "I", "", I, 0.0, 100.0*A);

                  /* Get type */

                  type = (long)RDB[rad + NUCLIDE_RAD_TYPE];

                  /* Check type */

                  if (type == PARTICLE_TYPE_NEUTRON)
                    {
                      /* Store to neutron source rates */

                      WDB[DATA_TOT_NEUTRON_DEC_SRC_RATE] =
                        RDB[DATA_TOT_NEUTRON_DEC_SRC_RATE] + I;

                      WDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE] =
                        RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE] + I;
                    }
                  else if (type == PARTICLE_TYPE_GAMMA)
                    {
                      /* Store to photon source rates */

                      WDB[DATA_TOT_PHOTON_DEC_SRC_RATE] =
                        RDB[DATA_TOT_PHOTON_DEC_SRC_RATE] + I;

                      WDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE] =
                        RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE] + I;
                    }
                  else if (type == PARTICLE_TYPE_ELECTRON)
                    {
                      /* Store to total source rate */

                      WDB[DATA_TOT_BETA_DEC_SRC_RATE] =
                        RDB[DATA_TOT_BETA_DEC_SRC_RATE] + I;
                    }
                  else if (type == PARTICLE_TYPE_ALPHA)
                    {
                      /* Store to total source rate */

                      WDB[DATA_TOT_ALPHA_DEC_SRC_RATE] =
                        RDB[DATA_TOT_ALPHA_DEC_SRC_RATE] + I;
                    }
                  else
                    Die(FUNCTION_NAME, "Invalid radiation type");

                  /* Next radiation */

                  rad = NextItem(rad);
                }

              /* Add to burnable material values */

              if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
                {
                  WDB[DATA_BURN_DECAY_HEAT] =
                    RDB[DATA_BURN_DECAY_HEAT] + H;
                  WDB[DATA_BURN_SFRATE] = RDB[DATA_BURN_SFRATE] + SF;
                }

              /* Add partials */

              if ((long)RDB[nuc + NUCLIDE_Z] > 89)
                {
                  WDB[DATA_ACT_ACTIVITY] = RDB[DATA_ACT_ACTIVITY] + A;
                  WDB[DATA_ACT_DECAY_HEAT] = RDB[DATA_ACT_DECAY_HEAT] + H;
                  WDB[DATA_ACT_INH_TOX] = RDB[DATA_ACT_INH_TOX] + HT;
                  WDB[DATA_ACT_ING_TOX] = RDB[DATA_ACT_ING_TOX] + GT;
                }
              else if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
                       NUCLIDE_FLAG_FP)
                {
                  WDB[DATA_FP_ACTIVITY] = RDB[DATA_FP_ACTIVITY] + A;
                  WDB[DATA_FP_DECAY_HEAT] = RDB[DATA_FP_DECAY_HEAT] + H;
                  WDB[DATA_FP_INH_TOX] = RDB[DATA_FP_INH_TOX] + HT;
                  WDB[DATA_FP_ING_TOX] = RDB[DATA_FP_ING_TOX] + GT;
                }

              /* Add isotopic */

              if ((long)RDB[nuc + NUCLIDE_ZAI] == 380900)
                WDB[DATA_SR90_ACTIVITY] = RDB[DATA_SR90_ACTIVITY] + A;
              else if ((long)RDB[nuc + NUCLIDE_ZAI] == 521320)
                WDB[DATA_TE132_ACTIVITY] = RDB[DATA_TE132_ACTIVITY] + A;
              else if ((long)RDB[nuc + NUCLIDE_ZAI] == 531310)
                WDB[DATA_I131_ACTIVITY] = RDB[DATA_I131_ACTIVITY] + A;
              else if ((long)RDB[nuc + NUCLIDE_ZAI] == 531320)
                WDB[DATA_I132_ACTIVITY] = RDB[DATA_I132_ACTIVITY] + A;
              else if ((long)RDB[nuc + NUCLIDE_ZAI] == 551340)
                WDB[DATA_CS134_ACTIVITY] = RDB[DATA_CS134_ACTIVITY] + A;
              else if ((long)RDB[nuc + NUCLIDE_ZAI] == 551370)
                WDB[DATA_CS137_ACTIVITY] = RDB[DATA_CS137_ACTIVITY] + A;

              /* Add to material data */

              WDB[mat + MATERIAL_ACTIVITY] =
                RDB[mat + MATERIAL_ACTIVITY] + A;

              WDB[mat + MATERIAL_DECAY_HEAT] =
                RDB[mat + MATERIAL_DECAY_HEAT] + H;

              WDB[mat + MATERIAL_SFRATE] =
                RDB[mat + MATERIAL_SFRATE] + SF;

              /* Next isotope */

              iso = NextItem(iso);
            }

          /* Compare to maximum intensities */

          if (RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE]/vol >
              RDB[DATA_PHOTON_DEC_SRC_MAX_I])
            WDB[DATA_PHOTON_DEC_SRC_MAX_I] =
              RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE]/vol;

          if (RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE]/vol >
              RDB[DATA_NEUTRON_DEC_SRC_MAX_I])
            WDB[DATA_NEUTRON_DEC_SRC_MAX_I] =
              RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE]/vol;

          /* Add to source volumes */

          if (RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE] > 0.0)
            WDB[DATA_PHOTON_DEC_SRC_VOL] =
              RDB[DATA_PHOTON_DEC_SRC_VOL] + vol;

          if (RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE] > 0.0)
            WDB[DATA_NEUTRON_DEC_SRC_VOL] =
              RDB[DATA_NEUTRON_DEC_SRC_VOL] + vol;
        }

      /* Add to parent (JLe: Tää lisättiin 1.10.2015 / 2.1.25 että  */
      /* decay gamma sourcen normeerauksen saa kiinnitettyä parent- */
      /* materiaaliin. Voi olla että tää sotkee jotain kokonais-    */
      /* arvoihin liittyviä laskuja) */

      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > 0.0)
        {
          /* Sum only photon source rate for now  */
          /* (needed in radioactive decay source) */

          WDB[ptr + MATERIAL_PHOTON_DEC_SRC_RATE] =
            RDB[ptr + MATERIAL_PHOTON_DEC_SRC_RATE] +
            RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE];

          WDB[ptr + MATERIAL_NEUTRON_DEC_SRC_RATE] =
            RDB[ptr + MATERIAL_NEUTRON_DEC_SRC_RATE] +
            RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE];
        }

      /* Next material */

      mat = NextItem(mat);
    }

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/

  /***** Set decay source rate used for normalization ************************/

  /* Check radioactive decay source */

  if ((long)RDB[DATA_USE_DECAY_SRC] == YES)
    {
      /* Avoid compiler warning */

      I = 0.0;

      /* Get material pointer */

      if ((mat = (long)RDB[DATA_NORM_PTR_RAD_SRC_MAT]) > VALID_PTR)
        {
          /* Check transport mode */

          if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
            I = RDB[mat + MATERIAL_NEUTRON_DEC_SRC_RATE];
          else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
            I = RDB[mat + MATERIAL_PHOTON_DEC_SRC_RATE];
          else
            Die(FUNCTION_NAME, "WTF?");
        }
      else
        {
          /* Check transport mode */

          if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
            I = RDB[DATA_TOT_NEUTRON_DEC_SRC_RATE];
          else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
            I = RDB[DATA_TOT_PHOTON_DEC_SRC_RATE];
          else
            Die(FUNCTION_NAME, "WTF?");
        }

      /* Check value */

      if (I <= 0.0)
        {
          /* Check if neutron source */

          if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
            Error(0, "No neutron decay source (may be due to lack of data)");
          else
            Error(0, "Zero emission rate in decay source mode");
        }

      /* Put value for normalization */

      WDB[DATA_NORM_DECAY_SRC_RATE] = I;
    }
  else
    WDB[DATA_NORM_DECAY_SRC_RATE] = 0.0;

  /***************************************************************************/

  /***** Xenon entropy *******************************************************/

  /* Calculate mean concentration */

  sum2 = 0.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check physical flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)
        {
          /* Check divisor type */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
            Die(FUNCTION_NAME, "Divided parent material");

          /* Add Xe-135 concentration */

          if ((iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
            sum2 = sum2 + RDB[iso + COMPOSITION_ADENS];
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Calculate entropies */

  sum1 = 0.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check physical flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)
        {
          /* Check divisor type */

          if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
            Die(FUNCTION_NAME, "Divided parent material");

          /* Add Xe-135 concentration */

          if (sum2 > 0.0)
            if ((iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
              {
                /* Get density */

                adens = RDB[iso + COMPOSITION_ADENS];

                /* Add to entropy */

                sum1 = sum1 - log2(adens/sum2)/RDB[DATA_N_BURN_MATERIALS];
              }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Put value */

  if (sum1 > 0.0)
    WDB[DATA_XENON_ENTROPY] = log2(RDB[DATA_N_BURN_MATERIALS])/sum1;

  /***************************************************************************/
}

/*****************************************************************************/
