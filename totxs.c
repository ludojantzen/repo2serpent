/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : totxs.c                                        */
/*                                                                           */
/* Created:       2011/07/21 (JLe)                                           */
/* Last modified: 2018/06/05 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns continuous-energy or coarse multi-group total xs for */
/*              neutrons and photons                                         */
/*                                                                           */
/* Comments: - Used in Tracking() and Collision()                            */
/*                                                                           */
/*           - Tota arvoa ei voi tallentaa sillä samasta rea-structuresta    */
/*             haetaan sekä CE että MG -vaikutusalat.                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TotXS:"

/*****************************************************************************/

double TotXS(long mat, long type, double E, long id)
{
  long rea, ptr, mt;
  double xs;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return 0.0;

  /* Avoid compiler warning */

  rea = -1;

  /* Check type and get reaction pointer */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Use temperature majorant or total */

      if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE)
        rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS];
      else
        rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      /* Use total photon cross section */

      rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
    }
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Avoid compiler warning */

  xs = 0.0;

  /* Check for coarse multi-group xs */

  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
    {
      /* Get cross section (NOTE: Nää lisättiin 9.2.2018 / 2.1.28 / JLE */
      /* vastaavat tehdään tuolla MacroXS():n puolella. */

      xs = MGXS(rea, E, -1);

      /* Add poison cross section */

      if (((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == YES) ||
          ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == YES))
        xs = xs + PoisonXS(mat, E, MT_MACRO_TOTXS, id);

      /* Add iterated cross section */

      if ((long)RDB[DATA_ITER_MODE] == ITER_MODE_NUCLIDE)
        xs = xs + IterNucXS(mat, E, MT_MACRO_TOTXS, -1, id);

      /* Add cross sections from data interface nuclides */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_EXT_ADENS_MAT)
        xs = xs + DataIFCXS(mat, E, MT_MACRO_TOTXS, -1, id);

      /* Add cross sections for on-the-fly burnup mode */

      if ((long)RDB[DATA_OTF_BURN_MODE] == YES)
        xs = xs + OTFBurnXS(mat, E, MT_MACRO_TOTXS, id);
    }
  else
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Check mt for type */

      if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_TMP_MAJORANTXS))
        xs = MacroXS(rea, E, id);
      else if (mt == MT_MACRO_TOTPHOTXS)
        xs = PhotonMacroXS(rea, E, id);
      else
        Die(FUNCTION_NAME, "Invalid reaction mode");
    }

  /* Return cross section */

  return xs;
}

/*****************************************************************************/
