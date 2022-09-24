/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : detresponse.c                                  */
/*                                                                           */
/* Created:       2011/07/08 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates value for detector response function              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DetResponse:"

/*****************************************************************************/

double DetResponse(long det, long loc0, long part, long mat, double E,
                   double g, long id)
{
  long mt, rea, nuc, type, ptr, loc1, loc2, np, n, i;
  double f, T, Er, E0, E1, f0, f1;

  /* Get mt and pointer to reaction */

  mt = (long)RDB[loc0 + DET_RBIN_MT];
  rea = (long)RDB[loc0 + DET_RBIN_PTR_REA];

  /* Get type */
  
  type = (long)RDB[det + DET_PARTICLE];

  /* Check type and get response */

  if (mt == 0)
    f = 1.0;
  else if (mt == MT_NEUTRON_DENSITY)
    {
      /* Check type */

      if (type == PARTICLE_TYPE_GAMMA)
        f = 0.0;
      else
        {
          /* Response function is 1/v */

          f = 1.0/Speed(type, E);
        }
    }
  else if (1 == 2)
    {
      /* Arbitrary response function (not used atm) */

      f = ResponseFunction(E, 0);
    }
  else if (mt == MT_MACRO_MAJORANT)
    {
      /* Get majorant cross section */

      f = DTMajorant(type, E, id);
    }
  else if ((np = (long)RDB[loc0 + DET_RBIN_ATTN_NP]) > 0)
    {
      /* Photon mass attennuation, get pointers */

      loc1 = (long)RDB[loc0 + DET_RBIN_ATTN_PTR_E];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      loc2 = (long)RDB[loc0 + DET_RBIN_ATTN_PTR_F];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Avoid compiler warning */

      f = 0.0;

      /* Grid search */

      if ((n = SearchArray(&RDB[loc1], E, np)) > -1)
        {
          /* Check value */

          CheckValue(FUNCTION_NAME, "n", "", n, 0, np - 2);

          /* Get values */

          E0 = RDB[loc1 + n];
          E1 = RDB[loc1 + n + 1];
          f0 = RDB[loc2 + n];
          f1 = RDB[loc2 + n + 1];

          /* Get interpolated value */

          f = ENDFInterp(5, E, E0, E1, f0, f1);
        }
    }
  else if (mt == MT_USER_DEFINED)
    {
      /* Pointer to function */

      ptr = (long)RDB[loc0 + DET_RBIN_PTR_FUN];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Check type */

      if ((long)RDB[ptr + FUN_TYPE] != 1)
        Die(FUNCTION_NAME, "Function type not supported");

      /* Get number of points */

      np = (long)RDB[ptr + FUN_NE];

      /* Get pointers */

      loc1 = (long)RDB[ptr + FUN_PTR_E];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      loc2 = (long)RDB[ptr + FUN_PTR_F];
      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

      /* Avoid compiler warning */

      f = 0.0;

      /* Grid search */

      if ((n = SearchArray(&RDB[loc1], E, np)) > -1)
        {
          /* Check value */

          CheckValue(FUNCTION_NAME, "n", "", n, 0, np - 2);

          /* Get values */

          E0 = RDB[loc1 + n];
          E1 = RDB[loc1 + n + 1];
          f0 = RDB[loc2 + n];
          f1 = RDB[loc2 + n + 1];

          /* Get interpolated value */

          f = ENDFInterp((long)RDB[ptr + FUN_INT], E, E0, E1, f0, f1);
        }
    }
  else if (mt == MT_PHOTON_DOSE)
    {
      /* Photon dose from elemental data and material composition */

      if ((long)RDB[loc0 + DET_RBIN_VOID_MODE] == NO)
        mat = (long)RDB[loc0 + DET_RBIN_PTR_MAT];

      /* Collision can be in void */

      if (mat > VALID_PTR)
        {
          /* Get number of points */

          np = (long)RDB[mat + MATERIAL_PHOTON_ATT_NE];
          CheckValue(FUNCTION_NAME, "np", "", np, 10, 2000);

          /* Get pointers */

          loc1 = (long)RDB[mat + MATERIAL_PTR_PHOTON_ATT_E];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          loc2 = (long)RDB[mat + MATERIAL_PTR_PHOTON_ATT_F];
          CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

          /* Avoid compiler warning */

          f = 0.0;

          /* Grid search */

          if ((n = SearchArray(&RDB[loc1], E, np)) > -1)
            {
              /* Check value */

              CheckValue(FUNCTION_NAME, "n", "", n, 0, np - 2);

              /* Get values */

              E0 = RDB[loc1 + n];
              E1 = RDB[loc1 + n + 1];
              f0 = RDB[loc2 + n];
              f1 = RDB[loc2 + n + 1];

              /* Get interpolated value */

              f = ENDFInterp(5, E, E0, E1, f0, f1);
            }

          /* Multiply by density factor */

          f = f*g;
        }
      else
        f = 0.0;
    }
  else if (mt < 0)
    {
      /* Check void mode */

      if ((long)RDB[loc0 + DET_RBIN_VOID_MODE] == YES)
        {
          if (mat < VALID_PTR)
            rea = -1;
          else if (mt == MT_MACRO_TOTXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
              else
                rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
            }
          else if (mt == MT_MACRO_ABSXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_ABSXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_ELAXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_ELAXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_INLPRODXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_INLPXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_FISSXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_FISSXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_HEATXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_HEATTXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_PHOTXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_PHOTPXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_PROTPXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_PROTPXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_DEUTPXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_DEUTPXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_TRITPXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_TRITPXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_HE3PXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_HE3PXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_HE4PXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_HE4PXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_FISSE)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_FISSE];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_NSF)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_NSF];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_TMP_MAJORANTXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS];
              else
                rea = -1;
            }
          else if (mt == MT_MACRO_TOTPHOTXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = -1;
              else
                rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
            }
          else if (mt == MT_MACRO_HEATPHOTXS)
            {
              if (type == PARTICLE_TYPE_NEUTRON)
                rea = -1;
              else if ((long)RDB[part + PARTICLE_PHOTON_TYPE] ==
                       PHOTON_TYPE_TTB)
                rea = -1;
              else
                rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS];
            }
          else
            Die(FUNCTION_NAME, "Invalid mt %ld", mt);
        }

      /* Get cross section */

      if (rea > VALID_PTR)
        {
          /* Avoid compiler warning */

          f = -1.0;

          /* Check type */

          if (type == PARTICLE_TYPE_NEUTRON)
            f = MacroXS(rea, E, id);
          else if (type == PARTICLE_TYPE_GAMMA)
            f = PhotonMacroXS(rea, E, id);
          else
            Die(FUNCTION_NAME, "Invalid mode");

          /* Renormalization for energy deposition */
          
          if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
              ((long)RDB[DATA_EDEP_KEFF_CORR]) &&
              ((long)RDB[DATA_EDEP_MODE] > EDEP_MODE_MT458) &&
              (mt == MT_MACRO_HEATXS))
            f = f*RDB[DATA_CYCLE_KEFF];

          /* If fission is not sampled, it is included in capture */
          /* (18.7.2013 / 2.1.15) */

          if ((long)RDB[DATA_NPHYS_SAMPLE_FISS] == NO)
            if (rea == (long)RDB[mat + MATERIAL_PTR_ABSXS])
              if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
                f = f - MacroXS(rea, E, id);

          /* Multiply by density factor */

          f = f*g;
        }
      else
        f = 0.0;
    }
  else
    {
      /* Check reaction pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get nuclide type */

      type = (long)RDB[nuc + NUCLIDE_TYPE];

      /* Get cross section (NOTE: tää tapa millä mikroskooppinen vaikutusala */
      /* lasketaan OTF-moodissa ei ole yhteneväinen muun käytännön kanssa.   */
      /* Tässä oletetaan että vaikutusala on törmäyspisteen lämpötilassa,    */
      /* ilman OTF:ää lämpötila on nuklidikohtainen vakio.)                  */

      if (type == NUCLIDE_TYPE_PHOTON)
        f = PhotonMicroXS(rea, E, id);
      else if ((long)RDB[DATA_TMS_MODE] == TMS_MODE_NONE)
        {
          /* Get microscopic cross section */

          f = MicroXS(rea, E, id);
        }
      else
        {
          if ((mat < VALID_PTR) || ((T = GetTemp(mat, id)) < ZERO))
            {
              /* Get microscopic cross section */

              f = MicroXS(rea, E, id);
            }
          else
            {
              /* Get microscopic cross section */

              f = DopMicroXS(mat, rea, E, &Er, T, id);
            }
        }
    }

  /* Check multiplying by temperature */

  if (RDB[det + DET_TYPE] == (double)DETECTOR_TYPE_T_MULTI)
    {
      /* Get local temperature */

      T = GetTemp(mat, id);

      /* Multiply response by temperature */

      f = f*T;
    }

  /* Check conversion */

  if ((long)RDB[loc0 + DET_RBIN_CONVERT] == DET_CONVERT_GDOSE)
    f = f*E*GRAYH;

  /* Check branching ratio multiplier */

  if ((rea > VALID_PTR) && ((i = (long)RDB[loc0 + DET_RBIN_RFS]) > -1))
    f = f*BranchFrac(rea, i, E, id);

  /* Return value */

  return f;
}

/*****************************************************************************/
