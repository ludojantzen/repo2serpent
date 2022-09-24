/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoredirectsensitivity.c                       */
/*                                                                           */
/* Created:       2018/06/28 (VVa)                                           */
/* Last modified: 2018/09/02 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores reaction rates required to calculate the direct part  */
/*              of detector sensitivity.                                     */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreDirectSensitivity:"

/*****************************************************************************/

void ScoreDirectSensitivity(long part, double wgt, long binidx, long rbin,
                            long mat, double E, double flx, double g, long id)
{
  long sens, label, ptr, lst;
  long imat, izai, irea, iene;
  long mymat;

  long stp, mt, rls, type, rea, nuc, reaType, isfission, idx, iso;
  double partval, XS, adens, T, Er, mult, dnu, pnu;

  /* No need to score sensitivities during inactive cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Get pointer to sensitivity block or return */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /******************************************************************/
  /** Score contribution of different sensitivity bins to detector **/
  /******************************************************************/

  /* Get mt */

  mt = (long)RDB[rbin + DET_RBIN_MT];

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Get response material (void means local material passed as argument) */

  if ((long)RDB[rbin + DET_RBIN_VOID_MODE] == NO)
    mymat = (long)RDB[rbin + DET_RBIN_PTR_MAT];
  else
    mymat = mat;

  if (mymat < VALID_PTR)
    rls = -1;
  else if (mt == MT_MACRO_TOTXS)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
        rls = (long)RDB[mymat + MATERIAL_PTR_TOT_REA_LIST];
      else
        rls = -1;
    }
  else if (mt == MT_MACRO_ABSXS)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
        rls = (long)RDB[mymat + MATERIAL_PTR_ABS_REA_LIST];
      else
        rls = -1;
    }
  else if (mt == MT_MACRO_ELAXS)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
        rls = (long)RDB[mymat + MATERIAL_PTR_ELA_REA_LIST];
      else
        rls = -1;
    }
  else if (mt == MT_MACRO_FISSXS)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
        rls = (long)RDB[mymat + MATERIAL_PTR_FISS_REA_LIST];
      else
        rls = -1;
    }
  else if (mt == MT_MACRO_NSF)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
        rls = (long)RDB[mymat + MATERIAL_PTR_FISS_REA_LIST];
      else
        rls = -1;
    }
  else
    rls = -1;

  /* Get pointer to statistical buffer */

  ptr = (long)RDB[rbin + DET_RBIN_PTR_SENS_DIRECT_STAT_ARRAY];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  stp = (long)RDB[ptr + binidx];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);

  if (rls > VALID_PTR)
    rls = (long)RDB[rls + RLS_PTR_REA0];
  else
    rls = -1;

  /* Loop over reaction list */

  while (rls > VALID_PTR)
    {
      /* Get pointer to reaction */

      rea = (long)RDB[rls + RLS_DATA_PTR_REA];

      /* Check reaction pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get the indices for the perturbation */

      FindSensIndices(mat, rea, E, &imat, &izai, &iene);

      /* Get reaction/event type */

      GetSensEventType(mat, rea, E, &irea, &reaType, &isfission);

      /* Handle all perturbations directly affecting current reaction */

      if ((imat > 0) && (irea > 0) && (izai > 0) && (iene > 0))
        {
          /* Get atomic density */

          if ((long)RDB[rbin + DET_RBIN_VOID_MODE] == YES)
            {
              if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_EXT_ADENS))
                {
                  /* Read atomic density directly from composition */

                  /* Composition index */

                  idx = (long)RDB[rls + RLS_DATA_COMP_IDX];

                  /* Pointer to composition */

                  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
                  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

                  /* Pointer to composition */

                  iso = ListPtr(iso, idx);
                  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

                  /* Atomic density */

                  adens = RDB[iso + COMPOSITION_ADENS];
                }
              else
                {
                  /* Get atomic density from data interface */

                  adens = DataIFCAdens(mat, nuc, id);
                  CheckValue(FUNCTION_NAME, "adens", "from data ifc", adens, 0, 1e5);
                }

              /* Multiply the local atomic density with density factor */

              adens = adens*g;
            }
          else
            adens = 1.0;

          /* Get cross section (NOTE: tää tapa millä mikroskooppinen vaikutusala */
          /* lasketaan OTF-moodissa ei ole yhteneväinen muun käytännön kanssa.   */
          /* Tässä oletetaan että vaikutusala on törmäyspisteen lämpötilassa,    */
          /* ilman OTF:ää lämpötila on nuklidikohtainen vakio.)                  */

          if (type == NUCLIDE_TYPE_PHOTON)
            XS = PhotonMicroXS(rea, E, id);
          else if ((long)RDB[DATA_TMS_MODE] == TMS_MODE_NONE)
            {
              /* Get microscopic cross section */

              XS = MicroXS(rea, E, id);
            }
          else
            {
              if ((mat < VALID_PTR) || ((T = GetTemp(mat, id)) < ZERO))
                {
                  /* Get microscopic cross section */

                  XS = MicroXS(rea, E, id);
                }
              else
                {
                  /* Get microscopic cross section */

                  XS = DopMicroXS(mat, rea, E, &Er, T, id);
                }
            }

          /* Reset multiplier */

          mult = 1.0;

          if (mt == MT_MACRO_NSF)
            {
              /* Get total nubar as multiplier or set it to zero */

              if ((ptr = (long)RDB[rea + REACTION_PTR_TNUBAR]) > VALID_PTR)
                mult = Nubar(ptr, E, id);
              else
                mult = 0.0;
            }

          /* Multiply with atomic density and flux */

          partval = XS*adens*flx*mult;

          /* Compress label */

          label = CompressSensLabel(imat, izai, irea, iene, 1.0);

          /* Score to detector statistic */

          AddBuf1D(partval, wgt, stp, id, label);

          /*************************************************************/
          /* Score direct terms for some additional perturbations here */
          /*************************************************************/

          /* Nubar direct contribution to NSF */

          if ((mt == MT_MACRO_NSF) &&
              ((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_NUBAR))
            {
              /* Get pointer to index list */

              lst = (long)RDB[sens + SENS_PTR_PERT_INDICES];

              /*****************/
              /* Delayed nubar */
              /*****************/

              /* Get index for delayed nubar perturbation */

              irea = (long)RDB[lst + NUBAR_DEL_IDX];

              /* Get delayed nubar as multiplier or set it to zero */

              if ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR)
                dnu = Nubar(ptr, E, id);
              else
                dnu = 0.0;

              /* Multiply with atomic density, flux and nubar */

              partval = XS*adens*flx*dnu;

              /* Compress label */

              label = CompressSensLabel(imat, izai, irea, iene, 1.0);

              /* Score to detector statistic */

              AddBuf1D(partval, wgt, stp, id, label);

              /****************/
              /* Prompt nubar */
              /****************/

              /* Get index for prompt nubar perturbation */

              irea = (long)RDB[lst + NUBAR_PRO_IDX];

              /* Get prompt nubar as multiplier or set it to zero */

              if ((ptr = (long)RDB[rea + REACTION_PTR_TNUBAR]) > VALID_PTR)
                pnu = Nubar(ptr, E, id) - dnu;
              else
                pnu = 0.0;

              /* Multiply with atomic density, flux and nubar*/

              partval = XS*adens*flx*pnu;

              /* Compress label */

              label = CompressSensLabel(imat, izai, irea, iene, 1.0);

              /* Score to detector statistic */

              AddBuf1D(partval, wgt, stp, id, label);
            }
        }

      /* Next item from reaction list */

      rls = NextItem(rls);
    }
}
