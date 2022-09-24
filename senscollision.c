/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : senscollision.c                                */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Samples between accepted and rejected collisions for         */
/*              sensitivity calculations. Scores and stores events.          */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SensCollision"

/*****************************************************************************/

long SensCollision(long mat, long rea, long part, double E, long id)
{
  long imat, izai, iene, isfission, is_na_coll, irea, label;
  long sens, reaType, hitmiss, ptr, ncol, nuc, dum, iene2, ptr1, msh, uni;
  double Et, k, T, kT, tmpdbl, dFEdT, T0, x, y, z, val;

  /* Get pointer to sensitivity block or return a hit */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return COLL_HIT;

  /* Check particle type */

  if ((long)RDB[part + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
    return COLL_HIT;

  /* Get the indices for the perturbation */

  FindSensIndices(mat, rea, E, &imat, &izai, &iene);

  is_na_coll = 0;

  /* Check if the material index correspons to void fraction material */

  if (mat == (long)RDB[sens + SENS_PTR_VOID_MAT])
    is_na_coll = 1;

  /* Reset reaction and event type */

  irea = -1;
  isfission = -1;

  /* Get reaction/event type */

  GetSensEventType(mat, rea, E, &irea, &reaType, &isfission);

  /* Sample rejection due to Sens */

  if (RandF(id) < 0.5)
    hitmiss = COLL_MISS;
  else
    hitmiss = COLL_HIT;

  /* Score custom sensitivity perturbations */

  ScorePerturbations(part, mat, rea, reaType, E, (double)hitmiss, 1.0, id);

  /* Modify number of collisions in void coef material */

  if (is_na_coll == 1)
    WDB[part + PARTICLE_NA_COLL] = RDB[part + PARTICLE_NA_COLL] + (double)hitmiss;

  /* If target nuclide is not included in perturbation list and not */
  /* scoring catch-all nuclide, do not create event */

  if (izai < 0)
    return hitmiss;

  /* If target material is not included in perturbation list and not */
  /* scoring catch-all material, do not create event */

  if (imat < 0)
    return hitmiss;

  /* Get target energy */

  if (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_TEMPERATURE) > 0)
    {

      /******************************************************************/
      /********** Handle cross section temperature effect ***************/
      /******************************************************************/

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];

      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);

      if ((Et = TestValuePair(nuc + NUCLIDE_PREV_COL_ET, (double)ncol, id))
          < 0.0)
        Die(FUNCTION_NAME, "Target energy hasn't been sampled for %ld, col %ld", (long)RDB[nuc + NUCLIDE_ZAI], ncol);
      else if ((Et < 1e9) && (Et > ZERO))
        {
          /* Get indices */

          FindSensIndices(mat, rea, Et, &dum, &dum, &iene2);

          /* Get pointer to index list */

          ptr1 = (long)RDB[sens + SENS_PTR_PERT_INDICES];
          CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);

          /* Override energy idx for temperature sensitivity */

          if ((msh = (long)RDB[sens + SENS_PTR_MESH]) > VALID_PTR)
            {
              /* Get pointer to root universe */

              uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

              /* Get coordinates */

              ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              x = GetPrivateData(ptr, id);

              ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              y = GetPrivateData(ptr, id);

              ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              z = GetPrivateData(ptr, id);

              /* Get mesh index */

              if ((iene2 = MeshIndex(msh, x, y, z, -1.0)) < 0)
                iene2 = 0;
              else
                iene2 += 1;
            }

          /* Calculate label for accepted */

          label = CompressSensLabel(imat, izai, (long)RDB[ptr1 + TEMPERATURE_IDX], iene2, hitmiss);

          /* Calculate derivative of distribution here */

          k = 1.38064852e-23;

          /* Get temperature used for sampling target energy */

          if ((T = TestValuePair(nuc + NUCLIDE_PREV_COL_DT, (double)ncol, id))
              < 0.0)
            Die(FUNCTION_NAME, "Could not obtain collision Maxwellian temperature "
                "for %ld in %s, col %ld", (long)RDB[nuc + NUCLIDE_ZAI],
                GetText(mat + MATERIAL_PTR_NAME), ncol);

          kT = k*T/MEV*1e6;

          Et = Et*1e6;

          /* Derivative of MB distribution wrt temperature */

          /*
            dFEdT = (2*pow(Et,(3.0/2.0))*exp(-Et/(kT))/(sqrt(PI)*pow(k,(5.0/2.0))*pow(T,(7.0/2.0))) -
            3*sqrt(Et)*exp(-Et/(kT))/(sqrt(PI)*pow(k,(3.0/2.0))*pow(T,(5.0/2.0))));
          */

          /* Get temperature at collision point */

          if ((T0 = TestValuePair(nuc + NUCLIDE_PREV_COL_T, (double)ncol, id))
              < 0.0)
            Die(FUNCTION_NAME, "Could not obtain collision material temperature "
                "for %ld in %s, col %ld", (long)RDB[nuc + NUCLIDE_ZAI],
                GetText(mat + MATERIAL_PTR_NAME), ncol);

          /* Calculate the multiplier to the score */
          /* val = (dFE/FE)/(dT/T) or actually (dFE/dT/FE)*(1/T) */

          val = (Et/(kT*T) - 1.5/T)/(1.0/T0);

          /* Create event */

          StoreSensEvent(part, label, Et/1e6, val, id);
        }

      /******************************************************************/
      /********** Handle elastic scattering temperature effect **********/
      /******************************************************************/

      /* This doesn't really work yet, doesn't take in account DBRC correctly? */

      if ((reaType == (long)ELA_SCATT_IDX) && (RDB[nuc + NUCLIDE_XS_TEMP] > 0))
        {

          /* Reaction is an elastic scattering */

          /* Get temperature for target velocity */

          if(RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
            T = GetTemp(mat, id);
          else
            T = RDB[nuc + NUCLIDE_TEMP];

          kT = T*KELVIN;

          /* Sample target velocity already here */
          /* (we don't need to store the stuff as it is stored to nuclide) */

          TargetVelocity(rea, E, &tmpdbl, &tmpdbl, &tmpdbl,
                         1, 0, 0, kT, id);

          /* Check if target energy has been sampled */

          if (RDB[nuc + NUCLIDE_XS_TEMP] == 0)
            {
              /* Sampled at dopmicroxs.c */

              ptr = nuc + NUCLIDE_PREV_COL_ET;
            }
          else
            {
              /* Sampled at targetvelocity.c */

              ptr = nuc + NUCLIDE_PREV_COL_TV_ET;
            }

          /* Check if targetvelocity did sample an energy */

          if ((Et = TestValuePair(ptr, (double)ncol, id))
              < 0.0)
            Die(FUNCTION_NAME, "Scattering target energy hasn't "
                "been sampled for %ld, col %ld", (long)RDB[nuc + NUCLIDE_ZAI], ncol);
          else if (Et > ZERO)
            {
              /* Get indices */

              FindSensIndices(mat, rea, Et, &dum, &dum, &iene2);

              /* Get pointer to index list */

              ptr1 = (long)RDB[sens + SENS_PTR_PERT_INDICES];
              CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);

              /* Override energy idx for temperature sensitivity */

              if ((msh = (long)RDB[sens + SENS_PTR_MESH]) > VALID_PTR)
                {
                  /* Get pointer to root universe */

                  uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
                  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

                  /* Get coordinates */

                  ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_X];
                  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
                  x = GetPrivateData(ptr, id);

                  ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Y];
                  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
                  y = GetPrivateData(ptr, id);

                  ptr = (long)RDB[uni + UNIVERSE_PTR_PRIVA_Z];
                  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
                  z = GetPrivateData(ptr, id);

                  /* Get mesh index */

                  if ((iene2 = MeshIndex(msh, x, y, z, -1.0)) < 0)
                    iene2 = 0;
                  else
                    iene2 += 1;
                }

              /* Calculate label for accepted */

              label = CompressSensLabel(imat, izai, (long)RDB[ptr1 + SCATT_TEMP_IDX], iene2, hitmiss);

              /* Calculate derivative of distribution here */

              k = 1.38064852e-23;

              /* Get temperature used for sampling target energy */

              kT = k*T/MEV*1e6;

              Et = Et*1e6;

              dFEdT = (Et/(kT*T) - 1.5/T)/(1.0/T);

              /* Create  event */

              StoreSensEvent(part, label, Et/1e6, dFEdT, id);
            }
        }
    }



  /* If not perturbing basic reaction cross sections, do not create event */

  if ((((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XS) == 0) &&
      (((long)RDB[sens + SENS_PERT_FLAGS] & SENS_PERT_FLAG_XSMT) == 0))
    return hitmiss;

  /* If the reaction is a (successfull) fission, event is created in */
  /* fission.c */

  if (isfission*hitmiss == 1)
    return hitmiss;

  /* Calculate Sens label */

  label = CompressSensLabel(imat, izai, irea, iene, hitmiss);

  /* Create a new event */

  StoreSensEvent(part, label, E, 1.0, id);

  return hitmiss;

}
