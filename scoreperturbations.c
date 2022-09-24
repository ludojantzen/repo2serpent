/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreperturbations.c                           */
/*                                                                           */
/* Created:       2017/05/22 (VVa)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores custom (xGPT) perturbations for sensitivity.          */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScorePerturbations:"

/*****************************************************************************/

void ScorePerturbations(long part, long mat, long rea, long reaType, double E,
                        double hitmiss, double mult, long id)
{
  long loc0, sens, ipert, label, ptr, zai, nuc, mtOK, mt, ene, type, idx;
  double val, E0, E1, val0, val1;

  /* Get pointer to sensitivity block or return */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Check incoming pointers */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get particle coordinates */
  /*
  x = RDB[part + PARTICLE_X];
  y = RDB[part + PARTICLE_Y];
  z = RDB[part + PARTICLE_Z];
  */
  /* Get pointer to first custom perturbation */

  loc0 = (long)RDB[sens + SENS_PTR_PERT0];

  /* Loop over custom perturbations */

  while (loc0 > VALID_PTR)
    {
      /* Check that the material is correct for this perturbation */

      if ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_MAT_LIST]) > VALID_PTR)
        {
          /* Loop over acceptable materials */

          while ((long)RDB[ptr] > VALID_PTR)
            {
              /* Compare material pointer to the one in the list */

              if ((long)RDB[ptr] == mat)
                break;

              /* Get next material pointer from material list */

              ptr++;
            }

          /* Check if material was not found */

          if ((long)RDB[ptr] < VALID_PTR)
            {
              /* Was not found */

              /* Next item */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }
        }

      /* Check that the ZAI is correct for this perturbation */

      if ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_ZAI_LIST]) > VALID_PTR)
        {
          /* Get reaction nuclide */

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get reaction zai */

          zai = (long)RDB[nuc + NUCLIDE_ZAI];

          /* Loop over acceptable zais */

          while ((long)RDB[ptr] > 10010)
            {
              /* Compare zai to the one in the list */

              if ((long)RDB[ptr] == zai)
                break;

              /* Get next zai from zai list */

              ptr++;
            }

          /* Check if zai was not found */

          if ((long)RDB[ptr] < 10010)
            {
              /* Was not found */

              /* Next item */

              loc0 = NextItem(loc0);

              /* Cycle loop */

              continue;
            }
        }

      /* Reset flag for included mt */

      mtOK = 0;

      /* Check if perturbation only affects some reactions */

      if (((long)RDB[loc0 + SENS_PERT_PTR_MT_LIST] < VALID_PTR)
          && ((long)RDB[loc0 + SENS_PERT_PTR_REA_LIST] < VALID_PTR))
        mtOK = 1;
      else
        {

          /* If reaction type list has been given */

          if ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_REA_LIST]) > VALID_PTR)
            {
              /* Loop over reaction list to compare reaction type */

              while ((long)RDB[ptr] >= 0)
                {
                  /* Compare reaction index to the one in the list */

                  if ((long)RDB[ptr] == reaType)
                    break;

                  /* Get next reaction index from list */

                  ptr++;

                }

              /* Check if reaction was found */

              if ((long)RDB[ptr] > 0)
                mtOK = 1;
            }

          /* If reaction mt list has been given and mt is not already OK */

          if ((mtOK == 0) &&
              ((ptr = (long)RDB[loc0 + SENS_PERT_PTR_MT_LIST]) > VALID_PTR))
            {
              /* Get reaction mt */

              mt = (long)RDB[rea + REACTION_MT];

              /* Loop over mt list to compare reaction mt */

              while ((long)RDB[ptr] > 0)
                {
                  /* Compare mt to the one in the list */

                  if ((long)RDB[ptr] == mt)
                    break;

                  /* Get next mt from list */

                  ptr++;
                }

              /* Check if MT was found */

              if ((long)RDB[ptr] > 0)
                mtOK = 1;
            }
        }

      /* Check if the reaction type is OK for perturbation */

      if (!mtOK)
        {
          /* Next item */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Get value at this point */

      if ((ene = (long)RDB[loc0 + SENS_PERT_PTR_EGRID]) > VALID_PTR)
        {
          /* Get pointer to grid data */

          ptr = (long)RDB[ene + ENE_PTR_GRID];

          /* Find index from grid */

          idx = GridSearch(ptr, E);

          /* Check if energy point was not found */

          if (idx < 0)
            {
              loc0 = NextItem(loc0);

              continue;
            }

          /* Get pointer to grid data */

          ptr = (long)RDB[ene + ENE_PTR_GRID];
          CheckPointer(FUNCTION_NAME, "(ptr, 1)", DATA_ARRAY, ptr);

          ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
          CheckPointer(FUNCTION_NAME, "(ptr, 1)", DATA_ARRAY, ptr);

          /* Get energies */

          E0 = RDB[ptr + idx];
          E1 = RDB[ptr + idx + 1];

          /* Get pointer to values */

          ptr = (long)RDB[loc0 + SENS_PERT_PTR_EGRID_VAL];

          /* Get values */

          val0 = RDB[ptr + idx];
          val1 = RDB[ptr + idx + 1];

          /* Get interpolation type */

          type = (long)RDB[loc0 + SENS_PERT_EGRID_INTERP];

          /* Interpolate to get current value */

          val = ENDFInterp(type, E, E0, E1, val0, val1);
        }
      else
        val = 1.0;

      /* No use scoring zero valued perturbations */

      if (val == 0.0)
        {
          loc0 = NextItem(loc0);

          continue;
        }

      /* Get index of this perturbation */

      ipert = (long)RDB[loc0 + SENS_PERT_INDEX];

      /* Calculate Sens label */

      label = CompressSensLabel(0, 0, ipert, 0, hitmiss);

      /* Create a new event */

      StoreSensEvent(part, label, E, val*mult, id);

      /* Get next perturbation */

      loc0 = NextItem(loc0);
    }

}
