/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoremicrodep.c                                */
/*                                                                           */
/* Created:       2017/01/01 (JLe)                                           */
/* Last modified: 2020/03/04 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Scores micro-depletion cross sections                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreMicroDep:"

/*****************************************************************************/

void ScoreMicroDep(double flx, long mat, double E, double wgt, long id)
{
  long rea0, rea, stp, nuc, loc0, loc1, ptr, ng, ntot, RFS, gcu, idx, iso, mt;
  double val, E0, E1, E2, f, g, T, Er, adens, mul;
  const double *adens0;

  /* Number of energy groups */

  ntot = (long)RDB[DATA_ERG_FG_NG];

  /* Get pointer to few-group structure */

  ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get few-group index */

  if ((ng = GridSearch(ptr, E)) > -1)
    {
      /* Convert index */

      ng = ntot - ng - 1;
      CheckValue(FUNCTION_NAME, "ng", "", ng, 0, ntot - 1);
    }
  else
    return;

  /* Loop over micro depletion data */

  loc0 = (long)RDB[DATA_PTR_MDEP0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over materials */

      loc1 = (long)RDB[loc0 + MDEP_PTR_MAT];
      while (loc1 > VALID_PTR)
        {
          /* Compare pointes */

          if ((long)RDB[loc1 + MDEP_MAT_PTR_MAT] == mat)
            break;

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Check pointer */

      if (loc1 < VALID_PTR)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Pointer to average densities */

      ptr = (double)RDB[loc0 + MDEP_PTR_AVG_ADENS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      adens0 = &RDB[ptr];

      /* Pointer to gcu */

      gcu = (long)RDB[loc0 + MDEP_PTR_GCU];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);

      /* Pointer to stats */

      stp = (long)RDB[gcu + GCU_MICRO_MICRO_DEP_XS];
      CheckPointer(FUNCTION_NAME, "(stp)", RES2_ARRAY, stp);

      /* Loop over reactions */

      loc1 = (long)RDB[loc1 + MDEP_MAT_PTR_REA];
      while (loc1 > VALID_PTR)
        {
          /* Compare to minimum energy */

          if (E < RDB[loc1 + MDEP_REA_EMIN])
            break;

          /* Pointer to reaction data */

          rea0 = (long)RDB[loc1 + MDEP_REA_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

          /* Get pointer to parent data if branching to isomeric state */

          if ((rea = (long)RDB[rea0 + REACTION_PTR_BRANCH_PARENT]) < VALID_PTR)
            rea = rea0;

          nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Check if reaction is associated with isomeric branching */

          if ((long)RDB[rea + REACTION_PTR_ISO_BRA] > VALID_PTR)
            {
              /* Get fraction */

              RFS = (long)RDB[rea0 + REACTION_RFS];
              f = BranchFrac(rea0, RFS, E, id);
            }
          else
            f = 1.0;

          /* Get cross section */

          if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_NONE)
            val = MicroXS(rea, E, id);
          else if ((T = GetTemp(mat, id)) < ZERO)
            val = MicroXS(rea, E, id);
          else
            val = DopMicroXS(mat, rea, E, &Er, T, id);

          /* Get MT */

          mt = (long)RDB[loc1 + MDEP_REA_MT];

          /* Check special multipliers for fission */

          if (mt == 452)
            {
               /* Get total nubar */

               if ((ptr = (long)RDB[rea + REACTION_PTR_TNUBAR]) > VALID_PTR)
                 mul = Nubar(ptr, E, id);
               else
                 mul = 0.0;

               /* Multiply */

               val = val*mul;
            }
          else if (mt == 455)
            {
               /* Get delayed nubar */

              if ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR)
                mul = Nubar(ptr, E, id);
              else
                mul = 0.0;

               /* Multiply */

               val = val*mul;
            }
          else if (mt == 456)
            {
               /* Get total nubar */

               if ((ptr = (long)RDB[rea + REACTION_PTR_TNUBAR]) > VALID_PTR)
                 mul = Nubar(ptr, E, id);
               else
                 mul = 0.0;

               /* Subtract delayed nubar */

              if ((ptr = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR)
                mul = mul - Nubar(ptr, E, id);

              /* Multiply */

              val = val*mul;
            }
          else if (mt == 458)
            {
              /* Get fission energy */

              mul = FissE(rea, E, id)/MEV;
              CheckValue(FUNCTION_NAME, "fE", "", mul, 0.0, 300.0);

              /* Multiply */

              val = val*mul;
            }

          /* Get fission yield interpolation energies */

          E0 = RDB[rea0 + REACTION_FISSY_IE0];
          E1 = RDB[rea0 + REACTION_FISSY_IE1];
          E2 = RDB[rea0 + REACTION_FISSY_IE2];

          /* Check regions (assume E0 == 0 if not fission) */

          if (E0 == 0.0)
            g = 1.0;
          else if ((E < E0) || (E > E2))
            g = 0.0;
          else
            {
              /* Interpolation factor */

              if (E < E1)
                {
                  /* Below interpolation energy */

                  if (E0 < 0.0)
                    g = 1.0;
                  else
                    g = (E - E0)/(E1 - E0);
                }
              else
                {
                  /* Above interpolation energy */

                  if (E2 > 1000.0)
                    g = 1.0;
                  else
                    g = (E2 - E)/(E2 - E1);
                }
            }

          /* Check factor */

          CheckValue(FUNCTION_NAME, "g", "", g, 0.0, 1.0);

          /* Get index */

          idx = (long)RDB[loc1 + MDEP_REA_IDX];
          CheckValue(FUNCTION_NAME, "idx", "", idx, 0, 100000);

          /* Get pointer to composition */

          iso = (long)RDB[loc1 + MDEP_REA_PTR_ISO];
          CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

          /* Get density */

          adens = RDB[iso + COMPOSITION_ADENS];
          CheckValue(FUNCTION_NAME, "(adens)", "", adens, 0.0, INFTY);

          /* Score */

          if (adens0[idx] > 0.0)
            AddPrivateRes(stp + idx*ntot + ng, flx*adens*val*f*g*wgt, id);
          else
            AddPrivateRes(stp + idx*ntot + ng, flx*val*f*g*wgt, id);

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
