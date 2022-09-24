/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : divideburnmat.c                                */
/*                                                                           */
/* Created:       2011/07/03 (JLe)                                           */
/* Last modified: 2019/11/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Divides burnable materials in multiple rings etc.            */
/*                                                                           */
/* Comments: - This subroutine divides burnable materials in pin-type        */
/*             regions into multiple zones and rings, as is done in          */
/*             Serpent 1. Division of more complicated 3D structures using   */
/*             the "div"-cards is carried out in makedepletionzones.c        */
/*                                                                           */
/*           - Tää nopeutui merkittävästi kun ton haun erotti omaksi         */
/*             osakseen (tehokkaampi välimuistin käyttö?)                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DivideBurnMat:"

/*****************************************************************************/

void DivideBurnMat()
{
  long nst, type, reg0, mat, nr, n, m, surf, loc0, loc1, loc2, ptr, nptr, reg;
  long reg1, n0;
  double r0, r1, r;
  char name[MAX_STR];

  /* Check if burnable materials are defined in input */

  if ((long)RDB[DATA_BURN_MATERIALS_FLAG] == NO)
    return;
  else
    fprintf(outp, "Dividing burnable materials...\n");

  /***************************************************************************/

  /***** Find pointers *******************************************************/

  /* Loop over nests */

  nst = (long)RDB[DATA_PTR_NST0];
  while (nst > VALID_PTR)
    {
      /* Get nest type */

      type = (long)RDB[nst + NEST_TYPE];

      /* Get pointer if simple types, otherwise set to null */

      if ((type == SURF_CYL) || (type == SURF_CYLX) || (type == SURF_CYLY) ||
          (type == SURF_CYLZ) || (type == SURF_SPH))
        reg0 = (long)RDB[nst + NEST_PTR_REGIONS];
      else
        reg0 = -1;

      /* Loop over regions */

      while (reg0 > VALID_PTR)
        {
          /* Find material (cell pointer is used as temporary storage space) */

          if ((long)RDB[reg0 + NEST_REG_PTR_CELL] > VALID_PTR)
            {
              mat = (long)RDB[DATA_PTR_M0];
              mat = SeekListStr(mat, MATERIAL_PTR_NAME,
                                GetText(reg0 + NEST_REG_PTR_CELL));
            }
          else
            mat = -1;

          /* Put pointer */

          WDB[reg0 + NEST_REG_TMP_PTR] = (double)mat;

          /* Next region */

          reg0 = NextItem(reg0);
        }

      /* Next nest */

      nst = NextItem(nst);
    }

  /***************************************************************************/

  /***** Divide **************************************************************/

  /* Loop over nests */

  nst = (long)RDB[DATA_PTR_NST0];
  while (nst > VALID_PTR)
    {
      /* Get nest type */

      type = (long)RDB[nst + NEST_TYPE];

      /* Get pointer if simple types, otherwise set to null */

      if ((type == SURF_CYL) || (type == SURF_CYLX) || (type == SURF_CYLY) ||
          (type == SURF_CYLZ) || (type == SURF_SPH))
        reg0 = (long)RDB[nst + NEST_PTR_REGIONS];
      else
        reg0 = -1;

      /* Reset radii */

      r0 = 0.0;
      r1 = 0.0;

      /* Reset number of extra regions */

      n0 = 0;

      /* Loop over regions */

      while (reg0 > VALID_PTR)
        {
          /* Pointer to surface */

          if ((surf = (long)RDB[reg0 + NEST_REG_PTR_SURF_IN]) > VALID_PTR)
            {
              /* Pointer to parameter list */

              loc0 = (long)RDB[surf + SURFACE_PTR_PARAMS];
              CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

              /* Get radius */

              if (type == SURF_SPH)
                r1 = RDB[loc0 + 3];
              else
                r1 = RDB[loc0 + 2];
            }
          else
            r1 = INFTY;

          /* Get material pointer */

          mat = (long)RDB[reg0 + NEST_REG_TMP_PTR];

          /* Get number of burnable zones */

          if (mat > VALID_PTR)
            nr = (long)RDB[mat + MATERIAL_BURN_RINGS];
          else
            nr = 0;

          /* Check if material is associated with a divisor */

          if ((mat < VALID_PTR) ||
              ((long)RDB[mat + MATERIAL_PTR_DIV] > VALID_PTR))
            {
              /* Next region */

              reg0 = NextItem(reg0);

              /* Cycle loop */

              continue;
            }

          /* Region pointer */

          reg = reg0;

          /* Loop over regions */

          for (n = 0; n < nr; n++)
            {
              /* Material name */

              sprintf(name, "%sp%sr%ld", GetText(mat + MATERIAL_PTR_NAME),
                      GetText(nst + NEST_PTR_NAME), n0 + n + 1);

              /* Find match in previous */

              reg1 = PrevItem(reg);
              while (reg1 > VALID_PTR)
                {
                  /* Compare material names */

                  if ((long)WDB[reg1 + NEST_REG_PTR_CELL] > VALID_PTR)
                    if (!strcmp(name, GetText(reg1 + NEST_REG_PTR_CELL)))
                      break;

                  /* Previous */

                  reg1 = PrevItem(reg1);
                }

              /* Check pointer */

              if (reg1 > VALID_PTR)
                {
                  /* Regions reg1 and reg contain same material, update */
                  /* number of extra regions */

                  n0++;

                  /* Rename */

                  sprintf(name, "%sp%sr%ld", GetText(mat + MATERIAL_PTR_NAME),
                          GetText(nst + NEST_PTR_NAME), n0 + n + 1);
                }

              /* Put name */

              nptr = PutText(name);

              /* Calculate radius (ei toimi annulaarisilla) */

              if (type == SURF_SPH)
                r = cbrt(((double)(n + 1))/((double)nr))*(r1 - r0) + r0;
              else
                r = sqrt(((double)(n + 1))/((double)nr))*(r1 - r0) + r0;

              /* Duplicate material */

              loc1 = DuplicateItem(mat);

              /* Set divider types */

              WDB[mat + MATERIAL_DIV_TYPE] = (double)MAT_DIV_TYPE_PARENT;
              WDB[loc1 + MATERIAL_DIV_TYPE] = (double)MAT_DIV_TYPE_S1;

              /* Put output flags */

              if ((long)RDB[DATA_BURN_MAT_OUTPUT] == BURN_OUT_MAT_DIV)
                {
                  WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)NO;
                  WDB[loc1 + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
                }
              else if ((long)RDB[DATA_BURN_MAT_OUTPUT] == BURN_OUT_MAT_PARENT)
                {
                  WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
                  WDB[loc1 + MATERIAL_BURN_PRINT_OUTPUT] = (double)NO;
                }
              else if ((long)RDB[DATA_BURN_MAT_OUTPUT] == BURN_OUT_MAT_BOTH)
                {
                  WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
                  WDB[loc1 + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
                }
              else
                Die(FUNCTION_NAME, "Invalid output mode");

              /* Update total number of zones */

              WDB[mat + MATERIAL_DIV_N_TOT_ZONES] =
                RDB[mat + MATERIAL_DIV_N_TOT_ZONES] + 1.0;

              /* Set zone index */

              WDB[loc1 + MATERIAL_DIV_ZONE_IDX] =
                WDB[mat + MATERIAL_DIV_N_TOT_ZONES];

              /* Set number of zones and sub-zones */

              WDB[mat + MATERIAL_DIV_N_SUB_ZONES] = (double)nr;
              WDB[mat + MATERIAL_DIV_N_ZONES] =
                RDB[mat + MATERIAL_DIV_N_TOT_ZONES]/((double)nr);

              /* Put parent pointer */

              WDB[loc1 + MATERIAL_DIV_PTR_PARENT] = (double)mat;

              /* Reset composition pointer and number of zones */

              WDB[loc1 + MATERIAL_PTR_COMP] = NULLPTR;
              WDB[loc1 + MATERIAL_DIV_N_ZONES] = 0.0;
              WDB[loc1 + MATERIAL_DIV_N_SUB_ZONES] = 0.0;
              WDB[loc1 + MATERIAL_DIV_N_TOT_ZONES] = 0.0;

              /* Put name pointer */

              WDB[loc1 + MATERIAL_PTR_NAME] = (double)nptr;

              /* Duplicate region */

              reg = DuplicateItem(reg);

              /* Check surface pointer */

              if (surf < VALID_PTR)
                WDB[reg + NEST_REG_PTR_SURF_IN] = NULLPTR;
              else
                {
                  /* Duplicate surface */

                  loc2 = DuplicateItem(surf);

                  /* Put pointer */

                  WDB[reg + NEST_REG_PTR_SURF_IN] = (double)loc2;

                  /* Get number of surface parameters */

                  m = (long)RDB[loc2 + SURFACE_N_PARAMS];

                  /* Create parameter list */

                  ptr = ReallocMem(DATA_ARRAY, m);

                  /* Put pointer */

                  WDB[loc2 + SURFACE_PTR_PARAMS] = (double)ptr;

                  /* Put radius (tää ei nyt sit kopioi noita muita) */

                  WDB[ptr + m - 1] = r;
                }

              /* Put name pointer */

              WDB[reg + NEST_REG_PTR_CELL] = (double)nptr;

              /* Put material pointer */

              WDB[reg + NEST_REG_TMP_PTR] = (double)loc1;
            }

          /* Remove original region */

          if (nr > 0)
            RemoveItem(reg0);

          /* Set burn flag */

          SetOption(nst + NEST_OPTIONS, OPT_BURN_MAT);

          /* Update radius */

          r0 = r1;

          /* Next region */

          reg0 = NextItem(reg);
        }

      /* Next nest */

      nst = NextItem(nst);
    }

  /***************************************************************************/

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
