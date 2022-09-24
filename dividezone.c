/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dividezone.c                                   */
/*                                                                           */
/* Created:       2012/05/10 (JLe)                                           */
/* Last modified: 2020/05/09 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Divides material zone based on the div-card                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DivideZone:"

/*****************************************************************************/

void DivideZone(long mat, long *idx, long uni, long j, double x0, double y0,
                double z0)
{
  long div, nx, ny, nz, nrad, nseg, ntot, n, ptr, new, out, lst;
  char tmpstr[MAX_STR];

  /* Check material and universe pointers */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get pointer to divisor */

  if ((div = (long)RDB[mat + MATERIAL_PTR_DIV]) < VALID_PTR)
    return;

  /* Check cell-based separation */

  if ((long)RDB[div + DIV_SEP] == NO)
    {
      /* Check if material is associated with another universe. */
      /* (NOTE: this is a limitation in the methodology) */

      lst = (long)RDB[div + DIV_PTR_MAT_LIST];
      while (lst > VALID_PTR)
        {
          /* Get universe pointer */

          ptr = (long)RDB[lst + DIV_MAT_LIST_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check */

          if (uni != ptr)
            Error(div,
                  "Material %s located in universes %s and %s cannot be subdivided\nwithout sep entry",
                  GetText(mat + MATERIAL_PTR_NAME),
                  GetText(ptr + UNIVERSE_PTR_NAME),
                  GetText(uni + UNIVERSE_PTR_NAME));

          /* Next */

          lst = NextItem(lst);
        }

      /* Check if already defined */

      if (RDB[div + DIV_PTR_MAT_LIST] > VALID_PTR)
        return;
    }

  /* Store maximum */

  if (j > (long)RDB[div + DIV_LVL_MAX])
    WDB[div + DIV_LVL_MAX] = (double)j;

  /* Adjust level */

  if ((j = j - (long)RDB[div + DIV_SEP_LVL]) < 0)
    Error(div, "Value in \"sep\" entry exceeds number of geometry levels");

  /* Check if zone already exists */

  if ((long)RDB[div + DIV_SEP_LVL] > 0)
    {
      ptr = (long)RDB[div + DIV_PTR_MAT_LIST];
      while (ptr > VALID_PTR)
        {
          if ((long)RDB[ptr + DIV_MAT_LIST_ZONE_IDX] == idx[j])
            return;

          ptr = NextItem(ptr);
        }
    }

  /* Get bin sizes */

  nx = (long)RDB[div + DIV_NX];
  ny = (long)RDB[div + DIV_NY];
  nz = (long)RDB[div + DIV_NZ];
  nrad = (long)RDB[div + DIV_NRAD];
  nseg = (long)RDB[div + DIV_NSEG];

  /* Total */

  ntot = nx*ny*nz*nrad*nseg;

  /* Add extra zone */

  if ((long)RDB[div + DIV_LIMS_CHECK] == NO)
    ntot = ntot + 1;

  /* Allocate memory for region */

  div = NewItem(div + DIV_PTR_MAT_LIST, DIV_MAT_LIST_BLOCK_SIZE);

  /* Put zone index and universe pointer */

  WDB[div + DIV_MAT_LIST_ZONE_IDX] = (double)idx[j];
  WDB[div + DIV_MAT_LIST_PTR_UNIV] = (double)uni;

  /* Update number of zones and sub-zones (total number is updated later) */

  WDB[mat + MATERIAL_DIV_N_ZONES] = RDB[mat + MATERIAL_DIV_N_ZONES] + 1.0;
  WDB[mat + MATERIAL_DIV_N_SUB_ZONES] = (double)ntot;

  /* Get output flag */

  if ((out = (long)RDB[div + DIV_OUTPUT_FLAG]) == 0)
    out = (long)RDB[DATA_BURN_MAT_OUTPUT];

  /* Check value */

  CheckValue(FUNCTION_NAME, "ntot", "", ntot, 1, 1000000000);

  /* Allocate memory for pointers */

  ptr = ReallocMem(DATA_ARRAY, ntot);
  WDB[div + DIV_MAT_LIST_PTR_REG] = (double)ptr;

  /* Loop over bins */

  for (n = 0; n < ntot; n++)
    {
      /* Duplicate material */

      new = DuplicateItem(mat);

      /* Put origin (used for domain decomposition) */

      WDB[new + MATERIAL_DD_X0] = x0;
      WDB[new + MATERIAL_DD_Y0] = y0;
      WDB[new + MATERIAL_DD_Z0] = z0;

      /* Set divider types */

      WDB[mat + MATERIAL_DIV_TYPE] = (double)MAT_DIV_TYPE_PARENT;
      WDB[new + MATERIAL_DIV_TYPE] = (double)MAT_DIV_TYPE_NEW;

      /* Put material-wise output flags */

      if (out == BURN_OUT_MAT_DIV)
        {
          WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)NO;
          WDB[new + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
        }
      else if (out == BURN_OUT_MAT_PARENT)
        {
          WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
          WDB[new + MATERIAL_BURN_PRINT_OUTPUT] = (double)NO;
        }
      else if (out == BURN_OUT_MAT_BOTH)
        {
          WDB[mat + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
          WDB[new + MATERIAL_BURN_PRINT_OUTPUT] = (double)YES;
        }
      else
        Die(FUNCTION_NAME, "Invalid output mode");

      /* Update total number of zones */

      WDB[mat + MATERIAL_DIV_N_TOT_ZONES] =
        RDB[mat + MATERIAL_DIV_N_TOT_ZONES] + 1.0;

      /* Set zone index */

      WDB[new + MATERIAL_DIV_ZONE_IDX] = RDB[mat + MATERIAL_DIV_N_TOT_ZONES];

      /* Reset given volume */

      WDB[new + MATERIAL_VOLUME_GIVEN] = -1.0;

      /* Put parent pointer */

      WDB[new + MATERIAL_DIV_PTR_PARENT] = (double)mat;

      /* Reset composition pointer and number of zones */

      WDB[new + MATERIAL_PTR_COMP] = NULLPTR;
      WDB[new + MATERIAL_DIV_N_ZONES] = 0.0;
      WDB[new + MATERIAL_DIV_N_SUB_ZONES] = 0.0;
      WDB[new + MATERIAL_DIV_N_TOT_ZONES] = 0.0;

      /* Reset physical and majorant flags on parent */

      ResetOption(mat + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);
      ResetOption(mat + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);
      SetOption(new + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);
      SetOption(new + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);

      /* Put name */

      sprintf(tmpstr, "%sz%ld", GetText(mat + MATERIAL_PTR_NAME),
              (long)RDB[new + MATERIAL_DIV_ZONE_IDX]);
      WDB[new + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);

      /* Force DT if more than one sub-region */

      if (ntot > 1)
        {
          /* Check existing block */

          if ((long)RDB[new + MATERIAL_DT_MODE] == DT_MAT_BLOCK)
            Error(0, "Delta-tracking cannot be blocked in material %s",
                  GetText(mat + MATERIAL_PTR_NAME));

          /* Check if not full surface-tracking mode */

          if ((RDB[DATA_DT_NTHRESH] < 1.0) || (nseg > 1))
            {
              /* Enforce delta-tracking (20/04/08 2.1.31 JLe: Tämä on      */
              /* testattu toimivaksi, ja todennäköisesti on myös aika     */
              /* tehokas. Ainoa syy olla käyttämättä on DD, missä tulee   */
              /* muuten joku selittämätön bias jos alialuejako on päällä. */

              WDB[new + MATERIAL_DT_MODE] = (double)DT_MAT_FORCE;
            }
        }

      /* Put pointer in list */

      WDB[ptr++] = (double)new;
    }
}

/*****************************************************************************/
