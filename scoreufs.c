/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreufs.c                                     */
/*                                                                           */
/* Created:       2012/04/22 (JLe)                                           */
/* Last modified: 2019/03/13 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Scores reaction rates needed for UFS weight mesh             */
/*                                                                           */
/* Comments: - Ton lattice-optionk kanssa ei voi käyttää suoraan tracking-   */
/*             rutiinin tallentamaa indeksiä koska se laskee symmetria-      */
/*             tapauksissa leikkaantuneet kopit väärin.                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreUFS:"

/*****************************************************************************/

void ScoreUFS(double flx, long mat, double wgt, double x, double y, double z,
              double E, double g, long id)
{
  long msh, lat, rea, mode, type, nx, ny, nz, i, j, k;
  double f, zmin, zmax, pr, pz;

  /* Check that mode is on */

  if (((mode = (long)RDB[DATA_UFS_MODE]) == UFS_MODE_NONE) ||
      (mode == UFS_MODE_RMX))
    return;

  /* Check active cycles */

  if (RDB[DATA_CYCLE_IDX] > RDB[DATA_CRIT_SKIP])
    return;

  /* Get pointer to mesh */

  msh = (long)RDB[DATA_UFS_PTR_SRC_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Reset value */

  f = -1.0;

  /* Check mode */

  if (mode == UFS_MODE_FLUX)
    {
      /* Flux */
      
      f = flx*wgt;      
    }
  else if (mat > VALID_PTR)
    {
      /* Reaction rate */
      
      if (mode == UFS_MODE_COL)
        {
          /* Total (collision) */      
          
          if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTXS]) > VALID_PTR)
            f = MacroXS(rea, E, id)*flx*wgt*g;
        }
      else if (mode == UFS_MODE_FISS)
        {
          /* Fission */      
          
          if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
            f = MacroXS(rea, E, id)*flx*wgt*g;
        }
      else
        Die(FUNCTION_NAME, "Invalid mode");
    }

  /* Check value */

  if (f < ZERO)
    return;

  /* Check pointer to lattice */

  if ((lat = (long)RDB[DATA_UFS_PTR_LAT]) > VALID_PTR)
    {
      /* Get lattice type */
      
      type = (long)RDB[lat + LAT_TYPE];

      /* Check */

      if ((type != LAT_TYPE_S) && (type != LAT_TYPE_HX) && 
          (type != LAT_TYPE_HY))
        Die(FUNCTION_NAME, "Lattice type %ld not supported with UFS", type);

      /* Get lattice size */
      
      nx = (long)RDB[lat + LAT_NX];
      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 10000);

      ny = (long)RDB[lat + LAT_NY];
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 10000);

      /* Get axial dimension */
      
      nz = (long)RDB[DATA_UFS_NZ];
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 10000);

      zmin = RDB[DATA_UFS_ZMIN];
      CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);

      zmax = RDB[DATA_UFS_ZMAX];
      CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

      /* Get radial and axial pitch */

      pr = RDB[lat + LAT_PITCH];
      pz = (zmax - zmin)/((double)nz);

      CheckValue(FUNCTION_NAME, "pr", "", pr, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "pz", "", pz, ZERO, INFTY);

      /* Transfer co-ordinates */
      
      x = x - RDB[lat + LAT_ORIG_X0];
      y = y - RDB[lat + LAT_ORIG_Y0];

      /* Get indexes */
      
      GetLatticeIndexes(pr, pr, pz, x, y, z, &i, &j, &k, type);

      /* Adjust */

      i = i + (long)((double)nx/2.0);
      j = j + (long)((double)ny/2.0);

      /* Check */

      if ((i < 0) || (i > nx - 1) || (j < 0) || (j > ny - 1) ||
          (k < 0) || (k > nz - 1))
        return;

      /* Score source point for UFS */
      
      AddMeshIdx(msh, f, i, j, k, id);
    }
  else
    {
      /* Score source point for UFS in Cartesian mesh */
      
      AddMesh(msh, f, x, y, z, id);
    }
}

/*****************************************************************************/
