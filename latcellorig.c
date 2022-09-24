/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : latcellorig.c                                  */
/*                                                                           */
/* Created:       2019/03/26 (JLe)                                           */
/* Last modified: 2019/03/26 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns coordinates of lattice cell origin.                  */
/*                                                                           */
/* Comments: - This subroutine was written to obtain the local origin        */
/*             in a multi-level geometry divided using the "div" card        */
/*             for domain decomposition.                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LatCellOrig:"

/*****************************************************************************/

void LatCellOrig(long lat, long idx, double *x, double *y, double *z)
{
  long ptr, nx, ny, nz, np, i, j, k, type;
  double pitch, x0, y0, z0, px, py, pz;

  /* Check lattice pointer */

  CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

  /* Remember initial co-ordinates */

  x0 = *x;
  y0 = *y;
  z0 = *z;

  /* Check type */

  switch (type = (long)RDB[lat + LAT_TYPE])
    {
    case LAT_TYPE_ZSTACK:
      {
        /*********************************************************************/

        /***** Vertical stack ************************************************/

        /* Transfer coordinates */

        *x = *x + RDB[lat + LAT_ORIG_X0];
        *y = *y + RDB[lat + LAT_ORIG_Y0];

        /* Get number of layers */

        np = (long)RDB[lat + LAT_NTOT];

        /* Get pointer to planes */

        ptr = (long)RDB[lat + LAT_PTR_Z];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Get local origin */

        if ((idx < 0) || (idx > np - 1))
          Die(FUNCTION_NAME, "Index error");
        else
          *z = *z + RDB[ptr + idx];

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_S:
    case LAT_TYPE_HX:
    case LAT_TYPE_HY:
      {
        /*********************************************************************/

        /***** Square or hex type lattice ************************************/

        /* Get parameters */

        nx = (long)RDB[lat + LAT_NX];
        ny = (long)RDB[lat + LAT_NY];
        pitch = RDB[lat + LAT_PITCH];

        /* Calculate size */

        np = nx*ny;

        /* Transfer co-ordinates */

        *x = *x + RDB[lat + LAT_ORIG_X0];
        *y = *y + RDB[lat + LAT_ORIG_Y0];

        /* If even number of cells, shift origin by pitch/2. This results */
        /* from the fact that the local origin is not in the centre of a  */
        /* lattice element. */

        *x = *x + (1 - (nx % 2))*0.5*pitch;
        *y = *y + (1 - (ny % 2))*0.5*pitch;

        /* Get indexes */

        j = (long)((double)idx/((double)nx));
        i = idx - j*nx;

        /* Shift */

        i = i - (long)(0.5*nx);
        j = j - (long)(0.5*ny);

        /* Check */

        if (idx != i + (long)(nx/2.0) + (j + (long)(ny/2.0))*nx)
          Die(FUNCTION_NAME, "Indexing error");

        /* Transfer to position */

        if (type == LAT_TYPE_S)
          {
            *x = *x + i*pitch;
            *y = *y + j*pitch;
          }
        else if (type == LAT_TYPE_HX)
          {
            *x = *x + (i + COS60*j)*pitch;
            *y = *y + j*SIN60*pitch;
          }
        else if (type == LAT_TYPE_HY)
          {
            *x = *x + j*SIN60*pitch;
            *y = *y + (i + COS60*j)*pitch;
          }
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_CUBOID:
    case LAT_TYPE_XPRISM:
    case LAT_TYPE_YPRISM:
      {
        /*********************************************************************/

        /***** 3D cuboidal or prismatic lattice ******************************/

        /* Get parameters */

        nx = (long)RDB[lat + LAT_NX];
        ny = (long)RDB[lat + LAT_NY];
        nz = (long)RDB[lat + LAT_NZ];

        px = RDB[lat + LAT_PITCHX];
        py = RDB[lat + LAT_PITCHY];
        pz = RDB[lat + LAT_PITCHZ];

        /* Calculate size */

        np = nx*ny*nz;

        /* Transfer co-ordinates */

        *x = *x + RDB[lat + LAT_ORIG_X0];
        *y = *y + RDB[lat + LAT_ORIG_Y0];
        *z = *z + RDB[lat + LAT_ORIG_Z0];

        /* If even number of cells, shift origin by pitch/2. This results */
        /* from the fact that the local origin is not in the centre of a  */
        /* lattice element. */

        *x = *x + (1 - (nx % 2))*0.5*px;
        *y = *y + (1 - (ny % 2))*0.5*py;
        *z = *z + (1 - (nz % 2))*0.5*pz;

        /* Get indexes */

        k = (long)((double)idx/((double)(nx*ny)));
        j = (long)((double)(idx - k*nx*ny))/((double)nx);
        i = idx - k*nx*ny - j*nx;

        /* Shift */

        i = i - (long)(0.5*nx);
        j = j - (long)(0.5*ny);
        k = k - (long)(0.5*nz);

        /* Check */

        if (idx != i + (long)(nx/2.0) + (j + (long)(ny/2.0))*nx
            + (k + (long)(nz/2.0))*nx*ny)
          Die(FUNCTION_NAME, "Indexing error");

        /* Transfer to position */

        if (type == LAT_TYPE_CUBOID)
          {
            *x = *x + i*px;
            *y = *y + j*py;
            *z = *z + j*pz;
          }
        else if (type == LAT_TYPE_XPRISM)
          {
            *x = *x + (i + COS60*j)*px;
            *y = *y + j*SIN60*py;
            *z = *z + j*pz;
          }
        else if (type == LAT_TYPE_YPRISM)
          {
            *x = *x + j*SIN60*pz;
            *y = *y + (i + COS60*j)*pz;
            *z = *z + j*pz;
          }
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Break case */

        break;

        /*********************************************************************/
      }
    default:
      {
        /*********************************************************************/

        /***** Unsupported types *********************************************/

        /* Use very negative values to indicate that the estimation */
        /* failed. */

        *x = -INFTY;
        *y = -INFTY;
        *z = -INFTY;

        /*********************************************************************/
      }
    }
}

/*****************************************************************************/
