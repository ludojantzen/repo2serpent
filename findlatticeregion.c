/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findlatticeregion.c                            */
/*                                                                           */
/* Created:       2010/10/08 (JLe)                                           */
/* Last modified: 2019/12/11 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Finds lattice region based on coordinates                    */
/*                                                                           */
/* Comments: - Toi lvl pointteri on nolla kun toi on super-imposed. Silloin  */
/*             ei pitäisi kirjoittaa mitään mihinkään ja regioni annetaan    */
/*             palautusarvona.                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindLatticeRegion:"

/*****************************************************************************/

long FindLatticeRegion(long lat, long lvl, double *x, double *y, double *z,
                       long *ridx, long id)
{
  long uni, ptr, reg, n, idx, nx, ny, nz, np, i, j, k, type, ncol;
  double pitch, x0, y0, z0, r2, R1, R2, phi, wdth, zmin, zmax, px, py, pz, h;

  /* Check lattice pointer */

  CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);

  /* Remember initial co-ordinates */

  x0 = *x;
  y0 = *y;
  z0 = *z;

  /* Avoid compiler warning */

  uni = -1;

  /* Reset index */

  idx = 0;

  /* Check type */

  switch (type = (long)RDB[lat + LAT_TYPE])
    {
    case LAT_TYPE_CLU:
      {
        /*********************************************************************/

        /***** Cluster type lattice ******************************************/

        /* Calculate square radius from centre */

        r2 = *x**x + *y**y;

        /* Reset radii */

        R1 = 0.0;
        R2 = 0.0;

        /* Get pointer to rings */

        ptr = (long)RDB[lat + LAT_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Loop over rings */

        n = 0;
        while ((reg = ListPtr(ptr, n)) > 0)
          {
            /* Update radii */

            R1 = R2;
            R2 = RDB[reg + RING_RLIM];

            /* Compare radii */

            if (r2 < R2*R2)
              break;

            /* Update index */

            idx = idx + (long)RDB[reg + RING_N_SEC];

            /* Next region */

            n++;
          }

        /* Check Values */

        CheckValue(FUNCTION_NAME, "r2", "", r2, R1*R1, R2*R2);
        CheckValue(FUNCTION_NAME, "idx", "", idx, 0, RDB[lat + LAT_NTOT] - 1);

        /* Calculate angle */

          phi = PolarAngle(*x, *y);

        /* Get number of sectors */

        np = (long)RDB[reg + RING_N_SEC];

        /* Calculate sector width */

        wdth = 2.0*PI/((double)np);

        /* Add tilt angle and half sector width */

        phi = phi - RDB[reg + RING_TILT] + wdth/2.0;

        /* Adjust */

        while (phi < 0.0)
          phi = phi + 2.0*PI;

        while (phi >= 2.0*PI)
          phi = phi - 2.0*PI;

        /* Get sector index */

        i = (long)((double)np*phi/(2.0*PI));
        CheckValue(FUNCTION_NAME, "i", "", i, 0, RDB[reg + RING_N_SEC] - 1);

        /* Add sector index to element index */

        idx = idx + i;
        CheckValue(FUNCTION_NAME, "idx", "", idx, 0, RDB[lat + LAT_NTOT] - 1);

        /* Sector center angle */

        phi = i*wdth + RDB[reg + RING_TILT];

        /* Adjust */

        while (phi < 0.0)
          phi = phi + 2.0*PI;

        while (phi >= 2.0*PI)
          phi = phi - 2.0*PI;

        /* Transfer co-ordinates to the lattice-element */

        *x = *x - RDB[reg + RING_RAD]*cos(phi);
        *y = *y - RDB[reg + RING_RAD]*sin(phi);

         /* Check super-imposed */

        if (lvl < VALID_PTR)
          return idx;

        /* Pointer to items */

        ptr = (long)RDB[reg + RING_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Pointer to universe */

        uni = (long)RDB[ptr + i];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        /* Put surface type */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_PAD, id);

        /* Put number of surface parameters */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 6.0, id);

        /* Put surface coordinates */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, 0.0, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, 0.0, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, R1, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, R2, id);

        if (wdth < 360.0)
          {
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4,
                           -180.0*(phi - wdth/2.0)/PI, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C5,
                           -180.0*(phi + wdth/2.0)/PI, id);
          }
        else
          {
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4, 0.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C5, 2.0*PI, id);
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_ZSTACK:
      {
        /***** Vertical stack ************************************************/

        /* Transfer co-ordinates */

        *x = *x - RDB[lat + LAT_ORIG_X0];
        *y = *y - RDB[lat + LAT_ORIG_Y0];

        /* Get number of layers */

        np = (long)RDB[lat + LAT_NTOT];

        /* Get pointer to planes */

        ptr = (long)RDB[lat + LAT_PTR_Z];
        CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

        /* Loop over stack */

        for (n = 0; n < np; n++)
          if (*z < RDB[ptr + n])
            break;

        /* Check lower boundary */

        if (n == 0)
          return -1;

        /* Get boundary values */

        zmin = RDB[ptr + n - 1];

        if (n < np)
          zmax = RDB[ptr + n];
        else
          zmax = INFTY;

        /* Check values */

        CheckValue(FUNCTION_NAME, "z", "(stack 1)", *z, zmin, zmax);

        /* Set index */

        idx = n - 1;

        /* Transfer co-ordinates to the stack-element */

        *z = *z - zmin;

        /* Check super-imposed */

        if (lvl < VALID_PTR)
          return idx;

        /* Get Universe pointer */

        uni = (long)RDB[lat + LAT_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        uni = (long)RDB[uni + idx];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        /* Put surface type */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_PZ, id);

        /* Put number of surface parameters */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 1.0, id);

        /* Put surface coordinates */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, zmin, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, zmax, id);

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_S:
    case LAT_TYPE_HX:
    case LAT_TYPE_HY:
      {
        /***** Square and hex type lattice ***********************************/

        /* Get parameters */

        nx = (long)RDB[lat + LAT_NX];
        ny = (long)RDB[lat + LAT_NY];
        pitch = RDB[lat + LAT_PITCH];

        /* Calculate size */

        np = nx*ny;

        /* Transfer co-ordinates */

        *x = *x - RDB[lat + LAT_ORIG_X0];
        *y = *y - RDB[lat + LAT_ORIG_Y0];

        /* If even number of cells, shift origin by pitch/2. This results */
        /* from the fact that the local origin is not in the centre of a  */
        /* lattice element. */

        *x = *x - (1 - (nx % 2))*0.5*pitch;
        *y = *y - (1 - (ny % 2))*0.5*pitch;

        /* Get position (i,j) in lattice */

        GetLatticeIndexes(pitch, pitch, pitch, *x, *y, 0.0, &i, &j, &k, type);

        if ((i < -nx/2) || (i > nx/2))
          return -1;
        else if ((j < -ny/2) || (j > ny/2))
          return -1;

        /* Index to lattice element */

        idx = i + (long)(nx/2.0) + (j + (long)(ny/2.0))*nx;

        /* Check value */

        if ((idx < 0) || (idx > np - 1))
          return -1;

        /* Transfer co-ordinates to lattice position */

        if (type == LAT_TYPE_S)
          {
            *x = *x - i*pitch;
            *y = *y - j*pitch;
          }
        else if (type == LAT_TYPE_HX)
          {
            *x = *x - (i + COS60*j)*pitch;
            *y = *y - j*SIN60*pitch;
          }
        else if (type == LAT_TYPE_HY)
          {
            *x = *x - j*SIN60*pitch;
            *y = *y - (i + COS60*j)*pitch;
          }
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Check super-imposed */

        if (lvl < VALID_PTR)
          return idx;

        /* Get Universe pointer */

        uni = (long)RDB[lat + LAT_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        /* Check content */

        if ((uni = (long)RDB[uni + idx]) < VALID_PTR)
          return -1;

        /* Set surface type */

        if (type == LAT_TYPE_S)
          PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_SQC, id);
        else if (type == LAT_TYPE_HX)
          PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_HEXXC, id);
        else if (type == LAT_TYPE_HY)
          PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_HEXYC, id);
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Put number of surface parameters */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 3.0, id);

        /* Put surface coordinates */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, pitch/2.0, id);

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_XTRIAG:
      {
        /***** Triangular lattice type ***************************************/

        /* Get parameters */

        nx = (long)RDB[lat + LAT_NX];
        ny = (long)RDB[lat + LAT_NY];
        h = 0.5*RDB[lat + LAT_PITCH];

        px = 3.0*h;
        py = 3.0*h/SQRT3;
        pz = INFTY;

        /* Calculate size */

        np = nx*ny;

        /* Transfer coordinates */

        *x = *x - RDB[lat + LAT_ORIG_X0] + 2*h;
        *y = *y - RDB[lat + LAT_ORIG_Y0];

        /* Check that number of elements is odd */

        if (!(nx % 2) || !(ny % 2))
          Error(lat, "Number of elements should be odd in both directions");

        /* Get position (i,j) in lattice */

        GetLatticeIndexes(px, py, pz, *x, *y, 0.0, &i, &j, &k, type);

        if ((i < -nx/2) || (i > nx/2))
          return -1;
        else if ((j < -ny/2) || (j > ny/2))
          return -1;

        /* Index to lattice element */

        idx = i + (long)(nx/2.0) + (j + (long)(ny/2.0))*nx;

        /* Check value */

        if ((idx < 0) || (idx > np - 1))
          return -1;

        /* Transfer co-ordinates to lattice position */

        if (type == LAT_TYPE_XTRIAG)
          {
            *x = *x - i*px;
            *y = *y - j*py;

            if ((i + j) % 2)
              *x = *x - h;
            else
              *x = *x  - 2*h;
          }
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Check super-imposed */

        if (lvl < VALID_PTR)
          return idx;

        /* Get Universe pointer */

        uni = (long)RDB[lat + LAT_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        /* Check content */

        if ((uni = (long)RDB[uni + idx]) < VALID_PTR)
          return -1;

        /* Set surface type */

        if (type == LAT_TYPE_XTRIAG)
          PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_TRIAG, id);
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Put number of surface parameters */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 4.0, id);

        /* Put surface coordinates */

        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, h, id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x , id);
        PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);

        /* Put type */

        if ((i + j) % 2)
          PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, 3.0, id);
        else
          PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, 1.0, id);

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_CUBOID:
    case LAT_TYPE_XPRISM:
    case LAT_TYPE_YPRISM:
      {
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

        *x = *x - RDB[lat + LAT_ORIG_X0];
        *y = *y - RDB[lat + LAT_ORIG_Y0];
        *z = *z - RDB[lat + LAT_ORIG_Z0];

        /* If even number of cells, shift origin by pitch/2. This results */
        /* from the fact that the local origin is not in the centre of a  */
        /* lattice element. */

        *x = *x - (1 - (nx % 2))*0.5*px;
        *y = *y - (1 - (ny % 2))*0.5*py;
        *z = *z - (1 - (nz % 2))*0.5*pz;

        /* Get position (i,j,k) in lattice */

        if (type == LAT_TYPE_CUBOID)
          GetLatticeIndexes(px, py, pz, *x, *y, *z, &i, &j, &k, LAT_TYPE_S);
        else if (type == LAT_TYPE_XPRISM)
          GetLatticeIndexes(px, py, pz, *x, *y, *z, &i, &j, &k, LAT_TYPE_HX);
        else if (type == LAT_TYPE_YPRISM)
          GetLatticeIndexes(px, py, pz, *x, *y, *z, &i, &j, &k, LAT_TYPE_HY);

        /* Override index if single value */

        if (nx == 1)
          i = 0;
        if (ny == 1)
          j = 0;
        if (nz == 1)
          k = 0;

        /* Index to lattice element */

        idx = i + (long)(nx/2.0) + (j + (long)(ny/2.0))*nx
          + (k + (long)(nz/2.0))*nx*ny;

        /* Check value */

        if ((idx < 0) || (idx > np - 1))
          return -1;

        /* Transfer co-ordinates to lattice position  */

        if (type == LAT_TYPE_CUBOID)
          {
            *x = *x - i*px;
            *y = *y - j*py;
            *z = *z - k*pz;
          }
        else if (type == LAT_TYPE_XPRISM)
          {
            *x = *x - (i + COS60*j)*px;
            *y = *y - j*SIN60*py;
            *z = *z - k*pz;
          }
        else if (type == LAT_TYPE_YPRISM)
          {
            *x = *x - j*SIN60*px;
            *y = *y - (i + COS60*j)*py;
            *z = *z - k*pz;
          }
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Check super-imposed */

        if (lvl < VALID_PTR)
          return idx;

        /* Get Universe pointer */

        uni = (long)RDB[lat + LAT_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        uni = (long)RDB[uni + idx];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        /* Set surface type and parameters */

        if (type == LAT_TYPE_CUBOID)
          {
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_CUBOID, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 6.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x - px/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, x0 - *x + px/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, y0 - *y - py/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, y0 - *y + py/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4, z0 - *z - pz/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C5, z0 - *z + pz/2.0, id);
          }
        else if (type == LAT_TYPE_XPRISM)
          {
            if (px != py)
              Die(FUNCTION_NAME, "px != py");

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_HEXXC, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 5.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, px/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, z0 - *z - pz/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4, z0 - *z + pz/2.0, id);
          }
        else if (type == LAT_TYPE_YPRISM)
          {
            if (px != py)
              Die(FUNCTION_NAME, "px != py");

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_HEXYC, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 5.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, py/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, z0 - *z - pz/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4, z0 - *z + pz/2.0, id);
          }
        else
          Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);

        /* Break case */

        break;

        /*********************************************************************/
      }
    case LAT_TYPE_INFS:
    case LAT_TYPE_INFHX:
    case LAT_TYPE_INFHY:
      {
        /***** Infinite 3D square or hex type lattice ************************/

        /* Put index */

        idx = 0;


        /* Get parameters */

        pitch = RDB[lat + LAT_PITCH];

        /* Transfer co-ordinates */

        *x = *x - RDB[lat + LAT_ORIG_X0];
        *y = *y - RDB[lat + LAT_ORIG_Y0];

        /* Get position (i,j,k) in lattice */

        GetLatticeIndexes(pitch, pitch, pitch, *x, *y, *z, &i, &j, &k, type);

        /* Transfer to lattice position*/

        if (type == LAT_TYPE_INFS)
          {
            *x = *x - i*pitch;
            *y = *y - j*pitch;
            *z = *z - k*pitch;
          }
        else if (type == LAT_TYPE_INFHX)
          {
            *x = *x - (i + COS60*j)*pitch;
            *y = *y - j*SIN60*pitch;
            *z = *z - k*pitch;
          }
        else if (type == LAT_TYPE_INFHY)
          {
            *x = *x - j*SIN60*pitch;
            *y = *y - (i + COS60*j)*pitch;
            *z = *z - k*pitch;
          }

        /* Check super-imposed (this does not return index to lattice-element */
        /* if super-imposed since all lattice elements are essentially        */
        /* the same) */

        if (lvl < VALID_PTR)
          return 0;

        /* Get Universe pointer */

        uni = (long)RDB[lat + LAT_PTR_FILL];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        uni = (long)RDB[uni];
        CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

        /* Set surface type and parameters */

        if (type == LAT_TYPE_INFS)
          {
            /* Put type */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_CUBE, id);

            /* Put number of surface parameters */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 4.0, id);

            /* Put surface coordinates */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, z0 - *z, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, pitch/2.0, id);
          }
        else if (type == LAT_TYPE_INFHX)
          {
            /* Put type */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_HEXXPRISM, id);

            /* Put number of surface parameters */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 5.0, id);

            /* Put surface coordinates */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, pitch/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, z0 - *z - pitch/2.0,id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4, z0 - *z + pitch/2.0,id);
          }
        else if (type == LAT_TYPE_INFHY)
          {
            /* Put type */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_TYPE, SURF_HEXYPRISM, id);

            /* Put number of surface parameters */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_NP, 5.0, id);

            /* Put surface coordinates */

            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C0, x0 - *x, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C1, y0 - *y, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C2, pitch/2.0, id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C3, z0 - *z - pitch/2.0,id);
            PutPrivateData(lvl + LVL_PRIV_LAT_SURF_C4, z0 - *z + pitch/2.0,id);
          }

        /* Break case */

        break;

        /*********************************************************************/
      }
    default:
      {
        Die (FUNCTION_NAME, "Unsupported lattice type %d.\n", type);
      }
    }

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Store index */

  StoreValuePair(lat + LAT_COL_CELL_IDX, (double)ncol, (double)idx, id);

  /* Use as region index */

  *ridx = idx;

  /* Return universe number */

  return uni;
}

/*****************************************************************************/
