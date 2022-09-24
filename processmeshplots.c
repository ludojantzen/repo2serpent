/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processmeshplots.c                             */
/*                                                                           */
/* Created:       2011/03/20 (JLe)                                           */
/* Last modified: 2019/04/03 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Processes mesh plots                                         */
/*                                                                           */
/* Comments: - Tässä varaillaan muistia vähän tarpeettomasti jos tyyppi on   */
/*             1-arvoinen jakauma.                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessMeshPlots:"

/*****************************************************************************/

void ProcessMeshPlots()
{
  long mpl, nx, ny, ptr, n, det, ax, type, ndis, div, sz;
  double xmin, xmax, ymin, ymax, zmin, zmax, lims[6];
  char tmpstr[MAX_STR];

  /* Check if mesh plots are defined */

  if ((long)RDB[DATA_PTR_MPL0] < VALID_PTR)
    return;

  fprintf(outp, "Processing mesh plots...\n");

  /* Loop over plots and calculate total size */

  sz = 0;

  mpl = (long)RDB[DATA_PTR_MPL0];
  while(mpl > VALID_PTR)
    {
      /* Get size */

      nx = (long)RDB[mpl + MPL_NX];
      ny = (long)RDB[mpl + MPL_NY];

      /* Get number of distributions and divider flag */

      ndis = (long)RDB[mpl + MPL_NDIST];
      div = (long)WDB[mpl + MPL_DIV];

      /* Add to size */

      if (div == NO)
        sz = sz + ndis*(nx*ny + 1);
      else
        sz = sz + 2*ndis*(nx*ny + 1);

      /* Next plot */

      mpl = NextItem(mpl);
    }

  /* Allocate memory for data */

  PreallocMem(sz + 1, RES2_ARRAY);

  /* Reset counter */

  n = 1;

  /* Loop over plots */

  mpl = (long)RDB[DATA_PTR_MPL0];
  while(mpl > VALID_PTR)
    {
      /* Get limits */

      if ((xmin = RDB[mpl + MPL_XMIN]) == -INFTY)
        xmin = RDB[DATA_GEOM_MINX];

      if ((xmax = RDB[mpl + MPL_XMAX]) == INFTY)
        xmax= RDB[DATA_GEOM_MAXX];

      if ((ymin = RDB[mpl + MPL_YMIN]) == -INFTY)
        ymin = RDB[DATA_GEOM_MINY];

      if ((ymax = RDB[mpl + MPL_YMAX]) == INFTY)
        ymax= RDB[DATA_GEOM_MAXY];

      if ((zmin = RDB[mpl + MPL_ZMIN]) == -INFTY)
        zmin = RDB[DATA_GEOM_MINZ];

      if ((zmax = RDB[mpl + MPL_ZMAX]) == INFTY)
        zmax= RDB[DATA_GEOM_MAXZ];

      /* Get type */

      type = (long)RDB[mpl + MPL_TYPE];

      /* Get number of distributions and divider flag */

      ndis = (long)RDB[mpl + MPL_NDIST];
      div = (long)WDB[mpl + MPL_DIV];

      /* Get axis and set dimensions */

      ax = (long)RDB[mpl + MPL_AX];

      /* Check is and set coordinates */

      if (ax == 1)
        {
          /* Distribution in yz-plane */

          WDB[mpl + MPL_XMIN] = ymin;
          WDB[mpl + MPL_XMAX] = ymax;

          WDB[mpl + MPL_YMIN] = zmin;
          WDB[mpl + MPL_YMAX] = zmax;

          WDB[mpl + MPL_ZMIN] = xmin;
          WDB[mpl + MPL_ZMAX] = xmax;
        }
      else if (ax == 2)
        {
          /* Distribution in xz-plane */

          WDB[mpl + MPL_XMIN] = xmin;
          WDB[mpl + MPL_XMAX] = xmax;

          WDB[mpl + MPL_YMIN] = zmin;
          WDB[mpl + MPL_YMAX] = zmax;

          WDB[mpl + MPL_ZMIN] = ymin;
          WDB[mpl + MPL_ZMAX] = ymax;
        }
      else if (ax == 3)
        {
          /* Distribution in xy-plane */

          WDB[mpl + MPL_XMIN] = xmin;
          WDB[mpl + MPL_XMAX] = xmax;

          WDB[mpl + MPL_YMIN] = ymin;
          WDB[mpl + MPL_YMAX] = ymax;

          WDB[mpl + MPL_ZMIN] = zmin;
          WDB[mpl + MPL_ZMAX] = zmax;
        }
      else if (ax == 4)
        {
          /* Cylindrical mesh (tää olettaa että luetaan 4 parametria) */

          if ((xmin < 0.0) && (xmin == RDB[DATA_GEOM_MINX]))
            Error(mpl, "Inner radius must be provided");
          else if (xmin < 0.0)
            Error(mpl, "Invalid value %1.5E for inner radius", xmin);
          else if ((xmax <= 0.0) && (xmax == RDB[DATA_GEOM_MAXX]))
            Error(mpl, "Outer radius must be provided");
          else if (xmax <= 0.0)
            Error(mpl, "Invalid value %1.5E for outer radius", xmax);
          else if ((zmin != RDB[DATA_GEOM_MINZ]) ||
                   (zmax != RDB[DATA_GEOM_MAXZ]))
            Error(mpl, "Only r- and z-boundies should be given");

          WDB[mpl + MPL_XMIN] = xmin;
          WDB[mpl + MPL_XMAX] = xmax;

          WDB[mpl + MPL_YMIN] = 0.0;
          WDB[mpl + MPL_YMAX] = 2.0*PI;

          WDB[mpl + MPL_ZMIN] = ymin;
          WDB[mpl + MPL_ZMAX] = ymax;

          /* Check that symmetry is not used */

          if ((long)RDB[mpl + MPL_SYM] > 0)
            Error(mpl, "Symmetry not allowed with cylindrical mesh");
        }

      /* Get size */

      nx = (long)RDB[mpl + MPL_NX];
      ny = (long)RDB[mpl + MPL_NY];

      /* Put mesh variables */

      lims[0] = RDB[mpl + MPL_XMIN];
      lims[1] = RDB[mpl + MPL_XMAX];
      lims[2] = RDB[mpl + MPL_YMIN];
      lims[3] = RDB[mpl + MPL_YMAX];
      lims[4] = RDB[mpl + MPL_ZMIN];
      lims[5] = RDB[mpl + MPL_ZMAX];

      /* Allocate memory */

      if (ax < 4)
        ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1, nx, ny, 1,
                         lims, -1);
      else
        ptr = CreateMesh(MESH_TYPE_CYLINDRICAL, MESH_CONTENT_RES, -1, nx, 1,
                         ny, lims, -1);

      /* Put pointer */

      WDB[mpl + MPL_PTR_VAL1] = (double)ptr;

      /* Check divider flag */

      if (div == YES)
        {
          /* Allocate memory */

          ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1, nx, ny,
                           1, lims, -1);

          /* Put pointer */

          WDB[mpl + MPL_PTR_DIV1] = (double)ptr;
        }

      /* Check double distribution */

      if (ndis == 2)
        {
          /* Allocate memory */

          ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1, nx, ny,
                           1, lims, -1);

          /* Put pointer */

          WDB[mpl + MPL_PTR_VAL2] = (double)ptr;

          /* Check divider flag */

          if (div == YES)
            {
              /* Allocate memory */

              ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1,
                               nx, ny, 1, lims, -1);

              /* Put pointer */

              WDB[mpl + MPL_PTR_DIV2] = (double)ptr;
            }
        }

      /* Set file name */

      sprintf(tmpstr, "%s_mesh%ld", GetText(DATA_PTR_INPUT_FNAME), n++);
      WDB[mpl + MPL_PTR_FNAME] = (double)PutText(tmpstr);

      /* Link detector */

      if ((type == MPL_TYPE_DET) ||
          (type == MPL_TYPE_DET_IMP))
        {
          /* Find detector */

          det = (long)RDB[DATA_PTR_DET0];
          CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

          if ((det = SeekListStr(det, DET_PTR_NAME,
                                 GetText(mpl + MPL_PTR_DET))) < VALID_PTR)
            Error(mpl, "Detector %s not defined", GetText(mpl + MPL_PTR_DET));

          /* Get pointer to response function */

          ptr = (long)RDB[det + DET_PTR_RBINS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check response function */

          if (((long)RDB[ptr + DET_RBIN_MT] == MT_PHOTON_PULSE_HEIGHT) ||
              ((long)RDB[ptr + DET_RBIN_MT] == MT_MACRO_HEATPHOTANA))
            Error(mpl, "Detector response %ld not allowed with mesh plots",
                  (long)RDB[ptr + DET_RBIN_MT]);

          /* Check for more responses */

          if (NextItem(ptr) > VALID_PTR)
            Error(mpl,
                  "Multiple response functions not allowed with mesh plots");

          /* Set pointer */

          WDB[mpl + MPL_PTR_DET] = (double)det;
        }

      /* Check importance detector */

      if (type == MPL_TYPE_DET_IMP)
        {
          /* Record events */

          SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_IMPORTANCE);
        }

      /* Next plot */

      mpl = NextItem(mpl);
    }

  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
