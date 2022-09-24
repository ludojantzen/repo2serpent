/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcfblims.c                                */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Reads fuel behavior multi-physics interfaces output limits   */
/*                                                                           */
/* Comments:   -Split from readifcfb.c for 2.1.25                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCFBLims:"

/*****************************************************************************/

void ReadIFCFBLims(FILE *fp, long loc1, long update)
{
  long loc2, ptr;
  long nz, na, nr, n;
  double zmin, zmax, amin, amax, rmin, rmax, emin, emax;

  /***********************/
  /* Read output meshing */
  /***********************/

  /* Allocate memory for power limits or get pointer from memory */

  if (!update)
    {
      /* Allocate memory */

      loc2 = ReallocMem(DATA_ARRAY, FUEP_LIM_SIZE);

      /* Store pointer */

      WDB[loc1 + IFC_FUEP_OUT_PTR_LIM] = (double)loc2;
    }
  else
    loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];

  /* Read number of axial zones */

  if (fscanf(fp, "%ld", &nz) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Store number of axial zones */

  WDB[loc2 + FUEP_NZ] = (double)nz;

  /* Read axial limits if equally spaced zones */

  if (nz > 0)
    {
      if(fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
      WDB[loc2 + FUEP_ZMIN] = zmin;

      CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);
      WDB[loc2 + FUEP_ZMAX] = zmax;

    }
  else if (nz < 0)
    {
      WDB[loc2 + FUEP_ZMIN] = 0.0;
      WDB[loc2 + FUEP_ZMAX] = 0.0;

      /* Get number of axial zones to be read */

      nz = -1*nz;

      /* Store number of axial zones */

      WDB[loc2 + FUEP_NZ] = (double)nz;

      /* Read axial meshing with a dedicated routine */

      ptr = ReadIFCBins(nz, fp, -INFTY/10.0, update);

      /* Store pointer or get it from memory */

      if (!update)
        WDB[loc1 + IFC_FUEP_OUT_PTR_Z] = (double)ptr;
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_Z];

      /* Store min/max coordinates */

      WDB[loc2 + FUEP_ZMIN] = RDB[ptr +  0];
      WDB[loc2 + FUEP_ZMAX] = RDB[ptr + nz];
    }
  else
    Die(FUNCTION_NAME, "Zero axial output zones");

  /* Read number of angular zones */

  if (fscanf(fp, "%ld", &na) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc2 + FUEP_NA] = (double)na;

  /* Read angular limits if equally spaced zones */

  if (na > 0)
    {
      if(fscanf(fp, "%lf %lf", &amin, &amax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      CheckValue(FUNCTION_NAME, "amin", "", amin, -360.0, 360.0);
      WDB[loc2 + FUEP_AMIN] = amin;

      CheckValue(FUNCTION_NAME, "amax", "", amax, amin, 360.0);
      WDB[loc2 + FUEP_AMAX] = amax;

    }
  else if (na < 0)
    {
      WDB[loc2 + FUEP_AMIN] = 0.0;
      WDB[loc2 + FUEP_AMAX] = 0.0;

      /* Get number of zones to be read */

      na = -1*na;

      /* Store number of zones */

      WDB[loc2 + FUEP_NA] = (double)na;

      /* Read zone limits with a dedicated routine */

      ptr = ReadIFCBins(na, fp, -360.0, update);

      /* Store pointer or get it from memory */

      if (!update)
        WDB[loc1 + IFC_FUEP_OUT_PTR_PHI] = (double)ptr;
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_PHI];

      /* Store min/max coordinates */

      WDB[loc2 + FUEP_AMIN] = RDB[ptr +  0];
      WDB[loc2 + FUEP_AMAX] = RDB[ptr + na];
    }
  else
    Die(FUNCTION_NAME, "Zero angular output zones");

  /* Read number of radial zones */

  if (fscanf(fp, "%ld", &nr) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc2 + FUEP_NR] = (double)nr;

  /* Read radial limits if equal area zones */

  if (nr > 0)
    {
      if(fscanf(fp, "%lf %lf", &rmin, &rmax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      CheckValue(FUNCTION_NAME, "rmin", "", rmin, 0.0, INFTY);
      WDB[loc2 + FUEP_RMIN] = rmin;

      CheckValue(FUNCTION_NAME, "amax", "", rmax, rmin, INFTY);
      WDB[loc2 + FUEP_RMAX] = rmax;

    }
  else if (nr < 0)
    {
      WDB[loc2 + FUEP_RMIN] = 0.0;
      WDB[loc2 + FUEP_RMAX] = 0.0;

      /* Get number of zones to be read */

      nr = -1*nr;

      /* Store number of zones */

      WDB[loc2 + FUEP_NR] = (double)nr;

      /* Read zone limits with a dedicated routine */

      ptr = ReadIFCBins(nr, fp, 0.0, update);

      /* Store pointer or get it from memory */

      if(!update)
        WDB[loc1 + IFC_FUEP_OUT_PTR_R2] = (double)ptr;
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_R2];

      /* Radial limits must be lifted to the second power */

      if (!update)
        for (n = 0; n < nr + 1; n++)
          WDB[ptr + n] = RDB[ptr + n]*RDB[ptr + n];

      /* Store min/max coordinates */

      WDB[loc2 + FUEP_RMIN] = RDB[ptr +  0];
      WDB[loc2 + FUEP_RMAX] = RDB[ptr + nr];

    }
  else
    Die(FUNCTION_NAME, "Zero radial output zones");

  /* Read output meshing (flux) */

  if(!update)
    {
      loc2 = ReallocMem(DATA_ARRAY, FUEP_LIM_SIZE);
      WDB[loc1 + IFC_FUEP_OUT_PTR_FLIM] = (double)loc2;
    }
  else
    loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];

  /* Read number of axial zones */

  if (fscanf(fp, "%ld", &nz) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc2 + FUEP_NZ] = (double)nz;

  /* Read axial limits if equally spaced zones */

  if (nz > 0)
    {
      if(fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
      WDB[loc2 + FUEP_ZMIN] = zmin;

      CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);
      WDB[loc2 + FUEP_ZMAX] = zmax;

    }
  else if (nz < 0)
    {
      WDB[loc2 + FUEP_ZMIN] = 0.0;
      WDB[loc2 + FUEP_ZMAX] = 0.0;

      /* Get number of axial zones to be read */

      nz = -1*nz;

      /* Store number of axial zones */

      WDB[loc2 + FUEP_NZ] = (double)nz;

      /* Read axial meshing with a dedicated routine */

      ptr = ReadIFCBins(nz, fp, -INFTY/10.0, update);

      /* Store pointer or get it from memory */

      if (!update)
        WDB[loc1 + IFC_FUEP_OUT_PTR_FZ] = (double)ptr;
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FZ];

      /* Store min/max coordinates */

      WDB[loc2 + FUEP_ZMIN] = RDB[ptr +  0];
      WDB[loc2 + FUEP_ZMAX] = RDB[ptr + nz];
    }
  else
    Die(FUNCTION_NAME, "Zero axial output zones");

  /* Read number of angular zones */

  if (fscanf(fp, "%ld", &na) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc2 + FUEP_NA] = (double)na;

  /* Read angular limits if equally spaced zones */

  if (na > 0)
    {
      if(fscanf(fp, "%lf %lf", &amin, &amax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      CheckValue(FUNCTION_NAME, "amin", "", amin, -360.0, 360.0);
      WDB[loc2 + FUEP_AMIN] = amin;

      CheckValue(FUNCTION_NAME, "amax", "", amax, amin, 360.0);
      WDB[loc2 + FUEP_AMAX] = amax;

    }
  else if (na < 0)
    {
      WDB[loc2 + FUEP_AMIN] = 0.0;
      WDB[loc2 + FUEP_AMAX] = 0.0;

      /* Get number of zones to be read */

      na = -1*na;

      /* Store number of zones */

      WDB[loc2 + FUEP_NA] = (double)na;

      /* Read zone limits with a dedicated routine */

      ptr = ReadIFCBins(na, fp, -360.0, update);

      /* Store pointer or get it from memory */

      if (!update)
        WDB[loc1 + IFC_FUEP_OUT_PTR_PHI] = (double)ptr;
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_PHI];

      /* Store min/max coordinates */

      WDB[loc2 + FUEP_AMIN] = RDB[ptr +  0];
      WDB[loc2 + FUEP_AMAX] = RDB[ptr + na];
    }
  else
    Die(FUNCTION_NAME, "Zero angular output zones");

  /* Read number of radial zones */

  if (fscanf(fp, "%ld", &nr) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc2 + FUEP_NR] = (double)nr;

  /* Read radial limits if equal area zones */

  if (nr > 0)
    {
      if(fscanf(fp, "%lf %lf", &rmin, &rmax) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      CheckValue(FUNCTION_NAME, "rmin", "", rmin, 0.0, INFTY);
      WDB[loc2 + FUEP_RMIN] = rmin;

      CheckValue(FUNCTION_NAME, "amax", "", rmax, rmin, INFTY);
      WDB[loc2 + FUEP_RMAX] = rmax;

    }
  else if (nr < 0)
    {
      WDB[loc2 + FUEP_RMIN] = 0.0;
      WDB[loc2 + FUEP_RMAX] = 0.0;

      /* Get number of radial bins to read */

      nr = -1*nr;

      /* Store number of radial bins */

      WDB[loc2 + FUEP_NR] = (double)nr;

      /* Read binning with dedicated routine */

      ptr = ReadIFCBins(nr, fp, 0.0, update);

      /* Store pointer to binning or get it from memory */

      if (!update)
        WDB[loc1 + IFC_FUEP_OUT_PTR_FR2] = (double)ptr;
      else
        ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FR2];

      /* Radial limits must be lifted to the second power */

      if (!update)
        for (n = 0; n < nr + 1; n++)
          WDB[ptr + n] = RDB[ptr + n]*RDB[ptr + n];

      /* Set radial limits */

      WDB[loc2 + FUEP_RMIN] = RDB[ptr +  0];
      WDB[loc2 + FUEP_RMAX] = RDB[ptr + nr];
    }
  else
    Die(FUNCTION_NAME, "Zero radial output zones");

  /* Read fast flux energy limits */

  if (fscanf(fp, "%lf %lf", &emin, &emax) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc2 + FUEP_EMIN] = emin;
  WDB[loc2 + FUEP_EMAX] = emax;

}
