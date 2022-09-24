/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcptavg.c                                 */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2018/11/08 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads point average multi-physics interfaces                 */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCPtAvg:"

/*****************************************************************************/

void ReadIFCPtAvg(long ifc, long update)
{
  long loc1, ptr, nmat;
  long n, nz, nr, np, dim, type;
  double xmin, xmax, ymin, ymax, zmin, zmax, Tmin, Tmax, dmax, x, y, z, d, T;
  double old, new, maxeps, maxdiff, L2abs, L2rel;
  double rad, ex;
  char tmpstr[MAX_STR];
  FILE *fp, *fout;

  if(update)
    Warn(FUNCTION_NAME, "Interface updating for type %ld not thoroughly tested", (long)RDB[ifc + IFC_TYPE]);

  /* Open file for reading */

  if ((fp = fopen(GetText(ifc + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(ifc, "Multi-physics interface file \"%s\" does not exist",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "Could not read interface type from \"%s\"",
        GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Read material name */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

  /* Allocate memory for the material array and store material name */

  if (update == NO)
    {
      /* IFC files always contain only one material name */

      nmat = 1;

      /* Store number of materials */

      WDB[ifc + IFC_N_MAT] = (double)nmat;

      /* Allocate memory for material array */

      ptr = ReallocMem(DATA_ARRAY, nmat + 1);

      /* Store pointer to material array */

      WDB[ifc + IFC_PTR_MAT_ARR] = (double)ptr;

      /* Store material name */

      WDB[ptr++] = (double)PutText(tmpstr);

      /* Null terminate the array */

      WDB[ptr] = NULLPTR;
    }

  /* Read output flag */

  if (fscanf(fp, "%ld", &n) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

  /* Store output flag */

  WDB[ifc + IFC_CALC_OUTPUT] = (double)n;

  /* Read output file name and axial binning for output */

  if (n == YES)
    {
      /* Read file name and binning for other types */

      if (fscanf(fp, "%s %ld %lf %lf %ld", tmpstr, &nz, &zmin,
                 &zmax, &nr) == EOF)
        Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

      if(!update)
        WDB[ifc + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

      /* Put number of axial zones */

      if (nz < 1)
        Error(ifc, "Error in number of axial zones");
      else
        WDB[ifc + IFC_NZ] = (double)nz;

      /* Put axial limits */

      if (zmin < zmax)
        {
          WDB[ifc + IFC_ZMIN] = zmin;
          WDB[ifc + IFC_ZMAX] = zmax;
        }
      else
        Error(ifc, "Error in axial boundaries");

      /* Put number of radial zones */

      if (nr < 1)
        Error(ifc, "Error in number of radial zones");
      else
        WDB[ifc + IFC_NR] = (double)nr;

      /* Loop over previous and check duplicate file names */

      loc1 = PrevItem(ifc);
      while (loc1 > VALID_PTR)
        {
          /* Compare file names */

          if ((long)RDB[loc1 + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
            if (CompareStr(ifc + IFC_PTR_OUTPUT_FNAME,
                           loc1 + IFC_PTR_OUTPUT_FNAME))
              Error(ifc,
                    "Duplicate output file name with distribution \"%s\"",
                    GetText(loc1 + IFC_PTR_INPUT_FNAME));

          /* Pointer to previous */

          loc1 = PrevItem(loc1);
        }
    }



  /***********************************************************************/

  /***** Average of points *******************************************/

  /* Reset limiting values */

  Tmax = -INFTY;
  Tmin = INFTY;
  dmax = 0.0;

  /* Reset mesh boundaries */

  xmin =  INFTY;
  xmax = -INFTY;
  ymin =  INFTY;
  ymax = -INFTY;
  zmin =  INFTY;
  zmax = -INFTY;

  /* Read dimension, exclusion radius and exponent (TODO: Tohon */
  /* ne muut tyypit: eksponentiaali, ym.) */

  if (fscanf(fp, "%ld %lf %lf", &dim, &rad, &ex) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

  /* Put dimension */

  if ((dim < 1) || (dim > 3))
    Error(ifc, "Error in dimension");
  else
    WDB[ifc + IFC_DIM] = (double)dim;

  /* Exclusion radius */

  if (rad > 0.0)
    WDB[ifc + IFC_EXCL_RAD] = rad;
  else
    Error(ifc, "Invalid exclusion radius entered");

  /* Exponent */

  WDB[ifc + IFC_EXP] = ex;

  /* Read number of points */

  if (fscanf(fp, "%ld", &np) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

  /* Check number of points */

  if (np > 0)
    WDB[ifc + IFC_NP] = (double)np;
  else
    Error(ifc, "Invalid number of points entered");


  /* Get pointer to points */

  if(update)
    loc1 = (long)RDB[ifc + IFC_PTR_POINTS];
  else
    loc1 = -1;

  /* Reset maximum of convergence criterion */

  maxeps  = 0.0;
  maxdiff = 0.0;
  L2abs   = 0.0;
  L2rel   = 0.0;

  /* Loop over points and read data */

  for (n = 0; n < np; n++)
    {
      /* Avoid compiler warning */

      x = 0.0;
      y = 0.0;
      z = 0.0;

      /* Read values */

      if (dim == 3)
        {
          if (fscanf(fp, "%lf %lf %lf %lf %lf",
                     &x, &y, &z, &d, &T) == EOF)
            Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);
        }
      else if (dim == 2)
        {
          if (fscanf(fp, "%lf %lf %lf %lf", &x, &y, &d, &T) == EOF)
            Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);
        }
      else if (dim == 1)
        {
          if (fscanf(fp, "%lf %lf %lf", &z, &d, &T) == EOF)
            Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);
        }
      else
        Die(FUNCTION_NAME, "error in dimension");

      /* Allocate memory for point */

      if (update == (long)NO)
        loc1 = NewItem(ifc + IFC_PTR_POINTS, IFC_PT_LIST_BLOCK_SIZE);

      /**********************************************************/
      /* Convergence criterions based on momentary distribution */
      /**********************************************************/

      if (update)
        {
          old = RDB[loc1 + IFC_PT_TMP];
          new = T;

          CalcConvCriteria(new, old,
                           &maxdiff, &maxeps, &L2abs, &L2rel);
        }

      /* Put data */

      WDB[loc1 + IFC_PT_DF] = d;
      WDB[loc1 + IFC_PT_TMP] = T;
      WDB[loc1 + IFC_PT_X] = x;
      WDB[loc1 + IFC_PT_Y] = y;
      WDB[loc1 + IFC_PT_Z] = z;

      /* Compare to limits */

      if (x - rad < xmin)
        xmin = x - rad;
      if (x + rad> xmax)
        xmax = x + rad;

      if (y - rad < ymin)
        ymin = y - rad;
      if (y + rad > ymax)
        ymax = y + rad;

      if (z - rad < zmin)
        zmin = z - rad;
      if (z + rad > zmax)
        zmax = z + rad;

      if (fabs(d) > fabs(dmax))
        dmax = d;

      if (T > Tmax)
        Tmax = T;

      if (T < Tmin)
        Tmin = T;

      /* Next point */
      if(update)
        loc1 = NextItem(loc1);
    }

  /* Put maximum density and temperature                 */
  /* For updates, these are checked in processifcptavg.c */

  WDB[ifc + IFC_MAX_DENSITY] = dmax;
  WDB[ifc + IFC_MAX_TEMP] = Tmax;
  WDB[ifc + IFC_MIN_TEMP] = Tmin;

  if(!update)
    {
      /* Put boundaries */

      WDB[ifc + IFC_MESH_XMIN] = xmin;
      WDB[ifc + IFC_MESH_XMAX] = xmax;
      WDB[ifc + IFC_MESH_YMIN] = ymin;
      WDB[ifc + IFC_MESH_YMAX] = ymax;
      WDB[ifc + IFC_MESH_ZMIN] = zmin;
      WDB[ifc + IFC_MESH_ZMAX] = zmax;

    }
  /* Set TMS on */

  if (Tmax > 0.0)
    if(!update)
      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;

  if (update)
    {

      /* Open output file for convergence */

      if (WDB[DATA_SOL_REL_ITER] == (double)0)
        {
          /* On first iteration open a new file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "w");

          /* Reset idx */

          fprintf(fout,"\nidx = 1;\n\n");
        }
      else
        {
          /* On subsequent iterations append to old file */

          sprintf(tmpstr, "%s_Tconv%ld.m", GetText(ifc + IFC_PTR_INPUT_FNAME),
                  (long)RDB[DATA_BURN_STEP]);

          fout = fopen(tmpstr, "a");

          /* Increment idx */

          fprintf(fout,"\nidx = idx + 1;\n\n");
        }

      /* Write out convergence criteria of momentary power distribution */

      fprintf(fout, "T_eps(idx) = %E;\n", maxeps);
      fprintf(fout, "T_delta(idx) = %E;\n", maxdiff);
      fprintf(fout, "T_L2_of_absolute(idx) = %E;\n", sqrt(L2abs));
      fprintf(fout, "T_L2_of_relative(idx) = %E;\n", sqrt(L2rel));

      /* Close output file */

      fclose(fout);
    }

  /*******************************************************************/

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
