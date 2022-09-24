/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcfunc.c                                  */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2018/11/08 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Reads function dependent multi-physics interfaces            */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCFunc:"

/*****************************************************************************/

void ReadIFCFunc(long ifc, long update)
{
  long loc1, ptr, nmat;
  long type, n, nz, nr, np;
  double zmin, zmax, d;
  char tmpstr[MAX_STR];
  FILE *fp;

  if(update == 1)
    {
      Warn(FUNCTION_NAME, "Interface updating for type %ld not yet implemented", RDB[ifc + IFC_TYPE]);
      return;
    }

  /* Open file for reading */

  if ((fp = fopen(GetText(ifc + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(ifc, "Multi-physics interface file \"%s\" does not exist",
          GetText(ifc + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "fscanf error at function line %ld", __LINE__);

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
        Die(FUNCTION_NAME, "fscanf error");

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
    }

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

  /*******************************************************************/

  /***** Functional density distribution *****************************/

  /* Get number of parameters */

  if (fscanf(fp, "%ld", &np) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Check */

  if (np < 1)
    Error(ifc, "Invalid number of parameters");
  else
    WDB[ifc + IFC_FUNC_NP] = (double)np;

  /* Allocate memory for parameters */

  ptr = ReallocMem(DATA_ARRAY, np);

  /* Put pointer */

  WDB[ifc + IFC_FUNC_PTR_PARAM] = (double)ptr;

  /* Loop over parameters */

  for (n = 0; n < np; n++)
    {
      /* Get value */

      if (fscanf(fp, "%lf", &d) == EOF)
        Die(FUNCTION_NAME, "fscanf error");

      /* Put data */

      WDB[ptr++] = d;
    }

  /* Reset maximum density to avoid error in processinterface.c */
  /* (value is not used for anything) */

  WDB[ifc + IFC_MAX_DENSITY] = 0.0;

  /*******************************************************************/

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
