/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readstlgeometry.c                              */
/*                                                                           */
/* Created:       2014/03/03 (JLe)                                           */
/* Last modified: 2017/11/39 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Reads STL type geometry data                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadSTLGeometry()"

/*****************************************************************************/

void ReadSTLGeometry()
{
  short count;
  int tmpi;
  float tmpf;
  long loc0, loc1, loc2, sld, nf, ptr, n, m, i, N, ***arr;
  char fname[MAX_STR], str[MAX_STR], sname[MAX_STR], hdr[81];
  double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax, val;
  FILE *fp;

  /* Check pointer */

  if ((loc0 = (long)RDB[DATA_PTR_STL0]) < VALID_PTR)
    return;

  fprintf(outp, "Reading stereolithography (STL) geometries...\n");

  /* Loop over definitions */

  while (loc0 > VALID_PTR)
    {
      /***********************************************************************/
  
      /***** First loop to get geometry boundaries ***************************/

      /* Reset boundaries */

      xmin = INFTY;
      xmax = -INFTY;
      ymin = INFTY;
      ymax = -INFTY;
      zmin = INFTY;
      zmax = -INFTY;

      /* Loop over files */

      loc1 = (long)RDB[loc0 + STL_PTR_FILES];
      while (loc1 > VALID_PTR)
        {
          /* Get file name */

          sprintf(fname, "%s", GetText(loc1 + STL_FILE_PTR_FNAME));

          /* Open file for reading */

          if ((fp = fopen(fname, "r")) == NULL)
            Error(loc0, "STL geometry file \"%s\" does not exist", fname);
          
          /* Check file type */

          if (TestASCIIFile(fname) == NO)
            {
              /****************************************************************/

              /***** Binary format ********************************************/

              /* Read header */

              if (fread(hdr, sizeof(char), 80, fp) == 0)
                Error(loc0, "Format error in STL file");

              /* Put EOF */

              hdr[80] = '\0';

              /* Read number of facets */
              
              if (fread(&tmpi, 4, 1, fp) == 0)
                Error(loc0, "Format error in STL file");

              /* Convert to long */

              if ((nf = (long)tmpi) < 1)
                Error(loc0, "No facets in STL file \"%s\"", fname);

              /* Loop over vertices */

              for (n = 0; n < nf; n++)
                {
                  /* Read direction vector */

                  for (m = 0; m < 3; m++)
                    {
                      /* Read value */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");
                    }

                  /* Read points */

                  for (m = 0; m < 3; m++)
                    {
                      /* Read x-coordinate */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      x = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);

                      /* Read y-coordinate */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      y = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);

                      /* Read z-coordinate */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      z = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);
                      
                      /* Scale and shift origin */

                      x = x*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_X0];
                      y = y*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Y0];
                      z = z*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Z0];

                      /* Compare to boundaries */

                      if (x < xmin)
                        xmin = x;
                      if (x > xmax)
                        xmax = x;

                      if (y < ymin)
                        ymin = y;
                      if (y > ymax)
                        ymax = y;

                      if (z < zmin)
                        zmin = z;
                      if (z > zmax)
                        zmax = z;
                    }

                  /* Read attribute byte count */

                  if (fread(&count, 2, 1, fp) == 0)
                    Error(loc0, "Format error in STL file");
                }

              /****************************************************************/
            }
          else
            {
              /****************************************************************/

              /***** ASCII format *********************************************/

              /* Check file format */

              TestDOSFile(fname);

              /* Read data */

              while(fscanf(fp, "%s", str) != EOF)
                {
                  /* Check input string */

                  if (!strcmp(str, "facet"))
                    {
                      /* Read normal vector */

                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (1)");
                      else if (strcmp(str, "normal"))
                        Error(loc0, "Format error in STL file (2)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (3)");

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (4)");
                      else if (strcmp(str, "outer"))
                        Error(loc0, "Format error in STL file (5)");

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (6)");
                      else if (strcmp(str, "loop"))
                        Error(loc0, "Format error in STL file (7)");

                      /* Read first point */
                      
                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (8)");
                      else if (strcmp(str, "vertex"))
                        Error(loc0, "Format error in STL file (9)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (10)");

                      /* Scale and shift origin */

                      x = x*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_X0];
                      y = y*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Y0];
                      z = z*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Z0];

                      /* Compare to boundaries */

                      if (x < xmin)
                        xmin = x;
                      if (x > xmax)
                        xmax = x;

                      if (y < ymin)
                        ymin = y;
                      if (y > ymax)
                        ymax = y;

                      if (z < zmin)
                        zmin = z;
                      if (z > zmax)
                        zmax = z;

                      /* Read second point */
                      
                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (11)");
                      else if (strcmp(str, "vertex"))
                        Error(loc0, "Format error in STL file (12)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (13)");

                      /* Scale and shift origin */

                      x = x*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_X0];
                      y = y*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Y0];
                      z = z*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Z0];

                      /* Compare to boundaries */

                      if (x < xmin)
                        xmin = x;
                      if (x > xmax)
                        xmax = x;

                      if (y < ymin)
                        ymin = y;
                      if (y > ymax)
                        ymax = y;

                      if (z < zmin)
                        zmin = z;
                      if (z > zmax)
                        zmax = z;

                      /* Read third point */
                      
                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (14)");
                      else if (strcmp(str, "vertex"))
                        Error(loc0, "Format error in STL file (15)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (16)");

                      /* Scale and shift origin */

                      x = x*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_X0];
                      y = y*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Y0];
                      z = z*RDB[loc1 + STL_FILE_SCALING] - 
                        RDB[loc1 + STL_FILE_Z0];

                      /* Compare to boundaries */

                      if (x < xmin)
                        xmin = x;
                      if (x > xmax)
                        xmax = x;

                      if (y < ymin)
                        ymin = y;
                      if (y > ymax)
                        ymax = y;

                      if (z < zmin)
                        zmin = z;
                      if (z > zmax)
                        zmax = z;

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (17)");
                      else if (strcmp(str, "endloop"))
                        Error(loc0, "Format error in STL file (18)");

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (19)");
                      else if (strcmp(str, "endfacet"))
                        Error(loc0, "Format error in STL file (20)");
                    }
                }

              /****************************************************************/
            }

          /* Close file */

          fclose(fp);

          /* Next solid */

          loc1 = NextItem(loc1);
        }

      /* Store boundaries */

      WDB[loc0 + STL_XMIN] = xmin;
      WDB[loc0 + STL_XMAX] = xmax;
      WDB[loc0 + STL_YMIN] = ymin;
      WDB[loc0 + STL_YMAX] = ymax;
      WDB[loc0 + STL_ZMIN] = zmin;
      WDB[loc0 + STL_ZMAX] = zmax;

      /* Allocate memory for pointer array */
          
      N = (long)RDB[DATA_STL_TEMP_ARRAY_SIZE];

      arr = (long ***)Mem(MEM_ALLOC, N, sizeof(long **));        
          
      for(n = 0; n < N; n++)
        {
          arr[n] = (long **)Mem(MEM_ALLOC, N, sizeof(long *));
          for(m = 0; m < N; m++)
            arr[n][m] = (long *)Mem(MEM_ALLOC, N, sizeof(long));
        }

      /* Reset pointers */

      for (n = 0; n < N; n++)
        for (m = 0; m < N; m++)
          for (i = 0; i < N; i++)
            arr[n][m][i] = NULLPTR;

      /***********************************************************************/

      /***** Second loop to read data ****************************************/

      /* Loop over files */

      loc1 = (long)RDB[loc0 + STL_PTR_FILES];
      while (loc1 > VALID_PTR)
        {
          /* Get file name */

          sprintf(fname, "%s", GetText(loc1 + STL_FILE_PTR_FNAME));

          /* Open file for reading */

          if ((fp = fopen(fname, "r")) == NULL)
            Die(FUNCTION_NAME, "Unable to open file");
          
          /* Check file type */

          if (TestASCIIFile(fname) == NO)
            {
              /****************************************************************/

              /***** Binary format ********************************************/

              /* Read header */

              if (fread(hdr, sizeof(char), 80, fp) == 0)
                Error(loc0, "Format error in STL file");

              /* Put EOF */

              hdr[80] = '\0';

              /* Create structure */

              sld = NewItem(loc0 + STL_PTR_SOLIDS, STL_SOLID_BLOCK_SIZE);

              /* Put name */
                      
              WDB[sld + STL_SOLID_PTR_STL_NAME] = 
                RDB[loc1 + STL_FILE_PTR_SNAME];

              /* Put file name */

              WDB[sld + STL_SOLID_PTR_FNAME] = 
                RDB[loc1 + STL_FILE_PTR_FNAME];

              /* Read number of facets */
              
              if (fread(&tmpi, 4, 1, fp) == 0)
                Error(loc0, "Format error in STL file");

              /* Convert to long */

              if ((nf = (long)tmpi) < 1)
                Error(loc0, "No facets in STL file \"%s\"", fname);

              /* Allocate memory for facets */

              ptr = ReallocMem(DATA_ARRAY, nf*STL_FACET_BLOCK_SIZE);
              WDB[sld + STL_SOLID_PTR_FACETS] = (double)ptr;

              /* Put number of facets */

              WDB[sld + STL_SOLID_N_FACETS] = (double)nf;

              /* Loop over vertices */

              for (n = 0; n < nf; n++)
                {
                  /* Read direction vector */

                  for (m = 0; m < 3; m++)
                    {
                      /* Read value */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      val = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "val", "", val, -1.0, 1.0);

                      /* Put data */

                      WDB[ptr + STL_FACET_NORM_U + m] = val;
                    }

                  /* Read points */

                  for (m = 0; m < 3; m++)
                    {
                      /* Read x-coordinate */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      x = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);

                      /* Read y-coordinate */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      y = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);

                      /* Read z-coordinate */
                      
                      if (fread(&tmpf, 4, 1, fp) == 0)
                        Error(loc0, "Format error in STL file");

                      /* Convert to double */

                      z = (double)tmpf;
                      CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);

                      /* Store point */

                      loc2 = AddSTLPoint(arr, loc0, loc1, sld, x, y, z);
                      WDB[ptr + STL_FACET_PTR_PT1 + m] = (double)loc2;
                    }

                  /* Read attribute byte count */

                  if (fread(&count, 2, 1, fp) == 0)
                    Error(loc0, "Format error in STL file");

                  /* Store pointer to solid */

                  WDB[ptr + STL_FACET_PTR_SOLID] = (double)sld;
                
                  /* Update pointer */

                  ptr = ptr + STL_FACET_BLOCK_SIZE;
                }

              /****************************************************************/
            }
          else
            {
              /****************************************************************/

              /***** ASCII format *********************************************/

              /* Check file format */

              TestDOSFile(fname);

              /* Loop over file and count facets */

              nf = 0;
              while(fscanf(fp, "%s", str) != EOF)
                if (!strcmp(str, "facet"))
                  nf++;

              /* Check */

              if (nf < 1)
                Error(loc0, "No facets in STL file \"%s\"", fname);

              /* Rewind file */

              rewind(fp);

              /* Allocate memory for facets */

              ptr = ReallocMem(DATA_ARRAY, nf*STL_FACET_BLOCK_SIZE);
              
              /* Reset solid name and pointer */
              
              sname[0] = '\0';
              sld = -1;

              /* Read data */

              while(fscanf(fp, "%s", str) != EOF)
                {
                  /* Check input string */

                  if (!strcmp(str, "solid"))
                    {
                      /* Create structure */

                      sld = NewItem(loc0 + STL_PTR_SOLIDS, 
                                    STL_SOLID_BLOCK_SIZE);

                      /* Read name */
                      
                      if (fscanf(fp, "%s", sname) == EOF)

                        Error(loc0, "Format error in STL file (21)");

                      /* Put name */

                      if (strcmp(GetText(loc1 + STL_FILE_PTR_SNAME), "-1"))
                        WDB[sld + STL_SOLID_PTR_STL_NAME] = 
                          RDB[loc1 + STL_FILE_PTR_SNAME];
                      else
                        WDB[sld + STL_SOLID_PTR_STL_NAME] = 
                          (double)PutText(sname);

                      /* Put file name */

                      WDB[sld + STL_SOLID_PTR_FNAME] = 
                        RDB[loc1 + STL_FILE_PTR_FNAME];
                      
                      /* Put facet pointer */

                      WDB[sld + STL_SOLID_PTR_FACETS] = (double)ptr;

                      /* Reset number of facets */

                      nf = 0;
                    }
                  else if (!strcmp(str, "endsolid"))
                    {
                      /* Check name */

                      if (sname[0] == '\0')
                        Error(loc0, "Format error in STL file (22)");

                      /* Put number of facets */

                      WDB[sld + STL_SOLID_N_FACETS] = (double)nf;
                    }
                  else if (!strcmp(str, "facet"))
                    {
                      /* Check pointer */

                      if (sld < VALID_PTR)
                        Error(loc0, "Format error in STL file (25)");

                      /* Reset variables */

                      x = -INFTY;
                      y = -INFTY;
                      z = -INFTY;

                      /* Read normal vector */

                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (26)");
                      else if (strcmp(str, "normal"))
                        Error(loc0, "Format error in STL file (27)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (28)");

                      /* Store */

                      WDB[ptr + STL_FACET_NORM_U] = x;
                      WDB[ptr + STL_FACET_NORM_V] = y;
                      WDB[ptr + STL_FACET_NORM_W] = z;

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (29)");
                      else if (strcmp(str, "outer"))
                        Error(loc0, "Format error in STL file (30)");

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (31)");
                      else if (strcmp(str, "loop"))
                        Error(loc0, "Format error in STL file (32)");

                      /* Reset variables */

                      x = -INFTY;
                      y = -INFTY;
                      z = -INFTY;

                      /* Read first point */
                      
                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (33)");
                      else if (strcmp(str, "vertex"))
                        Error(loc0, "Format error in STL file (34)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (35)");

                      /* Store */
                      
                      loc2 = AddSTLPoint(arr, loc0, loc1, sld, x, y, z);
                      WDB[ptr + STL_FACET_PTR_PT1] = (double)loc2;

                      /* Reset variables */

                      x = -INFTY;
                      y = -INFTY;
                      z = -INFTY;

                      /* Read second point */
                      
                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (36)");
                      else if (strcmp(str, "vertex"))
                        Error(loc0, "Format error in STL file (37)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (38)");

                      /* Store */
                      
                      loc2 = AddSTLPoint(arr, loc0, loc1, sld, x, y, z);
                      WDB[ptr + STL_FACET_PTR_PT2] = (double)loc2;

                      /* Reset variables */

                      x = -INFTY;
                      y = -INFTY;
                      z = -INFTY;

                      /* Read third point */
                      
                      if (fscanf(fp, "%s %lf %lf %lf", str, &x, &y, &z) == EOF)
                        Error(loc0, "Format error in STL file (39)");
                      else if (strcmp(str, "vertex"))
                        Error(loc0, "Format error in STL file (40)");
                      else if ((x == -INFTY) || (y == -INFTY) || (z == -INFTY))
                        Error(loc0, "Format error in STL file (41)");

                      /* Store */
                      
                      loc2 = AddSTLPoint(arr, loc0, loc1, sld, x, y, z);
                      WDB[ptr + STL_FACET_PTR_PT3] = (double)loc2;

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (42)");
                      else if (strcmp(str, "endloop"))
                        Error(loc0, "Format error in STL file (43)");

                      /* Check next word */

                      if (fscanf(fp, "%s", str) == EOF)
                        Error(loc0, "Format error in STL file (44)");
                      else if (strcmp(str, "endfacet"))
                        Error(loc0, "Format error in STL file (45)");

                      /* Add counter */

                      nf++;

                      /* Store pointer to solid */

                      WDB[ptr + STL_FACET_PTR_SOLID] = (double)sld;

                      /* Update pointer */

                      ptr = ptr + STL_FACET_BLOCK_SIZE;
                    }
                }

              /***************************************************************/
            }

          /* Close file */

          fclose(fp);

          /* Next solid */

          loc1 = NextItem(loc1);
        }

      /* Free pointer array */
          
      for(n = 0; n < N; n++)
        {
          for(m = 0; m < N; m++)
            Mem(MEM_FREE, arr[n][m]);
          Mem(MEM_FREE, arr[n]);
        } 

      Mem(MEM_FREE, arr);

      /***********************************************************************/

      /* Next */

      loc0 = NextItem(loc0);
    }
  
  /* Exit OK */
  
  fprintf(outp, "OK.\n\n");
}

/*****************************************************************************/
