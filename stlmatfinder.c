/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stlmatfinder.c                                 */
/*                                                                           */
/* Created:       2014/03/07 (JLe)                                           */
/* Last modified: 2016/12/07 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Used for reverse-engineering materials in STL geometries     */
/*              from corresponding CSG models                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "STLMatFinder:"

/*****************************************************************************/

void STLMatFinder()
{
  long stl, sld, msh, cell, N, n, i, j, m, id, mat;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w, d;
  unsigned long seed;
  char bname0[MAX_STR], bname[MAX_STR], fname[MAX_STR];
  signed int c;
  FILE *fp1, *fp2;

  /* Print notification */

  fprintf(outp, "\nReverse-engineering material STL compositions...\n\n");
  
  /* Set dummy OpenMP id */
  
  id = 0;

  /* Init random number sequence */
                  
  seed = ReInitRNG(0);
  SEED[id*RNG_SZ] = seed;

  /* Set plotter mode */
  
  WDB[DATA_PLOTTER_MODE] = (double)YES;

  /***************************************************************************/

  /***** Mode 1: write points in STL files into file *************************/

  if (1 == 2)
    {
      /* Set number of points */

      N = 1000;

      /* Open file for writing */

      fp1 = fopen("points.txt", "w");

      /* Loop over geometries */
      
      stl = (long)RDB[DATA_PTR_STL0];
      while (stl > VALID_PTR)
        {
          /* Get pointer to search mesh */
  
          msh = (long)RDB[stl + STL_PTR_FACET_MESH];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
          
          /* Loop over solids */

          sld = (long)RDB[stl + STL_PTR_SOLIDS];
          while (sld > VALID_PTR)
            {
              /* Print cell name */

              cell = (long)RDB[sld + STL_SOLID_PTR_CELL];
              CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
              fprintf(fp1, "%s\n", GetText(cell + CELL_PTR_NAME));

              /* Print file name and number of points */

              fprintf(fp1, "%s\n", GetText(sld + STL_SOLID_PTR_FNAME));
              fprintf(fp1, "%ld\n", N);

              /* Get boundaries */

              xmin = RDB[sld + STL_SOLID_XMIN];
              xmax = RDB[sld + STL_SOLID_XMAX];
              ymin = RDB[sld + STL_SOLID_YMIN];
              ymax = RDB[sld + STL_SOLID_YMAX];
              zmin = RDB[sld + STL_SOLID_ZMIN];
              zmax = RDB[sld + STL_SOLID_ZMAX];

              /* Adjust */

              d = 0.01*(xmax - xmin);
              xmin = xmin - d;
              xmax = xmax + d;

              d = 0.01*(ymax - ymin);
              ymin = ymin - d;
              ymax = ymax + d;

              d = 0.01*(zmax - zmin);
              zmin = zmin - d;
              zmax = zmax + d;
              
              /* Print (for progress) */
              
              fprintf(outp, "%s : %E %E : %E %E : %E %E\n", 
                      GetText(sld + STL_SOLID_PTR_FNAME), xmin, xmax, 
                      ymin, ymax, zmin, zmax);

              /* Sample points */

              n = 0;
              i = 0;
              j = 0;

              do
                {
                  /* Sample point */

                  x = RandF(id)*(xmax - xmin) + xmin;
                  y = RandF(id)*(ymax - ymin) + ymin;
                  z = RandF(id)*(zmax - zmin) + zmin;

                  /* Sample direction */
              
                  IsotropicDirection(&u, &v, &w, id);

                  /* Perform ray test */
                  
                  if (STLRayTest(sld, msh, x, y, z, u, v, w, 
                                 STL_SEARCH_MODE_FAST, id) == YES)
                    {
                      /* Print coordinates */

                      fprintf(fp1, "%11.5E %11.5E %11.5E\n", x, y, z);
                      
                      /* Update count */

                      n++;
                    }

#ifdef mmmmmmmmmmm

                  /* Get pointer to cell */

                  if ((ptr = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                    if (ptr == cell)
                      {
                        /* Print coordinates */
                        
                        fprintf(fp1, "%16.10E %16.10E %16.10E\n", x, y, z);
                        
                        /* Update count */
                        
                        n++;
                      }

#endif

                  /* Check if failures */

                  if (j++ == 1000000)
                    {
                      Die(FUNCTION_NAME, "Infinite loop");
                      
                      /* Adjust */

                      d = 0.20*(xmax - xmin);
                      xmin = xmin - d;
                      xmax = xmax + d;
                      
                      d = 0.20*(ymax - ymin);
                      ymin = ymin - d;
                      ymax = ymax + d;
                      
                      d = 0.20*(zmax - zmin);
                      zmin = zmin - d;
                      zmax = zmax + d;
                    }
                }
              while (n < N);
              
              /* Next solid */
              
              sld = NextItem(sld);
            }
          
          /* Next geometry */

          stl = NextItem(stl);
        }

      /* Close file */

      fclose(fp1);
    }
  
  /***************************************************************************/

  /***** Mode 2: find corresponding materials in another geometry ************/

  else
    {
      /* Open file for reading */

      fp1 = fopen("points.txt", "r");
      fp2 = fopen("geometry.txt", "w");

      /* Reset previous body name and count */

      bname[0] = '\0';
      m = 0;

      /* Loop over input */

      while (fscanf(fp1, "%s\n", bname) != EOF)
        {
          /* Compare body name to previous and reset count */

          if (strcmp(bname, bname0))
            m = 1;
          else
            m++;

          /* Remember previous */

          strcpy(bname0, bname);

          /* Read file name */
          
          n = 0;
          while ((c = fgetc(fp1)) != '\n')
            fname[n++] = (char)c;

          fname[n] = '\0';

          /* Read number of points */

          if (fscanf(fp1, "%ld", &N) == EOF)
            Die(FUNCTION_NAME, "fscanf error 2");

          /* Reset material-wise counters (use FLOW_IDX) */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Reset count */

              WDB[mat + MATERIAL_FLOW_IDX] = 0.0;

              /* Next material */

              mat = NextItem(mat);
            }

          /* Print (for progress) */
              
          fprintf(outp, "%s\n", fname);

          /* Loop over points */

          for (n = 0; n < N; n++)
            {
              /* Read point */

              if (fscanf(fp1, "%lf %lf %lf\n", &x, &y, &z) == EOF)
                Die(FUNCTION_NAME, "fscanf error 2");

              /* Sample direction */
              
              IsotropicDirection(&u, &v, &w, id);

              /* Get pointer to cell */

              if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
                {
                  /* Add to material count */

                  if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
                    WDB[mat + MATERIAL_FLOW_IDX] = 
                      RDB[mat + MATERIAL_FLOW_IDX] + 1.0;
                }
            }
          
          /* Reset count */

          n = 0.0;

          /* Loop over materials and print */

          mat = (long)RDB[DATA_PTR_M0];
          while (mat > VALID_PTR)
            {
              /* Check count */

              if (RDB[mat + MATERIAL_FLOW_IDX] > 0.0)
                {
                  /* Print */

                  fprintf(fp2, "body %s%ld %s%ld %s", bname, m, bname, m,
                          GetText(mat + MATERIAL_PTR_NAME));

                  /* Print fraction */
                  
                  fprintf(fp2, " %% %4.1f%%\n", 
                          100.0*(RDB[mat + MATERIAL_FLOW_IDX]/((double)N)));
                          
                  /* Add count */

                  n++;
                }

              /* Next material */

              mat = NextItem(mat);
            }

          /* Check count */
          
          if (n == 0)
            {
              
              /* Print */
              
              fprintf(fp2, "body %s%ld %s%ld %s\n", bname, m, bname, m, "void");
            }

          fprintf(fp2, "file %s%ld \"%s\" 1.0 0 0 0\n\n", bname, m, fname);
               }

      /* Close files */

      fclose(fp1);
      fclose(fp2);     
    }

  /***************************************************************************/

  /* Terminate run */

  exit(0);

  /***************************************************************************/
}
