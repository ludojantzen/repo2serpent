/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readpbgeometry.c                               */
/*                                                                           */
/* Created:       2010/10/22 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Reads explicit stochastic geometry                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadPBGeometry:"

/*****************************************************************************/

void ReadPBGeometry()
{  
  long pbl, loc0, msh, nx, ny, nz, nv;
  double x, y, z, r, xmin, xmax, ymin, ymax, zmin, zmax, rmax, lims[6];
  char uni[MAX_STR];
  FILE *fp;

  /* Check pointer */

  if ((loc0 = (long)RDB[DATA_PTR_PB0]) < 0)
    return;
  
  fprintf(outp, "Reading pebble bed type geometry:\n\n");
  
  /* Loop over definitions */
  
  while (loc0 > 0)
    {
      /* Test file format */
      
      TestDOSFile(GetText(loc0 + PBED_PTR_FNAME));
      
      /* Open file for reading */
      
      if ((fp = fopen(GetText(loc0 + PBED_PTR_FNAME), "r")) == NULL)
        Error(loc0, "PB geometry file \"%s\" does not exist",
              GetText(loc0 + PBED_PTR_FNAME));
      else
        fprintf(outp, "Reading data from file \"%s\"...\n",
                GetText(loc0 + PBED_PTR_FNAME));

      /* Reset limiting values */
      
      xmin =  INFTY;
      xmax = -INFTY;
      ymin =  INFTY;
      ymax = -INFTY;
      zmin =  INFTY;
      zmax = -INFTY;
      
      rmax = 0.0;
      
      /* Loop over file */
      
      while ((nv = fscanf(fp, "%lf %lf %lf %lf %s\n", &x, &y, &z, &r, uni)) 
             != EOF)
        {
          /* Check failure */

          if (nv != 5)
            Error(loc0, "Format error in file \"%s\"", 
                  GetText(loc0 + PBED_PTR_FNAME));

          /* Update counter */

          WDB[loc0 + PBED_N_PEBBLES] = RDB[loc0 + PBED_N_PEBBLES] + 1.0;

          /* Create new pebble */

          pbl = NewItem(loc0 + PBED_PTR_PEBBLES, PEBBLE_BLOCK_SIZE);
          
          /* Set values */

          WDB[pbl + PEBBLE_X0] =  x;
          WDB[pbl + PEBBLE_Y0] = y;
          WDB[pbl + PEBBLE_Z0] = z;

          WDB[pbl + PEBBLE_RAD] = r;
          WDB[pbl + PEBBLE_PTR_UNIV] = (double)PutText(uni);
          
          /* Compare to limiting values */
          
          if (x + r > xmax)
            xmax = x + r;
          if (x - r < xmin)
            xmin = x - r;
          if (y + r > ymax)
            ymax = y + r;
          if (y - r < ymin)
            ymin = y - r;
          if (z + r > zmax)
            zmax = z + r;
          if (z - r < zmin)
            zmin = z - r;
          if (r > rmax)
            rmax = r;
        }
      
      /* Calculate mesh size */
      
      nx = (long)((xmax - xmin)/(2.0*rmax));
      ny = (long)((ymax - ymin)/(2.0*rmax));
      nz = (long)((zmax - zmin)/(2.0*rmax));
      
      /* Check cut-offs */
      
      if (nx < 5)
        nx = 5;
      else if (nx > 500)
        nx = 500;
      
      if (ny < 5)
        ny = 5;
      else if (ny > 500)
        ny = 500;
      
      if (nz < 5)
        nz = 5;
      else if (nz > 500)
        nz = 500;

      /* Put mesh variables */
      
      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;
      
      /* Create mesh */
      
      msh = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR, -1, nx, ny, nz, 
                       lims, -1);

      /* Put pointer */

      WDB[loc0 + PBED_PTR_SEARCH_MESH] = (double)msh;

      /* Close file */
      
      fclose(fp);
    
      /* Next geometry */

      loc0 = NextItem(loc0);
    }

  /* Done */

  fprintf(outp, "\n");
}

/*****************************************************************************/
