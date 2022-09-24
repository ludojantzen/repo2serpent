/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : disperse.c                                     */
/*                                                                           */
/* Created:       2012/02/08 (JLe)                                           */
/* Last modified: 2015/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Generates random particle distributions inside simple        */
/*              volumes                                                      */
/*                                                                           */
/* Comments: - Grow and shake algorithm and parallelepiped volume            */
/*             implemented by Cole Gentry.                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Disperse:"

/* NOTE: pallo ja kuutio on origokeskeisiä. */

/*****************************************************************************/

struct particle {
  long idx;
  double x;
  double y;
  double z;
  double crp;
  double rp;
  long u;
  long fixed;
  struct particle *next;
};

struct smsh {
  long sz;
  long *n;
};

/*****************************************************************************/

void Disperse()
{
  long type, N, Nprev, N0, Nmax, u, n, m, i, j, k, nx, ny, nz;
  long imin, imax, jmin, jmax, kmin, kmax, ni, id, count, *temp;
  long  flag, small;
  unsigned long seed;
  double r, rmin, rpmax, xmin, xmax, ymin, ymax, zmin, zmax, rp, x, y, z;
  double xlen, ylen, zlen, La, Lb, Lc, ps, th, ph;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3, a4, b4, c4, d4, a5;
  double b5, c5, d5, a6, b6, c6, d6;
  double dx, dy, dz, dr, V0, V1, V2, pf, rpfracmin, rpfracprev, crp, shake;
  double growth;
  char answer[256], fname[256];
  struct particle *p, *first, *p1; 
  struct smsh ***f;
  FILE *fp, *outfp;

  /***************************************************************************/

  /***** Prompt input data ***************************************************/

  /* Check if input file is given */

  if ((long)RDB[DATA_PTR_INPUT_FNAME] > 0)
    {
      fp = fopen(GetText(DATA_PTR_INPUT_FNAME), "r");
      outfp = fopen("/dev/null", "w");
    }
  else
    {
      fp = stdin;
      outfp = stdout;
    }

  printf("\nRandom particle distribution file generator launched...\n");

  if (fp == NULL)
    {
      printf("\nUnable to open input file.\n\n");
      exit(-1);
    }
  else if (fp != stdin)
    printf("\nReading data from file...\n\n");

  /* Get volume type */

  fprintf(outfp, "\nEnter volume type: 1 = sphere\n");
  fprintf(outfp, "                   2 = cylinder\n");
  fprintf(outfp, "                   3 = cube\n");
  fprintf(outfp, "                   4 = annular cylinder\n");
  fprintf(outfp, "                   5 = cuboid\n");
  fprintf(outfp, "                   6 = parallelepiped\n\n");

  if (fscanf(fp, "%ld", &type) == EOF)
    Error(0, "Missing type");

  fprintf(outfp, "\n");

  /* This is needed to avoid compiler warnings */

  a1 = 0.0;
  b1 = 0.0;
  c1 = 0.0;
  d1 = 0.0;
  a2 = 0.0;
  b2 = 0.0;
  c2 = 0.0;
  d2 = 0.0;
  a3 = 0.0;
  b3 = 0.0;
  c3 = 0.0;
  d3 = 0.0;
  a4 = 0.0;
  b4 = 0.0;
  c4 = 0.0;
  d4 = 0.0;
  a5 = 0.0;
  b5 = 0.0;
  c5 = 0.0;
  d5 = 0.0;
  a6 = 0.0;
  b6 = 0.0;
  c6 = 0.0;
  d6 = 0.0;
  ymin = 0.0;
  ymax = 0.0;
  xmin = 0.0;
  xmax = 0.0;
  V0 = 0.0;  

  p = NULL;
  temp= NULL;

  /* Get volume dimensions */

  if (type == 1)
    {
      fprintf(outfp, "Enter sphere radius (cm): ");

      if (fscanf(fp, "%lf", &r) == EOF)
        Error(0, "Missing radius");

      V0 = (4.0/3.0)*PI*r*r*r;

      xmin = -r;
      xmax =  r;
      ymin = -r;
      ymax =  r;
      zmin = -r;
      zmax =  r;
    }
  else if (type == 2)
    {
      fprintf(outfp, "Enter cylinder radius (cm): ");

      if (fscanf(fp, "%lf", &r) == EOF)
        Error(0, "Missing radius");

      fprintf(outfp, "Enter cylinder bottom coordinate (cm): ");

      if (fscanf(fp, "%lf", &zmin) == EOF)
        Error(0, "Missing bottom coordinate");

      fprintf(outfp, "Enter cylinder top coordinate (cm): ");

      if (fscanf(fp, "%lf", &zmax) == EOF)
        Error(0, "Missing top coordinate");

      V0 = PI*r*r*(zmax - zmin);

      xmin = -r;
      xmax =  r;
      ymin = -r;
      ymax =  r;
    }
  else if (type == 3)
    {
      fprintf(outfp, "Enter cube half-width (cm): ");

      if (fscanf(fp, "%lf", &r) == EOF)
        Error(0, "Missing half width");

      V0 = 8*r*r*r;

      xmin = -r;
      xmax =  r;
      ymin = -r;
      ymax =  r;
      zmin = -r;
      zmax =  r;
    }
  else if (type == 4)
    {
      fprintf(outfp, "Enter cylinder radius (cm): ");

      if (fscanf(fp, "%lf", &r) == EOF)
        Error(0, "Missing radius");

      fprintf(outfp, "Enter cylinder bottom coordinate (cm): ");

      if (fscanf(fp, "%lf", &zmin) == EOF)
        Error(0, "Missing bottom coordinate");

      fprintf(outfp, "Enter cylinder top coordinate (cm): ");

      if (fscanf(fp, "%lf", &zmax) == EOF)
        Error(0, "Missing top coordinate");

      fprintf(outfp, "Enter central hole radius (cm): ");

      if (fscanf(fp, "%lf", &rmin) == EOF)
        Error(0, "Missing central hole radius");

      V0 = PI*(r*r - rmin*rmin)*(zmax - zmin);

      xmin = -r;
      xmax =  r;
      ymin = -r;
      ymax =  r;
    }
  else if (type == 5)
    {
      fprintf(outfp, "Enter length in x-direction (cm): ");

      if (fscanf(fp, "%lf", &xlen) == EOF)
        Error(0, "Missing x-length");

      fprintf(outfp, "Enter length in y-direction (cm): ");

      if (fscanf(fp, "%lf", &ylen) == EOF)
        Error(0, "Missing y-length");

      fprintf(outfp, "Enter length in z-direction (cm): ");

      if (fscanf(fp, "%lf", &zlen) == EOF)
        Error(0, "Missing z-length");

      V0 = xlen*ylen*zlen;

      xmin = -xlen/2;
      xmax =  xlen/2;
      ymin = -ylen/2;
      ymax =  ylen/2;
      zmin = -zlen/2;
      zmax =  zlen/2;
    }
  else if (type == 6)
    {
      fprintf(outfp, "Enter length of La (cm): ");

      if (fscanf(fp, "%lf", &La) == EOF)
        Error(0, "Missing La");

      fprintf(outfp, "Enter length of Lb (cm): ");

      if (fscanf(fp, "%lf", &Lb) == EOF)
        Error(0, "Missing Lb");

      fprintf(outfp, "Enter length of Lc (cm): ");

      if (fscanf(fp, "%lf", &Lc) == EOF)
        Error(0, "Missing Lc");

      fprintf(outfp, "Enter Psi angle (degrees): ");

      if (fscanf(fp, "%lf", &ps) == EOF)
        Error(0, "Missing psi");

      fprintf(outfp, "Enter Theta angle (degrees): ");

      if (fscanf(fp, "%lf", &th) == EOF)
        Error(0, "Missing theta");

      fprintf(outfp, "Enter Phi angle (degrees): ");

      if (fscanf(fp, "%lf", &ph) == EOF)
        Error(0, "Missing phi");

      ps = ps*PI/180.0;
      th = th*PI/180.0;
      ph = ph*PI/180.0;

      /*Calculate the parallelepiped volume*/

      V0 = La*cos(th)*(Lb*Lc*cos(ps));

      /*Determine the maximum possible x, y, z values */

      xmax =  La + (Lb + Lc*(sin(th)*sin(ph)/cos(ps)))*sin(ps) + 
        Lc*(sin(th)*cos(ph));
      ymax =  Lb*cos(ps) + Lc*(sin(th)*sin(ph));
      zmax =  Lc*cos(th);

      /*Calculate the plane equations for each surface of the parallelepiped*/

      x1 = 0.0;
      y1 = 0.0;
      z1 = 0.0;
      x2 = La;
      y2 = 0.0;
      z2 = 0.0;
      x3 = La + (Lc*(sin(th)*sin(ph)/cos(ps)))*sin(ps) + Lc*(sin(th)*cos(ph));
      y3 = Lc*(sin(th)*sin(ph));
      z3 = zmax;

      a1 = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      b1 = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      c1 = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      d1 = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);

      /* Check Surface 1 Coefficients*/

      if (fabs(a1*x1 + b1*y1 + c1*z1 + d1) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (1) surface (1): %E",
            fabs(a1*x1 + b1*y1 + c1*z1 + d1));

      if (fabs(a1*x2 + b1*y2 + c1*z2 + d1) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (2) surface (1): %E",
          fabs(a1*x2 + b1*y2 + c1*z2 + d1));

      if (fabs(a1*x3 + b1*y3 + c1*z3 + d1) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (3) surface (1): %E",
          fabs(a1*x3 + b1*y3 + c1*z3 + d1));

      x1 = Lb*sin(ps);
      y1 = Lb*cos(ps);
      z1 = 0.0;
      x2 = La + Lb*sin(ps);
      y2 = Lb*cos(ps);
      z2 = 0.0;
      x3 = xmax;
      y3 = ymax;
      z3 = zmax;

      a2 = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      b2 = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      c2 = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      d2 = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);

      /* Check Surface 2 Coefficients*/

      if (fabs(a2*x1 + b2*y1 + c2*z1 + d2) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (1) surface (1): %E",
            fabs(a2*x1 + b2*y1 + c2*z1 + d2));

      if (fabs(a2*x2 + b2*y2 + c2*z2 + d2) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (2) surface (1): %E",
          fabs(a2*x2 + b2*y2 + c2*z2 + d2));

      if (fabs(a2*x3 + b2*y3 + c2*z3 + d2) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (3) surface (1): %E",
          fabs(a2*x3 + b2*y3 + c2*z3 + d2));

      x1 = 0.0;
      y1 = 0.0;
      z1 = 0.0;
      x2 = (Lc*(sin(th)*sin(ph)/cos(ps)))*sin(ps)+Lc*(sin(th)*cos(ph));
      y2 = Lc*(sin(th)*sin(ph));
      z2 = zmax;
      x3 = Lb*sin(ps);
      y3 = Lb*cos(ps);
      z3 = 0.0;

      a3 = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      b3 = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      c3 = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      d3 = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);

      /* Check Surface 3 Coefficients*/

      if (fabs(a3*x1 + b3*y1 + c3*z1 + d3) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (1) surface (1): %E",
            fabs(a3*x1 + b3*y1 + c3*z1 + d3));

      if (fabs(a3*x2 + b3*y2 + c3*z2 + d3) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (2) surface (1): %E",
          fabs(a3*x2 + b3*y2 + c3*z2 + d3));

      if (fabs(a3*x3 + b3*y3 + c3*z3 + d3) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (3) surface (1): %E",
          fabs(a3*x3 + b3*y3 + c3*z3 + d3));

      x1 = La;
      y1 = 0.0;
      z1 = 0.0;
      x2 = La+(Lc*(sin(th)*sin(ph)/cos(ps)))*sin(ps)+Lc*(sin(th)*cos(ph));
      y2 = Lc*(sin(th)*sin(ph));
      z2 = zmax;
      x3 = La+Lb*sin(ps);
      y3 = Lb*cos(ps);
      z3 = 0.0;

      a4 = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      b4 = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      c4 = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      d4 = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);

      /* Check Surface 4 Coefficients*/

      if (fabs(a4*x1 + b4*y1 + c4*z1 + d4) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (1) surface (1): %E",
            fabs(a4*x1 + b4*y1 + c4*z1 + d4));

      if (fabs(a4*x2 + b4*y2 + c4*z2 + d4) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (2) surface (1): %E",
          fabs(a4*x2 + b4*y2 + c4*z2 + d4));

      if (fabs(a4*x3 + b4*y3 + c4*z3 + d4) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (3) surface (1): %E",
          fabs(a4*x3 + b4*y3 + c4*z3 + d4));

      x1 = 0.0;
      y1 = 0.0;
      z1 = 0.0;
      x2 = La;
      y2 = 0.0;
      z2 = 0.0;
      x3 = La+Lb*sin(ps);
      y3 = Lb*cos(ps);
      z3 = 0.0;

      a5 = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      b5 = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      c5 = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      d5 = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);

      /* Check Surface 5 Coefficients*/

      if (fabs(a5*x1 + b5*y1 + c5*z1 + d5) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (1) surface (1): %E",
            fabs(a5*x1 + b5*y1 + c5*z1 + d5));

      if (fabs(a5*x2 + b5*y2 + c5*z2 + d5) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (2) surface (1): %E",
          fabs(a5*x2 + b5*y2 + c5*z2 + d5));

      if (fabs(a5*x3 + b5*y3 + c5*z3 + d5) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (3) surface (1): %E",
          fabs(a5*x3 + b5*y3 + c5*z3 + d5));

      x1 = (Lc*(sin(th)*sin(ph)/cos(ps)))*sin(ps)+Lc*(sin(th)*cos(ph));
      y1 = Lc*(sin(th)*sin(ph));
      z1 = zmax;
      x2 = La+(Lc*(sin(th)*sin(ph)/cos(ps)))*sin(ps)+Lc*(sin(th)*cos(ph));
      y2 = Lc*(sin(th)*sin(ph));
      z2 = zmax;
      x3 = xmax;
      y3 = ymax;
      z3 = zmax;

      a6 = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
      b6 = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
      c6 = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
      d6 = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 - x3*z2) - z1*(x2*y3 - x3*y2);

      /* Check Surface 6 Coefficients*/

      if (fabs(a6*x1 + b6*y1 + c6*z1 + d6) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (1) surface (1): %E",
            fabs(a6*x1 + b6*y1 + c6*z1 + d6));

      if (fabs(a6*x2 + b6*y2 + c6*z2 + d6) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (2) surface (1): %E",
          fabs(a6*x2 + b6*y2 + c6*z2 + d6));

      if (fabs(a6*x3 + b6*y3 + c6*z3 + d6) > 1E-6)
        Die(FUNCTION_NAME, "Error in coefficients (3) surface (1): %E",
          fabs(a6*x3 + b6*y3 + c6*z3 + d6));
    }
  else
    Error(0, "Invalid type %ld", type);

  /* Loop over particle types */

  Nmax = 0;
  N0 = 0;

  V1 = 0;
  V2 = 0;

  rpmax = -1.0;

  do
    {
      /* Get number of particles, radius and universe */

      fprintf(outfp,
              "Enter number of particles (> 1) or packing fraction (< 1): ");

      if (fscanf(fp, "%lf", &pf) == EOF)
        Error(0, "Missing number of particles or packing fraction");

      fprintf(outfp, "Enter particle radius (cm): ");

      if (fscanf(fp, "%lf", &rp) == EOF)
        Error(0, "Missing particle radius");

      fprintf(outfp, "Enter particle universe: ");

      if (fscanf(fp, "%ld", &u) == EOF)
        Error(0, "Missing universe");

      /* Compare radius to maximum */

      if (rp > rpmax)
        rpmax = rp;

      /* Calculate number of particles from packing fraction */

      if (pf > 1.0)
        N = (long)pf;
      else
        N = (long)round(V0*pf/((4.0/3.0)*PI*rp*rp*rp));

      /* Calculate volume */

      V1 = V1 + (double)N*(4.0/3.0)*PI*rp*rp*rp;

      /* Check packing fraction */

      if (V1/V0 > 0.7)
        Error(0, "Total packing fraction %1.2f too high", V1/V0);

      /* Allocate memory for data */

      Nmax = Nmax + N;
      p = (struct particle *)Mem(MEM_REALLOC, p, Nmax*sizeof(struct particle));

      /* Check */

      if (p == NULL)
        Error(0, "Memory allocation failed, too many particles");

      /* Create particles */

      for (n = N0; n < Nmax; n++)
        {
          p[n].crp = 0.0;
          p[n].rp = rp;
          p[n].u = u;
          p[n].fixed = NO;
        }

      N0 = Nmax;

      /* Prompt more particles */

      fprintf(outfp, "\nMore particles? (y/n): ");

      if (fscanf(fp, "%s", answer) == EOF)
        Error(0, "Missing answer");

      fprintf(outfp, "\n");
    }
  while (answer[0] == 'y');

  /* Get file name */

  fprintf(outfp, "Enter file name: ");

  if (fscanf(fp, "%s", fname) == EOF)
    Error(0, "Missing file name");

  fprintf(outfp, "\nUse grow and shake algorithm? (y/n): ");

  if (fscanf(fp, "%s", answer) == EOF)
    Error(0, "Missing answer");
  
  fprintf(outfp, "\n");

  if (answer[0] == 'n')
    {
      /* Reset values */

      growth = -INFTY;
      shake = -INFTY;

      /* Set particle radii */

      for (n = 0; n < Nmax; n++)
        p[n].crp = p[n].rp;  
    }
  else
    {
      /* Get Shake Factor */

      fprintf(outfp, 
              "Enter particle shake factor (fraction of particle radius): ");
      
      if (fscanf(fp, "%lf", &shake) == EOF)
        Error(0, "Missing shake factor");
      
      /* Get Growth Rate */
      
      fprintf(outfp, 
              "Enter particle growth rate (fraction of particle radius): ");

      if (fscanf(fp, "%lf", &growth) == EOF)
        Error(0, "Missing growth rate");
      
      /* Check Growth Rate */
      
      if ((growth <= 0.0001) || (growth > 1.0))
        Error(0, "Growth rate must be between 0.0001 and 1.0");

      fprintf(outfp, "\n");
    }
  
  /* Close input file */
      
  if (fp != stdin)
    fclose(fp);

  /* Close outfput file */

  if (outfp != stdout)
    fclose(outfp);

  /* Outfput to stdout */
  
  outfp = stdout;

  /***************************************************************************/

  /***** Prepare data for distribution ***************************************/

  /* Calculate mesh size */

  nx = (long)((xmax - xmin)/(2.0*rpmax));
  ny = (long)((ymax - ymin)/(2.0*rpmax));
  nz = (long)((zmax - zmin)/(2.0*rpmax));

  /* Check cut-offs */

  if (nx < 5)
    nx = 5;
  else if (nx > 200)
    nx = 200;

  if (ny < 5)
    ny = 5;
  else if (ny > 200)
    ny = 200;

  if (nz < 5)
    nz = 5;
  else if (nz > 200)
    nz = 200;

  /* Allocate data for search mesh */

  f = (struct smsh ***)Mem(MEM_ALLOC, nx, sizeof(struct smsh **));

  for(i = 0; i < nx; i++)
    {
      f[i] = (struct smsh **)Mem(MEM_ALLOC, ny, sizeof(struct smsh *));

      for(j = 0; j < ny; j++)
        f[i][j] = (struct smsh *)Mem(MEM_ALLOC, nz, sizeof(struct smsh));
    }

  /* Reset data */

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
        {
          f[i][j][k].sz = 0;
          f[i][j][k].n = NULL;
        }

  /* Set pointers for linked list */

  first = &p[0];

  for (n = 0; n < Nmax - 1; n++)
    p[n].next = &p[n + 1];

  p[n].next = NULL;

  /* Set indexes */

  for (n = 0; n < Nmax; n++)
    p[n].idx = n;

  /***************************************************************************/

  /***** Generate initial distribution ***************************************/

  /* Reset printing counter */

  /* Set id */

  id = 0;

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  /* Init random number sequence */

  seed = ReInitRNG(1);
  SEED[0] = seed;

  /* Randomize */

  fprintf(outfp, "Randomizing %ld particles for initial sampling...\n\n", 
          Nmax);
  
  Nprev = -1;
  do
    { 
      /* Pointer to first unset particle */

      p1 = first;

      /* Randomize */

      if (type == 1)
        {
          while (p1 != NULL)
            {
              rp = p1->crp;

              do
                {
                  x = (2*RandF(id) - 1)*(r - rp);
                  y = (2*RandF(id) - 1)*(r - rp);
                  z = (2*RandF(id) - 1)*(r - rp);
                }
              while (x*x + y*y + z*z > (r - rp)*(r - rp));

              p1->x = x;
              p1->y = y;
              p1->z = z;

              p1 = p1->next;
            }
        }      
      else if (type == 2)
        {
          while (p1 != NULL)
            {
              rp = p1->crp;

              do
                {
                  x = (2*RandF(id) - 1)*(r - rp);
                  y = (2*RandF(id) - 1)*(r - rp);
                  z = RandF(id)*(zmax - zmin) + zmin;
                }
              while ((x*x + y*y  > (r - rp)*(r - rp)) ||
                     (z < zmin + rp) || (z > zmax - rp));

              p1->x = x;
              p1->y = y;
              p1->z = z;

              p1 = p1->next;
            }
        }
      else if (type == 4)
        {
          while (p1 != NULL)
            {
              rp = p1->crp;

              do
                {
                  x = (2*RandF(id) - 1)*(r - rp);
                  y = (2*RandF(id) - 1)*(r - rp);
                  z = RandF(id)*(zmax - zmin) + zmin;
                }
              while ((x*x + y*y  > (r - rp)*(r - rp)) ||
                     (x*x + y*y  < (rmin + rp)*(rmin + rp)) ||
                     (z < zmin + rp) || (z > zmax - rp));

              p1->x = x;
              p1->y = y;
              p1->z = z;

              p1 = p1->next;
            }
        }
      else if (type == 5)
        {
          while (p1 != NULL)
            {
              rp = p1->crp;

              x = (2*RandF(id) - 1)*(xmax - rp);
              y = (2*RandF(id) - 1)*(ymax - rp);
              z = (2*RandF(id) - 1)*(zmax - rp);

              p1->x = x;
              p1->y = y;
              p1->z = z;

              p1 = p1->next;
            }
        }
      else if (type == 6)
        {
          while (p1 != NULL)
            {
              rp = p1->crp;

              do
                {
                  x = (RandF(id))*(xmax-rp);
                  y = (RandF(id))*(ymax-rp);
                  z = (RandF(id))*(zmax-rp);
                }
              while ((rp > -(a1*x + b1*y + c1*z + d1)/
                      sqrt(pow(a1,2) + pow(b1,2) + pow(c1,2))) ||
                     (rp >  (a2*x + b2*y + c2*z + d2)/
                      sqrt(pow(a2,2) + pow(b2,2) + pow(c2,2))) ||
                     (rp > -(a3*x + b3*y + c3*z + d3)/
                      sqrt(pow(a3,2) + pow(b3,2) + pow(c3,2))) ||
                     (rp >  (a4*x + b4*y + c4*z + d4)/
                      sqrt(pow(a4,2) + pow(b4,2) + pow(c4,2))) ||
                     (rp >  (a5*x + b5*y + c5*z + d5)/
                      sqrt(pow(a5,2) + pow(b5,2) + pow(c5,2))) ||
                     (rp > -(a6*x + b6*y + c6*z + d6)/
                      sqrt(pow(a6,2) + pow(b6,2) + pow(c6,2))));
              
              p1->x = x;
              p1->y = y;
              p1->z = z;

              p1 = p1->next;
            }
        }
      else
        {
          while (p1 != NULL)
            {
              rp = p1->crp;

              x = (2*RandF(id) - 1)*(r - rp);
              y = (2*RandF(id) - 1)*(r - rp);
              z = (2*RandF(id) - 1)*(r - rp);

              p1->x = x;
              p1->y = y;
              p1->z = z;

              p1 = p1->next;
            }
        }

      /* Pointer to first unfixed particle */

      p1 = first;

      /* Loop over distribution */

      N = 0;

     while(p1 != NULL)
        {
          /* Get mesh indexes */

          imin = (long)((double)nx*(p1->x - p1->crp - xmin)/(xmax - xmin));
          imax = (long)((double)nx*(p1->x + p1->crp - xmin)/(xmax - xmin));
          jmin = (long)((double)ny*(p1->y - p1->crp - ymin)/(ymax - ymin));
          jmax = (long)((double)ny*(p1->y + p1->crp - ymin)/(ymax - ymin));
          kmin = (long)((double)nz*(p1->z - p1->crp - zmin)/(zmax - zmin));
          kmax = (long)((double)nz*(p1->z + p1->crp - zmin)/(zmax - zmin));

          /* Test boundaries */

          if (imin < 0)
            imin = 0;
          if (imax > nx - 1)
            imax = nx - 1;

          if (jmin < 0)
            jmin = 0;
          if (jmax > ny - 1)
            jmax = ny - 1;

          if (kmin < 0)
            kmin = 0;
          if (kmax > nz - 1)
            kmax = nz - 1;

          /* Set fixed-flag */

          p1->fixed = YES;

          /* Loop particles in mesh to check overlap */

          for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
              for (k = kmin; k <= kmax; k++)
                for (ni = 0; ni < f[i][j][k].sz; ni++)
                  {
                    /* Get index */

                    m = f[i][j][k].n[ni];

                    /* Check that particle is fixed */

                    if (p[m].fixed == NO)
                      Die(FUNCTION_NAME, "Particle not fixed");

                    /* Test */

                    dx = p1->x - p[m].x;
                    dy = p1->y - p[m].y;
                    dz = p1->z - p[m].z;
                    dr = p1->crp + p[m].crp;

                    if (dx*dx + dy*dy + dz*dz <= dr*dr)
                      {
                        /* Overlap, reset flag */

                        p1->fixed = NO;

                        /* Break loop */

                        i = imax + 1;
                        j = jmax + 1;
                        k = kmax + 1;

                        break;
                      }
                  }

          /* Check if flag is set */

          if (p1->fixed == YES)
            {
              /* Add particle in mesh */

              for (i = imin; i <= imax; i++)
                for (j = jmin; j <= jmax; j++)
                  for (k = kmin; k <= kmax; k++)
                    {
                      /* Get number of existing */

                      ni = f[i][j][k].sz;

                      /* Allocate memory */

                      f[i][j][k].n = (long *)Mem(MEM_REALLOC, f[i][j][k].n,
                                                 (ni + 1)*sizeof(long));

                      /* Add index to list */

                      f[i][j][k].n[ni] = p1->idx;

                      /* Add size */

                      f[i][j][k].sz++;
                    }
              
              /* Add volume */
              
              V2 = V2 + (4.0/3.0)*PI*p1->rp*p1->rp*p1->rp;
            }
          else
            {
              /* Add counter */

              N++;
            }

          /* Next particle */

          p1 = p1->next;
        }

      /* Check if status has changed */

      if ((N > 0) && (N != Nprev))
        {
          /* Regenerate list */

          first = NULL;
          p1 = NULL;

          for (n = 0; n < Nmax; n++)
            if (p[n].fixed == NO)
              {
                if (first == NULL)
                  first = &p[n];
                else
                  p1->next = &p[n];

                p1 = &p[n];
                p1->next = NULL;
              }

          /* Print number of particles left */

          fprintf(outfp, 
                  "Overlapping particles: %6ld / %ld pf = %1.5f / %1.5f\n", 
                  N, Nmax, V2/V0, V1/V0);
          Nprev = N;
        }
    }
  while (N > 0);

  /***************************************************************************/

  /***** Grow and Shake Particles ********************************************/

  /* Check if grow and shake algorithm is in use */
  
  if (growth > 0.0)
    {
      /* Reset Fixed Flag (Now referring to Radius rather than position) */

      for(n = 0; n < Nmax; n++)  
        p[n].fixed = NO;

      /* Loop over distribution */      

      rpfracprev = 1.0;
      do
        {
          for(n = 0; n < Nmax; n++)
            {
              /* Attempt to increase particle size */                    

              if (p[n].fixed == NO)
                {
                  /* Increase test Radius size */

                  if (p[n].crp + growth*p[n].rp < p[n].rp) 
                    crp = p[n].crp + growth*p[n].rp;
                  else 
                    crp = p[n].rp;
                  
                  x = p[n].x;
                  y = p[n].y;
                  z = p[n].z;

                  /* Initialize flag */

                  flag = 1;

                  /* Check Boundaries */

                  if (type == 1)
                    {
                      if (x*x + y*y + z*z > (r - crp)*(r - crp)) 
                        flag = 0;
                    }
                  else if (type == 2)
                    {
                      if ((x*x + y*y  > (r - crp)*(r - crp)) ||
                          (z < zmin + crp) || (z > zmax - crp)) 
                        flag = 0;
                    }
                  else if (type == 4)
                    {
                      if ((x*x + y*y  > (r - crp)*(r - crp)) ||
                          (x*x + y*y  < (rmin + crp)*(rmin + crp)) ||
                          (z < zmin + crp) || (z > zmax - crp)) 
                        flag = 0;
                    }
                  else if (type == 5)
                    {
                      if ((x + crp > xmax) || (x - crp < xmin) ||
                          (y + crp > ymax) || (y - crp < ymin) ||
                          (z + crp > zmax) || (z - crp < zmin)) 
                        flag = 0;
                    }
                  else if (type == 6)
                    {
                      if ((crp > -(a1*x + b1*y + c1*z + d1)/
                           sqrt(pow(a1,2) + pow(b1,2) + pow(c1,2))) ||
                          (crp >  (a2*x + b2*y + c2*z + d2)/
                           sqrt(pow(a2,2) + pow(b2,2) + pow(c2,2))) ||
                          (crp > -(a3*x + b3*y + c3*z + d3)/
                           sqrt(pow(a3,2) + pow(b3,2) + pow(c3,2))) ||
                          (crp >  (a4*x + b4*y + c4*z + d4)/
                           sqrt(pow(a4,2) + pow(b4,2) + pow(c4,2))) ||
                          (crp >  (a5*x + b5*y + c5*z + d5)/
                           sqrt(pow(a5,2) + pow(b5,2) + pow(c5,2))) ||
                          (crp > -(a6*x + b6*y + c6*z + d6)/
                           sqrt(pow(a6,2) + pow(b6,2) + pow(c6,2)))) 
                        flag = 0;
                    }
                  else
                    {
                      if ((x + crp > xmax) || (x - crp < xmin) ||
                          (y + crp > ymax) || (y - crp < ymin) ||
                          (z + crp > zmax) || (z - crp < zmin)) 
                        flag = 0;
                    }

                  
                  /* Check Overlap */
                  
                  /* Get mesh indexes */
                  
                  imin = (long)((double)nx*(x - crp - xmin)/(xmax - xmin));
                  imax = (long)((double)nx*(x + crp - xmin)/(xmax - xmin));
                  jmin = (long)((double)ny*(y - crp - ymin)/(ymax - ymin));
                  jmax = (long)((double)ny*(y + crp - ymin)/(ymax - ymin));
                  kmin = (long)((double)nz*(z - crp - zmin)/(zmax - zmin));
                  kmax = (long)((double)nz*(z + crp - zmin)/(zmax - zmin));
                  
                  /* Test boundaries */
                  
                  if (imin < 0)
                    imin = 0;
                  if (imax > nx - 1)
                    imax = nx - 1;
                  
                  if (jmin < 0)
                    jmin = 0;
                  if (jmax > ny - 1)
                    jmax = ny - 1;
                  
                  if (kmin < 0)
                    kmin = 0;
                  if (kmax > nz - 1)
                    kmax = nz - 1;
                  
                  /* Loop particles in mesh to check overlap */
                  
                  for (i = imin; i <= imax; i++)
                    for (j = jmin; j <= jmax; j++)
                      for (k = kmin; k <= kmax; k++)
                        for (ni = 0; ni < f[i][j][k].sz; ni++)
                          {
                            /* Get index */

                            m = f[i][j][k].n[ni];

                            if (p[n].idx != m)
                              {
                                /* Test */
                                
                                dx = x - p[m].x;
                                dy = y - p[m].y;
                                dz = z - p[m].z;
                                dr = crp + p[m].crp;
                                
                                if (dx*dx + dy*dy + dz*dz < dr*dr)
                                  {
                                    /* Overlap, reset flag */
                                    
                                    flag = 0;
                                    
                                    /* Break loop */
                                    
                                    i = imax + 1;
                                    j = jmax + 1;
                                    k = kmax + 1;
                                    
                                    break;
                                  }
                              }
                          }
                  
                  /* Increase particle size if appropriate  */
                  
                  if (flag == 1)
                    {
                      /* Remove particle from old mesh *************/

                      /* Get mesh indexes */
                      
                      imin = (long)((double)nx*(p[n].x - p[n].crp - xmin)/
                                    (xmax - xmin));
                      imax = (long)((double)nx*(p[n].x + p[n].crp - xmin)/
                                    (xmax - xmin));
                      jmin = (long)((double)ny*(p[n].y - p[n].crp - ymin)/
                                    (ymax - ymin));
                      jmax = (long)((double)ny*(p[n].y + p[n].crp - ymin)/
                                    (ymax - ymin));
                      kmin = (long)((double)nz*(p[n].z - p[n].crp - zmin)/
                                    (zmax - zmin));
                      kmax = (long)((double)nz*(p[n].z + p[n].crp - zmin)/
                                    (zmax - zmin));
                      
                      /* Test boundaries */
                      
                      if (imin < 0)
                        imin = 0;
                      if (imax > nx - 1)
                        imax = nx - 1;
                      
                      if (jmin < 0)
                        jmin = 0;
                      if (jmax > ny - 1)
                        jmax = ny - 1;
                      
                      if (kmin < 0)
                        kmin = 0;
                      if (kmax > nz - 1)
                        kmax = nz - 1;
                      
                      /* Remove particle from old mesh */
                      
                      for (i = imin; i <= imax; i++)
                        for (j = jmin; j <= jmax; j++)
                          for (k = kmin; k <= kmax; k++)
                            {
                              if (f[i][j][k].sz > 1)
                                {
                                  /* Initialize Counter */

                                  count = 0;

                                  /* Initialize temp mesh  */
                                  temp = (long *)Mem(MEM_REALLOC, temp,
                                                     (f[i][j][k].sz - 1)*
                                                     sizeof(long));

                                  /* Store old mesh particles that are not */
                                  /* current in temp */

                                  for (ni = 0; ni < f[i][j][k].sz; ni++)
                                    {
                                      if (p[n].idx != f[i][j][k].n[ni])
                                        {
                                          temp[count] = f[i][j][k].n[ni];
                                          count++;
                                        }
                                    }
                                }
                              
                              /* Resize Mesh */
                              
                              ni = f[i][j][k].sz;
                              
                              if (f[i][j][k].sz == 0)
                                f[i][j][k].n = NULL;

                              else if (f[i][j][k].sz == 1)
                                {
                                  Mem(MEM_FREE, f[i][j][k].n);
                                  f[i][j][k].n = NULL;
                                  f[i][j][k].sz--;
                                }
                              else
                                {
                                  f[i][j][k].n = (long *)Mem(MEM_REALLOC, 
                                                             f[i][j][k].n,
                                                             (ni-1)*
                                                             sizeof(long));
                                  f[i][j][k].sz--;

                                  for (ni = 0; ni < f[i][j][k].sz; ni++)
                                    f[i][j][k].n[ni] = temp[ni];
                                }
                            }

                      /* Add particle to new mesh *************/
                      /* Set new particle location  */

                      p[n].x = x;
                      p[n].y = y;
                      p[n].z = z;
                      
                      /* Get mesh indexes */

                      imin = (long)((double)nx*(p[n].x - crp - xmin)/
                                    (xmax - xmin));
                      imax = (long)((double)nx*(p[n].x + crp - xmin)/
                                (xmax - xmin));
                      jmin = (long)((double)ny*(p[n].y - crp - ymin)/
                                    (ymax - ymin));
                      jmax = (long)((double)ny*(p[n].y + crp - ymin)/
                                    (ymax - ymin));
                      kmin = (long)((double)nz*(p[n].z - crp - zmin)/
                                    (zmax - zmin));
                      kmax = (long)((double)nz*(p[n].z + crp - zmin)/
                                    (zmax - zmin));
                      
                      /* Test boundaries */
                      
                      if (imin < 0)
                        imin = 0;
                      if (imax > nx - 1)
                        imax = nx - 1;
                      
                      if (jmin < 0)
                        jmin = 0;
                      if (jmax > ny - 1)
                        jmax = ny - 1;
                      
                      if (kmin < 0)
                        kmin = 0;
                      if (kmax > nz - 1)
                        kmax = nz - 1;
                      
                      /* Add particle in mesh */
                      
                      for (i = imin; i <= imax; i++)
                        for (j = jmin; j <= jmax; j++)
                          for (k = kmin; k <= kmax; k++)
                            {
                              /* Get number of existing */
                              
                              ni = f[i][j][k].sz;
                              
                              /* Allocate memory */
                              
                              f[i][j][k].n = (long *)Mem(MEM_REALLOC, 
                                                         f[i][j][k].n,
                                                         (ni + 1)*sizeof(long));

                              /* Add index to list */
                              
                              f[i][j][k].n[ni] = p[n].idx;
                              
                              /* Add size */
                              
                              f[i][j][k].sz++;
                            }
                      
                      
                      /* Increase Current Radius */
                      
                      p[n].crp = crp;
                      
                      /* Set to full radius and set fixed status if needed */
                      
                      if (p[n].crp == p[n].rp) 
                        p[n].fixed = YES;
                    }
                }
              
              /* Attempt to move particle to new location */

              if (type == 1)
                {
                  rp = p[n].rp;
                  crp = p[n].crp;
                  
                  do
                    {
                      x = (2*RandF(id) - 1)*(shake*rp) + p[n].x;
                      y = (2*RandF(id) - 1)*(shake*rp) + p[n].y;
                      z = (2*RandF(id) - 1)*(shake*rp) + p[n].z;
                    }
                  while (x*x + y*y + z*z > (r - crp)*(r - crp));
                }
              else if (type == 2)
                {
                  rp = p[n].rp;
                  crp = p[n].crp;             
                  
                  do
                    {
                      x = (2*RandF(id) - 1)*(shake*rp) + p[n].x;
                      y = (2*RandF(id) - 1)*(shake*rp) + p[n].y;
                      z = (2*RandF(id) - 1)*(shake*rp) + p[n].z;
                    }
                  while ((x*x + y*y  > (r - crp)*(r - crp)) ||
                         (z < zmin + crp) || (z > zmax - crp));
                }
              else if (type == 4)
                {
                  rp = p[n].rp;
                  crp = p[n].crp;
                  
                  do
                    {
                      x = (2*RandF(id) - 1)*(shake*rp) + p[n].x;
                      y = (2*RandF(id) - 1)*(shake*rp) + p[n].y;
                      z = (2*RandF(id) - 1)*(shake*rp) + p[n].z;
                    }
                  while ((x*x + y*y  > (r - crp)*(r - crp)) ||
                         (x*x + y*y  < (rmin + crp)*(rmin + crp)) ||
                         (z < zmin + crp) || (z > zmax - crp));
                }
              else if (type == 5)
                {
                  rp = p[n].rp;
                  crp = p[n].crp;
                  
                  do
                    {
                      x = (2*RandF(id) - 1)*(shake*rp) + p[n].x;
                      y = (2*RandF(id) - 1)*(shake*rp) + p[n].y;
                      z = (2*RandF(id) - 1)*(shake*rp) + p[n].z;
                    }
                  while ((x + crp > xmax) || (x - crp < xmin) ||
                         (y + crp > ymax) || (y - crp < ymin) ||
                         (z + crp > zmax) || (z - crp < zmin));
                }
              else if (type == 6)
                {
                  rp = p[n].rp;
                  crp = p[n].crp;
                  
                  do
                    {
                      x = (2*RandF(id) - 1)*(shake*rp) + p[n].x;
                      y = (2*RandF(id) - 1)*(shake*rp) + p[n].y;
                      z = (2*RandF(id) - 1)*(shake*rp) + p[n].z;
                    }
                  while ((crp > -(a1*x + b1*y + c1*z + d1)/
                          sqrt(pow(a1,2) + pow(b1,2) + pow(c1,2)))||
                         (crp >  (a2*x + b2*y + c2*z + d2)/
                          sqrt(pow(a2,2) + pow(b2,2) + pow(c2,2)))||
                         (crp > -(a3*x + b3*y + c3*z + d3)/
                          sqrt(pow(a3,2) + pow(b3,2) + pow(c3,2)))||
                         (crp >  (a4*x + b4*y + c4*z + d4)/
                          sqrt(pow(a4,2) + pow(b4,2) + pow(c4,2)))||
                         (crp >  (a5*x + b5*y + c5*z + d5)/
                          sqrt(pow(a5,2) + pow(b5,2) + pow(c5,2)))||
                         (crp > -(a6*x + b6*y + c6*z + d6)/
                          sqrt(pow(a6,2) + pow(b6,2) + pow(c6,2))) );
                }
              else
                {
                  rp = p[n].rp;
                  crp = p[n].crp;
                  
                  do
                    {
                      x = (2*RandF(id) - 1)*(shake*rp) + p[n].x;
                      y = (2*RandF(id) - 1)*(shake*rp) + p[n].y;
                      z = (2*RandF(id) - 1)*(shake*rp) + p[n].z;
                    }
                  while ((x + crp > xmax) || (x - crp < xmin) ||
                         (y + crp > ymax) || (y - crp < ymin) ||
                         (z + crp > zmax) || (z - crp < zmin));
                }
              
              /* Get mesh indexes */
              
              imin = (long)((double)nx*(x - crp - xmin)/(xmax - xmin));
              imax = (long)((double)nx*(x + crp - xmin)/(xmax - xmin));
              jmin = (long)((double)ny*(y - crp - ymin)/(ymax - ymin));
              jmax = (long)((double)ny*(y + crp - ymin)/(ymax - ymin));
              kmin = (long)((double)nz*(z - crp - zmin)/(zmax - zmin));
              kmax = (long)((double)nz*(z + crp - zmin)/(zmax - zmin));
              
              /* Test boundaries */
              
              if (imin < 0)
                imin = 0;
              if (imax > nx - 1)
                imax = nx - 1;
              
              if (jmin < 0)
                jmin = 0;
              if (jmax > ny - 1)
                jmax = ny - 1;
              
              if (kmin < 0)
                kmin = 0;
              if (kmax > nz - 1)
                kmax = nz - 1;
              
              /* Set match-flag */
              
              flag = 1;
              
              /* Loop particles in mesh to check overlap */
              
              for (i = imin; i <= imax; i++)
                for (j = jmin; j <= jmax; j++)
                  for (k = kmin; k <= kmax; k++)
                    for (ni = 0; ni < f[i][j][k].sz; ni++)
                      {
                        /* Get index */
                        
                        m = f[i][j][k].n[ni];
                        
                        if (p[n].idx != m)
                          {
                            
                            /* Test */
                            
                            dx = x - p[m].x;
                            dy = y - p[m].y;
                            dz = z - p[m].z;
                            dr = p[n].crp + p[m].crp;
                            
                            if (dx*dx + dy*dy + dz*dz < dr*dr)
                              {
                                /* Overlap, reset flag */
                                
                                flag = 0;

                                /* Break loop */
                                
                                i = imax + 1;
                                j = jmax + 1;
                                k = kmax + 1;
                                
                                break;
                              }
                          }
                      }
              
              /* Check if flag is set */
              
              if (flag == 1)
                {
                  
                  /* Remove particle from old mesh */

                  /* Get mesh indexes */

                  imin = (long)((double)nx*(p[n].x - p[n].crp - xmin)/
                                (xmax - xmin));
                  imax = (long)((double)nx*(p[n].x + p[n].crp - xmin)/
                                (xmax - xmin));
                  jmin = (long)((double)ny*(p[n].y - p[n].crp - ymin)/
                                (ymax - ymin));
                  jmax = (long)((double)ny*(p[n].y + p[n].crp - ymin)/
                                (ymax - ymin));
                  kmin = (long)((double)nz*(p[n].z - p[n].crp - zmin)/
                                (zmax - zmin));
                  kmax = (long)((double)nz*(p[n].z + p[n].crp - zmin)/
                                (zmax - zmin));
                  
                  /* Test boundaries */
                  
                  if (imin < 0)
                    imin = 0;
                  if (imax > nx - 1)
                    imax = nx - 1;
                  
                  if (jmin < 0)
                    jmin = 0;
                  if (jmax > ny - 1)
                    jmax = ny - 1;
                  
                  if (kmin < 0)
                    kmin = 0;
                  if (kmax > nz - 1)
                    kmax = nz - 1;
                  
                  /* Remove particle from old mesh */
                  
                  for (i = imin; i <= imax; i++)
                    for (j = jmin; j <= jmax; j++)
                      for (k = kmin; k <= kmax; k++)
                        {
                          if (f[i][j][k].sz > 1)
                            {
                              /* Initialize Counter */
                              
                              count = 0;
                              
                              /* Initialize temp mesh  */
                              temp = (long *)Mem(MEM_REALLOC, temp,
                                                 (f[i][j][k].sz - 1)*
                                                 sizeof(long));
                              
                              /* Store old mesh particles that are not */
                              /* current in temp */

                              for (ni = 0; ni < f[i][j][k].sz; ni++)
                                {
                                  if(p[n].idx != f[i][j][k].n[ni])
                                    {
                                      temp[count] = f[i][j][k].n[ni];
                                      count++;
                                    }
                                }
                            }
                          
                          /* Resize Mesh */
                          
                          if (f[i][j][k].sz == 0)
                            f[i][j][k].n = NULL;
                          
                          else if (f[i][j][k].sz == 1)
                            {
                              Mem(MEM_FREE, f[i][j][k].n);
                              f[i][j][k].n = NULL;
                              f[i][j][k].sz--;
                            }
                          else
                            {
                              f[i][j][k].n = (long *)Mem(MEM_REALLOC, 
                                                         f[i][j][k].n,
                                                         (f[i][j][k].sz - 1)*
                                                         sizeof(long));
                              f[i][j][k].sz--;
                              
                              for (ni = 0; ni < f[i][j][k].sz; ni++)
                                f[i][j][k].n[ni] = temp[ni];
                            }
                        }
                  
                  /* Add particle to new mesh */
                  /* Set new particle location  */

                  p[n].x = x;
                  p[n].y = y;
                  p[n].z = z;
                  
                  /* Get mesh indexes */
                  
                  imin = (long)((double)nx*(p[n].x - p[n].crp - xmin)/
                                (xmax - xmin));
                  imax = (long)((double)nx*(p[n].x + p[n].crp - xmin)/
                                (xmax - xmin));
                  jmin = (long)((double)ny*(p[n].y - p[n].crp - ymin)/
                                (ymax - ymin));
                  jmax = (long)((double)ny*(p[n].y + p[n].crp - ymin)/
                                (ymax - ymin));
                  kmin = (long)((double)nz*(p[n].z - p[n].crp - zmin)/
                                (zmax - zmin));
                  kmax = (long)((double)nz*(p[n].z + p[n].crp - zmin)/
                                (zmax - zmin));
                  
                  /* Test boundaries */
                  
                  if (imin < 0)
                    imin = 0;
                  if (imax > nx - 1)
                    imax = nx - 1;
                  
                  if (jmin < 0)
                    jmin = 0;
                  if (jmax > ny - 1)
                    jmax = ny - 1;
                  
                  if (kmin < 0)
                    kmin = 0;
                  if (kmax > nz - 1)
                    kmax = nz - 1;
                  
                  /* Add particle in mesh */
                  
                  for (i = imin; i <= imax; i++)
                    for (j = jmin; j <= jmax; j++)
                      for (k = kmin; k <= kmax; k++)
                        {
                          /* Get number of existing */
                          
                          ni = f[i][j][k].sz;
                          
                          /* Allocate memory */
                          
                          f[i][j][k].n = (long *)Mem(MEM_REALLOC, f[i][j][k].n,
                                                     (ni + 1)*sizeof(long));
                          
                          /* Add index to list */
                          
                          f[i][j][k].n[ni] = p[n].idx;
                          
                          /* Add size */
                          
                          f[i][j][k].sz++;
                        }
                }
            }
          
          /* Print Smallest Radius */
          
          small = 0;
          rpfracmin = 1.0;

          for(n = 0; n < Nmax; n++)
            if (p[n].crp / p[n].rp < rpfracmin)
              {
                rpfracmin = p[n].crp / p[n].rp;
                small = n;
              }
                
          if (small < 0)
            Die(FUNCTION_NAME, "Error");
     
          if (rpfracprev != rpfracmin)
            {
              fprintf(outfp, "Smallest particle radius fraction: %1.4f\n",
                      rpfracmin);
              rpfracprev = rpfracmin;
            }

        }
      while (rpfracmin < 1.0);
    }

  /* Verify */

#ifdef DEBUG

  fprintf(outfp, "\n");

  for (n = 0; n < Nmax; n++)
    {
      fprintf(outfp, "Verifying distribution: %6ld / %ld \n", n + 1, Nmax);

      for (m = n + 1; m < Nmax; m++)
        {
          /* Check overlap */

          dx = p[n].x - p[m].x;
          dy = p[n].y - p[m].y;
          dz = p[n].z - p[m].z;
          dr = p[n].rp + p[m].rp;

          if (dx*dx + dy*dy + dz*dz < dr*dr)
            Die(FUNCTION_NAME, "Particles overlapping");
        }
    }

#endif

  /* Finished. Write data to file */

  fprintf(outfp, "\nWriting final distribution to file \"%s\"...\n\n", fname);

  /* Write file */

  fp = fopen(fname, "w");

  for (n = 0; n < Nmax; n++)
    fprintf(fp, "%12.5E %12.5E %12.5E %11.5E %ld\n", p[n].x, p[n].y, p[n].z,
            p[n].rp, p[n].u);

  fclose(fp);

  fprintf(outfp, "%ld particles, packing fraction = %1.5f\n\n", Nmax, V1/V0);

  /* Free memory */

  Mem(MEM_FREE, p);

  if (temp != NULL)
    Mem(MEM_FREE, temp);

  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
        {
          for (k = 0; k < nz; k++)
            if ((f[i][j][k].n) != NULL)
              Mem(MEM_FREE, f[i][j][k].n);

          Mem(MEM_FREE, f[i][j]);
        }
      Mem(MEM_FREE, f[i]);
    }

  Mem(MEM_FREE, f);
}

/*****************************************************************************/
