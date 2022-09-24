/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : commandlineplotter.c                           */
/*                                                                           */
/* Created:       2020/06/17 (JLe)                                           */
/* Last modified: 2020/06/17 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Command-line plotter routine                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CommandLinePlotter:"

/*****************************************************************************/

void CommandLinePlotter()
{
#ifndef NO_GFX_MODE

  long gpl, xp, yp, bou, par, ax, qp, ptr, scale, *mtx, stop, n, m, l, dim;
  long R[MAX_PLOT_COLORS], G[MAX_PLOT_COLORS], B[MAX_PLOT_COLORS];
  long refresh, mat, mat0, r, g, b, wwd;
  double xmin, xmax, ymin, ymax, zmin, zmax, tmp, pos, E, fmin, fmax, f;
  double x, y, z, xx, yy, zz, dx, dy, dz, min1, max1, min2, max2;
  char line[MAX_STR], word[MAX_STR], fname[MAX_STR], mname[MAX_STR];
  gdImagePtr im;
  long palette[MAX_PLOT_COLORS];
  FILE *fp;

  /* Check plotter mode */

  if ((long)RDB[DATA_INTERACTIVE_PLOT_MODE] == NO)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    Error(0, "Command-line plotter does not support multiple MPI tasks");

  /* Print */

  printf("Command line plotter started...\n\n");

  /***************************************************************************/

  /***** Initialize data *****************************************************/

  /* Get pointer to plot or create default */

  if ((gpl = (long)RDB[DATA_PTR_GPL0]) < VALID_PTR)
    {
      /* Create a default plot structure */

      gpl = NewItem(DATA_PTR_GPL0, GPL_BLOCK_SIZE);

      /* Init values */

      WDB[gpl + GPL_XMIN] = -INFTY;
      WDB[gpl + GPL_XMAX] =  INFTY;
      WDB[gpl + GPL_YMIN] = -INFTY;
      WDB[gpl + GPL_YMAX] =  INFTY;
      WDB[gpl + GPL_ZMIN] = -INFTY;
      WDB[gpl + GPL_ZMAX] =  INFTY;
      WDB[gpl + GPL_POS] =  -INFTY;
      WDB[gpl + GPL_PLOT_BOUND] = 2.0;
      WDB[gpl + GPL_GEOM] = 3.0;
      WDB[gpl + GPL_PIX_X] = 500.0;
      WDB[gpl + GPL_PIX_Y] = 500.0;
    }

  /* Set plotter flag on */

  WDB[DATA_PLOTTER_MODE] = (double)YES;

  /* Put RGB values for plot */

  PutPlotColors(gpl, R, G, B);

  /* Get geometry boundaries */

  xmin = RDB[DATA_GEOM_MINX];
  xmax = RDB[DATA_GEOM_MAXX];
  ymin = RDB[DATA_GEOM_MINY];
  ymax = RDB[DATA_GEOM_MAXY];
  zmin = RDB[DATA_GEOM_MINZ];
  zmax = RDB[DATA_GEOM_MAXZ];

  /* Check if boundaries are set by user */

  if (RDB[gpl + GPL_XMIN] != -INFTY)
    xmin = RDB[gpl + GPL_XMIN];

  if (RDB[gpl + GPL_XMAX] !=  INFTY)
    xmax = RDB[gpl + GPL_XMAX];

  if (RDB[gpl + GPL_YMIN] != -INFTY)
    ymin = RDB[gpl + GPL_YMIN];

  if (RDB[gpl + GPL_YMAX] !=  INFTY)
    ymax = RDB[gpl + GPL_YMAX];

  if (RDB[gpl + GPL_ZMIN] != -INFTY)
    zmin = RDB[gpl + GPL_ZMIN];

  if (RDB[gpl + GPL_ZMAX] !=  INFTY)
    zmax = RDB[gpl + GPL_ZMAX];

  /* Check boundaries and swap */

  if (xmin > xmax)
    {
      tmp = xmax;
      xmax = xmin;
      xmin = tmp;
    }

  if (ymin > ymax)
    {
      tmp = ymax;
      ymax = ymin;
      ymin = tmp;
    }

  if (zmin > zmax)
    {
      tmp = zmax;
      zmax = zmin;
      zmin = tmp;
    }

  /* Get dimensions */

  dim = (long)RDB[DATA_GEOM_DIM];

  /* Move limits away from boundaries to avoid numerical errors in */
  /* some systems. */

  xmin = xmin + EXTRAP_L;
  xmax = xmax - EXTRAP_L;
  ymin = ymin + EXTRAP_L;
  ymax = ymax - EXTRAP_L;
  zmin = zmin + EXTRAP_L;
  zmax = zmax - EXTRAP_L;

  /* Set image size */

  xp = (long)RDB[gpl + GPL_PIX_X];
  yp = (long)RDB[gpl + GPL_PIX_Y];

  /* Boundary mode (1 = plot cell boundaries, 2 = plot material */
  /*                boundaries, 3 = plot both) */

  if ((((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0) ||
       ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)) &&
      ((long)RDB[gpl + GPL_IMP_SCALE] == 0))
    bou = 0;
  else
    bou = (long)RDB[gpl + GPL_PLOT_BOUND];

  /* Get weight window particle type */

  if ((ptr = (long)RDB[DATA_PTR_WWD0]) > VALID_PTR)
    {
      /* Check for multiple */

      if (NextItem(ptr) < VALID_PTR)
        par = (long)RDB[ptr + WWD_PARTICLE_TYPE];
      else
        par = 0;
    }
  else
    par = 0;

  /* Axis type */

  ax = (long)RDB[gpl + GPL_GEOM];
  CheckValue(FUNCTION_NAME, "ax", "", ax, 1, 3);

  /* Get position */

  if ((pos = RDB[gpl + GPL_POS]) == -INFTY)
    {
      /* Not defined */

      if (ax == PLOT_MODE_YZ)
        pos = 0.5*(xmin + xmax);
      else if (ax == PLOT_MODE_XZ)
        pos = 0.5*(ymin + ymax);
      else if (ax == PLOT_MODE_XY)
        pos = 0.5*(zmin + zmax);
    }

  /* Check if importances */

  if ((wwd = (long)RDB[DATA_PTR_WWD0]) < VALID_PTR)
    {
      /* Reset data */

      scale = 0;
      E = 0.0;
      fmin = 0.0;
      fmax = 0.0;
    }
  else
    {
      /* Scale for importance plotter */

      scale = (long)RDB[gpl + GPL_IMP_SCALE];
      CheckValue(FUNCTION_NAME, "scale", "", scale, 0, 4);

      /* Get energy */

      E = RDB[gpl + GPL_IMP_E];
      CheckValue(FUNCTION_NAME, "E", "", E, -INFTY, INFTY);

      /* Get maximum and maximum importances */

      fmin = RDB[gpl + GPL_IMP_MIN];
      CheckValue(FUNCTION_NAME, "fmin", "", fmin, 0.0, 1E+200);

      fmax = RDB[gpl + GPL_IMP_MAX];
      CheckValue(FUNCTION_NAME, "fmax", "", fmax, fmin, 1E+200);

      /* Check if values are not set */

      if (E == 0.0)
        E = 1.0;
      if (fmin == 0.0)
        fmin = RDB[wwd + WWD_MESH_MIN];
      if (fmax == 0.0)
        fmax = RDB[wwd + WWD_MESH_MAX];
    }

  /* Get quick plotter option */

  qp = (long)RDB[DATA_QUICK_PLOT_MODE];

  /* Default file name */

  sprintf(fname, "%s_geom1.png", GetText(DATA_PTR_INPUT_FNAME));

  /***************************************************************************/

  /***** Main loop ***********************************************************/

  /* Loop until exit */

  while (1 != 2)
    {
      /* Reset exit flag */

      stop = NO;

      /***********************************************************************/

      /***** Command prompt **************************************************/

      /* Calculate origin */

      x = 0.5*(xmax + xmin);
      y = 0.5*(ymax + ymin);
      z = 0.5*(zmax + zmin);

      /* Calculate widths */

      dx = xmax - xmin;
      dy = ymax - ymin;
      dz = zmax - zmin;

      /* Get position */

      if (ax == PLOT_MODE_YZ)
        x = pos;
      else if (ax == PLOT_MODE_XZ)
        y = pos;
      else if (ax == PLOT_MODE_XY)
        z = pos;

      /* Reset line */

      for (n = 0; n < MAX_STR; n++)
        line[n] = '\0';

      /* Print prompt */

      if (dim == 3)
        printf("[ %1.5f %1.5f %1.5f ] ", x, y, z);
      else
        printf("[ %1.5f %1.5f ] ", x, y);

      if (ax == PLOT_MODE_YZ)
        printf("YZ > ");
      else if (ax == PLOT_MODE_XZ)
        printf("XZ > ");
      else if (ax == PLOT_MODE_XY)
        printf("XY > ");

      /* Read command line */

      n = 0;
      while ((line[n] = fgetc(stdin)) != '\n')
        {
          /* Convert EOF to exit */

          if (line[n] == EOF)
            {
              sprintf(&line[n], "\nexit\n");
              break;
            }
          else
            n++;
        }

      /* Reset refresh flag */

      refresh = NO;

      /* Get first command word */

      l = 0;
      if (sscanf(&line[l], "%s", word) > 0)
        {
          /* Check */

          if (!strcasecmp(word, "exit") || !strcasecmp(word, "xx") ||
              !strcasecmp(word, "end"))
            {
              /***************************************************************/

              /***** Exit plotter ********************************************/

              /* Set flag */

              stop = YES;

              /* Print */

              printf("\nExit command-line plotter\n\n");

              /***************************************************************/
            }

          else if (!strcasecmp(word, "help") || !strcasecmp(word, "h"))
            {
              /***************************************************************/

              /***** Print instructions **************************************/

              /* Print */

              printf("\nList of commands:\n\n");
              printf("exit : exit plotter\n");
              printf("file : set output file name\n");
              printf("qp   : toggle quick plot mode on/off [off by default]\n");
              printf("bou  : boundary mode [0 = none, 1 = cell, 2 = material, 3 = both]\n");
              printf("pix  : set image size in pixels\n");
              printf("orig : set origin\n");
              printf("rgb  : set material color [name and RGB values]\n");
              printf("max  : reset to maximum dimensions\n");
              printf("sq   : set square shape\n");
              printf("box  : set plot area [bounding box coordinates]\n");
              printf("ax   : set plot axis [xz, yz or xy]\n");
              printf("zo   : zoom out [optional factor]\n");
              printf("zi   : zoom in [optional factor]\n");
              printf("w    : move west [optional distance]\n");
              printf("s    : move south [optional distance]\n");
              printf("n    : move north [optional distance]\n");
              printf("e    : move east [optional distance]\n");
              printf("u    : move up [optional distance]\n");
              printf("d    : move down [optional distance]\n");
              printf("imp  : WW mesh [0 = off, 4 = lin cell, 5 = log cell, 6 = lin src, 7 = log src]\n");
              printf("par  : set particle type for WW mesh plot [n = neutron, p = photon]\n");
              printf("ene  : set particle energy for WW mesh plot [in MeV]\n\n");

              /***************************************************************/
            }
          else if (!strcasecmp(word, "file"))
            {
              /***************************************************************/

              /***** Set output file *****************************************/

              /* Read name */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                {
                  /* Set file */

                  strcpy(fname, word);

                  /* Print */

                  printf("\nOutput file name set to: \"%s\"\n\n", fname);
                }
              else
                printf("\nOutput file name: \"%s\"\n\n", fname);

              /***************************************************************/
            }
          else if (!strcasecmp(word, "par"))
            {
              /***************************************************************/

              /***** Set particle type for importance mesh plot **************/

              /* Read option */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                {
                  /* Check */

                  if (!strcasecmp(word, "1") || !strcasecmp(word, "p") ||
                      !strcasecmp(word, "g"))
                    {
                      /* Set particle type */

                      par = PARTICLE_TYPE_GAMMA;

                      /* Print */

                      printf("\nParticle type set to: photon\n\n");

                      /* Set refresh flag */

                      refresh = YES;
                    }
                  else if(!strcasecmp(word, "2") || !strcasecmp(word, "n"))
                    {
                      /* Set particle type */

                      par = PARTICLE_TYPE_NEUTRON;

                      /* Print */

                      printf("\nParticle type set to: neutron\n\n");

                      /* Set refresh flag */

                      refresh = YES;
                    }
                  else
                    printf("\nInvalid particle type: %s\n\n", word);
                }
              else
                {
                  /* Check */

                  if (par == PARTICLE_TYPE_GAMMA)
                    printf("\nParticle type: photon\n\n");
                  else if (par == PARTICLE_TYPE_NEUTRON)
                    printf("\nParticle type: neutron\n\n");
                  else
                    printf("\nParticle type: undefined\n\n");
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "ene"))
            {
              /***************************************************************/

              /***** Set energy for importance mesh plot *********************/

              /* Read option */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                {
                  /* Convert */

                  if ((xx = atof(word)) > 0.0)
                    {
                      /* Put energy */

                      E = xx;

                      /* Print */

                      printf("\nParticle energy set to: %1.5E MeV\n\n", E);

                      /* Refresh */

                      refresh = YES;
                    }
                  else
                    printf("\nInvalid energy: %s\n\n", word);
                }
              else
                printf("\nEnergy: %1.5E MeV\n\n", E);

              /***************************************************************/
            }
          else if (!strcasecmp(word, "qp"))
            {
              /***************************************************************/

              /***** Toggle quick plot mode **********************************/

              if (qp == NO)
                {
                  /* Set */

                  qp = YES;

                  /* Print */

                  printf("\nQuick plot mode on\n\n");
                }
              else
                {
                  /* Set */

                  qp = NO;

                  /* Print */

                  printf("\nQuick plot mode off\n\n");
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "imp"))
            {
              /***************************************************************/

              /***** Importance mesh plotte **********************************/

              if ((long)RDB[DATA_PTR_WWD0] < VALID_PTR)
                printf("\nImportances not available\n\n");
              else
                {
                  /* Reset mode */

                  n = -1000000000;

                  /* Read mode */

                  l = l + strlen(word) + 1;
                  if (sscanf(&line[l], "%s", word) > 0)
                    n = atol(word);

                  /* Read limits */

                  xx = -INFTY;
                  yy = -INFTY;

                  l = l + strlen(word) + 1;
                  if (sscanf(&line[l], "%s", word) > 0)
                    xx = atof(word);

                  l = l + strlen(word) + 1;
                  if (sscanf(&line[l], "%s", word) > 0)
                    yy = atof(word);

                  /* Check */

                  if (n == -1000000000)
                    {
                      if (scale == 0)
                        printf("\nImportances off: [ %1.5E %1.5E ]\n\n",
                               fmin, fmax);
                      else
                        {
                          if (scale == 1)
                            printf("\nCell importances, linear scale: ");
                          else if (scale == 2)
                            printf("\nCell importances, log scale: ");
                          else if (scale == 3)
                            printf("\nSource importances, log scale: ");
                          else if (scale == 4)
                            printf("\nSource importances, log scale: ");
                          else
                            Die(FUNCTION_NAME, "WTF?");

                          printf("[ %1.5E %1.5E ]\n\n", fmin, fmax);
                        }
                    }
                  else if (n == 0)
                    {
                      /* Switch importance plot off */

                      PutPlotColors(NO, R, G, B);
                      scale = 0;

                      /* Print */

                      printf("\nImportances off\n\n");

                      /* Refresh */

                      refresh = YES;
                    }
                  else if ((n < 4) || (n > 7))
                    printf("\nInvalid option\n\n");
                  else if ((xx == -INFTY) && (yy == -INFTY))
                    {
                      /* Switch importance plot on */

                      PutPlotColors(YES, R, G, B);
                      scale = n - 3;

                      /* Print */

                      if (n == 4)
                        printf("\nCell importances, linear scale: ");
                      else if (n == 5)
                        printf("\nCell importances, log scale: ");
                      else if (n == 6)
                        printf("\nSource importances, log scale: ");
                      else if (n == 7)
                        printf("\nSource importances, log scale: ");
                      else
                        Die(FUNCTION_NAME, "WTF?");

                      printf("[ %1.5E %1.5E ]\n\n", fmin, fmax);

                      /* Refresh */

                      refresh = YES;
                    }
                  else if ((xx == -INFTY) || (yy == -INFTY))
                    printf("\nMissing values\n\n");
                  else if ((xx < 0) || (yy < 0) || (yy < xx))
                    printf("\nInvalid values\n\n");
                  else
                    {
                      /* Switch importance plot on */

                      PutPlotColors(YES, R, G, B);
                      scale = n - 3;

                      /* Put limits */

                      fmin = xx;
                      fmax = yy;

                      /* Print */

                      if (n == 4)
                        printf("\nCell importances, linear scale: ");
                      else if (n == 5)
                        printf("\nCell importances, log scale: ");
                      else if (n == 6)
                        printf("\nSource importances, log scale: ");
                      else if (n == 7)
                        printf("\nSource importances, log scale: ");
                      else
                        Die(FUNCTION_NAME, "WTF?");

                      printf("[ %1.5E %1.5E ]\n\n", fmin, fmax);

                      /* Refresh */

                      refresh = YES;
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "bou"))
            {
              /***************************************************************/

              /***** Toggle bouders ******************************************/

              /* Reset mode */

              n = -1000000000;

              /* Read mode */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                n = atol(word);

              /* Check */

              if (((n > -1) && (n < 4)) || (n == -1000000000))
                {
                  /* Check if given */

                  if (n != -1000000000)
                    {
                      /* Put value */

                      bou = n;

                      /* Put refresh flag */

                      refresh = YES;
                    }

                  /* Print */

                  if (bou == 0)
                    printf("\nBoundaries off\n\n");
                  else if (bou == 1)
                    printf("\nPlot cell boundaries\n\n");
                  else if (bou == 2)
                    printf("\nPlot material boundaries\n\n");
                  else if (bou == 3)
                    printf("\nPlot cell and material boundaries\n\n");
                  else
                    Die(FUNCTION_NAME, "Invalid boundary mode");
                }
              else
                printf("\nInvalid option: %ld\n\n", n);

              /***************************************************************/
            }
          else if (!strcasecmp(word, "pix"))
            {
              /***************************************************************/

              /***** Image size in pixels ************************************/

              /* Reset y-size */

              n = -1;
              m = -1;

              /* Read size */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                n = atol(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                m = atol(word);

              /* Set y-size if not given */

              if (m < 0)
                {
                  /* Check plot axis */

                  if (ax == PLOT_MODE_YZ)
                    m = (long)((zmax - zmin)/(ymax - ymin)*n);
                  else if (ax == PLOT_MODE_XZ)
                    m = (long)((zmax - zmin)/(xmax - xmin)*n);
                  else if (ax == PLOT_MODE_XY)
                    m = (long)((ymax - ymin)/(xmax - xmin)*n);
                  else
                    Die(FUNCTION_NAME, "Error in ax");
                }

              /* Check size */

              if ((n == -1) && (m == -1))
                printf("\nImage size: %ld x %ld pixels\n\n", xp, yp);
              else if ((n < 10) || (n > 10000) || (m < 10) || (m > 10000))
                printf("\nInvalid image size\n\n");
              else
                {
                  /* Put values */

                  xp = n;
                  yp = m;

                  /* Print */

                  printf("\nImage size set to: %ld x %ld pixels\n\n", xp, yp);

                  /* Set refresh flag */

                  refresh = YES;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "orig"))
            {
              /***************************************************************/

              /***** Set origin **********************************************/

              /* Reset sizes */

              xx = -INFTY;
              yy = -INFTY;
              zz = -INFTY;

              /* Read coordinates */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                xx = atof(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                yy = atof(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                zz = atof(word);

              /* Check dimensions */

              if (dim == 3)
                {
                  /* Check all 3 coordinates */

                  if ((xx == -INFTY) && (yy == -INFTY) && (zz == -INFTY))
                    printf("\nOrigin: [ %1.5f %1.5f %1.5f ]\n\n", x, y, z);
                  else if ((xx == -INFTY) || (yy == -INFTY) || (zz == -INFTY))
                    printf("\nMissing coordinates\n\n");
                  else
                    {
                      /* Put values */

                      x = xx;
                      y = yy;
                      z = zz;

                      /* Put boundaries */

                      xmin = x - 0.5*dx;
                      xmax = x + 0.5*dx;
                      ymin = y - 0.5*dy;
                      ymax = y + 0.5*dy;
                      zmin = z - 0.5*dz;
                      zmax = z + 0.5*dz;

                      printf("\nOrigin set to: [ %1.5f %1.5f %1.5f ]\n\n",
                             x, y, z);

                      /* Put refresh flag */

                      refresh = YES;
                    }
                }
              else
                {
                  /* Check all x- and y-coordinate */

                  if ((xx == -INFTY) && (yy == -INFTY))
                    printf("\nOrigin: [ %1.5f %1.5f %1.5f ]\n\n", x, y, z);
                  else if ((xx == -INFTY) || (yy == -INFTY))
                    printf("\nMissing coordinates\n\n");
                  else
                    {
                      /* Put values */

                      x = xx;
                      y = yy;

                      /* Put boundaries */

                      xmin = x - 0.5*dx;
                      xmax = x + 0.5*dx;
                      ymin = y - 0.5*dy;
                      ymax = y + 0.5*dy;

                      printf("\nOrigin set to: [ %1.5f %1.5f ]\n\n", x, y);

                      /* Put refresh flag */

                      refresh = YES;
                    }
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "rgb"))
            {
              /***************************************************************/

              /***** Set color value *****************************************/

              /* Reset variables */

              *mname = '\0';
              r = -1000000;
              g = -1000000;
              b = -1000000;

              /* Read material name */

              l = l + strlen(word) + 1;
              sscanf(&line[l], "%s", mname);

              /* Read RGB values */

              l = l + strlen(mname) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                r = atol(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                g = atol(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                b = atol(word);

              /* Check */

              if ((mname[0] == '\0') || (r == -1000000) || (g == -1000000) ||
                  (b == -1000000))
                {
                  /* Print */

                  printf("\nMaterial RGB values:\n\n");

                  /* Loop over materials */

                  mat = (long)RDB[DATA_PTR_M0];
                  while (mat > VALID_PTR)
                    {
                      /* Get color index */

                      n = (long)RDB[mat + MATERIAL_COLOUR_IDX];

                      /* Check */

                      if ((n > 4) && (n < MAX_PLOT_COLORS))
                        {
                          /* Print */

                          printf("%-10s : %3ld %3ld %3ld\n",
                                 GetText(mat + MATERIAL_PTR_NAME),
                                 R[n], G[n], B[n]);
                        }

                      /* Next material */

                      mat = NextItem(mat);
                    }

                  /* Newline */

                  printf("\n");
                }
              else if ((r < 0) || (r > 255) || (g < 0) || (g > 255) ||
                       (b < 0) || (b > 255))
                printf("\nInvalid values\n\n");
              else
                {
                  /* Reset count */

                  m = 0;

                  /* Find material */

                  mat = (long)RDB[DATA_PTR_M0];
                  while (mat > VALID_PTR)
                    {
                      /* Get pointer to parent */

                      mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];

                      /* Check */

                      if ((!strcmp(mname, GetText(mat + MATERIAL_PTR_NAME))) ||
                          ((mat0 > VALID_PTR) &&
                           !strcmp(mname, GetText(mat0 + MATERIAL_PTR_NAME))))
                        {
                          /* Get color index */

                          n = (long)RDB[mat + MATERIAL_COLOUR_IDX];

                          /* Check */

                          if ((n > 4) && (n < MAX_PLOT_COLORS))
                            {
                              /* Adjust */

                              R[n] = r;
                              G[n] = g;
                              B[n] = b;

                              /* Add count */

                              m++;
                            }
                        }

                      /* Next */

                      mat = NextItem(mat);
                    }

                  /* Check count */

                  if (m > 0)
                    {
                      /* Print */

                      printf("\nSetting material %s color to: %ld %ld %ld\n\n",
                             mname, R[n], G[n], B[n]);

                      /* Refresh */

                      refresh = YES;
                    }
                  else
                    printf("\nMaterial %s not found\n\n", mname);
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "max"))
            {
              /***************************************************************/

              /***** Reset to outer boundaries *******************************/

              /* Get geometry boundaries */

              xmin = RDB[DATA_GEOM_MINX];
              xmax = RDB[DATA_GEOM_MAXX];
              ymin = RDB[DATA_GEOM_MINY];
              ymax = RDB[DATA_GEOM_MAXY];
              zmin = RDB[DATA_GEOM_MINZ];
              zmax = RDB[DATA_GEOM_MAXZ];

              /* Check plot axis */

              if (ax == PLOT_MODE_YZ)
                yp = (long)((zmax - zmin)/(ymax - ymin)*xp);
              else if (ax == PLOT_MODE_XZ)
                yp = (long)((zmax - zmin)/(xmax - xmin)*xp);
              else if (ax == PLOT_MODE_XY)
                yp = (long)((ymax - ymin)/(xmax - xmin)*xp);
              else
                Die(FUNCTION_NAME, "Error in ax");

              /* Print */

              printf("\nReset to maximum dimensions:\n\n");

              if (ax == PLOT_MODE_YZ)
                printf("Y-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n",
                       ymin, ymax, zmin, zmax);
              else if (ax == PLOT_MODE_XZ)
                printf("X-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n",
                       xmin, xmax, zmin, zmax);
              else
                printf("X-dim  : [ %1.5f %1.5f ]\nY-dim  : [ %1.5f %1.5f ]\n",
                       xmin, xmax, ymin, ymax);

              printf("\nImage size set to: %ld x %ld pixels\n\n", xp, yp);

              /* Set refresh flag */

              refresh = YES;

              /***************************************************************/
            }
          else if (!strcasecmp(word, "sq"))
            {
              /***************************************************************/

              /***** Set square boundaries ***********************************/

              /* Check plot axis */

              if (ax == PLOT_MODE_YZ)
                {
                  /* Check larger */

                  if (dy < dz)
                    dy = dz;
                  else
                    dz = dy;
                }
              else if (ax == PLOT_MODE_XZ)
                {
                  /* Check larger */

                  if (dx < dz)
                    dx = dz;
                  else
                    dz = dx;
                }
              else if (ax == PLOT_MODE_XY)
                {
                  /* Check larger */

                  if (dx < dy)
                    dx = dy;
                  else
                    dy = dx;
                }
              else
                Die(FUNCTION_NAME, "WTF?");

              /* Put boundaries */

              xmin = x - 0.5*dx;
              xmax = x + 0.5*dx;
              ymin = y - 0.5*dy;
              ymax = y + 0.5*dy;
              zmin = z - 0.5*dz;
              zmax = z + 0.5*dz;

              /* Check plot axis */

              if (ax == PLOT_MODE_YZ)
                yp = (long)((zmax - zmin)/(ymax - ymin)*xp);
              else if (ax == PLOT_MODE_XZ)
                yp = (long)((zmax - zmin)/(xmax - xmin)*xp);
              else if (ax == PLOT_MODE_XY)
                yp = (long)((ymax - ymin)/(xmax - xmin)*xp);
              else
                Die(FUNCTION_NAME, "Error in ax");

              /* Print */

              printf("\nSquare dimensions:\n\n");

              if (ax == PLOT_MODE_YZ)
                printf("Y-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n",
                       ymin, ymax, zmin, zmax);
              else if (ax == PLOT_MODE_XZ)
                printf("X-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n",
                       xmin, xmax, zmin, zmax);
              else
                printf("X-dim  : [ %1.5f %1.5f ]\nY-dim  : [ %1.5f %1.5f ]\n",
                       xmin, xmax, ymin, ymax);

              printf("\nImage size set to: %ld x %ld pixels\n\n", xp, yp);

              /* Set refresh flag */

              refresh = YES;

              /***************************************************************/
            }
          else if (!strcasecmp(word, "box"))
            {
              /***************************************************************/

              /***** Plot area ***********************************************/

              /* Reset sizes */

              min1 = INFTY;
              max1 = -INFTY;
              min2 = INFTY;
              max2 = -INFTY;

              /* Read coordinates */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                min1 = atof(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                max1 = atof(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                min2 = atof(word);

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                max2 = atof(word);

              /* Check missing */

              if ((min1 == INFTY) && (max1 == -INFTY) &&
                       (min2 == INFTY) && (max2 == -INFTY))
                {
                  printf("\nGeometry dimensions:\n\n");

                  if (ax == PLOT_MODE_YZ)
                    printf("Y-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n\n",
                           ymin, ymax, zmin, zmax);
                  else if (ax == PLOT_MODE_XZ)
                    printf("X-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n\n",
                           xmin, xmax, zmin, zmax);
                  else
                    printf("X-dim  : [ %1.5f %1.5f ]\nY-dim  : [ %1.5f %1.5f ]\n\n",
                           xmin, xmax, ymin, ymax);
                }
              else if ((min1 == INFTY) || (max1 == -INFTY) ||
                       (min2 == INFTY) || (max2 == -INFTY))
                printf("\nMissing boundaries\n\n");
              else
                {
                  /* Check axis */

                  if (ax == PLOT_MODE_YZ)
                    {
                      ymin = min1;
                      ymax = max1;
                      zmin = min2;
                      zmax = max2;
                    }
                  else if (ax == PLOT_MODE_XZ)
                    {
                      xmin = min1;
                      xmax = max1;
                      zmin = min2;
                      zmax = max2;
                    }
                  else
                    {
                      xmin = min1;
                      xmax = max1;
                      ymin = min2;
                      ymax = max2;
                    }

                  /* Check boundaries and swap */

                  if (xmin > xmax)
                    {
                      tmp = xmax;
                      xmax = xmin;
                      xmin = tmp;
                    }

                  if (ymin > ymax)
                    {
                      tmp = ymax;
                      ymax = ymin;
                      ymin = tmp;
                    }

                  if (zmin > zmax)
                    {
                      tmp = zmax;
                      zmax = zmin;
                      zmin = tmp;
                    }

                  /* Check plot axis */

                  if (ax == PLOT_MODE_YZ)
                    yp = (long)((zmax - zmin)/(ymax - ymin)*xp);
                  else if (ax == PLOT_MODE_XZ)
                    yp = (long)((zmax - zmin)/(xmax - xmin)*xp);
                  else if (ax == PLOT_MODE_XY)
                    yp = (long)((ymax - ymin)/(xmax - xmin)*xp);
                  else
                    Die(FUNCTION_NAME, "Error in ax");

                  printf("\nGeometry dimensions set to:\n\n");

                  if (ax == PLOT_MODE_YZ)
                    printf("Y-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n",
                           ymin, ymax, zmin, zmax);
                  else if (ax == PLOT_MODE_XZ)
                    printf("X-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n",
                           xmin, xmax, zmin, zmax);
                  else
                    printf("X-dim  : [ %1.5f %1.5f ]\nY-dim  : [ %1.5f %1.5f ]\n",
                           xmin, xmax, ymin, ymax);

                  printf("\nImage size set to: %ld x %ld pixels\n\n", xp, yp);

                  /* Put refresh flag */

                  refresh = YES;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "ax") || !strcasecmp(word, "axis") ||
                   !strcasecmp(word, "xz") || !strcasecmp(word, "zx") ||
                   !strcasecmp(word, "yz") || !strcasecmp(word, "zy") ||
                   !strcasecmp(word, "xy") || !strcasecmp(word, "yx"))
            {
              /***************************************************************/

              /***** Plot axis ***********************************************/

              /* Check short notation */

              if (!strcasecmp(word, "xz") || !strcasecmp(word, "zx"))
                {
                  n = l + strlen(word) + 1;
                  sprintf(&line[n], " xz");
                }
              else if (!strcasecmp(word, "yz") || !strcasecmp(word, "zy"))
                {
                  n = l + strlen(word) + 1;
                  sprintf(&line[n], " yz");
                }
              else if (!strcasecmp(word, "xy") || !strcasecmp(word, "yx"))
                {
                  n = l + strlen(word) + 1;
                  sprintf(&line[n], " xy");
                }

              /* Read axis */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%s", word) > 0)
                {
                  /* Check */

                  if (!strcasecmp(word, "yz"))
                    n = PLOT_MODE_YZ;
                  else if (!strcasecmp(word, "xz"))
                    n = PLOT_MODE_XZ;
                  else if (!strcasecmp(word, "xy"))
                    n = PLOT_MODE_XY;
                  else
                    n = atol(word);

                  /* Check */

                  if ((n < 1) || (n > 3))
                    printf("\nInvalid axis: %ld\n", n);
                  else if ((dim != 3) &&
                           (n != PLOT_MODE_XY))
                    printf("\nAxis type not allowed for 2D geometry\n");
                  else
                    {
                      /* Put value */

                      ax = n;

                      /* Check plot axis */

                      if (ax == PLOT_MODE_YZ)
                        yp = (long)((zmax - zmin)/(ymax - ymin)*xp);
                      else if (ax == PLOT_MODE_XZ)
                        yp = (long)((zmax - zmin)/(xmax - xmin)*xp);
                      else if (ax == PLOT_MODE_XY)
                        yp = (long)((ymax - ymin)/(xmax - xmin)*xp);
                      else
                        Die(FUNCTION_NAME, "Error in ax");

                      /* Print */

                      if (ax == PLOT_MODE_YZ)
                        printf("\nImage axis set to: YZ\n");
                      if (ax == PLOT_MODE_XZ)
                        printf("\nImage axis set to: XZ\n");
                      if (ax == PLOT_MODE_XY)
                        printf("\nImage axis set to: XY\n");

                      printf("\nImage size set to: %ld x %ld pixels\n\n",
                             xp, yp);

                      /* Set refresh flag */

                      refresh = YES;
                    }
                }
              else
                {
                  if (ax == PLOT_MODE_YZ)
                    printf("\nImage axis: YZ\n\n");
                  if (ax == PLOT_MODE_XZ)
                    printf("\nImage axis: XZ\n\n");
                  if (ax == PLOT_MODE_XY)
                    printf("\nImage axis: XY\n\n");
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "zi") || !strcasecmp(word, "zo"))
            {
              /***************************************************************/

              /***** Zoom in / out********************************************/

              /* Check mode */

              if (!strcasecmp(word, "zo"))
                {
                  /* Read factor */

                  l = l + strlen(word) + 1;
                  if (sscanf(&line[l], "%lf", &f) > 0)
                    {
                      /* Put factor */

                      if (f < 1.0)
                        printf("\nInvalid zoom-factor: %f\n", f);
                    }
                  else
                    {
                      /* Put default */

                      f = 1.1;
                    }
                }
              else
                {
                  /* Read factor */

                  l = l + strlen(word) + 1;
                  if (sscanf(&line[l], "%lf", &f) > 0)
                    {
                      /* Put factor */

                      if (f < 1.0)
                        printf("\nInvalid zoom-factor: %f\n", f);
                      else
                        f = 1.0/f;
                    }
                  else
                    {
                      /* Put default */

                      f = 1.0/1.1;
                    }
                }

              /* Check */

              if (f > 0.0)
                {
                  /* Print */

                  if (f > 1.0)
                    printf("\nZoom out by factor %1.1f\n\n", f);
                  else
                    printf("\nZoom in by factor %1.1f\n\n", 1.0/f);

                  /* Adjust widths */

                  dx = dx*f;
                  dy = dy*f;
                  dz = dz*f;

                  /* Put boundaries */

                  xmin = x - 0.5*dx;
                  xmax = x + 0.5*dx;
                  ymin = y - 0.5*dy;
                  ymax = y + 0.5*dy;
                  zmin = z - 0.5*dz;
                  zmax = z + 0.5*dz;

                  if (ax == PLOT_MODE_YZ)
                    printf("Y-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n\n",
                       ymin, ymax, zmin, zmax);
                  else if (ax == PLOT_MODE_XZ)
                    printf("X-dim  : [ %1.5f %1.5f ]\nZ-dim  : [ %1.5f %1.5f ]\n\n",
                           xmin, xmax, zmin, zmax);
                  else
                    printf("X-dim  : [ %1.5f %1.5f ]\nY-dim  : [ %1.5f %1.5f ]\n\n",
                           xmin, xmax, ymin, ymax);

                  /* Put refresh flag */

                  refresh = YES;
                }

              /***************************************************************/
            }
          else if (!strcasecmp(word, "n") || !strcasecmp(word, "s") ||
                   !strcasecmp(word, "e") || !strcasecmp(word, "w") ||
                   !strcasecmp(word, "u") || !strcasecmp(word, "d"))
            {
              /***************************************************************/

              /***** Move origin *********************************************/

              /* Avoid compiler warning */

              m = -1;

              /* Put direction */

              if (!strcasecmp(word, "w"))
                m = 1;
              else if (!strcasecmp(word, "s"))
                m = 2;
              else if (!strcasecmp(word, "e"))
                m = 3;
              else if (!strcasecmp(word, "n"))
                m = 4;
              else if (!strcasecmp(word, "u"))
                m = 5;
              else if (!strcasecmp(word, "d"))
                m = 6;
              else
                Die(FUNCTION_NAME, "WTF");

              /* Read factor */

              l = l + strlen(word) + 1;
              if (sscanf(&line[l], "%lf", &f) > 0)
                {
                  /* Move by given amount */

                  if (f < 0.0)
                    printf("\nInvalid move factor: %f\n", f);
                  else
                    {
                      /* Check plot axis */

                      if (ax == PLOT_MODE_YZ)
                        {
                          /* Move origin */

                          if (m == 1)
                            y = y - f;
                          else if (m == 2)
                            z = z - f;
                          else if (m == 3)
                            y = y + f;
                          else if (m == 4)
                            z = z + f;
                          else if (m == 5)
                            x = x + f;
                          else if (m == 6)
                            x = x - f;
                        }
                      else if (ax == PLOT_MODE_XZ)
                        {
                          /* Move origin */

                          if (m == 1)
                            x = x - f;
                          else if (m == 2)
                            z = z - f;
                          else if (m == 3)
                            x = x + f;
                          else if (m == 4)
                            z = z + f;
                          else if (m == 5)
                            y = y + f;
                          else if (m == 6)
                            y = y - f;
                        }
                      else if (ax == PLOT_MODE_XY)
                        {
                          /* Move origin */

                          if (m == 1)
                            x = x - f;
                          else if (m == 2)
                            y = y - f;
                          else if (m == 3)
                            x = x + f;
                          else if (m == 4)
                            y = y + f;
                          else if (m == 5)
                            z = z + f;
                          else if (m == 6)
                            z = z - f;
                        }
                    }
                }
              else
                {
                  /* Move by 10% of width */

                  f = 0.1;

                  /* Check plot axis */

                  if (ax == PLOT_MODE_YZ)
                    {
                      /* Move origin */

                      if (m == 1)
                        y = y - f*dy;
                      else if (m == 2)
                        z = z - f*dz;
                      else if (m == 3)
                        y = y + f*dy;
                      else if (m == 4)
                        z = z + f*dz;
                      else if (m == 5)
                        x = x + f*dx;
                      else if (m == 6)
                        x = x - f*dx;
                    }
                  else if (ax == PLOT_MODE_XZ)
                    {
                      /* Move origin */

                      if (m == 1)
                        x = x - f*dx;
                      else if (m == 2)
                        z = z - f*dz;
                      else if (m == 3)
                        x = x + f*dx;
                      else if (m == 4)
                        z = z + f*dz;
                      else if (m == 5)
                        y = y + f*dy;
                      else if (m == 6)
                        y = y - f*dy;
                    }
                  else if (ax == PLOT_MODE_XY)
                    {
                      /* Move origin */

                      if (m == 1)
                        x = x - f*dx;
                      else if (m == 2)
                        y = y - f*dy;
                      else if (m == 3)
                        x = x + f*dx;
                      else if (m == 4)
                        y = y + f*dy;
                      else if (m == 5)
                        z = z + f*dz;
                      else if (m == 6)
                        z = z - f*dz;
                    }
                }

              /* Check */

              if (f > 0.0)
                {
                  /* Put boundaries */

                  xmin = x - 0.5*dx;
                  xmax = x + 0.5*dx;
                  ymin = y - 0.5*dy;
                  ymax = y + 0.5*dy;
                  zmin = z - 0.5*dz;
                  zmax = z + 0.5*dz;

                  if (dim == 3)
                    printf("\nOrigin set to: [ %1.5f %1.5f %1.5f]\n\n",
                           x, y, z);
                  else
                    printf("\nOrigin set to: [ %1.5f %1.5f ]\n\n", x, y);

                  /* Put refresh flag */

                  refresh = YES;
                }

              /***************************************************************/
            }
          else
            {
              /***************************************************************/

              /***** Unknown command *****************************************/

              /* Print note */

              printf("\nvalid command: %s\n\n", word);

              /***************************************************************/
            }
        }

      /* Check exit flag */

      if (stop == YES)
        break;

      /* Check refresh-flag */

      if (refresh == NO)
        continue;

      /* Put position */

      if (ax == PLOT_MODE_YZ)
        pos = x;
      else if (ax == PLOT_MODE_XZ)
        pos = y;
      else if (ax == PLOT_MODE_XY)
        pos = z;

      /***********************************************************************/

      /***** Draw image ******************************************************/

      /* Allocate memory */

      mtx = (long *)Mem(MEM_ALLOC, xp*yp, sizeof(long));

      /* Plot image */

      n = PlotImage(gpl, mtx, xp, yp, xmin, xmax, ymin, ymax, zmin, zmax, pos,
                    bou, par, ax, scale, E, fmin, fmax, qp);

      /* Check for geometry errors */

      if (n < 0)
        {
          if (n == GEOM_ERROR_NO_CELL)
            printf("Geometry errors (no cell).\n\n");
          else if (n == GEOM_ERROR_MULTIPLE_CELLS)
            printf("Geometry errors (overlap).\n\n");
          else
            printf("Geometry errors (unknown).\n\n");
        }

      /* Check matrix */

      for (n = 0; n < xp; n++)
        for (m = 0; m < yp; m++)
          if ((mtx[n*yp + m] < 0) || (mtx[n*yp + m] > MAX_PLOT_COLORS - 1))
            Die(FUNCTION_NAME, "Invalid colour index %ld (%ld %ld)",
                mtx[n*yp + m], n, m);

      /* Create image */

      im = gdImageCreate(xp, yp);

      /* Generate palette */

      for (n = 0; n < MAX_PLOT_COLORS; n++)
        palette[n] = gdImageColorAllocate(im, R[n], G[n], B[n]);

      /* Draw image */

#ifdef OPEN_MP
#pragma omp parallel private(n, m)
#endif
      {

#ifdef OPEN_MP
#pragma omp for
#endif
        for (n = 0; n < xp; n++)
          for (m = 0; m < yp; m++)
            {
              /* Flip y-axis in xy-plot */

              if (ax == PLOT_MODE_XY)
                gdImageSetPixel(im, n, m, palette[mtx[n*yp + yp - m - 1]]);
              else
                gdImageSetPixel(im, n, m, palette[mtx[n*yp + m]]);
            }
      }

      /* Open file for writing */

      if ((fp = fopen(fname, "w")) == NULL)
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Write image (png format) */

      gdImagePng(im, fp);

      /* Free image */

      gdImageDestroy(im);

      /* Free matrix */

      Mem(MEM_FREE, mtx);

      /* Close file */

      fclose(fp);

      /***********************************************************************/
    }

  /***************************************************************************/


  /* Exit subroutine */

  exit(-1);

#endif
}

/*****************************************************************************/
