/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : disperse2.c                                    */
/*                                                                           */
/* Created:       2013/08/26 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Experimental version of the disperser routine                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Disperse2:"

/*****************************************************************************/

void Disperse2()
{
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w;
  long uni0, cell0, root, cell, id, n;
  unsigned long seed;
  char uname[MAX_STR], cname[MAX_STR];

  /***************************************************************************/

  /***** Prompt input data ***************************************************/
  return;
  /* Get universe */
  
  fprintf(outp, "Enter universe:\n");
    
  if (scanf("%s", uname) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Get cell */
  
  fprintf(outp, "Enter cell:\n");
    
  if (scanf("%s", cname) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /***************************************************************************/

  /***** Find cell and universe pointers *************************************/

  /* Loop over universes */

  uni0 = (long)RDB[DATA_PTR_U0];
  while (uni0 > VALID_PTR)
    {
      /* Compare name */

      if (!strcmp(GetText(uni0 + UNIVERSE_PTR_NAME), uname))
        break;

      /* Next */

      uni0 = NextItem(uni0);
    }

  /* Check pointer */

  if (uni0 < VALID_PTR)
    Error(0, "Universe %s not found in geometry", uname);

  /* Loop over cells */

  cell0 = (long)RDB[DATA_PTR_C0];
  while (cell0 > VALID_PTR)
    {
      /* Compare name */

      if (!strcmp(GetText(cell0 + CELL_PTR_NAME), cname))
        break;

      /* Next */

      cell0 = NextItem(cell0);
    }

  /* Check pointer */

  if (cell0 < VALID_PTR)
    Error(0, "Cell %s not found in geometry", cname);

  /* Set thread id */

  id = 0;

  /* Init random number sequence */
      
  seed = ReInitRNG(1);
  SEED[0] = seed;

  /* Set plotter mode to avoid termination in geometry error, and quick */
  /* plotter mode to disable overlap check */

  WDB[DATA_PLOTTER_MODE] = (double)YES;
  WDB[DATA_QUICK_PLOT_MODE] = (double)YES;

  /* Override root universe pointer */

  root = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
  WDB[DATA_PTR_ROOT_UNIVERSE] = (double)uni0;

  /***************************************************************************/

  /***** Sample code *********************************************************/

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(uni0)", DATA_ARRAY, uni0);
  CheckPointer(FUNCTION_NAME, "(cell0)", DATA_ARRAY, cell0);

  /* Get universe boundaries */

  xmin = RDB[uni0 + UNIVERSE_MINX];
  xmax = RDB[uni0 + UNIVERSE_MAXX];
  ymin = RDB[uni0 + UNIVERSE_MINY];
  ymax = RDB[uni0 + UNIVERSE_MAXY];
  zmin = RDB[uni0 + UNIVERSE_MINZ];
  zmax = RDB[uni0 + UNIVERSE_MAXZ];
  
  fprintf(outp, "\nboundaries : \nx = [%E,%E] \ny = [%E %E] \nz = [%E %E]\n\n",
          xmin, xmax, ymin, ymax, zmin, zmax);

  /* Set direction vector (needed but not used) */

  u = 1.0;
  v = 0.0;
  w = 0.0;

  /* Sample 1000 random points and print coordinates that are inside the */
  /* universe and cell */

  for (n = 0; n < 1000; n++)
    {
      /* Sample random point */

      x = RandF(id)*(xmax - xmin) + xmin;
      y = RandF(id)*(ymax - ymin) + ymin;
      z = RandF(id)*(zmax - zmin) + zmin;

      /* Get cell pointer at position */

      cell = WhereAmI(x, y, z, u, v, w, id);

      /* Check pointer and print coordinates */

      if (cell == cell0)
        fprintf(outp, "%E %E %E\n", x, y, z);
    }

  /***************************************************************************/

  /* Reset plotter mode */

  WDB[DATA_PLOTTER_MODE] = (double)YES;

  /* Put root universe pointer */

  WDB[DATA_PTR_ROOT_UNIVERSE] = (double)root;

  /* Terminate calculation */

  exit(0);
}

/*****************************************************************************/
