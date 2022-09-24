/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcells.c                                 */
/*                                                                           */
/* Created:       2010/10/11 (JLe)                                           */
/* Last modified: 2019/05/13 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Creates surface lists, etc.                                  */
/*                                                                           */
/* Comments: - Tota surface pointterien hakurutiinia muutettiin ottamalla    */
/*             toinen nimi loopin ulkopuolelle (30.8.2017 / 2.1.30 / JLe).   */
/*                                                                           */
/*           - Haku pinnan numeroarvolla lisättiin 6.8.2017 / 2.1.30 / JLe.  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessCells:"

/* Väliaikainen testifunktio */

void Testaa(long, long *, long);

/*****************************************************************************/

void ProcessCells()
{
  long cell, ptr, infix[10000], surf, loc0, nspec, nmax, uni, new, loc1, n;
  char *sname, tmpstr1[MAX_STR], tmpstr2[MAX_STR];
  double r;

  fprintf(outp, "Processing cells...\n");

  /* Avoid compiler warning */

  nmax = -1;

  /* Loop over cells and add default infinite surfaces */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Check if cell has no surface list */

      if ((long)RDB[cell + CELL_PTR_SURF_LIST] < VALID_PTR)
        {
          /* Create an infinite surface */

          surf = NewItem(DATA_PTR_S0, SURFACE_BLOCK_SIZE);

          /* Put name */

          sprintf(tmpstr1, "%s_inf", GetText(cell + CELL_PTR_NAME));
          WDB[surf + SURFACE_PTR_NAME] = (double)PutText(tmpstr1);

          /* Put number of parameters and type */

          WDB[surf + SURFACE_N_PARAMS] = 0.0;
          WDB[surf + SURFACE_TYPE] = (double)SURF_INF;

          /* Create surface list for cell */

          sprintf(tmpstr2, " - %s", tmpstr1);
          WDB[cell + CELL_PTR_SURF_LIST] = (double)PutText(tmpstr2);
        }

      /* Next cell */

      cell = NextItem(cell);
    }

  /* Close and sort surface list */

  if ((surf = (long)RDB[DATA_PTR_S0]) > VALID_PTR)
    {
      /* Close list */

      CloseList(surf);

      /* Sort in ascending order */

      SortList(surf, SURFACE_NAME_NUMERIC, SORT_MODE_ASCEND);
    }

  /***************************************************************************/

  /***** Process surface lists according to format ***************************/

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Add counter */

      WDB[DATA_N_TOT_CELLS] = RDB[DATA_N_TOT_CELLS] + 1.0;

      /* Get parameters in infix format */

      nmax = ReadInfix(cell, infix, &nspec);

      /* Check special count (unions and parenthesis) */

      if (nspec == 0)
        {
          /* Only intersections in list */

          IntersectionList(cell, infix, nmax);
        }
      else
        {
          /* Convert list to postfix */

          ShuntingYard(cell, infix, nmax);

          /* Add counter */

          WDB[DATA_N_UNION_CELLS] = RDB[DATA_N_UNION_CELLS] + 1.0;

          /* Test notation */
          /*
          Testaa(cell, infix, nmax);
          */
        }

#ifndef TEST

      /***********************************************************************/

      /***** Link surfaces ***************************************************/

      /* Allocate memory for surface list */

      loc0 = ReallocMem(DATA_ARRAY, nmax + 1);
      WDB[cell + CELL_PTR_SURF_LIST] = (double)loc0;

      /* Check type */

      if ((ptr = (long)RDB[cell + CELL_PTR_SURF_COMP]) > VALID_PTR)
        {
          /*******************************************************************/

          /***** Composition list given **************************************/

          /* Loop over list */

          while ((long)RDB[ptr] != 0)
            {
              /* Check type */

              if ((long)RDB[ptr] > 0)
                {
                  /* Get surface name (no memory allocated in loop) */

                  sname = GetText(ptr);

                  /* Check if numeric value */

                  if ((n = NumericStr(sname)) > 0)
                    {
                      /* Pointer to list */

                      surf = (long)RDB[DATA_PTR_S0];
                      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

                      /* Find by binary search */

                      surf = SeekList(surf, SURFACE_NAME_NUMERIC, (double)n,
                                      SORT_MODE_ASCEND);
                    }
                  else
                    surf = -1;

                  /* Check if found */

                  if (surf < VALID_PTR)
                    {
                      /* Loop over surfaces */

                      surf = (long)RDB[DATA_PTR_S0];
                      while (surf > VALID_PTR)
                        {
                          /* Compare */

                          if (!strcmp(sname, GetText(surf + SURFACE_PTR_NAME)))
                            {
                              /* Break loop */

                              break;
                            }

                          /* Next */

                          surf = NextItem(surf);
                        }
                    }

                  /* Check */

                  if (surf < 0)
                    Error(cell, "Surface %s is not defined", sname);
                  else
                    {
                      /* Put pointers */

                      WDB[ptr] = (double)surf;
                      WDB[loc0++] = (double)surf;
                    }
                }

              /* Next */

              ptr++;
            }

          /* Put null terminator for surface list */

          WDB[loc0] = NULLPTR;

          /*******************************************************************/
        }
      else if ((ptr = (long)RDB[cell + CELL_PTR_SURF_INSC]) > VALID_PTR)
        {
          /*******************************************************************/

          /***** Intersection list given *************************************/

          /* Close list */

          CloseList(ptr);

          /* Loop over items */

          while (ptr > VALID_PTR)
            {
              /* Get surface name (no memory allocated in loop) */

              sname = GetText(ptr + CELL_INSC_PTR_SURF);

              /* Check if numeric value */

              if ((n = NumericStr(sname)) > 0)
                {
                  /* Pointer to list */

                  surf = (long)RDB[DATA_PTR_S0];
                  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

                  /* Find by binary search */

                  surf = SeekList(surf, SURFACE_NAME_NUMERIC, (double)n,
                                  SORT_MODE_ASCEND);
                }
              else
                surf = -1;

              /* Check if found */

              if (surf < VALID_PTR)
                {
                  /* Loop over surfaces */

                  surf = (long)RDB[DATA_PTR_S0];
                  while (surf > VALID_PTR)
                    {
                      /* Compare */

                      if (!strcmp(sname, GetText(surf + SURFACE_PTR_NAME)))
                        {
                          /* Break loop */

                          break;
                        }

                      /* Next */

                      surf = NextItem(surf);
                    }
                }

              /* Check */

              if (surf < 0)
                Error(cell, "Surface %s is not defined", sname);
              else
                {
                  /* Put pointers */

                  WDB[ptr + CELL_INSC_PTR_SURF] = (double)surf;
                  WDB[loc0++] = (double)surf;
                }

              /* Next */

              ptr = NextItem(ptr);
            }

          /* Put null terminator for surface list */

          WDB[loc0] = (double)NULLPTR;

          /*******************************************************************/
        }
      else
        Die(FUNCTION_NAME, "No surface list");

      /***********************************************************************/

#endif

      /* Next cell */

      cell = NextItem(cell);
    }

  /* Re-open surface list */

  if ((surf = (long)RDB[DATA_PTR_S0]) > VALID_PTR)
    ReopenList(surf);

  /***************************************************************************/

  /***** Additional checks for some surface types ****************************/

  /* Loop over surfaces */

  surf = (long)RDB[DATA_PTR_S0];
  while (surf > VALID_PTR)
    {
      /* Check torus */

      if (((long)RDB[surf + SURFACE_TYPE] == SURF_TORX) ||
          ((long)RDB[surf + SURFACE_TYPE] == SURF_TORY) ||
          ((long)RDB[surf + SURFACE_TYPE] == SURF_TORZ))
        {
          /* Pointer to parameters */

          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check radii */

          if (RDB[ptr + 3] < 0.0)
            Error(surf, "Invalid radius %E", RDB[ptr + 3]);
          else if (RDB[ptr + 3] == 0.0)
            {
              /* Convert to ellipsoid (nyt sphere) */

              if (RDB[ptr + 4] == RDB[ptr + 5])
                {
                  WDB[surf + SURFACE_TYPE] = (double)SURF_SPH;
                  WDB[surf + SURFACE_N_PARAMS] = 4.0;
                  WDB[ptr + 3] = RDB[ptr + 4];
                }
            }
          else if (RDB[ptr + 4] < ZERO)
            Error(surf, "Invalid radius %E", RDB[ptr + 4]);
          else if (RDB[ptr + 5] < ZERO)
            Error(surf, "Invalid radius %E", RDB[ptr + 5]);
        }

      /* Cylider along vector */

      else if ((long)RDB[surf + SURFACE_TYPE] == SURF_CYLV)
        {
          /* Pointer to parameters */

          ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Calculate length of direction vector */

          if ((r = sqrt(RDB[ptr + 3]*RDB[ptr + 3] + RDB[ptr + 4]*RDB[ptr + 4] +
                        RDB[ptr + 5]*RDB[ptr + 5])) == 0.0)
            Error(surf, "Invalid direction vector");

          /* Normalize */

          WDB[ptr + 3] = RDB[ptr + 3]/r;
          WDB[ptr + 4] = RDB[ptr + 4]/r;
          WDB[ptr + 5] = RDB[ptr + 5]/r;
        }

      /* Next */

      surf = NextItem(surf);
    }

  /***************************************************************************/

  /***** Create search lists *************************************************/

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Skip this loop to enforce using universe-based lists      */
      /* (NOTE: the cell-based lists were added 8.1.2016 / 2.1.28) */

      /*
      break;
      */
      /* Pointer to cell list */

      if ((loc0 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST]) > VALID_PTR)
        {
          /* Check size and skip universes with large number of cells */
          /* to avoid CPU time wasted for sorting */

          if (ListSize(loc0) > (long)RDB[DATA_MAX_CELL_SEARCH_LIST])
            {
              /* Print warning */

              Note(0, "Search lists not used for universe %s (%ld cells)",
                   GetText(uni + UNIVERSE_PTR_NAME), ListSize(loc0));

              /* Next universe */

              uni = NextItem(uni);

              /* Cycle loop */

              continue;
            }
        }

      /* Loop over cell list */

      while (loc0 > VALID_PTR)
        {
          /* Pointer to cell */

          cell = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Second loop over universe cell list (myös cell itse lisätään */
          /* jotta myöhemmin ei tuu mystistä pointer erroria) */

          loc1 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
          while (loc1 > VALID_PTR)
            {
              /* Allocate memory */

              new = NewItem(cell + CELL_PTR_SEARCH_LIST, CELL_LIST_BLOCK_SIZE);

              /* Copy data */

              memcpy(&WDB[new + LIST_DATA_SIZE], &RDB[loc1 + LIST_DATA_SIZE],
                     (CELL_LIST_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

              /* Allocate memory for counter */

              ptr = AllocPrivateData(1, PRIVA_ARRAY);
              WDB[new + CELL_LIST_PTR_COUNT] = (double)ptr;

              /* Next */

              loc1 = NextItem(loc1);
            }

          /* Close list */

          if ((ptr = (long)RDB[cell + CELL_PTR_SEARCH_LIST]) < VALID_PTR)
            Die(FUNCTION_NAME, "List not created");
          else
            CloseList(ptr);

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Next universe */

      uni = NextItem(uni);
    }

  /***************************************************************************/

  /***** Create surface search lists *****************************************/

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* NOTE: This is disabled for now */

      break;

      /* Loop over cell list */

      loc0 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
      while (loc0 > VALID_PTR)
        {
          /* Pointer to cell */

          cell = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Loop over surfaces */

          if ((loc1 = (long)RDB[cell + CELL_PTR_SURF_LIST]) > VALID_PTR)
            while ((surf = (long)RDB[loc1++]) > VALID_PTR)
              {
                /* Allocate memory */

                new = NewItem(surf + SURFACE_PTR_CELL_LIST,
                              CELL_LIST_BLOCK_SIZE);

                /* Copy data */

                memcpy(&WDB[new + LIST_DATA_SIZE], &RDB[loc0 + LIST_DATA_SIZE],
                       (CELL_LIST_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

                /* Allocate memory for counter */

                ptr = AllocPrivateData(1, PRIVA_ARRAY);
                WDB[new + CELL_LIST_PTR_COUNT] = (double)ptr;
              }

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Next universe */

      uni = NextItem(uni);
    }

  /* Loop over surfaces */

  surf = (long)RDB[DATA_PTR_S0];
  while (surf > VALID_PTR)
    {
      /* Close list */

      if ((ptr = (long)RDB[surf + SURFACE_PTR_CELL_LIST]) > VALID_PTR)
        CloseList(ptr);

      /* Next surface */

      surf = NextItem(surf);
    }

  /***************************************************************************/

  /* Exit OK */

  fprintf(outp, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/

void Testaa(long cell, long *infix, long ni)
{
  long n, ptr, stack[10000], ns, a, b;

  fprintf(outp, "\ncell %s: \n", GetText(cell + CELL_PTR_NAME));

  /* Print infix */

  fprintf(outp, "infix:   ");

  for (n = 0; n < ni; n++)
    {
      if (infix[n] == SURF_OP_OR)
        fprintf(outp, " + ");
      else if (infix[n] == SURF_OP_AND)
        fprintf(outp, " * ");
      else if (infix[n] == SURF_OP_NOT)
        fprintf(outp, " - ");
      else if (infix[n] == SURF_OP_LEFT)
        fprintf(outp, " ( ");
      else if (infix[n] == SURF_OP_RIGHT)
        fprintf(outp, " ) ");
      else
        {
          WDB[DATA_DUMMY] = (double)infix[n];
          fprintf(outp, " %s ", GetText(DATA_DUMMY));
        }
    }

  /* Print postfix */

  fprintf(outp, "\n");

  fprintf(outp, "postfix: ");

  ptr = (long)RDB[cell + CELL_PTR_SURF_COMP];
  while ((long)RDB[ptr] != 0)
    {
      if ((long)RDB[ptr] == SURF_OP_OR)
        fprintf(outp, " + ");
      else if ((long)RDB[ptr] == SURF_OP_AND)
        fprintf(outp, " * ");
      else if ((long)RDB[ptr] == SURF_OP_NOT)
        fprintf(outp, " - ");
      else if ((long)RDB[ptr] == SURF_OP_LEFT)
        fprintf(outp, " ( ");
      else if ((long)RDB[ptr] == SURF_OP_RIGHT)
        fprintf(outp, " ) ");
      else
        fprintf(outp, " %s ", GetText(ptr));

      /* Next */

      ptr++;
    }

  fprintf(outp, "\n");

  /* Evaluate */

  ns = 0;

  ptr = (long)RDB[cell + CELL_PTR_SURF_COMP];
  while ((long)RDB[ptr] != 0)
    {
      if ((long)RDB[ptr] == SURF_OP_OR)
        {
          /* Pop last two values */

          a = stack[--ns];
          b = stack[--ns];

          /* Push union to stack */

          stack[ns++] = a + b;
        }
      else if ((long)RDB[ptr] == SURF_OP_NOT)
        {
          /* Pop last value */

          a = stack[--ns];

          /* Push complement to stack */

          stack[ns++] = -a;
        }
      else if ((long)RDB[ptr] == SURF_OP_AND)
        {
          /* Pop last two values */

          a = stack[--ns];
          b = stack[--ns];

          /* Push intersection to stack */

          stack[ns++] = a * b;
        }
      else
        {
          /* Push value to stack */

          stack[ns++] = atoi(GetText(ptr));
        }

      ptr++;
    }

  fprintf(outp, "result:   %ld\n", stack[0]);
}

/*****************************************************************************/
