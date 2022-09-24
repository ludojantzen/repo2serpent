/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : Processbc.c                                    */
/*                                                                           */
/* Created:       2013/03/07 (JLe)                                           */
/* Last modified: 2018/08/10 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Creates a list of outer boundaries and finds surface used    */
/*              with repeated boundary conditions.                           */
/*                                                                           */
/* Comments: - List is used in stopatboundary.c, surface pointer in          */
/*             boundaryconditions.c                                          */
/*                                                                           */
/*           - NOTE: symmetrioiden käyttö tarkastetaan myös muualla.         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessBC:"

/*****************************************************************************/

void ProcessBC()
{
  long uni, type, loc0, loc1, side, nst, reg, surf, cell, ptr, tra, ok;

  /* Stop track at outer boundary if any of the BC's is black or albedos */
  /* are set. */

  if (((long)RDB[DATA_GEOM_BC1] == BC_BLACK) ||
      ((long)RDB[DATA_GEOM_BC2] == BC_BLACK) ||
      ((long)RDB[DATA_GEOM_BC3] == BC_BLACK) ||
      (RDB[DATA_GEOM_ALBEDO1] != 1.0) ||
      (RDB[DATA_GEOM_ALBEDO2] != 1.0) ||
      (RDB[DATA_GEOM_ALBEDO3] != 1.0))
    WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

  /* Get pointer to root universe */

  uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Reset OK flag */

  ok = NO;

  /* Get universe type */

  type = (long)RDB[uni + UNIVERSE_TYPE];    

  /* Get pointer to surfaces */

  if (type == UNIVERSE_TYPE_NEST)
    {
      /* Pointer to nest */

      nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
      CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

      /* Pointer to regions */

      reg = (long)RDB[nst + NEST_PTR_REGIONS];
      CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

      /* Pointer to last */

      reg = LastItem(reg);
      CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

      /* Pointer to surface */

      if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT]) < VALID_PTR)
        Error(nst, "Root universe must have at least two regions");

      /* Get surface type */

      type = (long)RDB[surf + SURFACE_TYPE];
      
      /* Check allowed types for repeated boundary conditions */
      
      if ((type == SURF_SQC) || (type == SURF_CUBE) || 
          (type == SURF_CUBOID) || (type == SURF_HEXYC) || 
          (type == SURF_HEXXC) || (type == SURF_HEXYPRISM) || 
          (type == SURF_HEXXPRISM) || (type == SURF_RECT))
        {
          /* Check fill pointer */

          if ((long)RDB[reg + NEST_REG_PTR_FILL] > VALID_PTR)
            Error(nst, "Outermost region cannot be filled with a universe");

          /* Pointer to cell */

          cell = (long)RDB[reg + NEST_REG_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

          /* Put pointer to nest structure */
          
          WDB[nst + NEST_PTR_BC_SURF] = (double)surf;

          /* Put pointer to cell structure */
                      
          WDB[cell + CELL_PTR_BC_SURF] = (double)surf;          

          /* Set OK flag */

          ok = YES;
        }

      /* Add pointer to list of boundaries */

      ptr = NewItem(DATA_PTR_OUTER_BOUNDS, BOUNDS_BLOCK_SIZE);
      WDB[ptr + BOUNDS_PTR_SURF] = (double)surf;
    }
  else if (type == UNIVERSE_TYPE_CELL)
    {
      /* Pointer to cell list */

      loc0 = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);  

      /* Loop over cells */

      while (loc0 > VALID_PTR)
        {
          /* Pointer to cell */

          cell = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
          CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);  

          /* Check type */

          if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
            {
              /* Add surfaces to list of boundaries */

              loc1 = (long)RDB[cell + CELL_PTR_SURF_LIST];
              CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

              while ((surf = (long)RDB[loc1++] ) > VALID_PTR)
                {
                  /* Find match */

                  ptr = (long)RDB[DATA_PTR_OUTER_BOUNDS];
                  while (ptr > VALID_PTR)
                    {
                      /* Compare pointer */

                      if (surf == (long)RDB[ptr + BOUNDS_PTR_SURF])
                        break;

                      /* Next */

                      ptr = NextItem(ptr);
                    }

                  /* Check if found */

                  if (ptr < VALID_PTR)
                    {
                      /* Add new */

                      ptr = NewItem(DATA_PTR_OUTER_BOUNDS, BOUNDS_BLOCK_SIZE);
                      WDB[ptr + BOUNDS_PTR_SURF] = (double)surf;
                    }
                }
              
              /* Loop over cell intersection list */

              ptr = (long)RDB[cell + CELL_PTR_SURF_INSC];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

              while (ptr > VALID_PTR)
                {
                  /* Pointer to surface */

                  surf = (long)RDB[ptr + CELL_INSC_PTR_SURF];
                  CheckPointer(FUNCTION_NAME, "(surf2)", DATA_ARRAY, surf);

                  /* Surface side */

                  side = (long)RDB[ptr + CELL_INSC_SIDE];

                  /* Get surface type */

                  type = (long)RDB[surf + SURFACE_TYPE];
                  
                  /* Check allowed types for repeated boundary conditions */
      
                  if ((side > 0) && ((type == SURF_SQC) || 
                                     (type == SURF_CUBE) || 
                                     (type == SURF_CUBOID) || 
                                     (type == SURF_HEXYC) || 
                                     (type == SURF_HEXXC) || 
                                     (type == SURF_HEXYPRISM) || 
                                     (type == SURF_HEXXPRISM) ||
                                     (type == SURF_RECT)))
                    {
                      /* Put pointer to cell structure */
                      
                      WDB[cell + CELL_PTR_BC_SURF] = (double)surf;
                      
                      /* Set OK flag */
                      
                      ok = YES;

                      /* Break loop */

                      break;
                    }

                  /* Next surface */

                  ptr = NextItem(ptr);
                }
            }

          /* Next cell */

          loc0 = NextItem(loc0);
        }
    }
  else
    Error(0, "Invalid root universe type");

  /* Check that boundary list was created */

  if ((long)RDB[DATA_PTR_OUTER_BOUNDS] < VALID_PTR)
    Error(0, "Geometry has no boundaries");

  /* Check that universe transformations are not used with root universe. */
  /* (ei toimi uuden stopatboundary.c:n kanssa) */

  if ((ptr = (long)RDB[uni + UNIVERSE_PTR_TRANS]) > VALID_PTR)
    Error(ptr, "Transformation not allowed with root universe %s",
          GetText(uni + UNIVERSE_PTR_NAME));

  /* Check repeated boundary conditions */
  
  if (((long)RDB[DATA_GEOM_BC1] == BC_BLACK) && 
      ((long)RDB[DATA_GEOM_BC2] == BC_BLACK) &&
      ((long)RDB[DATA_GEOM_BC3] == BC_BLACK))
    return;

  /* Check that surface is ok for boundary conditions */

  if (ok == NO)
    Error(0, "Invalid outer boundary for repeated boundary conditions");

  /* Check that universe symmetries are not used with root universe when */
  /* repeated boundary conditions are applied. */

  if ((ptr = (long)RDB[uni + UNIVERSE_PTR_SYM]) > VALID_PTR)
    if ((long)RDB[ptr + SYMMETRY_COORD_TRANS] == NO)
      Error(0, "Symmetry not allowed with root universe %s (repeated bc)",
            GetText(uni + UNIVERSE_PTR_NAME));

  /* Check that surface transformations are not used with outer boundaries */

  ptr = (long)RDB[DATA_PTR_OUTER_BOUNDS];
  while (ptr > VALID_PTR)
    {
      /* Check pointer */
      
      surf = (long)RDB[ptr + BOUNDS_PTR_SURF];
      CheckPointer(FUNCTION_NAME, "(surf3)", DATA_ARRAY, surf);

      /* Check transformations */

      if ((tra = (long)RDB[surf + SURFACE_PTR_TRANS]) > VALID_PTR)
        {
          /* Check poiner to name (nesteille ei aseteta?) */

          if ((long)RDB[surf + SURFACE_PTR_NAME] < VALID_PTR)
            Error(tra, "Transformations not allowed with boundary surfaces");
          else
            Error(tra, "Transformation not allowed with boundary surface %s",
                  GetText(surf + SURFACE_PTR_NAME));
        }

      /* Next */

      ptr = NextItem(ptr);
    }
}

/*****************************************************************************/
