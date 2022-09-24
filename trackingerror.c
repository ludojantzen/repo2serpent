/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : trackingerror.c                                */
/*                                                                           */
/* Created:       2014/01/16 (JLe)                                           */
/* Last modified: 2018/06/30 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Attempt to print a useful error message in case of a         */
/*              lost particle or other mysterious tracking behaviour.        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TrackingError:"

/*****************************************************************************/

void TrackingError(long etyp, double E, long mat, long type, long id)
{
  long cell, lvl0, lvl, uni, ptr, ncol, idx;
  double x, y, z, xt, yt, zt, u, v, w, totxs, majorant, dummy;
  
  /* Get coordinates at level 0 */
  
  lvl0 = (long)RDB[DATA_PTR_LVL0];
  CheckPointer(FUNCTION_NAME, "(lvl0)", DATA_ARRAY, lvl0);

  lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
  CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

  x = GetPrivateData(lvl + LVL_PRIV_X, id);
  y = GetPrivateData(lvl + LVL_PRIV_Y, id);
  z = GetPrivateData(lvl + LVL_PRIV_Z, id);
  u = GetPrivateData(lvl + LVL_PRIV_U, id);
  v = GetPrivateData(lvl + LVL_PRIV_V, id);
  w = GetPrivateData(lvl + LVL_PRIV_W, id);

  /* Reset source coordinates */

  xt = 0.0;
  yt = 0.0;
  zt = 0.0;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);
  
  /* Print standard data */

  fprintf(stdout,"\n***** %s (seed = %lu, MPI task = %d, OMP thread = %d)\n\n", 
          TimeStamp(), parent_seed, mpiid, OMP_THREAD_NUM);

  fprintf(stdout, "Tracking error:");

  if (etyp == TRACK_ERR_INF_LOOP)
    {
      /***********************************************************************/
      
      /***** Infinite loop ***************************************************/
      
      if (type == PARTICLE_TYPE_NEUTRON)
        {
          fprintf(stdout, " infinite loop (after %ld cycles)\n\n",
                  (long)RDB[DATA_NEUTRON_MAX_TRACK_LOOP]);
          fprintf(stdout, "Particle type : neutron\n\n");
        }
      else
        {
          fprintf(stdout, " infinite loop (after %ld cycles)\n\n",
                  (long)RDB[DATA_PHOTON_MAX_TRACK_LOOP]);
          fprintf(stdout, "Particle type : photon\n\n");
        }

      fprintf(stdout, "Collision = %ld\n", ncol);

      fprintf(stdout, "(x,y,z)   = %12.5E %12.5E %12.5E\n", x, y, z);
      fprintf(stdout, "(u,v,w)   = %12.5E %12.5E %12.5E\n", u, v, w);

      /* Get total cross section and majorant */

      if (mat > VALID_PTR)
        totxs = TotXS(mat, type, E, id);
      else
        totxs = 0.0;

      majorant = DTMajorant(type, E, id);

      /* Print information possibly related to infinite loop */
      
      fprintf(stdout, "E         =  %1.5E\n", E);

      if (type == PARTICLE_TYPE_NEUTRON)
        fprintf(stdout, "Emin      =  %1.5E\n", RDB[DATA_NEUTRON_EMIN]);
      else
        fprintf(stdout, "Emin      =  %1.5E\n", RDB[DATA_PHOTON_EMIN]);

      if (type == PARTICLE_TYPE_NEUTRON)
        fprintf(stdout, "Emax      =  %1.5E\n", RDB[DATA_NEUTRON_EMAX]);
      else
        fprintf(stdout, "Emax      =  %1.5E\n", RDB[DATA_PHOTON_EMAX]);

      fprintf(stdout, "material  = ");
      
      if (mat > VALID_PTR)
        fprintf(stdout, " %s\n", GetText(mat + MATERIAL_PTR_NAME));
      else
        fprintf(stdout, " void / undefined\n");
      
      fprintf(stdout, "totxs     =  %1.5E\n", totxs);
      fprintf(stdout, "majorant  =  %1.5E\n", majorant);

      /************************************************************************/
    }
  else if (etyp == TRACK_ERR_LATTICE)
    {
      /***********************************************************************/
      
      /***** Error in lattice ************************************************/
      
      fprintf(stdout, " lattice error (possibly undefined element)\n\n");

      fprintf(stdout, "(x,y,z)  = %12.5E %12.5E %12.5E\n", x, y, z);
      fprintf(stdout, "(u,v,w)  = %12.5E %12.5E %12.5E\n", u, v, w);

      /************************************************************************/
    }
  else if (etyp == TRACK_ERR_CELL_SEARCH)
    {
      /***********************************************************************/
      
      /***** Error in lattice ************************************************/
      
      fprintf(stdout, " cell search error (possibly undefined region)\n\n");

      fprintf(stdout, "(x,y,z)  = %12.5E %12.5E %12.5E\n", x, y, z);
      fprintf(stdout, "(u,v,w)  = %12.5E %12.5E %12.5E\n", u, v, w);

      /************************************************************************/
    }
  else if (etyp == TRACK_ERR_NO_MATERIAL)
    {
      /***********************************************************************/
      
      /***** Material pointer lost *******************************************/
      
      fprintf(stdout, " material pointer lost\n\n");

      fprintf(stdout, "(x,y,z)  = %12.5E %12.5E %12.5E\n", x, y, z);
      fprintf(stdout, "(u,v,w)  = %12.5E %12.5E %12.5E\n", u, v, w);

      /************************************************************************/
    }
  else if (etyp == TRACK_ERR_OUTSIDE)
    {
      /***********************************************************************/
      
      /***** Problem with outer boundary *************************************/

      /* Get pointer to root universe */

      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      fprintf(stdout, " Error handling boundary conditions.\n");
      fprintf(stdout, "                Geometry may be re-entrant or over-\n");
      fprintf(stdout, "                lapping outside root universe %s\n\n", 
              GetText(uni + UNIVERSE_PTR_NAME));

      fprintf(stdout, "(x,y,z)  = %12.5E %12.5E %12.5E\n", x, y, z);
      fprintf(stdout, "(u,v,w)  = %12.5E %12.5E %12.5E\n", u, v, w);

      /************************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid error type");

  /* Print number of lost particles */

  fprintf(stdout, "\nLost particles: %ld\n", 
          (long)RDB[DATA_UNDEF_POS_COUNT]);

  /* Check geometry errors */

  fprintf(stdout, "\nChecking geometry errors (this may fail)...\n\n");

#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Set plotter mode to check overlaps */
 
    WDB[DATA_PLOTTER_MODE] = (double)YES;
  }

  /* Add to track counter */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  AddPrivateData(ptr, 1.0, id);

  /* Get cell */

  if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
    BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, &xt, &yt, &zt, 
                       &dummy, id);

  /* Check return value */

  if (cell == GEOM_ERROR_NO_CELL)
    fprintf(stdout, "Geometry error: undefined region.\n");
  else if (cell == GEOM_ERROR_MULTIPLE_CELLS)
    fprintf(stdout, "Geometry error: overlap.\n");
  else if (cell < VALID_PTR)
    fprintf(stdout, "Geometry error: unknown.\n");
  else if (etyp == TRACK_ERR_INF_LOOP)
    {
      fprintf(stdout, "No geometry errors, but infinite loop may occur if\n");
      fprintf(stdout, "particle energy is outside the range of most cross\n");
      fprintf(stdout, "section libraries and collision probability becomes\n");
      fprintf(stdout, "very low. This can be fixed by limiting the minimum\n");
      fprintf(stdout, "and maximum energy");
      
      if (type == PARTICLE_TYPE_NEUTRON)
        fprintf(stdout, " (set egrid <tol> <Emin> <Emax>).\n");
      else
        fprintf(stdout, ".\n");
    }
  else
    fprintf(stdout, "No geometry errors.\n");

  fprintf(stdout, "\nGeometry structure at last position is printed below.\n");
  fprintf(stdout, "Some levels may not be printed correctly. Look for the\n");
  fprintf(stdout, "first level that fails: undefined region, suspicious\n");
  fprintf(stdout, "local coordinates, etc., and check the corresponding\n");
  fprintf(stdout, "geometry definition in the input.\n");

  /* Loop over levels */
  
  lvl0 = (long)RDB[DATA_PTR_LVL0];
  CheckPointer(FUNCTION_NAME, "(lvl0)", DATA_ARRAY, lvl0);

  while (lvl0 > VALID_PTR)
    {
      /* Pointer to private data */
      
      lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
      CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);
      
      /* Pointer to universe */

      if ((uni = (long)GetPrivateData(lvl + LVL_PRIV_PTR_UNIV, id)) 
          < VALID_PTR)
        {
          /* Undefined universe (should not happen) */
          
          fprintf(stdout, "\nUndefined universe at level %ld\n",
                  (long)RDB[lvl0 + LVL_NUMBER]);
          
          /* Break loop */
          
          break;
        }
      
      fprintf(stdout, "\nLevel %ld / universe %s :\n", 
              (long)RDB[lvl0 + LVL_NUMBER],
              GetText(uni + UNIVERSE_PTR_NAME));
      
      /* Check type */

      if ((type = (long)RDB[uni + UNIVERSE_TYPE]) == UNIVERSE_TYPE_NEST)
        {
          /* Get pointer to regions */

          ptr = (long)GetPrivateData(lvl + LVL_PRIV_PTR_NEST_REG, id);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          fprintf(stdout, "Nest universe / region %ld\n",
                  (long)RDB[ptr + NEST_REG_IDX]);
        }
      else if (type == UNIVERSE_TYPE_CELL)
        {
          /* Get pointer to cell */

          if ((ptr = (long)GetPrivateData(lvl + LVL_PRIV_PTR_CELL, id)) 
              > VALID_PTR)
            fprintf(stdout, "Cell universe / cell %s\n",
                    GetText(ptr + CELL_PTR_NAME));
          else
            fprintf(stdout, "Cell universe / overlap or undefined cell\n");

        }
      else if (type == UNIVERSE_TYPE_LATTICE)
        {
          /* Get pointer to lattice */

          ptr = (long)GetPrivateData(lvl + LVL_PRIV_PTR_LAT, id);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
          /* Get region index */

          if ((idx = (long)TestValuePair(ptr + LAT_COL_CELL_IDX, (double)ncol,
                                         id)) > -1)
            fprintf(stdout, "Lattice universe / region %ld\n", idx);
          else
            fprintf(stdout, "Lattice universe / undefined region\n");
        }
      else if (type == UNIVERSE_TYPE_PBED)
        {
          /* Print type */

          fprintf(stdout, "Explicit stochastic geometry universe\n");
        }
      else if (type == UNIVERSE_TYPE_UMSH)
        {
          /* Print type */

          fprintf(stdout, "Unstructured mesh universe\n");
        }
        
      fprintf(stdout, "(x,y,z)  = %12.5E %12.5E %12.5E\n", 
              GetPrivateData(lvl + LVL_PRIV_X, id),
              GetPrivateData(lvl + LVL_PRIV_Y, id),
              GetPrivateData(lvl + LVL_PRIV_Z, id));
      fprintf(stdout, "(u,v,w)  = %12.5E %12.5E %12.5E\n", 
              GetPrivateData(lvl + LVL_PRIV_U, id),
              GetPrivateData(lvl + LVL_PRIV_V, id),
              GetPrivateData(lvl + LVL_PRIV_W, id));
      
      /* Check last flag */
      
      if ((long)GetPrivateData(lvl + LVL_PRIV_LAST, id) == YES)
        break;
      
      /* Next level */
      
      lvl0 = NextItem(lvl0);
    }
  
  /* Print note about options */

  if (etyp == TRACK_ERR_INF_LOOP)
    {
      fprintf(stdout, "\nGeometry errors caused by infinite loops can be\n");  
      fprintf(stdout, "ignored by using the \"set inftrk\" input option:\n");  
      fprintf(stdout, "\nhttp://serpent.vtt.fi/mediawiki/index.php/Input_syntax_manual#set_inftrk\n");  
    }
  else
    {
      fprintf(stdout, "\nGeometry errors caused by undefined regions can be\n");  
      fprintf(stdout, "ignored by using the \"set lost\" input option:\n");  
      fprintf(stdout, "\nhttp://serpent.vtt.fi/mediawiki/index.php/Input_syntax_manual#set_lost\n");  
    }

  fprintf(stdout, "\nSimulation aborted.\n\n");
  
  /* Exit with value -1 to terminate all MPI tasks */

  exit(-1);
}

/*****************************************************************************/
