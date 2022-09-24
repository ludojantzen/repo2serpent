/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findstlsolid.c                                 */
/*                                                                           */
/* Created:       2014/03/05 (JLe)                                           */
/* Last modified: 2018/11/15 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds STL solid located in position                          */
/*                                                                           */
/* Comments: - The direction vector defines the direction in which the       */
/*             nearest boundary or known cell is searched.                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindSTLSolid:"

void STLRayTestDir(long sld, long msh, double x, double y, double z, 
                   double *u, double *v, double *w);

/*****************************************************************************/

long FindSTLSolid(long stl, double x, double y, double z, 
                  double u, double v, double w, long pre, long id)
{
  long sld, msh, ptr, pts, lst, ok, i, mode;
  
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);

  CheckValue(FUNCTION_NAME, "u", "", u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", w, -1.0, 1.0);

  /***************************************************************************/

  /***** Check for preassigned information in facet mesh *********************/

  /* Pointer to facet search mesh */

  msh = (long)RDB[stl + STL_PTR_FACET_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get pointer to search mesh cell */

  if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
    return NULLPTR;

  /* Check preassigned information */      
  
  if ((lst = (long)RDB[ptr]) > VALID_PTR)
    if ((sld = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT]) < VALID_PTR)
      return -sld;

  /* Enforce Delta-tracking if requested or switch to STL mode in */
  /* surface-tracking */

  if ((long)RDB[DATA_STL_ENFORCE_DT] == YES)
    {
      ptr = (long)RDB[DATA_DT_ENFORCE_NEXT_TRACK];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      PutPrivateData(ptr, YES, id);
    }
  else
    {
      ptr = (long)RDB[DATA_ST_USE_STL_MODE];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      PutPrivateData(ptr, YES, id);
    }

  /***************************************************************************/
  
  /****** Short list of solids in facet mesh *********************************/

  /* Check that short list exists and loop over list */

  if (lst > VALID_PTR)
    {
      /* Check */

      if (pre == NO)
        Die(FUNCTION_NAME, "WTF?");

      while (lst > VALID_PTR)
        {
          /* Pointer to solid */
          
          if ((sld = (long)RDB[lst + SEARCH_MESH_PTR_CELL_COUNT]) < VALID_PTR)
            break;
          
          /* Get mode */
          
          mode = (long)RDB[stl + STL_SEARCH_MODE];
          
          /* Pointer to fail statistics */

          pts = (long)RDB[RES_STL_RAY_TEST];
          CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

          /* Get optimal direction */
          /*
          STLRayTestDir(sld, msh, x, y, z, &u, &v, &w);
          */
          /* Resampling loop */
          
          for (i = 0; i < 100; i++)
            {
              /* Score total */

              AddBuf1D(1.0, 1.0, pts, id, 0);

              /* Check solid */
              
              if ((ok = STLRayTest(sld, msh, x, y, z, u, v, w, mode, id)) 
                  == YES)
                return sld;
              else if (ok == NO)
                break;
              else if (ok == STL_FACET_OVERLAP)
                {
                  /* Print error */
                  
                  if (mode == STL_SEARCH_MODE_FAST)
                    Error(stl, "Overlapping facets, try ray test mode 2");
                  else
                    Die(FUNCTION_NAME, "Overlap in safe mode");
                }
              else if (ok <= STL_RAY_TEST_FAIL_PARA)
                {
                  /* Score failure */
                  
                  AddBuf1D(1.0, 1.0, pts, id, -ok/1000);

                  /* Resample direction */
                  
                  IsotropicDirection(&u, &v, &w, id);
                }
              else
                Die(FUNCTION_NAME, "WTF?");
            }
                  
          /* Check for infinite loop */
          
          if (i == 100)
            {
              /* Record error */

              AddBuf1D(100.0, 1.0, pts, id, -STL_RAY_TEST_FAIL_STUCK/1000);

              /* Print warning */

              if ((long)RDB[DATA_STL_ENFORCE_DT] == NO)              
                Note(stl, "Particle stuck on boundary, try enforcing DT");
              else
                Note(stl, "Particle stuck on boundary");

              /* Put point outside */

              return NULLPTR;
            }

          /* Next in content list */
          
          lst = NextItem(lst);
        }
    }

  /* No list or search failed --> toi failure pitäis pystyä jotenkin */
  /* tunnistamaan. */
  
  /***************************************************************************/

  /***** Loop over all cells *************************************************/

  /* Pointer to solid search mesh */

  ptr = (long)RDB[stl + STL_PTR_SOLID_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, msh);

  /* Get pointer to search mesh cell */

  if ((lst = MeshPtr(ptr, x, y, z)) < VALID_PTR)
    return NULLPTR;

  /* NOTE 13.11.2018: jos ei olla short-listillä, niin pitäisi olla */
  /* ulkona. Testaa tämä ottamalla tuo testi käyttöön alempaa. */

  /*
  if (pre == YES)
    return NULLPTR;
  */

  /* Loop over list */

  lst = (long)RDB[lst];
  while (lst > VALID_PTR)
    {
      /* Pointer to solid */
      
      sld = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(sld)", DATA_ARRAY, sld);

      /* Get mode */
      
      mode = (long)RDB[stl + STL_SEARCH_MODE];
      
      /* Pointer to fail statistics */
      
      pts = (long)RDB[RES_STL_RAY_TEST];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      /* Resampling loop */
      
      for (i = 0; i < 100; i++)
        {
          /* Score total */

          AddBuf1D(1.0, 1.0, pts, id, 0);

          /* Check solid */
          
          if ((ok = STLRayTest(sld, msh, x, y, z, u, v, w, mode, id)) == YES)
            {
              /* Check plotter mode */

              if ((long)RDB[DATA_PLOTTER_MODE] == NO)
                {
                  /* Pointer to counter */
              
                  ptr = (long)RDB[lst + SEARCH_MESH_PTR_CELL_COUNT];
                  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
                  
                  /* Add counter */
                  
                  AddPrivateData(ptr, 1, id);
                }
              /*
              if (pre == YES)
                Die(FUNCTION_NAME, "WTF?");
              */
              /* Return pointer */

              return sld;
            }
          else if (ok == NO)
            break;
          else if (ok == STL_FACET_OVERLAP)
            {
              /* Print error */
              
              if (mode == STL_SEARCH_MODE_FAST)
                Error(stl, "Overlapping facets, try ray test mode 2");
              else
                Die(FUNCTION_NAME, "Overlap in safe mode");
            }
          else if (ok <= STL_RAY_TEST_FAIL_PARA)
            {
              /* Score failure */
              
              AddBuf1D(1.0, 1.0, pts, id, -ok/1000);

              /* Resample direction */
              
              IsotropicDirection(&u, &v, &w, id);
            }
          else
            Die(FUNCTION_NAME, "WTF?");
        }
      
      /* Check for infinite loop */
      
      if (i == 100)
        {
          /* Record error */
          
          AddBuf1D(100.0, 1.0, pts, id, -STL_RAY_TEST_FAIL_STUCK/1000);
          
          /* Print warning */
          
          if ((long)RDB[DATA_STL_ENFORCE_DT] == NO)              
            Note(stl, "Particle stuck on boundary, try enforcing DT");
          else
            Note(stl, "Particle stuck on boundary");

          /* Put point outside */
          
          return NULLPTR;
        }
      
      /* Next item in list */

      lst = NextItem(lst);
    }

  /***************************************************************************/

  /* Point is outside */

  return NULLPTR;
}

/*****************************************************************************/

/***** Calculate optimal direction for test ray ******************************/

void STLRayTestDir(long sld, long msh, double x, double y, double z, 
                   double *u, double *v, double *w)
{
  long loc0, loc1, ptr;
  double x0, y0, z0, xm, ym, zm, d, min;

  /* Get pointer to search mesh cell */

  if ((loc0 = MeshPtr(msh, x, y, z)) < VALID_PTR)
    return;

  /* Pointer to content */

  if ((loc0 = (long)RDB[loc0]) < VALID_PTR)
    return;
  
  /* Check pre-assigned data */

  if ((long)RDB[loc0 + SEARCH_MESH_CELL_CONTENT] < VALID_PTR)
    return;

  /* Reset minimum distance and coordinates */

  min = INFTY;
  xm = 0.0;
  ym = 0.0;
  zm = 0.0;

  /* Loop over content list */

  while (loc0 > VALID_PTR)
    {
      /* Pointer to content */

      loc1 = (long)RDB[loc0 + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      
      /* Get pointer to solid */
      
      ptr = (long)RDB[loc1 + STL_FACET_PTR_SOLID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Check pointer */

      if (ptr == sld)
        {
          /* Get points */

          ptr = (long)RDB[loc1 + STL_FACET_PTR_PT1];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          x0 = RDB[ptr + STL_POINT_X];
          y0 = RDB[ptr + STL_POINT_Y];
          z0 = RDB[ptr + STL_POINT_Z];

          ptr = (long)RDB[loc1 + STL_FACET_PTR_PT2];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          x0 = x0 + RDB[ptr + STL_POINT_X];
          y0 = y0 + RDB[ptr + STL_POINT_Y];
          z0 = z0 + RDB[ptr + STL_POINT_Z];

          ptr = (long)RDB[loc1 + STL_FACET_PTR_PT3];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          x0 = x0 + RDB[ptr + STL_POINT_X];
          y0 = y0 + RDB[ptr + STL_POINT_Y];
          z0 = z0 + RDB[ptr + STL_POINT_Z];

          x0 = x0/3.0;
          y0 = y0/3.0;
          z0 = z0/3.0;

          /* Calculate distance */

          d = (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0);
          
          /* Compare to minimum */

          if (d < min)
            {
              /* Store values */

              min = d;
              xm = x0;
              ym = y0;
              zm = z0;
            }
        }
      
      /* Next */
      
      loc0 = NextItem(loc0);
    }
  
  /* Check */

  if (min == INFTY)
    Die(FUNCTION_NAME, "WTF???");

  /* Calculate vector */

  x0 = xm - x;
  y0 = ym - y;
  z0 = zm - z;

  /* Calculate vector length */

  d = sqrt(x0*x0 + y0*y0 + z0*z0);
  CheckValue(FUNCTION_NAME, "d", "", d, ZERO, INFTY);

  /* Calculate direction cosines */

  *u = x0/d;
  *v = y0/d;
  *w = z0/d;
}

/*****************************************************************************/
