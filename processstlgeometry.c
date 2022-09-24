/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processstlgeometry.c                           */
/*                                                                           */
/* Created:       2014/03/03 (JLe)                                           */
/* Last modified: 2020/06/08 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: General processing of STL geometry data                      */
/*                                                                           */
/* Comments: - HUOM: Tota split criterionia ei välitetä createmesh.c:lle,    */
/*             vaan sinne välitetään kakkonen.                               */
/*                                                                           */
/*           - Koska nää cellit luodaan vasta nyt, niin niitä ei voi käyttää */
/*             lähteiden (ja detektorien?) kanssa. Niiden linkitys tehdään   */
/*             createuniverse.c:ssä, mitä kutsutaan ennen.                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSTLGeometry()"

/*****************************************************************************/

void ProcessSTLGeometry()
{
  long loc0, loc1, nf, ptr, pts, n, msh, sld, cell, mfile;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3, l, xmin, xmax, ymin, ymax;
  double zmin, zmax, dx, dy, dz, lims[6], mem;

  /* Check pointer */

  if ((loc0 = (long)RDB[DATA_PTR_STL0]) < VALID_PTR)
    return;

  /* Allocate memory for STL ray test counters */

  ptr = NewStat("STL_RAY_TEST", 1, 6);
  WDB[RES_STL_RAY_TEST] = (double)ptr;

  /* Loop over definitions */

  while (loc0 > VALID_PTR)
    {
      /***********************************************************************/

      /***** Create cell structures ******************************************/

      /* Check duplicate bodies */

      loc1 = (long)RDB[loc0 + STL_PTR_BODIES];
      while (loc1 > VALID_PTR)
        {
          /* Loop over remaining */

          ptr = NextItem(loc1);
          while (ptr > VALID_PTR)
            {
              /* Compare names */

              if (CompareStr(ptr + STL_BODY_PTR_BNAME,
                             loc1 + STL_BODY_PTR_BNAME))
                Error(loc0, "Duplicate definition of body \"%s\"",
                      GetText(loc1 + STL_BODY_PTR_BNAME));

              /* Next */

              ptr = NextItem(ptr);
            }

          /* Next */

          loc1 = NextItem(loc1);
        }

      /* Loop over bodies and create cells */

      loc1 = (long)RDB[loc0 + STL_PTR_BODIES];
      while (loc1 > VALID_PTR)
        {
          /* Allocate memory for structure */

          cell = NewItem(DATA_PTR_C0, CELL_BLOCK_SIZE);

          /* Put pointer */

          WDB[loc1 + STL_BODY_PTR_CELL] = (double)cell;

          /* Copy header data from STL block */

          WDB[cell + PARAM_PTR_NAME] = RDB[loc0 + PARAM_PTR_NAME];
          WDB[cell + PARAM_PTR_FNAME] = RDB[loc0 + PARAM_PTR_FNAME];
          WDB[cell + PARAM_LINE] = RDB[loc0 + PARAM_LINE];

          /* Put name */

          WDB[cell + CELL_PTR_NAME] = RDB[loc1 + STL_BODY_PTR_CNAME];

          /* Put universe and material pointers */

          WDB[cell + CELL_PTR_UNI] = RDB[loc0 + STL_PTR_NAME];
          WDB[cell + CELL_PTR_MAT] = RDB[loc1 + STL_BODY_PTR_MNAME];
          WDB[cell + CELL_PTR_FILL] = RDB[loc1 + STL_BODY_PTR_FILL];

          /* Allocate memory for collision counter */

          AllocValuePair(cell + CELL_COL_COUNT);

          /* Next body */

          loc1 = NextItem(loc1);
        }

      /* Link solids to cells */

      sld = (long)RDB[loc0 + STL_PTR_SOLIDS];
      while (sld > VALID_PTR)
        {
          /* Loop over bodies and find match */

          loc1 = (long)RDB[loc0 + STL_PTR_BODIES];
          while (loc1 > VALID_PTR)
            {
              /* Compare */

              if (CompareStr(sld + STL_SOLID_PTR_STL_NAME,
                             loc1 + STL_BODY_PTR_BNAME))
                break;

              /* Next body */

              loc1 = NextItem(loc1);
            }

          /* Check pointer */

          if (loc1 > VALID_PTR)
            {
              /* Link body and cell */

              WDB[sld + STL_SOLID_PTR_BODY] = (double)loc1;
              WDB[sld + STL_SOLID_PTR_CELL] = RDB[loc1 + STL_BODY_PTR_CELL];

              /* Add to counters */

              WDB[loc1 + STL_BODY_N_PARTS] =
                RDB[loc1 + STL_BODY_N_PARTS] + 1.0;

              WDB[loc1 + STL_BODY_N_POINTS] = RDB[loc1 + STL_BODY_N_POINTS] +
                RDB[sld + STL_SOLID_N_POINTS];

              WDB[loc1 + STL_BODY_N_FACETS] = RDB[loc1 + STL_BODY_N_FACETS] +
                RDB[sld + STL_SOLID_N_FACETS];
            }
          else
            {
              /* Print warning */

              Note(loc0, "Solid body \"%s\" not defined",
                   GetText(sld + STL_SOLID_PTR_STL_NAME));

              /* Copy pointer */

              ptr = sld;
              sld = NextItem(sld);

              /* Remove solid from list */

              RemoveItem(ptr);

              /* Cycle loop */

              continue;
            }

          /* Next solid */

          sld = NextItem(sld);
        }

      /***********************************************************************/

      /***** Create common search mesh ***************************************/

      fprintf(outp, "Processing STL geometry in universe %s...\n",
              GetText(loc0 + STL_PTR_NAME));

      /* Get mesh boundaries */

      xmin = RDB[loc0 + STL_XMIN];
      xmax = RDB[loc0 + STL_XMAX];
      ymin = RDB[loc0 + STL_YMIN];
      ymax = RDB[loc0 + STL_YMAX];
      zmin = RDB[loc0 + STL_ZMIN];
      zmax = RDB[loc0 + STL_ZMAX];

      /* Check boundaries */

      if ((xmin >= xmax) || (ymin >= ymax) || (zmin >= zmax))
        Error(loc0, "Structure is not 3D");

      /* Adjust boundaries (NOTE: geometries often have facts coinciding */
      /* with x = 0, y = 0 or z = 0. Break symmetry to avoid problems).  */

      dx = xmax - xmin;
      dy = ymax - ymin;
      dz = zmax - zmin;

      xmin = xmin - 0.01*dx;
      xmax = xmax + 0.02*dx;
      ymin = ymin - 0.03*dy;
      ymax = ymax + 0.04*dy;
      zmin = zmin - 0.05*dz;
      zmax = zmax + 0.06*dz;

      /* Put mesh variables */

      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;

      /* Get memory count */

      MemCount();

      /* Get pointer to size vector */

      ptr = (long)RDB[loc0 + STL_SEARCH_MESH_PTR_SZ];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get adaptive mesh split criterion */

      n = (long)RDB[loc0 + STL_SEARCH_MESH_ADA_SPLIT];
      CheckValue(FUNCTION_NAME, "n", "", n, 1, INFTY);

      /* Create facet mesh structure */

      msh = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR,
                       MESH_CONTENT_DATA_STL, 2, 0, 0, lims, ptr);

      /* Put pointer */

      WDB[loc0 + STL_PTR_FACET_MESH] = (double)msh;

      /* Create solid mesh structure */

      msh = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR,
                       -1, 10, 10, 10, lims, -1);
      /* Put pointer */

      WDB[loc0 + STL_PTR_SOLID_MESH] = (double)msh;

      /* Read mesh from file */

      mfile = ReadSTLMesh(loc0);

      /***********************************************************************/

      /***** Process solids **************************************************/

      /* Check number of components */

      if ((loc1 = (long)RDB[loc0 + STL_PTR_SOLIDS]) < VALID_PTR)
        Error(loc0, "Universe has no components");

      /* Loop over solids */

      while (loc1 > VALID_PTR)
        {
          /*******************************************************************/

          /***** Misc. checks ************************************************/

          /* Get number of facets */

          nf = (long)RDB[loc1 + STL_SOLID_N_FACETS];

          /* Reset solid bounding box coordinates */

          WDB[loc1 + STL_SOLID_XMIN] = INFTY;
          WDB[loc1 + STL_SOLID_XMAX] = -INFTY;
          WDB[loc1 + STL_SOLID_YMIN] = INFTY;
          WDB[loc1 + STL_SOLID_YMAX] = -INFTY;
          WDB[loc1 + STL_SOLID_ZMIN] = INFTY;
          WDB[loc1 + STL_SOLID_ZMAX] = -INFTY;

          /* Pointer to facets */

          ptr = (long)RDB[loc1 + STL_SOLID_PTR_FACETS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over facets and re-calculate normals */

          for (n = 0; n < nf; n++)
            {
              /* Reset facet bounding box coordinates */

              WDB[ptr + STL_FACET_XMIN] = INFTY;
              WDB[ptr + STL_FACET_XMAX] = -INFTY;
              WDB[ptr + STL_FACET_YMIN] = INFTY;
              WDB[ptr + STL_FACET_YMAX] = -INFTY;
              WDB[ptr + STL_FACET_ZMIN] = INFTY;
              WDB[ptr + STL_FACET_ZMAX] = -INFTY;

              /* Get points */

              pts = (long)RDB[ptr + STL_FACET_PTR_PT1];
              CheckPointer(FUNCTION_NAME, "(pts1)", DATA_ARRAY, pts);

              x1 = RDB[pts + STL_POINT_X];
              y1 = RDB[pts + STL_POINT_Y];
              z1 = RDB[pts + STL_POINT_Z];

              pts = (long)RDB[ptr + STL_FACET_PTR_PT2];
              CheckPointer(FUNCTION_NAME, "(pts2)", DATA_ARRAY, pts);

              x2 = RDB[pts + STL_POINT_X];
              y2 = RDB[pts + STL_POINT_Y];
              z2 = RDB[pts + STL_POINT_Z];

              pts = (long)RDB[ptr + STL_FACET_PTR_PT3];
              CheckPointer(FUNCTION_NAME, "(pts3)", DATA_ARRAY, pts);

              x3 = RDB[pts + STL_POINT_X];
              y3 = RDB[pts + STL_POINT_Y];
              z3 = RDB[pts + STL_POINT_Z];

              /* Check if points are on same line (jostain syystä tämä */
              /* ei välttämättä aiheuta fataalia erroria) */
              /*
              if (fabs((x2 - x1)*(x3 - x1) + (y2 - y1)*(y3 - y1) +
                       (z2 - z1)*(z3 - z1)) < ZERO)
                Note(loc0, "Degenerate facet at (%1.5E, %1.5E, %1.5E)",
                      x1, y1, z1);
              */
              /* Set bounding box for solid */

              if (x1 < RDB[loc1 + STL_SOLID_XMIN])
                WDB[loc1 + STL_SOLID_XMIN] = x1 - 0.1;
              if (x1 > RDB[loc1 + STL_SOLID_XMAX])
                WDB[loc1 + STL_SOLID_XMAX] = x1 + 0.1;
              if (y1 < RDB[loc1 + STL_SOLID_YMIN])
                WDB[loc1 + STL_SOLID_YMIN] = y1 - 0.1;
              if (y1 > RDB[loc1 + STL_SOLID_YMAX])
                WDB[loc1 + STL_SOLID_YMAX] = y1 + 0.1;
              if (z1 < RDB[loc1 + STL_SOLID_ZMIN])
                WDB[loc1 + STL_SOLID_ZMIN] = z1 - 0.1;
              if (z1 > RDB[loc1 + STL_SOLID_ZMAX])
                WDB[loc1 + STL_SOLID_ZMAX] = z1 + 0.1;
              if (x2 < RDB[loc1 + STL_SOLID_XMIN])
                WDB[loc1 + STL_SOLID_XMIN] = x2 - 0.1;
              if (x2 > RDB[loc1 + STL_SOLID_XMAX])
                WDB[loc1 + STL_SOLID_XMAX] = x2 + 0.1;
              if (y2 < RDB[loc1 + STL_SOLID_YMIN])
                WDB[loc1 + STL_SOLID_YMIN] = y2 - 0.1;
              if (y2 > RDB[loc1 + STL_SOLID_YMAX])
                WDB[loc1 + STL_SOLID_YMAX] = y2 + 0.1;
              if (z2 < RDB[loc1 + STL_SOLID_ZMIN])
                WDB[loc1 + STL_SOLID_ZMIN] = z2 - 0.1;
              if (z2 > RDB[loc1 + STL_SOLID_ZMAX])
                WDB[loc1 + STL_SOLID_ZMAX] = z2 + 0.1;
              if (x3 < RDB[loc1 + STL_SOLID_XMIN])
                WDB[loc1 + STL_SOLID_XMIN] = x3 - 0.1;
              if (x3 > RDB[loc1 + STL_SOLID_XMAX])
                WDB[loc1 + STL_SOLID_XMAX] = x3 + 0.1;
              if (y3 < RDB[loc1 + STL_SOLID_YMIN])
                WDB[loc1 + STL_SOLID_YMIN] = y3 - 0.1;
              if (y3 > RDB[loc1 + STL_SOLID_YMAX])
                WDB[loc1 + STL_SOLID_YMAX] = y3 + 0.1;
              if (z3 < RDB[loc1 + STL_SOLID_ZMIN])
                WDB[loc1 + STL_SOLID_ZMIN] = z3 - 0.1;
              if (z3 > RDB[loc1 + STL_SOLID_ZMAX])
                WDB[loc1 + STL_SOLID_ZMAX] = z3 + 0.1;

              /* Set bounding box for facet */

              if (x1 < RDB[ptr + STL_FACET_XMIN])
                WDB[ptr + STL_FACET_XMIN] = x1 - 0.1;
              if (x1 > RDB[ptr + STL_FACET_XMAX])
                WDB[ptr + STL_FACET_XMAX] = x1 + 0.1;
              if (y1 < RDB[ptr + STL_FACET_YMIN])
                WDB[ptr + STL_FACET_YMIN] = y1 - 0.1;
              if (y1 > RDB[ptr + STL_FACET_YMAX])
                WDB[ptr + STL_FACET_YMAX] = y1 + 0.1;
              if (z1 < RDB[ptr + STL_FACET_ZMIN])
                WDB[ptr + STL_FACET_ZMIN] = z1 - 0.1;
              if (z1 > RDB[ptr + STL_FACET_ZMAX])
                WDB[ptr + STL_FACET_ZMAX] = z1 + 0.1;
              if (x2 < RDB[ptr + STL_FACET_XMIN])
                WDB[ptr + STL_FACET_XMIN] = x2 - 0.1;
              if (x2 > RDB[ptr + STL_FACET_XMAX])
                WDB[ptr + STL_FACET_XMAX] = x2 + 0.1;
              if (y2 < RDB[ptr + STL_FACET_YMIN])
                WDB[ptr + STL_FACET_YMIN] = y2 - 0.1;
              if (y2 > RDB[ptr + STL_FACET_YMAX])
                WDB[ptr + STL_FACET_YMAX] = y2 + 0.1;
              if (z2 < RDB[ptr + STL_FACET_ZMIN])
                WDB[ptr + STL_FACET_ZMIN] = z2 - 0.1;
              if (z2 > RDB[ptr + STL_FACET_ZMAX])
                WDB[ptr + STL_FACET_ZMAX] = z2 + 0.1;
              if (x3 < RDB[ptr + STL_FACET_XMIN])
                WDB[ptr + STL_FACET_XMIN] = x3 - 0.1;
              if (x3 > RDB[ptr + STL_FACET_XMAX])
                WDB[ptr + STL_FACET_XMAX] = x3 + 0.1;
              if (y3 < RDB[ptr + STL_FACET_YMIN])
                WDB[ptr + STL_FACET_YMIN] = y3 - 0.1;
              if (y3 > RDB[ptr + STL_FACET_YMAX])
                WDB[ptr + STL_FACET_YMAX] = y3 + 0.1;
              if (z3 < RDB[ptr + STL_FACET_ZMIN])
                WDB[ptr + STL_FACET_ZMIN] = z3 - 0.1;
              if (z3 > RDB[ptr + STL_FACET_ZMAX])
                WDB[ptr + STL_FACET_ZMAX] = z3 + 0.1;

              /* Calculate vectors */

              x1 = x2 - x1;
              y1 = y2 - y1;
              z1 = z2 - z1;

              x2 = x3 - x2;
              y2 = y3 - y2;
              z2 = z3 - z2;

              /* Cross product */

              x3 =  (y1*z2 - z1*y2);
              y3 = -(x1*z2 - z1*x2);
              z3 =  (x1*y2 - y1*x2);

              /* Normalize */

              if ((l = x3*x3 + y3*y3 + z3*z3) > 0.0)
                {
                  l = sqrt(l);
                  x3 = x3/l;
                  y3 = y3/l;
                  z3 = z3/l;
                }
              else if (RDB[loc0 + STL_MERGE_RAD] > 0.0)
                Error(loc0,
                      "Degenerate facet at (%1.5E, %1.5E, %1.5E), try decreasing the merge radius", RDB[pts + STL_POINT_X], RDB[pts + STL_POINT_Y], RDB[pts + STL_POINT_Z]);
              else if ((long)RDB[loc0 + STL_SEARCH_MODE] ==
                       STL_SEARCH_MODE_SAFE)
                Note(loc0, "Degenerate facet at (%1.5E, %1.5E, %1.5E)",
                     RDB[pts + STL_POINT_X], RDB[pts + STL_POINT_Y],
                     RDB[pts + STL_POINT_Z]);
              else
                Error(loc0,
                      "Degenerate facet, try safe mode to ignore");

              /* Check direction of normal vector */
              /*
              if (fabs(RDB[ptr + STL_FACET_NORM_U]*x3 +
                       RDB[ptr + STL_FACET_NORM_V]*y3 +
                       RDB[ptr + STL_FACET_NORM_W]*z3 - 1.0) > 1E-6)
                Note(loc0,
                     "Inconsistent normal vector at (%1.5E, %1.5E, %1.5E)",
                     RDB[pts + STL_POINT_X], RDB[pts + STL_POINT_Y],
                     RDB[pts + STL_POINT_Z]);
              */
              /* Put normal vector */

              WDB[ptr + STL_FACET_NORM_U] = x3;
              WDB[ptr + STL_FACET_NORM_V] = y3;
              WDB[ptr + STL_FACET_NORM_W] = z3;

              /* Pointer to next facet */

              ptr = ptr + STL_FACET_BLOCK_SIZE;
            }

          /* Pointer to solid search mesh */

          msh = (long)RDB[loc0 + STL_PTR_SOLID_MESH];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

          /* Get bounding box */

          xmin = RDB[loc1 + STL_SOLID_XMIN];
          xmax = RDB[loc1 + STL_SOLID_XMAX];
          ymin = RDB[loc1 + STL_SOLID_YMIN];
          ymax = RDB[loc1 + STL_SOLID_YMAX];
          zmin = RDB[loc1 + STL_SOLID_ZMIN];
          zmax = RDB[loc1 + STL_SOLID_ZMAX];

          /* Adjust boundaries */

          dx = xmax - xmin;
          dy = ymax - ymin;
          dz = zmax - zmin;

          xmin = xmin - 0.01*dx;
          xmax = xmax + 0.01*dx;
          ymin = ymin - 0.01*dy;
          ymax = ymax + 0.01*dy;
          zmin = zmin - 0.01*dz;
          zmax = zmax + 0.01*dz;

          /* Add to search mesh */

          if (mfile == NO)
            AddSearchMesh(msh, loc1, xmin, xmax, ymin, ymax, zmin, zmax);

          /*******************************************************************/

          /***** Put facets to search mesh ***********************************/

          /* Pointer to facet search mesh */

          msh = (long)RDB[loc0 + STL_PTR_FACET_MESH];
          CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

          /* Pointer to facets */

          ptr = (long)RDB[loc1 + STL_SOLID_PTR_FACETS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over facets */

          for (n = 0; n < nf; n++)
            {
              /* Reset boundaries */

              xmin = INFTY;
              xmax = -INFTY;
              ymin = INFTY;
              ymax = -INFTY;
              zmin = INFTY;
              zmax = -INFTY;

              /* Get points */

              pts = (long)RDB[ptr + STL_FACET_PTR_PT1];
              CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

              x1 = RDB[pts + STL_POINT_X];
              y1 = RDB[pts + STL_POINT_Y];
              z1 = RDB[pts + STL_POINT_Z];

              pts = (long)RDB[ptr + STL_FACET_PTR_PT2];
              CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

              x2 = RDB[pts + STL_POINT_X];
              y2 = RDB[pts + STL_POINT_Y];
              z2 = RDB[pts + STL_POINT_Z];

              pts = (long)RDB[ptr + STL_FACET_PTR_PT3];
              CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

              x3 = RDB[pts + STL_POINT_X];
              y3 = RDB[pts + STL_POINT_Y];
              z3 = RDB[pts + STL_POINT_Z];

              /* Find minimum and maximum */

              if (x1 < xmin)
                xmin = x1;
              if (x1 > xmax)
                xmax = x1;

              if (y1 < ymin)
                ymin = y1;
              if (y1 > ymax)
                ymax = y1;

              if (z1 < zmin)
                zmin = z1;
              if (z1 > zmax)
                zmax = z1;

              if (x2 < xmin)
                xmin = x2;
              if (x2 > xmax)
                xmax = x2;

              if (y2 < ymin)
                ymin = y2;
              if (y2 > ymax)
                ymax = y2;

              if (z2 < zmin)
                zmin = z2;
              if (z2 > zmax)
                zmax = z2;

              if (x3 < xmin)
                xmin = x3;
              if (x3 > xmax)
                xmax = x3;

              if (y3 < ymin)
                ymin = y3;
              if (y3 > ymax)
                ymax = y3;

              if (z3 < zmin)
                zmin = z3;
              if (z3 > zmax)
                zmax = z3;

              /* Adjust boundaries */

              dx = xmax - xmin;
              dy = ymax - ymin;
              dz = zmax - zmin;

              xmin = xmin - 0.01*dx;
              xmax = xmax + 0.01*dx;
              ymin = ymin - 0.01*dy;
              ymax = ymax + 0.01*dy;
              zmin = zmin - 0.01*dz;
              zmax = zmax + 0.01*dz;

              /* Add to search mesh */

              if (mfile == NO)
                AddSearchMesh(msh, ptr, xmin, xmax, ymin, ymax, zmin, zmax);

                /* Pointer to next facet */

              ptr = ptr + STL_FACET_BLOCK_SIZE;
            }

          /*******************************************************************/

          /* Next solid */

          loc1 = NextItem(loc1);
        }

      /* Loop over bodies and print */

      loc1 = (long)RDB[loc0 + STL_PTR_BODIES];
      while (loc1 > VALID_PTR)
        {
          fprintf(outp, "\nSTL body \"%s\" :\n\n - Cell: \"%s\"\n",
                  GetText(loc1 + STL_BODY_PTR_BNAME),
                  GetText(loc1 + STL_BODY_PTR_CNAME));

          if ((long)RDB[loc1 + STL_BODY_PTR_MNAME] > VALID_PTR)
            fprintf(outp, " - Material: \"%s\"\n",
                    GetText(loc1 + STL_BODY_PTR_MNAME));
          else if ((ptr = (long)RDB[loc1 + STL_BODY_PTR_FILL]) > VALID_PTR)
            fprintf(outp, " - Filled with universe: \"%s\"\n",
                    GetText(ptr + UNIVERSE_PTR_NAME));
          else
            Die(FUNCTION_NAME, "Material and fill pointers null");

          fprintf(outp, " - %ld components\n",
                  (long)RDB[loc1 + STL_BODY_N_PARTS]);
          fprintf(outp, " - %ld triangular facets\n",
                  (long)RDB[loc1 + STL_BODY_N_FACETS]);
          fprintf(outp, " - %ld (unique) points\n",
                  (long)RDB[loc1 + STL_BODY_N_POINTS]);

          /* Next body */

          loc1 = NextItem(loc1);
        }

      /* Pointer to search mesh */

      msh = (long)RDB[loc0 + STL_PTR_FACET_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Reconfigure pointers for faster solid search */

      SetSTLMeshPointers(msh);

      /* Allow access to private data blocks (tarvitaan nyt ton */
      /* viritelmän takia) */

      WDB[DATA_PRIVA_MEM_READY] = (double)YES;

      /* Test solids separately */

      if (1 == 2)
        TestSTLSolids(loc0, 100000, 10);

      /* Fill search mesh with pointers and stuff */

      if (mfile == NO)
        {
          fprintf(outp, "\nFilling empty search mesh cells...\n");

          FillSTLMesh(loc0, msh, xmin, ymin, zmin);

          fprintf(outp, "OK.\n");
        }

      /* Get mesh boundaries */

      xmin = RDB[msh + MESH_MIN0];
      xmax = RDB[msh + MESH_MAX0];
      ymin = RDB[msh + MESH_MIN1];
      ymax = RDB[msh + MESH_MAX1];
      zmin = RDB[msh + MESH_MIN2];
      zmax = RDB[msh + MESH_MAX2];

      /* Allow allocation of private data blocks (tarvitaan nyt ton */
      /* viritelmän takia) */

      WDB[DATA_PRIVA_MEM_READY] = (double)NO;

      /* Print search mesh stuff */

      fprintf(outp, "\nAdaptive search mesh:\n\n");

      if (mfile == YES)
        fprintf(outp, " - Read from pre-generated file\n");

      fprintf(outp, " - Dimensions: x = [%1.1f, %1.1f]; ",
              RDB[loc0 + STL_XMIN], RDB[loc0 + STL_XMAX]);
      fprintf(outp, "y = [%1.1f, %1.1f]; ",
              RDB[loc0 + STL_YMIN], RDB[loc0 + STL_YMAX]);
      fprintf(outp, "z = [%1.1f, %1.1f]\n",
              RDB[loc0 + STL_ZMIN], RDB[loc0 + STL_ZMAX]);

      fprintf(outp, " - Depth:");

      ptr = (long)RDB[loc0 + STL_SEARCH_MESH_PTR_SZ];
      while ((n = (long)RDB[ptr++]) > 0)
        fprintf(outp, " %ld", n);

      fprintf(outp, "\n - %ld mesh cells \n",
              (long)RDB[loc0 + STL_SEARCH_MESH_CELLS]);

      fprintf(outp, " - %1.1f%% of volume consists of pre-assigned data\n",
              100.0*RDB[loc0 + STL_SEARCH_MESH_V]/
              ((xmax - xmin)*(ymax - ymin)*(zmax - zmin)));

      mem = (double)MemCount();

      if (mem < MEGA)
        fprintf(outp, " - %1.2f kb of memory allocated for search mesh\n\n",
                mem/KILO);
      else if (mem < GIGA)
        fprintf(outp, " - %1.2f Mb of memory allocated for search mesh\n\n",
                mem/MEGA);
      else
        fprintf(outp, " - %1.2f Gb of memory allocated for search mesh\n\n",
                    mem/GIGA);

      /* Write mesh to file */

      if (mfile == NO)
        WriteSTLMesh(loc0);

      /* Next */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
