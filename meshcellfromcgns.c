/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshcellfromcgns.c                             */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2017/11/29 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Takes a IFC_HEX_MSH cgns cell and creates a point list       */
/*              out of it                                                    */
/*                                                                           */
/* Comments:    -For fixhexmesh.c                                            */
/*                                                                           */
/*                            3---5                                          */
/*             7-----6        |\ /|                                          */
/*            /|    /|        | 4 |           4            3                 */
/*           4-----5 |        | | |          /|\          /|\                */
/*           | 3---|-2        0-|-2        3-----2       0-|-2               */
/*           |/    |/          \|/        /     /         \|/                */
/*           0-----1            1        0-----1           1                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshCellFromCGNS:"

/*****************************************************************************/

void  MeshCellFromCGNS(long ifc, long cgns, long *V,
                       long *hexfaces, long (*initFaces)[4], long celltype)
{
  long loc2, loc3;
  long ptr, ptr1, pointlist, surflist, surf;
  long ownrlist, nbrlist;
  long np, i, j, n, k, maxf, oppose;
  long side, id, p0, p1, topfound;
  long pointfour, bot[4], top[4], left[4];

  /* Get pointer to parents pointlist */

  pointlist = (long)RDB[ifc + IFC_PTR_POINT_LIST_PARENTS];;
  CheckPointer(FUNCTION_NAME, "(pointlist)", DATA_ARRAY, pointlist);

  /* Get pointer to parents surfacelist */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to parents owner list */

  ownrlist = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(ownrlist)", DATA_ARRAY, ownrlist);

  /* Get pointer to parents neighbour list */

  nbrlist  = (long)RDB[ifc + IFC_PTR_NBR_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(nbrlist)", DATA_ARRAY, nbrlist);

  /* Get pointer to face list */

  ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, ptr);

  /* Get pointer to side list */

  ptr1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(PTR1)", DATA_ARRAY, ptr1);

  /* Get number of faces */

  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:
      maxf = 4;

      break;

    case MESH_CELL_TYPE_PYRAMID:
      maxf = 5;

      break;

    case MESH_CELL_TYPE_PRISM:
      maxf = 5;

      break;

    case MESH_CELL_TYPE_HEX:
      maxf = 6;

      break;
    default:

      maxf = -1;

      Die(FUNCTION_NAME, "Invalid cell type %ld", celltype);

    }

  /* Make list of faces */

  for (j = 0; j < maxf; j++)
    {
      hexfaces[j] = (long)RDB[ptr + j];
    }

  /* Reset bottom face */

  bot[0] = -1;
  bot[1] = -1;
  bot[2] = -1;
  bot[3] = -1;

  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:

      /* Get one face (will be bottom) */

      for (j = 0; j < 4; j++)
        {

          /* Get index of first face */

          n = (long)RDB[ptr + j];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get side of face */

          side = (long)RDB[ptr1 + j];

          /* Get pointer to surface parameters (point list) */

          loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

          /* Loop over points */

          for (i = 0; i < 3; i++)
            {
              /* Get point */
              loc3 = (long)RDB[loc2 + i];

              /* Calculate point id */

              id = (long)((loc3 - pointlist)/3);

              /* Store first face */
              /* This will be the bottom */

              if ((j == 0) && (side == -1))
                {
                  /* Owner */

                  /* Store index (change direction) */

                  bot[2-i] = id;

                }
              else if ((j == 0) && (side == 1))
                {
                  /* Store index (do not change) */

                  bot[i] = id;

                }
              else
                {
                  /* Store point that is not in bottom */

                  if (PolyPInF(id, bot, 3) == NO)
                    bot[3] = id;

                }

              initFaces[j][i] = id;

            }

        }

      /* Create the hexahedron */

      for (i = 0; i < 4; i++)
        V[i] = bot[i];

      break;

    case MESH_CELL_TYPE_PYRAMID:
      /* Reset pointfour */

      pointfour = -1;

      /* Get only 4 point face (will be bottom) */

      for (j = 0; j < 5; j++)
        {

          /* Get index of face */

          n = (long)RDB[ptr + j];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get number of points on the face */
          /* and continue if not 4 */

          if ((np = (long)RDB[surf + UMSH_SURF_N_POINTS]) != 4)
            continue;

          /* Get side of face */

          side = (long)RDB[ptr1 + j];

          /* Get pointer to surface parameters (point list) */

          loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

          /* Loop over points */

          for (i = 0; i < np; i++)
            {
              /* Get point */
              loc3 = (long)RDB[loc2 + i];

              /* Calculate point id */

              id = (long)((loc3 - pointlist)/3);

              /* Store face  */
              /* This will be the bottom */

              if (side == -1)
                {
                  /* Owner */

                  /* Store index (change direction) */

                  bot[3-i] = id;

                }
              else if (side == 1)
                {
                  /* Store index (do not change) */

                  bot[i] = id;

                }
              else
                Die(FUNCTION_NAME, "WTF?");

              initFaces[j][i] = id;

            }

          break;
        }

      /* Loop over other faces to find another side */

      for (i = 0; i < 5; i++)
        {
          /* Do not handle the bottom face again */

          if (i == j)
            continue;

          /* Get index of face */

          n = (long)RDB[ptr + i];

          /* Get side of face */

          side = (long)RDB[ptr1 + i];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get number of points on the face */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Get pointer to surface parameters (point list) */

          loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

          /* Loop over points */

          for (k = 0; k < np; k++)
            {
              /* Get point */
              loc3 = (long)RDB[loc2 + k];

              id = (long)((loc3 - pointlist)/3);

              /* Find the point that does not exist in bot face */

              if (PolyPInF(id, bot, 4) == NO)
                pointfour = id;

              /* Store face points to initFaces */

              initFaces[i][k] = id;

            }

        }

      /* Create the hexahedron */

      for (i = 0; i < 4; i++)
        {
          V[i] = bot[i];
        }

      V[4] = pointfour;

      break;

    case MESH_CELL_TYPE_PRISM:

      /* Get first triangular face (will be bottom) */

      for (j = 0; j < 5; j++)
        {

          /* Get index of face */

          n = (long)RDB[ptr + j];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get number of points on the face */
          /* and continue if not 3 */

          if ((np = (long)RDB[surf + UMSH_SURF_N_POINTS]) != 3)
            continue;

          /* Get side of face */

          side = (long)RDB[ptr1 + j];

          /* Get pointer to surface parameters (point list) */

          loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

          /* Loop over points */

          for (i = 0; i < np; i++)
            {
              /* Get point */
              loc3 = (long)RDB[loc2 + i];

              /* Calculate point id */

              id = (long)((loc3 - pointlist)/3);

              /* Store face  */
              /* This will be the bottom */

              if (side == -1)
                {
                  /* Owner */

                  /* Store index (change direction) */

                  bot[2-i] = id;

                }
              else if (side == 1)
                {
                  /* Store index (do not change) */

                  bot[i] = id;

                }
              else
                Die(FUNCTION_NAME, "WTF?");

              initFaces[j][i] = id;

            }

          break;
        }

      /* Loop over other faces to find the opposing side */

      for (i = 0; i < 5; i++)
        {
          /* Do not handle the bottom face again */

          if (i == j)
            continue;

          /* Get index of face */

          n = (long)RDB[ptr + i];

          /* Get side of face */

          side = (long)RDB[ptr1 + i];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get number of points on the face */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Get pointer to surface parameters (point list) */

          loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

          oppose = YES;

          /* Loop over points */

          for (k = 0; k < np; k++)
            {
              /* Get point */
              loc3 = (long)RDB[loc2 + k];

              id = (long)((loc3 - pointlist)/3);

              /* Check if this is opposing side to the one we stored */
              /* They do not share any points */
              if (np == 3)
                {
                  if (PolyPInF(id, bot, 3) == YES)
                    oppose = NO;
                }
              else
                oppose = NO;

              /* Store second face (Top) */

              if (side == -1)
                {
                  /* Owner */

                  /* Store index (Do not rotate) */

                  if (oppose == YES)
                    top[k] = id;
                  else if (np == 4)
                    left[k] = id;
                  else
                    Die(FUNCTION_NAME, "Something is wrong with this prism");

                }
              else if (side == 1)
                {
                  /* Store index (Rotate) */

                  if (oppose == YES)
                    top[2-k] = id;
                  else if (np == 4)
                    left[3-k] = id;
                  else
                    Die(FUNCTION_NAME, "Something is wrong with this prism");
                }
              else
                Die(FUNCTION_NAME, "WTF?");

              /* Store face points to initFaces */

              initFaces[i][k] = id;

            }

        }

      /* Find pair of points from side that are connected in top and bottom */

      for (i = 0; i < 3; i++)
        if (((PolyPInF(left[i], bot, 3)) && (PolyPInF(left[i+1], top, 3))) ||
            ((PolyPInF(left[i], top, 3)) && (PolyPInF(left[i+1], bot, 3))))
          break;

      if (i == 3)
        Die(FUNCTION_NAME, "Bad side");

      if (PolyPInF(left[i], bot, 3))
        {
          p0 = left[i];
          p1 = left[i+1];
        }
      else
        {
          p0 = left[i+1];
          p1 = left[i];
        }

      /* Now p0 is in bot and p1 in top */
      /* They connect these two opposite faces */
      /* Let's rotate top until p0 and p1 match */

      for (i = 0; i < 3; i++)
        if (bot[i] == p0)
          break;

      for (j = 0; j < 3; j++)
        if (top[j] == p1)
          break;

      /* Top will have to be rotated i-j steps */

      HexRotateFace(top, i-j, 3);

      /* Create the hexahedron */

      for (i = 0; i < 3; i++)
        {
          V[i] = bot[i];
          V[3+i] = top[i];
        }

      break;

    case MESH_CELL_TYPE_HEX:

      /* Get index of first face */

      n = (long)RDB[ptr + 0];

      /* Get pointer to surface */

      surf = (long)RDB[surflist + n];
      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

      /* Get side of first face */

      side = (long)RDB[ptr1 + 0];

      /* Get number of points on the face */

      np = (long)RDB[surf + UMSH_SURF_N_POINTS];

      /* Get pointer to surface parameters (point list) */

      loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

      /* Loop over points */

      for (i = 0; i < np; i++)
        {
          /* Get point */
          loc3 = (long)RDB[loc2 + i];

          /* Calculate point id */

          id = (long)((loc3 - pointlist)/3);

          /* Store first face        */
          /* This will be the bottom */

          if (side == -1)
            {
              /* Owner */

              /* Store index (change direction) */

              bot[3-i] = id;

            }
          else if (side == 1)
            {
              /* Store index (do not change) */

              bot[i] = id;

            }
          else
            Die(FUNCTION_NAME, "WTF?");

          initFaces[0][i] = id;

        }

      /* We still need to find an opposing face (top) */
      /* And one of the others (side/left)            */
      /* Reset found flags                            */

      topfound = 0;

      /* Loop over other faces to find the opposing side */

      for (i = 1; i < 6; i++)
        {
          /* Get index of face */

          n = (long)RDB[ptr + i];

          /* Get side of face */

          side = (long)RDB[ptr1 + i];

          /* Get pointer to surface */

          surf = (long)RDB[surflist + n];
          CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

          /* Get number of points on the face */

          np = (long)RDB[surf + UMSH_SURF_N_POINTS];

          /* Get pointer to surface parameters (point list) */

          loc2 = (long)RDB[surf + UMSH_SURF_PTR_POINTS];

          oppose = YES;

          /* Loop over points */

          for (j = 0; j < np; j++)
            {
              /* Get point */
              loc3 = (long)RDB[loc2 + j];

              id = (long)((loc3 - pointlist)/3);

              /* Check if this is opposing side to the one we stored */
              /* They do not share any points */

              if (PolyPInF(id, bot, 4) == YES)
                oppose = NO;

              /* Store second face (Top) */

              if (side == -1)
                {
                  /* Owner */

                  /* Store index (Do not rotate) */

                  if (topfound == 0)
                    top[j] = id;
                  else
                    left[j] = id;

                }
              else if (side == 1)
                {
                  /* Store index (Rotate) */

                  if (topfound == 0)
                    top[3-j] = id;
                  else
                    left[3-j] = id;

                }
              else
                Die(FUNCTION_NAME, "WTF?");

              /* Store face points to initFaces */

              initFaces[i][j] = id;

            }

          /* Check if this was opposing face */

          if (oppose == YES)
            {
              topfound = 1;
            }
          else
            {
              /* One of the sides */
              /* Store sides */

              if (topfound == 0)
                for (j = 0; j < 4; j++)
                  left[j] = top[j];

            }

        }

      /* Find pair of points from side that are connected in top and bottom */

      for (i = 0; i < 3; i++)
        if (((PolyPInF(left[i], bot, 4)) && (PolyPInF(left[i+1], top, 4))) ||
            ((PolyPInF(left[i], top, 4)) && (PolyPInF(left[i+1], bot, 4))))
          break;

      if (i == 3)
        Die(FUNCTION_NAME, "Bad side");

      if (PolyPInF(left[i], bot, 4))
        {
          p0 = left[i];
          p1 = left[i+1];
        }
      else
        {
          p0 = left[i+1];
          p1 = left[i];
        }

      /* Now p0 is in bot and p1 in top */
      /* They connect these two opposite faces */
      /* Let's rotate top until p0 and p1 match */

      for (i = 0; i < 4; i++)
        if (bot[i] == p0)
          break;

      for (j = 0; j < 4; j++)
        if (top[j] == p1)
          break;

      /* Top will have to be rotated i-j steps */

      HexRotateFace(top, i-j, 4);

      /* Create the hexahedron */

      for (i = 0; i < 4; i++)
        {
          V[i] = bot[i];
          V[4+i] = top[i];
        }


      break;
    default:

      maxf = -1;

      Die(FUNCTION_NAME, "Invalid cell type %ld", celltype);

    }

 }

/*****************************************************************************/
