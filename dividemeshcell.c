/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dividemeshcell.c                               */
/*                                                                           */
/* Created:       2015/08/31 (VVa)                                           */
/* Last modified: 2017/11/29 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Splits a mesh cell into tetrahedrons                         */
/*                                                                           */
/* Comments:   -Tetrahedrons are just copied to the list of                  */
/*              child tetrahedrons                                           */
/*             -Based on "J. Dompierre, P. Labbe, M. Vallet and R. Camarero, */
/*              How to Subdivide Tets, Prisms, and Hexahedra                 */
/*              into Tetrahedra, Proceedings,                                */
/*              8th International Meshing Roundtable,                        */
/*              South Lake Tahoe, CA, U.S.A., pp.195-204, October 1999       */
/*                                                                           */
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

#define FUNCTION_NAME "DivideMeshCell:"

/*****************************************************************************/

void DivideMeshCell(long ifc, long cgns, long hex[8], long (*diags)[8],
                    long hexfaces[6], long **newnbrs, long **newnbrfaces,
                    long *totidx, long celltype)
{
  long nd, nc, i, j, k, pts[6][4], directions[6][4], direction;
  long p0, p1, p2, p3;
  long pointlist, tetlist, tet;
  long io, otet, fidx;
  double u, u0, u1, v, v0, v1, w, w0, w1;
  long *face[4];
  long face1[4] = {0,2,1,3};
  long face2[4] = {1,2,3,0};
  long face3[4] = {0,1,3,2};
  long face4[4] = {0,3,2,1};

  /* Store point list for easy access */

  face[0] = face1;
  face[1] = face2;
  face[2] = face3;
  face[3] = face4;

  /* Reset number of diagonals */

  nd = 0;

  /* Get pointer to parents pointlist */

  pointlist = (long)RDB[ifc + IFC_PTR_POINT_LIST_PARENTS];
  CheckPointer(FUNCTION_NAME, "(pointlist)", DATA_ARRAY, pointlist);

  /* Get pointer to tet list */

  tetlist = (long)RDB[ifc + IFC_PTR_TET_LIST];
  CheckPointer(FUNCTION_NAME, "(tetlist)", DATA_ARRAY, tetlist);

  /* Get number of subcells based on celltype */

  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:
      nc = 1;

      break;
    case MESH_CELL_TYPE_PYRAMID:
      nc = 2;

      break;
    case MESH_CELL_TYPE_PRISM:
      nc = 3;

      break;
    case MESH_CELL_TYPE_HEX:

      /* Count number of diagonals going through 6 */

      nd = 0;

      for (i = 0; i < 8; i++)
        if (diags[6][i] == 1)
          nd++;

      /* Get number of subcells based on number of diagonals */

      if (nd == 0)
        nc = 5;
      else
        nc = 6;

      break;
    default:
      nc = 0;
      Die(FUNCTION_NAME, "Unknown mesh cell type %ld\n", celltype);

    }


  /* Create the new cells */

  for (i = 0; i < nc; i++)
    {
      /* Create new tet cell */

      tet = ReallocMem(DATA_ARRAY, TET_BLOCK_SIZE);
      CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

      /* Put pointer to parent */

      WDB[tet + TET_PTR_PARENT] = (double)cgns;

      /* Put pointer to tet-list */

      WDB[tetlist + *totidx] = (double)tet;

      /* Increment total tet index */

      *totidx = *totidx + 1;
    }

  /* Get points and directions based on cell type */

  switch (celltype)
    {
    case MESH_CELL_TYPE_TET:
      /* Get points and directions */

      pts[0][0] = hex[0];
      pts[0][1] = hex[1];
      pts[0][2] = hex[2];
      pts[0][3] = hex[3];

      directions[0][0] = MESH_CELL_FACE_BOTTOM;
      directions[0][1] = MESH_CELL_FACE_RIGHT;
      directions[0][2] = MESH_CELL_FACE_LEFT;
      directions[0][3] = MESH_CELL_FACE_BACK;

      break;

    case MESH_CELL_TYPE_PYRAMID:
      if ((hex[0] < hex[2] ? hex[0] : hex[2])
          < (hex[1] < hex[3] ? hex[1] : hex[3]))
        {
          /* Diagonal through 0-2 */

          pts[0][0] = hex[0];
          pts[0][1] = hex[1];
          pts[0][2] = hex[2];
          pts[0][3] = hex[4];

          pts[1][0] = hex[0];
          pts[1][1] = hex[2];
          pts[1][2] = hex[3];
          pts[1][3] = hex[4];

          directions[0][0] = MESH_CELL_FACE_BOTTOM;
          directions[0][1] = MESH_CELL_FACE_RIGHT;
          directions[0][2] = MESH_CELL_FACE_FRONT;
          directions[0][3] = -2;

          directions[1][0] = MESH_CELL_FACE_BOTTOM;
          directions[1][1] = MESH_CELL_FACE_BACK;
          directions[1][2] = -1;
          directions[1][3] = MESH_CELL_FACE_LEFT;

        }
      else
        {

          /* Diagonal through 1-3 */

          pts[0][0] = hex[1];
          pts[0][1] = hex[2];
          pts[0][2] = hex[3];
          pts[0][3] = hex[4];

          pts[1][0] = hex[1];
          pts[1][1] = hex[3];
          pts[1][2] = hex[0];
          pts[1][3] = hex[4];

          directions[0][0] = MESH_CELL_FACE_BOTTOM;
          directions[0][1] = MESH_CELL_FACE_BACK;
          directions[0][2] = MESH_CELL_FACE_RIGHT;
          directions[0][3] = -2;

          directions[1][0] = MESH_CELL_FACE_BOTTOM;
          directions[1][1] = MESH_CELL_FACE_LEFT;
          directions[1][2] = -1;
          directions[1][3] = MESH_CELL_FACE_FRONT;

        }

      break;

    case MESH_CELL_TYPE_PRISM:
      if ((hex[1] < hex[5] ? hex[1] : hex[5])
          < (hex[2] < hex[4] ? hex[2] : hex[4]))
        {
          /* Diagonal through 1-5 */

          pts[0][0] = hex[0];
          pts[0][1] = hex[1];
          pts[0][2] = hex[2];
          pts[0][3] = hex[5];

          pts[1][0] = hex[0];
          pts[1][1] = hex[1];
          pts[1][2] = hex[5];
          pts[1][3] = hex[4];

          pts[2][0] = hex[0];
          pts[2][1] = hex[4];
          pts[2][2] = hex[5];
          pts[2][3] = hex[3];

          directions[0][0] = MESH_CELL_FACE_BOTTOM;
          directions[0][1] = MESH_CELL_FACE_RIGHT;
          directions[0][2] = -2;
          directions[0][3] = MESH_CELL_FACE_BACK;

          directions[1][0] = -1;
          directions[1][1] = MESH_CELL_FACE_RIGHT;
          directions[1][2] = MESH_CELL_FACE_LEFT;
          directions[1][3] = -3;

          directions[2][0] = -2;
          directions[2][1] = MESH_CELL_FACE_TOP;
          directions[2][2] = MESH_CELL_FACE_LEFT;
          directions[2][3] = MESH_CELL_FACE_BACK;

        }
      else
        {
          pts[0][0] = hex[0];
          pts[0][1] = hex[1];
          pts[0][2] = hex[2];
          pts[0][3] = hex[4];

          pts[1][0] = hex[0];
          pts[1][1] = hex[4];
          pts[1][2] = hex[2];
          pts[1][3] = hex[5];

          pts[2][0] = hex[0];
          pts[2][1] = hex[4];
          pts[2][2] = hex[5];
          pts[2][3] = hex[3];

          directions[0][0] = MESH_CELL_FACE_BOTTOM;
          directions[0][1] = MESH_CELL_FACE_RIGHT;
          directions[0][2] = MESH_CELL_FACE_LEFT;
          directions[0][3] = -2;

          directions[1][0] = -1;
          directions[1][1] = MESH_CELL_FACE_RIGHT;
          directions[1][2] = -3;
          directions[1][3] = MESH_CELL_FACE_BACK;

          directions[2][0] = -2;
          directions[2][1] = MESH_CELL_FACE_TOP;
          directions[2][2] = MESH_CELL_FACE_LEFT;
          directions[2][3] = MESH_CELL_FACE_BACK;

        }

      break;

    case MESH_CELL_TYPE_HEX:

      if (nd == 0)
        {
          pts[0][0] = hex[0];
          pts[0][1] = hex[1];
          pts[0][2] = hex[2];
          pts[0][3] = hex[5];

          pts[1][0] = hex[0];
          pts[1][1] = hex[2];
          pts[1][2] = hex[7];
          pts[1][3] = hex[5];

          pts[2][0] = hex[0];
          pts[2][1] = hex[2];
          pts[2][2] = hex[3];
          pts[2][3] = hex[7];

          pts[3][0] = hex[0];
          pts[3][1] = hex[5];
          pts[3][2] = hex[7];
          pts[3][3] = hex[4];

          pts[4][0] = hex[2];
          pts[4][1] = hex[7];
          pts[4][2] = hex[5];
          pts[4][3] = hex[6];

          pts[5][0] = -1;
          pts[5][1] = -1;
          pts[5][2] = -1;
          pts[5][3] = -1;

          directions[0][0] = MESH_CELL_FACE_BOTTOM;
          directions[0][1] = MESH_CELL_FACE_RIGHT;
          directions[0][2] = MESH_CELL_FACE_FRONT;
          directions[0][3] = -2;

          directions[1][0] = -3;
          directions[1][1] = -5;
          directions[1][2] = -1;
          directions[1][3] = -4;

          directions[2][0] = MESH_CELL_FACE_BOTTOM;
          directions[2][1] = MESH_CELL_FACE_BACK;
          directions[2][2] = -2;
          directions[2][3] = MESH_CELL_FACE_LEFT;

          directions[3][0] = -2;
          directions[3][1] = MESH_CELL_FACE_TOP;
          directions[3][2] = MESH_CELL_FACE_FRONT;
          directions[3][3] = MESH_CELL_FACE_LEFT;

          directions[4][0] = -2;
          directions[4][1] = MESH_CELL_FACE_TOP;
          directions[4][2] = MESH_CELL_FACE_BACK;
          directions[4][3] = MESH_CELL_FACE_RIGHT;

          directions[5][0] = -100;
          directions[5][1] = -100;
          directions[5][2] = -100;
          directions[5][3] = -100;

        }
      else if (nd == 1)
        {
          pts[0][0] = hex[0];
          pts[0][1] = hex[5];
          pts[0][2] = hex[7];
          pts[0][3] = hex[4];

          pts[1][0] = hex[0];
          pts[1][1] = hex[1];
          pts[1][2] = hex[7];
          pts[1][3] = hex[5];

          pts[2][0] = hex[1];
          pts[2][1] = hex[6];
          pts[2][2] = hex[7];
          pts[2][3] = hex[5];

          pts[3][0] = hex[0];
          pts[3][1] = hex[7];
          pts[3][2] = hex[2];
          pts[3][3] = hex[3];

          pts[4][0] = hex[0];
          pts[4][1] = hex[7];
          pts[4][2] = hex[1];
          pts[4][3] = hex[2];

          pts[5][0] = hex[1];
          pts[5][1] = hex[7];
          pts[5][2] = hex[6];
          pts[5][3] = hex[2];

          directions[0][0] = -2;
          directions[0][1] = MESH_CELL_FACE_TOP;
          directions[0][2] = MESH_CELL_FACE_FRONT;
          directions[0][3] = MESH_CELL_FACE_LEFT;

          directions[1][0] = -5;
          directions[1][1] = -3;
          directions[1][2] = MESH_CELL_FACE_FRONT;
          directions[1][3] = -1;

          directions[2][0] = -6;
          directions[2][1] = MESH_CELL_FACE_TOP;
          directions[2][2] = MESH_CELL_FACE_RIGHT;
          directions[2][3] = -2;

          directions[3][0] = -5;
          directions[3][1] = MESH_CELL_FACE_BACK;
          directions[3][2] = MESH_CELL_FACE_LEFT;
          directions[3][3] = MESH_CELL_FACE_BOTTOM;

          directions[4][0] = -2;
          directions[4][1] = -6;
          directions[4][2] = -4;
          directions[4][3] = MESH_CELL_FACE_BOTTOM;

          directions[5][0] = -3;
          directions[5][1] = MESH_CELL_FACE_BACK;
          directions[5][2] = -5;
          directions[5][3] = MESH_CELL_FACE_RIGHT;

        }
      else if (nd == 2)
        {
          pts[0][0] = hex[0];
          pts[0][1] = hex[4];
          pts[0][2] = hex[5];
          pts[0][3] = hex[6];

          pts[1][0] = hex[0];
          pts[1][1] = hex[3];
          pts[1][2] = hex[7];
          pts[1][3] = hex[6];

          pts[2][0] = hex[0];
          pts[2][1] = hex[7];
          pts[2][2] = hex[4];
          pts[2][3] = hex[6];

          pts[3][0] = hex[0];
          pts[3][1] = hex[1];
          pts[3][2] = hex[2];
          pts[3][3] = hex[5];

          pts[4][0] = hex[0];
          pts[4][1] = hex[3];
          pts[4][2] = hex[6];
          pts[4][3] = hex[2];

          pts[5][0] = hex[0];
          pts[5][1] = hex[6];
          pts[5][2] = hex[5];
          pts[5][3] = hex[2];

          directions[0][0] = MESH_CELL_FACE_FRONT;
          directions[0][1] = MESH_CELL_FACE_TOP;
          directions[0][2] = -3;
          directions[0][3] = -6;

          directions[1][0] = MESH_CELL_FACE_LEFT;
          directions[1][1] = MESH_CELL_FACE_BACK;
          directions[1][2] = -5;
          directions[1][3] = -3;

          directions[2][0] = MESH_CELL_FACE_LEFT;
          directions[2][1] = MESH_CELL_FACE_TOP;
          directions[2][2] = -2;
          directions[2][3] = -1;

          directions[3][0] = MESH_CELL_FACE_BOTTOM;
          directions[3][1] = MESH_CELL_FACE_RIGHT;
          directions[3][2] = MESH_CELL_FACE_FRONT;
          directions[3][3] = -6;

          directions[4][0] = -2;
          directions[4][1] = MESH_CELL_FACE_BACK;
          directions[4][2] = MESH_CELL_FACE_BOTTOM;
          directions[4][3] = -6;

          directions[5][0] = -1;
          directions[5][1] = MESH_CELL_FACE_RIGHT;
          directions[5][2] = -5;
          directions[5][3] = -4;
        }
      else if (nd == 3)
        {
          pts[0][0] = hex[0];
          pts[0][1] = hex[2];
          pts[0][2] = hex[3];
          pts[0][3] = hex[6];

          pts[1][0] = hex[0];
          pts[1][1] = hex[3];
          pts[1][2] = hex[7];
          pts[1][3] = hex[6];

          pts[2][0] = hex[0];
          pts[2][1] = hex[7];
          pts[2][2] = hex[4];
          pts[2][3] = hex[6];

          pts[3][0] = hex[0];
          pts[3][1] = hex[5];
          pts[3][2] = hex[6];
          pts[3][3] = hex[4];

          pts[4][0] = hex[1];
          pts[4][1] = hex[5];
          pts[4][2] = hex[6];
          pts[4][3] = hex[0];

          pts[5][0] = hex[1];
          pts[5][1] = hex[6];
          pts[5][2] = hex[2];
          pts[5][3] = hex[0];

          directions[0][0] = MESH_CELL_FACE_BOTTOM;
          directions[0][1] = MESH_CELL_FACE_BACK;
          directions[0][2] = -6;
          directions[0][3] = -2;

          directions[1][0] = MESH_CELL_FACE_LEFT;
          directions[1][1] = MESH_CELL_FACE_BACK;
          directions[1][2] = -1;
          directions[1][3] = -3;

          directions[2][0] = MESH_CELL_FACE_LEFT;
          directions[2][1] = MESH_CELL_FACE_TOP;
          directions[2][2] = -2;
          directions[2][3] = -4;

          directions[3][0] = -5;
          directions[3][1] = MESH_CELL_FACE_TOP;
          directions[3][2] = MESH_CELL_FACE_FRONT;
          directions[3][3] = -3;

          directions[4][0] = MESH_CELL_FACE_RIGHT;
          directions[4][1] = -4;
          directions[4][2] = MESH_CELL_FACE_FRONT;
          directions[4][3] = -6;

          directions[5][0] = MESH_CELL_FACE_RIGHT;
          directions[5][1] = -1;
          directions[5][2] = -5;
          directions[5][3] = MESH_CELL_FACE_BOTTOM;
        }
      else
        {
          /* This is just to avoid compiler warnings */

          Die(FUNCTION_NAME, "Shouldn't be here");

          pts[0][0] = -1;
          pts[0][1] = -1;
          pts[0][2] = -1;
          pts[0][3] = -1;

          pts[1][0] = -1;
          pts[1][1] = -1;
          pts[1][2] = -1;
          pts[1][3] = -1;

          pts[2][0] = -1;
          pts[2][1] = -1;
          pts[2][2] = -1;
          pts[2][3] = -1;

          pts[3][0] = -1;
          pts[3][1] = -1;
          pts[3][2] = -1;
          pts[3][3] = -1;

          pts[4][0] = -1;
          pts[4][1] = -1;
          pts[4][2] = -1;
          pts[4][3] = -1;

          pts[5][0] = -1;
          pts[5][1] = -1;
          pts[5][2] = -1;
          pts[5][3] = -1;

          directions[0][0] = -1;
          directions[0][1] = -1;
          directions[0][2] = -1;
          directions[0][3] = -1;

          directions[1][0] = -1;
          directions[1][1] = -1;
          directions[1][2] = -1;
          directions[1][3] = -1;

          directions[2][0] = -1;
          directions[2][1] = -1;
          directions[2][2] = -1;
          directions[2][3] = -1;

          directions[3][0] = -1;
          directions[3][1] = -1;
          directions[3][2] = -1;
          directions[3][3] = -1;

          directions[4][0] = -1;
          directions[4][1] = -1;
          directions[4][2] = -1;
          directions[4][3] = -1;

          directions[5][0] = -1;
          directions[5][1] = -1;
          directions[5][2] = -1;
          directions[5][3] = -1;
        }

      break;

    default:

      Die(FUNCTION_NAME, "Invalid cell type %ld", celltype);

    }

  /* Create cells */

  for (i = 0; i < nc; i++)
    {
      /* Set bounding box */

      p0 = pts[i][0];
      p1 = pts[i][1];
      p2 = pts[i][2];
      p3 = pts[i][3];

      p0 = pointlist + 3*p0;
      p1 = pointlist + 3*p1;
      p2 = pointlist + 3*p2;
      p3 = pointlist + 3*p3;

      /* Test that tetrahedron will be like */

      /*         3                      3            */
      /*        /|\                    /|\           */
      /*       / | \                  / | \          */
      /*  OK: 0 -|- 2        Not OK: 0 -|- 1         */
      /*       \ | /                  \ | /          */
      /*        \|/                    \|/           */
      /*         1                      2            */
      /*                                             */

      u0 = RDB[p1 + 0] - RDB[p0 + 0];
      v0 = RDB[p1 + 1] - RDB[p0 + 1];
      w0 = RDB[p1 + 2] - RDB[p0 + 2];

      u1 = RDB[p2 + 0] - RDB[p1 + 0];
      v1 = RDB[p2 + 1] - RDB[p1 + 1];
      w1 = RDB[p2 + 2] - RDB[p1 + 2];

      /* From p0 to p3 */

      u = RDB[p3 + 0] - RDB[p0 + 0];
      v = RDB[p3 + 1] - RDB[p0 + 1];
      w = RDB[p3 + 2] - RDB[p0 + 2];

      /* Check scalar triple product */

      if (u*(v0*w1-w0*v1) + v*(w0*u1-u0*w1) + w*(u0*v1-v0*u1) <= 0.0)
        Die(FUNCTION_NAME, "Not OK %f", u*(v0*w1-w0*v1) + v*(w0*u1-u0*w1) + w*(u0*v1-v0*u1));

      /* Bounding boxes are now only in parent cells
      TetPutBoundingBox(ncgns, pointlist, tmpface);
      */
      tet = (long)RDB[tetlist + *totidx - nc + i];

      /* Put points to Tet */
      /*
      printf("Tet %ld consists of %ld %ld %ld %ld\n", *totidx - nc + i, p0, p1, p2, p3);
      */
      WDB[tet + TET_POINTS + 0] = (double)p0;
      WDB[tet + TET_POINTS + 1] = (double)p1;
      WDB[tet + TET_POINTS + 2] = (double)p2;
      WDB[tet + TET_POINTS + 3] = (double)p3;

    }

  /* Set neighbours */

  for (i = 0; i < nc; i++)
    {
      /* Get pointer to tet */

      tet = (long)RDB[tetlist + *totidx - nc + i];

      /* Loop over surfaces (mostly to set neighbours?) */

      for (j = 0; j < 4; j++)
        {

          /* Get pointers to points */

          p0 = (long)RDB[tet + TET_POINTS + face[j][0]];
          p1 = (long)RDB[tet + TET_POINTS + face[j][1]];
          p2 = (long)RDB[tet + TET_POINTS + face[j][2]];


          /* Get first flag and unset it for future */

          if (directions[i][j] >= 0)
            {
              /* Out of cell */

              /* Get direction */

              direction = directions[i][j];

              /* Get face index */

              fidx = hexfaces[direction];
              /*
              fprintf(outp , "%ld is neighbour through surface %ld (directions[%ld][%ld] = %ld)\n",tet, fidx, i, j, direction);
              */
              /* Put this as one of the neighbours of that surface */

              for (k = 0; k < 4; k++)
                {
                  /* If there is a neighbour already set for this k-value, skip it */

                  if (newnbrs[k][fidx] != 0)
                    {
                      /*
                      fprintf(outp, "Index %ld of face index %ld already contains %ld\n", k, fidx, newnbrs[k][fidx]);
                      */
                      continue;
                    }
                  /*
                  fprintf(outp, "Storing %ld to index %ld of face index %ld\n", tet, k, fidx);
                  */

                  /* Store this as one of the tets sharing initial surface fidx */

                  newnbrs[k][fidx] = tet;

                  /* Store this tets face (j) sharing part of surface fidx */

                  newnbrfaces[k][fidx] = j;

                  break;
                }

              /* Check that we were able to store it */

              if (k >= 4)
                Die(FUNCTION_NAME, "All tets sharing surface %ld were already taken", fidx);

              /* Put this faces initial surface index as neighbour */

              WDB[tet  + TET_NEIGHBOURS + j] = -(double)fidx - 1;
            }
          else
            {
              /* Internal face so the neighbour should be one of the tets */
              /* created here */

              io = -(directions[i][j] + 1);

              /* Get pointer to other tet */

              otet = (long)RDB[tetlist + *totidx - nc + io];

              /* Put other tet as neighbour through this face */

              WDB[tet  + TET_NEIGHBOURS + j] = (double)otet;
              /*
              fprintf(outp, "%ld is neighbour of %ld\n", tet, otet);
              */
            }
        }
    }
}

/*****************************************************************************/
