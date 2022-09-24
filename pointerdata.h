/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pointerdata.h                                  */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2010/09/15 (JLe)                                           */
/* Version:       2.0.1                                                      */
/*                                                                           */
/* Description: Contains pointers and global array used by block functions   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

/* List header structure */

#define  LIST_BLOCK_SIZE     10

#define  LIST_PTR_NAME        1
#define  LIST_OPTIONS         2
#define  LIST_OWNER           3
#define  LIST_ITEM_SIZE       4
#define  LIST_N_ITEMS         5
#define  LIST_PTR_FIRST_ITEM  6
#define  LIST_PTR_LAST_ITEM   7
#define  LIST_PTR_ITEM_ARRAY  8
#define  LIST_PTR_NEXT        9

/* List item structure */

#define  ITEM_BLOCK_SIZE      5

#define  ITEM_PTR_DATA        0
#define  ITEM_IDX             1
#define  ITEM_PTR_LIST        2
#define  ITEM_PTR_PREV        3
#define  ITEM_PTR_NEXT        4

/* Header data for pointer block */

#define POINTER_HEADER_SIZE        4

#define POINTER_DATA_SIZE          0
#define POINTER_REAL_DATA_SIZE     1
#define POINTER_N_LISTS            2
#define POINTER_PTR_FIRST_LIST     3

/* Buffer size for memory allocation */

#define POINTER_BUFF_SIZE  1000*ITEM_BLOCK_SIZE

/* Global pointer array */

long *POINTER;

/* Variables for Open MP */

long omp_tasks;
long omp_id;

/*****************************************************************************/
