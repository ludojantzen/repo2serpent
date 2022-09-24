/*****************************************************************************/

/***** External libraries ****************************************************/

#define _XOPEN_SOURCE 600
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <sys/types.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/timeb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <signal.h>
#include <unistd.h>
#include <netdb.h>
#include <errno.h>

#ifndef NO_GFX_MODE
#include <gd.h>
#endif

/*****************************************************************************/

/***** Parallel calculation **************************************************/

#ifdef MPI

/* MPI library header */

#include <mpi.h>

/* mpirun executable path */

#define MPIRUN_PATH "mpirun"

#endif

#ifdef OPEN_MP

/* OpenMP library header */

#include <omp.h>

#endif

/* MPI mode: 1 = divide source size, 2 = divide number of active cycles */

#define MPI_MODE1

#define OLD_IFP

/*****************************************************************************/

/***** Code name and version *************************************************/

#define CODE_NAME     "Serpent"
#define CODE_VERSION  "2.1.32"
#define CODE_DATE     "May 16, 2019"
#define CODE_AUTHOR   "serpent@vtt.fi"

/*****************************************************************************/

/***** Constants *************************************************************/

/* Zero and infinity */

#define INFTY       1e+37
#define ZERO        1e-37

/* Natural constants */

#define M_NEUTRON   1.0086649670000E+00  /* neutron mass in amu              */
#define N_AVOGADRO  6.0220434469282E-01  /* Avogadro's constant in 1E-24     */
#define SHAKE       1E-8                 /* ancient time unit used at LANL   */
#define KELVIN      8.6173E-11           /* MeV per kelvin                   */
#define E_RESTMASS  0.5109989            /* Electron rest mass in MeV        */
#define LOG_E_RESTMASS -0.67138784142263 /* Log of electron rest mass        */
#define MEV         1.6021765314E-13     /* J per MeV                        */
#define ETOV2       1.9131290731E+18     /* square speed (cm2/s2) per energy */
#define BARN        1E-24                /* barn in cm2                      */
#define FS_CONST    7.2973525698E-03     /* Fine-structure constant (~1/137) */
#define GRAYH       5.7678382E-07        /* MeV/s to Gy/h                    */

#ifdef wanhat

#define SPD_C       29979245800.0        /* Speed of light in cm/s           */
#define NEUTRON_E0  939.668854490302     /* Neutron rest mass in MeV         */

#else
#define SPD_C       29979250000.0        /* Speed of light in cm/s           */
#define NEUTRON_E0  9.3958000000000E+02  /* Neutron rest mass in MeV         */

#endif

/* U-235 fission Q-value and energy deposition per fission */

#define U235_FISSQ 1.9372E+02
#define U235_FISSE 202.27*MEV

/* Math constants */

#define LOG2    0.69314718055995
#define PI      3.14159265358979
#define SQRT2   1.41421356237309
#define SQRT3   1.73205080756888
#define SIN30   0.50000000000000
#define COS30   0.86602540378444
#define TAN30   0.57735026918963
#define SIN45   0.70710678118655
#define SIN60   0.86602540378444
#define COS60   0.50000000000000
#define SQRTPI  1.77245385090552

/* Random number stride between neutrons = 2^STRIDE (Fast implementations */
/* for STRIDE values 12-19 are implemented for use in serial mode. (TVi)) */

#define STRIDE 16

/* Extrapolation length for boundary distances */

#define EXTRAP_L 1E-6

/* very negative null pointer for data arrays */

#define NULLPTR -1E+6

/* yes and no */

#define YES 1
#define NO  0

/* Maximum array sizes */

#define MAX_STR                  256
#define MAX_FP_NUCLIDES          1500
#define MAX_ISOTOPES             1000
#define MAX_INPUT_PARAMS         10000000
#define MAX_CELL_SURFACES        10000
#define MAX_LATTICE_ITEMS        100000
#define MAX_SURFACE_PARAMS       300
#define MAX_EGRID_NE             100000000
#define MAX_PRECURSOR_GROUPS     8
#define MAX_GEOMETRY_LEVELS      10000
#define MAX_EXT_K_GEN            10
#define MAX_GENERATIONS          1000000000
#define MAX_RMX_BUFF             10000
#define MAX_PLOT_COLORS          256

/* Tracking errors */

#define TRACK_ERR_INF_LOOP     1
#define TRACK_ERR_CELL_SEARCH  2
#define TRACK_ERR_LATTICE      3
#define TRACK_ERR_NO_MATERIAL  4
#define TRACK_ERR_OUTSIDE      5

/* Maximum allowed number of OpenMP threads */

#define MAX_OMP_THREADS 1000

/* Limiting values for parameters (used for sanity checks only) */

#define MAX_XS       1E+12  /* Maximum microscopic cross section  */
#define E_CHECK_MAX  220.0  /* Maximum energy */

/* Energy grid types */

#define EG_INTERP_MODE_LIN 1
#define EG_INTERP_MODE_LOG 2

/* Run mode options */

#define MODE_NORMAL         0
#define MODE_REPLAY         2
#define MODE_MPI            4

/* Simulation modes */

#define SIMULATION_MODE_CRIT    1
#define SIMULATION_MODE_SRC     2
#define SIMULATION_MODE_DYN     3
#define SIMULATION_MODE_DELDYN  4

/* Indexes used for particle balance */

#define BALA_N_SRC_EXT    0
#define BALA_N_SRC_FISS   1
#define BALA_N_SRC_NXN    2
#define BALA_N_SRC_VR     3
#define BALA_N_SRC_TOT    4

#define BALA_N_LOSS_CAPT  0
#define BALA_N_LOSS_FISS  1
#define BALA_N_LOSS_LEAK  2
#define BALA_N_LOSS_CUT   3
#define BALA_N_LOSS_ERR   4
#define BALA_N_LOSS_TOT   5

#define BALA_G_SRC_EXT    0
#define BALA_G_SRC_TTB    1
#define BALA_G_SRC_FLUOR  2
#define BALA_G_SRC_ANNIH  3
#define BALA_G_SRC_NREA   4
#define BALA_G_SRC_VR     5
#define BALA_G_SRC_TOT    6

#define BALA_G_LOSS_CAPT  0
#define BALA_G_LOSS_LEAK  1
#define BALA_G_LOSS_CUT   2
#define BALA_G_LOSS_ERR   3
#define BALA_G_LOSS_TOT   4

/* UFS mode */

#define UFS_MODE_NONE  0
#define UFS_MODE_COL   1
#define UFS_MODE_FLUX  2
#define UFS_MODE_FISS  3
#define UFS_MODE_RMX   4

/* Particle types (indeksit pitää mennä noin että sorttaus tyypin mukaan */
/* menee oikein) */

#define PARTICLE_TYPE_DUMMY      0
#define PARTICLE_TYPE_GAMMA      1
#define PARTICLE_TYPE_NEUTRON    2
#define PARTICLE_TYPE_PRECURSOR  3
#define PARTICLE_TYPE_ELECTRON   4
#define PARTICLE_TYPE_POSITRON   5
#define PARTICLE_TYPE_ALPHA      6

/* K-eff iteration modes */

#define ITER_MODE_NONE    0
#define ITER_MODE_ALBEDO  1
#define ITER_MODE_NUCLIDE 2
#define ITER_MODE_USER    3

/* Domain decomposition */

#define DD_MODE_NONE    0
#define DD_MODE_SIMPLE  1
#define DD_MODE_SECTOR  2
#define DD_MODE_AUTO    3
#define DD_MODE_LAT     4

/* CMM modes for calculating removal xs */

#define CMM_REMXS_IMPL  1
#define CMM_REMXS_ANA   2

/* CMM methods */

#define CMM_METHOD_OLD 1
#define CMM_METHOD_NEW 2

/* Burnup modes (solution method for Bateman's equations) */

#define BUMODE_TTA  1
#define BUMODE_CRAM 2

/* Special run modes in coefficient calculation */

#define SPECIAL_COEF_MODE_HIS_ONLY 1
#define SPECIAL_COEF_MODE_COE_ONLY 2

/* Normalization in burnup mode */

#define BURN_NORM_ALL       1
#define BURN_NORM_BURN      2
#define BURN_NORM_NOT_BURN  3

/* Radioactive decay source sampling modes */

#define RAD_SRC_MODE_ANALOG    1
#define RAD_SRC_MODE_IMPLICIT  2

/* Weight window mesh types */

#define WWD_MESH_TYPE_SSS   1
#define WWD_MESH_TYPE_MCNP  2
#define WWD_MESH_TYPE_ITER  3

/* Weight window iteration types */

#define WWIN_ITER_USR  1
#define WWIN_ITER_GEO  2
#define WWIN_ITER_TRK  3
#define WWIN_ITER_CON  4

/* Maximum number of time importance intervals */

#define TME_IMP_MAX_N  100

/* Array types */

#define DATA_ARRAY     1
#define RES1_ARRAY     2
#define ACE_ARRAY      3
#define PRIVA_ARRAY    4
#define BUF_ARRAY      5
#define RES2_ARRAY     6
#define RES3_ARRAY     7

/* Solid 3D geometry types */

#define SOLID_GEO_TYPE_OF   1
#define SOLID_GEO_TYPE_STL  2
#define SOLID_GEO_TYPE_IFC  3

/* Cell search list types */

#define CELL_SEARCH_LIST_UNI   0
#define CELL_SEARCH_LIST_CELL  1
#define CELL_SEARCH_LIST_SRC   2
#define CELL_SEARCH_LIST_SURF  3

/* STL geometry cell search mode */

#define STL_SEARCH_MODE_NONE  0
#define STL_SEARCH_MODE_FAST  1
#define STL_SEARCH_MODE_SAFE  2

/* Detector flag options */

#define DET_FLAG_OPT_RESET       0
#define DET_FLAG_OPT_SET         1
#define DET_FLAG_OPT_TEST_SET    2
#define DET_FLAG_OPT_TEST_UNSET  3

/* Lattice bin types */

#define DET_LBIN_LAT   1
#define DET_LBIN_PBED  2

/* STL fail flags */

#define STL_RAY_TEST_FAIL_PARA   -1000
#define STL_RAY_TEST_FAIL_EDGE   -2000
#define STL_RAY_TEST_FAIL_MESH   -3000
#define STL_RAY_TEST_FAIL_LOOP   -4000
#define STL_RAY_TEST_FAIL_STUCK  -5000
#define STL_FACET_OVERLAP        -6000

/* XS data types */

#define XS_TYPE_SAB         3
/*
#define XS_TYPE_NONE        0
#define XS_TYPE_CONTINUOUS  1
#define XS_TYPE_DOSIMETRY   2
#define XS_TYPE_DECAY       4
*/

/* THERM interpolation types */

#define THERM_INTERP_NONE      0
#define THERM_INTERP_STOCHMIX  1
#define THERM_INTERP_MAKXSF    2
#define THERM_INTERP_OTF       3

/* Calculation of analog reaction rates */

#define ARR_MODE_NONE  0
#define ARR_MODE_BALA  1
#define ARR_MODE_ALL   2

/* Diffusion coefficient used for calculating homogeneous diffusion flux */

#define HOMOFLUX_DIFF_INF 1
#define HOMOFLUX_DIFF_TRC 2

/* Nuclide types */

#define NUCLIDE_TYPE_NONE        0
#define NUCLIDE_TYPE_TRANSPORT   1
#define NUCLIDE_TYPE_DOSIMETRY   2
#define NUCLIDE_TYPE_SAB         3
#define NUCLIDE_TYPE_DECAY       4
#define NUCLIDE_TYPE_PHOTON      5
#define NUCLIDE_TYPE_DBRC        6
#define NUCLIDE_TYPE_TRANSMUXS   7

/* Nuclide flags */

#define NUCLIDE_FLAG_INITIAL              1
#define NUCLIDE_FLAG_AP                   2
#define NUCLIDE_FLAG_DP                   4
#define NUCLIDE_FLAG_FP                   8
#define NUCLIDE_FLAG_BP                  16
#define NUCLIDE_FLAG_TRANSPORT_DATA      32
#define NUCLIDE_FLAG_DOSIMETRY_DATA      64
#define NUCLIDE_FLAG_DECAY_DATA         128
#define NUCLIDE_FLAG_NFY_DATA           256
#define NUCLIDE_FLAG_SFY_DATA           512
#define NUCLIDE_FLAG_BRA_DATA          1024
#define NUCLIDE_FLAG_URES_AVAIL        2048
#define NUCLIDE_FLAG_URES_USED         4096
#define NUCLIDE_FLAG_FISSILE           8192
#define NUCLIDE_FLAG_SAB_DATA         16384
#define NUCLIDE_FLAG_SRC              32768
#define NUCLIDE_FLAG_PHOTON_DATA      65536
#define NUCLIDE_FLAG_DBRC            131072
#define NUCLIDE_FLAG_DEP             262144
#define NUCLIDE_FLAG_DELNU_PREC      524288
#define NUCLIDE_FLAG_TRANSMU_DATA   1048576
#define NUCLIDE_FLAG_NEW_DAUGHTERS  2097152
#define NUCLIDE_FLAG_TMS            4194304
#define NUCLIDE_FLAG_CRIT_ITER      8388608
#define NUCLIDE_FLAG_EXT_ADENS     16777216
#define NUCLIDE_FLAG_BA            33554432

/* Photon types */

#define PHOTON_TYPE_SRC           0
#define PHOTON_TYPE_TTB           1
#define PHOTON_TYPE_FLUOR         2
#define PHOTON_TYPE_ANNIH         3
#define PHOTON_TYPE_NREA          4

/* Decay spectra */

#define DECAY_SPEC_LINE 1
#define DECAY_SPEC_CONT 2

/* Reaction types */

#define REACTION_TYPE_SUM         1
#define REACTION_TYPE_PARTIAL     2
#define REACTION_TYPE_SPECIAL     3
#define REACTION_TYPE_DECAY       4
#define REACTION_TYPE_TRA_BRANCH  5
#define REACTION_TYPE_DEC_BRANCH  6

/* Macroscopic reaction MT's (also used with detector response functions) */

#define MT_MACRO_TOTXS           -1
#define MT_MACRO_ABSXS           -2
#define MT_MACRO_ELAXS           -3
#define MT_MACRO_HEATXS          -4
#define MT_MACRO_PHOTXS          -5
#define MT_MACRO_FISSXS          -6
#define MT_MACRO_NSF             -7
#define MT_MACRO_FISSE           -8
#define MT_MACRO_MAJORANT        -9
#define MT_MACRO_RECOILE        -10
#define MT_SOURCE_RATE          -11
#define MT_MACRO_HEATPHOTANA    -12

#define MT_NEUTRON_DENSITY      -15
#define MT_MACRO_INLPRODXS      -16
#define MT_LEAK_RATE            -17

#define MT_MACRO_TOTPHOTXS      -25
#define MT_MACRO_HEATPHOTXS     -26
#define MT_PHOTON_PULSE_HEIGHT  -27

#define MT_MACRO_TMP_MAJORANTXS -30

/* Light element productions (only triton included for now) */

#define MT_MACRO_PROTPXS        -53
#define MT_MACRO_DEUTPXS        -54
#define MT_MACRO_TRITPXS        -55
#define MT_MACRO_HE3PXS         -56
#define MT_MACRO_HE4PXS         -57

#define MT_MACRO_HEATTOT        -80

#define MT_USER_DEFINED        -100

/* Time-dependent source rates for different things */

#define MT_PRIMARY_LIVE_SOURCE     -110
#define MT_PRIMARY_DN_SOURCE       -111

#define MT_PRIMARY_DN_SOURCE_G1    -121
#define MT_PRIMARY_DN_SOURCE_G2    -122
#define MT_PRIMARY_DN_SOURCE_G3    -123
#define MT_PRIMARY_DN_SOURCE_G4    -124
#define MT_PRIMARY_DN_SOURCE_G5    -125
#define MT_PRIMARY_DN_SOURCE_G6    -126
#define MT_PRIMARY_DN_SOURCE_G7    -127
#define MT_PRIMARY_DN_SOURCE_G8    -128

#define MT_SECONDARY_PROMPT_SOURCE -130
#define MT_SECONDARY_DN_SOURCE     -131

#define MT_SECONDARY_DN_SOURCE_G1  -131
#define MT_SECONDARY_DN_SOURCE_G2  -132
#define MT_SECONDARY_DN_SOURCE_G3  -133
#define MT_SECONDARY_DN_SOURCE_G4  -134
#define MT_SECONDARY_DN_SOURCE_G5  -135
#define MT_SECONDARY_DN_SOURCE_G6  -136
#define MT_SECONDARY_DN_SOURCE_G7  -137
#define MT_SECONDARY_DN_SOURCE_G8  -138

/* Pre-defined photon attenuation coefficients -201 ... -248 */

#define MT_PHOTON_DOSE         -200

/* Electron MT's */

#define MT_ELECTRON_PE         -301
#define MT_ELECTRON_COMPTON    -302
#define MT_ELECTRON_PP_EL      -303
#define MT_ELECTRON_PP_POS     -304
#define MT_ELECTRON_AUGER      -305

/* Detector output conversion */

#define DET_CONVERT_NONDE  0
#define DET_CONVERT_GDOSE  1

/* Energy grid types */

#define GRID_TYPE_LIN  1
#define GRID_TYPE_LOG  2

/* User energy grid types */

#define EG_TYPE_ARB     1
#define EG_TYPE_UNI_E   2
#define EG_TYPE_UNI_L   3
#define EG_TYPE_PREDEF  4

/* User time binning types */

#define TB_TYPE_ARB       1
#define TB_TYPE_UNI_T     2
#define TB_TYPE_UNI_LOGT  3

/* Angular distribution types */

#define ANG_TYPE_EQUIBIN 1
#define ANG_TYPE_TABULAR 2

/* General options */

#define OPT_USED               1
#define OPT_EXISTS             2
#define OPT_BURN_MAT           4
#define OPT_FISSILE_MAT        8
#define OPT_PHYSICAL_MAT      16
#define OPT_REPLACED_MAT      32
#define OPT_INCLUDE_MAJORANT  64
#define OPT_EXT_ADENS_MAT    128
#define OPT_UNIV_SUB         256

/* Sort modes */

#define SORT_MODE_ASCEND         1
#define SORT_MODE_DESCEND        2
#define SORT_MODE_ASCEND_PRIVA   3
#define SORT_MODE_DESCEND_PRIVA  4

/* Fission yield types (used to identify fission reactions in   */
/* transmu chains, these must be different from all ENDF MT's). */

#define FISSION_YIELD_TYPE_NFY  -1
#define FISSION_YIELD_TYPE_SFY  -2

/* Data sizes in bytes */

#define KILO        1024.0
#define MEGA     1048576.0
#define GIGA  1073741824.0

/* Input parameter types */

#define PTYPE_LOGICAL 1
#define PTYPE_REAL    2
#define PTYPE_INT     3

/* Plotter modes */

#define PLOT_MODE_YZ 1
#define PLOT_MODE_XZ 2
#define PLOT_MODE_XY 3

/* Track plot tail colours */

#define TRACK_PLOT_NCOL 16

/* Boundary conditions */

#define BC_BLACK       1
#define BC_REFLECTIVE  2
#define BC_PERIODIC    3

/* Transformation order */

#define TRANS_ORDER_ROT_TRANS  0
#define TRANS_ORDER_TRANS_ROT  1

/* Transformation types */

#define TRANS_TYPE_UNI   1
#define TRANS_TYPE_SURF  2
#define TRANS_TYPE_FILL  3
#define TRANS_TYPE_DET   4
#define TRANS_TYPE_LAT   5

/* Moving transformation types */

#define TRANS_MOVE_NONE   0
#define TRANS_MOVE_VEL    1
#define TRANS_MOVE_ACC    2

#define TRANS_TLIM_NONE     0
#define TRANS_TLIM_STAY     1
#define TRANS_TLIM_RESET    2
#define TRANS_TLIM_RELEASE  3

/* Geometry is static from time-interval beginning to end */

#define TDEP_TRANS_TYPE_INTERVAL   0

/* Geometry is static from beginning of tracking-loop for particle */
/* to the end of the tracking loop */

#define TDEP_TRANS_TYPE_PARTICLE   1

/* Discontinuity factor region types */

#define ADF_REG_TYPE_FUEL 1
#define ADF_REG_TYPE_REF  2

/* Relative width of ADF mid and corner regions */

#define ADF_MID_WIDTH  0.1
#define ADF_CORN_WIDTH 0.1

/* Surface list operators */

#define SURF_OP_NOT    -1
#define SURF_OP_AND    -2
#define SURF_OP_OR     -3
#define SURF_OP_LEFT   -4
#define SURF_OP_RIGHT  -5

/* Timers */

#define TOT_TIMERS                22

#define TIMER_TRANSPORT            1
#define TIMER_TRANSPORT_ACTIVE     2
#define TIMER_TRANSPORT_TOTAL      3
#define TIMER_TRANSPORT_CYCLE      4
#define TIMER_BURNUP               5
#define TIMER_BURNUP_TOTAL         6
#define TIMER_PROCESS              7
#define TIMER_PROCESS_TOTAL        8
#define TIMER_BATEMAN              9
#define TIMER_BATEMAN_TOTAL       10
#define TIMER_RUNTIME             11
#define TIMER_INIT                12
#define TIMER_INIT_TOTAL          13
#define TIMER_VOLUME_CALC         14
#define TIMER_OMP_PARA            15
#define TIMER_MPI_OVERHEAD        16
#define TIMER_MPI_OVERHEAD_TOTAL  17
#define TIMER_FINIX               18
#define TIMER_DD_OVERHEAD         19
#define TIMER_RMX                 20
#define TIMER_LEAKAGE_CORR        21
#define TIMER_MISC                22

/* Geometry errors */

#define GEOM_ERROR_NO_CELL         -1
#define GEOM_ERROR_MULTIPLE_CELLS  -2
/*

#define GEOM_ERROR_POINTER_ERROR   -3
*/

/* MPI methods */

#define MPI_METH_BC       1
#define MPI_METH_RED      2
#define MPI_METH_ALL_RED  3

/* Mesh types */

#define MESH_TYPE_CARTESIAN    1
#define MESH_TYPE_CYLINDRICAL  2
#define MESH_TYPE_SPHERICAL    3
#define MESH_TYPE_HEXX         4
#define MESH_TYPE_HEXY         5
#define MESH_TYPE_ORTHOGONAL   6
#define MESH_TYPE_ADAPTIVE     7
#define MESH_TYPE_ICYL         8

/* Mesh content */

#define MESH_CONTENT_RES   1
#define MESH_CONTENT_DAT   2
#define MESH_CONTENT_PTR   3
#define MESH_CONTENT_NONE  4

/* Mesh content data types */

#define MESH_CONTENT_DATA_TET  1
#define MESH_CONTENT_DATA_STL  2

/* Predictor and corrector steps */

#define PREDICTOR_STEP  0
#define CORRECTOR_STEP  1

/* Mesh plot types (numbers 1-3 are reserved for default) */

#define MPL_TYPE_FLUXPOW    4
#define MPL_TYPE_COLPT      5
#define MPL_TYPE_COLWGT     6
#define MPL_TYPE_GAMMAHEAT  7
#define MPL_TYPE_DET        8
#define MPL_TYPE_DET_IMP    9
#define MPL_TYPE_FLUXTEMP  10
#define MPL_TYPE_DENSITY   11
#define MPL_TYPE_DT_NEFF   12
#define MPL_TYPE_DT_GEFF   13

/* Color palettes for mesh plots */

#define PALETTE_HOT             1
#define PALETTE_COLD            2
#define PALETTE_HOTCOLD         3
#define PALETTE_JET             4
#define PALETTE_BW              5
#define PALETTE_HSV             6
#define PALETTE_SPRING          7
#define PALETTE_SUMMER          8
#define PALETTE_AUTUMN          9
#define PALETTE_WINTER         10
#define PALETTE_GREEN_PURPLE   11
#define PALETTE_PURPLE_ORANGE  12
#define PALETTE_BLUE_RED       13

/* Color scales */

#define COLOR_SCALE_LIN  1
#define COLOR_SCALE_LOG  2

/* Surface current detector type */

#define SURF_DET_TYPE_SURF 1
#define SURF_DET_TYPE_UNIV 2

/* TMS modes (järjestyksellä on väliä) */

#define TMS_MODE_NONE 0
#define TMS_MODE_MG   1
#define TMS_MODE_CE   2

/* Entropy calculation mode */

#define ENTROPY_CALC_NONE        0
#define ENTROPY_CALC_SRC         1
#define ENTROPY_CALC_ALL         2

/* Memory allocation */

#define MEM_ALLOC    1
#define MEM_REALLOC  2
#define MEM_FREE     3
#define MEM_ALLOW    4
#define MEM_DENY     5

/* Data types for data-interfaces */

#define DATA_TYPE_IFC_ADENS          1

/* Multi-physics interface types */

#define IFC_TYPE_PT_AVG              1
#define IFC_TYPE_REG_MESH            2
#define IFC_TYPE_REG_MESH_MULTIMAT   21
#define IFC_TYPE_REG_MESH_MULTILVL   22
#define IFC_TYPE_REG_MESH_GENERAL    25
#define IFC_TYPE_FUNC                3
#define IFC_TYPE_FET_TEMP            31
#define IFC_TYPE_FET_DENSITY         32
#define IFC_TYPE_TET_MESH            4
#define IFC_TYPE_FUEP                5
#define IFC_TYPE_FPIP                6
#define IFC_TYPE_OPENFOAM            7
#define IFC_TYPE_OF_MAT              8
#define IFC_TYPE_OF_SOLID            9
#define IFC_TYPE_MIN                 IFC_TYPE_PT_AVG
#define IFC_TYPE_MAX                 IFC_TYPE_FET_DENSITY

/* Flags for printing relative or absolute error */

#define IFC_PRINT_ERROR_REL          2
#define IFC_PRINT_ERROR_ABS          3
#define IFC_PRINT_ERROR_BOTH         4

/* Flag for tallying interface output on the same mesh as incoming data uses */

#define IFC_OUTPUT_SAME_MESH         2

/* OpenFOAM file types */

#define OF_FILE_POINTS    1
#define OF_FILE_FACES     2
#define OF_FILE_OWNER     3
#define OF_FILE_NEIGHBOUR 4
#define OF_FILE_TEMP      5
#define OF_FILE_DENSITY   6
#define OF_FILE_MATERIAL  7
#define OF_FILE_MAP       8

/* OpenFOAM internal field types */

#define OF_INTERNAL_FIELD_UNKNOWN     0
#define OF_INTERNAL_FIELD_UNIFORM     1
#define OF_INTERNAL_FIELD_NONUNIFORM  2

/* Material-wise delta-tracking modes */

#define DT_MAT_BLOCK  1
#define DT_MAT_FORCE  2

/* Wielandt iteration modes */

#define WIELANDT_MODE_NONE   0
#define WIELANDT_MODE_FIX_K  1
#define WIELANDT_MODE_FIX_P  2

/* Source files */

#define SRC_FILE_TYPE_SERPENT1       1
#define SRC_FILE_TYPE_S1_RENORM      2
#define SRC_FILE_TYPE_FUSION_PLASMA  3
#define SRC_FILE_TYPE_WGT_RENORM     4

/* Plot only mode */

#define STOP_AFTER_PLOT_NONE    0
#define STOP_AFTER_PLOT_GEOM    1
#define STOP_AFTER_PLOT_TRACKS  2

/* Restart options */

#define RESTART_CHECK     1
#define RESTART_OVERRIDE  2
#define RESTART_REPLACE   3

/* Material divisor flags (keksi noille paremmat nimet) */

#define MAT_DIV_TYPE_NONE    0
#define MAT_DIV_TYPE_PARENT  1
#define MAT_DIV_TYPE_S1      2
#define MAT_DIV_TYPE_NEW     3

/* Materials in burnup output */

#define BURN_OUT_MAT_DIV     1
#define BURN_OUT_MAT_PARENT  2
#define BURN_OUT_MAT_BOTH    3

/* Print interval burnup output */

#define BURN_OUT_PRINT_FINAL 0
#define BURN_OUT_PRINT_ALL   1
#define BURN_OUT_PRINT_NONE  2

/* Tracking modes */

#define TRACK_MODE_DT 1
#define TRACK_MODE_ST 2

/* Critical spectrum calculation */

#define CRIT_SPECTRUM_OLD  0
#define CRIT_SPECTRUM_B1   1
#define CRIT_SPECTRUM_P1   2
#define CRIT_SPECTRUM_FM   3

#define FM_DIFF_INF  1
#define FM_DIFF_TRC  2
#define FM_DIFF_CMM  3

/* Super-imposed detector types */

#define SUPERDET_TYPE_CURRENT 1
#define SUPERDET_TYPE_TLEFLUX 2

/* Detector type flags */

#define DETECTOR_TYPE_CUMU       -1
#define DETECTOR_TYPE_UNI_E      -2
#define DETECTOR_TYPE_UNI_L      -3
#define DETECTOR_TYPE_ADD_BINS   -4
#define DETECTOR_TYPE_IMP_WGT    -5
#define DETECTOR_TYPE_SUM_SCORES -6

#define DETECTOR_TYPE_MULTI       2
#define DETECTOR_TYPE_DIVI        3
#define DETECTOR_TYPE_T_MULTI     4

/* Track types */

#define TRACK_END_STRT   0
#define TRACK_END_SURF   1
#define TRACK_END_LEAK   2
#define TRACK_END_COLL   3
#define TRACK_END_VIRT   4
#define TRACK_END_TCUT   5
#define TRACK_END_WWIN   6

/* Additional types */

#define TRACK_END_CAPT   7
#define TRACK_END_FISS   8
#define TRACK_END_SCAT   9
#define TRACK_END_ECUT  10
#define TRACK_END_WCUT  11
#define TRACK_END_BC    12
#define TRACK_END_FLAG  13
#define TRACK_END_DD    14
#define TRACK_END_TB    15

/* Last track point for plotter */

#define TRACK_PLOT_LAST 1000

/* Fission matrix types */

#define FISSION_MATRIX_TYPE_MAT  1
#define FISSION_MATRIX_TYPE_UNI  2
#define FISSION_MATRIX_TYPE_LVL  3
#define FISSION_MATRIX_TYPE_XYZ  4

/* Event types (0-11 are from TRACK_END...) */

#define EVENT_TYPE_FISS  13

/* Event record flags */

#define RECORD_EVENT_PLOTTER     1
#define RECORD_EVENT_IFP         2
#define RECORD_EVENT_IMPORTANCE  4
#define RECORD_EVENT_SENS        8

/* Photon production in neutron reactions */

#define PHOTON_PROD_ANA  1
#define PHOTON_PROD_IMP  2

/* Photon constants */

#define PHOTON_ZMAX           99
#define PHOTON_NSS_MAX        40
#define PHOTON_NRADTR_MAX  10000

/* Signaling modes for coupled calculation */

#define SIG_MODE_NONE   0
#define SIG_MODE_POSIX  1
#define SIG_MODE_FILE   2
#define SIG_MODE_SOCKET 3

/* Face indices for mesh cells */

#define MESH_CELL_FACE_LEFT    0
#define MESH_CELL_FACE_RIGHT   1
#define MESH_CELL_FACE_FRONT   2
#define MESH_CELL_FACE_BACK    3
#define MESH_CELL_FACE_BOTTOM  4
#define MESH_CELL_FACE_TOP     5

/* Different types of mesh cells */

#define MESH_CELL_TYPE_TET     0
#define MESH_CELL_TYPE_PYRAMID 1
#define MESH_CELL_TYPE_PRISM   2
#define MESH_CELL_TYPE_HEX     3
#define MESH_CELL_TYPE_POLY    4

/* Majorant extra xs types */

#define MAJORANT_EXTRA_FP_POISON_ITER   1
#define MAJORANT_EXTRA_MIXTURE_ITER     2
#define MAJORANT_EXTRA_NUCLIDE_ITER     3
#define MAJORANT_EXTRA_OTF_BURN         4
#define MAJORANT_EXTRA_NUCLIDE_EXT      5

/* Different modes for handling delayed neutron emission */

#define DELNU_MODE_OFF   0
#define DELNU_MODE_ON    1
#define DELNU_MODE_PREC  2

/* Different modes for tracking precursor concentrations */

#define PREC_MODE_NONE  0
#define PREC_MODE_MESH  1
#define PREC_MODE_POINT 2

#define CSPLINE_FAST_SZ 1000

/* Default values for FETs */

#define FET_DEFAULT_ORDER             6

#define FET_TYPES_START               1
#define FET_TYPE_CARTESIAN            1
#define FET_TYPE_CYLINDRICAL          2
#define FET_TYPE_SPHERICAL            3
#define FET_TYPES_END                 3

#define FET_ORIENTATIONS_START        1
#define FET_ORIENTATION_X             1
#define FET_ORIENTATION_Y             2
#define FET_ORIENTATION_Z             3
#define FET_ORIENTATIONS_END          3

#define FET_CALCULATIONS_START        1
#define FET_CALCULATION_STANDARD      1
#define FET_CALCULATION_ORTHONORMAL   2
#define FET_CALCULATIONS_END          2

#define FET_GENERATE_COEF             0
#define FET_USE_COEF                  1

/* Weight window meshes */

#define WWMESH_SRC    1
#define WWMESH_BOUND  2
#define WWMESH_COLL   3

/* Handling of weight window mesh bounds */

#define WW_MESH_BOUNDS_AVG  1
#define WW_MESH_BOUNDS_SEP  2

/* Response matrix solution modes */

#define RMX_MODE_SINGLE  1
#define RMX_MODE_MULTI   2
#define RMX_MODE_GVR     3
#define RMX_MODE_WWG     4

/* Different ENDF interpolation schemes used in ENDFinterp.c */
/* and calls to it */

#define ENDF_INTERP_TYPE_HIST      1
#define ENDF_INTERP_TYPE_LIN_LIN   2
#define ENDF_INTERP_TYPE_LIN_LOG   3
#define ENDF_INTERP_TYPE_LOG_LIN   4
#define ENDF_INTERP_TYPE_LOG_LOG   5

/* Sensitivity calculation */

#define N_SENS_MAT           3

#define SENS_TOT_MAT_IDX     0
#define SENS_SUM_MAT_IDX     1
#define SENS_NON_MAT_IDX     2

#define N_SENS_ZAI           3

#define SENS_TOT_ZAI_IDX     0
#define SENS_SUM_ZAI_IDX     1
#define SENS_NON_ZAI_IDX     2

/* These are not direct indices for the reactions in results arrays  */
/* The real indices are stored in an array in the main data block at */
/* ptr = (long)RDB[SENS_PERT_INDICES]; */

#define N_SENS_PERT     24

#define TOT_REA_IDX     0
#define ELA_SCATT_IDX   1
#define SAB_SCATT_IDX   2
#define INL_SCATT_IDX   3
#define CAPT_IDX        4
#define FISS_IDX        5
#define NXN_SCATT_IDX   6
#define TEMPERATURE_IDX 7
#define SCATT_TEMP_IDX  8
#define NUBAR_TOT_IDX   9
#define NUBAR_PRO_IDX   10
#define NUBAR_DEL_IDX   11
#define CHI_TOT_IDX     12
#define CHI_PRO_IDX     13
#define CHI_DEL_IDX     14
#define ELA_MU_IDX      15
#define INL_MU_IDX      16
#define ELA_P1_IDX      17
#define ELA_P2_IDX      18
#define ELA_P3_IDX      19
#define ELA_P4_IDX      20
#define ELA_P5_IDX      21
#define ELA_P6_IDX      22
#define ELA_P7_IDX      23

#define TOT_ENE_IDX     0

#define COLL_HIT        1
#define COLL_MISS       (-1)

#define SENS_MODE_NONE               0
#define SENS_MODE_BASIC              1

#define SENS_PERT_FLAG_SCATT_MOM     1
#define SENS_PERT_FLAG_CHI           2
#define SENS_PERT_FLAG_NUBAR         4
#define SENS_PERT_FLAG_XS            8
#define SENS_PERT_FLAG_XSMT         16
#define SENS_PERT_FLAG_ELA_MU       32
#define SENS_PERT_FLAG_INL_MU       64
#define SENS_PERT_FLAG_TEMPERATURE 128

#define SENS_RESP_FLAG_KEFF              1
#define SENS_RESP_FLAG_LEFF              2
#define SENS_RESP_FLAG_BEFF              4
#define SENS_RESP_FLAG_VOID              8
#define SENS_RESP_FLAG_BEFF_GROUPS      16
#define SENS_RESP_FLAG_LAMBDA           32
#define SENS_RESP_FLAG_LAMBDA_GROUPS    64

#define SENS_SCORE_FLAG_MAT_TOT      1
#define SENS_SCORE_FLAG_MAT_SUM      2
#define SENS_SCORE_FLAG_ZAI_TOT      4
#define SENS_SCORE_FLAG_ZAI_SUM      8
#define SENS_SCORE_FLAG_ENE_TOT     16
#define SENS_SCORE_FLAG_HIS         32

#define SENS_RESP_TYPE_KEFF          1
#define SENS_RESP_TYPE_LEFF          2
#define SENS_RESP_TYPE_BEFF          3
#define SENS_RESP_TYPE_VOID          4
#define SENS_RESP_TYPE_RATIO         5
#define SENS_RESP_TYPE_LAMBDA        6

#define SENS_SCORE_TYPE_DIRECT       0
#define SENS_SCORE_TYPE_EVENT        1

#define SENS_MU_BIN_TYPE_EQCOS       0
#define SENS_MU_BIN_TYPE_EQANG       1

/* Energy deposition modes */

#define EDEP_MODE_CONSTANT           0
#define EDEP_MODE_MT458              1
#define EDEP_MODE_LOCAL_PHOTON       2
#define EDEP_MODE_NEUTRON_PHOTON     3

/* Different representation types for fission energy release components */

#define FISSE_TYPE_SHER_BECK         0
#define FISSE_TYPE_TAB               1
#define FISSE_TYPE_POLY              2

/*****************************************************************************/

/***** Structures ************************************************************/

/* Complex number */

typedef struct {
  double re;
  double im;
} complex;


/* Matrix in compressed sparse column format*/

struct ccsMatrix{

  long m;            /* rivit */
  long n;            /* sarakket */
  long nnz;          /* nollasta eroavien lkm */

  long *colptr;      /* osoittimet 1. nollasta eroavaan joka sarakkeessa */
                     /* pituus on (n+1) */
  long *rowind;      /* rivi-indeksit */
  long *colind;      /* sarake-ineksit */
  long *rowptr;
  long *next;
  complex *values;   /* nollasta eroavien arvot */

};

/* Data structure to store nuclide data in depletion files */

struct depnuc {
  long ZAI;
  long Z;
  double AW;
  double lambda;
  double dh;
  double sf;
  double gI;
  double ingtox;
  double inhtox;
};

/* Scaled bremsstrahlung cross section data for elements */

typedef struct {
  long sz;          /* Number of data */
  long nE;          /* Size of E and the other arrays */
  long nkappa;      /* Size of kappa */
  long *map;        /* Index map */
  double *E;        /* Electron kinetic energy */
  double *kappa;    /* Photon energy per electron kinetic energy */
  double ***SXSe;   /* Scaled cross sections of electrons (dim: [sz][nE][nkappa]) */
  double ***SXSeT;  /* Transpose of SXSe (dim: [sz][nkappa][nE]) */
  double ***SXSp;   /* Scaled cross sections of positrons (dim: [sz][nE][nkappa]) */
  double ***SXSpT;  /* Transpose of SXSe (dim: [sz][nkappa][nE]) */
} ElBrSXSData;

/* Stopping power data struct for electrons and positrons */

typedef struct  {
  long sz;          /* Number of data */
  long nE;          /* Size of E and the other arrays */
  long *map;        /* Index map */
  double *E;        /* Electron kinetic energy */
  double **SPcole;  /* Collision stopping power of electrons (dim: [sz][nE])*/
  double **SPrade;  /* Radiative stopping power of electrons */
  double **SPtote;  /* Total stopping power of electrons */
  double **SPcolp;  /* Collision stopping power of positrons */
  double **SPradp;  /* Radiative stopping power of positrons */
  double **SPtotp;  /* Total stopping power of positrons */
  double **delta;   /* Density effect correction */
} ElSPData;

/* Ground state configuration data for an element */

typedef struct {
  long Z;          /* Atomic number */
  double I;        /* Mean excitation energy */
  long nss;        /* Number of shells in the ground state */
  long *gsconf;    /* Ground state configuration */
  double *ebi;     /* Electron binding energies for the ground state */
} GSconfigElemData;

/* Ground state configuration data struct */

typedef struct {
  long sz;                    /* Number of data */
  long *map;                  /* Index map */
  GSconfigElemData *elemData; /* Element-wise data */
} GSconfigData;

/*****************************************************************************/

/***** Function prototypes ***************************************************/

/* NOTE: Use (void) instead of () for empty argument list */

#ifdef __cplusplus
extern "C" {
#endif

void AddBranching(long);

void AddBuf(double, double, long, long, long, ...);

void AddBuf1D(double, double, long, long, long);

void AddChains(long, long, long);

void AddDDRes(long, double);

void AddItem(long, long);

void AddFET(const double *const, long, long, double, double, double,
            long, long, long);

void AddMesh(long, double, double, double, double, long);

void AddMeshIdx(long, double, long, long, long, long);

long AddNuclide(char *, long, char *, double, long, long);

void AddPrivateRes(long, double, long);

double *AddPts(double *, long *, const double *, long);

void AddSabData(void);

void AddSearchMesh(long, long, double, double, double, double, double, double);

void AddSortItem(long, long, long, long, long);

void AddStableNuclides(void);

void AddStat(double, long, ...);

long AddSTLPoint(long ***, long, long, long, double, double, double);

void AddValuePair(long, double, double, long);

void AdjustEnergyGrid(long, long, const double *);

long AdjustRMX(long, long);

void AdjustSabData(long);

void AllocFETCache(long, long, long, long);

void AllocInterfaceStat(void);

void AllocMacroXS(void);

void AllocMicroXS(void);

void AllocParticleStack(long, long);

void AllocPrecDet(void);

long AllocPrivateData(long, long);

void AllocValuePair(long);

void AllocStatHistory(long);

void Alpha(double, double, double *);

double AlphaXS(double);

void ApplyGCSymmetries(long);

long AtoI(char *, char *, char *, long);

double AtoF(char *, char *,char *, long);

void AtomicRelaxation(long, long, long, long, double, double, double, double,
                        double, double *, long);

void AverageTransmuXS(long, double, double, long);

void AvgRMX(long);

void AziRot(double, double *, double *, double *, long);

double B1FluxCorr(long gcu, double E);

long B1Flux(long, long, const double *, const double *, const double *,
            const double *, const double *, const double *, const double *,
            double *, double *);

void BanksToStore(void);

long BoundaryConditions(long *, double *, double *, double *, double *,
                        double *, double *, double *, double *, double *,
                        double *, long);

double BranchFrac(long, long, double, long);

void BroadcastIFCData(void);

void BroadCrossSection(long, long, long, double, double, double, double *,
                       long *, long, long);

void BsfunN(long, double *, double *, complex **, complex **, complex **);

double BufMean(long, ...);

double BufN(long, ...);

double BufVal(long, ...);

double BufWgt(long, ...);

void BurnMatCompositions(void);

void BurnMaterials(long, long);

long BurnMatrixSize(long);

void BurnupCycle(void);

complex c_add (complex, complex);

complex c_con (complex);

complex c_div (complex, complex);

complex c_mul (complex, complex);

double  c_norm(complex);

complex c_sub (complex, complex);

complex c_sqrt(complex);

complex c_exp(complex);

void CacheXS(void);

void CalcConvCriteria(double, double, double*, double*, double *, double *);

void CalcMicroGroupXS(void);

void CalculateActivities(void);

void CalculateBytes(void);

void CalculateDTMajorants(void);

void CalculateEntropies(void);

void CalculateLegendreMoments(long , long , double [7]);

void CalculateMasses(void);

void CalculateMGXS(void);

void CalculateRelAlpha(void);

void CalculateRelPopSize(void);

void CalculateTetBoundingBox(long, double [6]);

void CalculateTetCenter(long, long, long, long);

void CalculateTransmuXS(long, long);

void CalculateUresMajorants(long);

void ccsMatrixColSpace(struct ccsMatrix *, long , long );

void ccsMatrixCopy(struct ccsMatrix *, struct ccsMatrix *);

void ccsMatrixIsort(struct ccsMatrix *);

struct ccsMatrix *ccsMatrixNew(long, long, long);

void ccsMatrixFree(struct ccsMatrix  *);

void ccsMatrixPrint(struct ccsMatrix *);

void CellCount(long, long, long, long);

void CellVolumes(void);

void CheckCoefCalc(void);

long CheckDDDomain(long);

void CheckDuplicates(void);

void CheckNuclideData(void);

void CheckPolyhedMesh(long);

void CheckReaListSum(long, long, double, long, long);

void CheckUnused(void);

void ClearBuf(void);

void ClearInterfaceStat(void);

void ClearMicroGroupXS(void);

void ClearPrivateData(long);

void ClearPrivateRes(void);

void ClearRelTransmuXS(void);

void ClearStat(long);

void ClearTransmuXS(void);

void ClearValuePairs(void);

void CloseList(long);

void CoefCycle(void);

void CoefOutput(void);

void CombineActinides(void);

void CombineFissionYields(void);

long CompareStr(long, long);

void ComplexRea(long, long, double *, double, double, double, double *,
                double *, double *, double, double *, double, double *, long);

void ColDet(long, long, double, double, double, double, double, double, double,
            double, double, double, double, long);

void CollectBuf(void);

void CollectBurnData(void);

void CollectDet(void);

void CollectFET(const double *const, long, long);

void CollectDynData(void);

void CollectParallelData(void);

void CollectPrecDet(void);

void CollectResults(void);

void CollectSensResults(void);

void CollectUncResults(void);

void CollectVRMeshData(void);

long Collision(long, long, double, double, double, double *, double *,
               double *, double *, double *, double, long);

void CommandLinePlotter(void);

long CompressSensLabel(long, long, long, long,
                       long);

void ComptonScattering(long, long, long, double *, double, double, double,
                       double *, double *, double *, double, double, long);

void ContribSplit(long);

void CoordExpans(long, double *, double *, double *, double, long);

void CoordTrans(long, double *, double *, double *, double *, double *,
                double *, long);

void CountDynSrc(void);

long CountSensParams(char **, long);

void CovMatrixFromBlock(long *, long **, long, long, long, long);

void CovMatrixFromSingle(long *, long **, long, long, long, long);

void CreateGeometry(void);

long CreateMesh(long, long, long, long, long, long, const double *, long);

void CreateUMSHCells(long);

long CreateUniverse(long, char *, long);

void CSplineConstruct(const double *, const double *, long, double, double,
                      double **);

void CSplineInterpolate(double **, const double *, long, const double *,
                        double *f, long);

void CSplineInterpolate0(const double *, const double *, long, double, double,
                         const double *, double *, long, long);

double CSplineIntegrate(double **, const double *, double *, long, double,
                        double);

double CSplineIntegrate0(const double *, double *, long, double, double,
                         double, double);

void CSplineCumIntegral(double **, const double *, long, const double *,
                        double *, long);

void CSplineCumIntegral0(const double *, const double *, double *, long, double,
                         double, long, const double *, double **);

double CylDis(double, double, double, double, double);

double DataIFCAdens(long, long, long);

double DataIFCXS(long, double, long, long, long);

void DecayMeshPrecDet(void);

void DecayPointPrecDet(void);

void DecomposeElements(void);

void DefaultBraData(void);

void DeinitSocket(void);

double DensityFactor(long, double, double, double, double, long);

void DepletionPolyFit(long, long);

long DetBin(long, long, long, double, double, double, double, double, long);

void DetectorOutput(void);

long DeterministicLeakage(long, long, long, const double *, const double *,
                          const double *, const double *, const double *,
                          const double *, double *, double *);

long DetIdx(long, long, long, long, long, long, long, long, long, long);

double DetResponse(long, long, long, long, double, double, long);

void DFPos(long, double, double, double, long*, long *, long *, long *,
           long *, long *, double *, double *, double *, double *);

double *dfSol(long, long, double **, double *, complex **, complex **,
              complex *);

void DFSolver(long, long, long, const double *, double, double **,
              double **, complex **, complex **, complex *, complex *);

int Die(char *, ...);

void DiffCoefED(long, double, double, double, double, double,double,
                double, double, double, long);

void Disperse(void);

void Disperse2(void);

void DistributeMaterialData(void);

void DistributeDDMajorants(void);

void DistributeDDMaxAdens(void);

void DivideBurnMat(void);

void DivideMeshCell(long, long, long [8], long (*)[8], long [6],
                     long **, long **, long *, long);

void DividePolyhedCell(long, long, long *, long *S);

void DividePolyhedFace(long, long, long, long, long, long *);

void DivideZone(long, long *, long, long, double, double, double);

double DivSubDistance(long, long, double, double, double, double, double, double, double);

double DopMicroXS(long, long, double, double *, double, long);

void DopplerBroad(void);

double DTMajorant(long, double, long);

long DuplicateItem(long);

long DuplicateParticle(long, long);

long EBlockFromBank(long, long);

void EBlockToBank(long, long);

void ElasticScattering(long, long, long, double *, double *, double *, double *, long);

void Element(long, char *, char *, char *);

void ElBrSXSDataAlloc(ElBrSXSData *, long, long, long);

void ElBrSXSDataFree(ElBrSXSData *);

long ElBrSXSDataGetIdx(ElBrSXSData *, long);

void ElBrSXSDataRead(ElBrSXSData *);

void ElSPDataAlloc(ElSPData *, long, long);

void ElSPDataFree(ElSPData *);

long ElSPDataGetIdx(ElSPData *, long);

void ElSPDataRead(ElSPData *);

double ENDFColF(long, char *);

long ENDFColI(long, char *);

double ENDFInterp(long, double, double, double, double, double);

long ENDFMat(char *line);

long ENDFMF(char *);

long ENDFMT(char *);

void ENDFNewLine(char *, char *, FILE *);

long EneImportance(long, double, double, double, double, double, double,
                   double, double, double *, double, long);

void Error(long, ...);

void EstimateRuntime(void);

long EventFromBank(long, long);

void EventsToSensitivity(long, double, long, double, long);

void EventToBank(long, long);

double ExpandFE(const double *const, const double *const, double, double,
                double, long);

void ExpandPrivateArrays(void);

void ExpoDecFit(void);

void ExtractSensLabel(long, long *, long *, long *, long *, double *);

void FEIFC(const double *const, const double *const, long, double *, double *,
           double, double, double, long, long);

void FETFinalize(const double *const, long, long, long, double *, double *);

long FETIdx(const double *const, long, long, long);

void FillSTLMesh(long, long, double, double, double);

void FinalizeCCPrecDet(void);

void FinalizeMPI(void);

void FinalizeUMSH(long);

void FindInterfaceRegions(long, long, long, long, double, double, double, long);

long FindLatticeRegion(long, long, double *, double *, double *, long *, long);

void FindMaterialPointers(void);

long FindNestRegion(long, long, double, double, double, long);

long FindNuclideData(char *, long, char *, double, long, long);

long FindPBRegion(long, long, double *, double *, double *, long *, long *,
                  long);

long FindRegMeshIFCIndex(long, double, double, double, long);

void FindRowIndexes(struct ccsMatrix *);

long FindSensEIndex(double);

long FindSensMuIndex(double);

void FindSensIndices(long, long, double, long *, long *, long *);

long FindSensZAIIndex(long);

long FindTetCell(long, double, double, double, long);

long FindUniverseCell(long, double, double, double, long *, long);

long FindSTLSolid(long, double, double, double, double, double, double,
                  long, long);

long *FindXSLimits(long, long, double, double, double, double *, long *);

long FinixPtrFromFpe(long);

long FirstItem(long);

double FissE(long, double, long);

double FissEComp(long, long, long, double, long);

long FissMtxIndex(long, long);

void FissMtxOutput(void);

void Fission(long, long, long, double *, double, double, double, double,
             double *, double *, double *, double, double *, long);

void FixHexMesh(long, long, long);

void FixPolyhedMesh(long);

void FlushBank(void);

void FlushPrecSource(void);

void FormTransmuPaths(long, long, double, double, long, long);

void FreeMem(void);

long FromBank(long);

long FromCommonQue(long);

long FromSrc(long);

long FromStack(long, long);

long FromStore(long, long);

long FromTrkBank(long);

long FromQue(long);

long GaussianSubst(struct ccsMatrix *, complex *, complex *, complex *);

void GetBankedPrecursors(void);

void GetBurnIDs(void);

void GetSensEventType(long, long, double, long*, long*, long*);

long GetSignalFromSupervisor(void);

double *GetImportantPts(long, long *, long *, long *);

long GetLatticeIndexes(double, double, double, double, double, double, long *,
                       long *, long *, long);

long GetLineNumber(char *, long);

char **GetParams(char *, char *, long *, long *, long, long, char *);

void GetRequiredCovMatrices(long *);

void GetSensEventType(long, long, double, long*, long*, long*);

long GetSignalFromSupervisor(void);

double GetTemp(long, long);

char *GetText(long);

long GeoImportance(long, long, long, double, double, double, double, double,
                   double, double, double *, double, long);

void GeometryPlotter(long);

void GetIFPAncestorData(long, double *, long *, double *, double *);

double SensBufVal(long, long, long, long, long, long);

long SensCollision(long, long, long, double, long);

double GridFactor(long, double, long);

long GridSearch(long, double);

void GSconfigDataFree(GSconfigData *);

long GSconfigDataGetIdx(GSconfigData *, long);

void GSconfigDataRead(GSconfigData *);

void hessFactorization(long, complex **, complex **, complex **);

void HexNewTetFace(long, long, long, long (*)[6], long [6], long [6], long,
                   long, long, long, long, long *, long *,
                   long *, long, long);

void HexRotateCell(long *, long (*)[8]);

void HexRotateFace(long *, long, long);

double HisMean(long, long, ...);

double HisRelErr(long, long, ...);

double HisVal(long, long, ...);

void HomoFlux(long);

long ICMIdx(long, double, double, double, double *, double *, double *,
            double *, double *, double *, double *, double *, double *);

char *IdxStr(long, long);

void IFCPoint(long, double *, double *, double, long);

void ImportanceSolver(void);

long InCell(long, double, double, double, long, long);

void InelasticScattering(long, double *, double *, double *, double *, long);

void InitData(void);

void InitHistories(void);

void InitialCritSrc(void);

void InitMPI(int, char **);

void InitOMP(void);

void InitPrecDet(void);

void InitPrecDetSource(void);

void InitSignal(void);

void InitSocket(void);

long InSuperCell(long, long, double, double, double, long);

long InterpolateData(const double *, double *, long, const double *,
                     const double *, long, long, long *, long *, long);

void InterpolateFissE(double *, long, long);

double InterpolateMuPDF(long, double);

void InterpolateNubar(double *, long);

long InterpolateSab(long);

void IntersectionList(long, long *, long);

long InTetCell(long, double, double, double, long);

long InUMSHCell(long, long, double, double, double, long);

void InvokeBranch(long);

void Involute(long, long, const double *, long *, double *, double, double,
              double, double, double, double);

void IsotopeFractions(long);

long IsotoZAI(char *);

void IsotropicDirection(double *, double *, double *, long);

void IterateCC(void);

void IterateExternal(void);

void IterateKeff(void);

double IterNucXS(long, double, long, long, long);

void KleinNishina(double, double *, double *, long);

double LagrangeInterpCubic(const double *, const double *, long, double, long);

long LastItem(long);

void LatCellOrig(long, long, double *, double *, double *);

void Leak(long, double, double, double, double, double, double, double,
          double, long);

void LeakDet(long, long, double, double, double, double,
             double, double, double, double, double, long);

void LevelScattering(long, double *, double *, double *, double *, long);

long LIFOListSize(long, long);

void LinkGCUMaterials(long, long);

void LinkInterfaceMaterials(long);

void LinkReactions(void);

void LinkSabData(void);

void LUdecomposition(long, complex **, complex *, complex *);

double MacroUresCorr(long, double, double, long);

double MacroXS(long, double, long);

double MajorantXS(long, double, long);

double *MakeArray(double, double, long, long);

struct ccsMatrix *MakeBurnMatrix(long, long);

struct ccsMatrix *MakeBurnMatrixMSR(long, long, double **, long);

void MakeDepletionZones(long, long, long, long, long, long *, double, double,
                        double, long);

long MakeEnergyGrid(long, long, long, long, const double *, long);

void MakePalette(long *, long *, long *, long, long);

void MakeRing(long);

void MaterialBurnup(long, double *, double *, double, double, long, long);

void MaterialTotals(void);

void MaterialVolumes(void);

void MatlabOutput(void);

void matProduct(long, long, long, complex **, complex **, complex **);

long MatPtr(long, long);

void MatPos(void);

double *MatrixExponential(struct ccsMatrix *, double *, double);

void MaxSrcImp(long, double,  double,  double,  double, long);

void MaxSurfDimensions(long, double *, double *, double *, double *,
                          double *, double *);

double MaxwellEnergy(double, long);

double Mean(long, ...);

void *Mem(long, ...);

long MemCount(void);

void MeshCellBounds(long, long, double *, double *, double *, double *,
                    double *, double *);

void MeshCellConnectDiags(long [8], long (*)[8], long);

void MeshCellFromCGNS(long, long, long *, long *, long (*)[4], long);

void MeshCellGetFace(long [8], long *, long, long);

void MeshCellGetFacePos(long *, long, long);

void MeshCellIndirection(long*, long);

void MeshCellRotateLists(long [8], long (*)[4], long *, long);

long MeshCellType(long, long);

double MeshCellVol(long, double, double, double);

long MeshIndex(long, double, double, double, double);

void MeshPlotter(void);

long MeshPtr(long, double, double, double);

double MeshTot(long);

double MeshVal(long, double, double, double);

double MGXS(long, double, long);

void MicroCalc(void);

void MicroDepOutput(void);

double MicroMajorantXS(long, double, long);

double MicroXS(long, double, long);

double MinXS(long, double, long);

void MORAOutput(void);

long MoveDT(long, double, double, long *, double *, double *, double *,
            double *, double *, double *, double *, double *, double *,
            double *, double *, double, long);
void MoveIFC(void);

void MoveItemFirst(long);

void MoveItemLast(long);

void MoveItemRight(long);

long MoveST(long, double, double, long *, double *, double *, double *,
            double *, double *, double, double, double, long);

void MoveStore(void);

void MPITransfer(double *, double *, long, long, long);

long MyParallelMat(long, long);

double NearestBoundary(long);

double NearestPBSurf(long, double, double, double, double, double, double,
                     long);

double NearestMeshBoundary(long, double, double, double, double, double,
                           double, long *);

double NearestSTLSurf(long, double, double, double, double, double, double,
                      long);

double NearestUMSHSurf(long, double, double, double, double, double, double,
                       long);

void NestVolumes(void);

long NewItem(long, long);

long NewLIFOItem(long, long);

void NewReaList(long, long);

long NewStat(char *, long, ...);

long NextReaction(long, long *, double *, double *, double *, long);

long NextWord(char *, char *);

void NFKERMA(long, double *);

long NoCorrector(void);

double NormCoef(long);

void NormalizeCompositions(void);

void NormalizeCritSrc(void);

long NormalizeDynSrc(void);

void NormalizeImp(long);

void NormalizePrecDet(void);

void Note(long, char *, ...);

double Nubar(long, double, long);

long NumericGauss(struct ccsMatrix *, complex *, complex, complex *);

long NumericStr(char *);

void Nxn(long, long, double *, double, double, double, double *, double *,
         double *, double, double *, double, double *, long);

void OMPResetComp(long);

void OMPSetComp(long);

long OMPTestComp(void);

FILE *OpenDataFile(long, char *);

void OTFBurnTTA(long);

double OTFBurnXS(long, double, long, long);

double OTFSabXS(long, double, double, long);

void OTFSabScattering(long, double *, double *, double *, double *, long);

void OverrideIDs(void);

void PairProduction(long, long, long, double, double, double, double, double,
                    double, double, double, double, long);

void parlett(long, complex *, complex **, complex **, complex **);

long ParseCommandLine(int, char **);

void ParticlesFromStore(void);

void Pause(char *);

double PhiDis(double, double, double, double, double);

void Photoelectric(long, long, long, double, double, double, double, double,
                   double, double, double, double, long);

double PhotonMacroXS(long, double, long);

double PhotonMicroXS(long, double, long);

void PhotonProd(long, long, long, double, double, double, double, double,
                double, double, double, double, long);

long PlotImage(long, long *, long, long, double, double, double, double,
               double, double, double, long, long, long, long, double, double,
               double, long);

#ifndef NO_GFX_MODE

long PlotTracks(long, gdImagePtr, long, double, double, double, double,
                double, double, long, long);

#endif

void PoisonEq(void);

double PoisonXS(long, double, long, long);

double PolarAngle(double, double);

long PolyPInF(long, long *, long);

void PolynomialLegendre(long, double, long, long, long);

void PolynomialZernike(long, double, double, long, long, long);

long PolySameFace(long *, long * , long);

void PosAnnihilation(long, long, double, double, double, double, double,
                     double, double, double, double, long);

double PotCorr(long, double, double);

void PreallocMem(long, long);

void PrecDet(long, long, double, double, double, double, double, double,
             double, double, long);

void PrecursorPopControl(void);

void PrecursorsToStore(void);

void PrepareCCIter(void);

void PrepareCellSearchMesh(void);

void PrepareTransportCycle(void);

void PreSort(void);

void PreTrans(void);

void PrintCellStats(void);

void PrintCoeVals(FILE *, long);

void PrintCompositions(long);

void PrintCoreDistr(void);

void PrintCycleOutput(void);

void PrintDepMatrix(long, struct ccsMatrix *, double t, double *, double *,
                    long);

void PrintDepOutput(void);

void PrintDepVals(FILE *, char *, struct depnuc *, long, long, double **,
                  double *, double *);

void PrintGammaSpectra(void);

void PrintGeometryData(void);

void PrintHistoryOutput(void);

void PrintInterfaceOutput(void);

void PrintMaterialData(void);

void PrintMeshCell(long [8], long);

void PrintMixtures(void);

void PrintMVar(FILE *, long);

void PrintNuclideData(long, long);

void PrintPBData(void);

void PrintPrecDet(void);

void PrintProgress(long, long);

void PrintReactionLists(void);

void PrintSrcImp(void);

void PrintTitle(void);

void PrintTMSDiagnostics(void);

void PrintValues(FILE *, char *, long, long, long, long, long, long);

void ProbeVRMesh(long);

void ProcessBC(void);

void ProcessBUInterface(void);

void ProcessBurnMat(void);

void ProcessBurnupEGroups(void);

void ProcessCells(void);

void ProcessCellMesh(void);

void ProcessCovarianceData(void);

void ProcessCPD(void);

void ProcessComplementCells(void);

void ProcessCompton(long, long);

void ProcessDataInterfaces(void);

void ProcessDecayData(long);

void ProcessDecaySrc(void);

void ProcessDepHis(void);

void ProcessDetectors(void);

void ProcessDivisors(void);

void ProcessDIX(long);

void ProcessEDistributions(long, long);

void ProcessElectrons(void);

void ProcessEntropy(void);

void ProcessEvents(void);

void ProcessFissionYields(long);

void ProcessFissMtx(void);

void ProcessGC(void);

void ProcessICM(void);

void ProcessIFCFB(long, long);

void ProcessIFCFETMaterial(long, long);

void ProcessIFCFunc(long, long);

void ProcessIFCPtAvg(long, long);

void ProcessIFCRegMesh(long, long, long);

void ProcessIFCTetMesh(long, long);

void ProcessInelaXS(long);

void ProcessInterface(long);

void ProcessInventory(void);

void ProcessIterNucs(void);

void ProcessLattices(void);

void ProcessMaterials(void);

void ProcessMeshPlots(void);

void ProcessMixture(long, long);

void ProcessMSR(void);

void ProcessMuDistributions(long);

void ProcessNests(void);

void ProcessNubarData(long);

void ProcessNuclides(void);

void ProcessOTFBurn(void);

void ProcessPairProduction(long, long);

void ProcessPBGeometry(void);

void ProcessPhotoelectric(long, long);

void ProcessPhotoelectricFluorescenceCDF(double *, long);

void ProcessPhotonAtt(void);

void ProcessPhotonEcut(void);

void ProcessPhotonProd(long);

void ProcessPhotonRea(long);

void ProcessPoisons(void);

void ProcessPrecDet(void);

void ProcessRayleigh(long, long);

void ProcessReactionLists(void);

void ProcessRelaxation(void);

void ProcessReprocessors(void);

void ProcessRMX(void);

void ProcessSensEBlocks(void);

void ProcessSensitivities(void);

void ProcessSensStats(void);

void ProcessSingleInterface(long, long, long);

void ProcessSources(void);

void ProcessStats(void);

void ProcessSTLGeometry(void);

void ProcessSurfaces(void);

void ProcessSymmetries(void);

void ProcessTmpData(void);

void ProcessTimeBins(void);

void ProcessTransformations(void);

void ProcessTranspCorr(void);

void ProcessTTB(ElBrSXSData *, ElSPData *);

void ProcessUMSHGeometry(void);

void ProcessUresData(long);

void ProcessUserEGrids(void);

void ProcessVR(void);

void ProcessXSData(void);

void PulseDet(long, long, double, double, double, double, double, long);

void PutCompositions(void);

void PutMeshIdx(long, double, long, long, long);

void PutOTFBurnConc(void);

void PutPlotColors(long, long *, long *, long *);

void PutPoisonConc(void);

long PutText(char *);

void QRfactorization(long, complex **, complex **, complex **);

long RadGammaSrc(long, long, double *, double *, long*, long);

double Rand64(unsigned long *);

double RandF(long);

double RandNorm(double, double, long);

void RayleighScattering(long, double, double *, double *, double *, long);

void ReactionCount(void);

void ReactionCutoff(void);

long ReactionTargetZAI(long);

char *ReactionMT(long, long);

void ReadACEFile(long);

void ReadBRAFile(void);

void ReadCOVERXFile(long);

void ReadDataInterfaces(void);

void ReadDecayFile(void);

void ReadDirectoryFile(void);

void ReadFissionYields(void);

long ReadIFCBins(long, FILE *, double, long);

void ReadIFCFB(long, long);

void ReadIFCFBLims(FILE *, long, long);

void ReadIFCFE(long, long);

void ReadIFCFunc(long, long);

void ReadIFCOFMesh(long, long);

void ReadIFCPtAvg(long, long);

void ReadIFCRegMesh(long, long);

void ReadIFCTetMesh(long, long);

long ReadInfix(long, long *, long *);

void ReadInput(char *);

void ReadInterface(long, long);

double ReadMesh(long, long, long, long);

long ReadMeshPtr(long, long, long, long);

void ReadPendfData(long, FILE *, long, long);

void ReadOFPatches(long);

char *ReadOFData(FILE *, long);

void ReadOFDensities(long, long);

void ReadOFFaces(long);

void ReadOFHeader(FILE *, long *, long *, long *, long *, double *);

void ReadOFMapping(long);

void ReadOFMaterials(long);

void ReadOFMesh(long);

void ReadOFOwnNbr(long);

void ReadOFPoints(long);

void ReadOFTemperatures(long, long);

void ReadPBGeometry(void);

void ReadPhotonData(void);

void ReadPlasmaSrc(long);

void ReadRestartFile(long);

void ReadSourceFile(long, double *, double *, double *, double *, double *,
                    double *, double *, double *, double *);

void ReadSTLGeometry(void);

long ReadSTLMesh(long);

char *ReadTextFile(char *);

void ReadUMSHGeometry(void);

void ReadWWDMesh(long);

long ReallocMem(long, long);

double ReaMulti(long, long, double, long);

void ReceiveIFCInputData(void);

void RecoilDet(long, double, double, double, double, double, double, double,
               double, double, double, long);

void RecycleRMXData(long, long, long, long, long);

long ReDistributeQues(void);

void ReDistributeStacks(void);

void ReduceBuffer(void);

void ReducePrivateRes(void);

void RefreshInventory(void);

unsigned long ReInitRNG(long);

void RelaxCCResults(void);

void RelaxInterfacePower(void);

void RelaxTransmuXS(long, long);

double RelErr(long, ...);

long RemoveFlaggedItems(long, long, long, long);

void RemoveItem(long);

void RemoveStoreFiles(void);

void RemoveVoidCells(void);

void ReopenList(long);

void ReplaceItem(long, long);

void ReplacePhotonData(void);

void Reprocess(long);

void ResetOption(long, long);

void ResetOTFBurn(void);

void ResetPoisonConc(void);

void ResetTimer(long);

void ResizeDynSrc(void);

void ResizeFissionSrc(void);

double ResponseFunction(double, long);

void RIACycle(void);

void RROutput(void);

void SabScattering(long, double *, double *, double *, double *, long);

void SampleENDFLaw(long, long, double, double *, double *, long);

void SampleDelnu(void);

void SampleMeshDelnu(long, long, long);

double SampleMu(long, long, double, double [7], double *, long);

long SampleNu(double, long, long);

void SampleIFCData(long);

void SamplePlasmaSrc(long, double *, double *, double *, double *, double *,
                     double *, double *, double *, double *, long);

void SamplePointDelnu(long, long, long);

long SamplePrecursorGroup(long, double, long);

double SamplePTable(long, double, long, long);

long SampleReaction(long, long, double, double, long);

long SampleSrcPoint(long, long, long);

double SampleTabular(const double *, const double *, const double *, long,
                     long, long);

double ScalarProd(const double *, const double *);

void schurFactorization(long, complex **, complex **, complex **);

void Score(long, long, double, double, double, double, double, double, double,
           double, double, double, double, double, long);

void ScoreAlb(long, double, double, double, double, double, double, double,
              double, double, long);

void ScoreActiDet(double, long, long, double, double, double, double,
                  double, double, long);

void ScoreCapture(long, long, double, double, long);

long ScoreCMM(double, double, double, double, double, double, double, double,
              double, double, long);

void ScoreCPD(double, double, double, long);

void ScoreDF(double, double, double, double, double, double, double, double,
             double, long);

void ScoreDirectSensitivity(long, double, long, long, long, double, double,
                            double, long);

void ScoreFET(const double *const, long, long, double, double, double, double,
              double, long, long);

void ScoreFission(long, long, double, double, double, long, double,
                  double, double, double, long, long, long);

void ScoreGC(double, double, double, double, double, double, long, long,
             double, double, double, double, double, double, long);

void ScoreICMTrk(long, double, double, double, double, double, double, double,
                 double, double, long);

void ScoreICMCol(double, double, double, double, double, double, long, long,
                 double, double, double, double, double, double, long);

void ScoreInterfaceFlux(double, double, double, double, double, double, double,
                        long);

void ScoreInterfacePower(double, double, double, double, double, double, long,
                         long);

void ScoreIterNucs(double, long, double, double, double, long);

void ScoreMesh(long, long, double, double, double, double, double, double,
               double, double, double, long);

void ScoreMicroDep(double, long, double, double, long);

void ScoreOTFBurn(double, long, double, double, long);

void ScorePB(double, double, long);

void ScorePerturbations(long, long, long, long, double, double, double, long);

void ScorePhotonHeat(long, long, double, double, double, double, double,
                     double, double, long);

void ScorePinPower(long, long, double, double, double, long);

void ScorePoison(double, long, double, double, double, long);

void ScoreRMXCurr(long, double, double, double, double, double, long);

void ScoreRMXMFP(long, long, double, double, double, double, long);

void ScoreRMXResp(long, long, double, long);

void ScoreRMXSrc(long, double, double, double, double, double, long);

void ScoreScattering(long, long, double, double, double, double, double, long);

void ScoreSCAResp(long, double, double, double, double, long);

void ScoreSensLabel(long, long, long, double, double, long);

void ScoreSurf(long, double *, double *, double *, double, double, double,
               double, double, double, double, double, double, long);

void ScoreTimeConstants(double, double, long, long, long);

void ScoreTimeSource(long, long, long, double, double, double, double, double,
                     double, double, long);

void ScoreTransmuXS(double, long, double, double, long);

void ScoreUFS(double, long, double, double, double, double, double, double,
              long);

long SearchArray(const double *, double, long);

long SeekArray(const double *, long, long);

long SeekList(long, long, double, long);

long SeekListStr(long, long, char *);

void SendIFCInputTemplates(void);

void SendIFCMesh(long);

void SendIFCOutputData(void);

void SendIFCOutputTemplates(void);

void SensitivityOutput(void);

void SeparateBANuc(void);

void SetCIPopulation(void);

long SetCoefCalc(long);

void SetDecompID(void);

void SetDepStepSize(long, long);

void SetDetFlags(long, long);

void SetDirectPointers(long);

void SetFissE(void);

long SetHisvCalc(long, long);

void SetIFCTMSLimits(long, long);

void SetNormalization(long);

void SetOption(long, long);

void SetOptimization(void);

void SetPathLevels(long, long);

void SetPrecursorGroups(void);

void SetSTLMeshPointers(long);

void ShareInputData(void);

void ShuntingYard(long, long *, long);

void SignalExternal(int);

void SignalHandler(int);

void SimplifyUMSHSurfaces(long);

double **MSRReaList(double **, long, long);

void SocketReceive(char *, long);

void SocketReceiveDouble(double *);

void SocketReceiveLong(long *);

void SocketReceiveString(char *);

void SocketSend(const char *, long);

void SocketSendDouble(double);

void SocketSendLong(long);

void SocketSendString(const char *);

void SolveOTFBurn(void);

void SolveRMX(long, long);

void SolveRMXCrit(long);

void SortAll(void);

void SortArray(double *, long);

void SortList(long, long, long);

double Speed(long, double);

void SplitList(long, long, double, long);

void SrcDet(long, long, double, double, double, double, double, double,
            double, double, double, long);

double SrcRoulette(long, double, double, double, double, double, double *,
                   long);

void StartTimer(long);

long StatBin(long, ...);

double StatSum(long, ...);

void StatTests(void);

void StdComp(char *, char *);

double StdDev(long, ...);

double STLFacetDistance(long, double, double, double, double, double, double,
                        long, long);

void STLMatFinder(void);

long STLRayTest(long, long, double, double, double, double, double, double,
                long, long);

void StopAtBoundary (long *, double *, double *, double *, double *, double,
                     double, double, long);

long StopAtWWBound(long, long, double *, double *, double *, double, double,
                   double, double, double, double, double, double *, double *,
                   long *, long);

void StopCCIter(void);

void StopCI(void);

long StopSurfDetFlg(long, double *, double *, double *, double, double, double,
                    double, double, double, double, double *, double *, long *,
                    long);

void StopTimer(long);

void StoreComposition(long, double, double);

void StoreSensEvent(long, long, double, double, long);

void StoreHistoryPoint(long, long, long, double, double, double, double,
                       double, double, double, double, double, double, long, long);

long StoreRMXEvent(long, long, long, double, long, long);

void StoreTransmuXS(long, long, long, long);

void StoreValuePair(long, double, double, long);

void SumDivCompositions(void);

double SumPrivateData(long);

double SumPrivateRes(long);

void SumTotXS(long);

void SuperDet(long, double, double, double, double, double, double, double,
              double, double, double, long);

void SupervisedCycle(void);

void SurfaceNormal(long, double, double, double, double *, double *, double *,
                   long);

double SurfaceDistance(long, const double *, long, long, double, double,
                       double, double, double, double, long);

void SurfaceSrc(long, long, double *, double *, double *, double *, double *,
                double *, long);

double SurfaceVol(long);

void SwapItems(long, long);

void SwapUniverses(long, long);

struct ccsMatrix *SymbolicLU(struct ccsMatrix *);

double SymmetryBoundary(long, double, double, double, double, double, double);

void SystemStat(void);

void TargetVelocity(long, double, double *, double *, double *, double, double,
                    double, double, long);

long TestASCIIFile(char *);

void TestDOSFile(char *);

long TestFacetOverlap(long, double, double, double, double, double, double);

long TestHisvBreak();

double TestParam(char *, char *, long, char *, long, ...);

void TestSTLGeometry(void);

void TestSTLSolids(long, long, long);

long TestSurface(long, double, double, double, long, long);

long TestUniSym(double, double, double);

double TestValuePair(long, double, long);

void TestXS(void);

void TetPutBoundingBox(long, long, long[4]);

double TetraVol(long);

double *ThinGrid(double *, long *, double);

long TimeCutoff(long, long, long *, double *, double *,        double *, double *,
                double, double, double, double,        double, double *, double,
                double, long, long);

char *TimeIntervalStr(double);

char *TimeStamp(void);

char *TimeStr(long);

double TimerCPUVal(long);

double TimerVal(long);

void TmpMajorants(void);

void ToBank(long, long);

void ToCommonQue(long);

void ToLimbo(long, long);

void ToQue(long, long);

double TorusDis(double, double, double, double, double, double, double,
                double, double);

void ToStack(long, long);

void ToStore(long, long, long);

double TotXS(long, long, double, long);

void TrackFile();

void Tracking(long);

void TrackingError(long, double, long, long, long);

long TrackMode(long, long, double, double, double, long, long);

double TransportCorrection(long, double, long);

void TransportCycle(void);

complex trapz(long, double *, complex *);

double TrapzReal(const double *, const double *, long, double *, long);

void TrapzRealCum(const double *, const double *, double *, long, long);

double Truncate(double, long);

double *TTA(struct ccsMatrix *, double *, double);

double TTAChain(long, double, double *, double *);

void TTALoop (long, double, long, double *, double, double, double, double,
              double *, struct ccsMatrix *,long);

void TTByield(long, double, long, double *, long, long *, long *, long);

void TTBenergy(long mat, long part, double Ee, double x, double y, double z,
          double u, double v, double w, double wgt, double t, long mt,
          long idx, long nbr, double lEe, double *Edep, long id);

void TTB(long, long, double, double, double, double, double, double, double,
         double, double, long, double *, long);

double UFSFactor(double, double, double);

void UnionizeGrid(void);

void UniSym(long, double *, double *, double *, double *, double *, double *);

void UniverseBoundaries(void);

void UpdateCIStop(long, double *, long);

void UpdateIFCDensMax(long);

void UpdateIFCTempMinMax(long);

void UpdateIterNucs(void);

void UpdateMicroDens(void);

double UpdateRMXWgt(long, long);

double UresDiluMicroXS(long, double, long);

double UresFactor(long, double, long);

void UserIFC(long, double *, double *, double, double, double, double,
             long, const double *);

void UserIter(long);

void UserSrc(long, double *, double *, double *, double *,  double *, double *,
             double *, double *, double *, long);

void UserSurf(long, long, const double *, double *, double *, double *,
              double *, double *, double *, long *, double *, double, double,
              double, double, double, double);

long UserTransmuCut(long, long);

double ValuePairIdx(long, long);

double ValuePairVal(long, long);

double vectorNorm(long, complex *);

void VectorProd(double *, const double *, const double *);

void VirtGCUColFlags(double, double, double, long);

void VolumesMC(void);

void VRCycle(void);

void WalkerAliasInit(const double *, long, double *, double *, double *);

long WalkerAliasSample(const double *, const double *, const double *, long,
                       long);

void Warn(char *, char *, ...);

long WeightWindow(long, long, long, double, double, double, double, double,
                  double, double, double *, double, long, long);

long WhereAmI(double, double, double, double, double, double, long);

double *WorkArray(long, long, long, long);

void WriteCIMomFluxes(long);

void WriteDynSrc(void);

void WriteTetMeshtoGeo(void);

void WriteDepFile(void);

void WriteICMData(void);

void WriteSensResult(FILE *, char *, long, long, long,
                     long, long, long, long);

void WriteSourceFile(long, double, double, double, double, double, double,
                     double, double, double, double, long);

void WriteSTLMesh(long);

void WriteUMSHtoSTL(void);

void WriteWWMesh(long);

double WWDis(long, double, double, double, double, double, double);

double WWImportance(long, double, double, double, double, double, double,
                    double, long);

void WWinSrc(long, long, double *, double *, double *, double, double *, long);

void XSPlotter(void);

char *ZAItoIso(long, long);

double ZDis(double, double, double);

void ZoneCount(long, long, long);

#ifdef __cplusplus
}
#endif

/*****************************************************************************/

/***** Function prototypes for interactive plotter ***************************/

#ifdef INTERACTIVE_PLOTTER

long getGeometryPlotMatrix(long, long, long, double, double, double, double,
                           double, double, double, long*, long*);

#endif

/*****************************************************************************/

/***** Function prototypes for FINIX coupling ********************************/

#ifdef FINIX

/* Function prototypes */

void CollectFinix(void);
void CreateFinixIFC(void);
void FreeFinix(void);
void IterateFinix(void);
void PrintFinix(void);
void ProcessFinix(void);
void ReadFinixIFC(void);
void RunFinix(long, long);
void UpdateFinixIFC(void);
void UpdateFinixPower(long, long);
void WriteFinixIFC(void);
void WriteFinixInputFile(long, long);

#else

/* Replace by dummy definitions */

#define CollectFinix(void);
#define CreateFinixIFC(void);
#define DistributeFinix(void);
#define IterateFinix(void);
#define PrintFinix(void);
#define ProcessFinix(void);
#define ReadFinixIFC(void);
#define RunFinix(a);
#define UpdateFinixIFC(void);
#define UpdateFinixPower(a,b)
#define WriteFinixIFC(a);
#define WriteFinixInputFile(a,b);

#endif

/*****************************************************************************/

/***** Function prototypes for internal coupling *****************************/

#ifdef SerpentINT

/* Function prototypes */

void InitInternal(void);
void InternalIFCPoint(void);
void IterateInternal(void);
void ScoreInternaPower(void);

#else

/* Replace by dummy definitions */

#define InitInternal(void);
#define InternalIFCPoint(void);
#define IterateInternal(void);
#define ScoreInternalPower(void);

#endif

/*****************************************************************************/

/***** Functions replaced by macros in non-OpenMP mode ***********************/

#ifdef OPEN_MP

#define OMP_THREAD_NUM omp_get_thread_num()

#ifdef DEBUG

void AddPrivateData(long, double, long);
double GetPrivateData(long, long);
void PutPrivateData(long, double, long);

#else

#define AddPrivateData(a,b,c)(PRIVA[a + c*(long)RDB[DATA_REAL_PRIVA_SIZE]] += (double)b)
#define GetPrivateData(a,b)(PRIVA[a + b*(long)RDB[DATA_REAL_PRIVA_SIZE]])
#define PutPrivateData(a,b,c)(PRIVA[a + c*(long)RDB[DATA_REAL_PRIVA_SIZE]] = (double)b)

#endif

#else

#define OMP_THREAD_NUM 0

#ifdef DEBUG

void AddPrivateData(long, double, long);
double GetPrivateData(long, long);
double *PrivateDataPtr(long, long);
void PutPrivateData(long, double, long);

#else

#define AddPrivateData(a,b,c)(PRIVA[a] += (double)b)
#define GetPrivateData(a,b)(PRIVA[a])
#define PrivateDataPtr(a,b)(&PRIVA[a])
#define PutPrivateData(a,b,c)(PRIVA[a] = (double)b)

#endif

#endif


/*****************************************************************************/

/***** Functions replaced by macros in non-debug mode ************************/

#ifdef DEBUG

void CheckPointer(char *, char *, long, long);

void CheckValue(char *, char *, char *, double, double, double);

long ListPtr(long, long);

long NextItem(long);

long PrevItem(long);

double GetPrivateRes(long);

#else

#define CheckPointer(a, b, c, d)

#define CheckValue(a, b, c, d, e, f)

#define ListPtr(a, b)((long)RDB[(long)RDB[a + LIST_PTR_DIRECT] + b + 1])

#define NextItem(a)((long)RDB[a + LIST_PTR_NEXT])

#define PrevItem(a)((long)RDB[a + LIST_PTR_PREV])

#define GetPrivateRes(a)(RES2[a])

#endif

/*****************************************************************************/

/***** Additional macros *****************************************************/

#define ListCommon(a)((long)RDB[a + LIST_PTR_COMMON])
#define ListSize(a)((long)RDB[ListCommon(a) + LIST_COMMON_N_ITEMS])
#define ListRoot(a)((long)RDB[ListCommon(a) + LIST_COMMON_PTR_ROOT])
#define OMPPtr(a, b)((long)RDB[a] + b)
#define LorentzFactor(a)(a/NEUTRON_E0 + 1.0)
#define GetDDRes(a)(RES3[a])

/*****************************************************************************/

/***** Additional structures and prototypes for DD ***************************/

/* MGa: Struct to manage asynchronous sends for DD communications */

#ifdef MPI

typedef struct {
  long i; /* Index of the current buffer */
  long j; /* Index of the current position in the current buffer */
  long n; /* Total number of allocated buffers */
  long size; /* Number of buffers currently allocated (total - freed) */
  long mpiidto; /* Task to send to */
  long tag; /* Tag for the messages */
  long buffsize; /* Maximum buffer size */
  double **buffs; /* Buffers (buffs[i][j] is element j of buffer i) */
  MPI_Request **reqs; /* Requests to check for send completion (MPI_Isend(void)) */
  MPI_Status **stats; /* Statuses to get message info (MPI_Isend(void)) */
} DDSend;

/* MGa: Struct to manage asynchronous receives for DD communications */

typedef struct {
  long mpiidfrom; /* Task to receive from */
  long tag; /* Tag for the messages */
  long buffsize; /* Maximum buffer size */
  double *buff; /* Buffer */
  MPI_Request req; /* Request to check for receive completion (MPI_Irecv(void)) */
  MPI_Status stat; /* Status to get message info (MPI_Irecv(void)) */
} DDRecv;

#else

/* Use dummy types to avoid errors and compiler warnings when not */
/* compiled with MPI options */

typedef struct {
  long dummy;
} DDSend;

typedef struct {
  long dummy;
} DDRecv;

#endif

/* Function prototypes */

void CheckDDRecv(DDRecv *, long, int *, int *);
void CheckFinishDDSynch(long *);
void CheckFinishDDAsynch(long *);
void CleanUpDDComms(void);
void CleanUpDDRecv(DDRecv *, long);
void CleanUpDDSend(DDSend *, long, long);
void SendDDParticle(DDSend *, long);
void DistributeDDLimbo(void);
void FlushDDParticles(void);
void FreeDDComms(void);
void InitDDComms(void);
void InitDDRecv(DDRecv *, long, long, long);
void PostDDRecv(DDRecv *);
void PostDDSend(DDSend *);
void ReallocDDSend(DDSend *);
void ReceiveDDParticleCounters(void);
void ReceiveDDParticles(long *);
void ResetDDComms(void);
void ResetDDSend(DDSend *, long, long, long);
void SendDDParticleCounter(void);
void SetDDIDSimple(void);

/*****************************************************************************/

/***** Global arrays and variables *******************************************/
/*                                                                           */
/* Serpent 2 stores data in the following arrays:                            */
/*                                                                           */
/* WDB    Writable main data block that contains nuclear data, geometry,     */
/*        pointters and calculation parameters. Must be OpenMP protected     */
/*        when written during threaded routines.                             */
/*                                                                           */
/* RDB    Pointer to WDB array, but type-cast to constant double. This is to */
/*        induce a compiler warning when trying to write in the main data    */
/*        at the wrong place.                                                */
/*                                                                           */
/* PRIVA  Data block for storing OpenMP private values. Divided into         */
/*        segments and can be accessed by special routines without atomic or */
/*        critical pragmas.                                                  */
/*                                                                           */
/* BUF    Buffer for storing cycle/batch wise data for statistics. Divided   */
/*        into segments or accessed with atomic pragmas by special routines. */
/*                                                                           */
/* RES1   First results array, used for storing statistics. Not accessed by  */
/*        OpenMP threads.                                                    */
/*                                                                           */
/* RES2   Second results array, containing large tables of integral data for */
/*        which the statistics is not needed. Divided into segments or       */
/*        accessed with atomic pragmas by special routines.                  */
/*                                                                           */
/* RES3   Same as RES2, but used for storing material-wise data for burnup   */
/*        calculation when domain decomposition is in use. Always shared     */
/*        between OpenMP threads.                                            */
/*                                                                           */
/* ASCII  Data array for storing text strings.                               */
/*                                                                           */
/* ACE    Data array for storing cross sections and nuclide data before      */
/*        processing. Freed before transport cycle.                          */
/*                                                                           */
/*****************************************************************************/

/* Arrays */

double *ACE;
double *WDB;
double *PRIVA;
double *BUF;
double *RES1;
double *RES2;
double *RES3;

const double *RDB;

char *ASCII;

/* This size must be larger than cache line width */

#define RNG_SZ 100

unsigned long *SEED;
unsigned long *SEED0;

/*****************************************************************************/

/* Output pointers */

FILE *errp;
FILE *outp;

/* Number of mpi tasks and id */

int mpitasks;
int mpiid;

#ifdef MPI
MPI_Comm my_comm;
#endif

/* Random number seed */

unsigned long parent_seed;

/* Collision counter */

/* Timers */

struct {
  double t0;
  double t;
  double cpu_t0;
  double cpu_t;
  int on;
} timer[TOT_TIMERS + 1];

/*****************************************************************************/

/***** Variables used for particle communications in DD mode *****************/
/*                                                                           */
/* Note:                                                                     */
/*                                                                           */
/* The particle-tracking termination algorithm, i. e. the method to find out */
/* if all particle histories have been completed in DD mode is basically the */
/* one described in "An efficient, robust, domain-decomposition algorithm    */
/* for particle Monte Carlo" by T. A. Brunner and P. S. Brantley.            */
/*                                                                           */
/*****************************************************************************/

/* Parent and children in the binary tree for asynchronous termination check */

long dd_mpiid_parent, dd_mpiid_children[2];

/* Particle counters */
/* Note: dd_part_count (particles sent - particles received) is the local */
/* balance of particle transfers. When the sum of this counter over all the */
/* tasks is 0 the tracking is over. This reduction is first done */
/* asynchronously using a binary tree structure to get an efficient */
/* estimation, and then it's done with a synchronous operation to get the */
/* real instantaneous value. dd_part_count_children[i] is the counter of */
/* child task i in the binary tree. */

long dd_part_count, dd_part_count_children[2];

/* Size of the data block for each particle */

long dd_part_size;

/* Size of the particle buffers (number of particles * particle size) */

long dd_buff_size;

/* DDSend and DDRecv structs to manage particle messages */
/* Note: there's one DDSend and one DDRecv for each domain/task. */

DDSend *dd_part_sends;
DDRecv *dd_part_recvs;

/* DDSend and DDRecv structs to manage particle counter messages */
/* Note: there's one DDSend for the parent and one DDRecv for each children. */

DDSend dd_part_count_send;
DDRecv dd_part_count_recvs[2];

/* DDSend and DDRecv structs to manage synchronization request messages */
/* Note: there's one DDSend for each children and one DDRecv for the parent. */

DDSend dd_synch_req_sends[2];
DDRecv dd_synch_req_recv;

/*****************************************************************************/
