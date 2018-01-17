#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          WENO3
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     21

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  Eta                     0
#define  M_star                  1
#define  R_star                  2
#define  L_star                  3
#define  T_star                  4
#define  CAK_alpha               5
#define  CAK_delta               6
#define  CAK_ifrc                7
#define  CAK3D_nyy               8
#define  CAK3D_npp               9
#define  Q_factor                10
#define  Velocity_exponent       11
#define  Rotation                12
#define  Mean_mol_weight         13
#define  Magnetic_incl           14
#define  Cs_p                    15
#define  R_max                   16
#define  D_separation            17
#define  R_secondary             18
#define  Vinf_secondary          19
#define  Mdot_secondary_ratio    20

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.0e-12
#define  UNIT_LENGTH             (CONST_Rsun*9.0)
#define  UNIT_VELOCITY           1.0e5
#define  UNIT_TIME               (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_MASS               (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_G                  (CONST_G*UNIT_DENSITY*pow(UNIT_TIME,2))
#define  UNIT_kB                 ((CONST_kB*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5)))
#define  UNIT_B                  (sqrt((4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)))
#define  UNIT_L                  (pow(UNIT_TIME,-3)*(UNIT_DENSITY*pow(UNIT_LENGTH,5)))
#define  UNIT_c                  (CONST_c/UNIT_VELOCITY)
#define  tyear                   3.15569e+7
#define  tday                    8.64e+4
#define  L_sun                   3.846e+33
#define  VTK_VECTOR_DUMP         YES
#define  GLM_EXTENDED            YES
#define  CAK                     YES
#define  CHOMBO_LOGR             YES
#define  CHOMBO_CONS_AM          NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             VANLEER_LIM
