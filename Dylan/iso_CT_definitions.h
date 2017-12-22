#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 NO
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     17

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        YES
#define  RESISTIVITY             NO
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
#define  Q_factor                6
#define  Velocity_exponent       7
#define  Rotation                8
#define  Mean_mol_weight         9
#define  Magnetic_incl           10
#define  Cs_p                    11
#define  R_max                   12
#define  D_separation            13
#define  R_secondary             14
#define  Vinf_secondary          15
#define  Mdot_secondary_ratio    16

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.0e-12
#define  UNIT_LENGTH             (CONST_Rsun*9.0)
#define  UNIT_VELOCITY           1.0e+5
#define  UNIT_MASS               (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_TIME               (UNIT_LENGTH/UNIT_VELOCITY)
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

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   VANLEER_LIM
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO
