/* ********************************************************************* */
/*
   Initilisation file for a radiativly driven stellar wind with a 
   non-rigid dipole configuration magnetic field.

   The boundary and initial conditions are taken 
   from Runacres and Owocki (2002)    

   The method for calculating the radiative acceleration 
   comes from CAK (1975)      

   The model only works with polar corrdinates in 2D, with the 
   MHD module. 1D, 3D, HD, RHD, RMHD and other geometries do not 
   work at the moment.

*/
/* ********************************************************************* */
#include "pluto.h"                                                                  

/* ********************************************************************* */
typedef struct STAR
/*!
 * Type defineition for the star. This structure contains all the 
 * variables that are specific to the star. A structure of this 
 * type is initialised by every finction that needs stellar paramerter.
 * In this way, Only one function need calculate the parameters.
 *
 *********************************************************************** */
{

  double eta;
  double mass;
  double radius;
  double Eddington;
  double luminosity;
  double temperature;
  double alpha;
  double q_fac;
  double mass_loss;
  double sound_speed;
  double escape_velocity;
  double terminal_velocity;
  double vel_law_exponent;
  double Bfield;
  double gravity;
  double rotational_velocity;
  double mean_mol;
  double Bfield_angle;
  double surface_rho_param;

} Star;
/* ********************************************************************* */

void InitStar1(Star *star1);

void InitMagneticField(double *magnetic_field, double x1, double x2, 
                       double x3, Star star1);

void CAKAcceleration(const Data *d, RBox *box, Grid *grid, Star star1, int i, int j, int k);

double FiniteDiskCorrection(double dvdr, double vx1, double x1, double alpha);

double VelocityGradientVectorCartesian(const Data *d, RBox *box, Grid *grid, 
                                     double *x1, double *x2, double *x3,
                                     int i, int j, int k);

double AccelVectorRadial(const Data *d, Grid *grid, 
                       double f, double dvdr, Star star1, 
                       int i, int j, int k);

void AccelVectorNonRadial(double *gline, const Data *d, Grid *grid, 
                          double f, double *gradV, Star star1, 
                          int i, int j, int k);

double TwoDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, 
                            double r, double ddr);
double ThreeDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, double *x3,
                            double r, double ddr);

double VelocityStellarSurface2D(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, 
                            double *x3, double r);

double VelocityStellarSurface3D(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, 
                            double *x3, double r);

/* ********************************************************************* */
void InitStar1(Star *star1)
/*!
 * Calculate components of the velocity gradient.
 *
 * \param [in]  star1  poniter to star type data container.
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  star1->eta = g_inputParam[Eta];

  star1->mass = g_inputParam[M_star]*CONST_Msun/UNIT_MASS;

  star1->radius = g_inputParam[R_star]*CONST_Rsun/UNIT_LENGTH;

  star1->Eddington = 2.6e-5*g_inputParam[L_star]
                    *(1.0/g_inputParam[M_star]);

  star1->luminosity = g_inputParam[L_star]*L_sun/UNIT_L;

  star1->temperature = g_inputParam[T_star];

  star1->alpha = g_inputParam[CAK_alpha];

  star1->q_fac = g_inputParam[Q_factor];

  star1->mass_loss = star1->luminosity/(UNIT_c*UNIT_c)*
                      star1->alpha/(1.0 - star1->alpha)
                    *pow(star1->q_fac*star1->Eddington
                         /(1.0 - star1->Eddington), 
                         (1.0 - star1->alpha)/star1->alpha)
                    *pow(1.0 + star1->alpha, -1.0/star1->alpha);

  star1->mean_mol = g_inputParam[Mean_mol_waight];

  star1->sound_speed = sqrt(UNIT_kB*star1->temperature
    /(star1->mean_mol*(CONST_AH/UNIT_MASS)*CONST_amu));

  star1->escape_velocity = sqrt(2.0*UNIT_G*star1->mass*(1.0 - star1->Eddington));                              

  star1->terminal_velocity = star1->escape_velocity
    *sqrt(star1->alpha/(1.0 - star1->alpha));

  star1->vel_law_exponent = g_inputParam[Velocity_exponent];

  star1->Bfield = sqrt(star1->eta
                 *star1->mass_loss*UNIT_MASS/UNIT_TIME
                 *star1->terminal_velocity*UNIT_VELOCITY
                 /pow(UNIT_LENGTH, 2))/UNIT_B;

  star1->gravity = -UNIT_G*star1->mass*(1.0 - star1->Eddington);

  star1->rotational_velocity = g_inputParam[Rotation]*sqrt(UNIT_G*star1->mass);

  star1->Bfield_angle = g_inputParam[Magnetic_incl];

  star1->surface_rho_param = g_inputParam[Cs_p];

  return;
}
/* ********************************************************************* */

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*!
 * Initilise the state vector with physical values accourding to 
 * user supplied expressions.
 *
 * \param [in]  v   pointer to the state vector array.
 * \param [in]  x1  x1 position
 * \param [in]  x2  x2 position
 * \param [in]  x3  x3 position
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  double velocity, rho0;
  double r, r2;
#if PHYSICS == MHD 
#if BACKGROUND_FIELD == NO
  double magnetic_field[3];
#endif
#endif

  Star star1;
  InitStar1(&star1);

  rho0 = star1.mass_loss/(4.0*CONST_PI
           *(star1.sound_speed/star1.surface_rho_param));

  // Set the minimum pressure via the stellar temperature.
  g_smallPressure = rho0*star1.temperature/(KELVIN*star1.mean_mol);

# if COOLING !=NO
  /* Dylan,Asif: this sets minimum T of the sim.
     there is also gmaxCoolingRate that can limit timestep.
     it will be worthwhile to use it as well. */
  g_minCoolingTemp = star1.temperature;
  g_maxCoolingRate = 0.5;
# endif 

#if EOS == IDEAL
  g_gamma = 1.05;
#endif

#if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = star1.sound_speed;
#endif

#if ROTATING_FRAME == YES                                                          
  g_OmegaZ = star1.rotational_velocity;
#endif

  r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  r  = sqrt(r2);

  if (r > star1.radius) {

    velocity = star1.terminal_velocity
               *pow(1.0 - star1.radius/r, star1.vel_law_exponent);

    v[RHO] = (star1.mass_loss/(4.0*CONST_PI*velocity*r2));

    EXPAND(v[VX1] = velocity*x1/r;,
           v[VX2] = velocity*x2/r;,                               
           v[VX3] = velocity*x3/r;)

    v[TRC] = 0.0;

  } else if (r > 0.5*star1.radius && r <= star1.radius) {

    v[RHO] = (star1.mass_loss
               /(4.0*CONST_PI*(star1.sound_speed
               /star1.surface_rho_param)));

    velocity = star1.terminal_velocity
                *pow(1.0 - (1.0/1.001), star1.vel_law_exponent);

    EXPAND(v[VX1] = velocity*x1/r;,
           v[VX2] = velocity*x2/r;,                               
           v[VX3] = velocity*x3/r;)

    v[TRC] = 0.0;

  } else if (r <= 0.5*star1.radius) {

    v[RHO] = (star1.mass_loss
               /(4.0*CONST_PI*(star1.sound_speed
               /star1.surface_rho_param)));

    EXPAND(v[VX1] = 0.0;,                                                 
           v[VX2] = 0.0;,                               
           v[VX3] = 0.0;)

    v[TRC] = 0.0;

  }

#if EOS == IDEAL
  v[PRS] = (v[RHO]*star1.temperature/(KELVIN*star1.mean_mol));
#endif                                                               

#if PHYSICS == MHD                                   
#if BACKGROUND_FIELD == NO
  InitMagneticField(magnetic_field, x1, x2, x3, star1);
  EXPAND(v[BX1] = magnetic_field[0];,
         v[BX2] = magnetic_field[1];,
         v[BX3] = magnetic_field[2];)
#endif
#if BACKGROUND_FIELD == YES
  EXPAND(v[BX1] = 0.0;,
         v[BX2] = 0.0;,
         v[BX3] = 0.0;)
  EXPAND(v[AX1] = 0.0;,
         v[AX2] = 0.0;,
         v[AX3] = 0.0;)
#endif
#endif


  return;
}                                                                          
/* ********************************************************************* */

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 *
 *
 * \param [in]  d     pointer to the main PLUTO data structure.
 * \param [in]  grid  pointer to an array of Grid structures.
 *
 * \return  This function has no return value.
 *
 * TODO None
 *********************************************************************** */
{
  return;
}
/* ********************************************************************* */

#if BACKGROUND_FIELD == YES
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 *
 *
 * \param [in]  d     pointer to the main PLUTO data structure.
 * \param [in]  grid  pointer to an array of Grid structures.
 *
 * \return  This function has no return value.
 *
 * TODO None
 *********************************************************************** */
{
  double magnetic_field[3], r, r2;

  Star star1;
  InitStar1(&star1);

  InitMagneticField(magnetic_field, x1, x2, x3, star1);
  EXPAND(B0[0] = magnetic_field[0];,
         B0[1] = magnetic_field[1];,
         B0[2] = magnetic_field[2];)

  return;
}
/* ********************************************************************* */
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*!
 *
 *
 * \param [in]  d     pointer to the main PLUTO data structure.
 * \param [in]  grid  pointer to an array of Grid structures.
 *
 * \return  This function has no return value.
 *
 * TODO Impliment dynamic boundary conditions.
 *********************************************************************** */
{        

  int i, j, k;

  double r, r2;
  //double dr2, dr, ddr;
  double velocity, temp;

#if PHYSICS == MHD 
#if BACKGROUND_FIELD == NO
  double magnetic_field[3];
#endif
#endif

  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;
  //double *dx1 = grid[IDIR].dx;
  //double *dx2 = grid[JDIR].dx;
  //double *dx3 = grid[KDIR].dx;
  double ***gLx1 = d->gL[0];
  double ***gLx2 = d->gL[1];
  double ***gLx3 = d->gL[2];;

  Star star1;
  InitStar1(&star1);

  if(g_stepNumber < 2){
    double Bcgs = star1.Bfield*UNIT_B;
    double M_dot_cgs = star1.mass_loss*UNIT_MASS/UNIT_TIME;
    printf("Bcgs=%e, M_dotcgs=%e, Edd_gam=%e , Omega=%e  \n", 
           Bcgs, M_dot_cgs/6.35e25, star1.Eddington, 
           star1.rotational_velocity*UNIT_VELOCITY);
  }


  if(side == 0){
    DOM_LOOP(k,j,i){

      r2  = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
      r   = sqrt(r2);

      if (r <= 0.5*star1.radius) {

        d->Vc[RHO][k][j][i] = star1.mass_loss
                 /(4.0*CONST_PI*(star1.sound_speed
                 /star1.surface_rho_param));

#if EOS == IDEAL                                                              
        d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*star1.temperature
                                /(KELVIN*star1.mean_mol);
#endif

        D_EXPAND(d->Vc[VX1][k][j][i] = 0.0;,
                 d->Vc[VX2][k][j][i] = 0.0;,                               
                 d->Vc[VX3][k][j][i] = 0.0;)

#if PHYSICS == MHD 
#if BACKGROUND_FIELD == NO
        InitMagneticField(magnetic_field, x1[i], x2[j], x3[k], star1);
        EXPAND(d->Vc[BX1][k][j][i] = magnetic_field[0];,
               d->Vc[BX2][k][j][i] = magnetic_field[1];,
               d->Vc[BX3][k][j][i] = magnetic_field[2];)
#endif
#endif

        EXPAND(gLx1[k][j][i] = 0.0;,
               gLx2[k][j][i] = 0.0;,
               gLx3[k][j][i] = 0.0;)


        d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

      } else if (r > 0.5*star1.radius && r <= star1.radius) {
      
        d->Vc[RHO][k][j][i] = (star1.mass_loss
                                /(4.0*CONST_PI*(star1.sound_speed
                                /star1.surface_rho_param)));

#if EOS == IDEAL                                                              
        d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i]*star1.temperature
                                /(KELVIN*star1.mean_mol);
#endif

/*
        dr2 = EXPAND(dx1[i]*dx1[i], + dx2[j]*dx2[j], + dx3[k]*dx3[k]);
        dr = sqrt(dr2);
        ddr = 0.9*dr;
        if (r < star1.radius && r + 5.0*ddr > star1.radius ) {
#if DIMENSIONS == 2
          velocity = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, dr);
#endif
#if DIMENSIONS == 3
          velocity = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, dr);
#endif
        } else {
          velocity = 0.0;
        }
*/

        velocity = star1.sound_speed/star1.surface_rho_param;
        D_EXPAND(d->Vc[VX1][k][j][i] = velocity*x1[i]/r;,
                 d->Vc[VX2][k][j][i] = velocity*x2[j]/r;,
                 d->Vc[VX3][k][j][i] = velocity*x3[k]/r;)

#if PHYSICS == MHD 
#if BACKGROUND_FIELD == NO      
        InitMagneticField(magnetic_field, x1[i], x2[j], x3[k], star1);
        EXPAND(d->Vc[BX1][k][j][i] = magnetic_field[0];,
               d->Vc[BX2][k][j][i] = magnetic_field[1];,
               d->Vc[BX3][k][j][i] = magnetic_field[2];)
#endif
#endif

        d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

        EXPAND(gLx1[k][j][i] = 0.0;,
               gLx2[k][j][i] = 0.0;,
               gLx3[k][j][i] = 0.0;)

      }

#if EOS == IDEAL
      temp = KELVIN*d->Vc[PRS][k][j][i]*star1.mean_mol/d->Vc[RHO][k][j][i];

      //int mark = 0;
      if (temp < star1.temperature) {
        //printf("old_temp=%e, prs=%e, t_* prs=%e \n", temp, d->Vc[PRS][k][j][i], star1.temperature*d->Vc[RHO][k][j][i]/(KELVIN*star1.mean_mol));
        d->Vc[PRS][k][j][i] = star1.temperature*d->Vc[RHO][k][j][i]
                                /(KELVIN*star1.mean_mol);
        //printf("new_temp=%e \n", KELVIN*d->Vc[PRS][k][j][i]*star1.mean_mol/d->Vc[RHO][k][j][i]);
        //mark = 1;
        //printf("\n");
      }

      // Need to set the temp floor before calculating the CAK accel.
      //if (d->Vc[PRS][k][j][i] < (d->Vc[RHO][k][j][i])*star1.temperature
      //                           /(KELVIN*star1.mean_mol)){
      //  d->Vc[PRS][k][j][i] = (d->Vc[RHO][k][j][i])*star1.temperature
      //                           /(KELVIN*star1.mean_mol);
      //}
#endif
#if CAK == YES
      CAKAcceleration(d, box, grid, star1, i, j, k);
#endif

/*      temp = KELVIN*d->Vc[PRS][k][j][i]*star1.mean_mol/d->Vc[RHO][k][j][i];

      if (mark == 1) {
        printf("second old_temp=%e, prs=%e, t_* prs=%e \n", temp, d->Vc[PRS][k][j][i], star1.temperature*d->Vc[RHO][k][j][i]/(KELVIN*star1.mean_mol));
        d->Vc[PRS][k][j][i] = star1.temperature*d->Vc[RHO][k][j][i]
                                /(KELVIN*star1.mean_mol);
        printf("new_temp=%e \n", KELVIN*d->Vc[PRS][k][j][i]*star1.mean_mol/d->Vc[RHO][k][j][i]);
        printf("\n");
      }*/

    }
  }
  return;
}                                                                          
/* ********************************************************************* */

#if BODY_FORCE != NO
#if CAK == YES
/* ********************************************************************* */
void BodyForceVector(double *v, double *gla, double *g, double x1, double x2, double x3)
/*!
 *
 *
 * \param [in]       x1              x1 position
 * \param [in]       x2              x2 position
 * \param [in]       x3              x3 position
 * \param [in]       star1           poniter to star type data container.
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  Star star1;
  InitStar1(&star1);

  double r2, r, Fin_x1, Fin_x2, gg, g_in;

  /* - Distance from star - */
  r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r2);

  /* - Gravity outside bodies - */
  gg = star1.gravity/r/r;

  /* - Gravity inside bodies - */
  g_in = -(4.0/3.0)*CONST_PI*UNIT_G*v[RHO];     

  /* - Coriolis and centrifugal forces - */
  #if DIMENSIONS == 2
  Fin_x1 = pow(star1.rotational_velocity, 2)*x1;
  Fin_x2 = 0.0;
  #endif
  #if DIMENSIONS == 3
  Fin_x1 = pow(star1.rotational_velocity, 2)*x1 
             + 2.0*star1.rotational_velocity*v[VX2];
  Fin_x2 = pow(star1.rotational_velocity, 2)*x2 
             - 2.0*star1.rotational_velocity*v[VX1];
  #endif

  //double gl = sqrt(gla[0]*gla[0] + gla[1]*gla[1] + gla[2]*gla[2]);
  //printf("gl=%e, gg=%e, gl-gg=%e r=%e \n", gl, gg, gl-fabs(gg), r);

  if (r >= star1.radius){ // - External gravity + centrifugal + coriolis 
    g[IDIR] = gg*x1/r + Fin_x1 + gla[0];
    g[JDIR] = gg*x2/r + Fin_x2 + gla[1];
    g[KDIR] = gg*x3/r + gla[2];
  } else { // - Star interal gravity + rotation 
    g[IDIR] = g_in*x1 + Fin_x1;
    g[JDIR] = g_in*x2 + Fin_x2;
    g[KDIR] = g_in*x3;
  }

  return;
}
/* ********************************************************************* */
#endif

#if CAK == NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 *
 *
 * \param [in]       x1              x1 position
 * \param [in]       x2              x2 position
 * \param [in]       x3              x3 position
 * \param [in]       star1           poniter to star type data container.
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{
  Star star1;
  InitStar1(&star1);

  double r2, r, Fin_x1, Fin_x2, gg, g_in;

  /* - Distance from star - */
  r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r2);

  /* - Gravity outside bodies - */
  gg = star1.gravity/r/r;

  /* - Gravity inside bodies - */
  g_in = -(4.0/3.0)*CONST_PI*UNIT_G*v[RHO];     

  /* - Coriolis and centrifugal forces - */
  #if DIMENSIONS == 2
  Fin_x1 = pow(star1.rotational_velocity, 2)*x1;
  Fin_x2 = 0.0;
  #endif
  #if DIMENSIONS == 3
  Fin_x1 = pow(star1.rotational_velocity, 2)*x1 
             + 2.0*star1.rotational_velocity*v[VX2];
  Fin_x2 = pow(star1.rotational_velocity, 2)*x2 
             - 2.0*star1.rotational_velocity*v[VX1];
  #endif

  if (r >= star1.radius){ // - External gravity + centrifugal + coriolis 
    g[IDIR] = gg*x1/r + Fin_x1;
    g[JDIR] = gg*x2/r + Fin_x2;
    g[KDIR] = gg*x3/r;
  } else { // - Star interal gravity + rotation 
    g[IDIR] = g_in*x1 + Fin_x1;
    g[JDIR] = g_in*x2 + Fin_x2;
    g[KDIR] = g_in*x3;
  }
  return;
}
/* ********************************************************************* */
#endif
#endif

/* ********************************************************************* */
void CAKAcceleration(const Data *d, RBox *box, Grid *grid, Star star1, int i, int j, int k)
/*!
 * Calculate CAK line force acceleration.
 *
 * \param [in]  d      pointer to the main PLUTO data structure
 * \param [in]  grid   pointer to an array of Grid structures.
 * \param [in]  star1  poniter to star type data container.
 * \param [in]  i      x1 direction index
 * \param [in]  j      x2 direction index
 * \param [in]  k      x3 direction index
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  double f=0.0;
  double dvdr, gline;

  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;
  double r=0.0, r2=0.0, vr=0.0;

  r2  = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
  r   = sqrt(r2); 

  dvdr = VelocityGradientVectorCartesian(d, box, grid, x1, x2, x3, i, j, k);

  vr = EXPAND(d->Vc[VX1][k][j][i]*x1[i]/r, 
            + d->Vc[VX2][k][j][i]*x2[j]/r, 
            + d->Vc[VX3][k][j][i]*x3[k]/r);

  f = FiniteDiskCorrection(dvdr, vr, r, star1.alpha);

  gline = AccelVectorRadial(d, grid, f, dvdr, star1, i, j, k);

  if (r < star1.radius) {
    EXPAND(d->gL[0][k][j][i] = 0.0;,
           d->gL[1][k][j][i] = 0.0;,
           d->gL[2][k][j][i] = 0.0;)
  }else {
    EXPAND(d->gL[0][k][j][i] = gline*x1[i]/r;,
           d->gL[1][k][j][i] = gline*x2[j]/r;,
           d->gL[2][k][j][i] = gline*x3[k]/r;)
  }

  return;
}
/* ********************************************************************* */

/* ********************************************************************* */
double VelocityGradientVectorCartesian(const Data *d, RBox *box, Grid *grid, double *x1, 
                                     double *x2, double *x3,
                                     int i, int j, int k)
/*!
 * Calculate components of the velocity gradient.
 *
 * \param [in]       d      pointer to the main PLUTO data structure
 * \param [in]       x1     x1 direction array
 * \param [in]       x2     x2 direction array
 * \param [in]       x3     x3 direction array
 * \param [in]       i      x1 direction index
 * \param [in]       j      x2 direction index
 * \param [in]       k      x3 direction index
 * \param [in, out]  dvdr   array containing velocity gradient components
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  double *dx1 = grid[IDIR].dx;
  double *dx2 = grid[JDIR].dx;
  double *dx3 = grid[KDIR].dx;

  //double dx=0.0, dy=0.0, dz=0.0;
  double dr2=0.0, dr=0.0, ddr=0.0;
  double vrI[2], vr=0.0, dvdr=0.0;
  double r=0.0, r2=0.0;
  double P00, P11, P22, P01, P12; 

  r2  = EXPAND(x1[i]*x1[i],+x2[j]*x2[j],+x3[k]*x3[k]);
  r   = sqrt(r2);

  //dx = x1[i+1] - x1[i];
  //dy = x2[j+1] - x2[j];
  //dz = x3[k+1] - x3[k];
  //dr2  = EXPAND(dx*dx, + dy*dy, + dz*dz);
  dr2  = EXPAND(dx1[i]*dx1[i], + dx2[j]*dx2[j], + dx3[k]*dx3[k]);
  dr   = sqrt(dr2);
  ddr = 0.9*dr;

  // Get interpolated velocity values.
  #if DIMENSIONS == 2
  vrI[0] = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, -ddr);
  vrI[1] = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
  #endif
  #if DIMENSIONS == 3
  vrI[0] = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, -ddr);
  vrI[1] = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
  #endif

  vr = EXPAND(d->Vc[VX1][k][j][i]*x1[i]/r, 
            + d->Vc[VX2][k][j][i]*x2[j]/r, 
            + d->Vc[VX3][k][j][i]*x3[k]/r);
  P00 = vrI[0];
  P11 = vr;
  P22 = vrI[1];
  P01 = (0.5*ddr*P00 + 0.5*ddr*P11)/ddr;
  P12 = (0.5*ddr*P11 + 0.5*ddr*P22)/ddr;

  dvdr = fabs((1.0/12.0)*P00 - (2.0/3.0)*P01 
           + (2.0/3.0)*P12 - (1.0/12.0)*P22)/(0.5*ddr);
  dvdr = fmax(dvdr, 1.0e-8);

  return dvdr;
}
/* ********************************************************************* */

/* ********************************************************************* */
double FiniteDiskCorrection(double dvdr, double vx1, double r, double alpha)
/*!
 * Calculate finite disk correction factor.
 *
 * \param [in]  dvdr   velocity gradient
 * \param [in]  vx1    velocity at the current cell
 * \param [in]  r      radius
 * \param [in]  alpha  CAK force multiplyer
 *
 * \return f
 *
 * \TODO None
 *********************************************************************** */
{

  double beta_op=0.0, opa=0.0, oma=0.0, f=0.0;

  opa = 1.0 + alpha;
  oma = 1.0 - alpha;

  if (fabs(dvdr) > 1.0e-8){

    beta_op = (1.0 - vx1/(dvdr*r)) * pow(1.0/r, 2);

    if (beta_op >= 1.0){
      f = 1.0/opa;
    }else if(beta_op < -1.0e10){
      f = pow(-beta_op, alpha)/opa;
    }else if(fabs(beta_op) > 1.0e-3){
      f = (1.0 - pow(1.0 - beta_op, opa))/(beta_op * opa);
    }else{
      f = 1.0 - 0.5*alpha*beta_op*(1.0 + 1.0/3.0*oma*beta_op);
    }

  }else{
    f = 1.0;
  }

  return f;
}
/* ********************************************************************* */

/* ********************************************************************* */
double AccelVectorRadial(const Data *d, Grid *grid, 
                       double f, double dvdr, Star star1, 
                       int i, int j, int k)
/*!
 * Calculate components of the velocity gradient.
 *
 * \param [in, out]  gline  array containing acceleration components
 * \param [in]       f      finite disk correction factor
 * \param [in]       A      prarmeter 
 * \param [in]       B      B factor
 * \param [in]       x1     position in radius
 * \param [in]       temp   temperature at x1
 * \param [in]       star1  poniter to star type data container
 *
 * \return void
 *
 * \TODO 
 *********************************************************************** */
{

  double B=0.0, A=0.0, ke=0.0, temp=0.0;
  double *x1 = grid[IDIR].x;
  double *x2 = grid[JDIR].x;
  double *x3 = grid[KDIR].x;
  double r=0.0, r2=0.0;
  double gline=0.0, factor=0.0;

  ke = 4.0*CONST_PI*UNIT_G*star1.mass*UNIT_c
         *star1.Eddington/star1.luminosity;

  B = d->Vc[RHO][k][j][i]*star1.q_fac*UNIT_c*ke;

  A = 1.0/(1.0 - star1.alpha)*ke*star1.luminosity
        *star1.q_fac/(4.0*CONST_PI*UNIT_c);

  r2  = EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k]*x3[k]);
  r   = sqrt(r2);
     
  gline = f*A*pow(r, -2)*pow(dvdr/B, star1.alpha);

#if EOS == IDEAL
  // Accounting for total ionisation at high temp 
  // and for recombination at low temp.
  temp = d->Vc[PRS][k][j][i]*KELVIN*star1.mean_mol/d->Vc[RHO][k][j][i];
  factor = exp(-4.0*log(2.0)*pow((2.0 
             - temp/star1.temperature - star1.temperature/temp), 2));
  factor = fmax(factor, 1.0e-8);
  //gline *= factor;
#endif

  /* This set of if statements checks and excludes 
     the ghost zones from being accelerated. Basically
     stops strangeness from happening. */
/*
#if DIMENSIONS == 3
  if (i < 2 || j < 2 || k < 2) {
    gline = 0.0;
  }
  if (i > NX1+1 || j > NX2+1 || k > NX3+1){
    gline = 0.0;
  }
#endif
#if DIMENSIONS == 2
  if (i < 2 || j < 2) {
    gline = 0.0;
  }
  if (i > NX1+1 || j > NX2+1){
    gline = 0.0;
  }
#endif
*/
  //print("gline=%e \n", gline);

  return gline;
}
/* ********************************************************************* */

/* ********************************************************************* */
void AccelVectorNonRadial(double *gline, const Data *d, Grid *grid, 
                          double f, double *gradV, Star star1, 
                          int i, int j, int k)
/*!
 * Calculate components of the velocity gradient.
 *
 * \param [in, out]  gline  array containing acceleration components
 * \param [in]       f      finite disk correction factor
 * \param [in]       A      prarmeter 
 * \param [in]       B      B factor
 * \param [in]       x1     position in radius
 * \param [in]       temp   temperature at x1
 * \param [in]       star1  poniter to star type data container
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  double B, A, ke, temp;
  double gradV_mag, g_mag, unit_vec[3];
  double *x1 = grid[IDIR].x;
  double *x2 = grid[JDIR].x;
  double *x3 = grid[KDIR].x;
  double r, r2;

  // Unit vector in the direction of the velocity gradient.
  gradV_mag = sqrt(gradV[0]*gradV[0] 
                + gradV[1]*gradV[1] 
                + gradV[2]*gradV[2]);
  unit_vec[0] = gradV[0]/gradV_mag; 
  unit_vec[1] = gradV[1]/gradV_mag; 
  unit_vec[2] = gradV[2]/gradV_mag;

  ke = 4.0*CONST_PI*UNIT_G*star1.mass*UNIT_c
         *star1.Eddington/star1.luminosity;

  B = d->Vc[RHO][k][j][i]*star1.q_fac*UNIT_c*ke;

  A = 1.0/(1.0 - star1.alpha)*ke*star1.luminosity
        *star1.q_fac/(4.0*CONST_PI*UNIT_c);

  r2  = EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k]*x3[k]);
  r   = sqrt(r2);

  g_mag = f*A*pow(r, -2)*pow(gradV_mag/B, star1.alpha);

#if EOS == IDEAL
  // Accounting for total ionisation at high temp 
  // and for recombination at low temp.
  temp = d->Vc[PRS][k][j][i]*KELVIN*star1.mean_mol/d->Vc[RHO][k][j][i];
  g_mag *= exp(-4.0*log(2.0)*pow((2.0 - temp/star1.temperature 
             - star1.temperature/temp), 2));
#endif

  gline[0] = g_mag*unit_vec[0];
  gline[1] = g_mag*unit_vec[1];
  gline[2] = g_mag*unit_vec[2];

  return;
}
/* ********************************************************************* */



/* ********************************************************************* */
void InitMagneticField(double *magnetic_field, 
                       double x1, double x2, double x3, 
                       Star star1)
/*!
 * Calculate components of the velocity gradient.
 *
 * \param [in, out]  magnetic_field  array containing 
 *                                   magnetic field components
 * \param [in]       x1              x1 position
 * \param [in]       x2              x2 position
 * \param [in]       x3              x3 position
 * \param [in]       star1           poniter to star type data container
 *
 * \return void
 *
 * \TODO Move the magnetic field expressions to their own function.
 *********************************************************************** */
{

  double x, y, z, xp, yp, zp, r, r2;
  double br, btheta, bphi, bx, by, bz, bxp, byp, bzp, rp, rp2;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;

  star1.Bfield_angle *= 0.0174532925;

  r2 = EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r2);


#if DIMENSIONS == 2	
  // Rotate coordiantes.
  xp = x1*cos(star1.Bfield_angle) - x2*sin(star1.Bfield_angle);
  yp = x1*sin(star1.Bfield_angle) + x2*cos(star1.Bfield_angle);
  zp = 0.0;
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);
  // Calculate b-field components in rotated frame.
  if (r <= 0.5*star1.radius) {
    bx = 0.0;
    by = 16.0*star1.Bfield;
    bz = 0.0;
  } else {
    bx = 3.0*xp*yp*star1.Bfield*pow(rp,-5);
    by = (3.0*pow(yp,2)-pow(rp,2))*star1.Bfield*pow(rp,-5);
    bz = 0.0;
  }
  // Rotate B-field vector componets.  
  bxp = bx*cos(star1.Bfield_angle) + by*sin(star1.Bfield_angle);
  byp = -bx*sin(star1.Bfield_angle) + by*cos(star1.Bfield_angle);
  bzp = 0.0;
#endif
#if DIMENSIONS == 3
  // Rotate coordiantes.
  xp = x1*cos(star1.Bfield_angle) - x3*sin(star1.Bfield_angle);
  yp = x2;
  zp = x1*sin(star1.Bfield_angle) + x3*cos(star1.Bfield_angle);
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);
  // Calculate b-field components in rotated frame.
  if (r <= 0.5*star1.radius) {
    bx = 0.0;
    by = 0.0;
    bz = 16.0*star1.Bfield;
  } else {
    bx = 3.0*xp*zp*star1.Bfield*pow(rp,-5);
    by = 3.0*yp*zp*star1.Bfield*pow(rp,-5);
    bz = (3.0*pow(zp,2)-pow(rp,2))*star1.Bfield*pow(rp,-5);
  }
  // Rotate B-field vector componets.  
  bxp = bx*cos(star1.Bfield_angle) + bz*sin(star1.Bfield_angle);
  byp = by;
  bzp = -bx*sin(star1.Bfield_angle) + bz*cos(star1.Bfield_angle);
#endif

  magnetic_field[0] = bxp;
  magnetic_field[1] = byp; 
  magnetic_field[2] = bzp;

  return;
}
/* ********************************************************************* */


/* ********************************************************************* */
double TwoDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, double *x1, double *x2, 
                            double r, double ddr)
{
  /*
   *       _____________________________
   *  j+1 |              |              | 
   *      |              |              |
   *      |              |     *        |
   *      |              |              |
   *      |              |              |
   *    j |______________|______________|
   *      |              |              | 
   *      |              |              |
   *      |              |              |
   *      |              |              |
   *      |              |              |
   *  j-1 |______________|______________|
   *  
   *    i - 1          i            i + 1
   *  
   *  
   *  yb 3               4
   *      |```|`````````|
   *      |   |         |
   *   yI |___|_________|
   *      |   |         |
   *      |   |         |
   *      |___|_________|
   *  ya 1    xI         2
   *      xa            xb
   *
   * The interpolation points are always between xa, xb and ya, yb.
   *       
   * ddr is the distance forwards (backwards) from the point 
   * that the velocity gradient is needed to give the place 
   * to perform the interpolation.
   *
   *  r -> r±ddr
   *
   * |'''''''''''''''''''''''''|
   * |                         |
   * |                         |
   * |      * r+ddr            |
   * |     /|                  |
   * |    / |                  |
   * |   /  | r+ddr*cos(theta) |
   * |  /   |                  |
   * | /    |                  |
   * |/_____|__________________|
   * r     r+ddr*sin(theta)
   *
   * xI and yI are the interpolation componets. 
   * 
   */

  int u=0, s=0;
  double Ntot=0.0; // total volume of interpolation "space".
  double vrI=0.0; // final interpolated radial-velocity. 
  double vI[2]; // interpolated velocity components.
  double N[4]; // Areas used to waight nearest neighbours (NN).
  double V[2][4]; // Array to hold velocity componants at NN.
  double xa=0.0, xb=0.0; // bracketing x values.
  double ya=0.0, yb=0.0; // bracketing y values.
  double theta = atan2(x1[i],x2[j]); // convert to polar from Cartesian.
  double xI = (r+ddr)*sin(theta);
  double yI = (r+ddr)*cos(theta);

  /*
   * Eath of the following if statments checks which quadrent 
   * the interpolation point is in and gets the components 
   * of the velocity and the bracketing x and y values.
   */
  if (xI > x1[i] && yI > x2[j]){
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j]; yb = x2[j+1];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i];
      V[u-1][1] = d->Vc[u][k][j][i+1];
      V[u-1][2] = d->Vc[u][k][j+1][i];
      V[u-1][3] = d->Vc[u][k][j+1][i+1];      
    }
  } else if (xI < x1[i] && yI > x2[j]){
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j]; yb = x2[j+1];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i-1];
      V[u-1][1] = d->Vc[u][k][j][i];
      V[u-1][2] = d->Vc[u][k][j+1][i-1];
      V[u-1][3] = d->Vc[u][k][j+1][i];
    }
  } else if (xI < x1[i] && yI < x2[j]){
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j-1]; yb = x2[j];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i-1];
      V[u-1][1] = d->Vc[u][k][j-1][i];
      V[u-1][2] = d->Vc[u][k][j][i-1];
      V[u-1][3] = d->Vc[u][k][j][i];
    }
  } else if (xI > x1[i] && yI < x2[j]){
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j-1]; yb = x2[j];
    for (u = VX1; u < VX2+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i];
      V[u-1][1] = d->Vc[u][k][j-1][i+1];
      V[u-1][2] = d->Vc[u][k][j][i];
      V[u-1][3] = d->Vc[u][k][j][i+1];
    }
  }

  // Find total volume.
  N[0] = (xb - xI)*(yb - yI);
  N[1] = (xI - xa)*(yb - yI);
  N[2] = (xb - xI)*(yI - ya);
  N[3] = (xI - xa)*(yI - ya);
  Ntot = N[0] + N[1] + N[2] + N[3];
  // Normalise volumes by total.
  N[0] /= Ntot; 
  N[1] /= Ntot; 
  N[2] /= Ntot; 
  N[3] /= Ntot;

  // ==========================================
  // vI contains the interpolated velovities 
  // 
  // vI[0] = x velocity 
  // vI[1] = y "
  //
  // V contains the velocity componants at 
  // each of the sample points.
  // 
  // V[0, :] = x velocity at each sample point.
  // V[1, :] = y velocity at each sample point.
  //
  // N contains the waighting volumes for the 
  // the corner velocities.
  // ==========================================

  vI[0] = 0.0; vI[1] = 0.0;
  for (s = 0; s < 4; s++) {
    vI[0] += V[0][s]*N[s];
    vI[1] += V[1][s]*N[s];
  }
  vrI = (vI[0]*x1[i] + vI[1]*x2[j])/r;
  return vrI;

}
/* ********************************************************************* */


/* ********************************************************************* */
double ThreeDimensionalInterp(const Data *d, RBox *box, Grid *grid, 
                            int i, int j, int k, 
                            double *x1, double *x2, double *x3, 
                            double r, double ddr)
/*!
 * Calculate the velocity at the interpolation point in 3D..
 *
 * \param [out]  vI  Interpolated velocity.
 *
 * return vI 
 *
 * \TODO incoperate ddr calculation in each if statment so that the 
 *       stretched grids can be done properly. 
 *
 *       For AMR, dx1 = grid[IDIR].dx etc.. can be used.
 *********************************************************************** */
{
  int u=0, s=0, tag_if=-1, miss_tag=0;
  double Ntot=0.0;
  double vrI=0.0;
  double vI[3];
  double N[8];
  double V[3][8];
  double xa=0.0, xb=0.0;
  double ya=0.0, yb=0.0;
  double za=0.0, zb=0.0;
  double phi = atan2(x2[j],x1[i]);
  double theta = acos(x3[k]/r);
  double xI = (r+ddr)*sin(theta)*cos(phi);
  double yI = (r+ddr)*sin(theta)*sin(phi);
  double zI = (r+ddr)*cos(theta);

  if (xI > x1[i] && yI > x2[j] && zI > x3[k]){
    tag_if = 0;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i];
      V[u-1][1] = d->Vc[u][k][j][i+1];
      V[u-1][2] = d->Vc[u][k][j+1][i];
      V[u-1][3] = d->Vc[u][k][j+1][i+1];
      V[u-1][4] = d->Vc[u][k+1][j][i];
      V[u-1][5] = d->Vc[u][k+1][j][i+1];
      V[u-1][6] = d->Vc[u][k+1][j+1][i];
      V[u-1][7] = d->Vc[u][k+1][j+1][i+1];
    }
  }
  else if (xI < x1[i] && yI > x2[j] && zI > x3[k]){
    tag_if = 1;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j][i-1];
      V[u-1][1] = d->Vc[u][k][j][i];
      V[u-1][2] = d->Vc[u][k][j+1][i-1];
      V[u-1][3] = d->Vc[u][k][j+1][i];
      V[u-1][4] = d->Vc[u][k+1][j][i-1];
      V[u-1][5] = d->Vc[u][k+1][j][i];
      V[u-1][6] = d->Vc[u][k+1][j+1][i-1];
      V[u-1][7] = d->Vc[u][k+1][j+1][i];
    }
  }
  else if (xI > x1[i] && yI < x2[j] && zI > x3[k]){
    tag_if = 2;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i];
      V[u-1][1] = d->Vc[u][k][j-1][i+1];
      V[u-1][2] = d->Vc[u][k][j][i];
      V[u-1][3] = d->Vc[u][k][j][i+1];
      V[u-1][4] = d->Vc[u][k+1][j-1][i];
      V[u-1][5] = d->Vc[u][k+1][j-1][i+1];
      V[u-1][6] = d->Vc[u][k+1][j][i];
      V[u-1][7] = d->Vc[u][k+1][j][i+1];
    }
  }
  else if (xI < x1[i] && yI < x2[j] && zI > x3[k]){
    tag_if = 3;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k]; zb = x3[k+1];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k][j-1][i-1];
      V[u-1][1] = d->Vc[u][k][j-1][i];
      V[u-1][2] = d->Vc[u][k][j][i-1];
      V[u-1][3] = d->Vc[u][k][j][i];
      V[u-1][4] = d->Vc[u][k+1][j-1][i-1];
      V[u-1][5] = d->Vc[u][k+1][j-1][i];
      V[u-1][6] = d->Vc[u][k+1][j][i-1];
      V[u-1][7] = d->Vc[u][k+1][j][i];
    }
  }

  // zI < zP 
  else if (xI > x1[i] && yI > x2[j] && zI < x3[k]){
    tag_if = 4;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j][i];
      V[u-1][1] = d->Vc[u][k-1][j][i+1];
      V[u-1][2] = d->Vc[u][k-1][j+1][i];
      V[u-1][3] = d->Vc[u][k-1][j+1][i+1];
      V[u-1][4] = d->Vc[u][k][j][i];
      V[u-1][5] = d->Vc[u][k][j][i+1];
      V[u-1][6] = d->Vc[u][k][j+1][i];
      V[u-1][7] = d->Vc[u][k][j+1][i+1];
    }
  }
  else if (xI < x1[i] && yI > x2[j] && zI < x3[k]){
    tag_if = 5;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j]; yb = x2[j+1];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j][i-1];
      V[u-1][1] = d->Vc[u][k-1][j][i];
      V[u-1][2] = d->Vc[u][k-1][j+1][i-1];
      V[u-1][3] = d->Vc[u][k-1][j+1][i];
      V[u-1][4] = d->Vc[u][k][j][i-1];
      V[u-1][5] = d->Vc[u][k][j][i];
      V[u-1][6] = d->Vc[u][k][j+1][i-1];
      V[u-1][7] = d->Vc[u][k][j+1][i];
    }
  }
  else if (xI > x1[i] && yI < x2[j] && zI < x3[k]){
    tag_if = 6;
    xa = x1[i]; xb = x1[i+1];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j-1][i];
      V[u-1][1] = d->Vc[u][k-1][j-1][i+1];
      V[u-1][2] = d->Vc[u][k-1][j][i];
      V[u-1][3] = d->Vc[u][k-1][j][i+1];
      V[u-1][4] = d->Vc[u][k][j-1][i];
      V[u-1][5] = d->Vc[u][k][j-1][i+1];
      V[u-1][6] = d->Vc[u][k][j][i];
      V[u-1][7] = d->Vc[u][k][j][i+1];
    }
  }
  else if (xI < x1[i] && yI < x2[j] && zI < x3[k]){
    tag_if = 7;
    xa = x1[i-1]; xb = x1[i];
    ya = x2[j-1]; yb = x2[j];
    za = x3[k-1]; zb = x3[k];
    for (u = VX1; u < VX3+1; u++) {
      V[u-1][0] = d->Vc[u][k-1][j-1][i-1];
      V[u-1][1] = d->Vc[u][k-1][j-1][i];
      V[u-1][2] = d->Vc[u][k-1][j][i-1];
      V[u-1][3] = d->Vc[u][k-1][j][i];
      V[u-1][4] = d->Vc[u][k][j-1][i-1];
      V[u-1][5] = d->Vc[u][k][j-1][i];
      V[u-1][6] = d->Vc[u][k][j][i-1];
      V[u-1][7] = d->Vc[u][k][j][i];
    }
  }
  else {
    printf("Error in 3D interpolation: \n");
    printf("interpolation point outside of brackets. \n");
    printf("Quadrent: %i \n",tag_if);
    printf("NX1=%li, NX2=%li, NX3=%li \n", NX1, NX2, NX3);
    printf("NX1_TOT=%li, NX2_TOT=%li, NX3_TOT=%li \n", NX1_TOT, NX2_TOT, NX3_TOT);
    printf("i=%i, j=%i, k=%i \n", i, j, k);
    printf("x1m=%f x1=%f, x1p=%f, xI=%f \n", x1[i-1], x1[i], x1[i+1], xI);
    printf("x2m=%f x2=%f, x2p=%f, yI=%f \n", x2[j-1], x2[j], x2[j+1], yI);
    printf("x3m=%f x3=%f, x3p=%f, zI=%f \n", x3[k-1], x3[k], x3[k+1], zI);
    printf("r=%f, ddr=%f \n", r, ddr);
    miss_tag = 1;
  }

  // Find total volume.
  N[0] = (xb - xI)*(yb - yI)*(zb - zI);
  N[1] = (xI - xa)*(yb - yI)*(zb - zI);
  N[2] = (xb - xI)*(yI - ya)*(zb - zI);
  N[3] = (xI - xa)*(yI - ya)*(zb - zI);
  N[4] = (xb - xI)*(yb - yI)*(zI - za);
  N[5] = (xI - xa)*(yb - yI)*(zI - za);
  N[6] = (xb - xI)*(yI - ya)*(zI - za);
  N[7] = (xI - xa)*(yI - ya)*(zI - za);
  Ntot = N[0] + N[1] + N[2] + N[3] + N[4] + N[5] + N[6] + N[7];
  // Normalise volumes by total.
  N[0] /= Ntot; 
  N[1] /= Ntot; 
  N[2] /= Ntot; 
  N[3] /= Ntot;
  N[4] /= Ntot; 
  N[5] /= Ntot; 
  N[6] /= Ntot; 
  N[7] /= Ntot;
         
  // ==========================================
  // vI contains the interpolated velovities 
  // 
  // vI[0] = x velocity 
  // vI[1] = y "
  // vI[2] = z "
  //
  // V contains the velocity componants at 
  // each of the sample points.
  // 
  // V[0][:] = x velocity at each sample point.
  // V[1][:] = y velocity at each sample point.
  // V[2][:] = z velocity at each sample point.
  //
  // N contains the waighting volumes for the 
  // the corner velocities.
  // ==========================================

  vI[0] = 0.0; vI[1] = 0.0; vI[2] = 0.0;
  for (s = 0; s < 8; s++) {
    vI[0] += V[0][s]*N[s];
    vI[1] += V[1][s]*N[s];
    vI[2] += V[2][s]*N[s];
  }
  vrI = (vI[0]*x1[i] + vI[1]*x2[j] + vI[2]*x3[k])/r;

  if(isnan(vrI) || miss_tag == 1){
    printf("vrI=%f, vI[0]=%f, vI[1]=%f, vI[2]=%f Ntot=%f \n", 
            vrI, vI[0], vI[1], vI[2], Ntot);
    printf("N0=%f, N1=%f, N2=%f, N3=%f, N4=%f, N5=%f, N6=%f, N7=%f \n", 
            N[0], N[1], N[2], N[3], N[4], N[5], N[6], N[7]);
    printf("xa=%f, xI=%f, xb=%f \n", xa, xI, xb);
    printf("ya=%f, yI=%f, yb=%f \n", ya, yI, yb);
    printf("za=%f, zI=%f, zb=%f \n", za, zI, zb);
    printf("\n");
  }
  return vrI;
}
/* ********************************************************************* */


/* ********************************************************************* */
double VelocityStellarSurface2D(const Data *d, RBox *box, Grid *grid, 
                                int i, int j, int k, 
                                double *x1, double *x2, double *x3, 
                                double r)
{

  int l;
  double r_testx=0.0, r_testy=0.0, vel_mag=0.0;
  double dx=0.0, dy=0.0, dz=0.0, dr2=0.0, dr=0.0, ddr=0.0;

  dx = x1[i+1] - x1[i];
  dy = x2[j+1] - x2[j];
  dz = x3[k+1] - x3[k];
  dr2  = EXPAND(dx*dx,+dy*dy,+dz*dz);
  dr   = sqrt(dr2);
  ddr = 0.9*dr;

  if (r < 1.0 && r + ddr > 1.0) {
    vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
  }


  for (l=1; l<2; l++){

    // Upper right.
    if (x1[i] > 0.0 && x2[j] > 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)   
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
    // Lower right.
    else if (x1[i] > 0.0 && x2[j] < 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i+l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
    // lower left.
    else if (x1[i] < 0.0 && x2[j] < 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j-l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    }
    // Upper left.
    else if (x1[i] < 0.0 && x2[j] > 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0)) {
        vel_mag = TwoDimensionalInterp(d, box, grid, i, j, k, x1, x2, r, ddr);
      }
/*
      if (r < 1.0 && r_testx > 1.0 && r_testy > 1.0) {
        // these are the cells at the boundary.
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx < 1.0 && r_testy > 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i];,
                 vel_x2 = d->Vc[VX2][k][j+l][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      } else if (r < 1.0 && r_testx > 1.0 && r_testy < 1.0) {
        D_EXPAND(vel_x1 = d->Vc[VX1][k][j][i-l];,
                 vel_x2 = d->Vc[VX2][k][j][i];,
                 vel_x3 = d->Vc[VX3][k][j][i];)          
        vel_mag = sqrt(vel_x1*vel_x1 + vel_x2*vel_x2);
      }
*/
    } 
    else {
      vel_mag = 0.0;
    }
  }

  return fabs(vel_mag);
}
/* ********************************************************************* */


/* ********************************************************************* */
double VelocityStellarSurface3D(const Data *d, RBox *box, Grid *grid, 
                                int i, int j, int k, 
                                double *x1, double *x2, double *x3, 
                                double r)
{

  int l=0;
  double r_testx=0.0, r_testy=0.0, r_testz=0.0, vel_mag=0.0;
  double dx=0.0, dy=0.0, dz=0.0, dr2=0.0, dr=0.0, ddr=0.0;

  dx = x1[i+1] - x1[i];
  dy = x2[j+1] - x2[j];
  dz = x3[k+1] - x3[k];
  dr2  = EXPAND(dx*dx,+dy*dy,+dz*dz);
  dr   = sqrt(dr2);
  ddr = 0.9*dr;

  for (l=1; l<2; l++){
    // Top of sphere.
    if (x1[i] > 0.0 && x2[j] > 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] > 0.0 && x2[j] < 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] < 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] > 0.0 && x3[k] > 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k+l]*x3[k+l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    // Bottom of sphere.
    if (x1[i] > 0.0 && x2[j] > 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] > 0.0 && x2[j] < 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i+l]*x1[i+l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] < 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j-l]*x2[j-l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    }
    else if (x1[i] < 0.0 && x2[j] > 0.0 && x3[k] < 0.0) {
      r_testx = D_EXPAND(x1[i-l]*x1[i-l], + x2[j]*x2[j], + x3[k]*x3[k]);
      r_testx = sqrt(r_testx);
      r_testy = D_EXPAND(x1[i]*x1[i], + x2[j+l]*x2[j+l], + x3[k]*x3[k]);
      r_testy = sqrt(r_testy);
      r_testz = D_EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k-l]*x3[k-l]);
      r_testz = sqrt(r_testz);
      if (r < 1.0 && (r_testx > 1.0 || r_testy > 1.0 || r_testz > 1.0)) {
        vel_mag = ThreeDimensionalInterp(d, box, grid, i, j, k, x1, x2, x3, r, ddr);
      }
    } 
    else {
      vel_mag = 0.0;
    }
  }
  return vel_mag;
}
/* ********************************************************************* */

