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
void InitMagneticField(double *magnetic_field, double x1, double x2, double x3, Star star1);
void CAKAcceleration(const Data *d, Grid *grid, Star star1, int i, int j, int k);
double FiniteDiskCorrection(double *gradV, double vx1, double x1, double alpha);
void VelocityGradientVector(const Data *d, double *x1, double *x2, double *x3,
                            int i, int j, int k, double *gradV);
void AccelVectorRadial(double *gline, const Data *d, double f, double x1, 
                       double *gradV, Star star1, int i, int j, int k);
void AccelVectorNonRadial(double *gline, const Data *d, double f, double x1, 
                          double *gradV, Star star1, int i, int j, int k);

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

  double velocity, magnetic_field[3];

  Star star1;
  InitStar1(&star1);

  // Set the minimum pressure via the stellar temperature.
  g_smallPressure = v[RHO]*star1.temperature/(KELVIN*star1.mean_mol);

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

  velocity = star1.terminal_velocity
             *pow(1.0 - star1.radius/x1, star1.vel_law_exponent);
  v[RHO] = (star1.mass_loss/(4.0*CONST_PI*velocity*x1*x1));
#if EOS == IDEAL
  v[PRS] = (v[RHO]*star1.temperature/(KELVIN*star1.mean_mol));
#endif                                                               
  EXPAND(v[VX1] = velocity;,
         v[VX2] = 0.0;,
         v[VX3] = 0.0;)
  v[TRC] = 0.0;

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
  double magnetic_field[3];
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
 * TODO None
 *********************************************************************** */
{        

  int i, j, k, ghost;

  double vradial, vtheta, vphi, magnetic_field[3];

  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;
  double *dx1 = grid[IDIR].dx;
  double *dx2 = grid[JDIR].dx;
  double *dx3 = grid[KDIR].dx;
  double ***vx1 = d->Vc[VX1];
  double ***vx2 = d->Vc[VX2];                                                  
  double ***vx3 = d->Vc[VX3];
  double ***gLx1 = d->gL[0];
  double ***gLx2 = d->gL[1];
  double ***gLx3 = d->gL[2];
  double ***rho = d->Vc[RHO];
#if EOS == IDEAL
  double ***prs = d->Vc[PRS];
#endif
#if PHYSICS == MHD
  double ***bx1 = d->Vc[BX1];
  double ***bx2 = d->Vc[BX2];
  double ***bx3 = d->Vc[BX3];
#endif

  Star star1;
  InitStar1(&star1);

  ghost = (NX1_TOT - NX1)/2;

  if(g_stepNumber < 2){
    double Bcgs = star1.Bfield*UNIT_B;
    double M_dot_cgs = star1.mass_loss*UNIT_MASS/UNIT_TIME;
    double v_inf_cgs = star1.terminal_velocity*UNIT_VELOCITY;
    printf("Bcgs=%e, M_dotcgs=%e, Edd_gam=%e , Omega=%e  \n", 
           Bcgs, M_dot_cgs/6.35e25, star1.Eddington, 
           star1.rotational_velocity*UNIT_VELOCITY);
  }


  if(side == X1_BEG){                          
    BOX_LOOP(box,k,j,i){ 
  
#if PHYSICS == MHD
      if (star1.eta < 5000.0) {
#endif
        EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
               vtheta = 0.0;,
               vphi = 0.0;)
#if PHYSICS == MHD
      } else if (star1.eta > 5000.0) {
        EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
               vtheta = 2.0*vx2[k][j][ghost] - vx2[k][j][ghost+1];,
               vphi = 0.0;)
      }
#endif

      if (vradial > star1.sound_speed){
        EXPAND(vradial = star1.sound_speed;,
               vtheta = vtheta;,
               vphi = vphi;)
      } else if (vradial < -star1.sound_speed){
        EXPAND(vradial = -star1.sound_speed;,
               vtheta = vtheta;,
               vphi = vphi;)
      }

      if (vtheta > star1.sound_speed){
        EXPAND(vradial = vradial;,
               vtheta = star1.sound_speed;,
               vphi = vphi;)
      } else if (vradial < -star1.sound_speed){
        EXPAND(vradial = vradial;,
               vtheta = -star1.sound_speed;,
               vphi = vphi;)
      }

      EXPAND(vx1[k][j][i] = vradial;,
             vx2[k][j][i] = vtheta;,
             vx3[k][j][i] = vphi;)
      rho[k][j][i] = star1.mass_loss
        /(4.0*CONST_PI*(star1.sound_speed/star1.surface_rho_param));          
#if EOS == IDEAL
      prs[k][j][i] = rho[k][j][i]*star1.temperature/(KELVIN*star1.mean_mol);          
#endif

#if PHYSICS == MHD   
#if BACKGROUND_FIELD == NO
      InitMagneticField(magnetic_field, x1[i], x2[j], x3[k], star1);
      EXPAND(bx1[k][j][i] = magnetic_field[0];,            
             bx2[k][j][i] = magnetic_field[1];,                 
             bx3[k][j][i] = magnetic_field[2];) 
#endif
#if BACKGROUND_FIELD == YES
      EXPAND(bx1[k][j][i] = 0.0;,
             bx1[k][j][i] = 0.0;,
             bx1[k][j][i] = 0.0;)
#endif
#endif                                                            
    }
  }

                                                                   
#if CAK == YES
  if(side == 0){
    DOM_LOOP(k,j,i){
      CAKAcceleration(d, grid, star1, i, j, k);
#if EOS == IDEAL
      if (d->Vc[PRS][k][j][i] < (rho[k][j][i])*star1.temperature
                                 /(KELVIN*star1.mean_mol)){
        d->Vc[PRS][k][j][i] = (rho[k][j][i])*star1.temperature
                               /(KELVIN*star1.mean_mol);
      }
#endif
    }
  }
#endif

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
  if (x1 > 1.0){
    g[IDIR] = star1.gravity/x1/x1 + gla[0];
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
  } else {
    g[IDIR] = star1.gravity/x1/x1;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
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
  g[IDIR] = star1.gravity/x1/x1;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
  return;
}
/* ********************************************************************* */
#endif
#endif

/* ********************************************************************* */
void CAKAcceleration(const Data *d, Grid *grid, Star star1, int i, int j, int k)
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

  double temp, f;
  double gradV[3], gline[3];

  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;

  VelocityGradientVector(d, x1, x2, x3, i, j, k, gradV);

  f = FiniteDiskCorrection(gradV, d->Vc[VX1][k][j][i], x1[i], star1.alpha);

  AccelVectorRadial(gline, d, f, x1[i], gradV, star1, i, j, k);

  EXPAND(d->gL[0][k][j][i] = gline[0];,
         d->gL[1][k][j][i] = gline[1];,
         d->gL[2][k][j][i] = gline[2];)


  return;
}
/* ********************************************************************* */

/* ********************************************************************* */
void VelocityGradientVector(const Data *d, double *x1, double *x2, double *x3,
                            int i, int j, int k, double *gradV)
/*!
 * Calculate components of the velocity gradient.
 *
 * \param [in]       d     pointer to the main PLUTO data structure
 * \param [in]       x1    x1 direction array
 * \param [in]       x2    x2 direction array
 * \param [in]       x3    x3 direction array
 * \param [in]       i     x1 direction index
 * \param [in]       j     x2 direction index
 * \param [in]       k     x3 direction index
 * \param [in, out]  dvdx  array containing velocity gradient components
 *
 * \return void
 *
 * \TODO None
 *********************************************************************** */
{

  int Nghost, mark=0;
  double dxi, dxim1, dvdx1, dvdx2, dvdx3;

  dxi = x1[i+1] - x1[i];
  dxim1 = x1[i] - x1[i-1];

  dvdx1 = fabs(-dxi*d->Vc[VX1][k][j][i-1]/(dxim1*(dxi + dxim1)) 
            + (dxi - dxim1)*d->Vc[VX1][k][j][i]/(dxi*dxim1) 
            + dxim1*d->Vc[VX1][k][j][i+1]/(dxi*(dxi + dxim1)));
  dvdx2 = 0.0;
  dvdx3 = 0.0;

  dvdx1 = fmax(dvdx1, 1.0e-8);
  dvdx2 = fmax(dvdx2, 1.0e-8);
  dvdx3 = fmax(dvdx3, 1.0e-8);

/*
  if (fabs(dvdx1) < 1.0e-8 || isnan(dvdx1)) {
    dvdx1 = 1.0e-8;
  } else if (fabs(dvdx2) < 1.0e-8 || isnan(dvdx2)) {
    dvdx2 = 1.0e-8;
  } else if (fabs(dvdx3) < 1.0e-8 || isnan(dvdx3)) {
    dvdx3 = 1.0e-8;
  }
*/
  gradV[0] = dvdx1;
  gradV[1] = dvdx2;
  gradV[2] = dvdx3;

  return;
}
/* ********************************************************************* */

/* ********************************************************************* */
double FiniteDiskCorrection(double *gradV, double vx1, double x1, double alpha)
/*!
 * Calculate finite disk correction factor.
 *
 * \param [in]  gradV  velocity gradient vector
 * \param [in]  vx1    velocity at the current cell
 * \param [in]  x1     x1 position 
 * \param [in]  alpha  CAK force multiplyer
 *
 * \return f
 *
 * \TODO None
 *********************************************************************** */
{

  double beta_op, opa, oma, sigma, f, nu2_c;

  opa = 1.0 + alpha;
  oma = 1.0 - alpha;

  if (fabs(gradV[0]) > 1.0e-8){

    beta_op = (1.0 - vx1/(gradV[0]*x1)) * pow(1.0/x1, 2);

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
void AccelVectorRadial(double *gline, const Data *d, double f, double x1, 
                       double *gradV, Star star1, int i, int j, int k)
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

  ke = 4.0*CONST_PI*UNIT_G*star1.mass*UNIT_c
         *star1.Eddington/star1.luminosity;

  B = d->Vc[RHO][k][j][i]*star1.q_fac*UNIT_c*ke;

  A = 1.0/(1.0 - star1.alpha)*ke*star1.luminosity
        *star1.q_fac/(4.0*CONST_PI*UNIT_c);

  gline[0] = f*A*pow(x1, -2)*pow(gradV[0]/B, star1.alpha);

  //if (isnan(gradV[0])) {
    gline[0] = 0.0;
  //}

#if EOS == IDEAL
  // Accounting for total ionisation at high temp 
  // and for recombination at low temp.
  temp = d->Vc[PRS][k][j][i]*KELVIN*star1.mean_mol/d->Vc[RHO][k][j][i];
  gline[0] *= exp(-4.0*log(2.0)*pow((2.0 - temp/star1.temperature 
             - star1.temperature/temp), 2));
#endif

  gline[1] = 0.0;
  gline[2] = 0.0;

  if (isnan(d->Vc[RHO][k][j][i]) || isnan(d->Vc[PRS][k][j][i]) || isnan(d->Vc[VX1][k][j][i]) || isnan(d->Vc[VX2][k][j][i]) || isnan(d->Vc[VX3][k][j][i])) {
    printf("d->Vc[RHO][k][j][i]=%e \n", d->Vc[RHO][k][j][i]);
    printf("d->Vc[PRS][k][j][i]=%e \n", d->Vc[PRS][k][j][i]);
    printf("d->Vc[VX1][k][j][i]=%e \n", d->Vc[VX1][k][j][i]);
    printf("d->Vc[VX2][k][j][i]=%e \n", d->Vc[VX2][k][j][i]);
    printf("d->Vc[VX3][k][j][i]=%e \n", d->Vc[VX3][k][j][i]);
    printf("gradV[0]=%e \n", gradV[0]);
    printf("gline[0]=%e \n", gline[0]);
  }

  return;
}
/* ********************************************************************* */

/* ********************************************************************* */
void AccelVectorNonRadial(double *gline, const Data *d, double f, double x1, 
                          double *gradV, Star star1, int i, int j, int k)
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

  g_mag = f*A*pow(x1, -2)*pow(gradV_mag/B, star1.alpha);

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

  double x, y, z, xp, yp, zp, r;
  double br, btheta, bphi, bx, by, bz, bxp, byp, bzp, rp, rp2;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;

  star1.Bfield_angle *= 0.0174532925;

  // Convert to Cartesian.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);

  // Rotate Cartesian coordiantes.
  xp = x*cos(star1.Bfield_angle) - z*sin(star1.Bfield_angle);
  yp = y;
  zp = x*sin(star1.Bfield_angle) + z*cos(star1.Bfield_angle);
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  // Calculate b-field components in rotated frame.
  bx = 3.0*xp*zp*star1.Bfield*pow(rp,-5);
  by = 3.0*yp*zp*star1.Bfield*pow(rp,-5);
  bz = (3.0*pow(zp,2)-pow(rp,2))*star1.Bfield*pow(rp,-5);

  // Rotate B-field vector componets.  
  bxp = bx*cos(star1.Bfield_angle) + bz*sin(star1.Bfield_angle);
  byp = by;
  bzp = -bx*sin(star1.Bfield_angle) + bz*cos(star1.Bfield_angle);

  // Define spherical basis vectors.
  a11 = sin(x2)*cos(x3); a12 = sin(x2)*sin(x3); a13 = cos(x2);
  a21 = cos(x2)*cos(x3); a22 = cos(x2)*sin(x3); a23 = -sin(x2);
  a31 = -sin(x3);        a32 = cos(x3);         a33 = 0.0;

  // Change basis back to spherical polar.
  br = bxp*a11 + byp*a12 + bzp*a13;
  btheta = bxp*a21 + byp*a22 + bzp*a23;
  bphi = bxp*a31 + byp*a32 + bzp*a33;

  magnetic_field[0] = br;
  magnetic_field[1] = btheta; 
  magnetic_field[2] = bphi;

  return;
}
/* ********************************************************************* */


/* ********************************************************************* */
void VectorFieldRotation()
/*!
 *
 *
 *
 *
 *
 *********************************************************************** */

{

  return;
}


