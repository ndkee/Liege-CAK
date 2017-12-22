/* ********************************************************************* */
/*!
 * Initilisation file for a radiativly driven stellar wind with a 
 * non-rigid dipole configuration magnetic field.
 *
 * The boundary and initial conditions are taken 
 * from Runacres and Owocki (2002)    
 *
 * The method for calculating the radiative acceleration 
 * comes from CAK (1975)      
 *
 * The model only works with polar coordinates in 2D, with the 
 * MHD module. 1D, 3D, HD, RHD, RMHD and other geometries do not 
 * work at the moment.
 *
 * Authors: Simon Daley-Yates & Asif ud-Doula and Dylan Kee. 
 *
 * Last update: 30/11/2017
 *
 * TODO: Cooling function seems to be limiting the time step to much. 
 *       This may require playing around with the g_maxCoolingRate 
 *       parameter. At the moment I feel comfortable using this code for 
 *       radial line acceleration with rotation and isothermal eos. ideal 
 *       works but suffers from timestep decay.
 *
 * Asif: added colliding wind outer boundary condition-> set x1end b.c. to userdef;
 * added vector potential
 * to run with CT options with no background field; need: rotation 
 * for ROTATION_FRAME==NO. Dec 4, 2017
 ********************************************************************** */
#include "pluto.h"                                                                  

/* ********************************************************************* */
typedef struct STAR
/*!
 * Type definition for the star. This structure contains all the 
 * variables that are specific to the star. A structure of this 
 * type is initialised by every function that needs stellar parameters.
 * In this way, Only one function needs to calculate the parameters.
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
  double rmax;
  double l_separation;
  double r_secondary;
  double vinf_secondary;
  double mdot_secondary;

} Star;
/* ********************************************************************* */

void InitStar1(Star *star1);
void InitMagneticField(double *magnetic_field, double *vector_potential, double x1, double x2, double x3, Star star1);
#if CAK == YES
void CAKAcceleration(const Data *d, Grid *grid, Star star1, int i, int j, int k);
double FiniteDiskCorrection(double *gradV, double vx1, double x1, double alpha);
void VelocityGradientVector(const Data *d, double *x1, double *x2, double *x3,
                            int i, int j, int k, double *gradV);
void AccelVectorRadial(double *gline, const Data *d, double f, double x1, 
                       double *gradV, Star star1, int i, int j, int k);
void AccelVectorNonRadial(double *gline, const Data *d, double f, double x1, 
                          double *gradV, Star star1, int i, int j, int k);
void   StartUpGCAK(Data *, Grid *, double);
void   ofdwts(Data *, Grid *, Star, double);
void   gauleg(double, double, double *, double *, int);
double Weights(double, double, double, double, double, double, double, double, double, double, double, double);
double RX(double, double);
void   LimbSearch(double, double, double, double, double, double *);
double TPfunc(double, double);
void   CoordTransfm(double, double, double, double, double, double, double *);
double brent (double, double, double, double, double, double);
double bisect0 (double, double, double, double, double, double);
double sign(double, double);
void   RminSolve (double, double, double, double, double, double, double *);
double FuncPrime (double, double, double, double, double, double);
void   findEqCross (Data *, Grid *);
void   gcak3d(const Data *, Grid *, Star, int, int, int, double *);
void   sobolev(const Data *, Grid *, Star, int, int, int, double *);
double rayOptDepth(const Data *, const Grid *, int, int, int, int, int);
#endif

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

  star1->mean_mol = g_inputParam[Mean_mol_weight];

  star1->sound_speed = sqrt(UNIT_kB*star1->temperature
    /(star1->mean_mol*(CONST_AH/UNIT_MASS)*CONST_amu));

  star1->escape_velocity = sqrt(2.0*UNIT_G*star1->mass*(1.0 - star1->Eddington));                              

  star1->terminal_velocity = star1->escape_velocity
    *sqrt(star1->alpha/(1.0 - star1->alpha));

  star1->vel_law_exponent = g_inputParam[Velocity_exponent];

// Asif: eta_* is defined in terms of B_equator, B_pole is 2xB_eq for dipole
  star1->Bfield = 2.0*sqrt(star1->eta
                 *star1->mass_loss*UNIT_MASS/UNIT_TIME
                 *star1->terminal_velocity*UNIT_VELOCITY
                 /pow(UNIT_LENGTH, 2))/UNIT_B;

  star1->gravity = -UNIT_G*star1->mass*(1.0 - star1->Eddington);

  star1->rotational_velocity = g_inputParam[Rotation]*sqrt(UNIT_G*star1->mass);

  star1->Bfield_angle = g_inputParam[Magnetic_incl];

  star1->surface_rho_param = g_inputParam[Cs_p];

  star1->rmax=g_inputParam[R_max];
  
  star1->l_separation=g_inputParam[D_separation];

  star1->r_secondary=g_inputParam[R_secondary];
  
  star1->vinf_secondary=g_inputParam[Vinf_secondary]/UNIT_VELOCITY;

  star1->mdot_secondary=g_inputParam[Mdot_secondary_ratio]*star1->mass_loss;
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

  #if PHYSICS == MHD
  #if BACKGROUND_FIELD == NO
  double magnetic_field[3]; 
  double vector_potential[3];
  #endif
  #endif

  Star star1;
  InitStar1(&star1);

  rho0 = star1.mass_loss/(4.0*CONST_PI
           *(star1.sound_speed/star1.surface_rho_param));

  // Set the minimum pressure via the stellar temperature.
  g_smallPressure = rho0*star1.temperature/(KELVIN*star1.mean_mol);
  g_smallPressure =g_smallPressure*1.e-12; 

  #if COOLING !=NO
  // Dylan,Asif: this sets minimum T of the sim.
  // there is also gmaxCoolingRate that can limit timestep.
  // it will be worthwhile to use it as well.
  g_minCoolingTemp = star1.temperature;
  g_maxCoolingRate = 0.2;
  #endif 

  #if EOS == IDEAL
  g_gamma = 1.666667;
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
  InitMagneticField(magnetic_field, vector_potential, x1, x2, x3, star1);
  v[BX1] = magnetic_field[0];
  v[BX2] = magnetic_field[1];
  v[BX3] = magnetic_field[2];
  #if ASSIGN_VECTOR_POTENTIAL == YES
  v[AX1] = vector_potential[0];
  v[AX2] = vector_potential[1];
  v[AX3] = vector_potential[2];
  #endif
  #endif
  #if BACKGROUND_FIELD == YES
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;
  #if ASSIGN_VECTOR_POTENTIAL == YES
  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
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
  double magnetic_field[3], vector_potential[3];
  Star star1;
  InitStar1(&star1);
  InitMagneticField(magnetic_field, vector_potential, x1, x2, x3, star1);
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

  int i, j, k, ghost, nv;

  double vradial, vtheta, vphi, magnetic_field[3], vector_potential[3], temp;
  double tmax, t2max, phi_secondary;
  double r_to_secondary, xd, yd, zd;
  double dB, *Area, *dV, eta_max;
  double temp1, temp2, temp3, temp4, temp5;

  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;
  double *dx1 = grid[IDIR].dx;
  double *dx2 = grid[JDIR].dx;
  double *dx3 = grid[KDIR].dx;
  double ***vx1 = d->Vc[VX1];
  double ***vx2 = d->Vc[VX2];                                                  
  double ***vx3 = d->Vc[VX3];
  #if CAK == YES
  double ***gLx1 = d->gL[0];
  double ***gLx2 = d->gL[1];
  double ***gLx3 = d->gL[2];
  #endif
  double ***rho = d->Vc[RHO];
  #if EOS == IDEAL
  double ***prs = d->Vc[PRS];
  #endif
  #if PHYSICS == MHD
  double ***bx1 = d->Vc[BX1];
  double ***bx2 = d->Vc[BX2];
  double ***bx3 = d->Vc[BX3];
  #ifdef STAGGERED_MHD
  double ***bx1s = d->Vs[BX1s];
  double ***bx2s = d->Vs[BX2s];
  double ***bx3s = d->Vs[BX3s];
  #endif
  #endif

  Star star1;
  InitStar1(&star1);

  Area  = grid->A;
  dV = grid->dV;
//	set maximum eta_* for linear extrapolation
  eta_max=1000.0 ;
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

    if(box->vpos == CENTER){

      BOX_LOOP(box,k,j,i){ 
          dB = (Area[ghost]*vx1[k][j][ghost] - Area[ghost+1]*vx1[k][j][ghost+1])/dV[ghost];

        #if PHYSICS == MHD
        if(star1.eta < eta_max){
        #endif
          //EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
          EXPAND(vradial =  (vx1[k][j][ghost]*Area[ghost] + dV[i]*dB)/Area[i];, 
                 vtheta = 0.0;,
                 vphi = 0.0;)
        #if PHYSICS == MHD
        }else if(star1.eta > eta_max){
          //EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];, 
          EXPAND(vradial =  (vx1[k][j][ghost]*Area[ghost] + dV[i]*dB)/Area[i];, 
                 vtheta =  2.0*vx2[k][j][ghost] - vx2[k][j][ghost+1];, 
                 vphi = 0.0;)
        }
        #endif

        if(vradial > star1.sound_speed){
          EXPAND(vradial = star1.sound_speed;,
                 vtheta = vtheta;,
                 vphi = vphi;)
        }else if(vradial < -star1.sound_speed){
          EXPAND(vradial = -star1.sound_speed;,
                 vtheta = vtheta;,
                 vphi = vphi;)
        }

        if(vtheta > star1.sound_speed){
          EXPAND(vradial = vradial;,
                 vtheta = star1.sound_speed;,
                 vphi = vphi;)
        }else if(vradial < -star1.sound_speed){
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
        #ifndef STAGGERED_MHD
        #if BACKGROUND_FIELD == NO
        if(star1.eta < eta_max){
        InitMagneticField(magnetic_field, vector_potential, x1[i], x2[j], x3[k], star1);
        bx1[k][j][i] = magnetic_field[0];
        bx2[k][j][i] = 2.0*bx2[k][j][ghost] - bx2[k][j][ghost+1];
        bx3[k][j][i] = 2.0*bx3[k][j][ghost] - bx3[k][j][ghost+1];
        }else if(star1.eta > eta_max){
        InitMagneticField(magnetic_field, vector_potential, x1[i], x2[j], x3[k], star1);
        bx1[k][j][i] = magnetic_field[0];
        bx2[k][j][i] = magnetic_field[1];
        bx3[k][j][i] = magnetic_field[2];
        }
        #endif
        #if BACKGROUND_FIELD == YES
        bx1[k][j][i] = 0.0;
        bx2[k][j][i] = 0.0;
        bx3[k][j][i] = 0.0;
        #endif
        #endif
        #endif
      }

    #ifdef STAGGERED_MHD
    }else if(box->vpos == X2FACE){  
      BOX_LOOP(box,k,j,i){
        #if BACKGROUND_FIELD == NO
        if(star1.eta < eta_max){
        bx2s[k][j][i] = 2.0*bx2s[k][j][ghost] - bx2s[k][j][ghost+1];
        }else if(star1.eta > eta_max){
        InitMagneticField(magnetic_field, vector_potential, x1[i], x2[j] + dx2[j]/2.0, x3[k], star1);
        bx2s[k][j][i] = magnetic_field[1];
        }
        #endif
        #if BACKGROUND_FIELD == YES
        bx2s[k][j][i] = 0.0;
        #endif
      }
    }else if(box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){
        #if BACKGROUND_FIELD == NO
        if(star1.eta < eta_max){
        bx3s[k][j][i] = 2.0*bx3s[k][j][ghost] - bx3s[k][j][ghost+1];
        }else if(star1.eta > eta_max){
        InitMagneticField(magnetic_field, vector_potential, x1[i], x2[j], x3[k] + dx3[k]/2.0, star1);
        bx3s[k][j][i] = magnetic_field[2];
        }
        #endif
        #if BACKGROUND_FIELD == YES
        bx3s[k][j][i] = 0.0;
        #endif
      }
    #endif // staggered
    } // outside #ifdef STAGGERED_MHD to close first if

  }

  /* - Define user boundary for outer Rmax. - */
  if(side == X1_END){  // X1_END boundary
    if(box->vpos == CENTER){
     for(nv = 0; nv < NVAR; nv++){
        BOX_LOOP(box,k,j,i){
          #if GEOMETRY == CARTESIAN
          d->Vc[nv][k][j][i] = 2.0*d->Vc[nv][k][j][i-1] - d->Vc[nv][k][j][i-2];
          #else
          dB = (Area[i-1]*d->Vc[nv][k][j][i-1] - Area[i-2]*d->Vc[nv][k][j][i-2])/dV[i-1];
          d->Vc[nv][k][j][i] = (d->Vc[nv][k][j][i-1]*Area[i-1] + dV[i]*dB)/Area[i];
          #endif
          // Make changes here
          //
          // Calculate maximum cos(cone_angle) sustained by the secondary star.
          t2max =(star1.rmax/star1.l_separation);
          // Allow for the secondary to orbit in XY plane, make phi(time)
          phi_secondary=0.0;
          tmax=sin(x2[j])*cos(x3[k]-phi_secondary);
          // Place holder for general coordinates for the secondary in Cartesian
          xd=star1.l_separation;
          yd=0.0;
          zd=0.0;
          //printf("tmax=%e, t2max=%e, Separation=%e , star1.vinf_secondary=%e  \n",
          //tmax, t2max, star1.l_separation,
          //star1.vinf_secondary*UNIT_VELOCITY);
          if(tmax > t2max) {
            // Compute the radius from the seondary to boundary. 
            r_to_secondary=sqrt(pow(star1.rmax,2)+pow(xd,2)+pow(yd,2)+pow(zd,2)
                          -2.*star1.rmax*(zd*cos(x2[j])+sin(x2[j])*(xd*cos(x3[k])+yd*sin(x3[k]))));
            // nv == VX1; PRS
            vradial=star1.vinf_secondary
                   *pow(1.0 - star1.r_secondary/r_to_secondary, star1.vel_law_exponent);

            if (nv == RHO){
              rho[k][j][i] = star1.mdot_secondary
                           /(4.0*CONST_PI*vradial*pow(r_to_secondary/star1.r_secondary,2));
            }
            if (nv == VX1){
              vx1[k][j][i] = vradial*(star1.rmax-zd*cos(x2[j])-sin(x2[j])*(xd*cos(x3[k])+yd*sin(x3[k])))
                           /r_to_secondary;  
            }
            if (nv == VX2){
              vx2[k][j][i] = vradial*(zd*sin(x2[j]) - cos(x2[j])* (xd*cos(x3[k]) + yd*sin(x3[k])))/r_to_secondary;
            }
            if (nv == VX3){
              vx3[k][j][i] = -vradial*(yd*cos(x3[k])-xd*sin(x3[k]))/r_to_secondary;       
            }
            #if EOS == IDEAL
            if (nv == PRS){
              temp = KELVIN*d->Vc[PRS][k][j][i]*star1.mean_mol/d->Vc[RHO][k][j][i];
              if (temp/star1.temperature < 10.){
              temp=star1.temperature; 
              }
              if (temp/star1.temperature > 10.0) {
              temp=star1.temperature*10.;
              }
              prs[k][j][i] = rho[k][j][i]*temp/(KELVIN*star1.mean_mol);
            }
            #endif
          }
        #if PHYSICS == MHD
        #ifndef STAGGERED_MHD
        #if BACKGROUND_FIELD == NO
            if (nv == BX1) {
            //bx1[k][j][i] = 0.0;
            }
            if (nv == BX2) {
            //bx2[k][j][i] = 0.0;
            }
            if (nv == BX3) {
            //bx3[k][j][i] = 0.0;
            }
        #endif
        #if BACKGROUND_FIELD == YES
            if (nv == BX1) {
            bx1[k][j][i] = 0.0;
            }
            if (nv == BX2) {
            bx2[k][j][i] = 0.0;
            }
            if (nv == BX3) {
            bx3[k][j][i] = 0.0;
            }
        #endif
        #endif
        #endif
        }
      }
    #ifdef STAGGERED_MHD
    }else if(box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){
        #if BACKGROUND_FIELD == NO
        dB = (Area[i-1]*bx2s[k][j][i-1] - Area[i-2]*bx2s[k][j][i-2])/dV[i-1];
        bx2s[k][j][i] = (bx2s[k][j][i-1]*Area[i-1] + dV[i]*dB)/Area[i];
        //bx2s[k][j][i] = 2.0*bx2s[k][j][ghost] - bx2s[k][j][ghost+1];
        #endif
        #if BACKGROUND_FIELD == YES
        bx2s[k][j][i] = 0.0;
        #endif
      }
    }else if(box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){
        #if BACKGROUND_FIELD == NO
        dB = (Area[i-1]*bx3s[k][j][i-1] - Area[i-2]*bx3s[k][j][i-2])/dV[i-1];
        bx3s[k][j][i] = (bx3s[k][j][i-1]*Area[i-1] + dV[i]*dB)/Area[i];
        //bx3s[k][j][i] = 2.0*bx3s[k][j][ghost] - bx3s[k][j][ghost+1];
        #endif
        #if BACKGROUND_FIELD == YES
        bx3s[k][j][i] = 0.0;
        #endif
      }
    #endif // staggered
    } // outside #ifdef STAGGERED_MHD to close first if

  } // end X1END
                                                                   
  #if CAK == YES
  if(side == 0){
    DOM_LOOP(k,j,i){
      #if EOS == IDEAL
      temp = KELVIN*d->Vc[PRS][k][j][i]*star1.mean_mol/d->Vc[RHO][k][j][i];
      if (temp < star1.temperature) {
        d->Vc[PRS][k][j][i] = star1.temperature*d->Vc[RHO][k][j][i]
                                /(KELVIN*star1.mean_mol);
      }
      // Avoid large T, especially along the axes, set Tmax=9.e9, 
      // close to max in tabualted cooling, Asif
      if (temp > 9.e9) {
        d->Vc[PRS][k][j][i] = 9.e9*d->Vc[RHO][k][j][i]
                               /(KELVIN*star1.mean_mol);
      }
      #endif
      double gline[3];
      if(g_inputParam[CAK_ifrc] == 0){
        sobolev(d, grid, star1, k, j, i, gline);
      }else{
        gcak3d(d, grid, star1, k, j, i, gline);
      }
      EXPAND(d->gL[0][k][j][i] = gline[0];,
             d->gL[1][k][j][i] = gline[1];,
             d->gL[2][k][j][i] = gline[2];)
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
    g[JDIR] = gla[1];
    g[KDIR] = gla[2];
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
void InitMagneticField(double *magnetic_field, double *vector_potential, 
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
  double ar, atheta, aphi, ax, ay, az, axp, ayp, azp;

  star1.Bfield_angle *= 0.0174532925;

  // Convert to Cartesian.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);

  // Rotate Cartesian coordinates.
  xp = x*cos(star1.Bfield_angle) - z*sin(star1.Bfield_angle);
  yp = y;
  zp = x*sin(star1.Bfield_angle) + z*cos(star1.Bfield_angle);
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  // Calculate b-field components in rotated frame.
  bx = 3.0*xp*zp*star1.Bfield*pow(rp,-5);
  by = 3.0*yp*zp*star1.Bfield*pow(rp,-5);
  bz = (3.0*pow(zp,2)-pow(rp,2))*star1.Bfield*pow(rp,-5);

  ax = -star1.Bfield/2.*yp*pow(rp,-3);
  ay =  star1.Bfield/2.*xp*pow(rp,-3);
  az =  0.0;

  // Rotate B-field vector componets.  
  bxp = bx*cos(star1.Bfield_angle) + bz*sin(star1.Bfield_angle);
  byp = by;
  bzp = -bx*sin(star1.Bfield_angle) + bz*cos(star1.Bfield_angle);

  axp = ax*cos(star1.Bfield_angle) + az*sin(star1.Bfield_angle);
  ayp = ay;
  azp = -ax*sin(star1.Bfield_angle) + az*cos(star1.Bfield_angle);

  // Define spherical basis vectors.
  a11 = sin(x2)*cos(x3); a12 = sin(x2)*sin(x3); a13 = cos(x2);
  a21 = cos(x2)*cos(x3); a22 = cos(x2)*sin(x3); a23 = -sin(x2);
  a31 = -sin(x3);        a32 = cos(x3);         a33 = 0.0;

  // Change basis back to spherical polar.
  br     = bxp*a11 + byp*a12 + bzp*a13;
  btheta = bxp*a21 + byp*a22 + bzp*a23;
  bphi   = bxp*a31 + byp*a32 + bzp*a33;

  ar     = axp*a11 + ayp*a12 + azp*a13;
  atheta = axp*a21 + ayp*a22 + azp*a23;
  aphi   = axp*a31 + ayp*a32 + azp*a33;

  magnetic_field[0] = br;
  magnetic_field[1] = btheta; 
  magnetic_field[2] = bphi;

  vector_potential[0] = ar;
  vector_potential[1] = atheta;
  vector_potential[2] = aphi;

  return;
}

#if CAK == YES
/* ********************************************************************* */
void gcak3d(const Data *data, Grid *grid, Star star, int kl, int jl, int il, double *gline){
    //
    //      multiray CAK method for 3-D radiation force:
    //  	Uses CAK escape probability along rays,
    //      allowing for full 3-D wind structure, but
    //       still assuming 2-D asymmetric stellar radiation field,
    //       though possibly with an Oblate Finite Disk (OFD)
    //
    //     Operation controlled by ifrc switch as follows:
    //
    //         .eq. 0 =>  gr=gp=gt=0.          //no radiative force, defaults to Sobolev
    //         .eq.-1 =>  gr=gdx,gp=gt=0.	  //pure radial      (w/ OFD)
    //         .eq.-2 =>  gr=gdx,gp=gdz,gt=0	  //phi, but no thet (w/ OFD)
    //         .eq.-3 =>  gr=gdx,gp=gt=0.	  //thet, but no phi (w/ OFD)
    //         .eq.-4 =>  gr=gdx,gp=gdz,gt=gdy //both thet+phi, in circ. disk approx.
    //         .eq.-5 =>  gr=gdx,gp=gdz,gt=gdy //both thet+phi, w/  OFD, but fw=1.
    //         .eq.-6 =>  gr=gdx,gp=gdz,gt=gdy //both thet+phi, w/  OFD, & fw~I_c, w/ gd&ld
    //         .eq.-7 =>  gr="    "       "    //w/ OFD, ld, but NO gd.
    //
    //     ny=MIN(abs(nyside),1) rays, with:
    // 		y,wy = Gauss-Leg quad (nyside.gt.0)
    //		y,wy = Simpsons rule  (nyside.le.0)
    //
    //     Revision history:
    //     6/07/99: add vtherm<0 option to include lbc profile effect
    //    ~5/30/98: add ifrc=-7 option (ld w/o gd)
    //     4/26/96: first full OFD test version, with y~mu (vs. mu^2)
    //     4/21/96: include infrastructure for OFD (Oblate Finite Disk) factors
    //     3/03/96: adapted from 1.5-D CAK routine, "gcak2d.f"
    //     9/28/95: adapted from gssf2d
    //     9/23/95: implement nray 1.5D options
    //

    
    int iy,ipp,ip,im;
    int ifrco;
    int i,j,k,imax,jmax,kmax;
    int kp,km,k1,k2,kdphi,jp,jm;
    int ny1,npp1,ny2,npp2;
    int iysgn,izsgn;
    double Rmin,Rmax;
    double pprot,ppmax,ppmin;
    double requator,xmustar,dilfac,delfac,delmus;
    double dphirot,wphi1,wphi2;
    double dphi,dtheta,theta,costo,sinto,cotto;
    double cosp,sinp,cospsq,sinpsq,cost,costsq,sint,sintsq;
    double y,wy,pp,wpp,dpp,wppy,wtot;
    double r,rn,rsbr,vr,vt,vp,rho,vrbr,vtbr,vpbr;
    double dr,dy,dz,dvrdr,dvrdy,dvrdz,dvtdr,dvtdy,dvtdz,dvpdr,dvpdy,dvpdz;
    double a1,a1sq,a2,a3,a4e,a4o,tee,teo,toe,too;
    double dvpp,dsym,tmaxp,tmaxm;
    double tmp,tmp1,tmp2,tmpp,tmpm;
    double cak,q0,qbar,delta,aCAK,oma,tau,vtherm;
    double gii, gjj, gkk;
    double *xi,*dxi,*yi,*dyi,*zi,*dzi;
    double drp, drm, dthetap, dthetam, dphip, dphim, dyp, dym, dzp, dzm;
    
    const double xhabun = 0.73;
    const double xmpro  = 1.67e-24 / UNIT_MASS;
    const double xmbe   = xmpro * 2.0 / (1.0 + xhabun);
    const double elkap  = 0.2 * (1.0 + xhabun) / (UNIT_MASS*pow(UNIT_LENGTH,-2));
    
    
    Rmin = grid[IDIR].xl_glob[grid[IDIR].gbeg]; // in code units
    Rmax = grid[IDIR].xr_glob[grid[IDIR].gend]; // in code units
    
    gii = 0.;
    gjj = 0.;
    gkk = 0.;
    
    xi   = grid[IDIR].x;
    dxi  = grid[IDIR].dx;
    imax = grid[IDIR].np_tot;
    yi   = grid[JDIR].x;
    dyi  = grid[JDIR].dx;
    jmax = grid[JDIR].np_tot;
    zi   = grid[KDIR].x;
    dzi  = grid[KDIR].dx;
    kmax = grid[KDIR].np_tot;
    
    // dkee 20Jun17 make requator always star.radius, HARDCODED NO OFD

    requator = star.radius;
    
    ifrco = g_inputParam[CAK_ifrc];
    
    if (ifrco == 0) return;
    if (ifrco > 0) ifrco=-5;
    
    ny1  = g_inputParam[CAK3D_nyy];
    npp1 = g_inputParam[CAK3D_npp];
    ny2  = (int) MAX(fabs(g_inputParam[CAK3D_nyy]),1);
    npp2 = (int) MAX(2*fabs(g_inputParam[CAK3D_npp]),1);
    ppmax =  M_PI;
    ppmin = -M_PI;
    if (kmax == 1) ppmin=0.;
    if (jmax == 1){
        ppmax=M_PI/2.;
        ppmin=-ppmax;
    }
    
    aCAK   = g_inputParam[CAK_alpha];
    oma    = 1.-aCAK;
    delta  = g_inputParam[CAK_delta];
    qbar   = g_inputParam[Q_factor];
    q0     = g_inputParam[Q_factor];
    
    //
    // Keep track of phi rotation with time...
    //

    dphirot = 0.0;
    kdphi   = (int) dphirot;
    wphi2   = dphirot-kdphi;
    wphi1   = 1.-wphi2;
    // Begin Angle and Ray loops:
    //
    kp = kl+1;
    km = kl-1;
    if (kp >= kmax){
        kp=0;
    }
    if (km < 0){
        km=kmax-1;
    }
    if (kmax == 1){
        dphip = 1.;
        dphim = 1.;
    }else{
        dphip = zi[kp]-zi[kl];
        dphim = zi[kl]-zi[km];
    }
    k1 = kl-kdphi;                  //Compute rotated phi indices...
    k1 = k1-((k1-kmax)/kmax)*kmax-1; //Truncate to keep in kmax range..
    k2 = k1-1;
    if (k1 == 0){
        k2 = kmax-1;          //Special case for endpoint
    }
    jp = MIN(jl+1,jmax-1);
    jm = MAX(jl-1,0);
    if (jmax == 1){
        dthetap = 1.;
        dthetam = 1.;
        costo  = 0.;
        sinto  = 1.;
    }else{
        dthetap = yi[jp]-yi[jl];
        dthetam = yi[jl]-yi[jm];
        theta   = yi[jl];
        sinto   = sign(MAX(fabs(sin(theta)),1.e-10),sin(theta));
        costo   = cos(theta);
    }
    cotto = costo/sinto;
    
    r      = xi[il];
    vr     = data->Vc[VX1][kl][jl][il];
    vt     = data->Vc[VX2][kl][jl][il];
    vp     = data->Vc[VX3][kl][jl][il];
    rho    = data->Vc[RHO][kl][jl][il];
    vrbr   = vr/r;
    vtbr   = vt/r;
    vpbr   = vp/r;
    
    ip     = MIN(il+1,imax-1);
    im     = MAX(il-1,0);
    
    drp    =  xi[ip]-xi[il];
    drm    =  xi[il]-xi[im];

    dvrdr  = -drp/(drm * (drp + drm)) * data->Vc[VX1][kl][jl][im] + (drp - drm)/(drp * drm) * data->Vc[VX1][kl][jl][il] + drm/(drp * (drp + drm)) * data->Vc[VX1][kl][jl][ip];
    dvtdr  = -drp/(drm * (drp + drm)) * data->Vc[VX2][kl][jl][im] + (drp - drm)/(drp * drm) * data->Vc[VX2][kl][jl][il] + drm/(drp * (drp + drm)) * data->Vc[VX2][kl][jl][ip];
    dvpdr  = -drp/(drm * (drp + drm)) * data->Vc[VX3][kl][jl][im] + (drp - drm)/(drp * drm) * data->Vc[VX3][kl][jl][il] + drm/(drp * (drp + drm)) * data->Vc[VX3][kl][jl][ip];
    dyp    = r * dthetap;
    dym    = r * dthetam;
    dvrdy  = -dyp/(dym * (dyp + dym)) * data->Vc[VX1][kl][jm][il] + (dyp - dym)/(dyp * dym) * data->Vc[VX1][kl][jl][il] + dym/(dyp * (dyp + dym)) * data->Vc[VX1][kl][jp][il];
    dvtdy  = -dyp/(dym * (dyp + dym)) * data->Vc[VX2][kl][jm][il] + (dyp - dym)/(dyp * dym) * data->Vc[VX2][kl][jl][il] + dym/(dyp * (dyp + dym)) * data->Vc[VX2][kl][jp][il];
    dvpdy  = -dyp/(dym * (dyp + dym)) * data->Vc[VX3][kl][jm][il] + (dyp - dym)/(dyp * dym) * data->Vc[VX3][kl][jl][il] + dym/(dyp * (dyp + dym)) * data->Vc[VX3][kl][jp][il];
    dzp    = r * sinto * dphip;
    dzm    = r * sinto * dphim;
    dvrdz  = -dzp/(dzm * (dzp + dzm)) * data->Vc[VX1][km][jl][il] + (dzp - dzm)/(dzp * dzm) * data->Vc[VX1][kl][jl][il] + dzm/(dzp * (dzp + dzm)) * data->Vc[VX1][kp][jl][il];
    dvtdz  = -dzp/(dzm * (dzp + dzm)) * data->Vc[VX2][km][jl][il] + (dzp - dzm)/(dzp * dzm) * data->Vc[VX2][kl][jl][il] + dzm/(dzp * (dzp + dzm)) * data->Vc[VX2][kp][jl][il];
    dvpdz  = -dzp/(dzm * (dzp + dzm)) * data->Vc[VX3][km][jl][il] + (dzp - dzm)/(dzp * dzm) * data->Vc[VX3][kl][jl][il] + dzm/(dzp * (dzp + dzm)) * data->Vc[VX3][kp][jl][il];

    for(iy = 0; iy < ny2; iy++){
        y    =  data->y1d[iy];
        wy   = data->wy1d[iy];
        for(ipp = 0; ipp < npp2; ipp++){
            dpp   = (ppmax-ppmin)/npp2;
            if((iy % 3) == 0){
                pprot = 0;
            }else if((iy % 3) == 1){
                pprot = 1.0 / 3.0;
            }else{
                pprot = -1.0 / 3.0;
            }
            pp = data->pp1d[ipp];
            if(tan(pp) > 0){
                pp += pprot;
            }else{
                pp -= pprot;
            }
            wpp   = data->wpp1d[ipp];
            wppy = wpp*wy;
            sinp  = sin(pp);
            cosp  = cos(pp);
            cospsq= cosp*cosp;
            sinpsq= sinp*sinp;
            //
            //  Sum force integrands at each radius
            //
            
            delmus = 1.-data->ctpmax[ipp][jl][il]; // make y quad linear in mu, not mu**2
            cost   = 1.-y*delmus;     // (this is better for gtheta, gphi
            costsq = cost*cost;       //  integrals; tho not for gr)
            sintsq = 1.-costsq;
            sint   = sqrt(MAX(0.,sintsq));
            
            wtot=wppy*delmus*(wphi1*data->fw[ipp*npp2+iy][k1][jl][il]+wphi2*data->fw[ipp*npp2+iy][k2][jl][il]);
            
            a1     = cost;
            a1sq   = costsq;
            a2     = sint*cosp;
            a3     = sint*sinp;
            a4e    = -cost*sinto;	//even in cosp
            a4o    = -sint*costo*cosp;	//odd  in cosp
            
            tee    = a1sq*dvrdr+a2*a2*(dvtdy+vrbr)+a3*a3*(vrbr+vtbr*cotto+dvpdz); 		//even in both sinp, cosp
            teo    =  a1*a2*(dvtdr-vtbr+dvrdy);		//even in sinp, odd in cosp
            toe    =  a3*a1*(dvpdr+dvrdz-vpbr);                //odd  in sinp, even in cosp
            too    =  a3*a2*(dvpdy+dvtdz-vpbr*costo/sinto);	//odd  in both sinp, cosp

            dvpp   = fabs(tee+teo+toe+too)/rho; //+sinp,+cosp
            if (q0 > 0){
                tmaxp  = MAX(1.e-6,q0*elkap*UNIT_c/dvpp);
                dvpp   = (pow(1.+tmaxp,oma)-1.)/tmaxp;  //correct for finite Kappa_max
            }else{
                dvpp   = pow(dvpp/(-q0*elkap*UNIT_c),aCAK);
            }
            //            }
            
            dsym   = 0.;
            iysgn  = +1;
            izsgn  = +1;
            if(kmax == 1){
                dsym=fabs(tee+teo-toe-too)/rho; //-sinp,+cosp,
                iysgn = +1;
                izsgn = -1;
                if (q0 > 0){
                    tmaxm  = MAX(1.e-6,q0*elkap*UNIT_c/dsym);
                    dsym   = (pow(1.+tmaxm,oma)-1.)/tmaxm; //correct for finite Kappa_max
                }else{
                    dsym   = pow(dsym/(-q0*elkap*UNIT_c),aCAK);
                }
            }
            if(jmax == 1){
                dsym=fabs(tee-teo+toe-too)/rho; //+sinp,-cosp
                iysgn = -1;
                izsgn = +1;
                if (q0 > 0){
                    tmaxm  = MAX(1.e-6,q0 * elkap * UNIT_c / dsym);
                    dsym   = (pow(1.+tmaxm,oma)-1.)/tmaxm; //correct for finite Kappa_max
                }else{
                    dsym   = pow(dsym/(-q0*elkap*UNIT_c),aCAK);
                }
            }
            //
            // Add up dv/dl**alpha along this ray.
            //
            
            gii = gii+wtot*(dvpp+dsym)*a1;  
            gjj = gjj+wtot*(dvpp+iysgn*dsym)*a2;
            gkk = gkk+wtot*(dvpp+izsgn*dsym)*a3;
        }
    }
    
    //
    // Normalize forces.
    //
    
    tmp = elkap * qbar * star.luminosity / (oma * 4.0*M_PI * UNIT_c);

    rn   = star.radius;
    if (ifrco == -4){
        rn = xi[data->iminv[jl]];
    }
    
    tmp = tmp/(M_PI*rn*rn);
    
    if (delta != 0.0){
        // dkee 18Jul17 Hardcoded no oblateness
        rsbr     = star.radius / r;
        xmustar  = sqrt(MAX(0.0, 1.0 - rsbr*rsbr));
        dilfac   = 0.5 * (1.0 - xmustar);
        delfac   = pow(data->Vc[RHO][kl][jl][il] / dilfac / (1.e11 / (UNIT_MASS*pow(UNIT_LENGTH,-3)) * UNIT_MASS) / xmbe, delta);
        tmp      = tmp * delfac;
    }
      
    gline[IDIR] = gii*tmp;
    gline[JDIR] = gjj*tmp;
    gline[KDIR] = gkk*tmp;
    
    if ((ifrco == -1)||(ifrco == -2)){
        gline[JDIR]=0.;
    }
    if ((ifrco == -1)||(ifrco == -3)){
        gline[KDIR]=0.;
    }
    return;
}

/* ********************************************************************* */
void sobolev(const Data *data, Grid *grid, Star star, int kl, int jl, int il, double *gline){
    
    int    ip, im, imax;
    double dr, dvrdr, tmp, r;
    double q0, qbar, aCAK, delta, oma, opa;
    double tq0, beta_op, fdfac, dilfac, delfac;
    double rsbr, xmustar;
    double *xi;
    
    const double xhabun = 0.74;
    const double xmpro  = 1.67e-24 / UNIT_MASS;
    const double xmbe   = xmpro * 2.0 / (1.0 + xhabun);
    const double elkap  = 0.2 * (1.0 + xhabun) / (UNIT_MASS*pow(UNIT_LENGTH,-2));
    
    xi   = grid[IDIR].x;
    imax = grid[IDIR].np_tot;
    r    = xi[il];
    
    ip     = MIN(il+1,imax-1);
    im     = MAX(il-1,0);
    dr     = xi[ip]-xi[im];
    dvrdr  = (data->Vc[VX1][kl][jl][ip]-data->Vc[VX1][kl][jl][im])/dr;
    
    aCAK   = g_inputParam[CAK_alpha];
    delta  = g_inputParam[CAK_delta];
    qbar   = g_inputParam[Q_factor];
    q0     = g_inputParam[Q_factor];
    
    oma    = 1.-aCAK;
    opa    = 1.+aCAK;
    
    tmp = elkap * qbar * star.luminosity / (oma * 4.0*M_PI * UNIT_c * r*r);
    
    tq0 = fabs(dvrdr) / (q0 * elkap * data->Vc[RHO][kl][jl][il] * UNIT_c);
    
    if (delta != 0.0){
        rsbr     = star.radius / r;
        xmustar  = sqrt(MAX(0.0, 1.0 - rsbr*rsbr));
        dilfac   = 0.5 * (1.0 - xmustar);
        delfac   = pow(data->Vc[RHO][kl][jl][il] / dilfac / (1.e11 * pow(UNIT_LENGTH,3)) / xmbe, delta);
        tmp      = tmp * delfac;
    }
    
    if (dvrdr != 0.0){
        beta_op = (1.-data->Vc[VX1][kl][jl][il]/(dvrdr*r)) * pow(star.radius/r,2);
        if (beta_op >= 1.){
            fdfac = 1./opa;
        }else if(beta_op < -1.e10){
            fdfac = pow(-beta_op, aCAK) / opa;
        }else if(fabs(beta_op) > 1.e-3){
            fdfac = (1.0 - pow(1.0 - beta_op, opa)) / (beta_op * opa);
        }else{
            fdfac = 1.0 - 0.5 * aCAK * beta_op * (1.0 + 1.0/3.0 * oma * beta_op);
        }
    }else{
        fdfac = 1.;
    }
    
    gline[IDIR] = tmp*pow(tq0,aCAK)*fdfac;
    gline[JDIR] = 0.;
    gline[KDIR] = 0.;

    
#if DEBUGGING > 0
    if(gline[IDIR] < 0){
        printf("ERROR: gline[IDIR] = %e < 0\n", gline[IDIR]);
        printf("ERROR:  see gcak3d.c\n");
    }
    if(isnan(gline[IDIR])){
        printf("ERROR: gline[IDIR] = nan\n");
        printf("ERROR:  see gcak3d.c\n");
    }
#endif
    
    return;
}

/* ********************************************************************* */

// gcak3d setup which only must be done once
//
// Note to Dylan: This was originally written for serial code
//                Currently, this is only doing the ray quadrature so it's fine, but
//                check in future to make sure global grid arrays not expected

void StartUpGCAK(Data *data, Grid *grid, double Omega){
    
    int i,j,k,imax,jmax,kmax,tmpj;
    int ipp,iy,npp1,ny1,npp2,ny2,ifrco;
    double y,dy,pp,dpp;
    double ppmin,ppmax;
    double theta,stheta,thetad;
    double omtmp,wo,arot;
    double ctpm,tmp,tmp2;
    double *tmpx,*tmpy,*tmpz;
    double *xi,*xi_glob,*dxi,*yi,*dyi,*zi,*dzi;

    Star star;
    InitStar1(&star);

    xi      = grid[IDIR].x;
    xi_glob = grid[IDIR].x_glob;
    dxi     = grid[IDIR].dx;
    imax    = grid[IDIR].np_tot;
    yi      = grid[JDIR].x;
    dyi     = grid[JDIR].dx;
    jmax    = grid[JDIR].np_tot;
    zi      = grid[KDIR].x;
    dzi     = grid[KDIR].dx;
    kmax    = grid[KDIR].np_tot;
    
    tmpx = (double*)malloc(imax*sizeof(double));
    tmpy = (double*)malloc(jmax*sizeof(double));
    tmpz = (double*)malloc(kmax*sizeof(double));
    
    ny1  = g_inputParam[CAK3D_nyy];
    npp1 = g_inputParam[CAK3D_npp];
    ny2  = (int) MAX(fabs(g_inputParam[CAK3D_nyy]),1);
    npp2 = (int) MAX(2*fabs(g_inputParam[CAK3D_npp]),1);
    ifrco=g_inputParam[CAK_ifrc];
    
    ppmax =  M_PI;
    ppmin = -M_PI;
    if (kmax == 1) ppmin=0.;
    if (jmax == 1){
        ppmax=M_PI/2.;
        ppmin=-ppmax;
    }
    
//    wo = 0.5*(pow(g_inputParam[winflo],2)*xi[IBEG]/(G_GravityConstant * M_X1_BEG));
//    wo = 0.5*(pow(g_inputParam[winflo],2)*g_inputParam[R_star_CAK]/(G_GravityConstant * M_X1_BEG));
//    arot = wo*pow(1.-wo,2);
    
    for(i = 0; i < imax; i++){
        data->jmaxv[i]=jmax;
    }
    
    for (j = jmax-1;j >= 0; j--){
        data->iminv[j] = 0;
        theta  = yi[j]+0.5*dyi[j];
        stheta = sin(theta);
        thetad = 180.*theta/M_PI;
//        for(i=0; i < imax; i++){
//            tmp = xi[i]/xi[IBEG];
//            tmp = xi[i]/star.radius;
//            tmp2= arot*pow(tmp*stheta,2) +1./tmp-1.;
//            if((tmp2 > 0.)&&(jmax != 1)){
//                data->iminv[j]=i+1;
//                data->jmaxv[i]=j-1;
//            }
//        }
    }
    
    if (ny1 > 0){   //Gauss-Legendre quad
        gauleg(0.,1.,data->y1d,data->wy1d,ny2);
    }else{             //Simpsons Rule  quad
        dy  = 2./ny2;
        y   = -1.+0.5*dy;
        for(iy = 0; iy < ny2; iy++){
            data->y1d[iy]  = y;
            data->wy1d[iy] = dy;
            y        = y+dy;
        }
    }
    if (npp1 > 0){
        gauleg(ppmin,ppmax,data->pp1d,data->wpp1d,npp2);
    }else{
        dpp   =  (ppmax-ppmin)/npp2;
        pp    = ppmin+0.5*dpp;
        for(ipp = 0; ipp < npp2; ipp++){
            data->pp1d[ipp]  = pp;
            data->wpp1d[ipp] = dpp;
            pp         = pp+dpp;
        }
    }
    print1("\n");
    print1("2D or 1.5D case with ifrc = %i & ny,np = %i,%i\n"
           ,ifrco,ny2,npp2);
    print1("y,wy/pp,wpp:\n");
    for(i = 0; i< ny2; i++){
        print1("%i %e %e\n",i,data->y1d[i],data->wy1d[i]);
    }
    for(i = 0; i < npp2; i++){
        print1("%i %e %e\n",i,data->pp1d[i],data->wpp1d[i]);
    }
    
    print1("\n");
    
    for(i = 0; i < imax; i++){
        tmpx[i]=xi[i]+0.5*dxi[i];
    }
    for(j = 0; j < jmax; j++){
        tmpy[j]=yi[j]+0.5*dyi[j];
    }
    for(k = 0; k < kmax; k++){
        tmpz[k]=zi[k]+0.5*dzi[k];
    }
    
    omtmp = Omega;
    if(jmax == 1){ 
        omtmp=0.;
    }
    if (ifrco != -4) {
        print1("ofdwts called\n");
        print1("WARNING: ofdwts not implemented fully yet\n");
        print1("WARNING: no stellar oblateness or gravity darkening\n");
        ofdwts(data, grid, star, omtmp);
    }
    else{
        for(i = 0; i < imax; i++){
            for(j = 0; j < jmax; j++){
                ctpm = sqrt(MAX(0.,1.-pow(star.radius/xi[i],2)));
                for(ipp = 0; ipp < npp2; ipp++){
                    if (ifrco == -4) data->ctpmax[ipp][j][i] = ctpm; // Default ThetaPrime_Max = circ disk approx
                    for(k = 0; k < kmax; k++){
                        for(iy = 0; iy < ny2; iy++){
                            if ((ifrco == -4)||(ifrco == -5)){
                                data->fw[(ipp*npp2)+iy][k][j][i]=1.; // Default flux wt. =uniform disk
                            }
                        }
                    }
                }
            }
        }
    }
}

void gauleg(double x1, double x2, double *x, double *w, int n){
//
//
// Given the lower and upper limits of integration x1 and x2, and given n,
// this routine returns arrays x and w of length n, containing the
// abscissas and weights of the Gauss-Legendre N-point quadrature formula.
// Originated by G. Rybicki; this version adapted from
// Numerical Recipes, Press et al., p. 123.
//

    int i,j,m;
    double xm,xl,z,z1;
    double p1,p2,p3,pp;
    const double errgoal=3.e-14;
//
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//
    m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);
    for(i = 1; i <= m; i++){
        z = cos(M_PI*(i-0.25)/(n+0.5));
        z1= 2.*z;
        while(fabs(z-z1) > errgoal){
            p1=1.;
            p2=0.;
            for(j = 1; j <= n; j++){
                p3=p2;
                p2=p1;
                p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.);
            z1=z;
            z =z1-p1/pp;
        }
        x[i-1] = xm-xl*z;
        x[n-i] = xm+xl*z;
        w[i-1] = 2.*xl/((1.-z*z)*pp*pp);
        w[n-i] = w[i-1];
    }
    return;
}

//----------------------------------------------------------------------------
void ofdwts (Data *data, Grid *grid, Star star, double omfrac){
    
    //
    //  25-JAN-2000:  modified for 3D with spots
    //  27-MAY-1998:  modified to include ifrc= -7 (no gd) and bld (ld coef)
    //  25-APR-1996:  written by SRC, tested standalone.
    //
    //  Computes geometrical extent (ctpmax) of oblate Roche-model star for
    //  field points at various locations in wind.  Also computes integration
    //  weights (fw) that take gravity- and limb-darkening into account.
    //
    //  Inputs:  zxa(imax)  :  radius array of wind points (cm)
    //           zya(jmax)  :  colatitude array of wind points (rad)
    //           zza(jmax)  :  longitude  array of wind points (rad)
    //           pp(npp)    :  phi-prime array of azimuthal field rays (rad)
    //           y(ny)      :  y array of polar field rays (0-1)
    //           iminv(jmax):  locus of radial points just outside star
    //           omfrac     :  Omega / Omega_crit  (fractional ang. velocity)
    //           irayin     :  (0) use actual tpmax, (1) max(tpmax) = pi/2
    //           ifrc       : -6, limb & grav dark; -7, ld, no gd
    //           bld        : limb darkening coef, I(mu)=1+bld*(mu-2/3)
    //
    //  Outputs: ctpmax     :  3D array of cos(max(theta-prime)), i.e. stellar limb
    //           fw         :  5D array of integration weights
    //
    //----------------------------------------------------------------------------
    
    int imax,kmax,jmax;
    int i,j,k,ipp,iy,npp1,ny1,npp2,ny2,ifrco;
    double phimin,phimax,delphi,tmpfw;
    double xeq,ww,w2,sig1,phiobs,tt,stt,ctt,x,xs;
    double tt0,pp0,tpmdum,rpmin0;
    double *xi,*dxi,*yi,*dyi,*zi,*dzi;
    double tmpLimb[2];
    
    const double delsurf=3.0e-4;
    
    xi   = grid[IDIR].x;
    dxi  = grid[IDIR].dx;
    imax = grid[IDIR].np_tot;
    yi   = grid[JDIR].x;
    dyi  = grid[JDIR].dx;
    jmax = grid[JDIR].np_tot;
    zi   = grid[KDIR].x;
    dzi  = grid[KDIR].dx;
    kmax = grid[KDIR].np_tot;
    
    ny1  = g_inputParam[CAK3D_nyy];
    npp1 = g_inputParam[CAK3D_npp];
    ny2  = (int) MAX(fabs(g_inputParam[CAK3D_nyy]),1);
    npp2 = (int) MAX(2*fabs(g_inputParam[CAK3D_npp]),1);
    ifrco = g_inputParam[CAK_ifrc];
    
    phimin = zi[0];
    if (kmax > 1){
        phimax = 2.*zi[kmax-1]-zi[kmax-2];
    }else{
        phimax=2.*M_PI;
    }
    
    delphi = phimax-phimin;
    
    
    //  Compute variables that need to be computed only once: e.g.,
    //  gravity-darkening constant.
    
    xeq  = RX(omfrac,M_PI/2.);
    ww   = omfrac;
    w2   = omfrac*omfrac;
    sig1 = 1.0 - 0.196962*(w2) - 0.0942915*pow(w2,2) + 0.338118*pow(w2,3) - 1.30661*pow(w2,4) + 1.82861*pow(w2,5)- 0.927139*pow(w2,6);
    

    
    //  Begin main loops to compute ctpmax and fw.
    
    
    for(k = 0; k < kmax; k++){
        phiobs = zi[k];
        for(j = 0; j < jmax; j++){
            tt  = yi[j];
            stt = sin(tt);
            ctt = cos(tt);
            xs  = RX(omfrac,tt);
            for(i = data->iminv[j]; i < imax; i++){
//                x = xi[i] / xi[IBEG];
                x = xi[i] / star.radius;
                if (x < xs){
                    for(ipp = 0; ipp < npp2; ipp++){
                        data->ctpmax[ipp][j][i] = 0.;
                        for(iy = 0; iy < ny2; iy++){
                            data->fw[ipp*npp2+iy][k][j][i] = 0.0;
                        }
                    }
                }else if (x == xs){
                    x = x * (1.0 + delsurf);
                }
                for(ipp = 0; ipp < npp2; ipp++){
                    pp0 = data->pp1d[ipp];
                    LimbSearch(pp0, x, tt, ww, xeq, tmpLimb);
                    tpmdum = tmpLimb[0];
                    rpmin0 = tmpLimb[1];

                    data->ctpmax[ipp][j][i] = cos(tpmdum);
                    for(iy = 0; iy < ny2; iy++){
                        if(ifrco > -6){
                            data->fw[ipp*npp2+iy][k][j][i]=1.0;
                        }else{
                            tt0 = TPfunc (tpmdum,data->y1d[iy]);
                            data->fw[ipp*npp2+iy][k][j][i] = Weights(pp0,tt0,tpmdum,rpmin0,x,tt,phiobs,xeq,ww,sig1,phimin,phimax);
                        }
                    }
                }
            }
        }
    }
    return;
}


double Weights (double php, double thp, double thpm, double rrrm, double x, double tt, double phiobs, double xeq, double ww, double sig1, double phimin, double phimax){
    //
    //  Compute limb- and gravity-darkened integration weights.
    //
    //-------------------------------------------------------------------------
    
    // dkee 18Jan16 added pass for tt
    
    double sphp,cphp,sthp,cthp;
    double sth0,cth0,sph0,cph0;
    double sth1,cth1,sph1,cph1;
    double sth2,cth2,sph2,cph2;
    double stt,ctt,w2,grav;
    double ph0sgnd,sph0sgnd;
    double xmupp,Dmupp,bld,ifrco;
    double delphi,ph2;
    double spot,spotwido,spotampo,spotbiaso,spotlato;
    double spotphio;
    double ffdum,ff0,fffirst,ffsecnd,ffnew;
    double rrfirst,rrsecnd,rrmid,rrnew;
    double dnx,dny,dnz,ddd;
    double TH0,PH0,RR0;
    double capR0,TH00,PH00;
    double xReqtol,xtol;
    double gravr,gravt,gravx,gravy,gravz;
    double tmp, tmpfw;
    double tmpCoord[4];
    
    const double ttol=1.0e-6;
    const double small=1.0e-7;
    
    
    ifrco     = g_inputParam[CAK_ifrc];
//    spotwido  = g_inputParam[spotwid];
//    spotampo  = g_inputParam[spotamp];
//    spotbiaso = g_inputParam[spotbias];
//    spotlato  = g_inputParam[spotlat];
//    spotphio  = g_inputParam[spotphi];
    
    bld = 0.;

    if(ifrco <= -6){
        bld = 0.75;
    }

    spotwido  = 0.0;
    spotampo  = 0.0;
    spotbiaso = 0.0;
    spotlato  = 0.0;
    spotphio  = 0.0;
    
    cth1 = cos(spotlato);
    cph1 = cos(spotphio);
    sth1 = sin(spotlato);
    sph1 = sin(spotphio);
    
    sphp  = sin(php);
    cphp  = cos(php);
    sthp  = sin(thp);
    cthp  = cos(thp);
    
    stt = sin(tt);
    ctt = cos(tt);
    w2  = pow(ww,2);
    
    delphi = phimax-phimin;
    
//First check for special cases!
        
    if (thp == 0.0){
        TH0 = tt;
        PH0 = 0.0;
        RR0 = RX(ww,TH0);
    }else if (thp >= thpm){
        CoordTransfm (thp,php,rrrm,x,tt,ww,tmpCoord);
        RR0 = tmpCoord[0];
        TH0 = tmpCoord[1];
        PH0 = tmpCoord[2];
        ffdum = tmpCoord[3];
    }else{
    //  For each ray, find the stellar surface.
                    
        xReqtol = xeq*(1.0+1.0e-4);
        xtol    = xeq*small;
        
        if (x <= xReqtol){
            rrfirst = 0.0;
        }else{
            rrfirst = x - xReqtol;
        }
        rrsecnd = rrrm;
        
        CoordTransfm (thp,php,rrfirst,x,tt,ww,tmpCoord);
        capR0 = tmpCoord[0];
        TH00 = tmpCoord[1];
        PH00 = tmpCoord[2];
        fffirst = tmpCoord[3];
        
        CoordTransfm (thp,php,rrsecnd,x,tt,ww,tmpCoord);
        capR0 = tmpCoord[0];
        TH00 = tmpCoord[1];
        PH00 = tmpCoord[2];
        ffsecnd = tmpCoord[3];
        
//  Now zoom in to find the first intersection point (i.e. stellar surface),
//  and find out what *star-centered* coordinates the point has.
                                
        rrmid = 0.5*(rrfirst+rrsecnd);
        CoordTransfm (thp,php,rrmid,x,tt,ww,tmpCoord);
        RR0 = tmpCoord[0];
        TH0 = tmpCoord[1];
        PH0 = tmpCoord[2];
        ff0 = tmpCoord[3];
        
        rrnew = brent(rrfirst,rrmid,rrsecnd,fffirst,ff0,ffsecnd);
        if ((rrnew < rrfirst)||(rrnew > rrsecnd)){
            rrnew = bisect0(rrfirst,rrmid,rrsecnd,fffirst,ff0,ffsecnd);
        }
        CoordTransfm (thp,php,rrnew,x,tt,ww,tmpCoord);
        RR0 = tmpCoord[0];
        TH0 = tmpCoord[1];
        PH0 = tmpCoord[2];
        ffnew = tmpCoord[3];
        
        while((fabs(rrmid-rrnew) > xtol)&&(ffnew != 0.0)){
                                                    
            if (rrnew < rrmid){
                rrsecnd = rrmid;
                ffsecnd = ff0;
            }else{
                rrfirst = rrmid;
                fffirst = ff0;
            }
            rrmid = rrnew;
            ff0   = ffnew;
            
            rrnew = brent(rrfirst,rrmid,rrsecnd,fffirst,ff0,ffsecnd);
            if ((rrnew < rrfirst)||(rrnew > rrsecnd)){
                rrnew = bisect0(rrfirst,rrmid,rrsecnd,fffirst,ff0,ffsecnd);
            }
            CoordTransfm (thp,php,rrnew,x,tt,ww,tmpCoord);
            RR0 = tmpCoord[0];
            TH0 = tmpCoord[1];
            PH0 = tmpCoord[2];
            ffnew = tmpCoord[3];
            
        }
    }
    
//  NEXT, calculate the projected velocity gradient onto this ray (unit vector n).
                                                                
    
                                                                
    sth0 = sin(TH0);
    cth0 = cos(TH0);
    sph0 = sin(PH0);
    cph0 = cos(PH0);
    
    dnx = x*stt - RR0*sth0*cph0;
    dny =       - RR0*sth0*sph0;
    dnz = x*ctt - RR0*cth0;
    ddd = sqrt(dnx*dnx + dny*dny + dnz*dnz);
    dnx = dnx / ddd;
    dny = dny / ddd;
    dnz = dnz / ddd;
    
//  Calculate atmospheric parameters:  gravity darkening, limb darkening.
                                                                
    gravr = -1.0/RR0/RR0 + 8.0/27.0*RR0*w2*sth0*sth0;
    gravt =                8.0/27.0*RR0*w2*sth0*cth0;
    grav  = sqrt(gravr*gravr + gravt*gravt);
    gravx = gravr*sth0*cph0 + gravt*cth0*cph0;
    gravy = gravr*sth0*sph0 + gravt*cth0*sph0;
    gravz = gravr*cth0      - gravt*sth0;
    
    xmupp = -(gravx*dnx + gravy*dny + gravz*dnz) / grav;
    if (xmupp < 0.0){
        xmupp = 0.0;
    }
    Dmupp = 1.;
    if (bld >= -3.){
        Dmupp = 1.+bld*(xmupp-0.666666667);
    }
    tmpfw    = Dmupp;
    if (ifrco >= -6){
        tmpfw = tmpfw*grav/sig1;
    }
    
    if (spotampo != 0.){
// Exact way to compute angular distance to spot...

        ph0sgnd = sign(fabs(PH0),php);
        sph0sgnd = sin(ph0sgnd);
        ph2  = ph0sgnd-phiobs; // Ensure proper sense of phi' direction...
        while (ph2 < phimin){
            ph2=ph2+delphi;
        }
        while (ph2 > phimax){
            ph2=ph2-delphi;
        }
        tmp = ph2-spotphio;
        if (fabs(tmp+delphi) < fabs(tmp)){
            ph2=ph2+delphi;
        }
        if (fabs(tmp-delphi) < fabs(tmp)){
            ph2=ph2-delphi;
        }
        sph2 = sin(ph2);
        cph2 = cos(ph2);
        sth2 = sth0;
        cth2 = cth0;
        tmp = cth2*cth1+cph2*cph1*sth2*sth1+sph2*sph1*sth2*sth1;  //Dot product
        tmp = acos(tmp);                                          //Convert to radian angle...

//
//Gaussian or Sinusoidal spot ...
//
        spot=0.;
        if(spotwido > 0){        // Gaussian spot
            tmp = tmp*tmp;
            tmp = tmp/(spotwido*spotwido);
            spot= spotampo*exp(-tmp);
            if (bld == -4){
                spot= spot/xmupp;  //Compensate for surface spot area foreshortening...
            }
    //For spotbias
            if(spotbiaso > 0.){
                // >0,prograde bias the spot brightness
                spot= spot*(1.+spotbiaso*sph0sgnd*sth0)/(1.+spotbiaso);
            }
            if((spotbiaso == -1.)&&(php < 0.)){
                // =-1, confine spot to prograde direction!
                spot=0.;
            }
            if (spotbiaso == -2.){   // =-2,  make spot shine only away from its center...
                tmp=phiobs-spotphio;
                if(tmp > delphi/2.){
                    tmp=tmp-delphi;
                }else if(tmp < delphi/2.){
                    tmp=tmp+delphi;
                }
                tmp = tmp*php;
                if(tmp < 0.){
                    spot=0.;
                }
            }
            if (spotbiaso == -3.){
                spot=(spot+1.)*(1.-xmupp)-1.; //  mimic dark filament prominence
            }
            tmpfw  = tmpfw*(1.+spot);
        }else if(spotwido < 0){      // "Sinuspot"
            tmp = -cos(tmp/spotwido);
            tmpfw  = Dmupp*exp(tmp);
        }
    }
    return tmpfw;
}

void LimbSearch(double php, double x, double tt, double ww, double xeq, double *tmpLimb){
//
//  Find the limb of the star for a given pencil of rays specified
//  by an observer point (r,tt) and an observer-centered azimuth (php).
//
//-------------------------------------------------------------------------

    double xsqr,xxmin,xmus,tmaxtry,tmintry,tavgtry,tnew;
    double ffpa,ffpb,ffpc,ffpnew,rrpa,rrpb,rrpc,rrpnew,rrnew;
    double sphp,cphp;
    double thpm,rrrm;
    double tmpRmin[2];
    
    const double ttol = 1.0e-6;
    const double small = 1.0e-7;

    sphp  = sin(php);
    cphp  = cos(php);

//  First bracket the max and min possible values of thp:

    xsqr    = 1.0/x/x;
    xxmin   = MIN(1.0,xsqr);
    xmus    = sqrt( 1.0 - xxmin );
    tmintry = acos(xmus) - small;
    RminSolve (tmintry,x,xeq,tt,php,ww,tmpRmin);
    ffpa = tmpRmin[0];
    rrpa = tmpRmin[1];

    xsqr    = xeq*xeq/x/x;
    xxmin   = MIN(1.0,xsqr);
    xmus    = sqrt( 1.0 - xxmin );
    tmaxtry = acos(xmus) + small;
    RminSolve (tmaxtry,x,xeq,tt,php,ww,tmpRmin);
    ffpc = tmpRmin[0];
    rrpc = tmpRmin[1];

    tavgtry = 0.5*(tmaxtry + tmintry);
    RminSolve (tavgtry,x,xeq,tt,php,ww,tmpRmin);
    ffpb = tmpRmin[0];
    rrpb = tmpRmin[1];

//  Now zoom in on the place where the function f=(r0-R0) just has one
//  solution (i.e. where its minimum is zero AT zero).

    tnew = brent(tmintry,tavgtry,tmaxtry,ffpa,ffpb,ffpc);
    if ((tnew < tmintry)||(tnew > tmaxtry)){
        tnew = bisect0(tmintry,tavgtry,tmaxtry,ffpa,ffpb,ffpc);
    }

    RminSolve (tnew,x,xeq,tt,php,ww,tmpRmin);
    ffpnew = tmpRmin[0];
    rrpnew = tmpRmin[1];

    while ((fabs(tavgtry-tnew) > ttol)&&(ffpnew != 0.0)){
        if (tnew < tavgtry){
            tmaxtry = tavgtry;
            ffpc    = ffpb;
            rrpc    = rrpb;
        }else{
            tmintry = tavgtry;
            ffpa    = ffpb;
            rrpa    = rrpb;
        }
        tavgtry = tnew;
        ffpb    = ffpnew;
        rrpb    = rrpnew;

        tnew = brent(tmintry,tavgtry,tmaxtry,ffpa,ffpb,ffpc);
        if ((tnew < tmintry)||(tnew > tmaxtry)){
            tnew = bisect0(tmintry,tavgtry,tmaxtry,ffpa,ffpb,ffpc);
        }
        RminSolve (tnew,x,xeq,tt,php,ww,tmpRmin);
        ffpnew = tmpRmin[0];
        rrpnew = tmpRmin[1];
    }
    
    tmpLimb[0] = tnew;
    tmpLimb[1] = rrpnew;
    
    return;
}


void RminSolve(double thtry, double x, double xeq, double tt, double php, double ww, double *tmpRmin){
//
//  Finds the minimum of the function (r0-R0), in order to see if our
//  test ray is inside or outside the star.  Returns fff, the value of
//  the function when the derivative is zero, and rrpmin, the radius
//  where this takes place.
//
//-------------------------------------------------------------------------

    double sthp,cthp,xeqtol,xtol;
    double rri,rrf,rrmid,ffpri,ffprf,ffprmid;
    double rrnew,ffprnew,capR0,th0,ph0;
    double fff, rrpmin;
    double tmpCoord[4];

    sthp    = sin(thtry);
    cthp    = cos(thtry);

//  First, determine the min and max bounds over which to search.

    xeqtol = xeq * (1.0+1.0e-4);
    xtol   = xeq * 1.0e-7;

    if (x <= xeqtol){
        rri = 0.0;
        rrf = 2.0*x;
    }else{
        rri = x - xeqtol;
        rrf = x + xeqtol;
    }
    

    ffpri = FuncPrime (thtry,php,rri,tt,ww,x);
    ffprf = FuncPrime (thtry,php,rrf,tt,ww,x);

//  Now zoom in to find the first intersection point (i.e. stellar surface),
//  and find out what *star-centered* coordinates the point has.

    rrmid = 0.5*(rri+rrf);
    ffprmid = FuncPrime (thtry,php,rrmid,tt,ww,x);

    rrnew = brent(rri,rrmid,rrf,ffpri,ffprmid,ffprf);
    if ((rrnew < rri)||(rrnew > rrf)){
        rrnew = bisect0(rri,rrmid,rrf,ffpri,ffprmid,ffprf);
    }

    ffprnew = FuncPrime (thtry,php,rrnew,tt,ww,x);

    while ((fabs(rrmid-rrnew) > xtol)&&(ffprnew != 0.0)){

        if (rrnew < rrmid){
            rrf   = rrmid;
            ffprf = ffprmid;
        }else{
            rri   = rrmid;
            ffpri = ffprmid;
        }
        rrmid   = rrnew;
        ffprmid = ffprnew;

        rrnew = brent(rri,rrmid,rrf,ffpri,ffprmid,ffprf);
        if ((rrnew < rri)||(rrnew > rrf)){
            rrnew = bisect0(rri,rrmid,rrf,ffpri,ffprmid,ffprf);
        }
        ffprnew = FuncPrime (thtry,php,rrnew,tt,ww,x);
    }

    rrpmin  = rrnew;
// CHECK php,x,tt
    CoordTransfm (thtry,php,rrpmin,x,tt,ww,tmpCoord);
    fff = tmpCoord[3];

    tmpRmin[0] = fff;
    tmpRmin[1] = rrpmin;
    
    return;
}

void CoordTransfm (double thp, double php, double rrrp, double x, double tt, double ww, double *tmpCoord){
    //
    //  Transforms from observer-centered coordinates (rrrp,thp,php) to
    //  star-centered coordinates (rr0,t0,p0).
    //
    //-------------------------------------------------------------------------
    
    //dkee 18Jan16 added php to inputs
    
    double xxp,yyp,zzp,rr0,zz0,tmp;
    double xx0,r00;
    double sthp,cthp,sphp,cphp;
    double stt,ctt;
    double capR0, t0, p0, fffp;
    
    sthp = sin(thp);
    cthp = cos(thp);
    sphp = sin(php);
    cphp = cos(php);
    stt  = sin(tt);
    ctt  = cos(tt);
    
    xxp   = rrrp * sthp * cphp;
    yyp   = rrrp * sthp * sphp;
    zzp   = rrrp * cthp;
    
    xx0   = -xxp*ctt - (zzp-x)*stt;
    //     yy0   =  yyp
    zz0   =  xxp*stt - (zzp-x)*ctt;
    r00   = sqrt (xx0*xx0 + yyp*yyp + zz0*zz0);
    tmp   = MIN(1.,zz0/r00);
    //     if(tmp.gt.1.) write(6,*) zz0,r00,tmp
    //     t0    = acos (zz0/r00)
    t0    = acos (tmp);
    if ((xx0 == 0.0)&&(yyp == 0.0)){
        p0  = 0.0;
    }else{
        p0  = acos (MIN(1.,xx0/sqrt(xx0*xx0 + yyp*yyp)));
    }
    capR0 = RX(ww,t0);
    
    fffp  = (r00 - capR0);
    
    tmpCoord[0] = capR0;
    tmpCoord[1] = t0;
    tmpCoord[2] = p0;
    tmpCoord[3] = fffp;
    
    return;
}



double FuncPrime (double thp, double php, double rrrp, double tt, double ww, double x){
//
//  Calculates the derivative with respect to r' of the function (r0-R0)
//
//-------------------------------------------------------------------------

    double xxp,yyp,zzp,xx0,yy0,zz0;
    double dx00,dy00,dz00,r00,dr00;
    double aaa,daa,C0,dC0;
    double sthp,cthp,sphp,cphp,ctt,stt;
    double fprime;
    
    sthp = sin(thp);
    cthp = cos(thp);
    sphp = sin(php);
    cphp = cos(php);
    stt  = sin(tt);
    ctt  = cos(tt);

//  First compute r0 and its derivatives

    xxp   = rrrp * sthp * cphp;
    yyp   = rrrp * sthp * sphp;
    zzp   = rrrp * cthp;
    xx0   = -xxp*ctt - (zzp-x)*stt;
    yy0   =  yyp;
    zz0   =  xxp*stt - (zzp-x)*ctt;
    r00   = sqrt (xx0*xx0 + yy0*yy0 + zz0*zz0);

    if (rrrp != 0.0){
        dz00  = (zz0-x*ctt)/rrrp;
        dr00  = (r00 - (x/r00)*(xx0*stt+zz0*ctt))/rrrp;
    }else{
        dx00  = -sthp*cphp*ctt - cthp*stt;
        dy00  =  sthp*sphp;
        dz00  =  sthp*cphp*stt - cthp*ctt;
        dr00  = (xx0*dx00+yy0*dy00+zz0*dz00)/r00;
    }

//  Next compute R0 and its derivatives
//  CBard 9Dec15, fixed bug in calculation,
//  Specifically changed aaa, daa, and dC0

    aaa   = ww*sqrt(1.0 - zz0*zz0/r00/r00);
    daa   = ww*ww*zz0/r00/r00 * (-dz00 + zz0*dr00/r00);

    if (aaa != 0.0){
        C0  = 3.0/aaa * cos( (M_PI+acos(aaa))/3.0 );
        dC0 = (-C0 + sin((M_PI+acos(aaa))/3.0)/sqrt(1.0-aaa*aaa))/aaa/aaa;
    }else{
        C0  = 1.0;
        dC0 = 0.0;
    }

    fprime = dr00 - (dC0*daa);

    return fprime;
}

    
double brent (double a, double b, double c, double fa, double fb, double fc){
//
//  Given three bounded points of a function, fit a quadratic to these
//  points, and find the point at which the quadratic crosses y=0.
//  (if it cannot find this point, then set xnew outside the input bounds
//  so that this method won't be used at all!)
//
//---------------------------------------------------------------------------
        
// RHDT --
    
    double P,Q,R,S,T;
    double xnew;
        
    if(fc != 0.){
        R = fb / fc;
        T = fa / fc;
    }else{
        R = 0.;
        T = 0.;
    }
        
    if(fa != 0.){
        S = fb / fa;
    }else{
        S = 0.;
    }
    
// -- RHDT
    
    P = S*(T*(R-T)*(c-b) - (1.0-R)*(b-a));
    Q = (T-1.0)*(R-1.0)*(S-1.0);
        
    if (Q != 0.0){
        xnew = b + (P/Q);
    }else{
        xnew = c + 100.0;
    }
        
    return xnew;
}
    

double bisect0 (double a, double b, double c, double fa, double fb, double fc){
//
//  Given three bounded points of a function, find the bounded sub-region,
//  and bisect that sub-region in the hopes of getting closer to the root.
//
//---------------------------------------------------------------------------
        
    double factor;
    double xnew;
        
    factor = fb * (fc-fa);
        
    if (factor > 0.0){
        xnew = 0.5 * (a + b);
    }else{
        xnew = 0.5 * (b + c);
    }
        
    return xnew;
}
    
    
    
double TPfunc (double tmax, double y0){
//
//  Given a theta-prime-max (tmax), and a y ordinate (ranging from 0
//  to 1), compute the intermediary value of theta-prime.
//
//---------------------------------------------------------------------------
        
    double ctmax,ct;
        
//  Map y into mu      (new way that just takes care of AREA)
    
    ctmax  = cos(tmax);
    ct     = 1.0 + (ctmax-1.0)*y0;
    return acos(ct);
    
}

double RX (double w, double theta){
//
//  Calculates normalized radius of a Roche surface given angle theta and
//  fractional angular velocity w.
//
//---------------------------------------------------------------------------
        
    double st,tmp;
        
    tmp = 1.0;
        
    if ((w <= 0.0)||(theta == 0.0)||(theta == M_PI)){
        return tmp;
    }
        
    st = fabs(sin(theta));
    tmp = 3.0 * cos((M_PI+acos(w*st))/3.0) / (w*st);
    
    return tmp;
}

double sign(double A,double B){
    if (B > 0) return A;
    if (B < 0) return -A;
    return 0;
}

#endif
