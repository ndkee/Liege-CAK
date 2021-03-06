/*================================================================================*/
/*
   Initilisation file for a radiativly driven stellar wind with a non-rigid dipole  
   configuration magnetic field.

   The boundary and initial conditions are taken from Runacres and Owocki (2002)    

   The method for calculating the radiative acceleration comes from CAK (1975)      

   The model only works with polar corrdinates in 2D, with the 
   MHD module. 1D, 3D, HD, RHD, RMHD and other geometries do not 
   work at the moment.

*/
/*================================================================================*/
#include "pluto.h"                                                                  

void Init (double *v, double x1, double x2, double x3){
/*================================================================================*/

  double Mratio, Lratio, Bcgs, T, mu, a, b, Q, a_eff, M_star, Edd, eta, Rratio;
  double L, c, M_dot, cs, Bq, v_esc, v_inf, vv, beta, M_dot_cgs, v_inf_cgs;
  double x, y, z, xp, yp, zp, r, theta, Rcgs, omega;
  double br, btheta, bphi, bx, by, bz,  bxp, byp, bzp, rp, rp2;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;

  eta = g_inputParam[Eta];
  Rratio = g_inputParam[R_RATIO];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  omega = g_inputParam[OMEGA];
  T = g_inputParam[TT];
  mu = g_inputParam[MU];
  a = g_inputParam[AA];
  b = g_inputParam[Bb];
  Q = g_inputParam[QQ];
  beta = g_inputParam[BB];
  a_eff = g_inputParam[aa_eff];

  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  L = (Lratio*L_sun/UNIT_L);
  c = 3.0e+5;

  M_dot= (L/(c*c))*(a/(1.0 - a))*pow(Q*Edd/(1.0 - Edd), (1.0 - a)/a);
  M_dot=M_dot*pow(1 + a, -1.0/a);
  M_dot_cgs = M_dot*UNIT_MASS/UNIT_TIME;

  cs = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  v_esc = sqrt(2.0*UNIT_G*M_star*(1.0-Edd));                              
  v_inf = v_esc * sqrt((a/(1.0 - a)));                                      
  v_inf_cgs = v_inf*UNIT_VELOCITY;
  vv = v_inf*pow(1.0 - 1.0/x1, b);                                      


  Bcgs = sqrt(eta*M_dot_cgs*v_inf_cgs/pow(UNIT_LENGTH, 2));
  Bq = Bcgs/UNIT_B;

  g_smallPressure = (v[RHO])*T/(KELVIN*mu);

# if COOLING !=NO
g_minCoolingTemp = g_inputParam[TT]; //Dylan,Asif: this sets minimum T of the sim.
                                     // there is also gmaxCoolingRate that can limit timestep.
                                     // it will be worthwhile to use it as well.
# endif 
#if EOS == IDEAL
  g_gamma = 1.05;
#endif

#if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
#endif

#if ROTATING_FRAME == YES                                                          
 // g_OmegaZ = omega*sqrt((8.0*UNIT_G*M_star)/27.0);                                  
  g_OmegaZ = omega*sqrt((8.0*UNIT_G*M_star)/27.0);                                  
#endif                                                                             

  if(x1 < 1.02 && x2 < CONST_PI/100. && x3 < CONST_PI/100.){
    printf("Bcgs=%e, M_dotcgs=%e, Edd_gam=%e , Omega=%e  \n",Bcgs,M_dot_cgs/6.35e25,Edd, g_OmegaZ*UNIT_VELOCITY);
  }

#if EOS == ISOTHERMAL                                                              
  v[RHO] = (M_dot/(4.0*CONST_PI*vv*x1*x1));                                         
#endif                                                                             

#if EOS == IDEAL                                                              
  v[RHO] = (M_dot/(4.0*CONST_PI*vv*x1*x1));                                         
  v[PRS] = (v[RHO]*T/(KELVIN*mu));                                                  
#endif                                                               

  EXPAND(v[VX1] = vv;,                                                 
         v[VX2] = 0.0;,                               
         v[VX3] = 0.0;)                 
  v[TRC] = 0.0;                                                                     

#if PHYSICS == MHD                                   
#if BACKGROUND_FIELD == NO

  beta *= 0.0174532925;

  // Convert to Cartesian.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);
  // Rotate Cartesian coordiantes.
  xp = x*cos(beta) - z*sin(beta);
  yp = y;
  zp = x*sin(beta) + z*cos(beta);
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  // Calculate b-field components in rotated frame.
  bx = 3.0*xp*zp*Bq*pow(rp,-5);
  by = 3.0*yp*zp*Bq*pow(rp,-5);
  bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);

  // Rotate B-field vector componets.  
  bxp = bx*cos(beta) + bz*sin(beta);
  byp = by;
  bzp = -bx*sin(beta) + bz*cos(beta);

  // Define spherical basis vectors.
  a11 = sin(x2)*cos(x3); a12 = sin(x2)*sin(x3); a13 = cos(x2);
  a21 = cos(x2)*cos(x3); a22 = cos(x2)*sin(x3); a23 = -sin(x2);
  a31 = -sin(x3);        a32 = cos(x3);         a33 = 0.0;

  // Change basis back to spherical polar.
  br = bxp*a11 + byp*a12 + bzp*a13;
  btheta = bxp*a21 + byp*a22 + bzp*a23;
  bphi = bxp*a31 + byp*a32 + bzp*a33;

  EXPAND(v[BX1] = br;,            
         v[BX2] = btheta;,                 
         v[BX3] = bphi;) 

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

}                                                                          

/*================================================================================*/
void Analysis (const Data *d, Grid *grid)
{
}
/*================================================================================*/

/*================================================================================*/
#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)                                                    
{                                                                                                                    
  double Rratio, Lratio, Mratio;
  double M_dot, v_inf, a, Q, Edd, a_ff, L, c;
  double eta, M_dot_cgs, v_inf_cgs, Rcgs;
  double Bq, Bcgs, beta, r, v_esc, a_eff, M_star;
  double x, y, z;
  double xp, yp, zp;
  double theta;
  double br, btheta, bphi, bx, by, bz,  bxp, byp, bzp, rp, rp2;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;

  beta = g_inputParam[BB];
  eta = g_inputParam[Eta];
  Rratio = g_inputParam[R_RATIO];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  a = g_inputParam[AA];
  Q = g_inputParam[QQ];
  a_eff = g_inputParam[aa_eff];

  Rcgs = Rratio*UNIT_LENGTH;
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  L = (Lratio*L_sun/UNIT_L);
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  c = 3.0e+5;
  v_esc = sqrt(2.0*UNIT_G*M_star*(1.0 - Edd));                              
  v_inf = v_esc*sqrt((a/(1.0 - a)));                                      
  v_inf_cgs = v_inf*UNIT_VELOCITY;

  M_dot= (L/(c*c))*(a/(1.0 - a))*pow(Q*Edd/(1.0 - Edd),((1.0 - a)/a));
  M_dot=M_dot*pow(1 + a, -1.0/a);
  M_dot_cgs = M_dot*UNIT_MASS/UNIT_TIME;

  Bcgs = sqrt(eta*M_dot_cgs*v_inf_cgs/pow(UNIT_LENGTH, 2));
  Bq = Bcgs/UNIT_B;

  beta *= 0.0174532925;

  // Convert to Cartesian.
  x = x1*sin(x2)*cos(x3);
  y = x1*sin(x2)*sin(x3);
  z = x1*cos(x2);
  // Rotate Cartesian coordiantes.
  xp = x*cos(beta) - z*sin(beta);
  yp = y;
  zp = x*sin(beta) + z*cos(beta);
  rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
  rp = sqrt(rp2);

  // Calculate b-field components in rotated frame.
  bx = 3.0*xp*zp*Bq*pow(rp,-5);
  by = 3.0*yp*zp*Bq*pow(rp,-5);
  bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);

  // Rotate B-field vector componets.  
  bxp = bx*cos(beta) + bz*sin(beta);
  byp = by;
  bzp = -bx*sin(beta) + bz*cos(beta);

  // Define spherical basis vectors.
  a11 = sin(x2)*cos(x3); a12 = sin(x2)*sin(x3); a13 = cos(x2);
  a21 = cos(x2)*cos(x3); a22 = cos(x2)*sin(x3); a23 = -sin(x2);
  a31 = -sin(x3);        a32 = cos(x3);         a33 = 0.0;

  // Change basis back to spherical polar.
  br = bxp*a11 + byp*a12 + bzp*a13;
  btheta = bxp*a21 + byp*a22 + bzp*a23;
  bphi = bxp*a31 + byp*a32 + bzp*a33;

  printf("br=%e, btheta=%e, bphi=%e \n", br, btheta, bphi);

  EXPAND(B0[0] = br;,
         B0[1] = btheta;,
         B0[2] = bphi;)

}
#endif
/*================================================================================*/

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {        
/*================================================================================*/

  int i, j, k, ip, kp, jp, ghost;

  double h, temp;
  double Idr, Iv1, Iv2,dxi,dxim1;

  double Cs_p, Mratio, Lratio, T, mu, a, b, Q, a_eff, M_star, Edd, Rratio;
  double L, c, M_dot, ke, Omega2, A, Bcgs, cs, eta, Bq;
  double nu2_c, B, sigma, f, gg, beta, Rcgs, vv;
  double x, y, z, xp, yp, zp, r, theta, v_inf, v_esc, v_inf_cgs, M_dot_cgs;
  double vradial, vtheta, vphi;
  double dvdx1, dvdx2, dvdx3;
  double beta_op, opa, oma;
  double br, btheta, bphi, bx, by, bz,  bxp, byp, bzp, rp, rp2;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;

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

  eta = g_inputParam[Eta];
  Rratio = g_inputParam[R_RATIO];
  Cs_p = g_inputParam[Cs_P];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  T = g_inputParam[TT];
  mu = g_inputParam[MU];
  a = g_inputParam[AA];
  opa = 1.0 + a;
  oma = 1.0 - a;
  b = g_inputParam[Bb];
  Q = g_inputParam[QQ];
  a_eff = g_inputParam[aa_eff];

  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));
  L = (Lratio*L_sun/UNIT_L);
  c = 3.0e+5;

  M_dot= (L/(c*c))*(a/(1.0 - a))*pow(Q*Edd/(1.0 - Edd), ((1.0 - a)/a));
  M_dot=M_dot*pow(1.0 + a, -1.0/a);
//  M_dot=M_dot*(1.0+4.*sqrt(1.-a)*cs/a/v_esc); // accounts for sound speed, expected mdot will be higher but makes 
//  no difference for the simulation which controls rho_* 
  M_dot_cgs = M_dot*UNIT_MASS/UNIT_TIME;

  ke = ((4.0*CONST_PI*UNIT_G*M_star*c*Edd)/L);
  Omega2 = pow(0.5,2)*(8.0/27.0)*UNIT_G*M_star;
  A = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c)));
  cs = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
  Bq = Bcgs/UNIT_B;
  v_esc = sqrt(2.0*UNIT_G*M_star*(1.0 - Edd));                              
  v_inf = v_esc*sqrt((a/(1.0 - a)));                                      
  v_inf_cgs = v_inf*UNIT_VELOCITY;

#if EOS == ISOTHERMAL                                                  
  g_isoSoundSpeed = sqrt(UNIT_kB*T/(mu*(CONST_AH/UNIT_MASS)*CONST_amu));
#endif       

  beta = 0.0174532925*g_inputParam[BB];

#if EOS == IDEAL
  g_gamma = 1.05;
#endif

  Bcgs = sqrt(eta*M_dot_cgs*v_inf_cgs/pow(UNIT_LENGTH, 2));
  Bq = Bcgs/UNIT_B;

  ghost = (NX1_TOT - NX1)/2;

  if(side == X1_BEG){                          
    if(box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){ 
  
        rho[k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));          

#if EOS == IDEAL
        prs[k][j][i] = ((rho[k][j][i])*T/(KELVIN*mu));          
#endif



//
#if PHYSICS == MHD
        if (eta < 5000.0) {
#endif
          EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
                 vtheta = 0.0;,
                 vphi = 0.0;)
#if PHYSICS == MHD
        } else if (eta > 5000.0) {
          EXPAND(vradial = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
                 vtheta = 2.0*vx2[k][j][ghost] - vx2[k][j][ghost+1];,
                 vphi = 0.0;)
        }
#endif

        if (vradial > cs){
          EXPAND(vradial = cs;,
                 vtheta = vtheta;,
                 vphi = vphi;)
        } else if (vradial < -cs){
          EXPAND(vradial = -cs;,
                 vtheta = vtheta;,
                 vphi = vphi;)
        }

        if (vtheta > cs){
          EXPAND(vradial = vradial;,
                 vtheta = cs;,
                 vphi = vphi;)
        } else if (vradial < -cs){
          EXPAND(vradial = vradial;,
                 vtheta = -cs;,
                 vphi = vphi;)
        }

        EXPAND(vx1[k][j][i] = vradial;,
               vx2[k][j][i] = vtheta;,
               vx3[k][j][i] = vphi;)
//

        rho[k][j][i] = (M_dot/(4.0*CONST_PI*(cs/Cs_p)));          

#if EOS == IDEAL
        prs[k][j][i] = ((rho[k][j][i])*T/(KELVIN*mu));          
#endif

#if PHYSICS == MHD   
#if BACKGROUND_FIELD == NO

        beta *= 0.0174532925;

        // Convert to Cartesian.
        x = x1[i]*sin(x2[j])*cos(x3[k]);
        y = x1[i]*sin(x2[j])*sin(x3[k]);
        z = x1[i]*cos(x2[j]);
  
        // Rotate Cartesian coordiantes.
        xp = x*cos(beta) - z*sin(beta);
        yp = y;
        zp = x*sin(beta) + z*cos(beta);
        rp2 = EXPAND(xp*xp, + yp*yp, + zp*zp);
        rp = sqrt(rp2);

        // Calculate b-field components in rotated frame.
        bx = 3.0*xp*zp*Bq*pow(rp,-5);
        by = 3.0*yp*zp*Bq*pow(rp,-5);
        bz = (3.0*pow(zp,2)-pow(rp,2))*Bq*pow(rp,-5);

        // Rotate B-field vector componets.  
        bxp = bx*cos(beta) + bz*sin(beta);
        byp = by;
        bzp = -bx*sin(beta) + bz*cos(beta);

        // Define spherical basis vectors.
        a11 = sin(x2[j])*cos(x3[k]); a12 = sin(x2[j])*sin(x3[k]); a13 = cos(x2[j]);
        a21 = cos(x2[j])*cos(x3[k]); a22 = cos(x2[j])*sin(x3[k]); a23 = -sin(x2[j]);
        a31 = -sin(x3[k]);           a32 = cos(x3[k]);            a33 = 0.0;

        // Change basis back to spherical polar.
        br = bxp*a11 + byp*a12 + bzp*a13;
        btheta = bxp*a21 + byp*a22 + bzp*a23;
        bphi = bxp*a31 + byp*a32 + bzp*a33;

        EXPAND(bx1[k][j][i] = br;,            
               bx2[k][j][i] = btheta;,                 
               bx3[k][j][i] = bphi;) 

#endif
#if BACKGROUND_FIELD == YES
        EXPAND(bx1[k][j][i] = 0.0;,
               bx1[k][j][i] = 0.0;,
               bx1[k][j][i] = 0.0;)
#endif
#endif                                                            
      }
    }
  }
                                                                   
  if(side == 0){
    DOM_LOOP(k,j,i){

#if CAK == YES
      dxi = x1[i+1] - x1[i];
      dxim1 = x1[i] - x1[i-1];
      dvdx1 = fabs(-dxi*vx1[k][j][i-1]/(dxim1*(dxi + dxim1)) + (dxi - dxim1)*
                vx1[k][j][i]/(dxi*dxim1) + dxim1*vx1[k][j][i+1]/(dxi*(dxi + dxim1)));

      nu2_c = 1.0 - 1.0/(x1[i]*x1[i]);
      ke = 4.0*CONST_PI*UNIT_G*M_star*c*Edd/L;
      B = rho[k][j][i]*Q*c*ke;
      sigma = (x1[i]/fabs(vx1[k][j][i]))*(dvdx1) - 1.0; 

//      f = ((pow(1.0 + sigma, 1.0 + a) - pow(1.0 + sigma*nu2_c, 1.0 + a))/
//            ((1.0 + a)*(1.0 - nu2_c)*sigma*pow(1.0 + sigma, a)));  

    if (dvdx1 != 0.0){
        beta_op = (1.-vx1[k][j][i]/(dvdx1*x1[i])) * pow(1.0/x1[i],2);
        if (beta_op >= 1.){
            f = 1./opa;
        }else if(beta_op < -1.e10){
            f = pow(-beta_op, a) / opa;
        }else if(fabs(beta_op) > 1.e-3){
            f = (1.0 - pow(1.0 - beta_op, opa)) / (beta_op * opa);
        }else{
            f = 1.0 - 0.5 * a * beta_op * (1.0 + 1.0/3.0 * oma * beta_op);
        }
    }else{
        f = 1.;
    }

      A = ((1.0/(1.0-a))*((ke*L*Q)/(4.0*CONST_PI*c)));

      EXPAND(gLx1[k][j][i] = f*A*pow(x1[i], -2)*pow(dvdx1/B, a);,
             gLx2[k][j][i] = 0.0;,
             gLx3[k][j][i] = 0.0;)

#if EOS == IDEAL
      // Accounting for total ionisation at high temp 
      // and for recombination at low temp.
      temp = prs[k][j][i]*KELVIN*mu/rho[k][j][i];
      gLx1[k][j][i] *= exp(-4.0*log(2.0)*pow((2.0 - temp/T - T/temp), 2));
#endif

      // This if statement acounts for when the gradient 
      // is zero, leading to the acceleration going to -nan 
      // due to sigma.
      if (fabs(dvdx1) < 1.0e-8){
        //printf("x1=%e, dvdx1=%e, f=%e, sigma=%e, nu2_c=%e, gLx1=%e \n", x1, dvdx1, f, sigma, nu2_c, gLx1);
        gLx1[k][j][i] = 0.0;
      }
#endif
      

#if EOS == IDEAL
      if (d->Vc[PRS][k][j][i] < (rho[k][j][i])*T/(KELVIN*mu)){
        d->Vc[PRS][k][j][i] = (rho[k][j][i])*T/(KELVIN*mu);
      }
#endif

    }
  }
}                                                                          
/*================================================================================*/
#if BODY_FORCE != NO
#if CAK == YES
void BodyForceVector(double *v, double *gla, double *g, double x1, double x2, double x3)
{
  double L, A, ke, a, M_star, gg, Edd, Mratio, Lratio, Q, T;
  double h, c, dvdx1, nu2_c, B, sigma, f, gLx1, mu, temp;
  double Idr, Iv1, Iv2,dxi,dxim1;

  c = 3.0e+5;
  a = g_inputParam[AA];
  Q = g_inputParam[QQ];               
  mu = g_inputParam[MU];
  Mratio = g_inputParam[M_RATIO];
  Lratio = g_inputParam[L_RATIO];
  T = g_inputParam[TT];

  M_star = Mratio*CONST_Msun/UNIT_MASS;
  L = Lratio*L_sun/UNIT_L;
  Edd = 2.6e-5*(Lratio)*(1.0/Mratio);
  gg = -UNIT_G*M_star*(1.0 - Edd)/x1/x1;

  if (x1 > 1.0){
    g[IDIR] = gg + gla[0];
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
  } else {
    g[IDIR] = gg;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
  }
  
}
#endif

#if CAK == NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  double M_star, gg, Edd, Mratio, Lratio;

  Mratio = g_inputParam[M_RATIO];
  M_star = (Mratio*CONST_Msun/UNIT_MASS);
  Lratio = g_inputParam[L_RATIO];
  Edd = (2.6e-5*(Lratio)*(1.0/Mratio));

  gg = -UNIT_G*M_star*(1.0 - Edd)/x1/x1;

  g[IDIR] = gg;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif
#endif
/*================================================================================*/

