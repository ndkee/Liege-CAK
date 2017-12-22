#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  
  double ***bx, ***by, ***bz ;
  double ***brtotal, ***bttotal, ***bptotal, ***vphi ;
  double ***p, ***rho, ***vx, ***vy, ***vz;
  double *dx, *dy, *dz;
  double *x, *y, *z;

// Define total magnetic field b0+b1
  
  vphi     = GetUserVar("vphi");
  brtotal = GetUserVar("brtotal");
  bttotal = GetUserVar("bttotal");
  bptotal = GetUserVar("bptotal");
  rho = d->Vc[RHO];  /* pointer shortcut to density    */
  p   = d->Vc[PRS];  /* pointer shortcut to pressure   */
  vx  = d->Vc[VX1];  /* pointer shortcut to x-velocity */
  vy  = d->Vc[VX2];  /* pointer shortcut to y-velocity */
  vz  = d->Vc[VX3];  /* pointer shortcut to y-velocity */
  bx  = d->Vc[BX1];  /* pointer shortcut to x-velocity */
  by  = d->Vc[BX2];  /* pointer shortcut to y-velocity */
  bz  = d->Vc[BX3];  /* pointer shortcut to y-velocity */
  dx = grid[IDIR].dx; /* shortcut to dx */
  dy = grid[JDIR].dx; /* shortcut to dy */
  dz = grid[KDIR].dx; /* shortcut to dz */
  x = grid[IDIR].x; /* shortcut to dx */
  y = grid[JDIR].x; /* shortcut to dy */
  z = grid[KDIR].x; /* shortcut to dz */

  DOM_LOOP(k,j,i){
//  Bx[k][j][i]=B0[0]+V[BX1] 
//  vphi[k][j][i] = (g_OmegaZ * x[i] + vz[k][j][i])*UNIT_VELOCITY; 
  }
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





