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
  double r2, r;
  double *x1 = grid[IDIR].x;                                                  
  double *x2 = grid[JDIR].x;                                                  
  double *x3 = grid[KDIR].x;

  double ***gline1, ***gline2, ***gline3;
  double ***grav;

  grav = GetUserVar("grav");
  gline1 = GetUserVar("gline1");
  gline2 = GetUserVar("gline2");
  gline3 = GetUserVar("gline3");

 /* - Distance from star - */
  r2 = EXPAND(x1[i]*x1[i], + x2[j]*x2[j], + x3[k]*x3[k]);
  r = sqrt(r2);

  /* - Gravity outside bodies - */
 
  DOM_LOOP(k,j,i){

    gline1[k][j][i] = d->gL[0][k][j][i];
    gline2[k][j][i] = d->gL[1][k][j][i];
    gline3[k][j][i] = 0.0;//d->gL[2][k][j][i];

    grav[k][j][i] = (g_inputParam[M_star]*CONST_Msun/UNIT_MASS)/r/r/r; 

 
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





