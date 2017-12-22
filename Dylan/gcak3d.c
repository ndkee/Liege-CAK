#include "pluto.h"
#include "CAKAcceleration.h"
#include "boundary_fluxes.h"
#include "PhysicalConstantsCGS.h"
#include "ReadRestrictionsConfiguration.h"
#include "ReferenceSystem.h"

double CAKa_M_STAR;
double CAKa_L_STAR;
double CAKa_R_STAR;
double CAKa_T_STAR;


//
//---------------------------------------------------------------------------
//
void gcak3d(const Data *data, Grid *grid, int kl, int jl, int il, double *gline){
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
    //
    //      include 'global.h'
    //      include 'zone.h'
    
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
    const double xmbe   = u_AtomicMassUnit * 2.0 / (1.0 + xhabun);
    const double elkap  = 0.2 * (1.0 + xhabun) / ReferenceOpacity;
    
//    print1("In gcak3d\n");
    
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
    
    
    // dkee 20Jun17 make requator always CAKa_R_STAR, HARDCODED NO OFD
    //    if(yi[jmax-1] > 1.1*M_PI/2.){
    //        requator = xi[data->iminv[jmax/2-1]]+0.5*dxi[data->iminv[jmax/2-1]];
    //    }else{
    //        requator = xi[data->iminv[jmax-1]]+0.5*dxi[data->iminv[jmax-1]];
    //    }
    
    requator = CAKa_R_STAR;
    
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
    qbar   = g_inputParam[CAK_qbar];
    q0     = g_inputParam[CAK_q0];
    
    //
    // Keep track of phi rotation with time...
    //
//    dphirot = (g_inputParam[winflo]+g_inputParam[spotdw])*g_time/(requator*dzi[0]); //ASSUMES DPHI=CONST.!!!!
    dphirot = 0.0;
    kdphi   = (int) dphirot;
    wphi2   = dphirot-kdphi;
    wphi1   = 1.-wphi2;
    // Begin Angle and Ray loops:
    //
    kp = kl+1;
    km = kl-1;
    //if(data->Voldv[VX3][kl][jl][il] > 0.){
    //    kp = kl;
    //}else{
    //    kp = kl + 1;
    //}
    //km = kp - 1;
    if (kp >= kmax){
        kp=0;
    }
    if (km < 0){
        km=kmax-1;
    }
    if (kmax == 1){
        //dphi = 1.;
        dphip = 1.;
        dphim = 1.;
    }else{
        //dphi = zi[kp]-zi[km];
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
    //if(data->Voldv[VX2][kl][jl][il] > 0.){
    //    jp = jl;
    //}else{
    //    jp = jl + 1;
    //}
    //jm = jp - 1;
    if (jmax == 1){
        //dtheta = 1.;
        dthetap = 1.;
        dthetam = 1.;
        costo  = 0.;
        sinto  = 1.;
    }else{
        //dtheta = yi[jp]-yi[jm];
        dthetap = yi[jp]-yi[jl];
        dthetam = yi[jl]-yi[jm];
        theta   = yi[jl];
        sinto   = sign(MAX(fabs(sin(theta)),1.e-10),sin(theta));
        costo   = cos(theta);
    }
    cotto = costo/sinto;
    
    //    print1("Before Force Loop\n");
    
    r      = xi[il];
    vr     = data->Voldv[VX1][kl][jl][il];
    vt     = data->Voldv[VX2][kl][jl][il];
    vp     = data->Voldv[VX3][kl][jl][il];
    rho    = data->Vc[RHO][kl][jl][il];
    vrbr   = vr/r;
    vtbr   = vt/r;
    vpbr   = vp/r;
    
    ip     = MIN(il+1,imax-1);
    im     = MAX(il-1,0);
    //if(data->Voldv[VX1][kl][jl][il] > 0.){
    //    ip = il;
    //}else{
    //    ip = il + 1;
    //}
    //im = ip - 1;
    //dr     =  xi[ip]-xi[im];
    //dvrdr  = (data->Voldv[VX1][kl][jl][ip]-data->Voldv[VX1][kl][jl][im])/dr;
    //dvtdr  = (data->Voldv[VX2][kl][jl][ip]-data->Voldv[VX2][kl][jl][im])/dr;
    //dvpdr  = (data->Voldv[VX3][kl][jl][ip]-data->Voldv[VX3][kl][jl][im])/dr;
    //dy     = r*dtheta;
    //dvrdy  = (data->Voldv[VX1][kl][jp][il]-data->Voldv[VX1][kl][jm][il])/dy;
    //dvtdy  = (data->Voldv[VX2][kl][jp][il]-data->Voldv[VX2][kl][jm][il])/dy;
    //dvpdy  = (data->Voldv[VX3][kl][jp][il]-data->Voldv[VX3][kl][jm][il])/dy;
    //dz     = r*dphi*sinto;
    //dvrdz  = (data->Voldv[VX1][kp][jl][il]-data->Voldv[VX1][km][jl][il])/dz;
    //dvtdz  = (data->Voldv[VX2][kp][jl][il]-data->Voldv[VX2][km][jl][il])/dz;
    //dvpdz  = (data->Voldv[VX3][kp][jl][il]-data->Voldv[VX3][km][jl][il])/dz;
    
    drp    =  xi[ip]-xi[il];
    drm    =  xi[il]-xi[im];

    dvrdr  = -drp/(drm * (drp + drm)) * data->Voldv[VX1][kl][jl][im] + (drp - drm)/(drp * drm) * data->Voldv[VX1][kl][jl][il] + drm/(drp * (drp + drm)) * data->Voldv[VX1][kl][jl][ip];
    dvtdr  = -drp/(drm * (drp + drm)) * data->Voldv[VX2][kl][jl][im] + (drp - drm)/(drp * drm) * data->Voldv[VX2][kl][jl][il] + drm/(drp * (drp + drm)) * data->Voldv[VX2][kl][jl][ip];
    dvpdr  = -drp/(drm * (drp + drm)) * data->Voldv[VX3][kl][jl][im] + (drp - drm)/(drp * drm) * data->Voldv[VX3][kl][jl][il] + drm/(drp * (drp + drm)) * data->Voldv[VX3][kl][jl][ip];
    dyp    = r * dthetap;
    dym    = r * dthetam;
    dvrdy  = -dyp/(dym * (dyp + dym)) * data->Voldv[VX1][kl][jm][il] + (dyp - dym)/(dyp * dym) * data->Voldv[VX1][kl][jl][il] + dym/(dyp * (dyp + dym)) * data->Voldv[VX1][kl][jp][il];
    dvtdy  = -dyp/(dym * (dyp + dym)) * data->Voldv[VX2][kl][jm][il] + (dyp - dym)/(dyp * dym) * data->Voldv[VX2][kl][jl][il] + dym/(dyp * (dyp + dym)) * data->Voldv[VX2][kl][jp][il];
    dvpdy  = -dyp/(dym * (dyp + dym)) * data->Voldv[VX3][kl][jm][il] + (dyp - dym)/(dyp * dym) * data->Voldv[VX3][kl][jl][il] + dym/(dyp * (dyp + dym)) * data->Voldv[VX3][kl][jp][il];
    dzp    = r * sinto * dphip;
    dzm    = r * sinto * dphim;
    dvrdz  = -dzp/(dzm * (dzp + dzm)) * data->Voldv[VX1][km][jl][il] + (dzp - dzm)/(dzp * dzm) * data->Voldv[VX1][kl][jl][il] + dzm/(dzp * (dzp + dzm)) * data->Voldv[VX1][kp][jl][il];
    dvtdz  = -dzp/(dzm * (dzp + dzm)) * data->Voldv[VX2][km][jl][il] + (dzp - dzm)/(dzp * dzm) * data->Voldv[VX2][kl][jl][il] + dzm/(dzp * (dzp + dzm)) * data->Voldv[VX2][kp][jl][il];
    dvpdz  = -dzp/(dzm * (dzp + dzm)) * data->Voldv[VX3][km][jl][il] + (dzp - dzm)/(dzp * dzm) * data->Voldv[VX3][kl][jl][il] + dzm/(dzp * (dzp + dzm)) * data->Voldv[VX3][kp][jl][il];

    for(iy = 0; iy < ny2; iy++){
        //        print1("iy = %i\n",iy);
        y    =  data->y1d[iy];
        wy   = data->wy1d[iy];
        for(ipp = 0; ipp < npp2; ipp++){
            //            print1("ipp = %i\n",ipp);
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
            
            //            dvpp = 0.;
            //            if ((tee+teo+toe+too) > 0.){ // no force when vel grad negative
            dvpp   = fabs(tee+teo+toe+too)/rho; //+sinp,+cosp
            //                tmaxp  = MAX(1.e-6,xkmvth/dvpp)
            //                dvpp   = ((1.+tmaxp)**oma-1.)/tmaxp //correct for finite Kappa_max
            if (q0 > 0){
                tmaxp  = MAX(1.e-6,q0*elkap*c_SpeedOfLight/dvpp);
                dvpp   = (pow(1.+tmaxp,oma)-1.)/tmaxp;  //correct for finite Kappa_max
            }else{
                dvpp   = pow(dvpp/(-q0*elkap*c_SpeedOfLight),aCAK);
            }
            //            }
            
            dsym   = 0.;
            iysgn  = +1;
            izsgn  = +1;
            if(kmax == 1){
                //              if((tee+teo-toe-too) > 0.){
                dsym=fabs(tee+teo-toe-too)/rho; //-sinp,+cosp, no force for neg vel grad
                iysgn = +1;
                izsgn = -1;
                if (q0 > 0){
                    tmaxm  = MAX(1.e-6,q0*elkap*c_SpeedOfLight/dsym);
                    dsym   = (pow(1.+tmaxm,oma)-1.)/tmaxm; //correct for finite Kappa_max
                }else{
                    dsym   = pow(dsym/(-q0*elkap*c_SpeedOfLight),aCAK);
                }
                //              }
            }
            if(jmax == 1){
                //              dsym = 0.;
                //             if((tee-teo+toe-too) > 0.){
                dsym=fabs(tee-teo+toe-too)/rho; //+sinp,-cosp
                iysgn = -1;
                izsgn = +1;
                if (q0 > 0){
                    tmaxm  = MAX(1.e-6,q0 * elkap * c_SpeedOfLight / dsym);
                    dsym   = (pow(1.+tmaxm,oma)-1.)/tmaxm; //correct for finite Kappa_max
                }else{
                    dsym   = pow(dsym/(-q0*elkap*c_SpeedOfLight),aCAK);
                }
                //              }
            }
            //
            // Add up dv/dl**alpha along this ray.
            //
            if ((g_inputParam[ABL_disk_tau] == 1) && (xi[il]/CAKa_R_STAR > 1.)){
                tau = rayOptDepth(data, grid, kl, jl, il, iy, ipp);
                //print1("%e %i %i\n",tau,il,jl);
            }else{
                tau=0.;
            }
            
            gii = gii+wtot*(dvpp+dsym)*a1*exp(-tau);  //take account of poss. backstreaming
            gjj = gjj+wtot*(dvpp+iysgn*dsym)*a2*exp(-tau);
            gkk = gkk+wtot*(dvpp+izsgn*dsym)*a3*exp(-tau);
        }
    }
    
    //
    // Normalize forces.
    //
    //       write(6,*) xkmvth,xkmvta
    
    tmp = elkap * qbar * CAKa_L_STAR / (oma * 4.0*M_PI * c_SpeedOfLight);

    rn   = CAKa_R_STAR;
    if (ifrco == -4){
        rn = xi[data->iminv[jl]];
    }
    
    tmp = tmp/(M_PI*rn*rn);
    
    if (delta != 0.0){
        // dkee 18Jul17 Hardcoded no oblateness
        rsbr     = CAKa_R_STAR / r;
        xmustar  = sqrt(MAX(0.0, 1.0 - rsbr*rsbr));
        dilfac   = 0.5 * (1.0 - xmustar);
        delfac   = pow(data->Vc[RHO][kl][jl][il] / dilfac / (1.e11 / ReferenceDensity * ReferenceMass) / xmbe, delta);
        tmp      = tmp * delfac;
    }


//    printf("tmp = %e\n", tmp);
//    printf(" elkap = %e\n", elkap);
//    printf(" qbar = %e\n", qbar);
//    printf(" CAKa_L_STAR = %e\n", CAKa_L_STAR);
//    printf(" oma = %e\n", oma);
//    printf(" r = %e\n", r);
//    printf(" oma = %e\n", oma);
//    printf(" delta = %e\n", delta);
//    printf(" CAKa_R_STAR = %e\n", CAKa_R_STAR);
//    printf(" rsbr = %e\n", rsbr);
//    printf(" xmustar = %e\n", xmustar);
//    printf(" dilfac = %e\n", dilfac);
//    printf(" delfac = %e\n", delfac);
          

    gline[IDIR] = gii*tmp;
    gline[JDIR] = gjj*tmp;
    gline[KDIR] = gkk*tmp;
    
    if ((ifrco == -1)||(ifrco == -2)){
        gline[JDIR]=0.;
    }
    if ((ifrco == -1)||(ifrco == -3)){
        gline[KDIR]=0.;
    }
    //    if((kmax == 1)&&(data->Voldv[VX1][kl][jl][il] < 0.)){
    //        g[KDIR]=0.;
    //    }
    return;
}


void sobolev(const Data *data, Grid *grid, int kl, int jl, int il, double *gline){
    
    int    ip, im, imax;
    double dr, dvrdr, tmp, r;
    double q0, qbar, aCAK, delta, oma, opa;
    double tq0, beta_op, fdfac, dilfac, delfac;
    double rsbr, xmustar;
    double *xi;
    
    const double xhabun = 0.74;
    //    const double xmpro  = 1.67e-24 / ReferenceMass;
    const double xmbe   = u_AtomicMassUnit * 2.0 / (1.0 + xhabun);
    const double elkap  = 0.2 * (1.0 + xhabun) / ReferenceOpacity;
    
    xi   = grid[IDIR].x;
    imax = grid[IDIR].np_tot;
    r    = xi[il];
    
    ip     = MIN(il+1,imax-1);
    im     = MAX(il-1,0);
    dr     = xi[ip]-xi[im];
    dvrdr  = (data->Voldv[VX1][kl][jl][ip]-data->Voldv[VX1][kl][jl][im])/dr;
    
    //#if STELLAREVOLUTION == NO
    aCAK   = g_inputParam[CAK_alpha];
    delta  = g_inputParam[CAK_delta];
    qbar   = g_inputParam[CAK_qbar];
    q0     = g_inputParam[CAK_q0];
    //#else
    //    cakParams(tstar,aCAK,delta,qbar,q0);
    //#endif
    
    oma    = 1.-aCAK;
    opa    = 1.+aCAK;
    
    tmp = elkap * qbar * CAKa_L_STAR / (oma * 4.0*M_PI * c_SpeedOfLight * r*r);
    
    tq0 = fabs(dvrdr) / (q0 * elkap * data->Vc[RHO][kl][jl][il] * c_SpeedOfLight);
    
    if (delta != 0.0){
        rsbr     = CAKa_R_STAR / r;
        xmustar  = sqrt(MAX(0.0, 1.0 - rsbr*rsbr));
        dilfac   = 0.5 * (1.0 - xmustar);
        delfac   = pow(data->Vc[RHO][kl][jl][il] / dilfac / (1.e11 / ReferenceDensity * ReferenceMass) / xmbe, delta);
        tmp      = tmp * delfac;
    }
    
    if (dvrdr != 0.0){
        beta_op = (1.-data->Voldv[VX1][kl][jl][il]/(dvrdr*r)) * pow(CAKa_R_STAR/r,2);
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

//     elkap*qbar*CAKa_L_STAR/(oma*4.0*M_PI*CONST_c/UNIT_VELOCITY*r*r);
//        printf("tmp = %e\n", tmp);
//        printf(" qbar = %e\n", qbar);
//        printf(" q0 = %e\n", q0);
//        printf(" aCAK = %e\n", aCAK);
//        printf(" oma = %e\n", oma);
//        printf(" elkap = %e\n", elkap*ReferenceOpacity);
//        printf(" c_SpeedOfLight = %e\n", c_SpeedOfLight*ReferenceVelocity);
//        printf(" CAKa_L_STAR = %e\n", CAKa_L_STAR*ReferenceLuminosity);
//        printf(" CAKa_R_STAR = %e\n", CAKa_R_STAR*ReferenceLength);
//        printf(" Voldv = %e\n", data->Voldv[VX1][kl][jl][il]*ReferenceVelocity);
//        printf(" RHO = %e\n", data->Vc[RHO][kl][jl][il]*ReferenceDensity);
//        printf(" r = %e\n", r*ReferenceLength);
//        printf(" rsbr = %e\n", rsbr);
//        printf(" xmustar = %e\n", xmustar);
//        printf(" dilfac = %e\n", dilfac);
//        printf(" delfac = %e\n", delfac);

    
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


//void cakParams(double tstar, double aCAK, double delta, double qbar, double q0){
//
//    // Allow for dynamic adjustment of CAK parameters with stellar evolution
//    // Qbar, Q0, alpha from Puls+ 2000
//    // delta from Pauldrach+ 1986
//
//    double q0low, q0high, qbarlow, qbarhigh, alplow, alphigh, dellow, delhigh;
//    double tlow, thigh;
//
//    if(tstar <= 2.e4){
//       q0low    = 14505.;
//       q0high   = 5171.;
//       qbarlow  = 915.;
//       qbarhigh = 1597.;
//       alplow   = 0.44;
//       alphigh  = 0.58;
//       dellow   = 0.02;
//       delhigh  = 0.02;
//       tlow     = 1.e4;
//       thigh    = 2.e4;
//    }else if(tstar <= 3.e4){
//       q0low    = 5171.;
//       q0high   = 3630.;
//       qbarlow  = 1597.;
//       qbarhigh = 2498.;
//       alplow   = 0.58;
//       alphigh  = 0.64;
//       dellow   = 0.02;
//       delhigh  = 0.09;
//       tlow     = 2.e4;
//       thigh    = 3.e4;
//    }else if(tstar < 4.e4){
//       q0low    = 3630.;
//       q0high   = 1778.;
//       qbarlow  = 2498.;
//       qbarhigh = 1954.;
//       alplow   = 0.64;
//       alphigh  = 0.67;
//       dellow   = 0.09;
//       delhigh  = 0.07;
//       tlow     = 3.e4;
//       thigh    = 4.e4;
//    }else{
//       q0low    = 1778.;
//       q0high   = 2260.;
//       qbarlow  = 1954.;
//       qbarhigh = 1939.;
//       alplow   = 1778.;
//       alphigh  = 2260.;
//       dellow   = 0.07;
//       delhigh  = 0.07;
//       tlow     = 4.e4;
//       thigh    = 5.e4;
//    }
//
//    aCAK  = alplow + (alphigh-alplow)/(thigh-tlow)*(tstar-tlow);
//    delta = dellow + (delhigh-dellow)/(thigh-tlow)*(tstar-tlow);
//    qbar  = qbarlow + (qbarhigh-qbarlow)/(thigh-tlow)*(tstar-tlow);
//    q0    = q0low + (q0high-q0low)/(thigh-tlow)*(tstar-tlow);
//
//    return;
//
//}


double rayOptDepth(const Data *data, const Grid *grid, int kl, int jl, int il, int iy, int ipp){
    
    int    npp, imax, kmax, ire, iphie;
    double xo, yo, zo, zs, musz, hr;
    double xe, ye, re, phie, rtmp, phitmp;
    double rho, prs, w, sumw;
    double dos, a;
    double ro, tho, phio, y, phip;
    double tau;
    const double xhabun=0.73;
    const double elkap=0.2*(1.+xhabun) / ReferenceOpacity;

    npp = (int) MAX(2*fabs(g_inputParam[CAK3D_npp]),1);

    ire   = grid->iEqCross[(ipp*npp)+iy][kl][jl][il];
    iphie = grid->kEqCross[(ipp*npp)+iy][kl][jl][il];

    re   = grid->rEqCross[(ipp*npp)+iy][kl][jl][il];
//    phie = grid->phiEqCross[(ipp*npp)+iy][kl][jl][il];
    zs   = grid->zsEqCross[(ipp*npp)+iy][kl][jl][il];
    musz = grid->musEqCross[(ipp*npp)+iy][kl][jl][il];

    imax = grid[IDIR].np_int_glob;
    kmax = grid[KDIR].np_int_glob;

    if(ire < 0){
       return 0.0;
    }

    ro   = grid[IDIR].x[il]/CAKa_R_STAR;
    tho  = grid[JDIR].x[jl];
//    phio = grid[KDIR].x[kl];
//    y    = data->y1d[iy];
//    phip = data->pp1d[ipp];
//    re   = grid[IDIR].x_glob[ire + grid[IDIR].nghost]/CAKa_R_STAR;
   
//    dos = sqrt(ro*ro - y) - sqrt(1.0 - y);
//    a   = dos*sqrt(y)/ro;

//    zs   = sqrt(1.0 - a*a)*cos(tho) + a*cos(phip)*sin(tho);
//    xo   = ro*sin(tho)*cos(phio);
//    yo   = ro*sin(tho)*sin(phio);
    zo   = ro*cos(tho);
//    musz = (zo - zs)/dos;
    
//    t = (re - xi[ire-1])/(xi[ire] - xi[ire-1]);

//    rho=(t*data->Vc[RHO][1][jmax/2][ire] + (1-t)*data->Vc[RHO][1][jmax/2][ire-1] + t*data->Vc[RHO][1][jmax/2+1][ire] + (1-t)*data->Vc[RHO][1][jmax/2+1][ire-1])/2.;

//    rho = 1e-10/ReferenceDensity;

//    rho  = 0.;
//    prs  = 0.;
//    sumw = 0.;

//    xe   = re*cos(phie);
//    ye   = re*sin(phie);

//    for(kl = MAX(0,iphie-3); kl <= MIN(kmax-1,iphie+3); kl++){
//#if DIMENSIONS == 3
//        phitmp = grid[KDIR].x_glob[kl + grid[KDIR].nghost];
//#else
//        phitmp = phie;
//#endif
//        for(il = MAX(0,ire-3); il <= MIN(imax-1,ire+3); il++){
//            rtmp = grid[IDIR].x_glob[il + grid[IDIR].nghost]/CAKa_R_STAR;
//            w    = exp(-4.*log(2.)*(pow(rtmp*cos(phitmp) - xe,2.) + pow(rtmp*sin(phitmp) - ye,2.))/pow(CAKa_R_STAR*0.1,2.));
//            rho += data->rhoEq[(kl*imax) + il]*w;
//#if EOS != ISOTHERMAL
//            prs += data->prsEq[(kl*imax) + il]*w;
//#endif
//            sumw += w;
            //printf("il, kl, rho, w = %i, %i, %e, %e\n",il, kl, data->rhoEq[(kl*imax) + il], w);
//        }
//    }

//    rho /= sumw;
//    prs /= sumw;

    rho = data->rhoEq[(iphie*imax) + ire];

//    printf("zo = %e, zs = %e\n",zo,zs);

    //if(ire == 0) printf("rho = %e\n",rho); 

#if EOS == ISOTHERMAL
    hr = g_isoSoundSpeed/sqrt(G_GravityConstant*CAKa_M_STAR/CAKa_R_STAR)*pow(re,1.5);
#else
//    prs=(t*data->Vc[PRS][1][jmax/2][ire] + (1-t)*data->Vc[PRS][1][jmax/2][ire-1] + t*data->Vc[PRS][1][jmax/2+1][ire] + (1-t)*data->Vc[PRS][1][jmax/2+1][ire-1])/2.;
    prs = data->prsEq[(iphie*imax) + ire];
        
    hr = sqrt(prs/rho)/sqrt(G_GravityConstant*CAKa_M_STAR/CAKa_R_STAR)*pow(re,1.5);
#endif
    
    tau=MAX(1.e-10,(sqrt(M_PI/2.)*rho*elkap*fabs((erfc(-zo/(sqrt(2.)*hr))-erfc(-zs/(sqrt(2.)*hr)))*hr*CAKa_R_STAR/musz)));
    
    return tau;
    
}
