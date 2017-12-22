// gcak3d setup which only must be done once
//
// Note to Dylan: This was originally written for serial code
//                Currently, this is only doing the ray quadrature so it's fine, but
//                check in future to make sure global grid arrays not expected

#include <float.h>
#include "pluto.h"
#include "CAKAcceleration.h"
#include "boundary_fluxes.h"
#include "PhysicalConstantsCGS.h"
#include "ReadRestrictionsConfiguration.h"
#include "StellarEvolution.h"

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
//        data->iminv[j]=IBEG;
        data->iminv[j] = 0;
        theta  = yi[j]+0.5*dyi[j];
        stheta = sin(theta);
        thetad = 180.*theta/M_PI;
//        for(i=0; i < imax; i++){
//            tmp = xi[i]/xi[IBEG];
//            tmp = xi[i]/CAKa_R_STAR;
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
        print1("WARNING: only allows circular disk approximation still\n");
        ofdwts(data, grid, omtmp);
    }
    else{
        for(i = 0; i < imax; i++){
            for(j = 0; j < jmax; j++){
                ctpm = sqrt(MAX(0.,1.-pow(CAKa_R_STAR/xi[i],2)));
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
    if (g_inputParam[ABL_disk_tau] == 1){
        findEqCross (data, grid);
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
void ofdwts (Data *data,Grid *grid, double omfrac){
    
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
                x = xi[i] / CAKa_R_STAR;
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
            tmp = -DBL_MAX*cos(tmp/spotwido);
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

void findEqCross (Data *data, Grid *grid){

    FILE   *fptr;
    double roty[3][3], rotz[3][3];
    double rvsapp[3], rvo[3], rvsy[3], rvs[3], musv[3];
    double a, dos, zs, muz, re, phie, dse, hr, phiMin, phiMax, delPhi, s, t;
    double ro, tho, phio, y, phip;
    int    ny, npp, imax, jmax, kmax;
    int    i, j, k, iy, ipp, ire, iphie;
    double *xi, *xi_glob, *yi, *zi, *zi_glob;

    ny  = (int) MAX(fabs(g_inputParam[CAK3D_nyy]),1);
    npp = (int) MAX(2*fabs(g_inputParam[CAK3D_npp]),1);

    xi       = grid[IDIR].x;
    xi_glob  = grid[IDIR].x_glob;
    yi       = grid[JDIR].x;
    zi       = grid[KDIR].x;
    zi_glob  = grid[KDIR].x_glob;

    imax      = grid[IDIR].np_tot;
    jmax      = grid[JDIR].np_tot;
    kmax      = grid[KDIR].np_tot;

#if DIMENSIONS > 2
    phiMin = zi_glob[grid[KDIR].gbeg];
    phiMax = zi_glob[grid[KDIR].gend];
    delPhi = phiMax - phiMin;
#endif

//    printf("grid[KDIR].gbeg = %i \n",grid[KDIR].gbeg);
//    printf("grid[KDIR].gend = %i \n",grid[KDIR].gend);

    for(k = 0; k < kmax; k++){

        phio = zi[k];

        rotz[0][0] = cos(phio);
        rotz[0][1] = -sin(phio);
        rotz[1][0] = sin(phio);
        rotz[1][1] = cos(phio);
        rotz[2][2] = 1.0;

        for(j = 0; j < jmax; j++){

            tho = yi[j];

            roty[0][0] = sin(tho);
            roty[0][2] = -cos(tho);
            roty[1][1] = 1.0;
            roty[2][0] = cos(tho);
            roty[2][2] = sin(tho);

            for(i = 0; i < imax; i++){
                ro = xi[i]/CAKa_R_STAR;
                for(iy = 0; iy < ny; iy++){
                    y = data->y1d[iy];
                    for(ipp = 0; ipp < npp; ipp++){
                        phip = data->pp1d[ipp];
                        if(ro < 1.0){
                            grid->iEqCross[(ipp*npp)+iy][k][j][i] = -1;
                            grid->kEqCross[(ipp*npp)+iy][k][j][i] = -1;
                        }else{
                            rvo[0] = ro*sin(tho)*cos(phio);
                            rvo[1] = ro*sin(tho)*sin(phio);
                            rvo[2] = ro*cos(tho);

                            dos = sqrt(ro*ro - y) - sqrt(1.0 - y);
                            a   = dos*sqrt(y)/ro;

                            rvsapp[0] = sqrt(1.0 - a*a);
                            rvsapp[1] = -a * sin(phip);
                            rvsapp[2] = a * cos(phip);

                            rvsy[0] = roty[0][0]*rvsapp[0] + roty[0][2]*rvsapp[2];
                            rvsy[1] = roty[1][1]*rvsapp[1];
                            rvsy[2] = roty[2][0]*rvsapp[0] + roty[2][2]*rvsapp[2];

                            rvs[0] = rotz[0][0]*rvsy[0] + rotz[0][1]*rvsy[1];
                            rvs[1] = rotz[1][0]*rvsy[0] + rotz[1][1]*rvsy[1];
                            rvs[2] = rotz[2][2]*rvsy[2];

                            musv[0] = (rvo[0]-rvs[0])/dos;
                            musv[1] = (rvo[1]-rvs[1])/dos;
                            musv[2] = (rvo[2]-rvs[2])/dos;

                            dse  = MAX(0.0,-rvs[2]/musv[2]);
                            re   = sqrt(pow(rvs[0] + dse*musv[0],2) + pow(rvs[1] + dse*musv[1],2));
                            phie = atan2((rvs[0] + dse*musv[0]),(rvs[1] + dse*musv[1]));

                            if(re > ro*sin(tho)){
                               re   = ro*sin(tho);
                               phie = phio;
                            }
                            if(re < 1.0){
                               re   = 1.0;
                               phie = atan2(rvs[0],rvs[1]);
                            }
#if DIMENSIONS > 2
                            while(phie < phiMin){
                               phie += delPhi;
                            }
                            while(phie > phiMax){
                               phie -= delPhi;
                            }
#endif
                            grid->rEqCross[ipp*npp+iy][k][j][i]   = re;
                            grid->phiEqCross[ipp*npp+iy][k][j][i] = phie;
                            grid->zsEqCross[ipp*npp+iy][k][j][i]  = rvs[2];
                            grid->musEqCross[ipp*npp+iy][k][j][i] = musv[2];

                            re *= CAKa_R_STAR;

                            ire = grid[IDIR].nghost;
                            while (xi_glob[ire] < re){
                               ire++;
                            }
                            if((xi_glob[ire] - re) < (re - xi_glob[ire-1])){
                               grid->iEqCross[ipp*npp+iy][k][j][i] = MAX(ire - grid[IDIR].nghost,0);
                            }else{
                               grid->iEqCross[ipp*npp+iy][k][j][i] = MAX(ire - 1 - grid[IDIR].nghost,0);
                            }
                            grid->kEqCross[ipp*npp+iy][k][j][i] = 0;
#if DIMENSIONS > 2
                            iphie = grid[KDIR].nghost;
                            while (zi_glob[iphie] < phie){
                               iphie++;
                            }
                            if((zi_glob[iphie] - phie) < (phie - zi_glob[iphie-1])){
                               grid->kEqCross[ipp*npp+iy][k][j][i] = MAX(iphie - grid[KDIR].nghost,0);
                            }else{
                               grid->kEqCross[ipp*npp+iy][k][j][i] = MAX(iphie - 1 - grid[KDIR].nghost,0);
                            }
#endif
                        }
                    }
                }
            }
        }
    }

    return;
}
