/*     ----------------------------------------------------------------      


                               UNIT SGP4UNIT;


    This file contains the SGP4 procedures. The code was originally          
    released in a 1980 paper. In 1997 and 1998, the updated and combined     
    code (SGP4 and SDP4) was released by NASA on the Internet. This version  
    follows the unrestricted Web version from                                
                 seawifs.gsfc.nasa.gov/~seawifsp/src/bobdays/                

                            Companion code for                               
               Fundamentals of Astrodyanmics and Applications                
                                     2001                                    
                              by David Vallado                               

       (H)               email valladodl@worldnet.att.net                    
       (W) 303-344-6037, email davallado@west.raytheon.com                   

       *****************************************************************     

    Current :                                                                
              14 May 01  David Vallado                                       
                           Original Baseline                                 
    Changes :                                                                
                     97  NASA                                                
                           Internet version                                  
                     80  NORAD                                               
                           Original baseline                                 

       ----------------------------------------------------------------      

                                IMPLEMENTATION

       ----------------------------------------------------------------      */
#include "asttime.h"
#include "ast2body.h"
#include "astutil.h"
#include "sgp4unit.h"

char  Help;
FILE *SGP4File;
char  common[25];

/*-----------------------------------------------------------------------------
*
*                           PROCEDURE DPPER
*
*  This PROCEDURE provides deep space long period periodic contributions
*    to the mean elements.  By design, these periodics are zero at epoch.
*    This used to be DSCOM, but it's really a recurring FUNCTION.
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2         -
*    se3         -
*    Sgh2, Sgh3, Sgh4        -
*    Sh2, Sh3
*    Si2, Si3         -
*    Sl2, Sl3, Sl4         -
*    T           -
*    xgh2        -
*    xgh3        -
*    xgh4        -
*    xh2, xh3    -
*    xi2, xi3    -
*    xl2, xl3, xl4         -
*    zmol        -
*    zmos        -
*    Ep          - Eccentricity                           0.0 - 1.0
*    Inclp       - Inclination
*    Omegap      - Longitude of ascending node
*    Argpp       - Argument of perigee
*    Mp          - Mean Anomaly
*
*  Outputs       :
*    Ep          - Eccentricity                           0.0 - 1.0
*    Inclp       - Inclination
*    Omegap      - Longitude of ascending node
*    Argpp       - Argument of perigee
*    Mp          - Mean Anomaly
*
*  Locals        :
*    alfdp       -
*    betdp       -
*    cosip       -
*    cosop       -
*    dalf        -
*    dbet        -
*    dls         -
*    f2          -
*    f3          -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel         -
*    ses         -
*    sghl        -
*    sghs        -
*    shl         -
*    shs         -
*    sil         -
*    sinip       -
*    sinop       -
*    sinzf       -
*    sis         -
*    sll         -
*    sls         -
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  Coupling      :
*    None.
*
*  References    :
*    NORAD Spacetrack Report #3
*
  ----------------------------------------------------------------------------*/
void DPPer
    (
      Real e3,   Real ee2, Real peo,  Real pgho, Real pho,  Real pinco,
      Real plo,  Real se2, Real se3,  Real sgh2, Real sgh3, Real sgh4,
      Real sh2,  Real sh3, Real si2,  Real si3,  Real sl2,  Real sl3,
      Real sl4,  Real T,   Real xgh2, Real xgh3, Real xgh4, Real xh2,
      Real xh3,  Real xi2, Real xi3,  Real xl2,  Real xl3,  Real xl4,
      Real zmol, Real zmos,
      SINT init,
      Real& Ep, Real& Inclp, Real& Omegap, Real& Argpp, Real& Mp
    )
{
  /* --------------------- Local Variables ------------------------ */
  const Real TwoPi = 2.0 * PI;
  Real alfdp, betdp, cosip, cosop, dalf, dbet, dls,
       f2,    f3,    pe,    pgh,   ph,   pinc, pl ,
       sel,   ses,   sghl,  sghs,  shll, shs,  sil,
       sinip, sinop, sinzf, sis,   sll,  sls,  xls,
       xnoh,  zf,    zm,    zel,   zes,  znl,  zns;

  /* ---------------------- Constants ----------------------------- */
  zns   = 1.19459E-5;
  zes   = 0.01675;
  znl   = 1.5835218E-4;
  zel   = 0.05490;

  /* --------------- CALCULATE TIME VARYING PERIODICS ----------- */
  zm    = zmos + zns * T;
  zf    = zm + 2.0 * zes * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  ses   = se2* f2 + se3 * f3;
  sis   = si2 * f2 + si3 * f3;
  sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
  sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
  zm    = zmol + znl * T;
  zf    = zm + 2.0 * zel * sin(zm);
  sinzf = sin(zf);
  f2    =  0.5 * sinzf * sinzf - 0.25;
  f3    = -0.5 * sinzf * cos(zf);
  sel   = ee2 * f2 + e3 * f3;
  sil   = xi2 * f2 + xi3 * f3;
  sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
  sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
  shll  = xh2 * f2 + xh3 * f3;
  pe    = ses + sel;
  pinc  = sis + sil;
  pl    = sls + sll;
  pgh   = sghs + sghl;
  ph    = shs + shll;

  if (init == 0)
  {
    pe    = pe - peo;
    pinc  = pinc - pinco;
    pl    = pl - plo;
    pgh   = pgh - pgho;
    ph    = ph - pho;
    Inclp = Inclp + pinc;
    Ep    = Ep + pe;
    sinip = sin(Inclp);
    cosip = cos(Inclp);

    /* ----------------- APPLY PERIODICS DIRECTLY ------------ */
    if (Inclp >= 0.2)
    {
      ph     = ph / sinip;
      pgh    = pgh - cosip * ph;
      Argpp  = Argpp + pgh;
      Omegap = Omegap + ph;
      Mp     = Mp + pl;
    }
    else
    {
      /* ---- APPLY PERIODICS WITH LYDDANE MODIFICATION ---- */
      sinop  = sin(Omegap);
      cosop  = cos(Omegap);
      alfdp  = sinip * sinop;
      betdp  = sinip * cosop;
      dalf   =  ph * cosop + pinc * cosip * sinop;
      dbet   = -ph * sinop + pinc * cosip * cosop;
      alfdp  = alfdp + dalf;
      betdp  = betdp + dbet;
      Omegap = Mod(Omegap, TwoPi);
      xls    = Mp + Argpp + cosip * Omegap;
      dls    = pl + pgh - pinc * Omegap * sinip;
      xls    = xls * dls;
      xnoh   = Omegap;
      Omegap = Atan2(alfdp, betdp);
      if (fabs(xnoh - Omegap) > PI)
        if (Omegap < xnoh)
          Omegap = Omegap + TwoPi;
        else
          Omegap = Omegap - TwoPi;
      Mp    = Mp + pl;
      Argpp = xls - Mp - cosip * Omegap;
    }
  }

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%84s\n", 
                        " ------------------After DPPER : ---------------- ");
      fprintf(SGP4File, "    Inputs : \n");
      fprintf(SGP4File, 
              "%7s%13.7f%7s%13.7f%7s%13.7%7s%13.7f%7s%7d%13.7f%7s%13.7f\n",
              "e3", e3, "ee2", ee2, "peo", peo, "pgho", pgho, 
              "pho", pho, "pinco", pinco);
      fprintf(SGP4File,
              "%7s%13.7f%7s%13.7f%7s%13.7%7s%13.7f%7s%7d%13.7f%7s%13.7f\n",
              "plo", plo, "se2", se2, "se3", se3, "sgh2", sgh2, 
              "sgh3", sgh3, "sgh4", sgh4);
      fprintf(SGP4File,
              "%7s%13.7f%7s%13.7f%7s%13.7%7s%13.7f%7s%7d%13.7f%7s%13.7f\n",
              "sh2", sh2, "sh3", sh3, "si2", si2, "si3", si3, 
              "sl2", sl2, "sl3", sl3);
      fprintf(SGP4File,
              "%7s%13.7f%7s%13.7f%7s%13.7%7s%13.7f%7s%7d%13.7f%7s%13.7f\n",
              "sl4", sl4, "T", T, "xgh2", xgh2, "xgh3", xgh3, 
              "xgh4", xgh4, "xh2", xh2);
      fprintf(SGP4File,
              "%7s%13.7f%7s%13.7f%7s%13.7%7s%13.7f%7s%7d%13.7f%7s%13.7f\n",
              "xh3", xh3, "xi2", xi2, "xi3", xi3, "xl2", xl2, 
              "xl3", xl3, "xl4", xl4);
      fprintf(SGP4File,
              "%7s%13.7f%7s%13.7f%7s%13.7\n",
              "zmol", zmol, "zmos", zmos, "Init", init);
      fprintf(SGP4File,  "    In/Out : n");
      fprintf(SGP4File,
              "%7s%13.7f%7s%13.7f%7s%13.7%7s%13.7%7s%13.7ff\n",
          "EP", Ep, "Inclp", Inclp, "Omegap", Omegap, "Argpp", Argpp, "Mp", Mp);
    }
}

/*-----------------------------------------------------------------------------
*
*                           PROCEDURE DSCOM
*
*  This PROCEDURE provides deep space common items used by both the secular
*    and periodics subroutines.  Input is provided as shown. This routine
*    used to be called DPPER, but the functions inside weren't correct.
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    Epoch       -
*    Ep          - Eccentricity
*    Argpp       - Argument of perigee
*    Tc          -
*    Inclp       - Inclination
*    Omegap      - Longitude of Ascending Node
*    Np          - Mean Motion
*
*  Outputs       :
*    Sinim       -
*    Cosim       -
*    Sinomm      -
*    Cosomm      -
*    Snodm       -
*    Cnodm       -
*    Day         -
*    E3          -
*    Ee2         -
*    EM          - Eccentricity
*    EMSq        - Eccentricity squared
*    Gam         -
*    Peo         -
*    PGho        -
*    Pho         -
*    PInco       -
*    Plo         -
*    RTemSq      -
*    Se2, Se3         -
*    Sgh2, Sgh3, Sgh4        -
*    Sh2, Sh3
*    Si2, Si3         -
*    Sl2, Sl3, Sl4         -
*    S1, S2, S3, S4, S5, S6, S7          -
*    SS1, SS2, SS3, SS4, SS5, SS6, SS7         -
*    SZ1, SZ2, SZ3         -
*    SZ11, SZ12, SZ13, SZ21, SZ22, SZ23, SZ31, SZ32, SZ33        -
*    Xgh2, Xgh3, Xgh4        -
*    Xh2, Xh3         -
*    Xi2, Xi3         -
*    Xl2, Xl3, Xl4         -
*    Nm          - Mean Motion
*    Z1, Z2, Z3          -
*    Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33         -
*    ZMol        -
*    ZMos        -
*
*  Locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem        -
*    stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg       -
*    zcosgl      -
*    zcosh       -
*    zcoshl      -
*    zcosi       -
*    zcosil      -
*    zsing       -
*    zsingl      -
*    zsinh       -
*    zsinhl      -
*    zsini       -
*    zsinil      -
*    zx          -
*    zy          -
*
*  Coupling      :
*    None.
*
*  References    :
*    NORAD Spacetrack Report #3
*
  ----------------------------------------------------------------------------*/
void DSCom
    (
      Real Epoch,   Real Ep,     Real Argpp, Real Tc, 
      Real Inclp,   Real Omegap, Real Np,
      Real& SNODM,  Real& CNODM, Real& SINIM,  Real& COSIM, Real& SINOMM, 
      Real& COSOMM, Real& Day,   Real& E3,     Real& Ee2,   Real& EM, 
      Real& EMSQ,   Real& GAM,   Real& Peo,    Real& Pgho,  Real& Pho, 
      Real& PInco,  Real& Plo,   Real& RTemSQ, Real& Se2,   Real& Se3, 
      Real& Sgh2,   Real& Sgh3,  Real& Sgh4,   Real& Sh2,   Real& Sh3, 
      Real& Si2,    Real& Si3,   Real& Sl2,    Real& Sl3,   Real& Sl4,
      Real& S1,     Real& S2,    Real& S3,     Real& S4,    Real& S5, 
      Real& S6,     Real& S7,    Real& SS1,    Real& SS2,   Real& SS3, 
      Real& SS4,    Real& SS5,   Real& SS6,    Real& SS7,   Real& SZ1, 
      Real& SZ2,    Real& SZ3,   Real& SZ11,   Real& SZ12,  Real& SZ13, 
      Real& SZ21,   Real& SZ22,  Real& SZ23,   Real& SZ31,  Real& SZ32, 
      Real& SZ33,   Real& Xgh2,  Real& Xgh3,   Real& Xgh4,  Real& Xh2, 
      Real& Xh3,    Real& Xi2,   Real& Xi3,    Real& Xl2,   Real& Xl3, 
      Real& Xl4,    Real& Nm,    Real& Z1,     Real& Z2,    Real& Z3, 
      Real& Z11,    Real& Z12,   Real& Z13,    Real& Z21,   Real& Z22, 
      Real& Z23,    Real& Z31,   Real& Z32,    Real& Z33,   Real& Zmol, 
      Real& Zmos
    )
{
  /* -------------------------- Constants ------------------------- */
  const Real ZES     =  0.01675;
  const Real ZEL     =  0.05490;
  const Real C1SS    =  2.9864797E-6;
  const Real C1L     =  4.7968065E-7;
  const Real ZSINIS  =  0.39785416;
  const Real ZCOSIS  =  0.91744867;
  const Real ZCOSGS  =  0.1945905;
  const Real ZSINGS  = -0.98088458;
  const Real TwoPi   =  2.0 * PI;

  /* --------------------- Local Variables ------------------------ */
  SINT LsFlg;
  Real A1    , A2    , A3    , A4    , A5    , A6    , A7    ,
       A8    , A9    , A10   , BETASQ, CC    , CTEM  , STEM  ,
       X1    , X2    , X3    , X4    , X5    , X6    , X7    ,
       X8    , XNODCE, XNOI  , ZCOSG , ZCOSGL, ZCOSH , ZCOSHL,
       ZCOSI , ZCOSIL, ZSING , ZSINGL, ZSINH , ZSINHL, ZSINI ,
       ZSINIL, ZX    , ZY;

  Nm     = Np;
  EM     = Ep;
  SNODM  = sin(Omegap);
  CNODM  = cos(Omegap);
  SINOMM = sin(Argpp);
  COSOMM = cos(Argpp);
  SINIM  = sin(Inclp);
  COSIM  = cos(Inclp);
  EMSQ   = EM * EM;
  BETASQ = 1.0 - EMSQ;
  RTemSQ = sqrt(BETASQ);

  /* ----------------- INITIALIZE LUNAR SOLAR TERMS --------------- */
  Peo    = 0.0;
  PInco  = 0.0;
  Plo    = 0.0;
  Pgho   = 0.0;
  Pho    = 0.0;
  Day    = Epoch + 18261.5 + Tc / 1440.0;
  XNODCE = Mod(4.5236020 - 9.2422029E-4 * Day, TwoPi);
  STEM   = sin(XNODCE);
  CTEM   = cos(XNODCE);
  ZCOSIL = 0.91375164 - 0.03568096 * CTEM;
  ZSINIL = sqrt(1.0 - ZCOSIL * ZCOSIL);
  ZSINHL = 0.089683511 * STEM / ZSINIL;
  ZCOSHL = sqrt(1.0 - ZSINHL * ZSINHL);
  GAM    = 5.8351514 + 0.0019443680 * Day;
  ZX     = 0.39785416 * STEM / ZSINIL;
  ZY     = ZCOSHL * CTEM + 0.91744867 * ZSINHL * STEM;
  ZX     = Atan2(ZX, ZY);
  ZX     = GAM + ZX - XNODCE;
  ZCOSGL = cos(ZX);
  ZSINGL = sin(ZX);

  /* ------------------------- DO SOLAR TERMS --------------------- */
  ZCOSG = ZCOSGS;
  ZSING = ZSINGS;
  ZCOSI = ZCOSIS;
  ZSINI = ZSINIS;
  ZCOSH = CNODM;
  ZSINH = SNODM;
  CC    = C1SS;
  XNOI  = 1.0 / Nm;

  for (LsFlg = 1; LsFlg <= 2; LsFlg++)
  {
    A1  =   ZCOSG * ZCOSH + ZSING * ZCOSI * ZSINH;
    A3  =  -ZSING * ZCOSH + ZCOSG * ZCOSI * ZSINH;
    A7  =  -ZCOSG * ZSINH + ZSING * ZCOSI * ZCOSH;
    A8  =   ZSING * ZSINI;
    A9  =   ZSING * ZSINH + ZCOSG * ZCOSI * ZCOSH;
    A10 =   ZCOSG * ZSINI;
    A2  =   COSIM * A7 + SINIM * A8;
    A4  =   COSIM * A9 + SINIM * A10;
    A5  =  -SINIM * A7 + COSIM * A8;
    A6  =  -SINIM * A9 + COSIM * A10;

    X1  =  A1 * COSOMM + A2 * SINOMM;
    X2  =  A3 * COSOMM + A4 * SINOMM;
    X3  = -A1 * SINOMM + A2 * COSOMM;
    X4  = -A3 * SINOMM + A4 * COSOMM;
    X5  =  A5 * SINOMM;
    X6  =  A6 * SINOMM;
    X7  =  A5 * COSOMM;
    X8  =  A6 * COSOMM;

    Z31 = 12.0 * X1 * X1 - 3.0 * X3 * X3;
    Z32 = 24.0 * X1 * X2 - 6.0 * X3 * X4;
    Z33 = 12.0 * X2 * X2 - 3.0 * X4 * X4;
    Z1  =  3.0 *  (A1 * A1 + A2 * A2) + Z31 * EMSQ;
    Z2  =  6.0 *  (A1 * A3 + A2 * A4) + Z32 * EMSQ;
    Z3  =  3.0 *  (A3 * A3 + A4 * A4) + Z33 * EMSQ;
    Z11 = -6.0 * A1 * A5 + EMSQ *  (-24.0 * X1 * X7-6.0 * X3 * X5);
    Z12 = -6.0 *  (A1 * A6 + A3 * A5) + EMSQ * 
           (-24.0 * (X2 * X7 + X1 * X8) - 6.0 * (X3 * X6 + X4 * X5));
    Z13 = -6.0 * A3 * A6 + EMSQ * (24.0 * X2 * X8 - 6.0 * X4 * X6);
    Z21 =  6.0 * A2 * A5 + EMSQ * (24.0 * X1 * X5 - 6.0 * X3 * X7);
    Z22 =  6.0 *  (A4 * A5 + A2 * A6) + EMSQ * 
           (24.0 * (X2 * X5 + X1 * X6) - 6.0 * (X4 * X7 + X3 * X8));
    Z23 =  6.0 * A4 * A6 + EMSQ * (24.0 * X2 * X6 - 6.0 * X4 * X8);
    Z1  = Z1 + Z1 + BETASQ * Z31;
    Z2  = Z2 + Z2 + BETASQ * Z32;
    Z3  = Z3 + Z3 + BETASQ * Z33;
    S3  = CC * XNOI;
    S2  = -0.5 * S3 / RTemSQ;
    S4  = S3 * RTemSQ;
    S1  = -15.0 * EM * S4;
    S5  = X1 * X3 + X2 * X4;
    S6  = X2 * X3 + X1 * X4;
    S7  = X2 * X4 - X1 * X3;

    /* ----------------------- DO LUNAR TERMS ------------------- */
    if (LsFlg == 1)
    {
      SS1   = S1;
      SS2   = S2;
      SS3   = S3;
      SS4   = S4;
      SS5   = S5;
      SS6   = S6;
      SS7   = S7;
      SZ1   = Z1;
      SZ2   = Z2;
      SZ3   = Z3;
      SZ11  = Z11;
      SZ12  = Z12;
      SZ13  = Z13;
      SZ21  = Z21;
      SZ22  = Z22;
      SZ23  = Z23;
      SZ31  = Z31;
      SZ32  = Z32;
      SZ33  = Z33;
      ZCOSG = ZCOSGL;
      ZSING = ZSINGL;
      ZCOSI = ZCOSIL;
      ZSINI = ZSINIL;
      ZCOSH = ZCOSHL * CNODM + ZSINHL * SNODM;
      ZSINH = SNODM * ZCOSHL - CNODM * ZSINHL;
      CC    = C1L;
   }
  }

  Zmol = Mod(4.7199672 + 0.22997150  * Day - GAM, TwoPi);
  Zmos = Mod(6.2565837 + 0.017201977 * Day, TwoPi);

  /* ------------------------ DO SOLAR TERMS ---------------------- */
  Se2  =  2.0 * SS1 * SS6;
  Se3  =  2.0 * SS1 * SS7;
  Si2  =  2.0 * SS2 * SZ12;
  Si3  =  2.0 * SS2 * (SZ13 - SZ11);
  Sl2  = -2.0 * SS3 * SZ2;
  Sl3  = -2.0 * SS3 * (SZ3 - SZ1);
  Sl4  = -2.0 * SS3 * (-21.0 - 9.0 * EMSQ) * ZES;
  Sgh2 =  2.0 * SS4 * SZ32;
  Sgh3 =  2.0 * SS4 * (SZ33 - SZ31);
  Sgh4 =-18.0 * SS4 * ZES;
  Sh2  = -2.0 * SS2 * SZ22;
  Sh3  = -2.0 * SS2 * (SZ23 - SZ21);

  /* ------------------------ DO LUNAR TERMS ---------------------- */
  Ee2  =  2.0 * S1 * S6;
  E3   =  2.0 * S1 * S7;
  Xi2  =  2.0 * S2 * Z12;
  Xi3  =  2.0 * S2 * (Z13 - Z11);
  Xl2  = -2.0 * S3 * Z2;
  Xl3  = -2.0 * S3 * (Z3 - Z1);
  Xl4  = -2.0 * S3 * (21.0 - 9.0 * EMSQ) * ZEL;
  Xgh2 =  2.0 * S4 * Z32;
  Xgh3 =  2.0 * S4 * (Z33 - Z31);
  Xgh4 =-18.0 * S4 * ZEL;
  Xh2  = -2.0 * S2 * Z22;
  Xh3  = -2.0 * S2 * (Z23 - Z21);

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%84s\n", 
          " ------------------ After DSCOM :-----------------");
      fprintf(SGP4File, "    Inputs :\n");
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "EPOCH", Epoch, "Ep", Ep, "Argpp", Argpp, "Tc", Tc, 
          "Inclp", Inclp, "Omegap", Omegap);
      fprintf(SGP4File, "%7s%13.7f\n", "Np", Np);
      fprintf(SGP4File, "    Outputs :\n");
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "SNODM", SNODM, "CNODM", CNODM, "SINIM", SINIM, "COSIM", COSIM, 
          "SINOMM", SINOMM, "COSOMM", COSOMM);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "DAY", Day, "E3", E3, "Ee2", Ee2, "EM", EM, 
          "EMSQ", EMSQ, "GAM", GAM);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Peo", Peo, "Pgho", Pgho, "Pho", Pho, "PInco", PInco, 
          "Plo", Plo);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "RTemSQ", RTemSQ, "Se2", Se2, "Se3", Se3, "Sgh2", Sgh2, 
          "Sgh3", Sgh3, "Sgh4", Sgh4);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Sh2", Sh2, "Sh3", Sh3, "Si2", Si2, "Si3", Si3, 
          "Sl2", Sl2, "Sl3", Sl3);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Sl4", Sl4, "S1", S1, "S2", S2, "S3", S3, 
          "S4", S4, "S5", S5);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "S6", S6, "S7", S7, "SS1", SS1, "SS2", SS2, 
          "SS3", SS3, "SS4", SS4);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "SS5", SS5, "SS6", SS6, "SS7", SS7, "SZ1", SZ1, 
          "SZ2", SZ2, "SZ3", SZ3);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "SZ11", SZ11, "SZ12", SZ12, "SZ13", SZ13, "SZ21", SZ21, 
          "SZ22", SZ22, "SZ23", SZ23);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "SZ31", SZ31, "SZ32", SZ32, "SZ33", SZ33, "Xgh2", Xgh2, 
          "Xgh3", Xgh3, "Xgh4", Xgh4);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Xh2", Xh2, "Xh3", Xh3, "Xi2", Xi2, "Xi3", Xi3, 
          "Xl2", Xl2, "Xl3", Xl3);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Xl4", Xl4, "Nm", Nm, "Z1", Z1, "Z2", Z2, 
          "Z3", Z3, "Z11", Z11);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Z12", Z12, "Z13", Z13, "Z21", Z21, "Z22", Z22, 
          "Z23", Z23, "Z31", Z31);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Z32", Z32, "Z33", Z33, "Zmol", Zmol, "Zmos", Zmos); 
    }
}

/*-----------------------------------------------------------------------------
*
*                           PROCEDURE DSINIT
*
*  This PROCEDURE provides Deep Space contributions to Mean Motion Dot due
*    to geopotential resonance with half day and one day orbits.
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    Cosim       -
*    Emsq        - Eccentricity squared
*    Argpo       - Argument of Perigee
*    S1, S2, S3, S4, S5          -
*    Sinim       -
*    Ss1, Ss2, Ss3, Ss4, Ss5         -
*    Sz1, Sz3
*    Sz11, Sz13, Sz21, Sz23, Sz31, Sz33        -
*    T           - Time
*    Tc          -
*    GSTo        -
*    Mo          - Mean Anomaly
*    MDot        - Mean Anomaly dot (rate)
*    No          - Mean Motion
*    Omegao      - Longitude of ascending node
*    OmegaDot    - Longitude of ascending node dot (rate)
*    XPIDOT      -
*    Z1, Z3      -
*    Z11, Z13, Z21, Z23, Z31, Z33         -
*    EM          - Eccentricity
*    Argpm       - Argument of perigee
*    Inclm       - Inclination
*    Mm          - Mean Anomaly
*    Nm          - Mean Motion
*    Omegam      - Longitude of ascending node
*
*  Outputs       :
*    EM          - Eccentricity
*    Argpm       - Argument of perigee
*    Inclm       - Inclination
*    Mm          - Mean Anomaly
*    Nm          - Mean motion
*    Omegam      - Longitude of ascending node
*    IRez        - Flag foir resonances               1-One day      - 2-Half day
*    Atime       -
*    D2201, D2211, D3210, D3222, D4410, D4422, D5220, D5232, D5421, D5433       -
*    Dedt        -
*    Didt        -
*    DMDT        -
*    DNDT        -
*    DNODT       -
*    DOMDT       -
*    Del1, Del2, Del3        -
*    Ses         -
*    Sghl        -
*    Sghs        -
*    Sgs         -
*    Shl         -
*    Shs         -
*    Sis         -
*    Sls         -
*    THETA       -
*    Xfact       -
*    Xlamo       -
*    Xli         -
*    Xni
*
*  Locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543        -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533        -
*    sini2       -
*    temp        -
*    temp1       -
*    Theta       -
*    xno2        -
*
*  Coupling      :
*    None.
*
*  References    :
*    NORAD Spacetrack Report #3
*
  ----------------------------------------------------------------------------*/
void DSInit
    (
      Real COSIM,  Real EMSQ, Real Argpo, Real S1,   Real S2,     Real S3, 
      Real S4,     Real S5,   Real SINIM, Real SS1,  Real SS2,    Real SS3,
      Real SS4,    Real SS5,  Real SZ1,   Real SZ3,  Real SZ11,   Real SZ13, 
      Real SZ21,   Real SZ23, Real SZ31,  Real SZ33, Real T,      Real TC, 
      Real GSTo,   Real Mo,   Real MDot,  Real No,   Real Omegao, Real OmegaDot,
      Real XPIDOT, Real Z1,   Real Z3,    Real Z11,  Real Z13,    Real Z21, 
      Real Z23,    Real Z31,  Real Z33,   Real Ecco, Real Eccsq,
      Real& EM,    Real& Argpm, Real& Inclm, Real& Mm, Real& Nm, Real& Omegam,
      SINT IRez,
      Real& Atime, Real& D2201, Real& D2211, Real& D3210, Real& D3222, 
      Real& D4410, Real& D4422, Real& D5220, Real& D5232, Real& D5421, 
      Real& D5433, Real& Dedt,  Real& Didt,  Real& Dmdt,  Real& Dndt,
      Real& Dnodt, Real& Domdt, Real& Del1,  Real& Del2,  Real& Del3,
      Real& Xfact, Real& Xlamo, Real& Xli,   Real& Xni
    )
{
  /* --------------------- Local Variables ------------------------ */
  const Real TwoPi = 2.0 * PI;

  Real Ainv2 , AONV  , COSISQ, EOC   , F220  , F221  , F311  ,
       F321  , F322  , F330  , F441  , F442  , F522  , F523  ,
       F542  , F543  , G200  , G201  , G211  , G300  , G310  ,
       G322  , G410  , G422  , G520  , G521  , G532  , G533  ,
       SES   , SGS   , Sghl  , SGHS  , SHS   , Shll  , SIS   ,
       SINI2 , SLS   , Temp  , Temp1 , Theta , Xno2  , Q22   , 
       Q31   , Q33   , Root22, Root44, Root54, RPTim , Root32, 
       Root52, X2o3  , XKe   , ZNL   , emo   , ZNS;

  Q22    = 1.7891679E-6;
  Q31    = 2.1460748E-6;
  Q33    = 2.2123015E-7;
  Root22 = 1.7891679E-6;
  Root44 = 7.3636953E-9;
  Root54 = 2.1765803E-9;
  RPTim  = 4.37526908801129966E-3;
  Root32 = 3.7393792E-7;
  Root52 = 1.1428639E-7;
  X2o3   = 2.0 / 3.0;
  XKe    = 7.43669161331734132E-2;
  ZNL    = 1.5835218E-4;
  ZNS    = 1.19459E-5;

  /* -------------------- Deep Space Initialization ------------ */
  IRez = 0;
  if ((Nm < 0.0052359877) && (Nm > 0.0034906585))
    IRez = 1;
  if ((Nm >= 8.26E-3) && (Nm <= 9.24E-3) && (EM >= 0.5))
    IRez = 2;

  /* ------------------------ DO SOLAR TERMS ------------------- */
  SES  =  SS1 * ZNS * SS5;
  SIS  =  SS2 * ZNS * (SZ11 + SZ13);
  SLS  = -ZNS * SS3 * (SZ1 + SZ3 - 14.0 - 6.0 * EMSQ);
  SGHS =  SS4 * ZNS * (SZ31 + SZ33 - 6.0);
  SHS  = -ZNS * SS2 * (SZ21 + SZ23);
  if (Inclm < 5.2359877E-2)
    SHS = 0.0;
  if (SINIM != 0.0)
    SHS = SHS / SINIM;
  SGS  = SGHS - COSIM * SHS;

  /* ------------------------- DO LUNAR TERMS ------------------ */
  Dedt = SES + S1 * ZNL * S5;
  Didt = SIS + S2 * ZNL * (Z11 + Z13);
  Dmdt = SLS - ZNL * S3 * (Z1 + Z3 - 14.0 - 6.0 * EMSQ);
  Sghl = S4 * ZNL * (Z31 + Z33 - 6.0);
  Shll = -ZNL * S2 * (Z21 + Z23);
  if (Inclm < 5.2359877E-2)
    Shll = 0.0;
  Domdt = SGS + Sghl;
  Dnodt = SHS;
  if (SINIM != 0.0)
  {
    Domdt = Domdt - COSIM / SINIM * Shll;
    Dnodt = Dnodt + Shll / SINIM;
  }

  /* ----------- CALCULATE DEEP SPACE RESONANCE EFFECTS -------- */
  Dndt   = 0.0;
  Theta  = Mod(GSTo + TC * RPTim, TwoPi);
  EM     = EM + Dedt * T;
  // shouldn't emsq be changed now?????
  Inclm  = Inclm + Didt * T;
  Argpm  = Argpm + Domdt * T;
  Omegam = Omegam + Dnodt * T;
  Mm     = Mm + Dmdt * T;
  if (Inclm < 0.0)
  {
    Inclm  = -Inclm;
    Argpm  = Argpm - PI;
    Omegam = Omegam + PI;
  }

  /* -------------- Initialize the resonance terms ------------- */
  if (IRez != 0)
    AONV = Power(Nm / XKe, X2o3);

  /* ---------- GEOPOTENTIAL RESONANCE FOR 12 HOUR ORBITS ------ */
  if (IRez == 2)
  {
    COSISQ = COSIM * COSIM;
    emo    = EM;
    EM     = Ecco;
    EMSQ   = Eccsq;
    EOC    = EM * EMSQ;
    G201   = -0.306 - (EM - 0.64) * 0.440;

    if (EM <= 0.65)
    {
      G211 =    3.616  -  13.2470 * EM +  16.2900 * EMSQ;
      G310 =  -19.302  + 117.3900 * EM - 228.4190 * EMSQ +  156.5910 * EOC;
      G322 =  -18.9068 + 109.7927 * EM - 214.6334 * EMSQ +  146.5816 * EOC;
      G410 =  -41.122  + 242.6940 * EM - 471.0940 * EMSQ +  313.9530 * EOC;
      G422 = -146.407  + 841.8800 * EM - 1629.014 * EMSQ + 1083.4350 * EOC;
      G520 = -532.114  + 3017.977 * EM - 5740.032 * EMSQ + 3708.2760 * EOC;
    }
    else
    {
      G211 =   -72.099 +   331.819 * EM -   508.738 * EMSQ +   266.724 * EOC;
      G310 =  -346.844 +  1582.851 * EM -  2415.925 * EMSQ +  1246.113 * EOC;
      G322 =  -342.585 +  1554.908 * EM -  2366.899 * EMSQ +  1215.972 * EOC;
      G410 = -1052.797 +  4758.686 * EM -  7193.992 * EMSQ +  3651.957 * EOC;
      G422 = -3581.690 + 16178.110 * EM - 24462.770 * EMSQ + 12422.520 * EOC;
      if (EM > 0.715)
        G520 =-5149.66 + 29936.92 * EM - 54087.36 * EMSQ + 31324.56 * EOC;
      else
        G520 = 1464.74 -  4664.75 * EM +  3763.64 * EMSQ;
    }
    if (EM < 0.7)
    {
      G533 = -919.22770 + 4988.6100 * EM - 9064.7700 * EMSQ + 5542.21  * EOC;
      G521 = -822.71072 + 4568.6173 * EM - 8491.4146 * EMSQ + 5337.524 * EOC;
      G532 = -853.66600 + 4690.2500 * EM - 8624.7700 * EMSQ + 5341.4  * EOC;
    }
    else
    {
      G533 =-37995.780 + 161616.52 * EM - 229838.20 * EMSQ + 109377.94 * EOC;
      G521 =-51752.104 + 218913.95 * EM - 309468.16 * EMSQ + 146349.42 * EOC;
      G532 =-40023.880 + 170470.89 * EM - 242699.48 * EMSQ + 115605.82 * EOC;
    }

    SINI2 =  SINIM * SINIM;
    F220 =  0.75 * (1.0 + 2.0 * COSIM+COSISQ);
    F221 =  1.5 * SINI2;
    F321 =  1.875 * SINIM  *  (1.0 - 2.0 * COSIM - 3.0 * COSISQ);
    F322 = -1.875 * SINIM  *  (1.0 + 2.0 * COSIM - 3.0 * COSISQ);
    F441 = 35.0 * SINI2 * F220;
    F442 = 39.3750 * SINI2 * SINI2;
    F522 =  9.84375 * SINIM * (SINI2 * (1.0 - 2.0 * COSIM- 5.0 * COSISQ) + 
            0.33333333 * (-2.0 + 4.0 * COSIM + 6.0 * COSISQ) );
    F523 = SINIM * (4.92187512 * SINI2 * (-2.0 - 4.0 * COSIM +
           10.0 * COSISQ) + 6.56250012 * (1.0+2.0 * COSIM - 3.0 * COSISQ));
    F542 = 29.53125 * SINIM * (2.0 - 8.0 * COSIM+COSISQ *
           (-12.0 + 8.0 * COSIM + 10.0 * COSISQ));
    F543 = 29.53125 * SINIM * (-2.0 - 8.0 * COSIM+COSISQ *
           (12.0 + 8.0 * COSIM - 10.0 * COSISQ));
    Xno2  =  Nm * Nm;
    Ainv2 =  AONV * AONV;
    Temp1 =  3.0 * Xno2 * Ainv2;
    Temp  =  Temp1 * Root22;
    D2201 =  Temp * F220 * G201;
    D2211 =  Temp * F221 * G211;
    Temp1 =  Temp1 * AONV;
    Temp  =  Temp1 * Root32;
    D3210 =  Temp * F321 * G310;
    D3222 =  Temp * F322 * G322;
    Temp1 =  Temp1 * AONV;
    Temp  =  2.0 * Temp1 * Root44;
    D4410 =  Temp * F441 * G410;
    D4422 =  Temp * F442 * G422;
    Temp1 =  Temp1 * AONV;
    Temp  =  Temp1 * Root52;
    D5220 =  Temp * F522 * G520;
    D5232 =  Temp * F523 * G532;
    Temp  =  2.0 * Temp1 * Root54;
    D5421 =  Temp * F542 * G521;
    D5433 =  Temp * F543 * G533;
    Xlamo =  Mod(Mo + Omegao + Omegao-Theta - Theta, TwoPi);
    Xfact =  MDot + Dmdt + 2.0 * (OmegaDot + Dnodt - RPTim) - No;
    EM    = emo;
  }

  /* ---------------- SYNCHRONOUS RESONANCE TERMS -------------- */
  if (IRez == 1)
  {
    G200  = 1.0 + EMSQ * (-2.5 + 0.8125 * EMSQ);
    G310  = 1.0 + 2.0 * EMSQ;
    G300  = 1.0 + EMSQ * (-6.0 + 6.60937 * EMSQ);
    F220  = 0.75 * (1.0 + COSIM) * (1.0 + COSIM);
    F311  = 0.9375 * SINIM * SINIM * (1.0 + 3.0 * COSIM) - 0.75 * (1.0 + COSIM);
    F330  = 1.0 + COSIM;
    F330  = 1.875 * F330 * F330 * F330;
    Del1  = 3.0 * Nm * Nm * AONV * AONV;
    Del2  = 2.0 * Del1 * F220 * G200 * Q22;
    Del3  = 3.0 * Del1 * F330 * G300 * Q33 * AONV;
    Del1  = Del1 * F311 * G310 * Q31 * AONV;
    Xlamo = Mod(Mo + Omegao + Argpo - Theta, TwoPi);
    Xfact = MDot + XPIDOT - RPTim + Dmdt + Domdt + Dnodt - No;
  }

  /* ------------ FOR SGP4, INITIALIZE THE INTEGRATOR ---------- */
  if (IRez != 0)
  {
    Xli   = Xlamo;
    Xni   = No;
    Atime = 0.0;
    Nm    = No + Dndt;
  }

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%84s\n",
          " ------------------After DSINIT : ---------------");
      fprintf(SGP4File, "    Inputs :\n");
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "Cosim", COSIM, "Emsq", EMSQ, "Argpo", Argpo, 
          "S1", S1, "S2", S2, "S3", S3);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "S4", S4, "S5", S5, "Sinim", SINIM, 
          "Ss1", SS1, "Ss2", SS2, "Ss3", SS3);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "Ss4",  SS4, "Ss5", SS5, "Sz1", SZ1, "Sz3", SZ3, "Sz11", SZ11, 
          "Sz13", SZ13);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "GSTo",     GSTo, "Mo", Mo, "MDot", MDot, "No", No, "Omegao", Omegao,
          "OmegaDot", OmegaDot);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "XPIDOT", XPIDOT, "Z1", Z1, "Z3", Z3, "Z11", Z11, "Z13", Z13, 
          "Z21", Z21);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Z23", Z23, "Z31", Z31, "Z33", Z33);
      fprintf(SGP4File, "    In / Out :\n");
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "EM", EM, "Argpm", Argpm, "Inclm", Inclm, "Mm", Mm, "Nm", Nm, 
          "Omegam", Omegam);
      fprintf(SGP4File, "    Outputs :\n");
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "IREZ", IRez, "Atime", Atime, "D2201", D2201, "D2211", D2211, 
          "D3210", D3210, "D3222", D3222);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "D4410", D4410, "D4422", D4422, "D5220", D5220, "D5232", D5232, 
          "D5421", D5421, "D5433", D5433);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "Dedt", Dedt, "Didt", Didt, "DMDT", Dmdt, "DNDT", Dndt, 
          "DNODT", Dnodt, "DOMDT", Domdt);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s&13.7f%7s%13.7f\n",
          "Del1", Del1, "Del2", Del2, "Del3", Del3, "Xfact", Xfact, 
          "Xlamo", Xlamo, "Xli", Xli);
      fprintf(SGP4File, "%7s%13.7f\n", "Xni", Xni);
    }
}

/*-----------------------------------------------------------------------------
*
*                           PROCEDURE DSPACE
*
*  This PROCEDURE provides deep space contributions to mean elements for
*    perturbing third body.  These effects ave been averaged over one
*    revolution of the sun and moon.  For Earth resonance effects, the
*    effects have been averaged over No revolutions of the satellite. (Mean motion)
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    D2201, D2211, D3210, D3222, D4410, D4422, D5220, D5232, D5421, D5433       -
*    Dedt        -
*    Del1        -
*    Del2        -
*    Del3        -
*    Didt        -
*    Dmdt        -
*    Dnodt       -
*    Domdt       -
*    IRez        -
*    Argpo       - Argument of perigee
*    ArgpDot     - Argument of perigee dot (rate)
*    T           - Time
*    TC          -
*    GSTo        - GST
*    Xfact       -
*    Xlamo       -
*    No          - Mean Motion
*    Atime       -
*    EM          - Eccentricity
*    Ft          -
*    Argpm       - Argument of perigee
*    Inclm       - Inclination
*    Xli         -
*    Mm          - Mean Anomaly
*    Xni         - Mean Motion
*    Omegam      - Longitude of ascending node
*
*  Outputs       :
*    Atime       -
*    EM          - Eccentricity
*    Argpm       - Argument of perigee
*    Inclm       - Inclination
*    xli         -
*    Mm          - Mean Anomaly
*    xni         -
*    Omegam      - Longitude of ascending node
*    Dndt        -
*    Nm          - Mean Motion
*
*  Locals        :
*    delt        -
*    Ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  Coupling      :
*    SREZ        -
*
*  References    :
*    NORAD Spacetrack Report #3
*
  ----------------------------------------------------------------------------*/
void DSpace
    (
      SINT IRez,
      Real D2201, Real D2211, Real D3210, Real D3222, Real D4410, Real D4422, 
      Real D5220, Real D5232, Real D5421, Real D5433, Real Dedt,  Real Del1,
      Real Del2,  Real Del3,  Real Didt,  Real Dmdt,  Real Dnodt, Real Domdt,
      Real Argpo, Real ArgpDot, Real T,   Real TC,    Real GSTo,  Real Xfact,
      Real Xlamo, Real No,
      Real& Atime, Real& EM,    Real& Argpm, Real& Inclm, Real& Xli, Real& Mm,
      Real& XNi,  Real& Omegam, Real& Dndt,  Real& Nm
    )
{
  const Real TwoPi = 2.0 * PI;

  SINT IRetn , IRet;
  Real Delt, Ft, Theta, X2li, X2omi, Xl, Xldot , Xnddt, Xndt, Xomi, G22, G32, 
       G44, G52, G54, FASX2, FASX4, FASX6, RPtim , Step2, Stepn , Stepp;

  FASX2 = 0.13130908;
  FASX4 = 2.8843198;
  FASX6 = 0.37448087;
  G22   = 5.7686396;
  G32   = 0.95240898;
  G44   = 1.8014998;
  G52   = 1.0508330;
  G54   = 4.4108898;
  RPtim = 4.37526908801129966E-3;
  Stepp =    720.0;
  Stepn =   -720.0;
  Step2 = 259200.0;

  /* ----------- CALCULATE DEEP SPACE RESONANCE EFFECTS ----------- */
  Dndt   = 0.0;
  Theta  = Mod(GSTo + TC * RPtim, TwoPi);
  EM     = EM + Dedt * T;
// shouldn't emsq be changed now?????
  Inclm  = Inclm + Didt * T;
  Argpm  = Argpm + Domdt * T;
  Omegam = Omegam + Dnodt * T;
  Mm     = Mm + Dmdt * T;
  if (Inclm < 0.0)
  {
    Inclm  = -Inclm;
    Argpm  = Argpm - PI;
    Omegam = Omegam + PI;
  }

  /* - UPDATE RESONANCES : NUMERICAL (EULER-MACLAURIN) INTEGRATION - */
  /* ------------------------- EPOCH RESTART ----------------------  */
  if (IRez != 0)
  {
    if ((Atime = 0.0) || ((T >= 0.0) && (Atime < 0.0)) ||
        ((T < 0.0) && (Atime >= 0.0)))
    {
      if (T >= 0.0)
        Delt = Stepp;
      else
        Delt = Stepn;
    }
    IRetn = 381; // added for do loop
    IRet  =   0; // added for loop
    while (IRetn == 381)
    {
      if ((fabs(T) < fabs(Atime)) || (IRet == 351))
      {
        if (T >= 0.0)
          Delt = Stepn;
        else
          Delt = Stepp;
        IRet  = 351;
        IRetn = 381;
      }
      else
      {
        if (T > 0.0)  // Error IF prev IF has Atime:=0.0 and t:=0.0 (ge)
          Delt = Stepp;
        else
          Delt = Stepn;
        if (fabs(T - Atime) >= Stepp)
        {
          IRet  = 0;
          IRetn = 381;
        }
        else
        {
          Ft    = T - Atime;
          IRetn = 0;
        }
      }

      /* ------------------- DOT TERMS CALCULATED ------------- */
      /* ----------- NEAR - SYNCHRONOUS RESONANCE TERMS ------- */
      if (IRez != 2)
      {
        Xndt  = Del1 * sin(Xli - FASX2) + Del2 * sin(2.0 * (Xli - FASX4)) +
                Del3 * sin(3.0 * (Xli - FASX6));
        Xldot = XNi + Xfact;
        Xnddt = Del1 * cos(Xli - FASX2) + 
                2.0 * Del2 * cos(2.0 * (Xli - FASX4)) +
                3.0 * Del3 * cos(3.0 * (Xli - FASX6));
        Xnddt = Xnddt * Xldot;
      }
      else
      {
        /* --------- NEAR - HALF-DAY RESONANCE TERMS -------- */
        Xomi  = Argpo + ArgpDot * Atime;
        X2omi = Xomi + Xomi;
        X2li  = Xli + Xli;
        Xndt  = D2201 * sin(X2omi + Xli - G22) + D2211 * sin(Xli - G22) +
              D3210 * sin(Xomi + Xli - G32)  + D3222 * sin(-Xomi + Xli - G32)+
              D4410 * sin(X2omi + X2li - G44)+ D4422 * sin(X2li - G44) +
              D5220 * sin(Xomi + Xli - G52)  + D5232 * sin(-Xomi + Xli - G52)+
              D5421 * sin(Xomi + X2li - G54) + D5433 * sin(-Xomi + X2li - G54);
        Xldot = XNi + Xfact;
        Xnddt = D2201 * cos(X2omi + Xli - G22) + D2211 * cos(Xli - G22) +
              D3210 * cos(Xomi + Xli - G32) + D3222 * cos(-Xomi + Xli - G32) +
              D5220 * cos(Xomi + Xli - G52) + D5232 * cos(-Xomi + Xli - G52) +
              2.0 * (D4410 * cos(X2omi + X2li - G44) +
              D4422 * cos(X2li - G44) + D5421 * cos(Xomi + X2li - G54) +
              D5433 * cos(-Xomi + X2li - G54));
        Xnddt = Xnddt * Xldot;
      }

      /* ----------------------- INTEGRATOR ------------------- */
      if (IRetn == 381)
      {
        Xli   = Xli + Xldot * Delt + Xndt * Step2;
        XNi   = XNi + Xndt * Delt + Xnddt * Step2;
        Atime = Atime + Delt;
      }
    }
  
    Nm = XNi + Xndt * Ft + Xnddt * Ft * Ft * 0.5;
    Xl = Xli + Xldot * Ft + Xndt * Ft * Ft * 0.5;
    if (IRez != 1)
    {
      Mm   = Xl - 2.0 * Omegam + 2.0 * Theta;
      Dndt = Nm - No;
    }
    else
    {
      Mm   = Xl - Omegam - Argpm+ Theta;
      Dndt = Nm - No;
    }

    Nm = No + Dndt;
  }

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%84s\n", 
                        " ------------------After DSPACE :--------------- ");
      fprintf(SGP4File, "    Inputs : \n");
      fprintf(SGP4File, 
              "%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
              "IRez", IRez, "D2201", D2201, "D2211", D2211,
              "D3210", D3210, "D3222", D3222, "D4410", D4410);
      fprintf(SGP4File, 
              "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
              "D4422", D4422, "D5220", D5220, "D5232", D5232,
              "D5421", D5421, "D5433", D5433, "Dedt", Dedt);
      fprintf(SGP4File, 
              "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
              "Del1", Del1, "Del2", Del2, "Del3", Del3,
              "Didt", Didt, "Dmdt", Dmdt, "Dnodt", Dnodt);
      fprintf(SGP4File, 
              "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
              "Domdt", Domdt, "Argpo", Argpo, "ArgpDot", ArgpDot,
              "T", T, "TC", TC, "GSTo", GSTo);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f\n",
                        "Xfact", Xfact, "Xlamo", Xlamo, "No", No);
      fprintf(SGP4File, "    In / Out : \n");
      fprintf(SGP4File, 
              "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
              "Atime", Atime, "EM", EM, "Argpm", Argpm,
              "Inclm", Inclm, "Xli", Xli, "Mm", Mm);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f\n", "Xni", XNi, "Omegam", Omegam);
      fprintf(SGP4File, "    Outputs : \n");
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f\n", "Dndt", Dndt, "Nm", Nm);
    }
}

/*-----------------------------------------------------------------------------
*
*                           PROCEDURE INITL
*
*  This PROCEDURE initializes the SPG4 propagator.
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    Ecco        - Eccentricity                           0.0 - 1.0
*    Epoch       - EPoch time in days from Jan 0, 1950. 0 hr
*    Inclo       - Inclination of satellite
*    No          - Mean motion of satellite
*    Satn        - Satellite number
*
*  Outputs       :
*    Ainv        - 1.0 / a
*    Ao          - Semi major axis
*    Con41       -
*    Con42       - 1.0 - 5.0 COS(i)
*    Cosio       - Cosine of inclination
*    Cosio2      - Cosio squared
*    Einv        - 1.0 / e
*    Eccsq       - Eccentricity squared
*    Method      -
*    Omeosq      - 1.0 - Ecco * Ecco
*    Posq        - Semi-parameter squared
*    rp          - Radius of perigee
*    RTeosq      - Square root of (1.0 - Ecco*Ecco)
*    Sinio       - Sine of inclination
*    GSTo        - GST at time of observation               rad
*    No          - Mean motion of satellite
*
*  Locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  Coupling      :
*    None.
*
*  References    :
*    NORAD Spacetrack Report #3
*
  ----------------------------------------------------------------------------*/
void InitL
    (
      LINT Satn, Real Ecco, Real Epoch, Real Inclo, Real& No, SINT& Method,
      Real& AINV, Real& AO, Real& CON41, Real& CON42, Real& COSIO, Real& COSIO2,
      Real& EINV, Real& EccSQ, Real& OMEOSQ, Real& POSQ, Real& rp, Real& RTEOSQ,
      Real& SINIO , Real& GSTo
    )
{
  /* --------------------- Local Variables ------------------------ */
  const Real TwoPi = 2.0 * PI;

  Real AK, D1, DEL, ADEL, PO, X2O3, J2, XKE;

  /* -------------------- WGS-72 EARTH CONSTANTS ----------------- */
  XKE  = 7.43669161331734132E-2;
  J2   = 1.082616E-3;
  X2O3 = 2.0 / 3.0;

  /* ------------- CALCULATE AUXILLARY EPOCH QUANTITIES ---------- */
  EccSQ  = Ecco * Ecco;
  OMEOSQ = 1.0 - EccSQ;
  RTEOSQ = sqrt(OMEOSQ);
  COSIO  = cos(Inclo);
  COSIO2 = COSIO * COSIO;

  /* ------------------ UN-KOZAI THE MEAN MOTION ----------------- */
  /* make sure .elm file has correct line endings!                 */
  AK    = Power(XKE / No, X2O3);
  D1    = 0.75 * J2 * (3.0 * COSIO2 - 1.0) / (RTEOSQ * OMEOSQ);
  DEL   = D1 / (AK * AK);
  ADEL  = AK * (1.0 - DEL * DEL - DEL *
          (1.0 / 3.0 + 134.0 * DEL * DEL / 81.0));
  DEL   = D1/(ADEL * ADEL);
  No    = No / (1.0 + DEL);

  AO    = Power(XKE / No, X2O3);
  SINIO = sin(Inclo);
  PO    = AO * OMEOSQ;
  CON42 = 1.0 - 5.0 * COSIO2;
  CON41 = -CON42-COSIO2-COSIO2;
  AINV  = 1.0 / AO;
  EINV  = 1.0 / Ecco;
  POSQ  = PO * PO;
  rp    = AO * (1.0 - Ecco);
  Method = 0;
  if (rp < 1.0)
    printf(" *** SATN%d EPOCH ELTS SUB-ORBITAL ***\n", Satn);

  GSTo = GSTime(Epoch + 2433281.5);

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%8fs\n", 
                        " ------------------After INITL  :---------------");
      fprintf(SGP4File, "    Inputs : \n");
      fprintf(SGP4File, "%7s%13d%7s%13s%7s%13.7f%7s%13.7f%7s%13.7f\n",
        "Satn", Satn, "Yr", " ", "Ecco", Ecco, "Epoch", Epoch, "Inclo", Inclo);
      fprintf(SGP4File, "    In/Out : \n");
      fprintf(SGP4File, "%7s%13.7f\n", "No", No);
      fprintf(SGP4File, "    Outputs : \n");
      fprintf(SGP4File, 
              "%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
               "Method", Method, "Ainv", AINV, "Ao", AO, "Con41", CON41, 
               "Con42", CON42, "Cosio", COSIO);
      fprintf(SGP4File, "%7s%13.7f\n", "Cosio2", COSIO2);
      fprintf(SGP4File, "%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13d%7s%13.7f\n",
                         "Einv", EINV, "Eccsq", EccSQ, "Omeosq", OMEOSQ, 
                         "posq", POSQ, "rp", rp, "Rteosq", RTEOSQ);
      fprintf(SGP4File, "%7s%13d%7s%13.7f\n", "Sinio", SINIO, "GSTo", GSTo);
    }
}

/*-----------------------------------------------------------------------------
*
*                             PROCEDURE SGP4
*
*  This PROCEDURE is the SGP4 SPADOC compatible prediction model from
*    Spacetrack report #3.  This is an updated version of SGP4 and SDP4, which
*    were originally published separately in Spacetrack Report #3. The version
*    follows the NASA release on the INternet. There are some fixes that are
*    added to the Em tolerances.
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    BSTAR       -
*    Ecco        - Initial eccentricity
*    Argpo       - Initial argument of perigee
*    Inclo       - Initial inclination
*    NEValues.T  - Time in minutes from the epoch of the satellite (set outside)
*    Mo          - Initial Mean Anomaly
*    No          - Initial Mean Motion
*    Omegao      - Initial Long of Ascending node
*    NeValues    - Near Earth common values for subsequent calls
*    DsValues    - Deep Space common values for calls
*
*  Outputs       :
*    X           - X component of position                          ER
*    Y           - Y component of position
*    Z           - Z component of position
*    XDoT        - X component of velocity                          ER/min
*    YDoT        - Y component of velocity
*    ZDoT        - Z component of velocity
*
*  Locals        :
*    AM          -
*    Axnl, Aynl        -
*    Betal       -
*    COSIM       -
*    COSOMM      -
*    Cnod        -
*    Cos2u       -
*    Coseo1      -
*    Cosi        -
*    Cosip       -
*    Cosisq      -
*    Cossu       -
*    Cosu        -
*    Delm        -
*    Delomg      -
*    Dndt        -
*    EM          -
*    EMSQ        -
*    Ecose       -
*    El2         -
*    Eo1         -
*    Ep          -
*    Esine       -
*    Argpm       -
*    Argpp       -
*    ArgpDF      -
*    Pl          -
*    R           -
*    RTEMSQ      -
*    Rdot        -
*    Rdotl       -
*    Rl          -
*    Rvdot       -
*    Rvdotl      -
*    SINIM       -
*    SINOMM      -
*    Sin2u       -
*    Sineo1      -
*    Sini        -
*    Sinip       -
*    Sinu        -
*    Snod        -
*    Su          -
*    T2, T3, T4          -
*    Tc          -
*    Tem5        -
*    Temp, Temp1, Temp2, Tempa, Tempe, Templ
*    U           -
*    Ux, Uy, Uz          -
*    Vx, Vy, Vz          -
*    Inclm       - Inclination
*    Mm          - Mean Anomaly
*    Nm          - Mean Motion
*    Omegam      - Longi of Ascending node
*    Xinc        -
*    Xincp       -
*    Xl          -
*    Xlm         -
*    Mp          -
*    Xmdf        -
*    Xmx         -
*    Xmy         -
*    OmegaDF     -
*    Xnode       -
*    Omegap      -
*    Np          -
*
*  Coupling      :
*    DPPER
*    DPSPACE
*
*  References    :
*    NORAD Spacetrack Report #3
*
  ----------------------------------------------------------------------------*/
void SGP4(ElSetRec SatRec, Vector& R, Vector& V, SINT& Error)
{
  Real AM   , Axnl  , Aynl , Betal ,  Cosim , Cnod  ,
       Cos2u, Coseo1, Cosi , Cosip ,  Cosisq, Cossu , Cosu,
       Delm , Delomg, EM   , EMSQ  ,  Ecose , El2   , EO1 ,
       Ep   , Esine , Argpm, Argpp ,  ArgpDF, Pl,     MRt ,
       MVt  , Rdotl, Rl    , Rvdot, Rvdotl,  Sinim ,
       Sin2u, Sineo1, Sini , Sinip ,  Sinsu , Sinu  ,
       Snod , Su    , T2   , T3    ,  T4    , Tem5  , Temp,
       Temp1, Temp2 , Tempa, Tempe ,  Templ , U     , Ux  ,
       Uy   , Uz    , Vx   , Vy    ,  Vz    , Inclm , Mm  ,
       Nm   , Omegam, Xinc , Xincp ,  Xl    , Xlm   , Mp  ,
       Xmdf , Xmx   , Xmy  , OmegaDF, Xnode , Omegap, Tc  , Dndt,
       TwoPi, X2O3,   J2,    J3,      XKE,    J3OJ2;

  /* -------------------- WGS-72 EARTH CONSTANTS ----------------- */
  /* ------------------ SET MATHEMATICAL CONSTANTS --------------- */
  TwoPi = 2.0 * PI;
  X2O3  = 2.0 / 3.0;
  XKE   = 7.43669161331734132E-2;
  J2    = 1.082616E-3;
  J3    = -2.53881E-6;
  J3OJ2 = J3 / J2;

  /* --------------------- CLEAR SGP4 ERROR FLAG ----------------- */
  Error = 0;

  /* ------- UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG ----- */
  Xmdf    = SatRec.Mo + SatRec.NEValues.MDot * SatRec.NEValues.T;
  ArgpDF  = SatRec.Argpo + SatRec.NEValues.ArgpDot * SatRec.NEValues.T;
  OmegaDF = SatRec.Omegao + SatRec.NEValues.OmegaDot * SatRec.NEValues.T;
  Argpm   = ArgpDF;
  Mm      = Xmdf;
  T2      = SatRec.NEValues.T * SatRec.NEValues.T;
  Omegam  = OmegaDF + SatRec.NEValues.OmegaCF * T2;
  Tempa   = 1.0 - SatRec.NEValues.Cc1 * SatRec.NEValues.T;
  Tempe   = SatRec.BStar * SatRec.NEValues.Cc4 * SatRec.NEValues.T;
  Templ   = SatRec.NEValues.T2cof * T2;
  
  if (SatRec.NEValues.Isimp != 1)
  {
    Delomg = SatRec.NEValues.Omgcof * SatRec.NEValues.T;
    Delm   = SatRec.NEValues.Xmcof * 
             (Power((1.0 + SatRec.NEValues.Eta * cos(Xmdf)), 3) - 
             SatRec.NEValues.Delmo);
    Temp   = Delomg + Delm;
    Mm     = Xmdf + Temp;
    Argpm  = ArgpDF - Temp;
    T3     = T2 * SatRec.NEValues.T;
    T4     = T3 * SatRec.NEValues.T;
    Tempa  = Tempa - SatRec.NEValues.D2 * T2 - SatRec.NEValues.D3 * T3 - 
                     SatRec.NEValues.D4 * T4;
    Tempe  = Tempe + SatRec.BStar * SatRec.NEValues.Cc5 * (sin(Mm) - 
                     SatRec.NEValues.Sinmao);
    Templ  = Templ + SatRec.NEValues.T3cof * T3 + T4 * (SatRec.NEValues.T4cof + 
                     SatRec.NEValues.T * SatRec.NEValues.T5cof);
  }

  Nm    = SatRec.No;
  EM    = SatRec.Ecco;
  Inclm = SatRec.Inclo;
  if (SatRec.NEValues.Method == 2)
  {
    Tc = SatRec.NEValues.T;
    DSpace
        (
          SatRec.DSValues.IRez,
          SatRec.DSValues.D2201, SatRec.DSValues.D2211, SatRec.DSValues.D3210,
          SatRec.DSValues.D3222, SatRec.DSValues.D4410, SatRec.DSValues.D4422,
          SatRec.DSValues.D5220, SatRec.DSValues.D5232, SatRec.DSValues.D5421,
          SatRec.DSValues.D5433, SatRec.DSValues.Dedt,  SatRec.DSValues.Del1,
          SatRec.DSValues.Del2,  SatRec.DSValues.Del3,  SatRec.DSValues.Didt,
          SatRec.DSValues.Dmdt,  SatRec.DSValues.Dnodt, SatRec.DSValues.Domdt,
          SatRec.Argpo,
          SatRec.NEValues.ArgpDot, SatRec.NEValues.T,
          Tc,
          SatRec.DSValues.GSTo, SatRec.DSValues.Xfact, SatRec.DSValues.Xlamo,
          SatRec.No,
          SatRec.DSValues.Atime,
          EM,
          Argpm,
          Inclm,
          SatRec.DSValues.Xli,
          Mm,
          SatRec.DSValues.Xni,
          Omegam,
          Dndt,
          Nm
        );
  }
 
  if (Nm <= 0.0)
  {
    printf("Error Nm %f\n", Nm);
    Error = 2;
  }
  AM = Power((XKE / Nm),X2O3) * Tempa * Tempa;
  Nm = XKE / Power(AM, 1.5);
  EM = EM - Tempe;
  if ((EM >= 1.0) || (EM < -0.001))
  {
    printf("Error Em %f\n", EM);
    Error = 1;
  }
  if (EM < 0.0)
    EM  = 1.0E-6;
  Mm     = Mm + SatRec.No * Templ;
  Xlm    = Mm + Argpm + Omegam;
  EMSQ   = EM * EM;
  Temp   = 1.0 - EMSQ;
  // RTEMSQ = sqrt(Temp);
  Omegam = Mod(Omegam, TwoPi);
  Argpm  = Mod(Argpm, TwoPi);
  Xlm    = Mod(Xlm, TwoPi);
  Mm     = Mod(Xlm - Argpm - Omegam, TwoPi);
  if (Mm < 0.0)
    Mm = Mm + TwoPi;

  /* ----------------- COMPUTE EXTRA MEAN QUANTITIES ------------- */
  Sinim = sin(Inclm);
  Cosim = cos(Inclm);
  // Sinim = sin(Argpm);
  // Cosim = cos(Argpm);

  /* -------------------- ADD LUNAR-SOLAR PERIODICS -------------- */
  // Np    = Nm;
  Ep     = EM;
  Xincp  = Inclm;
  Argpp  = Argpm;
  Omegap = Omegam;
  Mp     = Sinim;
  Cosip  = Cosim;
  if (SatRec.NEValues.Method == 2)
  {
    DPPer
        (
          SatRec.DSValues.E3,   SatRec.DSValues.Ee2,  SatRec.DSValues.Peo, 
          SatRec.DSValues.Pgho, SatRec.DSValues.Pho,  SatRec.DSValues.Pinco,
          SatRec.DSValues.Plo,  SatRec.DSValues.Se2,  SatRec.DSValues.Se3, 
          SatRec.DSValues.Sgh2, SatRec.DSValues.Sgh3, SatRec.DSValues.Sgh4,
          SatRec.DSValues.Sh2,  SatRec.DSValues.Sh3,  SatRec.DSValues.Si2, 
          SatRec.DSValues.Si3,  SatRec.DSValues.Sl2,  SatRec.DSValues.Sl3,
          SatRec.DSValues.Sl4,  SatRec.NEValues.T,    SatRec.DSValues.Xgh2, 
          SatRec.DSValues.Xgh3, SatRec.DSValues.Xgh4, SatRec.DSValues.Xh2,
          SatRec.DSValues.Xh3,  SatRec.DSValues.Xi2,  SatRec.DSValues.Xi3, 
          SatRec.DSValues.Xl2,  SatRec.DSValues.Xl3,  SatRec.DSValues.Xl4,
          SatRec.DSValues.Zmol, SatRec.DSValues.Zmos, 
          0, Ep, Xincp, Omegap, Argpp, Mp
        );
    if (Xincp < 0.0)
    {
      Xincp  = -Xincp;
      Omegap = Omegap + PI;
      Argpp  = Argpp - PI;
    }
    if ((Ep < 0.0 ) || ( Ep > 1.0))
    {
      printf("Error Ep %f\n", Ep);
      Error = 3;
    }
  }

  /* -------------------- LONG PERIOD PERIODICS ------------------ */
  if (SatRec.NEValues.Method == 2)
  {
    Sinip =  sin(Xincp);
    Cosip =  cos(Xincp);
    SatRec.NEValues.Aycof = -0.5*J3OJ2*Sinip;
    SatRec.NEValues.Xlcof = -0.25 * J3OJ2 * Sinip * (3.0 + 5.0 * Cosip) / 
                            (1.0+Cosip);
  }
  Axnl = Ep * cos(Argpp);
  Temp = 1.0 / (AM * (1.0 - Ep * Ep));
  Aynl = Ep* sin(Argpp) + Temp * SatRec.NEValues.Aycof;
  Xl   = Mp + Argpp + Omegap + Temp * SatRec.NEValues.Xlcof * Axnl;

  /* --------------------- SOLVE KEPLER'S EQUATION --------------- */
  U    = Mod(Xl - Omegap, TwoPi);
  EO1  = U;
  Tem5 = 9999.9;
  while (fabs(Tem5) >= 1.0E-12)
  {
    Sineo1 = sin(EO1);
    Coseo1 = cos(EO1);
    Tem5   = 1.0 - Coseo1 * Axnl - Sineo1 * Aynl;
    Tem5   = (U - Aynl * Coseo1 + Axnl * Sineo1 - EO1) / Tem5;
    EO1    = EO1 + Tem5;
  }

  /* ------------- SHORT PERIOD PRELIMINARY QUANTITIES ----------- */
  Ecose = Axnl*Coseo1+Aynl*Sineo1;
  Esine = Axnl*Sineo1-Aynl*Coseo1;
  El2   = Axnl*Axnl+Aynl*Aynl;
  Pl    = AM*(1.0-El2);
  if (Pl < 0.0)
  {
    printf("Error Pl %f\n", Pl);
    Error = 4;
  }
  else
  {
    Rl     = AM * (1.0 - Ecose);
    Rdotl  = sqrt(AM) * Esine/Rl;
    Rvdotl = sqrt(Pl) / Rl;
    Betal  = sqrt(1.0 - El2);
    Temp   = Esine / (1.0 + Betal);
    Sinu   = AM / Rl * (Sineo1 - Aynl - Axnl * Temp);
    Cosu   = AM / Rl * (Coseo1 - Axnl + Aynl * Temp);
    Su     = Atan2(Sinu, Cosu);
    Sin2u  = (Cosu + Cosu) * Sinu;
    Cos2u  = 1.0 - 2.0 * Sinu * Sinu;
    Temp   = 1.0 / Pl;
    Temp1  = 0.5 * J2 * Temp;
    Temp2  = Temp1 * Temp;

    /* -------------- UPDATE FOR SHORT PERIOD PERIODICS ------------ */
    if (SatRec.NEValues.Method == 2)
    {
      Cosisq                 = Cosip * Cosip;
      SatRec.NEValues.CON41  = 3.0*Cosisq - 1.0;
      SatRec.NEValues.X1mth2 = 1.0 - Cosisq;
      SatRec.NEValues.X7thm1 = 7.0*Cosisq - 1.0;
    }
    MRt   = Rl * (1.0 - 1.5 * Temp2 * Betal * SatRec.NEValues.CON41) + 
            0.5 * Temp1 * SatRec.NEValues.X1mth2 * Cos2u;
    Su    = Su - 0.25 * Temp2 * SatRec.NEValues.X7thm1 * Sin2u;
    Xnode = Omegap + 1.5 * Temp2 * Cosip * Sin2u;
    Xinc  = Xincp + 1.5 * Temp2 * Cosip * Sinip * Cos2u;
    MVt   = Rdotl - Nm * Temp1 * SatRec.NEValues.X1mth2 * Sin2u / XKE;
    Rvdot = Rvdotl + Nm * Temp1 * (SatRec.NEValues.X1mth2 * Cos2u + 
            1.5 * SatRec.NEValues.CON41) / XKE;

    /* --------------------- ORIENTATION VECTORS ------------------- */
    Sinsu =  sin(Su);
    Cossu =  cos(Su);
    Snod  =  sin(Xnode);
    Cnod  =  cos(Xnode);
    Sini  =  sin(Xinc);
    Cosi  =  cos(Xinc);
    Xmx   = -Snod * Cosi;
    Xmy   =  Cnod * Cosi;
    Ux    =  Xmx * Sinsu + Cnod * Cossu;
    Uy    =  Xmy * Sinsu + Snod * Cossu;
    Uz    =  Sini * Sinsu;
    Vx    =  Xmx * Cossu - Cnod * Sinsu;
    Vy    =  Xmy * Cossu - Snod * Sinsu;
    Vz    =  Sini * Cossu;

    /* ------------------- POSITION AND VELOCITY ------------------- */
    R.Set(MRt * Ux, 1);
    R.Set(MRt * Uy, 2);
    R.Set(MRt * Uz, 3);
    V.Set(MVt * Ux * Rvdot * Vx, 1);
    V.Set(MVt * Uy * Rvdot * Vy, 2);
    V.Set(MVt * Uz * Rvdot * Vz, 3);
  }

  if (Error > 0)
    printf("*** Error: t:= %f *** code = %3d\n", SatRec.NEValues.T, Error);

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%84s\n", 
                        "------------------After SGP4   :---------------");
      fprintf(SGP4File, "    Inputs : \n");
      fprintf(SGP4File, "%7s%13d%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Isimp", SatRec.NEValues.Isimp, "Method", SatRec.NEValues.Method, 
          "Aycof", SatRec.NEValues.Aycof, "BSTAR",  SatRec.BStar,
          "CON41", SatRec.NEValues.CON41, "Cc1",    SatRec.NEValues.Cc1);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Cc4", SatRec.NEValues.Cc4, "Cc5",   SatRec.NEValues.Cc5,
          "D2",  SatRec.NEValues.D2,  "D3",    SatRec.NEValues.D3,
          "D4",  SatRec.NEValues.D4,  "Delmo", SatRec.NEValues.Delmo);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Ecco",   SatRec.Ecco,            "Eta",     SatRec.NEValues.Eta,
          "Argpo",  SatRec.Argpo,           "ArgpDot", SatRec.NEValues.ArgpDot,
          "Omgcof", SatRec.NEValues.Omgcof, "Sinmao",  SatRec.NEValues.Sinmao);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "T",     SatRec.NEValues.T,     "T2cof",  SatRec.NEValues.T2cof,
          "T3cof", SatRec.NEValues.T3cof, "T4cof",  SatRec.NEValues.T4cof,
          "T5cof", SatRec.NEValues.T5cof, "X1mth2", SatRec.NEValues.X1mth2);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "X7thm1", SatRec.NEValues.X7thm1, "Inclo",  SatRec.Inclo,
          "Mo",     SatRec.Mo,              "MDot",   SatRec.NEValues.MDot,
          "XNO",    SatRec.No,              "Omegao", SatRec.Omegao);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "OmegaDot", SatRec.NEValues.OmegaDot, "Xlcof", SatRec.NEValues.Xlcof,
          "Xmcof",    SatRec.NEValues.Xmcof,    
          "OmegaCF",  SatRec.NEValues.OmegaCF);
      fprintf(SGP4File, "    Outputs : \n");
      fprintf(SGP4File, 
          "%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Error", Error, "X", R.Get(1), "Y", R.Get(2), "Z", R.Get(3),
          "XDOT", V.Get(1), "YDOT", V.Get(2), "ZDOT", V.Get(3));
      fprintf(SGP4File, "    Extra Inputs for DS : \n");
      fprintf(SGP4File, "%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "IRez",  SatRec.DSValues.IRez,  "D2201", SatRec.DSValues.D2201,
          "D2211", SatRec.DSValues.D2211, "D3210", SatRec.DSValues.D3210,
          "D3222", SatRec.DSValues.D3222);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "D4410", SatRec.DSValues.D4410, "D4422", SatRec.DSValues.D4422,
          "D5220", SatRec.DSValues.D5220, "D5232", SatRec.DSValues.D5232,
          "D5421", SatRec.DSValues.D5421, "D5433", SatRec.DSValues.D5433);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Dedt", SatRec.DSValues.Dedt, "Del1", SatRec.DSValues.Del1,
          "Del2", SatRec.DSValues.Del2, "Del3", SatRec.DSValues.Del3,
          "Didt", SatRec.DSValues.Didt, "Dmdt", SatRec.DSValues.Dmdt);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Dnodt", SatRec.DSValues.Dnodt, "Domdt", SatRec.DSValues.Domdt,
          "E3",    SatRec.DSValues.E3,    "Ee2",   SatRec.DSValues.Ee2,
          "Peo",   SatRec.DSValues.Peo,   "Pgho",  SatRec.DSValues.Pgho);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Pho", SatRec.DSValues.Pho, "Pinco", SatRec.DSValues.Pinco,
          "Plo", SatRec.DSValues.Plo, "Se2",   SatRec.DSValues.Se2,
          "Se3", SatRec.DSValues.Se3, "Sgh2",  SatRec.DSValues.Sgh2);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Sgh3", SatRec.DSValues.Sgh3, "Sgh4", SatRec.DSValues.Sgh4,
          "Sh2",  SatRec.DSValues.Sh2,  "Sh3",  SatRec.DSValues.Sh3,
          "Si2",  SatRec.DSValues.Si2,  "Si3",  SatRec.DSValues.Si3);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Sl2",   SatRec.DSValues.Sl2,   "Sl3",  SatRec.DSValues.Sl3,
          "Sl4",   SatRec.DSValues.Sl4,   "GSTo", SatRec.DSValues.GSTo,
          "Xfact", SatRec.DSValues.Xfact, "Xgh2", SatRec.DSValues.Xgh2);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Xgh3", SatRec.DSValues.Xgh3, "Xgh4", SatRec.DSValues.Xgh4,
          "Xh2",  SatRec.DSValues.Xh2,  "Xh3",  SatRec.DSValues.Xh3,
          "Xi2",  SatRec.DSValues.Xi2, "Xi3",   SatRec.DSValues.Xi3);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Xl2",  SatRec.DSValues.Xl2,  "Xl3",   SatRec.DSValues.Xl3,
          "Xl4",  SatRec.DSValues.Xl4,  "Xlamo", SatRec.DSValues.Xlamo,
          "Zmol", SatRec.DSValues.Zmol, "Zmos",  SatRec.DSValues.Zmos);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Atime", SatRec.DSValues.Atime, "Xli", SatRec.DSValues.Xli,
          "Xni",   SatRec.DSValues.Xni);
    }
}

/*-----------------------------------------------------------------------------
*
*                             PROCEDURE SGP4INIT
*
*  This PROCEDURE initializes variables for SGP4. The fix to conform to
*    Ver 3.01 consisted only of adding a VAR so the 'o elements would be
*    passed back for deep space initialization cases.
*
*  Author        : David Vallado                  303-344-6037    1 Mar 2001
*
*  Inputs        :
*    Satn        - Satellite Number
*    Year        - Year of observation                    1950-2049
*    BSTAR       - SGP4 type drag Coefficient              kg/m2ER
*    Ecco        - Eccentricity
*    Epoch       - EPoch time in days from Jan 0, 1950. 0 hr
*    Argpo       - Argument of perigee (output IF ds)
*    Inclo       - Inclination
*    Mo          - Mean Anomaly (output IF ds)
*    No          - Mean motion
*    Omegao      - Longitude of ascending node
*
*  Outputs       :
*    No          - Mean motion
*    Omegao      - Longitude of ascending node
*    Init        - Flag for first pass 1-yes, 0-not first pass
*    NeValues    - Near Earth common values for subsequent calls
*    DsValues    - Deep Space common values for calls
*
*  Locals        :
*    CNODM       -
*    SNODM       -
*    COSIM       -
*    SINIM       -
*    COSOMM      -
*    SINOMM      -
*    Cc1sq       -
*    Cc2         -
*    Cc3         -
*    Coef        -
*    Coef1       -
*    Cosio4      -
*    DAY         -
*    Dndt        -
*    EM          -
*    EMSQ        -
*    Eeta        -
*    Etasq       -
*    GAM         -
*    Argpm       -
*    Ndem        -
*    Inclm       - Inclination
*    Mm          - Mean Anomaly
*    Nm          - Mean Motion
*    Perige      -
*    Pinvsq      -
*    Psisq       -
*    Qzms24      -
*    RTEMSQ      -
*    S1, S2, S3, S4, S5, S6, S7          -
*    SFour       -
*    SS1, SS2, SS3, SS4, SS5, SS6, SS7         -
*    SZ1, SZ2, SZ3
*    SZ11, SZ12, SZ13, SZ21, SZ22, SZ23, SZ31, SZ32, SZ33        -
*    Tc          -
*    Temp        -
*    Temp1, Temp2, Temp3       -
*    Tsi         -
*    XPIDOT      -
*    Xhdot1      -
*    Z1, Z2, Z3          -
*    Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33         -
*
*  Coupling      :
*    INITL       -
*    DSCom       -
*    DPPER       -
*    DSINIT      -
*
*  References    :
*    NORAD SpaceTrack Report No.3 "Models For Propagation of Norad Element Sets"
*              1 NOV 1980:  Felix R. Hoots, Ronald L. Roehrich
*
  ----------------------------------------------------------------------------*/
void SGP4Init
    (
      LINT Satn, SINT Year, Real Epoch, Real& BStar, Real& Ecco, Real& Argpo,
      Real& Inclo, Real& Mo, Real& No, Real& Omegao,
      SINT& Init, NearEarthType& NEValues, DeepSpaceType& DSValues
    )
{
  /* --------------------- Local Variables ------------------------ */
  Real AO, AINV, EINV, CON42, COSIO, SINIO, COSIO2, EccSQ, 
       OMEOSQ, POSQ,   rp,     RTEOSQ,
       CNODM , SNODM , COSIM , SINIM , COSOMM, SINOMM, Cc1sq ,
       Cc2   , Cc3   , Coef  , Coef1 , COSIO4, DAY   , Dndt  ,
       EM    , EMSQ  , Eeta  , EtaSQ , GAM   , Argpm , Omegam,
       Inclm , Mm    , Nm    , Perige, PinvSQ, PsiSQ , Qzms24,
       RTEMSQ, S1    , S2    , S3    , S4    , S5    , S6    ,
       S7    , SFour , SS1   , SS2   , SS3   , SS4   , SS5   ,
       SS6   , SS7   , SZ1   , SZ2   , SZ3   , SZ11  , SZ12  ,
       SZ13  , SZ21  , SZ22  , SZ23  , SZ31  , SZ32  , SZ33  ,
       Tc    , Temp  , Temp1 , Temp2 , Temp3 , TSI   , XPIDOT,
       Xhdot1, Z1    , Z2    , Z3    , Z11   , Z12   , Z13   ,
       Z21   , Z22   , Z23   , Z31   , Z32   , Z33,
       Qzms2t, SS, RadiusEarthKm, twopi, J2,J3oJ2,J4,X2o3;

  /* ------------------------ INITIALIZATION --------------------- */
  /* ----------- set all Near Earth variables to zero ------------ */
  NEValues.Isimp   = 0;   NEValues.Method = 0;   NEValues.Aycof    = 0.0; 
  NEValues.CON41   = 0.0; NEValues.Cc1    = 0.0; NEValues.Cc4      = 0.0; 
  NEValues.Cc5     = 0.0; NEValues.D2     = 0.0; NEValues.D3       = 0.0; 
  NEValues.D4      = 0.0; NEValues.Delmo  = 0.0; NEValues.Eta      = 0.0; 
  NEValues.ArgpDot = 0.0; NEValues.Omgcof = 0.0; NEValues.Sinmao   = 0.0;
  NEValues.T       = 0.0; NEValues.T2cof  = 0.0; NEValues.T3cof    = 0.0; 
  NEValues.T4cof   = 0.0; NEValues.T5cof  = 0.0; NEValues.X1mth2   = 0.0; 
  NEValues.X7thm1  = 0.0; NEValues.MDot   = 0.0; NEValues.OmegaDot = 0.0; 
  NEValues.Xlcof   = 0.0; NEValues.Xmcof  = 0.0; NEValues.OmegaCF  = 0.0;

  /* ----------- set all Deep Space variables to zero ------------ */
  DSValues.IRez  = 0;   DSValues.D2201 = 0.0; DSValues.D2211 = 0.0; 
  DSValues.D3210 = 0.0; DSValues.D3222 = 0.0; DSValues.D4410 = 0.0; 
  DSValues.D4422 = 0.0; DSValues.D5220 = 0.0; DSValues.D5232 = 0.0; 
  DSValues.D5421 = 0.0; DSValues.D5433 = 0.0; DSValues.Dedt  = 0.0; 
  DSValues.Del1  = 0.0; DSValues.Del2  = 0.0; DSValues.Del3  = 0.0;
  DSValues.Didt  = 0.0; DSValues.Dmdt  = 0.0; DSValues.Dnodt = 0.0; 
  DSValues.Domdt = 0.0; DSValues.E3    = 0.0; DSValues.Ee2   = 0.0; 
  DSValues.Peo   = 0.0; DSValues.Pgho  = 0.0; DSValues.Pho   = 0.0; 
  DSValues.Pinco = 0.0; DSValues.Plo   = 0.0; DSValues.Se2   = 0.0; 
  DSValues.Se3   = 0.0; DSValues.Sgh2  = 0.0; DSValues.Sgh3  = 0.0;
  DSValues.Sgh4  = 0.0; DSValues.Sh2   = 0.0; DSValues.Sh3   = 0.0; 
  DSValues.Si2   = 0.0; DSValues.Si3   = 0.0; DSValues.Sl2   = 0.0; 
  DSValues.Sl3   = 0.0; DSValues.Sl4   = 0.0; DSValues.GSTo  = 0.0; 
  DSValues.Xfact = 0.0; DSValues.Xgh2  = 0.0; DSValues.Xgh3  = 0.0; 
  DSValues.Xgh4  = 0.0; DSValues.Xh2   = 0.0; DSValues.Xh3   = 0.0;
  DSValues.Xi2   = 0.0; DSValues.Xi3   = 0.0; DSValues.Xl2   = 0.0; 
  DSValues.Xl3   = 0.0; DSValues.Xl4   = 0.0; DSValues.Xlamo = 0.0; 
  DSValues.Zmol  = 0.0; DSValues.Zmos  = 0.0; DSValues.Atime = 0.0; 
  DSValues.Xli   = 0.0; DSValues.Xni   = 0.0;

  RadiusEarthKm = 6378.135;
  SS     = 78.0 / RadiusEarthKm + 1.0;
  Qzms24 = Power(((120.0 - 78.0) / RadiusEarthKm), 4);
  J2     =  1.082616E-3;
  J3oJ2  = -2.53881E-6 / J2;
  X2o3   =  2.0 / 3.0;
  J4     = -1.65597E-6;

  InitL
      (
        Satn, Ecco, Epoch, Inclo, No, NEValues.Method, AINV, AO, 
        NEValues.CON41, CON42, COSIO, COSIO2, EINV, EccSQ, OMEOSQ, 
        POSQ, rp, RTEOSQ, SINIO, DSValues.GSTo
      );

  if ((OMEOSQ >= 0.0 ) || ( No >= 0.0))
  {
    NEValues.Isimp = 0;
    if (rp < (220.0 / RadiusEarthKm + 1.0))
      NEValues.Isimp = 1;
    SFour  = SS;
    Qzms24 = Qzms2t;
    Perige = (rp - 1.0) * RadiusEarthKm;
    
    /* - For perigees below 156 km, S and Qoms2t are altered - */
    if (Perige < 156.0)
    {
      SFour = Perige - 78.0;
      if (Perige < 98.0)
        SFour = 20.0;
      Qzms24 = Power(((120.0 - SFour) / RadiusEarthKm), 4.0);
      SFour  = SFour / RadiusEarthKm + 1.0;
    }
    PinvSQ = 1.0 / POSQ;
    
    TSI  = 1.0 / (AO - SFour);
    NEValues.Eta  = AO * Ecco * TSI;
    EtaSQ = NEValues.Eta * NEValues.Eta;
    Eeta  = Ecco * NEValues.Eta;
    PsiSQ = fabs(1.0 - EtaSQ);
    Coef  = Qzms24 * Power(TSI, 4.0);
    Coef1 = Coef / Power(PsiSQ, 3.5);
    Cc2   = Coef1 * No * (AO * (1.0 + 1.5 * EtaSQ + Eeta * 
                   (4.0 + EtaSQ)) + 0.375 * J2 * TSI / PsiSQ * NEValues.CON41 * 
                   (8.0 + 3.0 * EtaSQ * (8.0 + EtaSQ)));
    NEValues.Cc1   = BStar * Cc2;
    Cc3   = 0.0;
    if (Ecco > 1.0E-4)
      Cc3 = -2.0 * Coef * TSI * J3oJ2 * No * SINIO / Ecco;
    NEValues.X1mth2 = 1.0 - COSIO2;
    NEValues.Cc4    = 2.0* No * Coef1 * AO * OMEOSQ * 
                      (NEValues.Eta * (2.0 + 0.5 * EtaSQ) + Ecco * 
                      (0.5 + 2.0 * EtaSQ) - J2 * TSI / (AO * PsiSQ) * 
                      (-3.0 * NEValues.CON41 * (1.0 - 2.0 * Eeta + EtaSQ * 
                      (1.5 - 0.5 * Eeta)) + 0.75 * NEValues.X1mth2 * 
                      (2.0 * EtaSQ - Eeta * (1.0 + EtaSQ)) * cos(2.0 * Argpo)));
    NEValues.Cc5 = 2.0 * Coef1 * AO * OMEOSQ * (1.0 + 2.75 *
                   (EtaSQ + Eeta) + Eeta * EtaSQ);
    COSIO4 = COSIO2 * COSIO2;
    Temp1  = 1.5 * J2 * PinvSQ * No;
    Temp2  = 0.5 * Temp1 * J2 * PinvSQ;
    Temp3  = -0.46875 * J4 * PinvSQ * PinvSQ * No;
    NEValues.MDot     = No + 0.5 * Temp1 * RTEOSQ * NEValues.CON41 + 0.0625 * 
                       Temp2 * RTEOSQ * (13.0 - 78.0 * COSIO2 + 137.0 * COSIO4);
    NEValues.ArgpDot  = -0.5 * Temp1 * CON42 + 0.0625 * Temp2 * 
                        (7.0 - 114.0 * COSIO2 + 395.0 * COSIO4) + 
                        Temp3 * (3.0 - 36.0 * COSIO2 + 49.0 * COSIO4);
    Xhdot1            = -Temp1 * COSIO;
    NEValues.OmegaDot = Xhdot1 + (0.5 * Temp2 * (4.0 - 19.0 * COSIO2) +
                         2.0 * Temp3 * (3.0 - 7.0 * COSIO2)) * COSIO;
    XPIDOT            =  NEValues.ArgpDot+ NEValues.OmegaDot;
    NEValues.Omgcof   = BStar * Cc3 * cos(Argpo);
    NEValues.Xmcof    = 0.0;
    if (Ecco > 1.0E-4)
      NEValues.Xmcof = -X2o3 * Coef * BStar / Eeta;
    NEValues.OmegaCF = 3.5 * OMEOSQ * Xhdot1 * NEValues.Cc1;
    NEValues.T2cof   = 1.5 * NEValues.Cc1;
    NEValues.Xlcof   = -0.25 * J3oJ2 * SINIO * 
                       (3.0 + 5.0 * COSIO) / (1.0 + COSIO);
    NEValues.Aycof   = -0.5 * J3oJ2 * SINIO;
    NEValues.Delmo   = Power((1.0 + NEValues.Eta * cos(Mo)), 3);
    NEValues.Sinmao  = sin(Mo);
    NEValues.X7thm1  = 7.0 * COSIO2 - 1.0;

    Init = 0;

    /* --------------- Deep Space Initialization ------------- */
    if ((TwoPi / No) >= 225.0)
    {
      NEValues.Method = 2;
      NEValues.Isimp  = 1;
      Tc    =  0.0;
      Inclm = Inclo;

      DSCom
          (
            Epoch, Ecco,   Argpo, Tc,       Inclo,        Omegao, No,
            SNODM, CNODM,  SINIM, COSIM,    SINOMM,       COSOMM,
            DAY, DSValues.E3, DSValues.Ee2, EM,           EMSQ, GAM,
            DSValues.Peo,  DSValues.Pgho,   DSValues.Pho, DSValues.Pinco, 
            DSValues.Plo,  RTEMSQ,          DSValues.Se2, DSValues.Se3,
            DSValues.Sgh2, DSValues.Sgh3,   DSValues.Sgh4,
            DSValues.Sh2,  DSValues.Sh3,    DSValues.Si2, DSValues.Si3, 
            DSValues.Sl2,  DSValues.Sl3,    DSValues.Sl4, S1, S2, S3, S4, S5,
            S6,   S7,   SS1,  SS2,  SS3,  SS4,  SS5,  SS6,  SS7, SZ1, SZ2, SZ3,
            SZ11, SZ12, SZ13, SZ21, SZ22, SZ23, SZ31, SZ32, SZ33, 
            DSValues.Xgh2, DSValues.Xgh3,   DSValues.Xgh4, DSValues.Xh2, 
            DSValues.Xh3,  DSValues.Xi2,    DSValues.Xi3,  DSValues.Xl2, 
            DSValues.Xl3,  DSValues.Xl4,    Nm, Z1, Z2, Z3, Z11,
            Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, 
            DSValues.Zmol, DSValues.Zmos
          );
      DPPer
          (
            DSValues.E3, DSValues.Ee2, DSValues.Peo, DSValues.Pgho, 
            DSValues.Pho, DSValues.Pinco, DSValues.Plo, DSValues.Se2, 
            DSValues.Se3, DSValues.Sgh2, DSValues.Sgh3, DSValues.Sgh4,
            DSValues.Sh2, DSValues.Sh3, DSValues.Si2, DSValues.Si3, 
            DSValues.Sl2, DSValues.Sl3, DSValues.Sl4, NEValues.T, 
            DSValues.Xgh2, DSValues.Xgh3, DSValues.Xgh4, DSValues.Xh2,
            DSValues.Xh3, DSValues.Xi2, DSValues.Xi3, DSValues.Xl2, 
            DSValues.Xl3, DSValues.Xl4, DSValues.Zmol, DSValues.Zmos, 1,
            Ecco, Inclo, Omegao, Argpo, Mo
          );
      DSInit
          (
            COSIM, EMSQ, Argpo, S1, S2, S3, S4, S5, SINIM, SS1, SS2, SS3, SS4,
            SS5, SZ1, SZ3, SZ11, SZ13, SZ21, SZ23, SZ31, SZ33, NEValues.T, Tc,
            DSValues.GSTo, Mo, NEValues.MDot, No, Omegao, NEValues.OmegaDot,
            XPIDOT, Z1, Z3, Z11, Z13, Z21, Z23, Z31, Z33,
            Ecco, EccSQ, EM, Argpm, Inclm, Mm, Nm, Omegam,
            DSValues.IRez,  DSValues.Atime, 
            DSValues.D2201, DSValues.D2211, DSValues.D3210, DSValues.D3222 ,
            DSValues.D4410, DSValues.D4422, DSValues.D5220, DSValues.D5232, 
            DSValues.D5421, DSValues.D5433, DSValues.Dedt,  DSValues.Didt,
            DSValues.Dmdt,  Dndt,           DSValues.Dnodt, DSValues.Domdt ,
            DSValues.Del1,  DSValues.Del2,  DSValues.Del3,  DSValues.Xfact, 
            DSValues.Xlamo, DSValues.Xli,   DSValues.Xni
          );
    }
 
    /* ----------- Set variables IF not deep space ----------- */
    if (NEValues.Isimp != 1)
    {
      Cc1sq          = NEValues.Cc1 * NEValues.Cc1;
      NEValues.D2    = 4.0 * AO * TSI * Cc1sq;
      Temp           = NEValues.D2 * TSI * NEValues.Cc1 / 3.0;
      NEValues.D3    = (17.0 * AO + SFour) * Temp;
      NEValues.D4    = 0.5 * Temp * AO * TSI * (221.0 * AO + 31.0 * SFour) * 
                       NEValues.Cc1;
      NEValues.T3cof = NEValues.D2 + 2.0 * Cc1sq;
      NEValues.T4cof = 0.25 * (3.0 * NEValues.D3 + NEValues.Cc1 * 
                       (12.0 * NEValues.D2 + 10.0 * Cc1sq));
      NEValues.T5cof = 0.2 * (3.0 * NEValues.D4 + 
                       12.0 * NEValues.Cc1 * NEValues.D3 + 
                       6.0 * NEValues.D2 * NEValues.D2 + 
                       15.0 * Cc1sq * (2.0 * NEValues.D2 + Cc1sq));
    }
  }

  if (Help == 'Y')
    if (SGP4File != NULL)
    {
      fprintf(SGP4File, "%7s\n", 
                        " ------------------After SGP4Init :-------------");
      fprintf(SGP4File, "    Inputs  : \n");
      fprintf(SGP4File,
          "%7s%13d%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Satn", Satn, "Yr", Year, "BStar", BStar, "Ecco", Ecco, 
          "Epoch", Epoch, "Argpo", Argpo);
      fprintf(SGP4File, "f%7s%13.7f%7s%13.7f\n", "Inclo", Inclo, "Mo", Mo);
      fprintf(SGP4File, " In and Out variables \n");
      fprintf(SGP4File, "f%7s%13.7f\n", "No", No);
      fprintf(SGP4File, "    Outputs  :\n");
      fprintf(SGP4File, "%7s%13d%7s%13d%7s%13d%7s%13.7f\n",
          "INIT", Init, "Isimp", NEValues.Isimp, 
          "Method", NEValues.Method, "Aycof", NEValues.Aycof);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "CON41", NEValues.CON41, "Cc1", NEValues.Cc1, "Cc4", NEValues.Cc4);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Cc5", NEValues.Cc5, "D2",    NEValues.D2, "D3", NEValues.D3,
          "D4",  NEValues.D4,  "Delmo", NEValues.Delmo);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Eta", NEValues.Eta, "ArgpDot", NEValues.ArgpDot,
          "Omgcof", NEValues.Omgcof);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Sinmao", NEValues.Sinmao, "T2cof", NEValues.T2cof,
          "T3cof", NEValues.T3cof);
      fprintf(SGP4File,
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "T4cof",  NEValues.T4cof,  "T5cof",  NEValues.T5cof, 
          "GSTo",   DSValues.GSTo,   "X1mth2", NEValues.X1mth2, 
          "X7thm1", NEValues.X7thm1, "Xlcof",  NEValues.Xlcof);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Xmcof",   NEValues.Xmcof,   "MDot",     NEValues.MDot, 
          "OmegaCF", NEValues.OmegaCF, "OmegaDot", NEValues.OmegaDot);
      fprintf(SGP4File, "   IN and Outputs from Deep space satellites :\n");
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f\n",
          "T", NEValues.T, "Omegao", Omegao);
      fprintf(SGP4File, "%7s%13d%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Irez",  DSValues.IRez,  "Atime", DSValues.Atime, 
          "D2201", DSValues.D2201, "D2211", DSValues.D2211, 
          "D3210", DSValues.D3210);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "D3222", DSValues.D3222, "D4410", DSValues.D4410, 
          "D4422", DSValues.D4422, "D5220", DSValues.D5220, 
          "D5232", DSValues.D5232, "D5421", DSValues.D5421);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%f\n",
          "D5433", DSValues.D5433, "Dedt", DSValues.Dedt, 
          "Del1",  DSValues.Del1,  "Del2", DSValues.Del2, 
          "Del3",  DSValues.Del3);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%f\n",
          "Didt",  DSValues.Didt,  "Dmdt",  DSValues.Dmdt, 
          "Dnodt", DSValues.Dnodt, "Domdt", DSValues.Domdt, 
          "E3",    DSValues.E3);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Ee2", DSValues.Ee2, "Peo",   DSValues.Peo,   "Pgho", DSValues.Pgho, 
          "Pho", DSValues.Pho, "Pinco", DSValues.Pinco, "Plo",  DSValues.Plo);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Se2",  DSValues.Se2,  "Se3",  DSValues.Se3,  "Sgh2", DSValues.Sgh2,
          "Sgh3", DSValues.Sgh3, "Sgh4", DSValues.Sgh4, "Sh2",  DSValues.Sh2);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Sh3", DSValues.Sh3, "Si2", DSValues.Si2, "Si3", DSValues.Si3,
          "Sl2", DSValues.Sl2, "Sl3", DSValues.Sl3, "Sl4", DSValues.Sl4);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%f\n",
          "Xfact", DSValues.Xfact, "Xgh2", DSValues.Xgh2, "Xgh3", DSValues.Xgh3,
          "Xgh4",  DSValues.Xgh4,  "Xh2",  DSValues.Xh2);
      fprintf(SGP4File, "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%f\n",
          "Xh3",  DSValues.Xh3, "Xi2", DSValues.Xi2, "Xi3", DSValues.Xi3,
          "Xl2",  DSValues.Xl2, "Xl3", DSValues.Xl3);
      fprintf(SGP4File, 
          "%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f%7s%13.7f\n",
          "Xl4", DSValues.Xl4, "Xli",  DSValues.Xli,  "Xlamo", DSValues.Xlamo,
          "Xni", DSValues.Xni, "Zmol", DSValues.Zmol, "Zmos",  DSValues.Zmos);
    }
}
void TwoLine2RV
    (
      char* LongStr1, char* LongStr2, char Show, ElSetRec SatRec, 
      Vector& Ro, Vector& Vo, FILE* OutFile, char* SatName, char* LDate,
      char* DDate, char *SatCatFileName
    )
{
}

void TwoLine2RVSGP4
    (
      char* LongStr1, char* LongStr2, char Show, ElSetRec SatRec, 
      Vector& Ro, Vector& Vo, FILE* OutFile, char* SatName, char* LDate,
      char* DDate, char *SatCatFileName
    )
{
}
