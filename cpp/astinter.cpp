/*----------------------------------------------------------------      


                               UNIT ASTINTER;


    This file contains fundamental Astrodynamic procedures and functions     
    relating to interplanetary calculations.                                 

                            Companion code for                               
               Fundamentals of Astrodyanmics and Applications                
                                     2001                                    
                              by David Vallado                               

       (H)               email valladodl@worldnet.att.net                    
       (W) 303-344-6037, email davallado@west.raytheon.com                   

       *****************************************************************     

    Current :                                                                
              14 May 01  David Vallado                                       
                           2nd edition baseline                              
    Changes :                                                                
              23 Nov 87  David Vallado                                       
                           Original baseline                                 

       ----------------------------------------------------------------      

                                IMPLEMENTATION

       -----------------------------------------------------------------      */
#include <math.h>
#include <stdio.h>

#include "astmath.h"
#include "asttime.h"
#include "ast2body.h"

/*----------------------------------------------------------------------------
|
|                           PROCEDURE INTERPLANETARY
|
|  This PROCEDURE calculates the delta v's for an interplanetary mission using a
|    patched conic approximation.  The transfer assumes circular orbits for each
|    of the planets.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units        Ex Value
|    R1          - Radius of planet 1 from sun    km
|    R2          - Radius of planet 2 from sun    km
|    Rbo         - Radius at burnout at planet 1  km
|    Rimpact     - Radius at impact on planet 2   km
|    Mu1         - Grav parameter of planet 1     km3/s2
|    Mut         - Grav parameter of planet Sun   km3/s2
|    Mu2         - Grav parameter of planet 2     km3/s2
|
|  OutPuts       :
|    DeltaV1     - Hyperb Exc vel at planet 1 SOI km/s
|    DeltaV2     - Hyperb Exc vel at planet 2 SOI km/s
|    Vbo         - Burnout vel at planet 1        km/s
|    Vretro      - Retro vel at surface planet 2  km/s
|
|  Locals        :
|    SME1        - Spec Mech Energy of 1st orbit  Km2/s
|    SMEt        - Spec Mech Energy of trans orbitKm2/s
|    SME2        - Spec Mech Energy of 2nd orbit  Km2/s
|    Vcs1        - Vel of 1st orbit at dv 1 point Km/s
|    Vcs2        - Vel of 2nd orbit at dv 2 point Km/s
|    Vt1         - Vel of Trans orbit at dv 1 pnt Km/s
|    Vt2         - Vel of Trans orbit at dv 2 pnt Km/s
|    A           - Semimajor Axis of Trans orbit  Km
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001,
|
 ---------------------------------------------------------------------------*/
void Interplanetary
    (
      Real R1, Real R2, Real Rbo, Real Rimpact, Real Mu1, Real Mut, Real Mu2,
      Real& DeltaV1, Real& DeltaV2, Real& Vbo, Real& Vretro
    )
{
  Real SME1, SME2, SMEt, Vcs1, Vcs2, Vt1, Vt2, A, TP;

  /* - Find a, SME, apogee and perigee velocities of trans orbit -- */
  A    = (R1 + R2) * 0.5;
  SMEt = -Mut / (2.0 * A);
  Vt1  = sqrt(2.0 * ((Mut / R1) + SMEt));
  Vt2  = sqrt(2.0 * ((Mut / R2) + SMEt));
  
  /* ---- Find circular velocities of launch and target planet ---- */
  Vcs1 = sqrt(Mut / R1);
  Vcs2 = sqrt(Mut / R2);

  /* ---- Find delta velocities for Hohmann transfer portion  ----- */
  DeltaV1 = fabs(Vt1  - Vcs1);
  DeltaV2 = fabs(Vcs2 - Vt2);

  /* - Find SME and burnout/impact vel of launch / target planets -*/
  SME1   = DeltaV1 * DeltaV1 * 0.5;
  SME2   = DeltaV2 * DeltaV2 * 0.5;
  Vbo    = sqrt(2.0 * ((Mu1 / Rbo) + SME1));
  Vretro = sqrt(2.0 * ((Mu2 / Rimpact) + SME2));

  if ((Show == 'Y') || (Show == 'S'))
  {
    TP = PI * sqrt(A * A * A / Mut);  // Transfer Period in secs
    if (FileOut != NULL)
    {
      fprintf(FileOut, "     Transfer Period = %8.3f yrs  or  %8.3f days\n",
                       TP / 3.1536E07, TP / 86400.0);
    }
  }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE PLANETRV
|
|  This PROCEDURE calculate the planetary ephemerides using the Epoch J2000.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|                                               
|  Inputs          Description                    Range / Units
|    PlanetNum   - Number of planet               1..9
|    JD          - Julian Date                    days from 4713 BC
|
|  OutPuts       :
|    R           - XYZ position vector            AU
|    V           - XYZ velocity vector            AU / day
|
|  Locals        :
|    ArgLat      -
|    TrueLon     -
|    LonPer      -
|    TUUT1       - Julian cenuries from Jan 1, 2000
|    TU2         - Tu squared
|    TU3         - TU Cubed
|    N           -
|    obliquity   - angle between ecliptic and
|                    Earth equator                rad
|    a           - Semi or axis
|    ecc         - eccentricity
|    p           - semi-parameter
|    incl        - inclination
|    omega       - ecliptic long of ascending node
|    argp        - ecliptic arg of perihelion
|    nu          - ecliptic true anomaly
|    m           - ecliptic mean anomaly
|    LLong       - True longitude
|    LongP       - longitude of perihelion
|    e0          -
|
|  Coupling      :
|    RealMOD     -
|    LnCom1      -
|    NewtonM     -
|    RandV       -
|    ROT1        -
|
|  References    :
|    Vallado       2001, 281-285, Alg 33, Ex 5-5
|
 ---------------------------------------------------------------------------*/
void PlanetRV(Integer PlanetNum, char *WhichEpoch, char *WhichCoord, Real JD,
              Vector& R, Vector& V)
{
  const Real TwoPi = 2.0 * PI;
  const Real Rad   = 180.0 / PI;

  Real TUDaySun, ArgLat, TrueLon, LonPer, Tut1, Tut12, Tut13, n, Eps, 
       a, ecc, p, incl, Omega, ArgP, Nu, LLong, LongP, M, E0;
  Vector Re(3), Ve(3);

  /* ----------------- Approximate TTDB with TUT1 ----------------- */
  Tut1  = (JD - 2451545.0) / 36525.0;
  Tut12 = Tut1 * Tut1;
  Tut13 = Tut12 * Tut1;

  if (strcmp(WhichEpoch, "J20") == 0)
    switch (PlanetNum)
    {
      case 1:   /* -----------Mercury --------- */
        a     =   0.387098310;
        ecc   =   0.20563175 + 0.000020406 * Tut1 - 
                  0.0000000284 * Tut12 -0.00000000017 * Tut13;
        incl  =   7.004986   - 0.0059516 * Tut1   + 
                  0.00000081 * Tut12 + 0.000000041 * Tut13;
        LongP =  77.456119   + 0.1588643 * Tut1   - 
                  0.00001343 * Tut12 + 0.000000039 * Tut13;
        Omega =  48.330893   - 0.1254229 * Tut1   - 
                  0.00008833 * Tut12 - 0.000000196 * Tut13;
        LLong = 252.250906   + 149472.6746358 * Tut1 - 
                  0.00000535 * Tut12  + 0.000000002 * Tut13;
        break;
      case 2:   /* -----------Venus  ---------- */
        a     =   0.723329820;
        ecc   =   0.00677188 - 0.000047766 * Tut1 + 
                  0.0000000975 * Tut12 + 0.00000000044 * Tut13;
        incl  =   3.394662 - 0.0008568 * Tut1 - 
                  0.00003244 * Tut12  + 0.000000010 * Tut13;
        LongP = 131.563707 + 0.0048646 * Tut1 - 
                  0.00138232 * Tut12 - 0.000005332 * Tut13;
        Omega =  76.679920 - 0.2780080 * Tut1 - 
                  0.00014256 * Tut12  - 0.000000198 * Tut13;
        LLong = 181.979801 + 58517.8156760 * Tut1 + 
                  0.00000165 * Tut12 - 0.000000002 * Tut13;
        break;
      case 3:   /* -----------Earth  ---------- */
        a     =   1.000001018;
        ecc   =   0.01670862 - 0.000042037 * Tut1 - 
                  0.0000001236 * Tut12 + 0.00000000004 * Tut13;
        incl  =   0.0000000 + 0.0130546 * Tut1 - 
                  0.00000931 * Tut12  - 0.000000034 * Tut13;
        LongP = 102.937348 + 0.3225557 * Tut1 + 
                  0.00015026 * Tut12  + 0.000000478 * Tut13;
        Omega =   0.0;
        LLong = 100.466449 + 35999.3728519 * Tut1 - 
                  0.00000568 * Tut12 + 0.000000000 * Tut13;
        break;
      case 4:   /* -----------Mars   ---------- */
        a     =   1.523679342;
        ecc   =   0.09340062 + 0.000090483 * Tut1 - 
                  0.0000000806 * Tut12  - 0.00000000035 * Tut13;
        incl  =   1.849726 - 0.0081479 * Tut1 - 
                  0.00002255 * Tut12  - 0.000000027 * Tut13;
        LongP = 336.060234 + 0.4438898 * Tut1 - 
                  0.00017321 * Tut12 + 0.000000300 * Tut13;
        Omega =  49.558093 - 0.2949846 * Tut1 - 
                  0.00063993 * Tut12 - 0.000002143 * Tut13;
        LLong = 355.433275 + 19140.2993313 * Tut1 + 
                  0.00000261 * Tut12 - 0.000000003 * Tut13;
        break;
      case 5:   /* -----------Jupiter --------- */
        a     =   5.202603191  + 0.0000001913 * Tut1;
        ecc   =   0.04849485 + 0.000163244 * Tut1 - 
                  0.0000004719 * Tut12 - 0.00000000197 * Tut13;
        incl  =   1.303270 - 0.0019872 * Tut1 + 
                  0.00003318 * Tut12  + 0.000000092 * Tut13;
        LongP =  14.331309 + 0.2155525 * Tut1 + 
                  0.00072252 * Tut12 - 0.000004590 * Tut13;
        Omega = 100.464441 + 0.1766828 * Tut1 + 
                  0.00090387 * Tut12 - 0.000007032 * Tut13;
        LLong =  34.351484 + 3034.9056746 * Tut1 - 
                  0.00008501 * Tut12  + 0.000000004 * Tut13;
        break;
      case 6:   /* -----------Saturn  --------- */
        a     =   9.554909596 - 0.0000021389 * Tut1;
        ecc   =   0.05550862 - 0.000346818 * Tut1 - 
                  0.0000006456 * Tut12 + 0.00000000338 * Tut13;
        incl  =   2.488878 + 0.0025515 * Tut1 - 
                  0.00004903 * Tut12 + 0.000000018 * Tut13;
        LongP =  93.056787 + 0.5665496 * Tut1 + 
                  0.00052809 * Tut12 + 0.000004882 * Tut13;
        Omega = 113.665524 - 0.2566649 * Tut1 - 
                  0.00018345 * Tut12 + 0.000000357 * Tut13;
        LLong =  50.077471 + 1222.1137943 * Tut1 + 
                  0.00021004 * Tut12 - 0.000000019 * Tut13;
        break;
      case 7:   /* -----------Uranus  --------- */
        a     =  19.218446062 - 0.0000000372 * Tut1 + 
                  0.00000000098 * Tut12;
        ecc   =   0.04629590 - 0.000027337 * Tut1 + 
                  0.0000000790 * Tut12 + 0.00000000025 * Tut13;
        incl  =   0.773196 - 0.0016869 * Tut1 + 
                  0.00000349 * Tut12 + 0.000000016 * Tut13;
        LongP = 173.005159 + 0.0893206 * Tut1 - 
                  0.00009470 * Tut12 + 0.000000413 * Tut13;
        Omega =  74.005947 + 0.0741461 * Tut1 + 
                  0.00040540 * Tut12 + 0.000000104 * Tut13;
        LLong = 314.055005 + 428.4669983 * Tut1 - 
                  0.00000486 * Tut12  + 0.000000006 * Tut13;
        break;
      case 8:   /* -----------Neptune --------- */
        a     =  30.110386869 - 0.0000001663 * Tut1 + 
                  0.00000000069 * Tut12;
        ecc   =   0.00898809 + 0.000006408 * Tut1 - 
                  0.0000000008 * Tut12;
        incl  =   1.769952 + 0.0002257 * Tut1 + 
                  0.00000023 * Tut12  - 0.000000000 * Tut13;
        LongP =  48.123691 + 0.0291587 * Tut1 + 
                  0.00007051 * Tut12 - 0.000000000 * Tut13;
        Omega = 131.784057 - 0.0061651 * Tut1 - 
                  0.00000219 * Tut12 - 0.000000078 * Tut13;
        LLong = 304.348665 + 218.4862002 * Tut1 + 
                  0.00000059 * Tut12 - 0.000000002 * Tut13;
        break;
      case 9:   /* -----------Pluto  ---------- */
        a     =  39.53758;
        ecc   =   0.250877;
        incl  =  17.13233;
        LongP = 224.6148;
        Omega = 110.4065;
        LLong = 218.88735;
        break;
      default:  /* ----------bad index ---------*/
        a     = 0.0;
        ecc   = 0.0;
        incl  = 0.0;
        LongP = 0.0;
        Omega = 0.0;
        LLong = 0.0;
        break;
    }
  if (strcmp(WhichEpoch, "ODA") == 0)
    /* ------------ Mean equinox of date in degrees (XYZ) ----------- */
    switch (PlanetNum)
    {
      case 1:   /* -----------Mercury --------- */
        a     =   0.387098310;
        ecc   =   0.20563175 +    0.000020406 * Tut1 - 
                  0.0000000284 * Tut12 - 0.00000000017 * Tut13;
        incl  =   7.004986 +      0.0018215 * Tut1 - 
                  0.00001809 * Tut12 + 0.000000053 * Tut13;
        LongP =  77.456119 +      1.5564775 * Tut1 + 
                  0.00029589 * Tut12 + 0.000000056 * Tut13;
        Omega =  48.330893 +      1.1861890 * Tut1 + 
                  0.00017587 * Tut12 + 0.000000211 * Tut13;
        LLong = 252.250906 + 149474.0722491 * Tut1 + 
                  0.00030397 * Tut12  + 0.000000018 * Tut13;
        break;
      case 2:   /* -----------Venus  ---------- */
        a     =   0.723329820;
        ecc   =   0.00677188 -   0.000047766 * Tut1 + 
                  0.0000000975 * Tut12 + 0.00000000044 * Tut13;
        incl  =   3.394662 +     0.0010037 * Tut1 - 
                  0.00000088 * Tut12  - 0.000000007 * Tut13;
        LongP = 131.563707 +     1.4022188 * Tut1 - 
                  0.00107337 * Tut12 - 0.000005315 * Tut13;
        Omega =  76.679920 +     0.9011190 * Tut1 + 
                  0.00040665 * Tut12  - 0.000000080 * Tut13;
        LLong = 181.979801 + 58519.2130302 * Tut1 + 
                  0.00031060 * Tut12 + 0.000000015 * Tut13;
        break;
      case 3:   /* -----------Earth  ---------- */
        a     =   1.000001018;
        ecc   =   0.01670862 -   0.000042037 * Tut1 - 
                  0.0000001236 * Tut12 + 0.00000000004 * Tut13;
        incl  =   0.0;
        LongP = 102.937348 +     1.7195269 * Tut1 + 
                  0.00045962 * Tut12  + 0.000000499 * Tut13;
        LLong = 100.466449 + 36000.7698231 * Tut1 + 
                  0.00030368 * Tut12 + 0.000000021 * Tut13;
        break;
      case 4:   /* -----------Mars   ---------- */
        a     =   1.523679342;
        ecc   =   0.09340062 +   0.000090483 * Tut1 - 
                  0.0000000806 * Tut12 - 0.00000000035 * Tut13;
        incl  =   1.849726 -     0.0006010 * Tut1 + 
                  0.00001276 * Tut12  - 0.000000006 * Tut13;
        LongP = 336.060234 +     1.8410331 * Tut1 + 
                  0.00013515 * Tut12 + 0.000000318 * Tut13;
        Omega =  49.558093 +     0.7720923 * Tut1 + 
                  0.00001605 * Tut12 + 0.000002325 * Tut13;
        LLong = 355.433275 + 19141.6964746 * Tut1 + 
                  0.00031097 * Tut12 + 0.000000015 * Tut13;
        break;
      case 5:   /* -----------Jupiter --------- */
        a     =   5.202603191 + 0.0000001913 * Tut1;
        ecc   =   0.04849485 +  0.000163244 * Tut1 - 
                  0.0000004719 * Tut12 - 0.00000000197 * Tut13;
        incl  =   1.303270 -    0.0054966 * Tut1 + 
                  0.00000465 * Tut12  - 0.000000004 * Tut13;
        LongP =  14.331309 +    1.6126668 * Tut1 + 
                  0.00103127 * Tut12 - 0.000004569 * Tut13;
        Omega = 100.464441 +    1.0209550 * Tut1 + 
                  0.00040117 * Tut12 + 0.000000569 * Tut13;
        LLong =  34.351484 + 3036.3027889 * Tut1 + 
                  0.00022374 * Tut12  + 0.000000025 * Tut13;
        break;
      case 6:   /* -----------Saturn  --------- */
        a     =   9.554909596 - 0.0000021389 * Tut1;
        ecc   =   0.05550862 -  0.000346818 * Tut1 - 
                  0.0000006456 * Tut12 + 0.00000000338 * Tut13;
        incl  =   2.488878 -    0.0037363 * Tut1 - 
                  0.00001516 * Tut12 + 0.000000089 * Tut13;
        LongP =  93.056787 +    1.9637694 * Tut1 + 
                  0.00083757 * Tut12 + 0.000004899 * Tut13;
        Omega = 113.665524 +    0.8770979 * Tut1 - 
                  0.00012067 * Tut12 - 0.000002380 * Tut13;
        LLong =  50.077471 + 1223.5110141 * Tut1 + 
                  0.00051952 * Tut12 - 0.000000003 * Tut13;
        break;
      case 7:   /* -----------uranus  --------- */
        a     =  19.218446062 - 0.0000000372 * Tut1  + 
                  0.00000000098 * Tut12;
        ecc   =   0.04629590 -  0.000027337 * Tut1 + 
                  0.0000000790 * Tut12 + 0.00000000025 * Tut13;
        incl  =   0.773196 +   0.0007744 * Tut1 + 
                  0.00003749 * Tut12 - 0.000000092 * Tut13;
        LongP = 173.005159 +   1.4863784 * Tut1 + 
                  0.00021450 * Tut12 + 0.000000433 * Tut13;
        Omega =  74.005947 +   0.5211258 * Tut1 + 
                  0.00133982 * Tut12 + 0.000018516 * Tut13;
        LLong = 314.055005 + 429.8640561 * Tut1 + 
                  0.00030434 * Tut12  + 0.000000026 * Tut13;
        break;
      case 8:   /* -----------Neptune --------- */
        a     =  30.110386869 - 0.0000001663 * Tut1 + 
                  0.00000000069 * Tut12;
        ecc   =   0.00898809 +  0.000006408 * Tut1 -  
                  0.0000000008 * Tut12;
        incl  =   1.769952 -   0.0093082 * Tut1 - 
                  0.00000708 * Tut12  + 0.000000028 * Tut13;
        LongP =  48.123691 +   1.4262677 * Tut1 + 
                  0.00037918 * Tut12 - 0.000000003 * Tut13;
        Omega = 131.784057 +   1.1022057 * Tut1 + 
                  0.00026006 * Tut12 - 0.000000636 * Tut13;
        LLong = 304.348665 + 219.8833092 * Tut1 + 
                  0.00030926 * Tut12 + 0.000000018 * Tut13;
        break;
      case 9:   /* -----------Pluto  ---------- */
        a     =  39.53758;
        ecc   =   0.250877;
        incl  =  17.13233;
        LongP = 224.6148;
        Omega = 110.4065;
        LLong = 218.88735;
        break;
      default:
        a     = 0.0;
        ecc   = 0.0;
        incl  = 0.0; 
        LongP = 0.0;
        Omega = 0.0;
        LLong = 0.0;
        break;
    }
  
  incl  = incl / Rad;  // Convert to radians
  LongP = LongP / Rad;
  Omega = Omega / Rad;
  LLong = LLong / Rad;

  LLong = Mod(LLong, TwoPi );
  LongP = Mod(LongP, TwoPi );
  Omega = Mod(Omega, TwoPi );
  ArgP  = LongP - Omega;
  M     = LLong - LongP;

  NewtonM(ecc, M, E0, Nu);

  /* ------------- Find Heliocentric ecliptic r and v ------------- */
  p       = a * (1.0- ecc* ecc);
  ArgLat  = ArgP + Nu;
  TrueLon = Omega + ArgP + Nu;
  LonPer  = Omega + ArgP + PI;
  RandV(p, ecc, incl, Omega, ArgP, Nu, ArgLat, TrueLon, LonPer, R,V );

  /* ---- Correct the velocity because we used TTdb - days!! ------ */
  TUDaySun = 1.0 / 58.1324409;  // 1.0 / days per sun TU
  V = TUDaySun * V;

  if (strcmp(WhichCoord, "GEO") == 0)
  {
    /* ------------ Find obliquity of the ecliptic angle -------- */
    Eps = 23.439291 - 0.0130042*Tut1 - 0.000000164*Tut12 +  0.000000504*Tut13;
    Eps = Mod(Eps, 360.0);
    Eps = Eps / Rad;

    /* -------------- Rotate to Geocentric coordinates ---------- */
    R = R.Rot1(-Eps);
    V = V.Rot1(-Eps);
  }

  if (Show == 'S')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "New Case, tut1 =%11.7f  Planet %3d %s %s----------- \n",
                       Tut1, PlanetNum, WhichCoord, WhichEpoch);
      fprintf(FileOut, "r helio %11.6f %11.6f %11.6f %11.6f\n",
                        R.Get(1), R.Get(2), R.Get(3), R.Mag());
      fprintf(FileOut, "%11.6f %11.6f %11.6f %11.6f\n",
                        V.Get(1), V.Get(2), V.Get(3), V.Mag());
      fprintf(FileOut, "%11.1f %11.1f %11.1f km",
                       R.Get(1) * 149597870.0, R.Get(2) * 149597870.0, 
                       R.Get(3) * 149597870.0, R.Mag() * 149597870.0);
      fprintf(FileOut, "%11.6f %11.6f %11.6f",
                        V.Get(1) * 29.784691674, V.Get(2) * 29.784691674, 
                        V.Get(3) * 29.784691674, V.Mag() * 29.784691674);
    }
  if (Show == 'Y')
    if (FileOut != NULL)
      fprintf(FileOut, "%5d %10.6f %13.6f %10.6f %12.6f %14.8f\n",
                      WhichCoord, a, ecc, incl * Rad, Omega * Rad, LongP * Rad);
}
