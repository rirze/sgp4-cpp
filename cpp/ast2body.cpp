/*     ----------------------------------------------------------------      


                               UNIT AST2BODY;


    This file contains fundamental Astrodynamic procedures and functions     
    using 2-body dynamics. The routines span a wide range of material, and   
    they come from chapters 2, 3, 5, and 11.                                 

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
                           Original Baseline                                 

       ----------------------------------------------------------------      

                                  IMPLEMENTATION

       ----------------------------------------------------------------      */
#include <math.h>
#include <stdio.h>

#include "ast2body.h"

extern char Show;
extern FILE *FileOut;

/*------------------------------------------------------------------------------
|
|                           PROCEDURE CHECKHITEARTH
|
|  This PROCEDURE checks to see IF the trajectory hits the earth during an 
|    orbital transfer.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInt        - Initial Position vector of Int ER
|    V1t         - Initial Velocity vector of trnsER/TU
|    RTgt        - Initial Position vector of Tgt ER
|    V2t         - Final Velocity vector of trns  ER/TU
|
|  Outputs       :
|    HitEarth    - Is Earth was impacted          'Y' 'N'
|
|  Locals        :
|    SME         - Specific mechanical energy
|    rp          - Radius of Perigee              ER
|    TransA      - Semi-or axis of transfer       ER
|    TransE      - Eccentricity of transfer
|    TransP      - Semi-paramater of transfer     ER
|    HBar        - Angular momentum vector of
|                  transfer orbit
|
|  Coupling      :
|    DOT         - DOT product of vectors
|    MAG         - Magnitude of a vector
|    CROSS       - CROSS product of vectors
|
|  References    :
|    Vallado       2001, 472-474, Alg 57
|
  ----------------------------------------------------------------------------*/
void CheckHitEarth
    (
      Vector Rint, Vector V1t, Vector Rtgt, Vector V2t, char& HitEarth
    )
{
  Vector HBar(3);
  Real   SME, rp, TransP, TransA, TransE;

  HitEarth = 'N';

  /* ----------- Find IF trajectory intersects Earth -------------- */
  if ((Rint.Dot(V1t) < 0.0) && (Rtgt.Dot(V2t) > 0.0))
  {
    /* ----------------  Find H N and E vectors   ---------------- */
    HBar = Rint.Cross(V1t);

    if (Rint.Mag() > 0.00001)
    {
      SME    = V1t.Mag() * V1t.Mag() * 0.5 - (1.0 / Rint.Mag());
      TransP = Rint.Mag() * Rint.Mag();
      TransE = 1.0;
      if (fabs(SME) > 0.00001)
      {
        TransA = -1.0 / (2.0 * SME);
        TransE = sqrt((TransA - TransP) / TransA );
        rp     = TransA * (1.0 - TransE);
      }
      else
        rp = TransP * 0.5;  // Parabola
    }
    else
      printf("The orbit does not exist \n");

    if (fabs(rp) < 1.0)
      HitEarth = 'Y';
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE ELORB
|
|  This PROCEDURE finds the classical orbital elements given the Geocentric
|    Equatorial Position and Velocity vectors.  Special cases for equatorial
|    and circular orbits are also handled.  IF the elements are Infinite, they
|    are set to 999999.9. IF elements are Undefined, they are set to 999999.1.
|    Be sure to check for these during output or subsequent use!!
|
|  Algorithm     : Initialze variables
|                  IF the HBar magnitude exists, contiNue, otherwise exit and
|                       assign undefined values
|                    Find vectors and values
|                    Determine the type of orbit with IF statements
|                    Find angles depending on the orbit type
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - IJK Position vector            ER
|    V           - IJK Velocity vector            ER / TU
|
|  Outputs       :
|    P           - SemiLatus rectum               ER
|    A           - semimajor axis                 ER
|    Ecc         - Eccentricity
|    Incl        - inclination                    0.0 to Pi rad
|    Omaga       - Longitude of Ascending Node    0.0 to 2Pi rad
|    Argp        - Argument of Perigee            0.0 to 2Pi rad
|    Nu          - True anomaly                   0.0 to 2Pi rad
|    M           - Mean anomaly                   0.0 to 2Pi rad
|    ArgLat      - Argument of Latitude      (CI) 0.0 to 2Pi rad
|    LamTrue     - True Longitude            (CE) 0.0 to 2Pi rad
|    LonPer      - Longitude of Periapsis    (EE) 0.0 to 2Pi rad
|
|  Locals        :
|    HBar        - Angular Momentum H Vector      ER2 / TU
|    EBar        - Eccentricity     E Vector
|    NBar        - Line of Nodes    N Vector
|    c1          - V**2 - u/R
|    RDotV       - R DOT V
|    Hk          - Hk unit vector
|    SME         - Specfic Mechanical Energy      ER2 / TU2
|    i           - index
|    E           - Eccentric, Parabolic,
|                  Hyperbolic Anomaly             rad
|    Temp        - Temporary variable
|    TypeOrbit   - Type of orbit                  EE, EI, CE, CI
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    CROSS       - CROSS product of two vectors
|    DOT         - DOT product of two vectors
|    ARCCOS      - Arc Cosine FUNCTION
|    ANGLE       - Find the ANGLE between two vectors
|    NEWTONNU    - Find the mean anomaly
|    REALMOD     - MOD FUNCTION for REAL variables
|
|  References    :
|    Vallado       2001, 118-122, Alg 9, Ex 2-5
|
 -----------------------------------------------------------------------------*/
void ElOrb
    (
      Vector R, Vector V,
      Real& P,  Real& A, Real& Ecc,    Real& Incl,    Real& Omega, Real& Argp,
      Real& Nu, Real& M, Real& ArgLat, Real& TrueLon, Real& LonPer
    )
{
  const Real TwoPi     = 2.0 * PI;
  const Real HalfPi    = 0.5 * PI;
  const Real Small     = 0.00000001;    // Small value for tolerances
  const Real Infinite  = 999999.9;      // Infinite value
  const Real Undefined = 999999.1;      // Undefined value

  char TypeOrbit[3];

  Vector EBar(3), HBar(3), NBar(3);

  Real   E, c1, RDotV, Hk, SME, Temp;

  /* -------------------  Find H N and E vectors   ---------------- */
  HBar = R.Cross(V);
  if (HBar.Mag() > Small)
  {
    NBar.Set(-HBar.Get(2), 1);
    NBar.Set( HBar.Get(1), 2);
    NBar.Set(0.0, 3);
    c1 = V.Mag() * V.Mag() - 1.0 / R.Mag();
    RDotV = R.Dot(V);
    for (UINT i = 1; i <= 3; i++)
      EBar.Set(c1 * R.Get(i) - RDotV * V.Get(i), i);

    /* -------------  Find a e and semi-Latus rectum   ---------- */
    SME = V.Mag() * V.Mag() * 0.5 - 1.0 / R.Mag();
    if (fabs(SME) > Small)
      A = -1.0 / (2.0 * SME);
    else
      A = Infinite;       // Parabola
    Ecc = EBar.Mag();
    P = HBar.Mag() * HBar.Mag();

    /* ------------------  Find inclination   ------------------- */
    Hk = HBar.Get(3) / HBar.Mag();
    if (fabs(fabs(Hk) - 1.0) < Small)
      /* -------  Equatorial Orbits   ---------- */
      if (fabs(HBar.Get(3)) > 0.0)
        Hk = Sgn(HBar.Get(3)) * 1.0;
    Incl = acos(Hk);

    /* ---------  Determine type of orbit for Later use  -------- */
    /* ------- Elliptical, Parabolic, Hyperbolic Inclined ------- */
    strcpy(TypeOrbit, "EI");
    if (Ecc < Small)
      /* --------------  Circular Equatorial -------------- */
      if ((Incl < Small) || (fabs(Incl - PI) < Small))
        strcpy(TypeOrbit, "CE");
      else
        /* ------------  Circular Inclined -------------- */
        strcpy(TypeOrbit, "CI");
    else
      /* -- Elliptical, Parabolic, Hyperbolic Equatorial -- */
      if ((Incl < Small) || (fabs(Incl - PI) < Small))
        strcpy(TypeOrbit, "EE");

    /* -----------  Find Longitude of Ascending Node ------------ */
    if (NBar.Mag() > Small)
    {
      Temp = NBar.Get(1) / NBar.Mag();
      if (fabs(Temp) > 1.0)
        Temp = Sgn(Temp) * 1.0;
      Omega = acos(Temp);
      if (NBar.Get(2) < 0.0)
        Omega = TwoPi - Omega;
    }
    else
      Omega = Undefined;

    /* ----------------- Find Argument of perigee --------------- */
    if (strcmp(TypeOrbit, "EI") == 0)
    {
      Argp = NBar.Angle(EBar);
      if (EBar.Get(3) < 0.0)
        Argp = TwoPi - Argp;
    }
    else
      Argp = Undefined;

    /* -------------  Find True Anomaly at Epoch    ------------- */
    if (TypeOrbit[0] == 'E')
    {
      Nu = EBar.Angle(R);
      if (RDotV < 0.0)
        Nu = TwoPi - Nu;
    }
    else
      Nu = Undefined;

    /* -----  Find Argument of Latitude - Circular Inclined ----- */
    if (strcmp(TypeOrbit, "CI") == 0)
    {
      ArgLat = NBar.Angle(R);
      if (R.Get(3) < 0.0)
        ArgLat = TwoPi - ArgLat;
    }
    else
      ArgLat = Undefined;

    /* --- Find Longitude of Perigee - Elliptical Equatorial ---- */
    if ((EBar.Mag() > Small) && (strcmp(TypeOrbit, "EE") == 0))
    {
      Temp = EBar.Get(1) / EBar.Mag();
      if (fabs(Temp) > 1.0)
        Temp = Sgn(Temp) * 1.0;
      LonPer = acos(Temp);
      if (EBar.Get(2) < 0.0)
        LonPer = TwoPi - LonPer;
      if (Incl > HalfPi)
        LonPer = TwoPi - LonPer;
    }
    else
      LonPer = TwoPi - LonPer;

    /* ------- Find True Longitude - Circular Equatorial -------- */
    if ((R.Mag() > Small) && (strcmp(TypeOrbit, "OE") == 0))
    {
      Temp = R.Get(1) / R.Mag();
      if (fabs(Temp) > 1.0)
        Temp = Sgn(Temp) * 1.0;
      TrueLon = acos(Temp);
      if (R.Get(2) < 0.0)
        TrueLon = TwoPi - TrueLon;
      if (Incl > HalfPi)
        TrueLon = TwoPi - TrueLon;
    }
    else
      TrueLon = Undefined;

    /* ------------- Find Mean Anomaly for all orbits ----------- */
    NewtonU(Ecc, Nu, E, M);

    if (Show == 'Y')
    {
      if (FileOut != NULL)
      {
        fprintf(FileOut, "H =   %13.7f  %14.7f %14.7f %14.7f\n",
                          HBar.Get(1), HBar.Get(2), HBar.Get(3), HBar.Mag());
        fprintf(FileOut, "N =   %13.7f  %14.7f %14.7f %14.7f\n",
                          NBar.Get(1), NBar.Get(2), NBar.Get(3), NBar.Mag());
        fprintf(FileOut, "E =   %13.7f  %14.7f %14.7f %14.7f\n",
                          EBar.Get(1), EBar.Get(2), EBar.Get(3), EBar.Mag());
        fprintf(FileOut, "SME = %13.7f ER2/TU2\n", SME);
        fprintf(FileOut, "Anomaly %14.7f %s\n", E * 57.29578, TypeOrbit);
      }
    }
  }
  else
  {
    P       = Undefined;
    A       = Undefined;
    Ecc     = Undefined;
    Incl    = Undefined;
    Omega   = Undefined;
    Argp    = Undefined;
    Nu      = Undefined;
    M       = Undefined;
    ArgLat  = Undefined;
    TrueLon = Undefined;
    LonPer  = Undefined;
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE FINDC2C3
|
|  This PROCEDURE calcuLates the C2 and C3 functions for use in the Universal
|    Variable calcuLation of z.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|                                                                4 Feb 1992
|  Inputs          Description                    Range / Units
|    ZNew        - Z variable                     rad2
|
|  Outputs       :
|    C2New       - C2 FUNCTION value
|    C3New       - C3 FUNCTION value
|
|  Locals        :
|    SqrtZ       - Square root of ZNew
|
|  Coupling      :
|    SINH        - Hyperbolic Sine
|    COSH        - Hyperbolic Cosine
|
|  References    :
|    Vallado       2001, 70-71, Alg 1
|
 -----------------------------------------------------------------------------*/
void FindC2C3(Real ZNew, Real& C2New, Real& C3New)
{
  const Real Small = 0.00000001;       // Small value for tolerances
  Real SqrtZ;

  if (ZNew > Small)
  {
    SqrtZ  = sqrt(ZNew);
    C2New  = (1.0 - cos(SqrtZ)) / ZNew;
    C3New  = (SqrtZ - sin(SqrtZ)) / (SqrtZ * SqrtZ * SqrtZ);
  }
  else
    if (ZNew < -Small)
    {
      SqrtZ  = sqrt(-ZNew);
      C2New  = (1.0 - cosh(SqrtZ)) / ZNew;
      C3New  = (sinh(SqrtZ) - SqrtZ) / (SqrtZ * SqrtZ * SqrtZ);
    }
    else
    {
      C2New = 0.5;
      C3New = 1.0 / 6.0;
    }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE FINDTOF
|
|  This PROCEDURE finds the time of flight given the initial position vectors,
|    Semi-parameter, and the sine and cosine values for the change in true
|    anomaly.  The result uses p-iteration theory to analytically find the result.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ro          - Interceptor position vector    ER
|    R           - TARGET position vector         ER
|    p           - Semiparameter                  ER
|
|  Outputs       :
|    Tof         - Time for transfer              TU
|
|  Locals        :
|    SinDNu      - Sine of change in Nu           rad
|    CosDNu      - Cosine of change in Nu         rad
|    DeltaE      -
|    DeltaH      -
|    k           -
|    l           -
|    m           -
|    a           -
|    f           -
|    g           -
|    FDot        -
|    SinDeltaE   - Sine value
|    CosDeltaE   - Cosine value
|    RcrossR     - CROSS product of two positions
|
|  Coupling      :
|    CROSS       - CROSS product of two vectors
|    SINH        - Hyperbolic Sine
|    ARCCOSH     - Arc hyperbolic cosine
|    ATAN2       - Arc tangent FUNCTION which also resloves quadrants
|    POWER       - Raise a Number to some POWER
|
|  References    :
|    Vallado       2001, 130-133, Alg 11
|
 -----------------------------------------------------------------------------*/
void FindTOF(Vector Ro, Vector R, Real p, Real& Tof)
{
  const  Real Small = 0.00001;
  Vector RCrossR(3);
  Real   CosDNu, SinDNu, c , s, alpha, DeltaE, DeltaH, DNu,
         k, l, m, a, f, g, FDot, SinDeltaE, CosDeltaE;

  CosDNu  = Ro.Dot(R) / Ro.Mag() * Ro.Mag();
  RCrossR = Ro.Cross(R);
  SinDNu  = RCrossR.Mag() / (Ro.Mag() * Ro.Mag());

  k = Ro.Mag() * R.Mag() * (1.0 - CosDNu);
  l = Ro.Mag() + R.Mag();
  m = Ro.Mag() * R.Mag() * (1.0 + CosDNu);
  a = (m * k * p) / ((2.0 * m - l * l) * p * p + 2.0 * k * l * p - k * k);

  /* -------  Use F and G series to find Velocity Vectors  -------- */
  f = 1.0 - (R.Mag() / p) * (1.0 - CosDNu);
  g = Ro.Mag() * R.Mag() * SinDNu / sqrt(p);
  alpha = 1.0 / a;
  if (alpha > Small)
  {
    /* -------------------- Elliptical ------------------ */
    DNu  = atan2(SinDNu, CosDNu);
    FDot = sqrt(1.0 / p) * tan(DNu * 0.5) * (((1.0 - CosDNu) / p)-
           (1.0 /Ro.Mag()) - (1.0 / R.Mag()));
    CosDeltaE = 1.0 - (Ro.Mag() * R.Mag() * FDot) / sqrt(a);
    DeltaE    = atan2( SinDeltaE, CosDeltaE);
    Tof       = g + sqrt(a * a * a) * (DeltaE - SinDeltaE);
    if (Show == 'Y')
    {
      if (FileOut != NULL)
      {
        fprintf(FileOut, "f g fgdot %11.7f %11.7f %11.7f\n", f, g, FDot);
        fprintf(FileOut, "DNu  %11.7f %11.7f\n", 
                          DNu * 57.29577951, DeltaE * 57.2957795);
      }
    }
  }
  else
  {
    /* -------------------- Hyperbolic ------------------ */
    if (alpha < -Small)
    {
      DeltaH  = acosh(1.0 - (Ro.Mag() / a) * (1.0 - f));
      Tof     = g + sqrt(-a * a * a) * (sinh(DeltaH) - DeltaH);
    }
    else
    {
      /* ----------------- Parabolic ------------------ */
      DNu = atan2(SinDNu, CosDNu);
      c   = sqrt(R.Mag() * R.Mag() + Ro.Mag() * Ro.Mag() - 
                 2.0 * R.Mag() * Ro.Mag() * cos(DNu));
      s   = (Ro.Mag() + R.Mag() + c) * 0.5;
      Tof = (2.0 / 3.0) * sqrt(s * s * s * 0.5) * 
            (1.0 - Power((s - c) / s, 1.5));
      printf("Inside parabolic %11.7f\n", Tof * 13.44685);
    }
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE GEOC_GEOD
|
|  This PROCEDURE converts from Geodetic to Geocentric Latitude for positions
|    on the surface of the Earth.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Latgd       - Geodetic Latitude              -Pi to Pi rad
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    Latgc       - Geocentric Latitude            -Pi to Pi rad
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 146, Eq 3-11
|
 -----------------------------------------------------------------------------*/
void GeocGeod(Real& Latgc, Direction dir, Real& Latgd)
{
  const Real EESqrd = 0.006694385000;   // Eccentricity of Earth Sqrd

  if (dir == FROM)
    Latgc = atan((1.0 - EESqrd) * tan(Latgd));
  else
    Latgd = atan(tan(Latgc) / (1.0 - EESqrd));
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE IJKtoLATLON/A/B/E
|
|  These PROCEDUREs convert a Geocentric Equatorial (IJK) position vector into
|    latitude and longitude.  Geodetic and Geocentric latitude are found.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - IJK position vector            ER
|    JD          - Julian Date                    days from 4713 BC
|
|  OutPuts       :
|    Latgc       - Geocentric Latitude            -Pi to Pi rad
|    Latgd       - Geodetic Latitude              -Pi to Pi rad
|    Lon         - Longitude (WEST -)             -2Pi to 2Pi rad
|    Hellp       - Height above the ellipsoid     ER
|
|  Locals        :
|
|  Escobal:
|    Rc          - Range of SITE wrt earth center ER
|    Height      - Height above earth wrt SITE    ER
|    Alpha       - ANGLE from Iaxis to point, LST rad
|    OldDelta    - Previous value of DeltaLat     rad
|    DeltaLat    - Diff between Delta and
|                  Geocentric lat                 rad
|    Delta       - Declination ANGLE of R in IJK  rad
|    RSqrd       - Magnitude of r squared         ER2
|    SinTemp     - Sine of Temp                   rad
|
|  Almanac:
|    Temp        - Diff between Geocentric/
|                  Geodetic lat                   rad
|    GST         - Greenwich SIDEREAL time        rad
|    SinTemp     - Sine of Temp                   rad
|    OldDelta    - Previous value of DeltaLat     rad
|    RtAsc       - Right ascension                rad
|    Decl        - Declination                    rad
|    c           -
|    i           - index
|
|  Borkowski:
|
|
|
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ATAN2       - Arc Tangent which also resolves quadrant
|    ARCSIN      - Arc Sine of a value
|    GSTIME      - Greenwich SIDEREAL Time
|    GEOCGEOD    - Converts between geocentric and geodetic latitude
|    SGN         - Sign of a number +1 or -1
|
|  References    :
|    Vallado       2001, 174-179, Alg 12 and Alg 13, Ex 3-3
|
 -----------------------------------------------------------------------------*/
void IJKtoLatLonA
    (
      Vector R, Real JD, Real& Latgc, Real& Latgd, Real& Lon, Real& Hellp
    )
{
  const Real TwoPi  = 2.0 * PI;
  const Real Small  = 0.00000001;       // Small value for tolerances
  const Real EESqrd = 0.006694385000;   // Eccentricity of Earth Sqrd
  Integer i;
  Real RtAsc, OldDelta, c, Decl, Temp, GST, SinTemp;

  /* ---------------------  Initialize values   ------------------- */
  /* ------------------ Find Longitude value  --------------------- */
  Temp = sqrt(R.Get(1) * R.Get(1) + R.Get(2) * R.Get(2));
  if (fabs(Temp) < Small)
    RtAsc = Sgn(R.Get(3)) * PI * 0.5;
  else
    RtAsc = atan2(R.Get(2) / Temp, R.Get(1) / Temp);
  GST = GSTime(JD);
  Lon = RtAsc - GST;
  if (fabs(Lon) >= PI)   // Mod it ?
    if (Lon < 0.0)
      Lon = TwoPi + Lon;
    else
      Lon = Lon - TwoPi;
  Decl = asin(R.Get(3) / R.Mag());
  Latgd = Decl;
  if (Show == 'Y')
    printf("GST %14.7f %14.7f %14.7f %11.7\n", 
           GST * 57.29577, RtAsc * 57.29578, Decl * 57.29578, Temp);

  /* -------------- Iterate to find Geodetic Latitude ------------- */
  i = 1;
  while (1 == 1)
  {
    OldDelta = Latgd;
    SinTemp  = sin(Latgd);
    c        = 1.0 / (sqrt(1.0 - EESqrd * SinTemp * SinTemp));
    Latgd    = atan((R.Get(3) + c * EESqrd * SinTemp) / Temp);
    if (Show == 'Y')
      printf("A loops %3d  gd %16.10f %11.7f\n", i, Latgd * 57.29578, c);
    i++;
    if ((fabs(OldDelta - Latgd) < Small) || (i >= 10))
      break;
  }

  Hellp = (Temp / cos(Latgd)) - c;

  GeocGeod(Latgc, FROM, Latgd);
  
  if (i >= 10)
    printf("IJKtoLatLon did NOT converge in c2\n");

}

void IJKtoLatLonB
    (
      Vector R, Real JD, Real& Latgc, Real& Latgd, Real& Lon, Real& Hellp
    )
{
  const Real TwoPi = 2.0 * PI;
  const Real Small = 0.00000001;       // Small value for tolerances

  Real a, b, RtAsc, SqrtP, third, e, f, p, q, d, nu, g, t, c, aTemp, Temp, GST;
  Integer i;

  /* ----------------- Find Longitude value  ---------------------- */
  Temp = sqrt(R.Get(1) * R.Get(1) + R.Get(2) * R.Get(2));
  if (fabs(Temp) < Small)
    RtAsc = Sgn(R.Get(3)) * PI * 0.5;
  else
    RtAsc = atan2(R.Get(2) / Temp, R.Get(1) / Temp);
  GST = GSTime(JD);
  Lon = RtAsc - GST;
  if (fabs(Lon) >= PI)
    if (Lon < 0.0)
      Lon = TwoPi + Lon;
    else
      Lon = Lon - TwoPi;

  a = 1.0;
  b = Sgn(R.Get(3)) * 6356.75160056 / 6378.1363;
  /* --------------- Set up initial latitude value  --------------- */
  aTemp = 1.0 / (a * Temp);
  e = (b * R.Get(3) - a * a + b * b) * aTemp;
  f = (b * R.Get(3) + a * a - b * b) * aTemp;
  third = 1.0 / 3.0;
  p = 4.0 * third * (e * f + 1.0);
  q = 2.0 * (e * e - f * f);
  d = p * p * p + q * q;

  if (d > 0.0)
    nu = Power(sqrt(d) - q, third) - Power(sqrt(d) + q, third);
  else
  {
    SqrtP = sqrt(-p);
    nu = 2.0 * SqrtP * cos(third * acos(q / (p * SqrtP)));
  }
  g = 0.5 * (sqrt(e * e + nu) + e);
  t = sqrt(g * g + (f - nu * g) / (2.0 * g - e)) - g;
  
  Latgd = atan(a * (1.0 - t * t) / (2.0 * b * t));
  Hellp = (Temp - a * t) * cos(Latgd) + (R.Get(3) - b) * sin(Latgd);

  GeocGeod(Latgc, FROM, Latgd);
}

void IJKtoLatLonE
    (
      Vector R, Real JD, Real& Latgc, Real& Latgd, Real& Lon, Real& Hellp
    )
{
  const Real TwoPi  = 2.0 * PI;
  const Real Small  = 0.00000001;      // Small value for tolerances
  const Real EESqrd = 0.006694385000;  // Eccentricity of Earth Sqrd

  Real rsite, RtAsc, OldDelta, DeltaLat, RSqrd,
       OneMinuse2, Decl, SinTemp,  Temp, GST;
  Integer i;

  /* --------------------  Initialize values   -------------------- */
  OneMinuse2 = 1.0 - EESqrd;

  /* ----------------- Find Longitude value  ---------------------- */
  Temp = sqrt(R.Get(1) * R.Get(1) + R.Get(2) * R.Get(2));
  if (fabs(Temp) < Small)
    RtAsc = Sgn(R.Get(3) * PI * 0.5);
  else
    RtAsc = atan2( R.Get(2) / Temp, R.Get(1) / Temp);
  GST = GSTime(JD);
  Lon = RtAsc - GST;

  if (fabs(Lon) >= PI)
    if (Lon < 0.0)
      Lon = TwoPi + Lon;
    else
      Lon = Lon - TwoPi;

  /* --------------- Set up initial latitude value  --------------- */
  Decl     = asin(R.Get(3) / R.Mag());
  Latgc    = Decl;
  DeltaLat = 100.0;
  RSqrd    = R.Mag() * R.Mag();

  /* ----- Iterate to find Geocentric and Geodetic Latitude  ----- */
  i = 1;
  while (1 == 1)
  {
    OldDelta = DeltaLat;
    rsite    = sqrt(OneMinuse2 / (1.0 - EESqrd * sqrt(cos(Latgc))));
    Latgd    = atan(tan(Latgc) / OneMinuse2);
    Temp     = Latgd - Latgc;
    SinTemp  = sin(Temp);
    Hellp    = sqrt(RSqrd - rsite * rsite * SinTemp * SinTemp) - 
                    rsite * cos(Temp);
    DeltaLat = asin(Hellp * SinTemp / R.Mag());
    Latgc    = Decl - DeltaLat;
    i++;
    if (Show == 'Y')
      printf("E loops gc gd %12.6f 16.10f\n", 
              Latgc * 57.29578, Latgd * 57.29578);

    if ((fabs(OldDelta - DeltaLat) < Small) || (i >= 10))
      break;
  }

  if ( i >= 10)
    printf("IJKtoLatLon did NOT converge\n");
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE LIGHT
|
|  This PROCEDURE determines IF a spacecraft is sunlit or in the dark at a
|    particular time.  An obLate Earth and cylindrical shadow is assumed.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - Position vector of sat         ER
|    JD          - Julian Date at desired time    Days from 4713 BC
|    WhichKind   - Spherical or Ellipsoidal Earth 'S', 'E'*default
|
|  OutPuts       :
|    Vis         - Visibility Flag                'YES','NO '
|
|  Locals        :
|    RtAsc       - Suns Right ascension           rad
|    Decl        - Suns Declination               rad
|    RSun        - SUN vector                     AU
|    AUER        - Conversion from AU to ER
|
|  Coupling      :
|    SUN         - Position vector of SUN
|    LNCOM1      - Multiple a vector by a constant
|    SIGHT       - Does Line-of-SIGHT exist beteen vectors
|
|  References    :
|    Vallado       2001, 291-295, Alg 35, Ex 5-6
|
 -----------------------------------------------------------------------------*/
void Light(Vector R, Real JD, char WhichKind, char *Lit)
{
  Vector RSun;
  Real   AUER, RtAsc, Decl;

  AUER = 149597870.0 / 6378.1363;
  
  Sun(JD,RSun, RtAsc, Decl);
  RSun = AUER * RSun;

  /* ------------- Is the satellite in the shadow? ---------------- */
  Sight(RSun, R, WhichKind,  Lit);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE KEPLER
|
|  This PROCEDURE solves Keplers problem for orbit determination and returns a
|    future Geocentric Equatorial (IJK) position and velocity vector.  The
|    solution uses Universal variables.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ro          - IJK Position vector - initial  ER
|    Vo          - IJK Velocity vector - initial  ER / TU
|    DtTU        - Length of time to propagate    TU
|
|  OutPuts       :
|    R           - IJK Position vector            ER
|    V           - IJK Velocity vector            ER / TU
|    Error       - Error flag                     'ok', ...
|
|  Locals        :
|    F           - f expression
|    G           - g expression
|    FDot        - f DOT expression
|    GDot        - g DOT expression
|    XOld        - Old Universal Variable X
|    XOldSqrd    - XOld squared
|    XNew        - New Universal Variable X
|    XNewSqrd    - XNew squared
|    ZNew        - New value of z
|    C2New       - C2(psi) FUNCTION
|    C3New       - C3(psi) FUNCTION
|    DtTU        - change in time                 TU
|    TimeNew     - New time                       TU
|    RDotV       - Result of Ro DOT Vo
|    A           - Semi or axis                   ER
|    Alpha       - Reciprocol  1/a
|    SME         - Specific Mech Energy           ER2 / TU2
|    Period      - Time period for satellite      TU
|    S           - Variable for parabolic case
|    W           - Variable for parabolic case
|    H           - Angular momentum vector
|    Temp        - Temporary EXTENDED value
|    i           - index
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    DOT         - DOT product of two vectors
|    REALMOD     - MOD function for Real variables
|    COT         - Cotangent FUNCTION
|    POWER       - Raise a Number to some POWER
|    SGN         - Sign of a Number +1 or -1
|    FINDC2C3    - Find C2 and C3 functions
|    CROSS       - CROSS product of two vectors
|
|  References    :
|    Vallado       2001, 95-103, Alg 8, Ex 2-4
|
 -----------------------------------------------------------------------------*/
void Kepler(Vector Ro, Vector Vo, Real DtTU, Vector& R, Vector& V, char *Error)
{
  const Real HalfPi = PI / 2.0;
  const Real TwoPi  = 2.0 * PI;
  const Real Small  = 0.00000001;     // Small value for tolerances, 2 less 0's
  const Real Infinite = 999999.9;     // Infinite value

  const Integer NumIter = 35;

  Vector  h;
  Integer Ktr, i;
  Real    F, G, FDot, GDot, Rval, XOld, XOldSqrd, XNew, XNewSqrd, ZNew, p,
          C2New, C3New, DtNew, RDotV, A, Alpha, SME, Period, S, W, Temp;

  /* ---------------------  Initialize values   ------------------- */
  Ktr = 0;
  XOld = 0.0;
  ZNew = 0.0;
  strcpy(Error, "ok");

  if (fabs(DtTU) > Small)
  {
    RDotV = Ro.Dot(Vo);

    /* --------------  Find SME, Alpha, and A  ------------------ */
    SME = (Vo.Mag() * Vo.Mag() * 0.5) - (1.0 / Ro.Mag());
    Alpha = -SME * 2.0;

    if (fabs(SME) > Small)   // circle, ellipse, hyperbola
      A = -1.0 / (2.0 * SME);
    else
      A = Infinite;
    if (fabs(Alpha) < Small)  // Parabola
      Alpha = 0.0;

    /* -------------   Setup initial guess for x  --------------- */
    /* ---------------  Circle and Ellipse ---------------------- */
    if (Alpha >= Small)
    {
      Period = TwoPi * sqrt(Power(fabs(A), 3.0));
      /* --- Next 2 lines needed for 2body multi-rev -- */
      if (fabs(DtTU) > fabs(Period))
        DtTU = Mod(DtTU, Period);
      if (fabs(Alpha - 1.0) > Small)
        XOld = DtTU * Alpha;
      else
        /* ---- 1st guess can't be too close. ie a circle, r=a ---- */
        XOld = DtTU * Alpha * 0.97;
    }
    else
    {
      /* ------------------  Parabola  ---------------- */
      if (fabs(Alpha) < Small)
      {
        h = Ro.Cross(Vo);
        p = h.Mag() * h.Mag();
        S = 0.5 * (HalfPi - atan(3.0 * sqrt(1.0 / (p * p * p)) * DtTU));

        W = atan(Power(tan(S) ,1.0 / 3.0));
        XOld = sqrt(p) * (2.0 * Cot(2.0 * W));
        Alpha = 0.0;
      }
      else
      {
        /* --------------  Hyperbola  --------------- */
        Temp = -2.0*DtTU /
               (A * (RDotV + Sgn(DtTU) * sqrt(-A) * (1.0 - Ro.Mag() * Alpha)));
        XOld = Sgn(DtTU) * sqrt( -A ) * log(Temp);
      }
    }

    Ktr = 1;
    while (1 == 1)
    {
      XOldSqrd = XOld * XOld;
      ZNew     = XOldSqrd * Alpha;

      /* -------------- Find C2 and C3 functions -------------- */
      FindC2C3(ZNew, C2New, C3New);

      /* ------- Use a Newton iteration for New values -------- */
      DtNew = XOldSqrd * XOld * C3New + RDotV * XOldSqrd * C2New +
              Ro.Mag() * XOld * (1.0 - ZNew * C3New);
      Rval  = XOldSqrd * C2New + RDotV * XOld * (1.0 - ZNew * C3New) +
              Ro.Mag() * (1.0 - ZNew * C2New);

      /* ------------- CalcuLate New value for x -------------- */
      XNew = XOld + (DtTU - DtNew) / Rval;

      /* -------------------------------------------------------------
        Check IF the orbit is an ellipse and xNew > 2pi SQRT(a), the step
        size must be changed.  This is accomplished by multiplying Rval
        by 10.0.  NOTE !! 10.0 is arbitrary, but seems to produce good
        results.  The idea is to keep XNew from increasing too rapidily.
      -------------------------------------------------------------- */
/* including this doesn't work IF you don't MOD the DtTU */
/*
      if ((A > 0.0) && (fabs(XNew) > TwoPi * sqrt(A)) && (SME < 0.0))
      {
        dx   = (DtTU - DtNew) / Rval; // *7.0  * 10.0 
        XNew = XOld + Dx / 7.0;       // /(1.0 + Dx) 
        if (FileOut != NULL)
          fprintf(FileOut, "dx orig %11.7f  d10 %11.7f  d7 %11.7f %11.7f\n",
                           dx, dx * 0.1, dx / 7.0, dx / (1.0 + dx));
// Alternate method to test various values of change
printf("chgamt %2d  xn %18.8f rval %18.8f a %18.8f\n", Ktr, xNew, RVal, a);
XNew := XOld + (DtTU - DtNew) / (RVal * 10 chgamt);
      }
*/

      if ((Show == 'Y') || (Show == 'S'))
        if (FileOut != NULL)
        fprintf(FileOut, "%2d %9.6f %9.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
                         i, XOld, DtNew, Rval, XNew, C3New, C2New, ZNew);

      Ktr++;
      XOld = XNew;
      if ((fabs(DtNew - DtTU) < Small) || (Ktr >= NumIter))
        break;
    }

    if (Ktr >= NumIter)
    {
      strcpy(Error, "KNotConv");
      printf("Not converged in %2d iterations\n", NumIter);
      V.Clear();
      R = V;
    }
    else
    {
      /* ---- Find position and velocity vectors at New time ---- */
      XNewSqrd = XNew * XNew;
      F = 1.0 - (XNewSqrd * C2New / Ro.Mag());
      G = DtTU - XNewSqrd * XNew * C3New;
      R = F * Ro + G * Vo;
      GDot = 1.0 - (XNewSqrd * C2New / R.Mag());
      FDot = (XNew / (Ro.Mag() * R.Mag())) * (ZNew * C3New - 1.0);
      V = FDot * Ro + GDot * Vo;
      Temp = F * GDot - FDot * G;
      if (fabs(Temp - 1.0) > 0.00001)
      {
        strcpy(Error, "FandG");
        printf(" Error in f and g in KEPLER %14.10f\n", Temp);
      }
      if ((Show == 'Y') || (Show == 'S'))
        if (FileOut != NULL)
          fprintf(FileOut, "f g fgdot %11.7f %11.7f %11.7f %11.7f %11.7f\n",
                            F, G, FDot, GDot, (F * GDot - G * FDot));
    }
  }
  else
  {
    /* ------------ Set vectors to incoming since 0 time ---------- */
    R = Ro;
    V = Vo;
  }
  
  if (FileOut != NULL)
    fprintf(FileOut, "%12.5f %11.4f 12.5f %6d\n", 
                     ZNew, DtTU * 13.44685108, XOld, Ktr);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE MOON
|
|  This PROCEDURE calcuLates the Geocentric Equatorial (IJK) position vector
|    for the MOON given the Julian Date.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian Date                    days from 4713 BC
|
|  Outputs       :
|    RMoon       - IJK Position vector of MOON    ER
|    RtAsc       - Right Ascension                rad
|    Decl        - Declination                    rad
|
|  Locals        :
|    EclpLong    - Ecliptic Longitude
|    EclpLat     - Eclpitic Latitude
|    HzParal     - Horizontal Parallax
|    l           - Geocentric Direction Cosines
|    m           -             "     "
|    n           -             "     "
|    TTDB        - Julian Centuries of TDB from
|                  Jan 1, 2000 12h
|    Hr          - Hours                          0 .. 24
|    MIN         - Minutes                        0 .. 59
|    SEC         - Seconds                        0.0 .. 59.99
|    deg         - Degrees
|
|  Coupling      :
|    REALMOD       MOD FUNCTION for REAL variables
|    ARCSIN        Arc Sine FUNCTION
|    ATAN2         Arc Tangent formula which resolves quadrants
|
|  References    :
|    Vallado       2001, 272-275, Alg 31, Ex 5-3
|
 -----------------------------------------------------------------------------*/
void Moon(Real JD, Vector& RMoon, Real& RtAsc, Real& Decl)
{
  const Real TwoPi   = 2.0 * PI;
  const Real Deg2Rad = PI / 180.0;

  Real SEC, TTDB, l, m, n, Temp, Obliquity, EclpLong, EclpLat, HzParal, RMoon4;
  Integer Hr, MIN, Deg;

  /* --------------------  Initialize values   -------------------- */
  TTDB = (JD - 2451545.0) / 36525.0;
  
  EclpLong = 218.32 + 481267.883 * TTDB + 
               6.29 * sin((134.9  + 477198.85 * TTDB) * Deg2Rad) - 
               1.27 * sin((259.2  - 413335.38 * TTDB) * Deg2Rad) + 
               0.66 * sin((235.7  + 890534.23 * TTDB) * Deg2Rad) + 
               0.21 * sin((269.9  + 954397.70 * TTDB) * Deg2Rad) - 
               0.19 * sin((357.5  +  35999.05 * TTDB) * Deg2Rad) - 
               0.11 * sin((186.6  + 966404.05 * TTDB) * Deg2Rad);     // Deg

  EclpLat  =   5.13 * sin(( 93.3 + 483202.03 * TTDB) * Deg2Rad) + 
               0.28 * sin((228.2 + 960400.87 * TTDB) * Deg2Rad) - 
               0.28 * sin((318.3 +   6003.18 * TTDB) * Deg2Rad) - 
               0.17 * sin((217.6 - 407332.20 * TTDB) * Deg2Rad);     // Deg

  HzParal  =   0.9508 + 0.0518 * cos((134.9 + 477198.85 * TTDB) * Deg2Rad) + 
               0.0095 * cos((259.2 - 413335.38 * TTDB) * Deg2Rad) + 
               0.0078 * cos((235.7 + 890534.23 * TTDB) * Deg2Rad) + 
               0.0028 * cos((269.9 + 954397.70 * TTDB) * Deg2Rad);   // Deg

  EclpLong = Mod(EclpLong * Deg2Rad, TwoPi);
  EclpLat  = Mod(EclpLat * Deg2Rad, TwoPi);
  HzParal  = Mod(HzParal * Deg2Rad, TwoPi);

  Obliquity = 23.439291 - 0.0130042 * TTDB; // deg
  Obliquity = Obliquity * Deg2Rad;

  /* ------------- Find the geocentric direction cosines ---------- */
  l = cos(EclpLat)   * cos(EclpLong);
  m = cos(Obliquity) * cos(EclpLat) * sin(EclpLong) - 
      sin(Obliquity) * sin(EclpLat);
  n = sin(Obliquity) * cos(EclpLat) * sin(EclpLong) +
      cos(Obliquity) * sin(EclpLat);

  /* -------------- CalcuLate MOON position vector ---------------- */
  RMoon4 = 1.0 / sin(HzParal);
  RMoon.Set(RMoon4 * l, 1);
  RMoon.Set(RMoon4 * m, 2);
  RMoon.Set(RMoon4 * n, 3);

  /* --------------- Find Rt Ascension and Declination ------------ */
  RtAsc = atan2(m, l);
  Decl  = asin(n);

  if (Show == 'Y')
  {
    printf("MOON Test Case \n");
    printf("TTDB     %18.12f\n", TTDB);
    DMS_Rad(Deg, MIN, SEC, FROM, EclpLong);
    printf("Ecl Lon  %15.9f %3d %3d %7.3f", EclpLong / Deg2Rad, Deg, MIN, SEC);
    DMS_Rad(Deg, MIN, SEC, FROM, EclpLat);
    printf("Ecl Lat  %15.9f  %3d %3d %7.3f", EclpLat / Deg2Rad, Deg, MIN, SEC);
    printf("Parallax %15.9f\n", HzParal / Deg2Rad);
    DMS_Rad(Deg, MIN, SEC, FROM, Obliquity);
    printf("Obliquit %15.9f %3d %3d %7.3f", 
                     Obliquity/Deg2Rad, Deg2Rad, Deg, MIN, SEC);
    HMS_Rad(Hr, MIN, SEC, FROM, RtAsc);
    printf("RtAsc    %15.9f %3d %3d %7.3f HMS",
                     RtAsc / Deg2Rad, Deg2Rad, Deg, MIN, SEC);
    DMS_Rad(Deg, MIN, SEC, FROM, Obliquity);
    printf("decl     %15.9f %3d %3d %7.3f", 
                     Decl / Deg2Rad, Deg2Rad, Deg, MIN, SEC);
    printf("%14.7f %14.7f %14.7f %14.7f ER", 
            RMoon.Get(1), RMoon.Get(2), RMoon.Get(3), RMoon.Mag());
    Temp = 6378.1363;
    printf("%18.6f %18.6f %18.6f %18.6f ER", 
            RMoon.Get(1) * Temp, RMoon.Get(2) * Temp, 
            RMoon.Get(3) * Temp, RMoon.Mag()  * Temp);
  }
}

void MoonIll(Real MoonEl, Real f, Real& MoonIll)
{
  Real x, g, l0, l1, l2, l3, HzParal;

  x = MoonEl / 90.0;  // fractional part using deg
  g = 1.0;

  if (MoonEl >= 20.0)
  {
    l0 = -1.95;
    l1 =  4.06;
    l2 = -4.24;
    l3 =  1.56;
  }
  else
  {
    if ((MoonEl >= 5.0) && (MoonEl < 20.0))
    {
      l0 =  -2.58;
      l1 =  12.58;
      l2 = -42.58;
      l3 =  59.06;
    }
    else
    {
      if ((MoonEl > -0.8) && (MoonEl < 5.0))
      {
        l0 =   -2.79;
        l1 =   24.27;
        l2 = -252.95;
        l3 = 1321.29;
      }
      else
      {
        l0 = 0.0;
        l1 = 0.0;
        l2 = 0.0;
        l3 = 0.0;
        f  = 0.0;
        g  = 0.0;
      }
    }
  }
  
  l1 = l0 + l1*x + l2*x*x + l3*x*x*x; // deg finds ill of full moon with MoonEl
  l2 = (-0.00868*f - 2.2E-9*f*f*f*f); // deg finds correct for phase

/*
  HzParal = 0.9508 + 0.0518*cos( (134.9+477198.85*TTDB)*Deg2Rad ) + 
            0.0095*cos( (259.2-413335.38*TTDB)*Deg2Rad ) + 
            0.0078*cos( (235.7+890534.23*TTDB)*Deg2Rad ) + 
            0.0028*COS( (269.9+954397.70*TTDB)*Deg2Rad );   // Deg }
  HzParal = Mod( HzParal*Deg2Rad, TwoPi );
  l3:= (2.0* POWER(10.0,(HzParal*rad / 0.951))*g ); // use g to eliminate neg el passes
*/

  MoonIll = Power(10.0,(l1 + l2 /* + l3 */));
  if ((MoonIll < -1.0E+36) || (MoonIll > 0.999))
    MoonIll = 0.0;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE NEWTONE
|
|  This PROCEDURE solves Keplers equation when the Eccentric, paraboic, or
|    Hyperbolic anomalies are known. The Mean anomaly and true anomaly are
|    calcuLated.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ecc         - Eccentricity                   0.0 to
|    E0          - Eccentric Anomaly              -2Pi to 2Pi rad
|
|  Outputs       :
|    M           - Mean Anomaly                   0.0 to 2Pi rad
|    Nu          - True Anomaly                   0.0 to 2Pi rad
|
|  Locals        :
|    Sinv        - Sine of Nu
|    Cosv        - Cosine of Nu
|
|  Coupling      :
|    ATAN2       - Arc tangent FUNCTION which also resloves quadrants
|    SINH        - Hyperbolic Sine
|    COSH        - Hyperbolic Cosine
|
|  References    :
|    Vallado       2001, 85, Alg 6
|
 -----------------------------------------------------------------------------*/
void NewtonE(Real Ecc, Real E0, Real& M, Real& Nu)
{
  const Real Small = 0.00000001;       // Small value for tolerances

  Real Sinv, Cosv;

  /* -------------------------- Circular -------------------------- */
  if (fabs(Ecc) < Small)
  {
    M  = E0;
    Nu = E0;
  }
  else
  {
    /* ------------------------ Elliptical ---------------------- */
    if (Ecc < 0.999)
    {
      M    = E0 - Ecc * sin(E0);
      Sinv = (sqrt(1.0 - Ecc * Ecc) * sin(E0)) / (1.0-Ecc * cos(E0));
      Cosv = (cos(E0) - Ecc) / (1.0 - Ecc * cos(E0));
      Nu   = atan2(Sinv, Cosv);
    }
    else
    {
      /* ----------------------- Hyperbolic  ------------------ */
      if (Ecc > 1.0001)
      {
        M    = Ecc * sinh(E0) - E0;
        Sinv = (sqrt(Ecc * Ecc - 1.0) * sinh(E0)) / (1.0 - Ecc * cosh(E0));
        Cosv = (cosh(E0) - Ecc) / (1.0 - Ecc * cosh(E0));
        Nu   = atan2( Sinv, Cosv);
      }
      else
      {
        /* --------------------- Parabolic -------------------- */
        M  = E0 + (1.0 / 3.0) * E0 * E0 * E0;
        Nu = 2.0 * atan(E0);
      }
    }
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE NEWTONM
|
|  This PROCEDURE performs the Newton Rhapson iteration to find the
|    Eccentric Anomaly given the Mean anomaly.  The True Anomaly is also
|    calcuLated.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ecc         - Eccentricity                   0.0 to
|    M           - Mean Anomaly                   -2Pi to 2Pi rad
|
|  Outputs       :
|    E0          - Eccentric Anomaly              0.0 to 2Pi rad
|    Nu          - True Anomaly                   0.0 to 2Pi rad
|
|  Locals        :
|    E1          - Eccentric Anomaly, next value  rad
|    Sinv        - Sine of Nu
|    Cosv        - Cosine of Nu
|    Ktr         - Index
|    R1r         - CUBIC roots - 1 to 3
|    R1i         - imaginary component
|    R2r         -
|    R2i         -
|    R3r         -
|    R3i         -
|    S           - Variables for parabolic solution
|    W           - Variables for parabolic solution
|
|  Coupling      :
|    ATAN2       - Arc tangent FUNCTION which also resloves quadrants
|    CUBIC       - Solves a CUBIC polynomial
|    POWER       - Raises a base Number to an arbitrary POWER
|    SINH        - Hyperbolic Sine
|    COSH        - Hyperbolic Cosine
|    SGN         - Returns the sign of an argument
|
|  References    :
|    Vallado       2001, 72-75, Alg 2, Ex 2-1
|
 -----------------------------------------------------------------------------*/
void NewtonM(Real Ecc, Real M, Real& E0, Real& Nu)
{
  const Real NumIter = 50;
  const Real Small   = 0.00000001;       // Small value for tolerances

  Real s, w, E1, Sinv, Cosv, R1r, R1i, R2r, R2i, R3r, R3i;
  Integer Ktr;

  /* -------------------------- Hyperbolic  ----------------------- */
  if ((Ecc - 1.0) > Small)
  {
    /* ------------  Initial Guess ------------- */
    if (Ecc < 1.6)
      if (((M < 0.0) && (M > -PI)) || (M > PI))
        E0 = M - Ecc;
      else
        E0 = M + Ecc;
    else
      if ((Ecc < 3.6) && (fabs(M) > PI)) // just edges)
        E0 = M - Sgn(M) * Ecc;
      else
        E0 = M / (Ecc - 1.0); // best over 1.8 in middle
    Ktr = 1;
    E1 = E0 + ((M - Ecc * sinh(E0) + E0) / (Ecc * cosh(E0) - 1.0));
    while ((fabs(E1 - E0) > Small) && (Ktr <= NumIter))
    {
      E0 = E1;
      E1 = E0 + ((M - Ecc * sinh(E0) + E0) / (Ecc * cosh(E0) - 1.0));
      Ktr++;
    }
    /* ---------  Find True Anomaly  ----------- */
    Sinv = -(sqrt(Ecc * Ecc - 1.0) * sinh(E1)) / (1.0 - Ecc * cosh(E1));
    Cosv = (cosh(E1) - Ecc) / (1.0 - Ecc * cosh(E1));
    Nu   = atan2(Sinv, Cosv);
  }
  else
  {
    /* ---------------------- Parabolic ------------------------- */
    if (fabs(Ecc - 1.0) < Small)
    {
      Cubic(1.0 / 3.0, 0.0, 1.0, -M, R1r, R1i, R2r, R2i, R3r, R3i);
      E0 = R1r;
      if (FileOut != NULL)
        fprintf(FileOut, "roots %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n",
                          R1r, R1i, R2r, R2i, R3r, R3i);
/*
      S  = 0.5 * (HalfPi - atan(1.5 * M));
      W  = atan(Power(tan(S), 1.0 / 3.0));
      E0 = 2.0 * Cot(2.0* W );
*/
      Ktr = 1;
      Nu  = 2.0 * atan(E0);
    }
    else
    {
      /* --------------------- Elliptical --------------------- */
      if (Ecc > Small)
      {
        /* ------------  Initial Guess ------------- */
        if (((M < 0.0) && (M > -PI)) || (M > PI))
          E0 = M - Ecc;
        else
          E0 = M + Ecc;
        Ktr = 1;
        E1  = E0 + (M - E0 + Ecc * sin(E0)) / (1.0 - Ecc * cos(E0));
        while ((fabs(E1 - E0) > Small) && (Ktr <= NumIter))
        {
          Ktr++;
          E0 = E1;
          E1 = E0 + (M - E0 + Ecc * sin(E0)) / (1.0 - Ecc * cos(E0));
        }
        /* ---------  Find True Anomaly  ----------- */
        Sinv = (sqrt(1.0 - Ecc * Ecc) * sin(E1)) / (1.0-Ecc * cos(E1));
        Cosv = (cos(E1) - Ecc) / (1.0 - Ecc * cos(E1));
        Nu   = atan2( Sinv, Cosv);
      }
      else
      {
        /* --------------------- Circular --------------------- */
        Ktr = 0;
        Nu  = M;
        E0  = M;
      }
    }
  }
  if (Ktr > NumIter)
    printf("NewtonRhapson not converged in %3d Iterations\n", NumIter);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE NEWTONNU
|
|  This PROCEDURE solvoes Keplers equation when the true anomaly is known.
|    The Mean and Eccentric, parabolic, or hyperbolic anomalies are also found.
|    The parabolic limit at 168ø is arbitrary.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ecc         - Eccentricity                   0.0 to
|    Nu          - True Anomaly                   -2Pi to 2Pi rad
|
|  Outputs       :
|    E0          - Eccentric Anomaly              0.0 to 2Pi rad
|    M           - Mean Anomaly                   0.0 to 2Pi rad
|
|  Locals        :
|    SinE        - Sine of E
|    CosE        - Cosine of E
|    Ktr         - Index
|
|  Coupling      :
|    ATAN2       - Arc tangent FUNCTION which also resloves quadrants
|    ARCSINH     - Arc hyperbolic sine
|    SINH        - Hyperbolic Sine
|
|  References    :
|    Vallado       2001, 85, Alg 5
|
 -----------------------------------------------------------------------------*/
void NewtonU(Real Ecc, Real Nu, Real& E0, Real& M)
{
  Real SinE, CosE;

  E0 = 999999.9;
  M  = 999999.9;
  /* ---------------------------- Circular ------------------------ */
  if (fabs(Ecc) < 0.000001)
  {
    M  = Nu;
    E0 = Nu;
  }
  else
  {
    /* ----------------------- Elliptical ----------------------- */
    if (Ecc < 0.999)
    {
      SinE = (sqrt(1.0 - Ecc * Ecc) * sin(Nu)) / (1.0 + Ecc * cos(Nu));
      CosE = (Ecc + cos(Nu)) / (1.0 + Ecc * cos(Nu));
      E0   = atan2(SinE, CosE);
      M    = E0 - Ecc * sin(E0);
    }
    else
    {
      /* --------------------- Hyperbolic  -------------------- */
      if (Ecc > 1.0001)
      {
        if (((Ecc > 1.0) && (fabs(Nu) + 0.00001 < PI - acos(1.0 / Ecc))))
        {
          SinE = (sqrt(Ecc * Ecc-1.0) * sin(Nu)) / (1.0 + Ecc * cos(Nu));
          E0   = asinh(SinE);
          M    = Ecc * sinh(E0) - E0;
        }
      }
      else
      {
        /* ------------------ Parabolic --------------------- */
        if (fabs(Nu) < 168.0 / 57.29578)   // TAN undef at 180
        {
          E0 = tan(Nu * 0.5);
          M  = E0 + (E0 * E0 * E0) / 3.0;
        }
      }
    }
  }

  if (Ecc < 1.0)   // Since B and H are areas
  {
    M = Mod(M, 2.0 * PI);
    if (M < 0.0)
      M = M + 2.0 * PI;
    E0 = Mod(E0, 2.0 * PI);
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE PATH
|
|  This PROCEDURE determines the END position for a given range and azimuth
|    from a given point.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    LLat        - Start Geocentric Latitude      -Pi/2 to  Pi/2 rad
|    LLon        - Start Longitude (WEST -)       0.0 to 2Pi rad
|    Range       - Range between points           ER
|    Az          - Azimuth                        0.0 to 2Pi rad
|
|  OutPuts       :
|    TLat        - END Geocentric Latitude        -Pi/2 to  Pi/2 rad
|    TLon        - END Longitude (WEST -)         0.0 to 2Pi rad
|
|  Locals        :
|    SinDeltaN   - Sine of Delta N                rad
|    CosDeltaN   - Cosine of Delta N              rad
|    DeltaN      - ANGLE bteween the two points   rad
|
|  Coupling      :
|    ARCSIN      - Arc sine FUNCTION
|    REALMOD     - MOD FUNCTION for REAL variables
|    ATAN2       - Arc tangent FUNCTION which also resolves quadrants
|
|  References    :
|    Vallado       2001, 774-776, Eq 11-6, Eq 11-7
|
 -----------------------------------------------------------------------------*/
void Path(Real LLat, Real LLon, Real Range, Real Az, Real& TLat, Real& TLon)
{
  const Real TwoPi = 2.0 * PI;
  const Real Small = 0.00000001;       // Small value for tolerances

  Real  SinDN, CosDN, DeltaN;

  Az = Mod(Az, TwoPi);
  if (LLon < 0.0)
    LLon = TwoPi + LLon;
  if (Range > TwoPi)
    Range = Mod(Range, TwoPi);

  /* ------------------ Find Geocentric Latitude  ----------------- */
  TLat = asin(sin(LLat) * cos(Range) + cos(LLat) * sin(Range) * cos(Az));

  /* ----- Find Delta N, the ANGLE between the points ------------- */
  if ((fabs(cos(TLat)) > Small) && (fabs(cos(LLat)) > Small))
  {
    SinDN = sin(Az) * sin(Range) / cos(TLat);
    CosDN = (cos(Range) - sin(TLat) * sin(LLat)) / (cos(TLat) * cos(LLat));
    DeltaN = atan2(SinDN, CosDN);
  }
  else
  {
    /* ------- Case where launch is within 3nm of a Pole -------- */
    if (fabs(cos(LLat)) <= Small)
      if ((Range > PI) && (Range < TwoPi))
        DeltaN = Az + PI;
      else
        DeltaN = Az;
    /* ------ Case where END point is within 3nm of a pole ------ */
    if (fabs(cos(TLat)) <= Small)
      DeltaN = 0.0;
  }

  TLon = LLon + DeltaN;
  if (fabs(TLon) > TwoPi)
    TLon = Mod(TLon, TwoPi);
  if (TLon < 0.0)
    TLon = TwoPi + TLon;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RANDV
|
|  This PROCEDURE finds the position and velocity vectors in Geocentric
|    Equatorial (IJK) system given the classical orbit elements. Notice that
|    p is used for calcuLations and that semi-major axis a, is not.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    P           - SemiLatus rectum               ER
|    Ecc         - Eccentricity
|    Incl        - inclination                    0.0 to Pi rad
|    Omaga       - Longitude of Ascending Node    0.0 to 2Pi rad
|    Argp        - Argument of Perigee            0.0 to 2Pi rad
|    Nu          - True anomaly                   0.0 to 2Pi rad
|    ArgLat      - Argument of Latitude      (CI) 0.0 to 2Pi rad
|    LamTrue     - True Longitude            (CE) 0.0 to 2Pi rad
|    LonPer      - Longitude of Periapsis    (EE) 0.0 to 2Pi rad
|
|  Outputs       :
|    R           - IJK Position vector            ER
|    V           - IJK Velocity vector            ER / TU
|
|  Locals        :
|    Temp        - Temporary EXTENDED value
|    Rpqw        - PQW Position vector            ER
|    Vpqw        - PQW Velocity vector            ER / TU
|    SinNu       - Sine of Nu
|    CosNu       - Cosine of Nu
|    TempVec     - PQW Velocity vector
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ROT3        - Rotation about the 3rd axis
|    ROT1        - Rotation about the 1st axis
|
|  References    :
|    Vallado       2001, 122-126, Alg 10, Ex 2-5
|
 -----------------------------------------------------------------------------*/
void RandV
    (
      Real P, Real Ecc, Real Incl, Real Omega, Real Argp, 
      Real Nu, Real ArgLat, Real TrueLon, Real LonPer,
      Vector& R, Vector& V
    )
{
  const Real Small = 0.00000001;       // Small value for tolerances

  Vector Rpqw(3), Vpqw(3), TempVec(3);
  Real   Temp, SinNu, CosNu;

  /* --------------------------------------------------------------
  |   Determine what type of orbit is involved and set up the
  |   set up angles for the special cases.
   --------------------------------------------------------------*/
  if (Ecc < Small)
    /* ------------  Circular Equatorial  ------------ */
    if ((Incl < Small) || (fabs(Incl - PI) < Small))
    {
      Argp  = 0.0;
      Omega = 0.0;
      Nu    = TrueLon;
    }
    else
    {
      /* ----------  Circular Inclined  -------------- */
      Argp  = 0.0;
      Nu    = ArgLat;
    }
  else
    /* ----------  Elliptical Equatorial  ------------ */
    if ((Incl < Small) || (fabs(Incl - PI) < Small))
    {
      Argp  = LonPer;
      Omega = 0.0;
    }

  /* -----------  Form PQW position and velocity vectors ---------- */
  CosNu = cos(Nu);
  SinNu = sin(Nu);
  Temp  = P / (1.0 + Ecc * CosNu);
  Rpqw.Set(Temp * CosNu, 1);
  Rpqw.Set(Temp * SinNu, 2);
  Rpqw.Set(0.0, 3);
  if (fabs(P) < 0.0001)   // IF input is 0.0
    P = 0.0001;
  Vpqw.Set(-SinNu / sqrt(P), 1);
  Vpqw.Set((Ecc + CosNu) / sqrt(P), 2);
  Vpqw.Set(0.0, 3);

  /* ---------------  Perform transformation to IJK  -------------- */
  TempVec = Rpqw.Rot3(-Argp);
  TempVec = TempVec.Rot1(-Incl);
  V       = TempVec.Rot3(-Omega);

  TempVec = Vpqw.Rot3(-Argp);
  TempVec = TempVec.Rot1(-Incl);
  V       = TempVec.Rot3(-Omega);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RNGAZ
|
|  This PROCEDURE calcuLates the Range and Azimuth between two specified
|    ground points on a spherical Earth.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    LLat        - Start Geocentric Latitude      -Pi/2 to  Pi/2 rad
|    LLon        - Start Longitude (WEST -)       0.0 to 2Pi rad
|    TLat        - END Geocentric Latitude        -Pi/2 to  Pi/2 rad
|    TLon        - END Longitude (WEST -)         0.0 to 2Pi rad
|    Tof         - Time of Flight IF ICBM, or 0.0 TU
|
|  OutPuts       :
|    Range       - Range between points           ER
|    Az          - Azimuth                        0.0 to 2Pi rad
|
|  Locals        :
|    None.
|
|  Coupling      :
|    ARCCOS      - Arc Cosine FUNCTION
|
|  References    :
|    Vallado       2001, 774-775, Eq 11-3, Eq 11-4, Eq 11-5
|
 -----------------------------------------------------------------------------*/
void RngAz
    (
      Real LLat, Real LLon, Real TLat, Real TLon, Real Tof, 
      Real& Range, Real& Az
    )
{
  const Real TwoPi = 2.0 * PI;
  const Real Small = 0.00000001;          // Small value for tolerances

  Real OmegaEarth = 0.05883359221938136;  // Earth Rot rad/TU

  Range = acos(sin(LLat) * sin(TLat) +
               cos(LLat) * cos(TLat) * cos(TLon - LLon + OmegaEarth * Tof));

  /* ------- Check IF the Range is 0 or half the earth  ----------- */
  if (fabs(sin(Range) * cos(LLat)) < Small)
    if (fabs(Range - PI) < Small)
      Az = PI;
    else
      Az = 0.0;
  else
    Az = acos((sin(TLat) - cos(Range) * sin(LLat)) / (sin(Range) * cos(LLat)));

  /* ------- Check IF the Azimuth is grt than Pi ( 180deg ) ------- */
  if (sin(TLon - LLon + OmegaEarth * Tof) < 0.0)
    Az = TwoPi - Az;
}

/*------------------------------------------------------------------------------|
|                           PROCEDURE SIGHT
|
|  This PROCEDURE takes the position vectors of two satellites and determines
|    IF there is line-of-SIGHT between the two satellites.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R1          - Position vector of the 1st sat ER
|    R2          - Position vector of the 2nd sat ER
|    WhichKind   - Spherical or Ellipsoidal Earth 'S', 'E'*default
|
|  Outputs       :
|    LOS         - Line of SIGHT                  'YES','NO '
|
|  Locals        :
|    TR1         - Scaled R1 vector               ER
|    TR2         - Scaled R2 vector               ER
|    ADotB       - DOT product of a DOT b
|    TMin        - Minimum value of t from a to b
|    DistSqrd    - MIN Distance squared to Earth  ER
|    ASqrd       - Magnitude of A squared
|    BSqrd       - Magnitude of B squared
|
|  Coupling:
|    DOT         - DOT product of two vectors
|
|  References    :
|    Vallado       2001, 291-295, Alg 35, Ex 5-3
|
 -----------------------------------------------------------------------------*/
void Sight(Vector R1, Vector R2, char WhichKind, char *LOS)
{
  const Real EESqrd = 0.006694385000;   // Eccentricity of Earth Sqrd

  Vector TR1(3), TR2(3);
  Real   ADotB, TMin, DistSqrd, ASqrd, BSqrd, Temp;

  TR1 = R1;
  TR2 = R2;
  /* --------- Scale z component -------- */
  if (WhichKind == 'E')
    Temp = 1.0 / sqrt(1.0 - EESqrd);
  else
    Temp = 1.0;
  TR1.Set(TR1.Get(3) * Temp, 3);
  TR2.Set(TR2.Get(3) * Temp, 3);

  BSqrd = TR2.Mag() * TR2.Mag();
  ASqrd = TR1.Mag() * TR1.Mag();
  ADotB = TR1.Dot(TR2);
  /* ----------- Find TMin -------------- */
  DistSqrd = 0.0;
  if (fabs(ASqrd + BSqrd - 2.0 * ADotB) < 0.0001)
    TMin = 0.0;
  else
    TMin = (ASqrd - ADotB) / (ASqrd + BSqrd - 2.0 * ADotB);
  /* ------------ Check LOS ------------- */
  if ((TMin < 0.0) || (TMin > 1.0))
    strcpy(LOS, "YES");
  else
  {
    DistSqrd = (1.0 - TMin) * ASqrd + ADotB * TMin;
    if (DistSqrd > 1.0)
      strcpy(LOS, "YES");
    else
      strcpy(LOS, "NO");
  }

  if (Show == 'Y')
    printf("Denom %11.7f TMin %11.7f DistS %11.7f\n", 
           (ASqrd + BSqrd - 2.0*ADotB), TMin, DistSqrd);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE SUN
|
|  This PROCEDURE calcuLates the Geocentric Equatorial position vector for
|    the SUN given the Julian Date.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian Date                    days from 4713 BC
|
|  Outputs       :
|    RSun        - IJK Position vector of the SUN AU
|    RtAsc       - Right Ascension                rad
|    Decl        - Declination                    rad
|
|  Locals        :
|    MeanLong    - Mean Longitude
|    MeanAnomaly - Mean anomaly
|    EclpLong    - Ecliptic Longitude
|    Obliquity   - Mean Obliquity of the Ecliptic
|    TUT1        - Julian Centuries of UT1 from
|                  Jan 1, 2000 12h
|    TTDB        - Julian Centuries of TDB from
|                  Jan 1, 2000 12h
|    Hr          - Hours                          0 .. 24
|    MIN         - Minutes                        0 .. 59
|    SEC         - Seconds                        0.0 .. 59.99
|    Temp        - Temporary variable
|    deg         - Degrees
|
|  Coupling      :
|    REALMOD     - MOD FUNCTION for REAL variables
|    ARCSIN      - Arc sine FUNCTION
|
|  References    :
|    Vallado       2001, 263-267, Alg 29, Ex 5-1
|
 -----------------------------------------------------------------------------*/
void Sun(Real JD, Vector& RSun, Real& RtAsc, Real& Decl)
{
  const Real TwoPi = 2.0 * PI;
  const Real Deg2Rad = PI / 180.0;

  Real RSun4, SEC, temp, MeanLong, MeanAnomaly, EclpLong, Obliquity, TUT1, TTDB;
  Integer hr, MIN, deg;

  /* --------------------  Initialize values   -------------------- */
  TUT1 = (JD - 2451545.0) / 36525.0;

  MeanLong = 280.4606184 + 36000.77005361 * TUT1;
  MeanLong = Mod(MeanLong, 360.0);  // deg

  TTDB        = TUT1;
  MeanAnomaly = 357.5277233 + 35999.05034 * TTDB;
  MeanAnomaly = Mod(MeanAnomaly * Deg2Rad, TwoPi);  // rad
  if (MeanAnomaly < 0.0)
    MeanAnomaly = TwoPi + MeanAnomaly;

  EclpLong = MeanLong + 1.914666471 * sin(MeanAnomaly) + 
             0.019994643 * sin(2.0 * MeanAnomaly); // deg

  Obliquity = 23.439291 - 0.0130042 * TTDB;  // deg

  MeanLong = MeanLong * Deg2Rad;
  if (MeanLong < 0.0)
    MeanLong = TwoPi + MeanLong;
  EclpLong   = EclpLong * Deg2Rad;
  Obliquity  = Obliquity * Deg2Rad;

  /* --------- Find magnitude of SUN vector, THEN components ------ */
  RSun4 = 1.000140612 - 0.016708617 * cos(MeanAnomaly) - 
          0.000139589 * cos(2.0 * MeanAnomaly);   // in AU's
  RSun.Set(RSun4 * cos(EclpLong), 1);
  RSun.Set(RSun4 * cos(Obliquity) * sin(EclpLong), 2);
  RSun.Set(RSun4 * sin(Obliquity) * sin(EclpLong), 3);

  RtAsc = atan(cos(Obliquity) * tan(EclpLong));
  /* ---- Check that RtAsc is in the same quadrant as EclpLong ---- */
  if (EclpLong < 0.0)
    EclpLong = EclpLong + TwoPi;   // make sure it's in 0 to 2pi range
  if (fabs(EclpLong - RtAsc) > PI * 0.5)
    RtAsc = RtAsc + 0.5 * PI * Round((EclpLong - RtAsc) / (0.5 * PI));
  Decl = asin(sin(Obliquity) * sin(EclpLong));

  if (Show == 'Y')
  {
    printf("SUN Test Case\n");
    printf("TUT1     %18.12f\n", TUT1);
    printf("Mean Lon %15.9f\n", MeanLong / Deg2Rad);
    printf("M        %15.9f\n", MeanAnomaly / Deg2Rad);
    DMS_Rad(deg, MIN, SEC, FROM, Obliquity);
    printf("Obliquit %15.9f %3d %3d %7.3f\n", Obliquity/Deg2Rad, deg, MIN, SEC);
    HMS_Rad(hr, MIN, SEC, FROM, RtAsc);
    printf("RtAsc    %15.9f %3d %3d %7.3f HMS\n", RtAsc/Deg2Rad, deg, MIN, SEC);
    DMS_Rad(deg, MIN, SEC, FROM, Decl);
    printf("decl     %15.9f %3d %3d %7.3f\n", Decl / Deg2Rad, deg, MIN, SEC);
    printf("%14.7f %14.7f %14.7f %18.7f AU\n", 
            RSun.Get(1), RSun.Get(2), RSun.Get(3), RSun.Mag());
    temp = 149597870.0;
    printf("%18.6f %18.6f %18.6f km\n", 
            RSun.Get(1) * temp , RSun.Get(2) * temp, RSun.Get(3) * temp);
  }
}

void SunIll(Real JD, Real Lat, Real Lon, Real& SunIll, Real& SunAz, Real& SunEl)
{
  const Real AUER  = 149597870.0 / 6378.135;
  const Real TwoPi =          6.28318530717959;
  const Real Rad   = 180.0 / PI;
  const Real Deg2Rad = PI / 180.0;

  Vector RSun(3), RV(3), RhoSat(3);
  Real LST, GST, x, g, LHA, Sinv, Cosv, l0, l1, l2, l3, sRtAsc, sDecl;

  Sun(JD, RSun, sRtAsc, sDecl);  // AU's needed for Sun ill

  LSTime(Lon, JD,  LST, GST);

  LHA = LST - sRtAsc;

  SunEl = asin(sin(sDecl) * sin(Lat) + cos(sDecl) * cos(Lat) * cos(LHA));

  Sinv  = -sin(LHA)   * cos(sDecl) * cos(Lat) / (cos(SunEl) * cos(Lat));
  Cosv  = (sin(sDecl) - sin(SunEl) * sin(Lat)) / (cos(SunEl) * cos(Lat));
  SunAz = atan2(Sinv, Cosv);

  SunEl = SunEl * Rad;  // deg

  if (SunEl > -18.01)    // only do when Sun above astron horizon
  {
    x = SunEl / 90.0;  // fractional part using deg
    g = 1.0;

    if (SunEl >= 20)
    {
      l0 =  3.74;
      l1 = 3.97;
      l2 = -4.07;
      l3 = 1.47;
    }
    else
    {
      if ((SunEl >= -0.8) && (SunEl < 5.0))
      {
        l0 =    2.88;
        l1 =   22.26;
        l2 = -207.64;
        l3 = 1034.30;
      }
      else
      {
        if ((SunEl >= -5.0) && (SunEl < -0.8))
        {
          l0 =    2.88;
          l1 =   21.81;
          l2 = -258.11;
          l3 = -858.36;
        }
        else
        {
          if ((SunEl >= -12.0) && (SunEl < -5.0))
          {
            l0 =    2.70;
            l1 =   12.17;
            l2 = -431.69;
            l3 =-1899.83;
          }
          else
         {
           if ((SunEl >= -18.0) && (SunEl < -12.0))
           {
             l0 =   13.84;
             l1 =  262.72;
             l2 = 1447.42;
             l3 = 2797.93;
           }
           else
           {
             l0 = 0.0;
             l1 = 0.0;
             l2 = 0.0;
             l3 = 0.0;
           }
         }
        }
      }
    }

    l1 = l0 + l1 * x + l2 * x * x + 
         l3 * x * x * x; // deg finds illumination of full Sun with SunEl

    SunIll = Power(10.0, l1);
    if ((SunIll < -1.0E+36) || (SunIll > 999.999))
      SunIll = 0.0;
  }
  else
    SunIll = 0.0;
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE SATFOV
|
|  This PROCEDURE finds parameters reLating to a satellite's FOV.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Incl        - Inclination                    rad
|    Az          - Azimuth                        rad
|    SLatgd      - Geodetic Latitude of sat       rad
|    SLon        - Longitude of sat               rad
|    SAlt        - Altitudeof satellite           ER
|    TFOV        - Total field of view            rad
|    EtaCtr      - Ctr where sensor looks         rad
|
|  Outputs       :
|    FovMax      - Maximum field of view          rad
|    TotalRng    -
|    RhoMax      -
|    RhoMin      -
|    TgtLat      -
|    TgtLon      -
|
|  Locals        :
|    r           -
|    etaHopriz   -
|    RhoHoriz    -
|    gamma       -
|    rho         -
|    FovMin      -
|    Lat         -
|    Lon         -
|    MaxLat      -
|    MinLKat     -
|    i           - Index
|
|  Coupling      :
|    ARCSIN      - Arc Sine FUNCTION
|    PATH        - Finds tgt location given initial location, range, and az
|
|  References    :
|    Vallado       2001, 776-781, Eq 11-8 to Eq 11-13, Ex 11-1
|
  ----------------------------------------------------------------------------*/
void SatFOV
    (
      Real Incl, Real Az, Real SLatgd, Real SLon, 
      Real SAlt, Real tFOV, Real EtaCtr,
      Real& FovMax, Real& TotalRng, Real& RhoMax, 
      Real& RhoMin, Real& TgtLat, Real& TgtLon
    )
{
  const Real Rad2Deg = 180.0 / PI;

  Real r, EtaHoriz, RhoHoriz, GammaM, RhoM, Gamma, Rho, 
       FovMin, Lat, Lon, MaxLat, MinLat;

  /* -------- Find satellite parameters and limiting cases -------- */
  r        = 1.0 + SAlt;
  EtaHoriz = asin(1.0 / r);
  RhoHoriz = r* cos(EtaHoriz);

  /* ----------------- Find Ground range ANGLE -------------------- */
  FovMax = tFOV*0.5 + EtaCtr;
  GammaM = PI - asin(r* sin(FovMax));  // must use larger ANGLE
  RhoM   = cos(GammaM) + r * cos(FovMax);
  RhoMax = asin(RhoM * sin(FovMax));

  /* ----- Do minimum, IF the sensor looks off axis ----- */
  if (fabs(EtaCtr) > 0.00001)
  {
    FovMin   = EtaCtr - tFOV * 0.5;
    Gamma    = PI - asin(  r * sin(FovMin)); // use larger
    Rho      = asin(Rho * sin(FovMin));
    TotalRng = RhoMax - RhoMin;
  }
  else
  {
    Gamma    = 0.0;
    Rho      = 0.0;
    FovMin   = 0.0;
    RhoMin   = 0.0;
    TotalRng = 2.0 * RhoMax; // equal sided
  }

  /* --------------- Find location of center of FOV --------------- */
  if (fabs(EtaCtr) > 0.00001)
    Path(SLatgd, SLon , RhoMin + TotalRng * 0.5, Az,  Lat, Lon);
  else
  {
    Lat = SLatgd;
    Lon = SLon;
  }

  /* ------ Loop around the New circle with the sensor range ------ */
  for (UINT i = 0; i <= 72; i++)
  {
    Az = i * 5.0 / Rad2Deg;
    Path(Lat, Lon, TotalRng * 0.5,Az, TgtLat, TgtLon);
    if (i == 0)
    {
      MaxLat = TgtLat;
      printf("MAX %11.7f %14.7f %14.7f\n", 
              Az * Rad2Deg, TgtLat * Rad2Deg, TgtLon * Rad2Deg);
    }
    if (i == 36)
    {
      MinLat = TgtLat;
      printf("MIN %11.7f %14.7f %14.7f\n",
              Az * Rad2Deg, TgtLat * Rad2Deg, TgtLon * Rad2Deg);
    }
    if (i < 5)
      printf("%11.7f %14.7f %14.7f\n",
              Az * Rad2Deg, TgtLat * Rad2Deg, TgtLon * Rad2Deg);
  }

  if (Show == 'Y')
  {
    printf("%11.7f %14.7f %11.7f %14.7f\n", 
           r, EtaHoriz * Rad2Deg, RhoHoriz, RhoHoriz * 6378.1363);
    printf("%11.7f %14.7f %12.f %14.7f\n",
           FovMax * Rad2Deg, GammaM * Rad2Deg, RhoM, RhoM * 6378.1363);
    printf("%11.7f %14.7f %12.7f %14.7f\n",
           FovMin * Rad2Deg, Gamma * Rad2Deg, Rho, Rho * 6378.1363);
    printf("total rng %14.7f %14.7f\n", 
           TotalRng * Rad2Deg, TotalRng * 6378.1363);
  }
}
