/*     ----------------------------------------------------------------      


                               UNIT ASTPERT


    This file contains fundamental Astrodynamic procedures and functions     
    allowing the calcualtions of perturbations. Most of these routines are   
    discussed in Ch 8.                                                       

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

       ----------------------------------------------------------------      */
#include <math.h>
#include <stdio.h>

#include "ast2body.h"
#include "astiod.h"
#include "astpert.h"
#include "asttime.h"

Matrix C(3, 3), S(3, 3);

/*------------------------------------------------------------------------------
|
|                           PROCEDURE ATMOS
|
|  This PROCEDURE finds the atmospheric density at an altitude above an
|    oblate earth given the position vector in the Geocentric Equatorial
|    frame.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - IJK Position vector            ER
|
|  Outputs       :
|    Rho         - Density                        kg/m**3
|
|  Locals        :
|    Hellp       - Height above ellipsoid         ER
|    OldDelta    - Previous value of DeltaLat     rad
|    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
|    SinTemp     - Sine of Temp
|    RhoNom      - Nominal density at particular alt      gm/cm**3
|    NextBaseAlt - Next Base Altitude
|    LastBaseAlt - Last Base Altitude
|    H           - Scale Height                   km
|    i           - index
|    AtmosFile   - File of data for the
|                    exponential atmosphere
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ARCSIN      - Arc sine FUNCTION
|
|  References    :
|    Vallado       2001, 532-534, Ex 8-4
|
 -----------------------------------------------------------------------------*/
void Atmos(Vector R, Real& Rho)
{
  const Real Small   = 0.0000001;
  const Real EE2Sqrd = 0.006694385000;

  UINT i;
  FILE *AtmosFile;
  const char *AtmosFileName = "atmosexp.dat";
  Real Hellp, OldDelta, Latgd, SinTemp, c, Decl, Temp, H, RhoNom, 
       NextBaseAlt, LastBaseAlt;

  /* --------------------  Initialize values   -------------------- */
  if ((AtmosFile = fopen(AtmosFileName, "r")) == NULL)
  {
    printf("Can't open atmosphere density file %s\n", AtmosFileName);
    exit(0);
  }
  Decl  = asin(R.Get(3) / R.Mag());
  Latgd = Decl;

  /* ----- Iterate to find Geocentric and Geodetic Latitude  ------ */
  Temp = sqrt(Power(R.Get(1), 2) + Power(R.Get(2), 2));
  i = 1;
  while (1 == 1)
  {
    OldDelta = Latgd;
    SinTemp  = sin(Latgd);
    c        = 1.0 / (sqrt(1.0 - EE2Sqrd * SinTemp * SinTemp));
    Latgd    = atan((R.Get(3) + c * EE2Sqrd * SinTemp) / Temp);
    i++;
    if ((fabs( OldDelta - Latgd) < Small) || (i >= 10))
      break;
  }
  Hellp = ((Temp / cos(Latgd)) - c) * 6378.1363;

  if (i >= 10)
    printf("IJKtoLatLon did NOT converge\n");

  /* ----------- Determine density based on altitude -------------- */
  fscanf(AtmosFile, "%f", &NextBaseAlt);
  while (1 == 1)
  {
    LastBaseAlt = NextBaseAlt;
    fscanf(AtmosFile, "%f %f\n", &RhoNom, &H);
    fscanf(AtmosFile, "%f", &NextBaseAlt);
    if ((NextBaseAlt >= Hellp) || (feof(AtmosFile)))
      break;
  }

  Rho = RhoNom * exp((LastBaseAlt - Hellp) / H);
  fclose(AtmosFile);
}

/*------------------------------------------------------------------------------
|
|                                PROCEDURE COWELL
|
|  This PROCEDURE uses a fourth order Runge-Kutta integrator on a 6 dimension
|    First Order differential equation.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - Initial position vector        ER
|    V           - Initial velocity vector        ER/TU
|    ITime       - Initial Time (Julian Date)     Days from 4713 BC
|    FTime       - Final Time (Julian Date)       Days from 4713 BC
|    DtDay       - Step size                      Day
|    DerivType   - String of which perts to incl  'Y' and 'N'
|                  Options are in order : J2, J3,
|                  J4, SUN, MOON, Drag, Solarrad
|    BC          - Ballistic Coefficient          kg/m2
|
|  Outputs       :
|    R1          - Final position vector          ER
|    V1          - Final velocity vector          ER/TU
|
|  Locals        :
|    Time        - Current time during the loop   Days from 4713 BC
|    X           - State vector at each time      ER, ER/TU
|
|  Coupling      :
|    RK4         - Runge-Kutta algorithm
|    MAG         - Magnitude of a vector
|    INITMATRIX  - Initialize a matrix and fil with 0.0's
|    ASSIGNVAL   - Assign a value to a matrix
|    GETVAL      - Get a value from a matrix
|    DELMATRIX   - Delete a matrix
|
 -----------------------------------------------------------------------------*/
void Cowell
    (
      Vector R, Vector V, Real ITime , Real FTime, Real DtDay, char *DerivType,
      SINT Order, Real BC, Vector& R1, Vector& V1
    )
{
  const Real RadiusEarthKm = 6378.1363;
  const Real VKmPerSec     =    7.905366149846;

  Matrix X(6, 1), XDot(6, 1);
  Real   Time;

  for (UINT i = 1; i <= 6; i++)
    if (i <= 3)
      X.Set(R.Get(i), i, 1);
    else
      X.Set(V.Get(i-3), i, 1);

  /* ---------- Loop through the time interval desired ------------ */
  Time = ITime;
  if (Show == 'F')
    if (FileOut != NULL)
      fprintf(FileOut, "%8.3f %12.5f %12.5f %12.5f %11.7f %11.7f %11.7f\n",
                        (Time-ITime)*1440.0,
                        X.Get(1, 1) * RadiusEarthKm,
                        X.Get(2, 1) * RadiusEarthKm,
                        X.Get(3, 1) * RadiusEarthKm,
                        X.Get(4, 1) * VKmPerSec,
                        X.Get(5, 1) * VKmPerSec,
                        X.Get(6, 1) * VKmPerSec);

  /* ---- This would be more accurate with ------------------------
               Julian Date,
               Secs from beginning of day
               Decimal seconds
  ---------------------------------------------------------------*/

  while (Time <= FTime)
  {
    if (Time+DtDay > FTime)
    {
      DtDay = FTime - Time;
      FTime = FTime - 1.0;
      printf("Fixing end time dt is now  %11.6f min\n", DtDay*1440.0);
    }

    RK4(Time, DtDay, XDot, DerivType, Order, BC, X);
/*    RKF45( Time,DtDay,XDot,DerivType,Order,BC, X ); */

    Time = Time + DtDay;
    if (Show == 'F')
      if (FileOut != NULL)
        fprintf(FileOut, "%8.3f %12.5f %12.5f %12.5f %11.7f %11.7f %11.7f\n",
                          (Time-ITime)*1440.0,
                          X.Get(1, 1) * RadiusEarthKm,
                          X.Get(2, 1) * RadiusEarthKm,
                          X.Get(3, 1) * RadiusEarthKm,
                          X.Get(4, 1) * VKmPerSec,
                          X.Get(5, 1) * VKmPerSec,
                          X.Get(6, 1) * VKmPerSec);
  }
                        
  /* ------------------ Update the state vector ------------------- */
  for (UINT i = 1; i <= 6; i++)
    if (i <= 3)
      R1.Set(X.Get(i, 1), i);
    else
      V1.Set(X.Get(i, 1), i-3);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE DERIV
|
|  This PROCEDURE calculates the derivative of the two-body state vector for
|    use with the Runge-Kutta algorithm.  Note time is not needed.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    X           - State Vector                   ER, ER/TU
|
|  Outputs       :
|    XDot        - Derivative of State Vector     ER/TU,  ER/TU2
|
|  Locals        :
|    RCubed      - Cube of R
|
|  Coupling      :
|    ASSIGNVAL   - Assign a value to a matrix
|    GETVAL      - Get a value from a matrix
|
 -----------------------------------------------------------------------------*/
void Deriv(Matrix X, Matrix& XDot)
{
  Real R, RCubed;

  R = sqrt(Power(X.Get(1,1), 2) + Power(X.Get(2,1), 2) + Power(X.Get(3,1), 2));
  RCubed = R * R * R;

  /* --------------------  Velocity Terms  ------------------------ */
  XDot.Set(X.Get(4, 1), 1, 1);
  XDot.Set(X.Get(5, 1), 2, 1);
  XDot.Set(X.Get(6, 1), 3, 1);

  /* ------------------  Acceleration Terms  ---------------------- */
  XDot.Set(X.Get(1, 1), 4, 1);
  XDot.Set(X.Get(2, 1), 5, 1);
  XDot.Set(X.Get(3, 1), 6, 1);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE FullGeop
|
|  This PROCEDURE finds the Legendre polynomial value for the gravity field
|    for a given order.
|
|  Algorithm     : Find the answer
|                  Step up recursion
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Order       - Size of gravity field          1..70
|
|  Outputs       :
|    C           - Gravitational Coefficients
|    S           - Gravitational Coefficients
|
|  Locals        :
|    RCubed      - Cube of R
|
|  Coupling      :
|    IJKTOLATLONA- Find sub satellite point
|
|  References    :
|    Vallado       2001,
|
 -----------------------------------------------------------------------------*/
void FullGeop
    (
      Vector R, Vector V, Real ITime, SINT WhichOne, SINT Order, Real BC,
      Matrix C, Matrix S, Vector& APert
    )
{
  Matrix LArr(70, 70);
  Real   OORDelta, Temp, OOr, SumM1, SumM2, SumM3, DistPartR,
         DistPartPhi, DistPartLon, RDelta, Latgc, Latgd, hellp, Lon;

  /* -------------------- find latgc and lon ---------------------- */
  IJKtoLatLonA(R, ITime, Latgc, Latgd, Lon, hellp);

  /* --------------------- Find Legendre polynomials -------------- */
  LegPoly(Latgc, LArr);

  /* ---------- Partial derivatives of disturbing potential ------- */
  OOr         = 1.0 / R.Mag();
  DistPartR   = 0.0;
  DistPartPhi = 0.0;
  DistPartLon = 0.0;
  Temp        = OOr;

  for (UINT L = 3; L <= Order; L++)
  {
    Temp  = Temp * OOr;  // will do the power as each L is indexed
    SumM1 = 0.0;
    SumM2 = 0.0;
    SumM3 = 0.0;

    for (UINT m = 1; m <= L; m++)
    {
      SumM1 = SumM1 + LArr.Get(L, m) * (C.Get(L, m) * cos(m * Lon) +
              S.Get(L, m) * sin(m * Lon));
      SumM2 = SumM2 + (C.Get(L, m) * cos(m * Lon) + 
              S.Get(L, m) * sin(m * Lon)) *
              (LArr.Get(L, m + 1) - m * LArr.Get(L, m) * tan(Latgc));
      SumM3 = SumM3 + m * LArr.Get(L, m) * (S.Get(L, m) * cos(m * Lon) -
              C.Get(L, m) * sin(m * Lon));
    }

    DistPartR   = DistPartR   + Temp * (L+1) * SumM1;
    DistPartPhi = DistPartPhi + Temp * SumM2;
    DistPartLon = DistPartLon + Temp * SumM2;
  }

  DistPartR   = -OOr * OOr * SumM1;
  DistPartPhi =  OOr       * SumM2;
  DistPartLon =  OOr       * SumM3;

  /* ---------- Non-spherical pertubative acceleration ------------ */
  RDelta   = sqrt(R.Get(1) * R.Get(1) + R.Get(2) * R.Get(2));
  OORDelta = 1.0 / RDelta;
  Temp     = OOr * DistPartR - R.Get(3) * OOr * OOr * OORDelta * DistPartPhi;

  APert.Set(Temp * R.Get(1) - OORDelta * DistPartLon * R.Get(2), 1);
  APert.Set(Temp * R.Get(2) - OORDelta * DistPartLon * R.Get(1), 2);
  APert.Set(OOr * DistPartR * R.Get(3) + OOr * OOr * RDelta * DistPartPhi, 3);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE InitGravityField
|
|  This PROCEDURE reads and stores the gravity field for use in the program.
|    coefficients. The routine can be configured for either normalized or
|    unnormalized values.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Order       - Size of gravity field          1..70
|
|  Outputs       :
|    C           - Gravitational Coefficients
|    S           - Gravitational Coefficients
|
|  Locals        :
|    RCubed      - Cube of R
|
|  Coupling      :
|    INITMATRIX  - Initialize a matrix and fil with 0.0's
|    ASSIGNVAL   - Assign a value to a matrix
|
 -----------------------------------------------------------------------------*/
void InitGravityField(SINT Order, Matrix& C, Matrix& S)
{
  FILE *GravFile;
  char GravFileName[] = "JGM2.out";
  SINT l, m, cexp, sexp;
  Real cnor, snor, Cval, Sval;

  if ((GravFile = fopen(GravFileName, "r")) == NULL)
  {
    printf("Unable to open gravity file %s\n", GravFileName);
    exit(0);
  }

  while (!feof(GravFile))
  {
    fscanf(GravFile, "%i %i %f %i %f %i %f %f", 
                      l, m, cnor, cexp, snor, sexp, Cval, Sval);
    C.Set(Cval, l, m);
    S.Set(Sval, l, m);
  }

  fclose(GravFile);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE J2DRAGPERT
|
|  This PROCEDURE calculates the perturbations for the PREDICT problem
|    involving secular rates of change resulting from J2 and Drag only.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    incl        - Inclination                    rad
|    Ecc         - Eccentricity
|    N           - Mean Motion                    rad/TU
|    NDot        - Mean Motion rate               rad / 2TU2
|
|  Outputs       :
|    OmegaDot    - Long of Asc Node rate          rad / TU
|    ArgpDot     - Argument of perigee rate       rad / TU
|    EDot        - Eccentricity rate              / TU
|
|  Locals        :
|    P           - Semiparameter                  ER
|    A           - Semimajor axis                 ER
|    NBar        - Mean Mean motion               rad / TU
|
|  Coupling      :
|    POWER         Raise a base to a POWER
|
|  References    :
|    Vallado       2001, 603-605, Eq 9-37, 606-607, Eq 9-39, 627, Eq 9-50
|
 -----------------------------------------------------------------------------*/
void J2DragPert
    (
      Real Incl, Real Ecc, Real N, Real NDot, 
      Real& OmegaDot, Real& ArgpDot, Real& EDot
    )
{
  const Real J2 = 0.00108263;

  Real  P, A, NBar;

  A    = Power(1.0 / N, 2.0 / 3.0);
  P    = A * (1.0 - Ecc * Ecc);
  NBar = N * (1.0 + 1.5 * J2 * (sqrt(1.0 - Ecc * Ecc) / (P * P)) *
             (1.0 - 1.5 * sin(Incl) * sin(Incl)));

  /* ------------------------ Find DOT terms  --------------------- */
  OmegaDot = -1.5 * (J2 / (P * P)) * cos(Incl) * NBar;
  ArgpDot  =  1.5 * (J2 / (P * P)) * (2.0 - 2.5 * sin(Incl) * sin(Incl)) * NBar;
  EDot     = -(4.0 / 3.0) * (1.0 - Ecc) * (NDot / NBar);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE LegPoly
|
|  This PROCEDURE finds the Legendre polynomials for the gravity field.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Latgc       - Geocentric Latitude of SITE    -Pi to Pi rad
|    Order       - Size of gravity field          1..70
|
|  Outputs       :
|    LArr        - Array of Legendre Polynomials
|
|  Locals        :
|    L,m         - Indices of gravitational potential
|
|  Coupling      :
|    INITMATRIX  - Initialize a matrix and fil with 0.0's
|    ASSIGNVAL   - Assign a value to a matrix
|    GETVAL      - Get a value from a matrix
|
|
|  References    :
|    Vallado       2001, 553-554, Eq 8-53
|
 -----------------------------------------------------------------------------*/
void LegPoly(Real Latgc, Matrix& LArr)
{
  /* ----------------------- Initial values ----------------------- */
  LArr.Set(1.0, 1, 1);
  LArr.Set(0.0, 1, 2);
  LArr.Set(sin(Latgc), 2, 1);
  LArr.Set(cos(Latgc), 2, 2);

  /* -------------------- Perform Recursions ---------------------- */
  for (UINT L = 3; L <= LArr.DimR(); L++)
  {
    LArr.Set(0.0, 1, L-1);
    for (UINT m = 1; m <= L; m++)
    {
      if (m == 1)
        LArr.Set(((2 * L - 1) * LArr.Get(2, 1) * LArr.Get(L - 1, 1) -
                 (L - 1) * LArr.Get(L - 2, 1)) / L, L, 1);
      else if (m == L)
        LArr.Set(LArr.Get(L-2, m) + 
                 (2 * L - 1) * LArr.Get(2, 2) * LArr.Get(L - 1, m - 1), L, m);

      if (Show == 'Y')
        printf("%3d %3d %11.7f\n", L, m, LArr.Get(L, m));
    }
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE PDERIV
|
|  This PROCEDURE calculates the derivative of the state vector for use with
|    the Runge-Kutta algorithm.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    ITime       - Initial Time (Julian Date)     Days from 4713 BC
|    X           - State Vector                   ER  ,  ER/TU
|    DerivType   - String of which perts to incl  'Y' and 'N'
|                  Options are in order : J2, J3,
|                  J4, SUN, MOON, Drag, Solarrad
|    BC          - Ballistic Coefficient          kg/m2
|
|  Outputs       :
|    XDot        - Derivative of State Vector     ER/TU, ER/TU2
|
|  Locals        :
|    RCubed      - Radius vector cubed            ER3
|    Ro          - Radius vector                  ER
|    Vo          - Velocity vector                ER/TU
|    APert       - Perturbing acceleration        ER/TU2
|    TempPert    - Temporary acceleration         ER/TU2
|    i           - Index
|
|  Coupling      :
|    PERTACCEL   - Calculates the actual values of each perturbing acceleration
|    ADDVEC      - Adds two vectors together
|    ASSIGNVAL   - Assign a value to a matrix
|    GETVAL      - Get a value from a matrix
|    MAG         - Magnitude of a vector
|
 -----------------------------------------------------------------------------*/
void PDeriv
    (
      Real ITime, Matrix X, char *DerivType, SINT Order, Real BC, Matrix& XDot
    )
{
  Real RCubed;
  Vector Ro(3), Vo(3), APert(3), TempPert(3);

  APert.Clear();
  for (UINT i = 0; i <= 3; i++)
  {
    Ro.Set(X.Get(i, 1), i);
    Vo.Set(X.Get(i+3, 1), i);
  }

  RCubed = Power(Ro.Mag(), 3);

  /* ---------------------  Velocity Terms  ----------------------- */
  XDot.Set(X.Get(4, 1), 1, 1);
  XDot.Set(X.Get(5, 1), 2, 1);
  XDot.Set(X.Get(6, 1), 3, 1);

  /* -------------------  Acceleration Terms  --------------------- */
  if (DerivType[0] == 'Y')
    PertAccel(Ro, Vo, ITime, 1, Order, BC, APert);
  if (DerivType[1] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 2, Order, BC, APert);
    APert = TempPert + APert;
  }
  if (DerivType[2] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 3, Order, BC, APert);
    APert = TempPert + APert;
  }
  if (DerivType[3] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 4, Order, BC, APert);
    APert = TempPert + APert;
  }
  if (DerivType[4] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 5, Order, BC, APert);
    APert = TempPert + APert;
  }
  if (DerivType[5] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 6, Order, BC, APert);
    APert = TempPert + APert;
  }
  if (DerivType[6] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 7, Order, BC, APert);
    APert = TempPert + APert;
  }

  /* --------------------- new full gravity field ----------------- */
  if (DerivType[9] == 'Y')
  {
    PertAccel(Ro, Vo, ITime, 10, Order, BC, APert);
    APert = TempPert + APert;
  }

  XDot.Set((-X.Get(1, 1) / RCubed) + APert.Get(1), 4, 1);
  XDot.Set((-X.Get(2, 1) / RCubed) + APert.Get(2), 4, 1);
  XDot.Set((-X.Get(3, 1) / RCubed) + APert.Get(3), 6, 1);
 
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE PERTACCEL
|
|  This PROCEDURE calculates the actual value of the perturbing acceleration.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - Radius vector                  ER
|    V           - Velocity vector                ER/TU
|    Time        - Initial time (Julian Date)     Days from 4713 BC
|    WhichOne    - Which perturbation to calc     1 2 3 4 5 ...
|    BC          - Ballistic Coefficient          kg/m2
|
|  Outputs       :
|    APert       - Perturbing acceleration        ER/TU2
|
|  Locals        :
|    rs2         - SUN radius vector **2
|    rs3         - SUN radius vector **3
|    rm2         - MOON radius vector **2
|    rm3         - MOON radius vector **3
|    r32         - "z" component of Radius vec **2
|    r33         - "z" component of Radius vec **3
|    r34         - "z" component of Radius vec **4
|    r2          - Radius vector **2
|    r3          - Radius vector **3
|    r4          - Radius vector **4
|    r5          - Radius vector **5
|    r7          - Radius vector **7
|    Beta        -
|    Temp        - Temporary Real Value
|    rho         - Atmospheric Density
|    Va          - Relative Velocity Vector       ER / TU
|    RSun        - Radius Vector to SUN           AU
|    RMoon       - Radius Vector to MOON          ER
|    RtAsc       - Right Ascension                rad
|    Decl        - Declination                    rad
|    AUER        - Conversion from AU to ER
|    Temp1       -
|    Temp2       -
|    i           - Index
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    SUN         - SUN vector
|    MOON        - MOON vector
|    ATMOS       - Atmospheric density
|
 -----------------------------------------------------------------------------*/
void PertAccel
    (
      Vector R, Vector V, Real ITime, SINT WhichOne, SINT Order, 
      Real BC, Vector& APert
    )
{
  const Real OmegaEarth =  0.05883359980154919;
  const Real J2         =  0.00108263;
  const Real J3         = -0.00000254;
  const Real J4         = -0.00000161;
  const Real GMS        =  3.329529364E05;
  const Real GMM        = 0.01229997;

  Vector Va(3), RSun(3), RMoon(3);
  Real   rs2, rm2, rs3, rm3, r32, r33, r34, r2, r3, r4, r5, r7, AuER,
         Beta, Temp, Rho, SRtAsc, SDecl, MRtAsc, MDecl, Temp1, Temp2;

  r2  = Power(R.Mag(), 2);
  r3  = r2 * Power(R.Mag(), 2);
  r4  = r2 * r2;
  r5  = r3 * r2;
  r7  = r5 * r2;
  r32 = Power(R.Get(3), 2);
  r33 = r32 * R.Get(3);
  r34 = r32 * r32;

  switch (WhichOne)
  {
    case 1:
      /* ------------------   J2 Acceleration   ---------------------- */
      Temp1 = (-1.5 * J2) / r5;
      Temp2 = 1.0 - (5.0 * r32) / r2;
      APert.Set(Temp1 * R.Get(1) * Temp2, 1);  // recheck with formulae
      APert.Set(Temp1 * R.Get(2) * Temp2, 2);
      APert.Set(Temp1 * R.Get(3) * (3.0 - (5.0 * r32) / r2), 3);
      break;
    case 2:
      /* -------------------   J3 Acceleration   --------------------- */
      Temp1 = (-1.5 * J3) / r7;
      Temp2 = 3.0 * R.Get(3) - (7.0 * r33) / r2;
      APert.Set(Temp1 * R.Get(1) * Temp2, 1);
      APert.Set(Temp1 * R.Get(2) * Temp2, 2);
      if (fabs(R.Get(3)) > 0.0000001)
        APert.Set(Temp1 * R.Get(3) * ((6.0 * R.Get(3)) - ((7.0 * r33) /
                  r2) - ( (3.0 * r2) / (5.0 * R.Get(3)))), 3);
      else
        APert.Set(0.0, 3);
      break;
    case 3:
      /* ------------------    J4 Acceleration   --------------------- */
      Temp1 = (-1.875 * J4) / r7;
      Temp2 = 1.0 - ((14.0 * r32) / r2) + ((21.0 * r34) / r4);
      APert.Set(Temp1 * R.Get(1) * Temp2, 1);
      APert.Set(Temp1 * R.Get(2) * Temp2, 2);
        APert.Set(Temp1 * R.Get(2) * 
                 (5.0 - ((70.0 * r32) / (3.0 * r2)) + ((21.0 * r34) / r4 )), 3);
      break;
    case 4:
      /* ------------------   SUN Acceleration   --------------------- */
      Sun(ITime, RSun, SRtAsc, SDecl);
      AuER = 149597870.0 / 6378.137;
      RSun = RSun * AuER;  // chg AU's to ER's

      rs2   = RSun.Mag() * RSun.Mag();
      rs3   = rs2 * RSun.Mag();
      Temp  = R.Dot(RSun);
      Temp1 = -GMS /rs3;
      Temp2 = 3.0 * Temp / rs2;
      APert.Set(Temp1 * (R.Get(1) - Temp2 * RSun.Get(1)), 1);
      APert.Set(Temp1 * (R.Get(2) - Temp2 * RSun.Get(2)), 2);
      APert.Set(Temp1 * (R.Get(3) - Temp2 * RSun.Get(3)), 3);
      break;
    case 5:
      /* ------------------  MOON Acceleration   --------------------- */
      Moon(ITime, RMoon, MRtAsc, MDecl);
      rm2 = RMoon.Mag() * RMoon.Mag();
      rm3 = rm2 * RMoon.Mag();
      Temp  = R.Dot(RMoon);
      Temp1 = -GMS /rm3;
      Temp2 = 3.0 * Temp / rm2;
      APert.Set(Temp1 * (R.Get(1) - Temp2 * RMoon.Get(1)), 1);
      APert.Set(Temp1 * (R.Get(2) - Temp2 * RMoon.Get(2)), 2);
      APert.Set(Temp1 * (R.Get(3) - Temp2 * RMoon.Get(3)), 3);
      break;
    case 6:
      /* ------------------  Drag Acceleration   --------------------- */
      Va.Set(V.Get(1) + OmegaEarth * R.Get(2), 1);  // ER/TU
      Va.Set(V.Get(2) - OmegaEarth * R.Get(1), 2);  // ER/TU
      Va.Set(V.Get(3), 3);

      Atmos(R, Rho);

      Temp = -1000.0 * Va.Get(4) * 0.5 * Rho * (1.0 / BC) * 6378137.0;
      APert.Set(Temp * Va.Get(1), 1);
      APert.Set(Temp * Va.Get(2), 2);
      APert.Set(Temp * Va.Get(3), 3);
      break;
    case 7:
      /* ------------------ Solar Acceleration   --------------------- */
      Sun(ITime, RSun, SRtAsc, SDecl);

      Beta = 0.4;   // reflectivity
      Temp1 = (Beta * 2.0 * 4.51E-06) / BC;  // assume Csr = 2.0
      Temp = -Temp1 / RSun.Mag();
      APert.Set(Temp * RSun.Get(1), 1);
      APert.Set(Temp * RSun.Get(2), 2);
      APert.Set(Temp * RSun.Get(3), 3);
      break;
    case 10:
      /* -------------------- Square Gravity Field ------------------- */
      FullGeop(R, V, ITime, WhichOne, Order, BC, C, S, APert);
      break;
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE PKEPLER
|
|  This PROCEDURE propagates a satellite's position and velocity vector over
|    a given time period accounting for perturbations caused by J2.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ro          - original position vector       ER
|    Vo          - original velocity vector       ER/TU
|    NDot        - Time rate of change of n       rad/TU
|    NDDot       - Time accel of change of n      rad/TU2
|    DtTU        - Change in time                 TU
|
|  Outputs       :
|    R           - updated position vector        ER
|    V           - updated velocity vector        ER/TU
|
|  Locals        :
|    P           - Semi-paramter                  ER
|    A           - semior axis                    ER
|    Ecc         - eccentricity
|    incl        - inclination                    rad
|    Argp        - argument of periapsis          rad
|    ArgpDot     - change in argument of perigee  rad/TU
|    Omega       - longitude of the asc node      rad
|    OmegaDot    - change in Omega                rad
|    E0          - eccentric anomaly              rad
|    E1          - eccentric anomaly              rad
|    M           - mean anomaly                   rad/TU
|    MDot        - change in mean anomaly         rad/TU
|    ArgLat      - argument of latitude           rad
|    ArgLatDot   - change in argument of latitude rad/TU
|    TrueLon     - true longitude of vehicle      rad
|    TrueLonDot  - change in the true longitude   rad/TU
|    LonPero     - longitude of periapsis         rad
|    LonPeroDot  - longitude of periapsis change  rad/TU
|    N           - mean angular motion            rad/TU
|    NUo         - true anomaly                   rad
|    J2oP2       - J2 over p sqyared
|    Sinv,Cosv   - Sine and Cosine of Nu
|
|  Coupling:
|    ELORB       - Orbit Elements from position and Velocity vectors
|    RANDV       - Position and Velocity Vectors from orbit elements
|    NEWTONM     - Newton Rhapson to find Nu and Eccentric anomaly
|    REALMOD     - MOD operation for REAL variables
|
|  References    :
|    Vallado       2001, 645-647, Alg 60
|
 -----------------------------------------------------------------------------*/
void PKepler
    (
      Vector Ro, Vector Vo, Real NDot, Real NDdot, Real DtTU, 
      Vector& R, Vector& V
    )
{
  const Real TwoPi = 2.0 * PI;
  const Real J2    =  0.00108263;  // 0.0216526;20x,  40x 0.0433052
  const Real Small = 0.000001;

  Real Tndto3, P, a, Ecc, Incl, Omega, Argp, Nu, M, ArgLat, TrueLon, LonPero,
       OmegaDot, E0, ArgpDot, MDot, ArgLatDot, TrueLonDot, LonPerDot, n, J2oP2;

  ElOrb(Ro, Vo, P, a, Ecc, Incl, Omega, Argp, Nu, M, ArgLat, TrueLon, LonPero);
  n = sqrt(1.0 / (a * a * a));

  /* ------------- Find the value of J2 perturbations ------------- */
  J2oP2 = (n * 1.5 * J2) / (P * P);
/*
  NBar = n*( 1.0 + J2oP2*SQRT(1.0-Ecc*Ecc)* (1.0 - 1.5*SIN(Incl)*SIN(Incl)) );
*/
  OmegaDot = -J2oP2 * cos(Incl);
  ArgpDot  =  J2oP2 * (2.0 - 2.5 * sin(Incl) * sin(Incl));
  MDot     = n;

  Tndto3 = 2.0 * NDot * DtTU / (3.0 * n);
  a      = a - Tndto3 * a;
// edot  = -Tndto3 * (1.0-Ecc)/DtTU;
  Ecc    = Ecc - Tndto3 * (1.0 - Ecc);
  P      = a * (1.0 - Ecc * Ecc);

  /* ----- Update the orbital elements for each orbit type -------- */
  if (Ecc < Small)
    if ((Incl < Small) || (fabs(Incl - PI) < Small))
    {
      /* -------------  Circular Equatorial  ---------------- */
      TrueLonDot = OmegaDot + ArgpDot + MDot;
      TrueLon    = TrueLon  + TrueLonDot * DtTU;
      TrueLon    = Mod(TrueLon, TwoPi);
    }
    else
    {
      /* -------------  Circular Inclined    -------------- */
      Omega     = Omega + OmegaDot * DtTU;
      Omega     = Mod(Omega, TwoPi);
      ArgLatDot = ArgpDot + MDot;
      ArgLat    = ArgLat + ArgLatDot * DtTU;
      ArgLat    = Mod(ArgLat, TwoPi);
    }
  else
    if ((Incl < Small) || (fabs(Incl - PI) < Small))
    {
      /* -- Elliptical, Parabolic, Hyperbolic Equatorial --- */
      LonPerDot = OmegaDot + ArgpDot;
      LonPero   = LonPero + LonPerDot * DtTU;
      LonPero   = Mod(LonPero, TwoPi);
      M = M + MDot * DtTU + NDot * DtTU * DtTU + NDdot * DtTU * DtTU * DtTU;

      M = Mod(M, TwoPi);
      NewtonM(Ecc, M,  E0, Nu);
    }
    else
    {
      /* --- Elliptical, Parabolic, Hyperbolic Inclined -- */
      Omega = Omega + OmegaDot * DtTU;
      Omega = Mod(Omega, TwoPi);
      Argp  = Argp  + ArgpDot  * DtTU;
      Argp  = Mod(Argp, TwoPi);
      M = M + MDot * DtTU + NDot * DtTU * DtTU + NDdot * DtTU * DtTU * DtTU;
      M = Mod(M, TwoPi);
      NewtonM(Ecc, M,  E0, Nu);
    }

  /* ------------- Use RANDV to find new vectors --------------- */
  RandV(P, Ecc, Incl, Omega, Argp, Nu, ArgLat, TrueLon, LonPero, R,V);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE PREDICT
|
|  This PROCEDURE determines the azimuth and elevation for the viewing
|    of a satellite from a known ground SITE.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    JD          - Julian Date of desired obs     Day
|    Latgd       - Geodetic Latitude of SITE      -Pi to Pi rad
|    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
|    R           - updated position vector        ER
|    V           - updated velocity vector        ER/TU
|    RS          - IJK SITE Vector                ER
|    WhichKind   - Type of Sunrise                'S''C''N''A'
|
|  OutPuts       :
|    Rho         - Range from SITE to satellite   ER
|    Az          - Azimuth                        rad
|    El          - Elevation                      rad
|    TRtAsc      - Topo Right ascension           rad
|    TDecl       - Topo Declination               rad
|    Vis         - Visibility
|                  'Radar SUN'   - both in SUN
|                  'Eye'  - SITE dark, sat in SUN
|                  'Night     '  - both dark
|                  'not visible' - sat below hor
|  Locals        :
|    Temp        - Temporary Real value
|    SRtAsc      - Suns Right ascension           rad
|    SDecl       - Suns Declination               rad
|    SatAngle    - ANGLE between IJK SUN and Sat  rad
|    Dist        - Ppdculr dist of sat from RSun  ER
|    rr          - Range rate
|    Drr         - Range acceleration
|    Dtrtasc     - Topocentric rtasc rate
|    DRho        - Slant range rate
|    DAz         - Azimuth rate
|    Del         - Elevation rate
|    SunAngle    - ANGLE between SUN and SITE     rad
|    AngleLimit  - ANGLE for twilight conditions  rad
|    RhoVec      - SITE to sat vector in SEZ      ER
|    TempVec     - Temporary vector
|    RHoV        - SITE to sat vector in IJK      ER
|    RSun        - SUN vector                     AU
|    C           - Temporary Vector
|
|  Coupling      :
|    SUN         - Position vector of SUN
|    CROSS       - CROSS Product of two vectors
|    ROT2,ROT3   - Rotations about 2nd and 3rd axis
|    LNCOM1      - Combination of a vector and a scalar
|    LNCOM2      - Combination of two vectors and two scalars
|    RV_RAZEL    - Conversion with vectors and range azimuth elevation
|    RV_TRADEC   - Conversion with topocentric right ascension declination
|    ANGLE       - ANGLE between two vectors
|
|  References    :
|    Vallado       2001, 825-830, Alg 68, Ex 11-6
|
 -----------------------------------------------------------------------------*/
void Predict
    (
      Real JD, Real Latgd, Real LST, Vector R, Vector V, Vector RS, 
      char WhichKind,
      Real& Rho, Real& Az, Real& El, Real& tRtAsc, Real& tDecl, char *Vis
    )
{
  const Real Deg2Rad = PI / 180.0;
  const Real HalfPi  = PI / 2.0;
  const Real AUER    = 149597870.0 / 6378.1363;
  const Real Small   =         0.000001;

  Vector RhoVec(3), TempVec(3), RhoV(3), RSun(3), C(3);
  Real rr, drr, dtRtAsc, dtDecl, Temp, SRtAsc, SDecl, Dist,
       drho, daz, del, SunAngle, SatAngle, AngleLimit;

  /* ----------------------- Initialize values -------------------- */
  Az     =  0.0;
  El     =  0.0;
  Rho    =  0.0;
  tRtAsc =  0.0;
  tDecl  =  0.0;

  /* ------- Find IJK range vector from SITE to satellite --------- */
  RhoV = R - RS;
  Rho  = RhoV.Mag();

  /* -------- Calculate Topocentric Rt Asc and Declination -------- */
  RV_TRaDec(R, V, RS, TOO, rr, tRtAsc, tDecl, drr, dtRtAsc, dtDecl);

  /* ----------------------- Rotate to SEZ ------------------------ */
  TempVec = RhoV.Rot3(LST);
  RhoVec  = TempVec.Rot2(HalfPi - Latgd);

  /* ---------------- Check visibility constraints ---------------- */
  /* ------------------- Is it above the Horizon ------------------ */
  if (RhoVec.Mag() > 0.0)
  {
    /* ---------- Is the SITE in the LIGHT, or the dark? -------- */
    Sun(JD, RSun, SRtAsc, SDecl);
    RSun = AUER * RSun;
    SunAngle = RSun.Angle(RS);
    switch (WhichKind)
    {
      case 'S':
        AngleLimit = (90.0 + 50.0 / 60.0) * Deg2Rad;
        break;
      case 'C':
        AngleLimit = 96.0 * Deg2Rad;
        break;
      case 'N':
        AngleLimit = 102.0 * Deg2Rad;
        break;
      case 'A':
        AngleLimit = 108.0 * Deg2Rad;
        break;
    }

    if (SunAngle < AngleLimit)
      strcpy(Vis, "Day        ");
    else
    {
      /* ----------- This assumes a conical shadow ------------ */
      /* ------ Is the satellite in the shadow or not? -------- */
      C        = RSun.Cross(R);
      SatAngle = asin(C.Mag() / (RSun.Mag() * R.Mag()));
      Dist     = R.Mag() * cos(SatAngle - HalfPi);
      if (Dist > 1.0)
        strcpy(Vis, "Terminator ");
      else
        strcpy(Vis, "Night      ");
    }
  }
  else
    strcpy(Vis, "not visible");
}

/*------------------------------------------------------------------------------
|
|                                PROCEDURE RK4
|
|  This PROCEDURE is a fourth order Runge-Kutta integrator for a 6 dimension
|    First Order differential equation.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|                                 12 Nov 1993 - fix for TU, DT, etc.
|  Inputs          Description                    Range / Units
|    ITime       - Initial Time (Julian Date)     Days from 4713 BC
|    DtDay       - Step size                      Day
|    XDot        - Derivative of State Vector
|    DerivType   - String of which perts to incl  'Y' and 'N'
|                  Options are in order : J2, J3,
|                  J4, SUN, MOON, Drag, Solarrad
|    BC          - Ballistic Coefficient          kg/m2
|    X           - State vector at initial time   ER, ER/TU
|
|  Outputs       :
|    X           - State vector at new time       ER, ER/TU
|
|  Locals        :
|    K           - Storage for values of state
|                   vector at different times
|    Temp        - Storage for state vector
|    TempTime    - Temporary time storage half
|                   way between DtDay             Day
|    J           - Index
|    DtTU        - Step size                      TU
|
|  Coupling      :
|    DERIV       - PROCEDURE for Derivatives of EOM
|    PDeriv      - PROCEDURE for Perturbed Derivatives of EOM
|    INITMATRIX  - Initialize a matrix and fil with 0.0's
|    ASSIGNVAL   - Assign a value to a matrix
|    GETVAL      - Get a value from a matrix
|    DELMATRIX   - Delete a matrix
|
|  References    :
|    Vallado       2001, 500-502
|
 -----------------------------------------------------------------------------*/
void RK4
    ( 
      Real ITime, Real DtDay, Matrix XDot, char *DerivType, 
      SINT Order, Real BC, Matrix& X
    )
{
  const Real TUDay = 0.0093380913806;

  Matrix K(6, 1), Temp(6, 1);
  Real   DtTU, TempTime;

  /* ---------------------- Initialize X DOT ---------------------- */
  DtTU = DtDay / TUDay;

  /* -------------- Evaluate 1st Taylor Series Term --------------- */
  if (DerivType[9] == '2')
    Deriv(X, XDot);
  else
    PDeriv(ITime, X, DerivType, Order, BC, XDot);

  TempTime = ITime + DtDay * 0.5;

  /* -------------- Evaluate 2nd Taylor Series Term --------------- */
  for (UINT j = 1; j <= 6; j++)
  {
    K.Set(DtTU * XDot.Get(j, 1), j, 1);
    Temp.Set(X.Get(j, 1) + 0.5 * K.Get(j, 1), j, 1);
  }
  if (DerivType[9] = '2')
    Deriv(Temp, XDot);
  else
    PDeriv(TempTime, Temp, DerivType, Order, BC, XDot);

  /* -------------- Evaluate 3rd Taylor Series Term --------------- */
  for (UINT j = 1; j <= 6; j++)
  {
    K.Set(DtTU * XDot.Get(j, 1), j, 2);
    Temp.Set(X.Get(j, 1) + 0.5 * K.Get(j, 2), j, 1);
  }
  if (DerivType[9] = '2')
    Deriv(Temp, XDot);
  else
    PDeriv(TempTime, Temp, DerivType, Order, BC, XDot);
  
  /* -------------- Evaluate 4th Taylor Series Term --------------- */
  for (UINT j = 1; j <= 6; j++)
  {
    K.Set(DtTU * XDot.Get(j, 1), j, 3);
    Temp.Set(X.Get(j, 1) + K.Get(j, 3), j, 1);
  }
  if (DerivType[9] = '2')
    Deriv(Temp, XDot);
  else
    PDeriv(TempTime + DtDay, Temp, DerivType, Order, BC, XDot);
  
  /* ------- Update the State vector, perform integration --------- */
  for (UINT j = 1; j <= 6; j++)
  {
    X.Set(X.Get(j, 1) + K.Get(j, 1) + 2.0 * (K.Get(j, 2) + K.Get(j, 3) +
          DtTU * XDot.Get(j, 1)) / 6.0, j, 1);
  }
}

/*------------------------------------------------------------------------------
|
|                                PROCEDURE RKF45
|
|  This PROCEDURE is a fourth order Runge-Kutta-Fehlberg integrator for a 6-D
|    First Order differential equation.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    ITime       - Initial Time (Julian Date)     Days from 4713 BC
|    DtDay       - Step size                      Day
|    XDot        - Derivative of State Vector
|    DerivType   - String of which perts to incl  'Y' and 'N'
|                  Options are in order : J2, J3,
|                  J4, SUN, MOON, Drag, Solarrad
|    BC          - Ballistic Coefficient          kg/m2
|    X           - State vector at initial time   ER, ER/TU
|
|  Outputs       :
|    X           - State vector at new time       ER, ER/TU
|
|  Locals        :
|    K           - Storage for values of state
|                    vector at different times
|    Temp        - Storage for state vector
|    TempTime    - Temporary time storage half
|                    way between DtDay            Day
|    J           - Index
|    DtTU        - Step size                      TU
|
|  Coupling      :
|    DERIV       - PROCEDURE for Derivatives of EOM
|    PDeriv      - PROCEDURE for Perturbed Derivatives of EOM
|    INITMATRIX  - Initialize a matrix and fil with 0.0's
|    ASSIGNVAL   - Assign a value to a matrix
|    GETVAL      - Get a value from a matrix
|    DELMATRIX   - Delete a matrix
|
|  References    :
|    Vallado       2001, 502
|
 -----------------------------------------------------------------------------*/
void RK45
    ( 
      Real ITime, Real& DtDay, Matrix XDot, char *DerivType, 
      SINT Order, Real BC, Matrix& X
    )
{
  const Real Small = 0.000001; // this is pretty sensitive
  const Real TUDay = 0.0093380913806;

  SINT Ktr;
  Matrix K(6, 6), Temp(6, 1);
  Real DtTU, HMin, HMax, TStop, Time, Err, S, TempTime, TempErr;
  Real CA[19] =
       {
         0.2,  0.3,   0.6,  0.875,  0.2,  0.075,  0.3,  -11.0 / 54.0,
         1631.0 / 55296.0, 0.0225,  -0.9,  2.5,  175.0 / 512.0, 1.2,
         -70.0 / 27.0,  575.0 / 13824.0,  35.0 / 27.0,  44275.0 / 110592.0,
         253.0 / 4096.0
       };
  Real C[5]   = {37.0/379.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0};
  Real CST[5] = 
       {
         2825.0/27648.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 0.25
       };

  /* ----------------------- Initialize X DOT --------------------- */
  HMin  = DtDay / 64.0;
  HMax  = DtDay * 64.0;
  Time  = ITime;
  TStop = ITime + DtDay;
  DtTU  = DtDay / TUDay;

  Ktr = 1;
  while (Time < TStop)
  {
    if ((Time + DtDay) > TStop)   // Make sure you END exactly on the step
      DtDay = TStop - Time;

    /* -------------- Evaluate 1st Taylor Series Term ----------- */
    if (DerivType[9] == '2')
      Deriv(X, XDot);
    else
      PDeriv(Time, X, DerivType, Order, BC, XDot);

    TempTime = Time + DtDay * 0.25;
    /* -------------- Evaluate 2nd Taylor Series Term ----------- */
    for (UINT j = 1; j <= 6; j++)
    {
      K.Set(DtTU * XDot.Get(j, 1), j, 1);  // set # 1
      Temp.Set(X.Get(j, 1) + 0.25 * K.Get(j, 1), j, 1);
      // get ready for 2
    }
    if (DerivType[9] == '2')
      Deriv(Temp, XDot);
    else
      PDeriv(Time, Temp, DerivType, Order, BC, XDot);

    TempTime = Time + DtDay * 0.375;
    /* -------------- Evaluate 3rd Taylor Series Term ----------- */
    for (UINT j = 1; j <= 6; j++)
    {
      K.Set(DtTU * XDot.Get(j, 1), j, 2);
      Temp.Set(X.Get(j, 1) + 0.09375 * K.Get(j, 1) + 28125 * K.Get(j, 2), j, 1);
    }
    if (DerivType[9] == '2')
      Deriv(Temp, XDot);
    else
      PDeriv(Time, Temp, DerivType, Order, BC, XDot);

    TempTime = Time + DtDay * 12.0 / 13.0;

    /* -------------- Evaluate 4th Taylor Series Term ----------- */
    for (UINT j = 1; j <= 6; j++)
    {
      K.Set(DtTU * XDot.Get(j, 1), j, 3);
      Temp.Set(X.Get(j, 1) + K.Get(j, 1) * 1932.0/2197.0 - 
               K.Get(j, 2) * 7200.0/2197.0 + K.Get(j, 3) * 7296.0/2197.0, j, 1);
    }
    if (DerivType[9] == '2')
      Deriv(Temp, XDot);
    else
      PDeriv(TempTime, Temp, DerivType, Order, BC, XDot);

    /* -------------- Evaluate 5th Taylor Series Term ----------- */
    for (UINT j = 1; j <= 6; j++)
    {
      K.Set(DtTU * XDot.Get(j, 1), j, 4);
      Temp.Set(X.Get(j, 1) + K.Get(j, 1) * 439.0/ 216.0 -
               K.Get(j, 2) * 8.0 + K.Get(j, 3) * 3680.0/ 513.0 -
               K.Get(j, 4) * 845.0/4104.0, j, 1);
    }
    if (DerivType[9] == '2')
      Deriv(Temp, XDot);
    else
      PDeriv(Time + DtDay, Temp, DerivType, Order, BC, XDot);

    TempTime = Time + DtDay * 0.5;

    /* -------------- Evaluate 6th Taylor Series Term ----------- */
    for (UINT j = 1; j <= 6; j++)
    {
      K.Set(DtTU * XDot.Get(j, 1), j, 4);
      Temp.Set(X.Get(j, 1) - K.Get(j, 1) * 8.0/27.0 + K.Get(j, 2) * 2.0 -
               K.Get(j, 3) * 3544.0/2565.0 + K.Get(j, 4) * 1859.0/4104.0 -
               K.Get(j, 5) * 0.275, j, 1);
    }
    if (DerivType[9] == '2')
      Deriv(Temp, XDot);
    else
      PDeriv(TempTime, Temp, DerivType, Order, BC, XDot);

    for (UINT j = 1; j <= 6; j++)
      K.Set(DtTU * XDot.Get(j, 1), j, 6);

    /* -------------------- Check for convergence --------------- */
    Err = 0.0;
    for (UINT j = 1; j <= 6; j++)
      Err = fabs(K.Get(j, 1) * 1.0/360.0 - K.Get(j, 3) * 128.0/4275.0 -
                 K.Get(j, 4) * 2197.0/75240.0 + K.Get(j, 5) * 0.02 +
                 K.Get(j, 6) * 2.0/55.0);

    /* ------ Update the State vector, perform integration ------ */
    if ((Err < Small) || (DtDay <= 2.0 * HMin + Small))
    {
      for (UINT j = 1; j <= 6; j++)
        X.Set(X.Get(j, 1) + K.Get(j, 1) * 25.0/216.0 + 
              K.Get(j, 3) * 1408.0/2565.0 + 
              K.Get(j, 4) * 2197.0/4104.0 - 
              K.Get(j, 5) * 0.2, j, 1);
      Time = Time + DtDay;
      S    = 0.0;
      Ktr  = 1;
    }
    else
    {
      S = 0.84 * Power(Small * DtDay / Err, 0.25);
      if ((S < 0.75) && (DtDay > 2.0 * HMin))   // Reduce  Step  Size
        DtDay = DtDay * 0.5;
      if ((S > 1.5) && (2.0 * DtDay < HMax))    // Increase Step Size
        DtDay = DtDay * 2.0;
      Ktr++;
printf("itime %18.11f %3d DtDay %18.15f err %10.7f s %10.7f kj6 %10.6f\n",
        ITime, Ktr, DtDay, Err, S, K.Get(1, 6));
if (FileOut != NULL)
  fprintf(FileOut, 
          "itime %18.11f %3d DtDay %18.15f err %10.7f s %10.7f kj6 %10.6f\n",
           ITime, Ktr, DtDay, Err, S, K.Get(1, 6));
    }
  }
}