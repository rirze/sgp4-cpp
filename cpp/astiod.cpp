/*----------------------------------------------------------------      


                               UNIT ASTIOD;


    This file contains fundamental Astrodynamic procedures and functions     
    relating to the Initial Orbit Determination techniques. See Ch 7 for     
    a complete discussion of these routines.                                 

                            Companion code for                               
               Fundamentals of Astrodyanmics and Applications                
                                     2001                                    
                              by David Vallado                               

       (H)               email valladodl@worldnet.att.net                    
       (W) 303-344-6037, email davallado@west.raytheon.com                   

       *****************************************************************     


                             MICROCOSM
                    email: bookproject@smad.com

                                 or

              Contact: Donna Klungle at (310) 726-4100


       *****************************************************************     

    Current :                                                                
              14 May 01  David Vallado                                       
                           2nd edition baseline                              

       ----------------------------------------------------------------     

                                  IMPLEMENTATION

       ----------------------------------------------------------------      */
#include "ast2body.h"
#include "astiod.h"

#include <math.h>

/* --------- only for testgau test ----------- */
Real lambda, BigT, Testamt;

/* Utility functions for LambertBattin, etc */
static Real k(Real v);
static Real See(Real v);

/*---------------------------------------------------------------------------
|
|                           PROCEDURE SITE
|
|  This PROCEDURE finds the position and velocity vectors for a SITE.  The
|    answer is returned in the Geocentric Equatorial (IJK) coordinate system.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
|    Alt         - Altitude                       ER
|    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
|
|  OutPuts       :
|    RSijk       - IJK SITE position vector       ER
|    VSijk       - IJK SITE velocity vector       ER/TU
|
|  Locals        :
|    EarthRate   - IJK Earth's rotation rate      rad/TU
|    SinLat      - Variable containing  SIN(Lat)  rad
|    Temp        - Temporary Real value
|    Rdel        - Rdel component of SITE vector  ER
|    Rk          - Rk component of SITE vector    ER
|    CEarth      -
|
|  Coupling      :
|    MAG           Magnitude of a vector
|    CROSS         CROSS product of two vectors
|
|  References    :
|    Vallado       2001, 404-407, Alg 47, Ex 7-1
|
 ----------------------------------------------------------------------------*/
void Site(Real Latgd, Real Alt, Real LST, Vector& RSijk, Vector& VSijk)
{
  const Real EESqrd     = 0.00669437999013;
  const Real OmegaEarth = 0.05883359980154919;

  Vector EarthRate(3);
  Real   SinLat, CEarth, Rdel, Rk;

  /* ---------------------  Initialize values   ------------------- */
  SinLat = sin(Latgd);
  EarthRate.Set(0.0, 1);
  EarthRate.Set(0.0, 2);
  EarthRate.Set(OmegaEarth, 2);

  /* -------  Find Rdel and Rk components of SITE vector  --------- */
  CEarth = 1.0 / sqrt(1.0 - (EESqrd * SinLat * SinLat));
  Rdel   = (CEarth + Alt) * cos(Latgd);
  Rk     = ((1.0 - EESqrd) * CEarth + Alt) * SinLat;

  /* ----------------  Find SITE position vector  ----------------- */
  RSijk.Set(Rdel * cos(LST), 1);
  RSijk.Set(Rdel * sin(LST), 2);
  RSijk.Set(Rk, 3);

  /* ----------------  Find SITE velocity vector  ----------------- */
  VSijk = EarthRate.Cross(RSijk);

  if (Show == 'Y')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "RdelRk %11.7f %11.7f CEarth %11.7f %11.7f\n",
                        Rdel, Rk, CEarth, (1.0 - EESqrd) * CEarth);
    }
}

/* ------------------- Angles-only techniques --------------------- */
/*------------------------------------------------------------------------------
|
|                           PROCEDURE ANGLESGAUSS
|
|  This PROCEDURE solves the problem of orbit determination using three
|    optical sightings.  The solution PROCEDURE uses the Gaussian technique.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    TRtAsc1      - Right Ascension #1            rad
|    TRtAsc2      - Right Ascension #2            rad
|    TRtAsc3      - Right Ascension #3            rad
|    TDecl1       - Declination #1                rad
|    TDecl2       - Declination #2                rad
|    TDecl3       - Declination #3                rad
|    JD1          - Julian Date of 1st sighting   Days from 4713 BC
|    JD2          - Julian Date of 2nd sighting   Days from 4713 BC
|    JD3          - Julian Date of 3rd sighting   Days from 4713 BC
|    RSijk        - IJK SITE position vector      ER
|
|  OutPuts        :
|    R            - IJK position vector at t2     ER
|    V            - IJK velocity vector at t2     ER / TU
|
|  Locals         :
|    L1           - Line of SIGHT vector for 1st
|    L2           - Line of SIGHT vector for 2nd
|    L3           - Line of SIGHT vector for 3rd
|    Tau          - Taylor expansion series about
|                   Tau ( t - to )
|    TauSqr       - Tau squared
|    t21t23       - (t2-t1) * (t2-t3)
|    t31t32       - (t3-t1) * (t3-t2)
|    i            - index
|    D            -
|    Rho          - Range from SITE to sat at t2  ER
|    RhoDot       -
|    DMat         -
|    RS1          - SITE vectors
|    RS2          -
|    RS3          -
|    EarthRate    - Velocity of Earth rotation
|    P            -
|    Q            -
|    OldR         -
|    OldV         -
|    F1           - F coefficient
|    G1           -
|    F3           -
|    G3           -
|    L2DotRS      -
|
|  Coupling       :
|    MAG          - Magnitude of a vector
|    Detrminant   - Evaluate the determinant of a matrix
|    FACTOR       - Find roots of a polynomial
|    MATMULT      - Multiply two matrices together
|    ASSIGNVAL    - Assign a value to a matrix
|    GETVAL       - Get a value from a matrix
|    INITMATRIX   - Initialize a matrix and fil with 0.0's
|    DELMATRIX    - Delete a matrix
|    GIBBS        - GIBBS method of orbit determination
|    HGIBBS       - Herrick GIBBS method of orbit determination
|    ANGLE        - ANGLE between two vectors
|
|  References     :
|    Vallado       2001, 417-421, Alg 49, Ex 7-2 (425-427)
|
 -----------------------------------------------------------------------------*/
void AnglesGauss
    (
      Real TDecl1,  Real TDecl2,  Real TDecl3, 
      Real TRtAsc1, Real TRtAsc2, Real TRtAsc3, 
      Real JD1,     Real JD2,     Real JD3,
      Vector RS1, Vector RS2, Vector RS3, Vector& R2, Vector& V2
    )
{
  const Real TUDay = 0.00933809017716;
  const Real Mu    = 1.0;
  const Real Small = 0.0000001;
  const Real Rad   = 180.0 / PI;

  char Error[12];
  Real Poly[16];
  Real Roots[15][2];
  Vector R1(3), R3(3), L1(3), L2(3), L3(3);
  Matrix LMatIi(3, 3), CMat(3, 3), RhoMat(3, 3), 
         LMatI(3, 3), RSMat(3, 3), LIR(3,3);
  Real rDot, Tau1, Tau3, u, uDot, p, f1, g1, f3, g3, a, ecc, incl, omega, argp,
       Nu, m, l, ArgPer, BigR2, a1, a1u, a3, a3u, D, D1, D2, c1, c3, L2DotRS,
       Rhoold1, Rhoold2, Rhoold3, Temproot, theta, theta1, copa, TauSqr;
  char tc, chg;
  Byte ScHi;


  /* ----------------------   Initialize   ------------------------ */
  CMat.Clear();
  LIR.Clear();
  LMatI.Clear();
  LMatIi.Clear();
  RhoMat.Clear();
  RSMat.Clear();

  JD1 = JD1 / TUDay;   // Convert days to TU
  JD2 = JD2 / TUDay;
  JD3 = JD3 / TUDay;
  /* ---- set middle to 0, deltas to other times ---- */
  Tau1 = JD1 - JD2;
  Tau3 = JD3 - JD2;

  /* ----------------  Find Line of SIGHT vectors  ---------------- */
  L1.Set(cos(TDecl1) * cos(TRtAsc1), 1);
  L1.Set(cos(TDecl1) * sin(TRtAsc1), 2);
  L1.Set(sin(TDecl1), 3);

  L2.Set(cos(TDecl2) * cos(TRtAsc2), 1);
  L2.Set(cos(TDecl2) * sin(TRtAsc2), 2);
  L2.Set(sin(TDecl2), 3);

  L3.Set(cos(TDecl3) * cos(TRtAsc3), 1);
  L3.Set(cos(TDecl3) * sin(TRtAsc3), 2);
  L3.Set(sin(TDecl3), 3);

  /* -------------- Find L matrix and determinant ----------------- */
  /* -------- Called LMatI since it is only used for determ ------- */
  for (UINT i = 1; i <= 3; i++)
  {
    L1.Set(LMatIi.Get(i, 1), i);
    L2.Set(LMatIi.Get(i, 2), i);
    L3.Set(LMatIi.Get(i, 3), i);
    RS1.Set(RSMat.Get(i, 1), i);
    RS2.Set(RSMat.Get(i, 2), i);
    RS3.Set(RSMat.Get(i, 3), i);
  }

  D = LMatIi.Determinant();

  /* ------------------- Now assign the inverse ------------------- */
  LMatI.Set(( L2.Get(2) * L3.Get(3) - L2.Get(3) * L3.Get(2)) / D, 1, 1);
  LMatI.Set((-L1.Get(2) * L3.Get(3) + L1.Get(3) * L3.Get(2)) / D, 2, 1);
  LMatI.Set(( L1.Get(2) * L2.Get(3) - L1.Get(3) * L2.Get(2)) / D, 3, 1);
  LMatI.Set((-L2.Get(1) * L3.Get(3) + L2.Get(3) * L3.Get(1)) / D, 1, 2);
  LMatI.Set(( L1.Get(1) * L3.Get(3) - L1.Get(3) * L3.Get(1)) / D, 2, 2);
  LMatI.Set((-L1.Get(1) * L2.Get(3) + L1.Get(3) * L2.Get(1)) / D, 3, 2);
  LMatI.Set(( L2.Get(1) * L3.Get(2) - L2.Get(2) * L3.Get(1)) / D, 1, 3);
  LMatI.Set((-L1.Get(1) * L3.Get(2) + L1.Get(2) * L3.Get(1)) / D, 2, 3);
  LMatI.Set(( L1.Get(1) * L2.Get(2) - L1.Get(2) * L2.Get(1)) / D, 3, 3);

  LIR = LMatI * RSMat;

  /* ------------- Find f and g series at 1st and 3rd obs --------- */
  /* speed by assuming circ sat vel for uDot here ??                */
  /* some similartities in 1/6t3t1 ...                              */
  /* ---- keep separated this time ----                             */
  a1  = Tau3 / (Tau3 - Tau1);
  a1u = (Tau3 * ((Tau3 - Tau1) * (Tau3 - Tau1) - Tau3 * Tau3)) / 
        (6.0 * (Tau3 - Tau1));
  a3  = -Tau1 / (Tau3 - Tau1);
  a3u = -(Tau1 * ((Tau3 - Tau1) * (Tau3 - Tau1) - Tau1 * Tau1)) / 
        (6.0 * (Tau3 - Tau1));

  /* ---- Form initial guess of r2 ---- */
  D1 = LIR.Get(2, 1) * a1  - LIR.Get(2, 2) + LIR.Get(2, 3) * a3;
  D2 = LIR.Get(2, 1) * a1u                 + LIR.Get(2, 3) * a3u;

  /* -------- Solve eighth order poly NOT same as LAPLACE --------- */
  L2DotRS = L2.Dot(RS2);
  Poly[ 0] = 1.0;  // r2
  Poly[ 1] = 0.0;
  Poly[ 2] = -(D1 * D1 + 2.0 * D1 * L2DotRS + RS2.Mag() * RS2.Mag());
  Poly[ 3] = 0.0;
  Poly[ 4] = 0.0;
  Poly[ 5] = -2.0* Mu * (L2DotRS * D2 + D1 * D2);
  Poly[ 6] = 0.0;
  Poly[ 7] = 0.0;
  Poly[ 8] = -Mu * Mu * D2 * D2;
  Poly[ 9] = 0.0;
  Poly[10] = 0.0;
  Poly[11] = 0.0;
  Poly[12] = 0.0;
  Poly[13] = 0.0;
  Poly[14] = 0.0;
  Poly[15] = 0.0;
  Factor(Poly, 8, (Real **) Roots);

  /* ------------------- Select the correct root ------------------ */
  BigR2 = 0.0;
  for (UINT j = 0; j < 8; j++)
  {
    if (fabs(Roots[j][1]) < Small)
    {
      Temproot = Roots[j][0] * Roots[j][0];
      Temproot = Temproot * Temproot * Temproot * Temproot +
                 Poly[2]  * Temproot * Temproot * Temproot +
                 Poly[5]  * Roots[j][0] * Temproot + Poly[8];
      if (FileOut != NULL)
      {
        fprintf(FileOut, "Root %d %0.7f + %0.7f j  value = %0.7f\n",
                          j, Roots[j][0], Roots[j][1], Temproot);
        fprintf(FileOut, "Root %d %0.7f + %0.7f j  value = %0.7f\n",
                          j, Roots[j][0], Roots[j][1], Temproot);
      }
      if (Roots[j][0] > BigR2)
        BigR2 = Roots[j][0];
    }
  }
  printf("input r2 ");
  scanf("%f\n", &BigR2);

  /* ------------- Solve matrix with u2 better known -------------- */
  u = Mu / (BigR2 * BigR2 * BigR2);

  c1 = a1 + a1u * u;
  c3 = a3 + a3u * u;
  CMat.Set(-c1, 1, 1);
  CMat.Set(1.0, 2, 1);
  CMat.Set(-c3, 3, 1);
  RhoMat = LIR * CMat;

  Rhoold1 =  RhoMat.Get(1, 1) / c1;
  Rhoold2 = -RhoMat.Get(2, 1);
  Rhoold3 =  RhoMat.Get(3, 1) / c3;

  if (FileOut != NULL)
    fprintf(FileOut, " Now start refining the guess\n");

  /* --------- Loop through the refining process ------------ */
  for (UINT ll = 1; ll <= 3; ll++)
  {
    if (FileOut != NULL)
      fprintf(FileOut, " Iteration # %2d\n", ll);
    /* ---- Now form the three position vectors ----- */
    for (UINT i = 1; i <= 3; i++)
    {
      R1.Set( RhoMat.Get(1, 1) * L1.Get(i) / c1 + RS1.Get(i), i);
      R2.Set(-RhoMat.Get(2, 1) * L2.Get(i)      + RS2.Get(i), i);
      R3.Set( RhoMat.Get(3, 1) * L3.Get(i) / c3 + RS3.Get(i), i);
    }

    Gibbs(R1,R2,R3, V2, theta, theta1, copa, Error);

    if ((strcmp(Error, "ok") == 0) && (copa < 1.0 / Rad))
    {
      /* HGibbs to get middle vector ---- */
      JD1 = JD1 * TUDay;   // Convert TU to days
      JD2 = JD2 * TUDay;
      JD3 = JD3 * TUDay;

      HerrGibbs(R1,R2,R3,JD1,JD2,JD3, V2,theta,theta1,copa,Error);

      if (FileOut != NULL)
        fprintf(FileOut, "hgibbs\n");
    }

    ElOrb(R2, V2, p, a, ecc, incl, omega, argp, Nu, m, u, l, ArgPer);

    if (FileOut != NULL)
    {
      fprintf(FileOut, "t123 %18.7f %18.7f %18.7f TU\n", JD1, JD2, JD3);
      fprintf(FileOut, "t123 %18.7f %18.7f %18.7f days\n", 
                        JD1 * TUDay, JD2 * TUDay, JD3 * TUDay);
      fprintf(FileOut, "los 1    %12.6f %12.6f %12.6f\n",
                        L1.Get(1), L1.Get(2), L1.Get(3));
      fprintf(FileOut, "los 2    %12.6f %12.6f %12.6f\n",
                        L2.Get(1), L2.Get(2), L2.Get(3));
      fprintf(FileOut, "los 3    %12.6f %12.6f %12.6f\n",
                        L3.Get(1), L3.Get(2), L3.Get(3));
      FilePrint(RSMat, " RSMat ", 3, FileOut);
    }
    LMatIi = LMatI.Inverse();
    LIR    = LMatI * LMatIi;
    LMatI.Display(" LMatI Matrix\n", 3);
    LIR.Display(" I Matrix\n", 6);
    printf("%14.7f\n", D);
    if (FileOut != NULL)
      FilePrint(LMatI,  " Lmat MatInv " ,3, FileOut);
    LIR = LMatI * LMatIi;
    LIR.Display(" should be I Matrix ", 6);
    LIR.Display(" LIR Matrix ", 6);
    if (FileOut != NULL)
    {
      fprintf(FileOut, "tau  %11.7f %11.7f TU %14.7f\n",Tau1, Tau3, u);
      fprintf(FileOut, "a13, u %11.7f %11.7f %11.7f%11.7f\n", a1, a3, a1u, a3u);
      fprintf(FileOut, "d1, d2 %11.7f %11.7f lDotr %11.7f\n", D1, D2, L2DotRS);
      fprintf(FileOut, "coeff %11.7f %11.7f %11.7f %11.7f\n",
                        Poly[0], Poly[2], Poly[5], Poly[8]);
      FilePrint(CMat, " c Matrix ", 3, FileOut);
      FilePrint(RhoMat, " Rho Matrix ", 3, FileOut);
      fprintf(FileOut, "r1 %11.7f %11.7f %11.7f %11.7f\n",
                        R1.Get(1), R1.Get(2), R1.Get(3), R1.Mag());
      fprintf(FileOut, "r2 %11.7f %11.7f %11.7f %11.7f\n",
                        R2.Get(1), R2.Get(2), R2.Get(3), R2.Mag());
      fprintf(FileOut, "r3 %11.7f %11.7f %11.7f %11.7f\n",
                        R3.Get(1), R3.Get(2), R3.Get(3), R3.Mag());
      fprintf(FileOut, "hgGau r2 %12.6 %12.6 %12.6 %s %11.7f ",
                 R2.Get(1), R2.Get(2), R2.Get(3), R2.Mag(), Error, copa * Rad);
      fprintf(FileOut, "%12.6 %12.6 %12.6\n", V2.Get(1), V2.Get(2), V2.Get(3));
      fprintf(FileOut, 
              "p=%11.7f  a%11.7f  e %11.7f i %11.7f Omeg %10.6f argp %10.6f\n",
                        p, a, ecc, incl * Rad, omega * Rad, argp * Rad);
      printf("p=%11.7f  a%11.7f  e %11.7f i %11.7f W%10.6f w%10.6f\n",
                        p, a, ecc, incl * Rad, omega * Rad, argp * Rad);
    }
    if (ll <= 2)
    {
      /* ---- Now get an improved estimate of the f and g series ---- */
      /* or can the analytic functions be found now??                 */
      u    = Mu / (R2.Mag() * R2.Mag() * R2.Mag());
      rDot = R2.Dot(V2) / R2.Mag();
      uDot = (-3.0 * Mu * rDot) / sqrt(R2.Mag()* R2.Mag());

      TauSqr = Tau1 * Tau1;
      f1 =  1.0 - 0.5 * u * TauSqr - 
                  (1.0 / 6.0) * uDot * TauSqr * Tau1 +
                  (1.0 / 24.0) * u * u * TauSqr * TauSqr +
                  (1.0 / 30.0) * u * uDot * TauSqr * TauSqr * Tau1;
      g1 = Tau1 - (1.0 / 6.0) * u * Tau1 * TauSqr - 
                  (1.0 / 12.0) * uDot * TauSqr * TauSqr +
                  (1.0 / 120.0) * u * u * TauSqr * TauSqr * Tau1 +
                  (1.0 / 120.0) * u * uDot * TauSqr * TauSqr * TauSqr;
      TauSqr = Tau3 * Tau3;
      f3 =  1.0 - 0.5 * u * TauSqr - 
                  (1.0 / 6.0) * uDot* TauSqr * Tau3 +
                  (1.0 / 24.0) * u * u * TauSqr * TauSqr +
                  (1.0 / 30.0) * u * uDot * TauSqr * TauSqr * Tau3;
      g3 = Tau3 - (1.0 / 6.0) * u * Tau3 * TauSqr -
                  (1.0 / 12.0) * uDot * TauSqr * TauSqr +
                  (1.0 / 120.0) * u * u * TauSqr * TauSqr * Tau3 +
                  (1.0 / 120.0) * u * uDot * TauSqr * TauSqr * TauSqr;
      if (FileOut != NULL)
        fprintf(FileOut, "tau1 %11.7f tau3 %11.7f u %14.7f %14.7f\n",
                          Tau1, Tau3, u, uDot);
    }
    else
    {
      /* --- Now use exact method to find f and g --- */
      theta  = R1.Angle(R2);
      theta1 = R2.Angle(R3);

      f1 = 1.0 - ((R1.Mag() * (1.0 - cos(theta)) / p));
      g1 = (R1.Mag()*R2.Mag() * sin(-theta) ) / sqrt(p); // - ANGLE because backwards!!
      f3 = 1.0 - ((R3.Mag() * cos(1.0 - cos(theta1))/p));
      g3 = ( R3.Mag()*R2.Mag() * sin(theta1)) / sqrt(p);
      if (FileOut != NULL)
        fprintf(FileOut, "f1n%11.7f %11.7f f3 %11.7f %11.7f\n", f1, g1, f3, g3);
      c1 =  g3 / (f1 * g3 - f3 * g1);
      c3 = -g1 / (f1 * g3 - f3 * g1);
    }
    /* ---- Solve for all three ranges via matrix equation ---- */
    CMat.Set(-c1, 1, 1);
    CMat.Set(1.0, 2, 1);
    CMat.Set(-c3, 3, 1);
    RhoMat = LIR * CMat;
    if (FileOut != NULL)
    {
      fprintf(FileOut, "c1, c3 %11.7f %11.7f\n", c1, c3);
      FilePrint(RhoMat, " Rho Matrix ", 3, FileOut);
    }
    /* ---- Check for convergence ---- */
  }
  /* ---- Find all three vectors ri ---- */
  for (UINT i = 1; i <= 3; i++)
  {
    R1.Set( RhoMat.Get(1, 1) * L1.Get(i) / c1 + RS1.Get(i), i);
    R2.Set(-RhoMat.Get(2, 1) * L2.Get(i)      + RS2.Get(i), i);
    R3.Set( RhoMat.Get(3, 1) * L3.Get(i) / c3 + RS3.Get(i), i);
  }
  if (FileOut != NULL)
  {
    fprintf(FileOut, "r1 %11.7f %11.7f %11.7f %11.7f",
                      R1.Get(1), R1.Get(2), R1.Get(3), R1.Mag());
    fprintf(FileOut, "r2 %11.7f %11.7f %11.7f %11.7f",
                      R2.Get(1), R2.Get(2), R2.Get(3), R2.Mag());
    fprintf(FileOut, "r3 %11.7f %11.7f %11.7f %11.7f",
                      R3.Get(1), R3.Get(2), R3.Get(3), R3.Mag());
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE ANGLESLAPLACE
|
|  This PROCEDURE solves the problem of orbit determination using three
|    optical sightings and the method of Laplace.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    TRtAsc1      - Right Ascension #1            rad
|    TRtAsc2      - Right Ascension #2            rad
|    TRtAsc3      - Right Ascension #3            rad
|    TDecl1       - Declination #1                rad
|    TDecl2       - Declination #2                rad
|    TDecl3       - Declination #3                rad
|    JD1          - Julian Date of 1st sighting   Days from 4713 BC
|    JD2          - Julian Date of 2nd sighting   Days from 4713 BC
|    JD3          - Julian Date of 3rd sighting   Days from 4713 BC
|    RS1          - IJK SITE position vector #1   ER
|    RS2          - IJK SITE position vector #2   ER
|    RS3          - IJK SITE position vector #3   ER
|
|  OutPuts        :
|    R            - IJK position vector           ER
|    V            - IJK velocity vector           ER / TU
|
|  Locals         :
|    L1           - Line of SIGHT vector for 1st
|    L2           - Line of SIGHT vector for 2nd
|    L3           - Line of SIGHT vector for 3rd
|    LDot         - 1st derivative of L2
|    LDDot        - 2nd derivative of L2
|    RS2Dot       - 1st Derivative of RS2 - vel
|    RS2DDot      - 2nd Derivative of RS2
|    t12t13       - (t1-t2) * (t1-t3)
|    t21t23       - (t2-t1) * (t2-t3)
|    t31t32       - (t3-t1) * (t3-t2)
|    i            - index
|    D            -
|    D1           -
|    D2           -
|    D3           -
|    D4           -
|    OldR         - Previous iteration on r
|    Rho          - Range from SITE to satellite at t2
|    RhoDot       -
|    DMat         -
|    D1Mat        -
|    D2Mat        -
|    D3Mat        -
|    D4Mat        -
|    EarthRate    - Angular rotation of the earth
|    L2DotRS      - Vector L2 Dotted with RSijk
|    Temp         - Temporary vector
|    Temp1        - Temporary vector
|    Small        - Tolerance
|    Roots        -
|
|  Coupling       :
|    MAG          - Magnitude of a vector
|    DETERMINANT  - Evaluate the determinant of a matrix
|    CROSS        - CROSS product of two vectors
|    NORM         - Normlize a matrix
|    ASSIGNVAL    - Assign a value to a matrix
|    GETVAL       - Get a value from a matrix
|    INITMATRIX   - Initialize a matrix and fil with 0.0's

|    DELMATRIX    - Delete a matrix
|    FACTOR       - Find the roots of a polynomial
|
|  References     :
|    Vallado       2001, 413-417
|
 -----------------------------------------------------------------------------*/
void AnglesLaplace
    (
      Real TDecl1, Real TDecl2, Real TDecl3, 
      Real TRtAsc1, Real TRtAsc2, Real TRtAsc3, 
      Real JD1, Real JD2, Real JD3,
      Vector RS1, Vector RS2, Vector RS3, Vector& r2, Vector& v2
    )
{
  const Real OmegaEarth = 0.05883359980154919;
  const Real TUDay      = 0.00933809017716;   // TUDay:= 58.132440906;
  const Real Mu         = 1.0;
  const Real Small      = 0.0000001;

  Real Poly[16];
  Real Roots[15][2];
  Matrix DMat(3, 3), DMat1(3, 3), DMat2(3, 3), DMat3(3, 3), DMat4(3, 3);
  Vector L1(3), L2(3), L3(3), LDot(3), LDDot(3), RS2Dot(3), RS2DDot(3), 
         EarthRate(3), Temp(3), Temp1(3);
  Real   D, D1, D2, D3, D4, Rho, RhoDot, t1t13, t1t3, t31t3, TempRoot,
         tau1, tau3, BigR2, L2DotRS;
  char   tc, chg;
  Byte   ScHi;

  /* ----------------------   Initialize   ------------------------ */
  EarthRate.Set(0.0, 1);
  EarthRate.Set(0.0, 2);
  EarthRate.Set(OmegaEarth, 3);

  JD1 = JD1 / TUDay;   // days to TU
  JD2 = JD2 / TUDay;   // days to TU
  JD3 = JD3 / TUDay;   // days to TU

  /* ---- set middle to 0, deltas to other times ---- */
  tau1 = JD1 - JD2;
  tau3 = JD3 - JD2;

  /* --------------- Find Line of SIGHT vectors ------------------- */
  L1.Set(cos(TDecl1) * cos(TRtAsc1), 1);
  L2.Set(cos(TDecl1) * sin(TRtAsc1), 2);
  L2.Set(sin(TDecl1), 3);
  L1.Set(cos(TDecl2) * cos(TRtAsc2), 1);
  L2.Set(cos(TDecl2) * sin(TRtAsc2), 2);
  L2.Set(sin(TDecl2), 3);
  L1.Set(cos(TDecl3) * cos(TRtAsc3), 1);
  L2.Set(cos(TDecl3) * sin(TRtAsc3), 2);
  L2.Set(sin(TDecl3), 3);

  /* -------------------------------------------------------------- */
  // Using Lagrange Interpolation formula to derive an expression
  // for L(t), substitute t=t2 and differentiate to obtain the
  // derivatives of L.
  /* --------------------------------------------------------------- */
  t1t13 = 1.0 / (tau1 * (tau1 - tau3));
  t1t3  = 1.0 / (tau1 * tau3);
  t31t3 = 1.0 / ((tau3 - tau1) * tau3);
  for (UINT i = 1; i <= 3; i++)
  {
    LDot.Set((-tau3 * t1t13) * L1.Get(i) + 
             ((-tau1 - tau3) * t1t3) * L2.Get(i) + 
             (-tau1 * t31t3) * L3.Get(i), i);
    LDDot.Set((2.0 * t1t13) * L1.Get(i) +
              (2.0 * t1t3)  * L2.Get(i) +
              (2.0 * t31t3) * L3.Get(i), i);
  }

  /* -------------------- Find 2nd derivative of RSijk --------------- */
  Temp  = RS1.Cross(RS2);
  Temp1 = RS2.Cross(RS3);

  // needs a different test xxxx!!
  if ((fabs(Temp.Mag()) > Small) && (fabs(Temp1.Mag()) > Small))
  {
    /* ------- All sightings from one SITE --------- */
    // fix this testhere
    for (UINT i = 1; i <= 3; i++)
    {
      // eSc pg 268  doesn't seem to work!!! xx
      RS2Dot.Set((-tau3 * t1t13) * RS1.Get(i) +
                 ((-tau1 - tau3) * t1t3) * RS2.Get(i) +
                 (-tau1 * t31t3) * RS3.Get(i), i);
      RS2DDot.Set((2.0 * t1t13) * RS1.Get(i) +
                  (2.0 * t1t3 ) * RS2.Get(i) +
                  (2.0 * t31t3) * RS3.Get(i), i);
    }

    RS2Dot = EarthRate.Cross(RS2);
    RS2DDot = EarthRate.Cross(RS2Dot);
  }
  else
  {
    /* ---- Each sighting from a different SITE ---- */
    for (UINT i = 1; i <= 3; i++)
    {
      RS2Dot.Set((-tau3 * t1t13) * RS1.Get(i) +
                 ((-tau1 - tau3) * t1t3) * RS2.Get(i) +
                 (-tau1 * t31t3) * RS3.Get(i), i);
      RS2DDot.Set((2.0 * t1t13) * RS1.Get(i) +
                  (2.0 * t1t3 ) * RS2.Get(i) +
                  (2.0 * t31t3) * RS3.Get(i), i);
    }
  }
  for (UINT i = 1; i <= 3; i++)
  {
    DMat.Set(2.0 * L2.Get(i), i, 1);
    DMat.Set(2.0 * LDot.Get(i), i, 2);
    DMat.Set(2.0 * LDot.Get(i), i, 3);

    /* ----  Position determinants ---- */
    DMat1.Set(L2.Get(i), i, 1);
    DMat1.Set(LDot.Get(i), i, 2);
    DMat1.Set(RS2DDot.Get(i), i, 3);
    DMat2.Set(L2.Get(i), i, 1);
    DMat2.Set(LDot.Get(i), i, 2);
    DMat2.Set(RS2.Get(i), i, 3);

    /* ----  Velocity determinants ---- */
    DMat3.Set(L2.Get(i), i, 1);
    DMat3.Set(RS2DDot.Get(i), i, 2);
    DMat3.Set(LDDot.Get(i), i, 3);
    DMat4.Set(L2.Get(i), i, 1);
    DMat4.Set(RS2.Get(i), i, 2);
    DMat4.Set(LDDot.Get(i), i, 3);
  }

  D  = DMat.Determinant();
  D1 = DMat1.Determinant();
  D2 = DMat2.Determinant();
  D3 = DMat3.Determinant();
  D4 = DMat4.Determinant();

  /* --------------  Iterate to find Rho magnitude --------------- */
/*
|     r[4]:= 1.5;  { First Guess }
|     WriteLn( 'Input initial guess for r[4] ' );
|     Readln( r[4] );
|     i:= 1;
|     REPEAT
|         OldR:= r[4];
|         Rho:= -2.0*D1/D - 2.0*D2/(r[4]*r[4]*r[4]*D);
|         r[4]:= SQRT( Rho*Rho + 2.0*Rho*L2DotRS + RS2[4]*RS2[4] );
|         INC(i);
|         r[4]:= (OldR - r[4] ) / 2.0;            // Simple bissection }
|         WriteLn( FileOut,'Rho guesses ',i:2,'Rho ',Rho:14:7,' r[4] ',r[4]:14:7,oldr:14:7 );
|// seems to converge, but wrong Numbers 
|         INC(i);
|     UNTIL ( ABS( OldR-R[4] ) < Small ) or ( i >= 30 );
*/

  if (fabs(D) > 0.000001)
  {
    /* ---------------- Solve eighth order poly ----------------- */
    L2DotRS  = L2.Dot(RS2);
    Poly[ 0] =  1.0; // r2^8th variable!!!!!!!!!!!!!!
    Poly[ 1] =  0.0; 
    Poly[ 2] =  (L2DotRS * 4.0 * D1 / D - 4.0 * D1 * D1 / (D * D) - 
               RS2.Mag() * RS2.Mag());
    Poly[ 3] =  0.0; 
    Poly[ 4] =  0.0; 
    Poly[ 5] =  Mu * (L2DotRS * 4.0 * D2 / D - 8.0 * D1 * D2 / (D * D));
    Poly[ 6] =  0.0; 
    Poly[ 7] =  0.0; 
    Poly[ 8] = -4.0 * Mu * D2 * D2 / (D * D);
    Poly[ 9] =  0.0; 
    Poly[10] =  0.0; 
    Poly[11] =  0.0; 
    Poly[12] =  0.0; 
    Poly[13] =  0.0; 
    Poly[14] =  0.0; 
    Poly[15] =  0.0; 
    Factor(Poly, 8, (Real **)Roots);

    /* ---- Find correct (xx) root ---- */
    BigR2 = 0.0;
    for (UINT j = 0; j < 8; j++)
      if (fabs(Roots[j][2]) < Small)
      {
        printf("Root %d %f + %f\n", j+1, Roots[j][1], Roots[j][2]);
        TempRoot = Roots[j][0] * Roots[j][0];
        TempRoot = TempRoot * TempRoot * TempRoot * TempRoot +
                   Poly[2]  * TempRoot * TempRoot * TempRoot + 
                   Poly[5]  * Roots[j][0] * TempRoot + Poly[8];
        if (FileOut != NULL)
          fprintf(FileOut, "Root %d %f + %f j  value = %f\n",
                            j, Roots[j][0], Roots[j][1], TempRoot);
        if (Roots[j][0] > BigR2)
          BigR2 = Roots[j][0];
      }
    printf("BigR2 %14.7f\n", BigR2);
    if (FileOut != NULL)
      fprintf(FileOut, "BigR2 %14.7f\n", BigR2);
    printf("Keep this root ? ");
    scanf("%f\n", &BigR2);
    
    Rho = -2.0 * D1 / D - 2.0 * Mu * D2 / (BigR2 * BigR2 * BigR2 * D);

    /* ---- Find the middle position vector ---- */
    for (UINT k = 1; k <= 3; k++)
      r2.Set(Rho * L2.Get(k) + RS2.Get(k), k);

    /* -------- Find RhoDot magnitude -------- */
    RhoDot = -D3 / D - Mu * D4 / (r2.Mag() * r2.Mag() * r2.Mag() * D);
    if (FileOut != NULL)
    {
      fprintf(FileOut, "Rho %14.7f\n", Rho);
      fprintf(FileOut, "RhoDot %14.7f\n", RhoDot);
    }

    /* ----- Find middle velocity vector ----- */
    for (UINT i = 1; i <= 3; i++)
      v2.Set(RhoDot * L2.Get(i) + Rho * LDot.Get(i) + RS2Dot.Get(i), i);
  }
  else
    printf("Determinant value was zero %f\n", D);

  if (FileOut != NULL)
  {
    fprintf(FileOut, "t123 %18.7f %18.7f %18.7f TU\n", JD1, JD2, JD3);
    fprintf(FileOut, "t123 %18.7f %18.7f %18.7f Days\n", 
                      JD1 * TUDay, JD2 * TUDay, JD3 * TUDay);
    fprintf(FileOut, "tau  %11.7f %11.7f TU\n", tau1, tau3);
    fprintf(FileOut, "tau  %11.7f %11.7f MIN\n", 
                      tau1 * 13.446849, tau3 * 13.446849);
    fprintf(FileOut, "delta123 %12.6f %12.6f %12.6f\n",
                      TDecl1 * 57.2957, TDecl2 * 57.2957, TDecl3 * 57.2957);
    fprintf(FileOut, "RtAsc123 %12.6f %12.6f %12.6f\n",
                     TRtAsc1 * 57.2957, TRtAsc2 * 57.2957, TRtAsc3 * 57.2957);
    fprintf(FileOut, "RtAsc1   %12.6f %12.6f %12.6f\n",
                      RS1.Get(1), RS1.Get(2), RS1.Get(3));
    fprintf(FileOut, "RtAsc 2  %12.6f %12.6f %12.6f\n",
                      RS2.Get(1), RS2.Get(2), RS2.Get(3));
    fprintf(FileOut, "RtAsc  3 %12.6f %12.6f %12.6f\n",
                      RS3.Get(1), RS3.Get(2), RS3.Get(3));
    fprintf(FileOut, "los 1    %12.6f %12.6f %12.6f\n",
                      L1.Get(1), L1.Get(2), L1.Get(3));
    fprintf(FileOut, "los 2    %12.6f %12.6f %12.6f\n",
                      L2.Get(1), L2.Get(2), L2.Get(3));
    fprintf(FileOut, "los 3    %12.6f %12.6f %12.6f\n",
                      L3.Get(1), L3.Get(2), L3.Get(3));
    fprintf(FileOut, "LDot     %12.6f %12.6f %12.6f\n",
                      LDot.Get(1), LDot.Get(2), LDot.Get(3));
    fprintf(FileOut, "LDDot    %12.6f %12.6f %12.6f\n",
                      LDDot.Get(1), LDDot.Get(2), LDDot.Get(3));
    fprintf(FileOut, "RS2     %12.6f %12.6f %12.6f\n",
                      RS2.Get(1), RS2.Get(2), RS2.Get(3));
    fprintf(FileOut, "RS2Dot  %12.6f %12.6f %12.6f\n",
                      RS2Dot.Get(1), RS2Dot.Get(2), RS2Dot.Get(3));
    fprintf(FileOut, "RS2DDot %12.6f %12.6f %12.6f\n",
                      RS2DDot.Get(1), RS2DDot.Get(2), RS2DDot.Get(3));
    fprintf(FileOut, "D 01234 = %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                      D, D1, D2, D3, D4);
  }
  DMat.Display(" D Matrix ", 6);
  DMat1.Display(" D1 Matrix ", 6);
  DMat2.Display(" D2 Matrix ", 6);
  DMat3.Display(" D3 Matrix ", 6);
  DMat4.Display(" D4 Matrix ", 6);
}

/* -------------------- Conversion techniques --------------------- */
/*------------------------------------------------------------------------------
|
|                           PROCEDURE RADEC_AZEL
|
| This PROCEDURE converts right ascension declination values with
|   azimuth, and elevation.  Notice the range is not defined because
|   Right ascension declination only allows a unit vector to be formed.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RtAsc       - Right Ascension                0.0 to 2Pi rad
|    Decl        - Declination                    -Pi/2 to Pi/2 rad
|    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
|    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    Az          - Azimuth                        0.0 to 2Pi rad
|    El          - Elevation                      -Pi/2 to Pi/2 rad
|
|  Locals        :
|    LHA         - Local Hour ANGLE               -2Pi to 2Pi rad
|    Sinv        - Sine value
|    Cosv        - Cosine value
|
|  Coupling      :
|    ARCSIN      - Arc sine FUNCTION
|    ATAN2       - Arc Tangent FUNCTION that resolves quadrant ambiguites
|
|  References    :
|    Vallado       2001, 255-257, Alg 28
|
 -----------------------------------------------------------------------------*/
void RaDec_AzEl
    (
      Real& RtAsc, Real& Decl, Real& LST, Real& Latgd, Direction dir, 
      Real& Az, Real& El
    )
{
  Real Sinv, Cosv, LHA;
  if (dir == FROM)
  {
    Decl  = asin(sin(El) * sin(Latgd) + cos(El) * cos(Latgd) * cos(Az));

    Sinv  = -(sin(Az) * cos(El) * cos(Latgd)) / (cos(Latgd) * cos(Decl));
    Cosv  = (sin(El) - sin(Latgd) * sin(Decl)) / (cos(Latgd) * cos(Decl));
    LHA   = Atan2(Sinv, Cosv);
    RtAsc = LST - LHA;
  }
  else
  {
    LHA  = LST - RtAsc;
    El   = asin(sin(Decl) * sin(Latgd) + cos(Decl) * cos(Latgd) * cos(LHA));
    Sinv = -sin(LHA) * cos(Decl) * cos(Latgd) / (cos(El) * cos(Latgd));
    Cosv = (sin(Decl) - sin(El) * sin(Latgd)) / (cos(El) * cos(Latgd));
    Az   = Atan2(Sinv, Cosv);
  }

  if (Show == 'Y')
    if (FileOut != NULL)
      fprintf(FileOut, "%f\n", LHA * 180.0 / PI);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RADEC_ELATLON
|
|  This PROCEDURE converts right-ascension declination values with ecliptic
|    latitude and longitude values.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RtAsc       - Right Ascension                rad
|    Decl        - Declination                    rad
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    EclLat      - Ecliptic Latitude              -Pi/2 to Pi/2 rad
|    EclLon      - Ecliptic Longitude             -Pi/2 to Pi/2 rad
|
|  Locals        :
|    Obliquity   - Obliquity of the ecliptic      rad
|    Sinv        -
|    Cosv        -
|
|  Coupling      :
|    ARCSIN      - Arc sine FUNCTION
|    ATAN2       - Arc tangent FUNCTION that resolves quadrant ambiguites
|
|  References    :
|    Vallado       2001, 259, Eq 4-19, Eq 4-20
|
 -----------------------------------------------------------------------------*/
void RaDec_ELatLon
    (
      Real& RtAsc, Real& Decl, Direction dir, Real& EclLat, Real& EclLon
    )
{
  Real Sinv, Cosv, Obliquity;

  Obliquity = 0.40909280; // 23.439291/rad
  if (dir == FROM)
  {
    Decl  = asin(sin(EclLat) * cos(Obliquity) + 
            cos(EclLat) * sin(Obliquity) * sin(EclLon));
    Sinv  = (-sin(EclLat) * sin(Obliquity) +
            cos(EclLat) * cos(Obliquity) * sin(EclLon)) / cos(Decl);
    Cosv  = cos(EclLat) * cos(EclLon) / cos(Decl);
    RtAsc = Atan2(Sinv, Cosv);
  }
  else
  {
    EclLat = asin(-cos(Decl) * sin(RtAsc) * sin(Obliquity) +
             sin(Decl) * cos(Obliquity));
    Sinv   = (cos(Decl) * sin(RtAsc) * cos(Obliquity) +
             sin(Decl) * sin(Obliquity)) / cos(EclLat);
    Cosv   = cos(Decl) * cos(RtAsc) / cos(EclLat);
    EclLon = Atan2(Sinv, Cosv);
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RV_ELATLON
|
|  This PROCEDURE converts ecliptic latitude and longitude with position and
|    velocity vectors. Uses velocity vector to find the solution of singular
|    cases.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rijk        - IJK position vector            ER
|    Vijk        - IJK velocity vector            ER/TU
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    Rr          - Radius of the sat              ER
|    EclLat      - Ecliptic Latitude              -Pi/2 to Pi/2 rad
|    EclLon      - Ecliptic Longitude             -Pi/2 to Pi/2 rad
|    DRr         - Radius of the sat rate         ER/TU
|    DEclLat     - Ecliptic Latitude rate         -Pi/2 to Pi/2 rad
|    EEclLon     - Ecliptic Longitude rate        -Pi/2 to Pi/2 rad
|
|  Locals        :
|    Obliquity   - Obliquity of the ecliptic      rad
|    Temp        -
|    Temp1       -
|    Re          - Position vec in eclitpic frame
|    Ve          - Velocity vec in ecliptic frame
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ROT1        - Rotation about 1st axis
|    DOT         - DOT product
|    ARCSIN      - Arc sine FUNCTION
|    ATAN2       - Arc tangent FUNCTION that resolves quadrant ambiguites
|
|  References    :
|    Vallado       2001, 257-259, Eq 4-15
|
 -----------------------------------------------------------------------------*/
void RV_ELatLon
    (
      Vector& Rijk, Vector& Vijk, Direction dir,
      Real& rr, Real& EclLat, Real& EclLon, 
      Real& DRr, Real& DEclLat, Real& DEclLon
    )
{
  const Real Small = 0.00000001;

  Vector Re(3), Ve(3);
  Real   Obliquity, Temp, Temp1;

  Obliquity = 0.40909280; // 23.439291/rad
  if (dir == FROM)
  {
    Re.Set(rr * cos(EclLat) * cos(EclLon), 1);
    Re.Set(rr * cos(EclLat) * sin(EclLon), 2);
    Re.Set(rr               * sin(EclLon), 3);
    Ve.Set(DRr * cos(EclLat) * cos(EclLon) -
            rr * sin(EclLat) * cos(EclLon) * DEclLat -
            rr * cos(EclLat) * sin(EclLon) * DEclLon, 1);
    Ve.Set(DRr * cos(EclLat) * sin(EclLon) -
            rr * sin(EclLat) * cos(EclLon) * DEclLat +
            rr * cos(EclLat) * cos(EclLon) * DEclLon, 2);
    Ve.Set(DRr * sin(EclLat) + rr * cos(EclLat) * DEclLat, 3);

    Rijk = Re.Rot1(-Obliquity);
    Vijk = Ve.Rot1(-Obliquity);
  }
  else
  {
    /* -------------- Calculate Angles and Rates ---------------- */
    rr = Re.Mag();
    Temp = sqrt(Re.Get(1) * Re.Get(1) + Re.Get(2) * Re.Get(2));
    if (Temp < Small)
    {
      Temp1 = sqrt(Ve.Get(1) * Ve.Get(1) + Ve.Get(2) * Ve.Get(2));
      if (fabs(Temp1) > Small)
        EclLon = Atan2(Ve.Get(2) / Temp1, Ve.Get(1) / Temp1);
      else
        EclLon = 0.0;
    }
    else
      EclLon = Atan2(Re.Get(2) / Temp, Re.Get(1) / Temp);
    EclLat = asin(Re.Get(3) / Re.Mag());

    Temp1 = -Re.Get(2) * Re.Get(2) - Re.Get(1) * Re.Get(1); // different now
    DRr   = Re.Dot(Ve) / rr;
    if (fabs(Temp1) > Small)
      DEclLon = (Ve.Get(1) * Re.Get(2) - Ve.Get(2) * Re.Get(1)) / Temp1;
    else
      DEclLon = 0.0;
    if (fabs(Temp) > Small)
      DEclLat = (Ve.Get(3) - DRr * sin(EclLat)) / Temp;
    else
      DEclLat = 0.0;
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RV_RADEC
|
|  This PROCEDURE converts the right ascension and declination values with
|    position and velocity vectors of a satellite. Uses velocity vector to
|    find the solution of singular cases.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rijk        - IJK position vector            ER
|    Vijk        - IJK velocity vector            ER/TU
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    Rr          - Radius of the satellite        ER
|    RtAsc       - Right Ascension                rad
|    Decl        - Declination                    rad
|    DRr         - Radius of the satellite rate   ER/TU
|    DRtAsc      - Right Ascension rate           rad/TU
|    DDecl       - Declination rate               rad/TU
|
|  Locals        :
|    Temp        - Temporary position vector
|    Temp1       - Temporary variable
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ATAN2       - ARCTAN FUNCTION that resolves the quadrant ambiguities
|    DOT         - DOT product of two vectors
|    ARCSIN      - Arc sine FUNCTION
|
|  References    :
|    Vallado       2001, 246-248, Alg 25
|
 -----------------------------------------------------------------------------*/
void RV_RaDec
    (
      Vector& Rijk, Vector& Vijk, Direction dir,
      Real& rr, Real& RtAsc, Real& Decl, Real& DRr, Real& DRtAsc, Real& DDecl
    )
{
  const Real Small = 0.00000001;

  Real Temp, Temp1;

  if (dir == FROM)
  {
    Rijk.Set(rr * cos(Decl) * cos(RtAsc), 1);
    Rijk.Set(rr * cos(Decl) * sin(RtAsc), 2);
    Rijk.Set(rr * sin(Decl), 3);
    Vijk.Set(DRr * cos(Decl) * cos(RtAsc) - 
              rr * sin(Decl) * cos(RtAsc) * DDecl -
              rr * cos(Decl) * sin(RtAsc) * DRtAsc, 1);
    Vijk.Set(DRr * cos(Decl) * sin(RtAsc) - 
              rr * sin(Decl) * sin(RtAsc) * DDecl +
              rr * cos(Decl) * cos(RtAsc) * DRtAsc, 2);
    Vijk.Set(DRr * sin(Decl) + rr * cos(Decl) * DDecl, 3);
  }
  else
  {
    /* -------------- Calculate Angles and Rates ---------------- */
    rr   = Rijk.Mag();
    Temp = sqrt(Rijk.Get(1) * Rijk.Get(1) + Rijk.Get(2) * Rijk.Get(2));
    if (Temp < Small)
    {
      Temp1 = sqrt(Vijk.Get(1) * Vijk.Get(1) + Vijk.Get(2) * Vijk.Get(2));
      if (fabs(Temp1) > Small)
        RtAsc = Atan2(Vijk.Get(2) / Temp1, Vijk.Get(1) / Temp1);
      else
        RtAsc = 0.0;
    }
    else
      RtAsc = Atan2(Rijk.Get(2) / Temp, Rijk.Get(1) / Temp);
    Decl    = asin(Rijk.Get(3) / Rijk.Mag());

    Temp1 = -Rijk.Get(2) * Rijk.Get(2) - Rijk.Get(1) * Rijk.Get(1);
    DRr   = Rijk.Dot(Vijk) / rr;
    if (fabs(Temp1) > Small)
      DRtAsc = (Vijk.Get(1) * Rijk.Get(2) - Vijk.Get(2) * Rijk.Get(1)) / Temp1;
    else
      DRtAsc = 0.0;
    if (fabs(Temp) > Small)
      DDecl = (Vijk.Get(3) - DRr * sin(Decl)) / Temp;
    else
      DDecl = 0.0;
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RV_RAZEL
|
|  This PROCEDURE converts Range, Azimuth, and Elevation and their rates with
|    the Geocentric Equatorial (IJK) Position and Velocity vectors.  Notice the
|    value of small as it can affect rate term calculations. Uses velocity
|    vector to find the solution of singular cases.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rijk        - IJK position vector            ER
|    Vijk        - IJK velocity vector            ER/TU
|    RSijk       - IJK SITE Position Vector       ER
|    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
|    LST         - Local SIDEREAL Time            -2Pi to Pi rad
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    Rho         - Satellite Range from SITE      ER
|    Az          - Azimuth                        0.0 to 2Pi rad
|    El          - Elevation                      -Pi/2 to Pi/2 rad
|    DRho        - Range Rate                     ER / TU
|    DAz         - Azimuth Rate                   rad / TU
|    DEl         - Elevation rate                 rad / TU
|
|  Locals        :
|    RhoVijk     - IJK Range Vector from SITE     ER
|    DRhoVijk    - IJK Velocity Vector from SITE  ER / TU
|    Rhosez      - SEZ Range vector from SITE     ER
|    DRhosez     - SEZ Velocity vector from SITE  ER
|    WCrossR     - CROSS product result           ER / TU
|    EarthRate   - IJK Earth's rotation rate vec  rad / TU
|    TempVec     - Temporary vector
|    Temp        - Temporary EXTENDED value
|    Temp1       - Temporary EXTENDED value
|    i           - Index
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ADDVEC      - Add two vectors
|    CROSS       - CROSS product of two vectors
|    ROT3        - Rotation about the 3rd axis
|    ROT2        - Rotation about the 2nd axis
|    ATAN2       - Arc tangent FUNCTION which also resloves quadrants
|    DOT         - DOT product of two vectors
|    RVSEZ_RAZEL - Find R and V from SITE in Topocentric Horizon (SEZ) system
|    LNCOM2      - Combine two vectors and constants
|    ARCSIN      - Arc sine FUNCTION
|    SGN         - Returns the sign of a variable
|
|  References    :
|    Vallado       2001, 250-255, Alg 27
|
 -----------------------------------------------------------------------------*/
void RV_RAzEl
    (
      Vector& Rijk, Vector& Vijk, Vector& RSijk, Real Latgd, Real LST,
      Direction dir, 
      Real& Rho, Real& Az, Real& El, Real& DRho, Real& DAz, Real& DEl
    )
{
  const Real HalfPi = PI / 2.0;
  const Real OmegaEarth = 0.05883359980154919;
  const Real Small = 0.0000001;
 
  Real Temp, Temp1;
  Vector Rhoijk(3), DRhoijk(3), Rhosez(3), DRhosez(3), 
         WCrossR(3), EarthRate(3), TempVec(3);

  /* --------------------  Initialize values   -------------------- */
  EarthRate.Set(0.0, 1);
  EarthRate.Set(0.0, 2);
  EarthRate.Set(OmegaEarth, 3);

  if (dir == FROM)
  {
    /* ---------  Find SEZ range and velocity vectors ----------- */
    RVSez_RAzEl(Rhosez, DRhosez, FROM, Rho, Az, El, DRho, DAz, DEl);

    /* ----------  Perform SEZ to IJK transformation ------------ */
    TempVec = Rhosez.Rot2(Latgd - HalfPi);
    Rhoijk  = TempVec.Rot3(-LST);
    TempVec = DRhosez.Rot2(Latgd - HalfPi);
    DRhoijk = TempVec.Rot3(-LST);

    /* ---------  Find IJK range and velocity vectors -----------*/
    Rijk    = Rhoijk + RSijk;
    WCrossR = EarthRate.Cross(Rijk);
    Vijk    = DRhoijk + WCrossR;

    if (Show == 'Y')
      if (FileOut != NULL)
      {
        fprintf(FileOut, "rseb %18.7f %18.7f %18.7f %18.7f ER\n",
                Rhosez.Get(1), Rhosez.Get(2), Rhosez.Get(3), Rhosez.Mag());
        fprintf(FileOut, "vseb %18.7f %18.7f %18.7f %18.7f\n",
                DRhosez.Get(1), DRhosez.Get(2), DRhosez.Get(3), DRhosez.Mag());
        fprintf(FileOut, "RhoV ijk %117.f %11.7f %11.7f\n",
                          Rhoijk.Get(1), Rhoijk.Get(2), Rhoijk.Get(3));
        fprintf(FileOut, "DRhoV ijk %117.f %11.7f %11.7f\n",
                          DRhoijk.Get(1), DRhoijk.Get(2), DRhoijk.Get(3));
        fprintf(FileOut, " Mat r1 %11.7f %11.7f %11.7f\n",
                sin(Latgd) * cos(LST), -sin(LST), cos(Latgd) * cos(LST));
        fprintf(FileOut, " Mat r2 %11.7f %11.7f %11.7f\n",
                sin(Latgd) * sin(LST), cos(LST), cos(Latgd) * sin(LST));
        fprintf(FileOut, " Mat r3 %11.7f %11s %11.7f\n",
                          -cos(Latgd), "0.00", sin(Latgd));
      }
  }
  else
  {
    /* ------- Find IJK range vector from SITE to satellite ----- */
    WCrossR = EarthRate.Cross(Rijk);
    Rhoijk  = Rijk - RSijk;
    DRhoijk = Vijk - WCrossR;
    Rho     = Rhoijk.Mag();

    /* ------------ Convert to SEZ for calculations ------------- */
    TempVec = Rhoijk.Rot3(LST);
    Rhosez  = TempVec.Rot2(HalfPi - Latgd);
    TempVec = DRhoijk.Rot3(LST);
    DRhosez = TempVec.Rot2(HalfPi - Latgd);

    /* ------------ Calculate Azimuth and Elevation ------------- */
    Temp = sqrt(Rhosez.Get(1) * Rhosez.Get(1) + Rhosez.Get(2) * Rhosez.Get(2));
    if (fabs(Rhosez.Get(2)) < Small)
      if (Temp < Small)
      {
        Temp1 = sqrt(DRhosez.Get(1) * DRhosez.Get(1) + 
                     DRhosez.Get(2) * DRhosez.Get(2));
        Az    = Atan2(DRhosez.Get(2) / Temp1, -DRhosez.Get(1) / Temp1);
      }
      else
        if (Rhosez.Get(1) > 0.0)
          Az = PI;
        else
          Az = 0.0;
    else
      Az = Atan2(Rhosez.Get(2) / Temp, -Rhosez.Get(1) / Temp);

    if (Temp < Small)  // directly over the north pole
      El = Sgn(Rhosez.Get(3)) * HalfPi; // +- 90
    else
      El = asin(Rhosez.Get(3) / Rhosez.Mag());

    /* ----- Calculate Range, Azimuth and Elevation rates ------- */
    DRho = Rhosez.Dot(DRhosez) / Rho;
    if (fabs(Temp * Temp) > Small)
      DAz = (DRhosez.Get(1) * Rhosez.Get(2) - DRhosez.Get(2) * Rhosez.Get(1)) /
            (Temp * Temp);
    else
      DAz = 0.0;

    if (fabs(Temp) > 0.00000001)
      DEl = (DRhosez.Get(3) - DRho * sin(El)) / Temp;
    else
      DEl = 0.0;

    if (Show == 'Y')
      if (FileOut != NULL)
      {
        fprintf(FileOut, "rsez %18.7f %18.7f %18.7f %18.7f ER\n",
                Rhosez.Get(1), Rhosez.Get(2), Rhosez.Get(3), Rhosez.Mag());
        fprintf(FileOut, "vsez %18.7f %18.7f %18.7f %18.7f\n",
                DRhosez.Get(1), DRhosez.Get(2), DRhosez.Get(3), DRhosez.Mag());
      }
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RV_TRADEC
|
|  This PROCEDURE converts topocentric right-ascension declination with
|    position and velocity vectors. Uses velocity vector to find the
|    solution of singular cases.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rijk        - IJK position vector            ER
|    Vijk        - IJK velocity vector            ER/TU
|    RSijk       - IJK SITE position vector       ER
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    Rho         - Top Radius of the sat          ER
|    TRtAsc      - Top Right Ascension            rad
|    TDecl       - Top Declination                rad
|    DRho        - Top Radius of the sat rate     ER/TU
|    TDRtAsc     - Top Right Ascension rate       rad/TU
|    TDDecl      - Top Declination rate           rad/TU
|
|  Locals        :
|    RhoV        - IJK Range Vector from SITE     ER
|    DRhoV       - IJK Velocity Vector from SITE  ER / TU
|    Temp        - Temporary EXTENDED value
|    Temp1       - Temporary EXTENDED value
|    i           - Index
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    ATAN2       - Arc tangent FUNCTION that resolves the quadrant ambiguities
|    ARCSIN      - Arc sine FUNCTION
|    LNCOM2      - Linear combination of 2 vectors
|    ADDVEC      - Add two vectors
|    DOT         - DOT product of two vectors
|
|  References    :
|    Vallado       2001, 248-250, Alg 26
|
 -----------------------------------------------------------------------------*/
void RV_TRaDec
    (
      Vector& Rijk, Vector& Vijk, Vector& RSijk, Direction dir,
      Real& Rho, Real& TRtAsc, Real& TDecl, 
      Real& DRho, Real& DTRtAsc, Real& DTDecl
    )
{
  const Real Small = 0.00000001;
  const Real OmegaEarth = 0.05883359221938136;  // Earth Rot rad/TU

  Vector EarthRate(3), RhoV(3), DRhoV(3), VSijk(3);
  Real   Latgc, Temp, Temp1;

  Latgc = asin(RSijk.Get(3) / RSijk.Mag());
  EarthRate.Set(0.0, 1);
  EarthRate.Set(0.0, 2);
  EarthRate.Set(OmegaEarth, 2);
  VSijk = EarthRate.Cross(RSijk);

  if (dir == FROM)
  {
    /* --------  Calculate Topocentric Vectors ------------------ */
    RhoV.Set(Rho * cos(TDecl) * cos(TRtAsc), 1);
    RhoV.Set(Rho * cos(TDecl) * sin(TRtAsc), 2);
    RhoV.Set(Rho * sin(TDecl), 3);

    DRhoV.Set(DRho * cos(TDecl) * cos(TRtAsc) -
               Rho * sin(TDecl) * cos(TRtAsc) * DTDecl -
               Rho * cos(TDecl) * sin(TRtAsc) * DTRtAsc, 1);
    DRhoV.Set(DRho * cos(TDecl) * sin(TRtAsc) -
               Rho * sin(TDecl) * sin(TRtAsc) * DTDecl +
               Rho * cos(TDecl) * cos(TRtAsc) * DTRtAsc, 2);
    DRhoV.Set(DRho * sin(TDecl) + Rho * cos(TDecl) * DTDecl, 3);

    /* ------ Find IJK range vector from SITE to satellite ------ */
    Rijk = RhoV + RSijk;
    Vijk = DRhoV + cos(Latgc) * VSijk;

    if (Show == 'Y')
      if (FileOut != NULL)
      {
        fprintf(FileOut, "rtb %18.7f %18.7f %18.7f %18.7f ER\n",
                RhoV.Get(1), RhoV.Get(2), RhoV.Get(3), RhoV.Mag());
        fprintf(FileOut, "vtb %18.7f %18.7f %18.7f %18.7f\n",
                DRhoV.Get(1), DRhoV.Get(2), DRhoV.Get(3), DRhoV.Mag());
      }
  }
  else
  {
    /* ------ Find IJK range vector from SITE to satellite ------ */
    RhoV = Rijk - RSijk;
    DRhoV = Vijk - cos(Latgc) * VSijk;

    /* -------- Calculate Topocentric ANGLE and Rate Values ----- */
    Rho = RhoV.Mag();
    Temp = sqrt(RhoV.Get(1) * RhoV.Get(1) + RhoV.Get(2) * RhoV.Get(2));
    if (Temp < Small)
    {
      Temp1 = sqrt(DRhoV.Get(1) * DRhoV.Get(1) + DRhoV.Get(2) * DRhoV.Get(2));
      TRtAsc = Atan2(DRhoV.Get(2) / Temp1, DRhoV.Get(1) / Temp1);
    }
    else
      TRtAsc = Atan2(RhoV.Get(2) / Temp, RhoV.Get(1) / Temp);

    TDecl = asin(RhoV.Get(3) / RhoV.Mag());

    Temp1 = -RhoV.Get(2) * RhoV.Get(2) - RhoV.Get(1) * RhoV.Get(1);
    DRho = RhoV.Dot(DRhoV) / Rho;
    if (fabs(Temp1) > Small)
      DTRtAsc = (DRhoV.Get(1) * RhoV.Get(2) - DRhoV.Get(2) * RhoV.Get(1))/Temp1;
    else
      DTRtAsc = 0.0;
    if (fabs(Temp) > Small)
      DTDecl = (DRhoV.Get(3) - DRho * sin(TDecl)) / Temp;
    else
      DTDecl = 0.0;

    if (Show == 'Y')
      if (FileOut != NULL)
      {
        fprintf(FileOut, "rta %18.7f %18.7f %18.7f %18.7f ER\n",
                RhoV.Get(1), RhoV.Get(3), RhoV.Get(3), RhoV.Mag());
        fprintf(FileOut, "vta %18.7f %18.7f %18.7f %18.7f ER\n",
                DRhoV.Get(1), DRhoV.Get(3), DRhoV.Get(3), DRhoV.Mag());
      }
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RVSEZ_RAZEL
|
|  This PROCEDURE converts range, azimuth, and elevation values with slant
|    range and velocity vectors for a satellite from a radar SITE in the
|    Topocentric Horizon (SEZ) system.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RhoVec      - SEZ Satellite range vector     ER
|    DRhoVec     - SEZ Satellite velocity vector  ER / TU
|    Direction   - Which set of vars to output    FROM  TOO
|
|  OutPuts       :
|    Rho         - Satellite range from SITE      ER
|    Az          - Azimuth                        0.0 to 2Pi rad
|    El          - Elevation                      -Pi/2 to Pi/2 rad
|    DRho        - Range Rate                     ER / TU
|    DAz         - Azimuth Rate                   rad / TU
|    DEl         - Elevation rate                 rad / TU
|
|  Locals        :
|    SinEl       - Variable for SIN( El )
|    CosEl       - Variable for COS( El )
|    SinAz       - Variable for SIN( Az )
|    CosAz       - Variable for COS( Az )
|    Temp        -
|    Temp1       -
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    SGN         - Returns the sign of a variable
|    DOT         - DOT product
|    ARCSIN      - Arc sine FUNCTION
|    ATAN2       - Arc tangent FUNCTION that resolves quadrant ambiguites
|
|  References    :
|    Vallado       2001, 250-251, Eq 4-4, Eq 4-5
|
 -----------------------------------------------------------------------------*/
void RVSez_RAzEl
    (
      Vector& Rhosez, Vector& DRhosez, Direction dir,
      Real& Rho, Real& Az, Real& El, Real& DRho, Real& DAz, Real& DEl
    )
{
  const Real Small = 0.00000001;
  const Real HalfPi = PI / 2.0;

  Real Temp1, Temp, SinEl, CosEl, SinAz, CosAz;

  if (dir == FROM)
  {
    SinEl = sin(El);
    CosEl = cos(El);
    SinAz = sin(Az);
    CosAz = cos(Az);

    /* ----------------- Form SEZ range vector ------------------ */
    Rhosez.Set(-Rho * CosEl * CosAz, 1);
    Rhosez.Set( Rho * CosEl * SinAz, 2);
    Rhosez.Set( Rho * SinEl, 3);

    /* --------------- Form SEZ velocity vector ----------------- */
    DRhosez.Set(-DRho * CosEl * CosAz + 
                 Rhosez.Get(3) * DEl * CosAz + Rhosez.Get(2) * DAz, 1);
    DRhosez.Set(DRho * CosEl * SinAz - 
                 Rhosez.Get(3) * DEl * SinAz - Rhosez.Get(1) * DAz, 2);
    DRhosez.Set(DRho * SinEl + Rho * DEl * CosEl, 3);
  }
  else
  {
    /* ------------ Calculate Azimuth and Elevation ------------- */
    Temp = sqrt(Rhosez.Get(1) * Rhosez.Get(1) + Rhosez.Get(2) * Rhosez.Get(2));
    if (fabs(Rhosez.Get(2)) < Small)
      if (Temp < Small)
      {
        Temp1 = sqrt(DRhosez.Get(1) * DRhosez.Get(1) + 
                     DRhosez.Get(2) * DRhosez.Get(2));
        Az    = Atan2(DRhosez.Get(2) / Temp1, DRhosez.Get(1) / Temp1);
      }
      else
        if (DRhosez.Get(1)  > 0.0)
          Az = PI;
        else
          Az = 0.0;
    else
      Az = Atan2(Rhosez.Get(2) / Temp, Rhosez.Get(1) / Temp);

    if (Temp < Small)   // directly over the north pole
      El = Sgn(Rhosez.Get(3)) * HalfPi;  // +- 90
    else
      El = asin(Rhosez.Get(3) / Rhosez.Mag());

    /* ------  Calculate Range, Azimuth and Elevation rates ----- */
    DRho = Rhosez.Dot(DRhosez) / Rho;
    if (fabs(Temp * Temp) > Small)
      DAz = (DRhosez.Get(1) * Rhosez.Get(2) - DRhosez.Get(2) * Rhosez.Get(1)) /
             (Temp * Temp);
    else
      DAz = 0.0;

    if (fabs(Temp) > Small)
      DEl = (DRhosez.Get(3) - DRho * sin(El)) / Temp;
    else
      DEl = 0.0;
    if (Show == 'Y')
      if (FileOut != NULL)
      {
        fprintf(FileOut, "rsez %18.7f %18.7f %18.7f %18.7f ER\n",
                Rhosez.Get(1), Rhosez.Get(2), Rhosez.Get(3), Rhosez.Mag());
        fprintf(FileOut, "vsez %18.7f %18.7f %18.7f %18.7f ER\n",
                DRhosez.Get(1), DRhosez.Get(2), DRhosez.Get(3), DRhosez.Mag());
      }
  }
}

/* -------------------- Three vector techniques ------------------- */
/*------------------------------------------------------------------------------
|
|                           PROCEDURE GIBBS
|
|  This PROCEDURE performs the GIBBS method of orbit determination.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R1          - IJK Position vector #1         ER
|    R2          - IJK Position vector #2         ER
|    R3          - IJK Position vector #3         ER
|
|  OutPuts       :
|    V2          - IJK Velocity Vector for R2     ER / TU
|    Theta       - ANGLE between vectors          rad
|    Error       - Flag indicating success        'ok',...
|
|  Locals        :
|    tover2      -
|    l           -
|    Small       - Tolerance for roundoff errors
|    r1mr2       - Magnitude of r1 - r2
|    r3mr1       - Magnitude of r3 - r1
|    r2mr3       - Magnitude of r2 - r3
|    p           - P Vector     r2 x r3
|    q           - Q Vector     r3 x r1
|    w           - W Vector     r1 x r2
|    d           - D Vector     p + q + w
|    n           - N Vector (r1)p + (r2)q + (r3)w
|    s           - S Vector
|                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
|    b           - B Vector     d x r2
|    Theta1      - Temp ANGLE between the vectors rad
|    Pn          - P Unit Vector
|    R1N         - R1 Unit Vector
|    dn          - D Unit Vector
|    Nn          - N Unit Vector
|    i           - index
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    CROSS       - CROSS product of two vectors
|    DOT         - DOT product of two vectors
|    ADD3VEC     - Add three vectors
|    LNCOM2      - Multiply two vectors by two constants
|    LNCOM3      - Add three vectors each multiplied by a constant
|    NORM        - Creates a Unit Vector
|    ANGLE       - ANGLE between two vectors
|
|  References    :
|    Vallado       2001, 432-445, Alg 52, Ex 7-5
|
 -----------------------------------------------------------------------------*/
void Gibbs
    (
      Vector R1, Vector R2, Vector R3, 
      Vector& V2, Real& Theta, Real& Theta1, Real& Copa, char* Error
    )
{
  const Real Small = 0.000001;
  
  Real Tover2, L, R1mr2, R3mr1, R2mr3;
  Vector P(3), Q(3), W(3), D(3), N(3), S(3), B(3), Pn(3), R1N(3), Dn(3), Nn(3);

  /* --------------------  Initialize values   -------------------- */
  strcpy(Error, "ok");
  Theta = 0.0;
  Theta1 = 0.0;
  Copa   = 0.0;
  V2.Clear();

  /* ----------------------------------------------------------------
  |  Determine IF the vectors are coplanar.
  ----------------------------------------------------------------- */
  P   = R2.Cross(R3);
  Q   = R3.Cross(R1);
  W   = R1.Cross(R2);
  Pn  = P.Norm();
  R1N = R1.Norm();
  Copa = asin(Pn.Dot(R1N));

  if (fabs(R1N.Dot(Pn)) > 0.017452406)
    strcpy(Error, "not coplanar");

  /* ---------------- or don't continue processing ---------------- */
  D = P + Q + W;
  N = R1.Mag() * P + R2.Mag() * Q + R3.Mag() * W;
  Nn = N.Norm();
  Dn = D.Norm();

  /* ----------------------------------------------------------------
  |  Determine IF the orbit is possible.  Both D and N must be in
  |    the same direction, and non-zero.
  ----------------------------------------------------------------- */
  if ((fabs(D.Mag()) < Small) || (fabs(N.Mag()) < Small) || 
       (Nn.Dot(Dn) < Small))
    strcpy(Error, "impossible");
  else
  {
    Theta  = R1.Angle(R2);
    Theta1 = R2.Angle(R3);

    /* ------------ Perform GIBBS method to find V2 ----------- */
    R1mr2 = R1.Mag() - R2.Mag();
    R3mr1 = R3.Mag() - R1.Mag();
    R2mr3 = R3.Mag() - R3.Mag();
    S     = R1mr2 * R3 + R3mr1 * R2 + R2mr3 * R1;
    B = D.Cross(R2);
    L = 1.0 / sqrt(D.Mag() * N.Mag());
    Tover2 = L / R2.Mag();
    V2 = Tover2 * B + L * S;
  }

  if (((Show == 'Y') || (Show == 'S')) && (strcmp(Error, "ok") == 0))
    if (FileOut != NULL)
    {
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "P vector = ", P.Get(1), P.Get(2), P.Get(3));
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "Q vector = ", Q.Get(1), Q.Get(2), Q.Get(3));
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "W vector = ", W.Get(1), W.Get(2), W.Get(3));
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "N vector = ", N.Get(1), N.Get(2), N.Get(3));
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "D vector = ", D.Get(1), D.Get(2), D.Get(3));
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "S vector = ", S.Get(1), S.Get(2), S.Get(3));
      fprintf(FileOut, "%16s %9.3f %9.3f %9.3f\n",
                        "B vector = ", B.Get(1), B.Get(2), B.Get(3));
    }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HERRGIBBS
|
|  This PROCEDURE implements the Herrick-GIBBS approximation for orbit
|    determination, and finds the middle velocity vector for the 3 given
|    position vectors.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R1          - IJK Position vector #1         ER
|    R2          - IJK Position vector #2         ER
|    R3          - IJK Position vector #3         ER
|    JD1         - Julian Date of 1st sighting    days from 4713 BC
|    JD2         - Julian Date of 2nd sighting    days from 4713 BC
|    JD3         - Julian Date of 3rd sighting    days from 4713 BC
|
|  OutPuts       :
|    V2          - IJK Velocity Vector for R2     ER / TU
|    Theta       - ANGLE between vectors          rad
|    Error       - Flag indicating success        'ok',...
|
|  Locals        :
|    Dt21        - time delta between r1 and r2   TU
|    Dt31        - time delta between r3 and r1   TU
|    Dt32        - time delta between r3 and r2   TU
|    p           - P vector    r2 x r3
|    Pn          - P Unit Vector
|    R1N         - R1 Unit Vector
|    Theta1      - temporary ANGLE between vec    rad
|    TolAngle    - Tolerance ANGLE  (1 deg)       rad
|    Term1       - 1st Term for HGibbs expansion
|    Term2       - 2nd Term for HGibbs expansion
|    Term3       - 3rd Term for HGibbs expansion
|    i           - Index
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|    CROSS       - CROSS product of two vectors
|    DOT         - DOT product of two vectors
|    ARCSIN      - Arc sine FUNCTION
|    NORM        - Creates a Unit Vector
|    LNCOM3      - Combination of three scalars and three vectors
|    ANGLE       - ANGLE between two vectors
|
|  References    :
|    Vallado       2001, 439-445, Alg 52, Ex 7-4
|
 -----------------------------------------------------------------------------*/
void HerrGibbs
    (
      Vector R1, Vector R2, Vector R3, Real JD1, Real JD2, Real JD3,
      Vector& V2, Real& Theta, Real& Theta1, Real& Copa, char* Error
    )
{
  const Real TUDay = 0.00933809017716;

  Vector P, Pn, R1n;
  Real   Dt21, Dt31, Dt32, Term1, Term2, Term3, TolAngle;

  /* --------------------  Initialize values   -------------------- */
  strcpy(Error, "ok");
  Theta  = 0.0;
  Theta1 = 0.0;
  V2.Clear();
  Dt21 = (JD2 - JD1) / TUDay;
  Dt31 = (JD3 - JD1) / TUDay;   // differences in times
  Dt32 = (JD3 - JD2) / TUDay;

  /* ----------------------------------------------------------------
  |  Determine IF the vectors are coplanar.
  ---------------------------------------------------------------- */
  P   = R2.Cross(R3);
  Pn  = P.Norm();
  R1n = R1.Norm();
  Copa = asin(Pn.Dot(R1n));
  if (fabs(Pn.Dot(R1n)) > 0.017452406)
    strcpy(Error, "not coplanar");

  /* ----------------------------------------------------------------
  | Check the size of the angles between the three position vectors.
  |   Herrick GIBBS only gives "reasonable" answers when the
  |   position vectors are reasonably close.  10 deg is only an estimate.
  ---------------------------------------------------------------- */
  Theta  = R1.Angle(R2);
  Theta1 = R2.Angle(R3);
  if ((Theta > TolAngle) || (Theta1 > TolAngle))
    strcpy(Error, "ANGLE > 1ø");

  /* ------------ Perform Herrick-GIBBS method to find V2 --------- */
  Term1 = -Dt32 * 
          (1.0 / (Dt21 * Dt31) + 1.0 / (12.0 * R1.Mag() * R1.Mag() * R1.Mag()));
  Term2 = (Dt32 - Dt21) * 
          (1.0 / (Dt21 * Dt32) + 1.0 / (12.0 * R2.Mag() * R2.Mag() * R2.Mag()));
  Term3 = Dt21 * 
          (1.0 / (Dt32 * Dt31) + 1.0 / (12.0 * R3.Mag() * R3.Mag() * R3.Mag()));
  V2    = Term1 * R1 + Term2 * R2 + Term3 * R3;
}

/* ----------------------- Lambert techniques -------------------- */
/*------------------------------------------------------------------------------
|
|                           PROCEDURE LAMBERBATTIN
|
|  This PROCEDURE solves Lambert's problem using Battins method. The method is
|    developed in Battin (1987).
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Ro          - IJK Position vector 1          ER
|    R           - IJK Position vector 2          ER
|    DM          - direction of motion            'L','S'
|    DtTU        - Time between R1 and R2         TU
|
|  OutPuts       :
|    Vo          - IJK Velocity vector            ER / TU
|    V           - IJK Velocity vector            ER / TU
|    Error       - Error flag                     'ok',...
|
|  Locals        :
|    i           - Index
|    Loops       -
|    u           -
|    b           -
|    Sinv        -
|    Cosv        -
|    rp          -
|    x           -
|    xn          -
|    y           -
|    l           -
|    m           -
|    CosDeltaNu  -
|    SinDeltaNu  -
|    DNu         -
|    a           -
|    Tan2w       -
|    RoR         -
|    h1          -
|    h2          -
|    Tempx       -
|    eps         -
|    denom       -
|    chord       -
|    k2          -
|    s           -
|    f           -
|    g           -
|    fDot        -
|    am          -
|    ae          -
|    be          -
|    tm          -
|    gDot        -
|    arg1        -
|    arg2        -
|    tc          -
|    AlpE        -
|    BetE        -
|    AlpH        -
|    BetH        -
|    DE          -
|    DH          -
|
|  Coupling      :
|    ARCSIN      - Arc sine FUNCTION
|    ARCCOS      - Arc cosine FUNCTION
|    MAG         - Magnitude of a vector
|    ARCSINH     - Inverse hyperbolic sine
|    ARCCOSH     - Inverse hyperbolic cosine
|    SINH        - Hyperbolic sine
|    POWER       - Raise a base to a POWER
|    ATAN2       - Arc tangent FUNCTION that resolves quadrants
|
|  References    :
|    Vallado       2001, 464-467, Ex 7-5
|
 -----------------------------------------------------------------------------*/
void LambertBattin
    (
      Vector Ro, Vector R, char dm, char OverRev, Real DtTU, 
      Vector& Vo, Vector& V, char* Error
    )
{
  const Real Small = 0.000001;

  Vector RCrossR;
  SINT   i, Loops;
  Real   u, b, Sinv,Cosv, rp, x, xn, y, L, m, CosDeltaNu, SinDeltaNu,DNu, a,
         tan2w, RoR, h1, h2, Tempx, eps, Denom, chord, k2, s, f, g, FDot, am,
         ae, be, tm, GDot, arg1, arg2, tc, AlpE, BetE, AlpH, BetH, DE, DH;

  strcpy(Error, "ok");
  CosDeltaNu = Ro.Dot(R) / (Ro.Mag() * R.Mag());
  RCrossR    = Ro.Cross(R);
  SinDeltaNu = RCrossR.Mag() / (Ro.Mag() * R.Mag());
  DNu        = Atan2(SinDeltaNu, CosDeltaNu);

  RoR   = R.Mag() / Ro.Mag();
  eps   = RoR - 1.0;
  tan2w = 0.25 * eps * eps / (sqrt(RoR) + RoR *(2.0 + sqrt(RoR)));
  rp    = sqrt(Ro.Mag()*R.Mag()) * (Power(cos(DNu * 0.25), 2) + tan2w);

  if (DNu < PI)
    L = (Power(sin(DNu * 0.25), 2) + tan2w ) /
        (Power(sin(DNu * 0.25), 2) + tan2w + cos(DNu * 0.5));
  else
    L = (Power(cos(DNu * 0.25), 2) + tan2w - cos(DNu * 0.5)) /
        (Power(cos(DNu * 0.25), 2) + tan2w);

  m     = DtTU * DtTU / (8.0 * rp * rp * rp);
  xn    = 0.0;   // 0 for par and hyp
  chord = sqrt(Ro.Mag() * Ro.Mag() + R.Mag() * R.Mag() - 
               2.0 * Ro.Mag() * R.Mag() * cos(DNu));
  s     = (Ro.Mag() + R.Mag() + chord) * 0.5;

  Loops = 1;
  while (1 == 1)
  {
    x     = xn;
    Tempx = See(x);
    Denom = 1.0 / ((1.0 + 2.0 * x + L) * (3.0 + x * (1.0 + 4.0 * Tempx)));
    h1    = Power(L + x, 2) * ( 1.0 + (1.0 + 3.0 * x) * Tempx) * Denom;
    h2    = m * ( 1.0 + (x - L) * Tempx) * Denom;

    /* ------------------------ Evaluate CUBIC ------------------ */
    b  = 0.25 * 27.0 * h2 / Power(1.0 + h1, 3);
    u  = -0.5 * b / (1.0 + sqrt(1.0 + b));
    k2 = k(u);

    y  = ((1.0 + h1) / 3.0) * (2.0 + sqrt(1.0 + b) /
         (1.0 - 2.0 * u * k2 * k2));
    xn = sqrt(Power((1.0 - L) * 0.5, 2) + m / (y * y)) - (1.0 + L) * 0.5;

    Loops++;

    if ((fabs(xn - x) < Small) || (Loops > 30))
      break;
  }

  a =  DtTU * DtTU / (16.0 * rp * rp * xn * y * y);
  
  /* -------------------- Find Eccentric anomalies ---------------- */
  /* -------------------------- Hyperbolic ------------------------ */
  if (a < -Small)
  {
    arg1 = sqrt(s / (-2.0 * a));
    arg2 = sqrt((s - chord) / (-2.0 * a));
    /* -------- Evaluate f and g functions -------- */
    AlpH = 2.0 * asinh(arg1);
    BetH = 2.0 * asinh(arg2);
    DH   = AlpH - BetH;
    f    = 1.0 - (a / Ro.Mag()) * (1.0 - cosh(DH));
    GDot = 1.0 - (a / R.Mag()) * (1.0 - cosh(DH));
    FDot = -sqrt(-a) * sinh(DH) / (Ro.Mag() * R.Mag());
  }
  else
  {
    /* ------------------------- Elliptical --------------------- */
    if (a > Small)
    {
      arg1 = sqrt(s / (2.0 * a));
      arg2 = sqrt((s - chord) / (2.0 * a));
      Sinv = arg2;
      Cosv = sqrt(1.0 - (Ro.Mag()+R.Mag() - chord) / (4.0 * a));
      BetE = 2.0 * acos(Cosv);
      BetE = 2.0 * asin(Sinv);
      if (DNu > PI)
        BetE = -BetE;

      Cosv = sqrt(1.0 - s / (2.0 * a));
      Sinv = arg1;

      am = s * 0.5;
      ae = PI;
      be = 2.0 * asin(sqrt((s - chord) / s));
      tm = sqrt(am * am * am) * (ae - (be - sin(be)));
      if (DtTU > tm)
        AlpE = 2.0 * PI - 2.0 * asin(Sinv);
      else
        AlpE = 2.0 * asin(Sinv);
      DE   = AlpE - BetE;
      f    = 1.0 - (a / Ro.Mag()) * (1.0 - cos(DE));
      GDot = 1.0 - (a / R.Mag()) * (1.0 - cos(DE));
      g    = DtTU - sqrt(a * a * a) * (DE - sin(DE));
      FDot = -sqrt(a) * sin(DE) / (Ro.Mag() * R.Mag());
    }
    else
    {
      /* ------------------------- Parabolic -------------------- */
      arg1 = 0.0;
      arg2 = 0.0;
      strcpy(Error, "a = 0 ");
      if (FileOut != NULL)
        fprintf(FileOut, " a parabolic orbit \n");
    }
  }

  for (UINT i = 1; i <= 3; i++)
  {
    Vo.Set((R.Get(i) - f * Ro.Get(i))/ g, i);
    V.Set((GDot * R.Get(i) - Ro.Get(i))/ g, i);
  }

  if (strcmp(Error, "ok") == 0)
    Testamt = f * GDot - FDot * g;
  else
    Testamt = 2.0;

  if (FileOut != NULL)
    fprintf(FileOut, "%8.5f %3d\n", Testamt, Loops);

  BigT = sqrt(8.0 / (s * s * s)) * DtTU;
}
/*------------------------------------------------------------------------------
|
|                           PROCEDURE LAMBERTUNIV
|
|  This PROCEDURE solves the Lambert problem for orbit determination and returns
|    the velocity vectors at each of two given position vectors.  The solution
|    uses Universal Variables for calculation and a bissection technique for
|    updating psi.
|
|  Algorithm     : Setting the initial bounds:
|                  Using -8Pi and 4Pi2 will allow single rev solutions
|                  Using -4Pi2 and 8Pi2 will allow multi-rev solutions
|                  The farther apart the initial guess, the more iterations
|                    because of the iteration
|                  Inner loop is for special cases. Must be sure to exit both!
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R1          - IJK Position vector 1          ER
|    R2          - IJK Position vector 2          ER
|    DM          - direction of motion            'L','S'
|    DtTU        - Time between R1 and R2         TU
|
|  OutPuts       :
|    V1          - IJK Velocity vector            ER / TU
|    V2          - IJK Velocity vector            ER / TU
|    Error       - Error flag                     'ok', ...
|
|  Locals        :
|    VarA        - Variable of the iteration,
|                  NOT the semi or axis!
|    Y           - Area between position vectors
|    Upper       - Upper bound for Z
|    Lower       - Lower bound for Z
|    CosDeltaNu  - Cosine of true anomaly change  rad
|    F           - f expression
|    G           - g expression
|    GDot        - g DOT expression
|    XOld        - Old Universal Variable X
|    XOldCubed   - XOld cubed
|    ZOld        - Old value of z
|    ZNew        - New value of z
|    C2New       - C2(z) FUNCTION
|    C3New       - C3(z) FUNCTION
|    TimeNew     - New time                       TU
|    Small       - Tolerance for roundoff errors
|    i, j        - index
|
|  Coupling      
|    MAG         - Magnitude of a vector
|    DOT         - DOT product of two vectors
|    FINDC2C3    - Find C2 and C3 functions
|
|  References    :
|    Vallado       2001, 459-464, Alg 55, Ex 7-5
|
 -----------------------------------------------------------------------------*/
void LambertUniv
    (
      Vector Ro, Vector R, char Dm, char OverRev, Real DtTU, 
      Vector& Vo, Vector& V, char* Error)
{
  const Real TwoPi   = 2.0 * PI;
  const Real Small   = 0.0000001;
  const UINT NumIter = 40;

  UINT Loops, i, YNegKtr;
  Real VarA, Y, Upper, Lower, CosDeltaNu, F, G, GDot, XOld, XOldCubed, FDot,
       PsiOld, PsiNew, C2New, C3New, dtNew;

  /* --------------------  Initialize values   -------------------- */
  strcpy(Error, "ok");
  PsiNew = 0.0;
  Vo.Clear();
  V.Clear();

  CosDeltaNu = Ro.Dot(R) / (Ro.Mag() * R.Mag());
  if (Dm == 'L')
    VarA = -sqrt(Ro.Mag() * R.Mag() * (1.0 + CosDeltaNu));
  else
    VarA =  sqrt(Ro.Mag() * R.Mag() * (1.0 + CosDeltaNu));

  /* ----------------  Form Initial guesses   --------------------- */
  PsiOld = 0.0;
  PsiNew = 0.0;
  XOld   = 0.0;
  C2New  = 0.5;
  C3New  = 1.0 / 6.0;

  /* -------- Set up initial bounds for the bissection ------------ */
  if (OverRev == 'N')
  {
    Upper = TwoPi * TwoPi;
    Lower = -4.0 * TwoPi;
  }
  else
  {
    Upper = -0.001 + 4.0 * TwoPi * TwoPi; // at 4, not alw work, 2.0, makes
    Lower =  0.001+TwoPi*TwoPi;           // orbit bigger, how about 2 revs??xx
  }

  /* --------  Determine IF the orbit is possible at all ---------- */
  if (fabs(VarA) > Small)
  {
    Loops   = 0;
    YNegKtr = 1; // y neg ktr
    while (1 == 1)
    {
      if (fabs(C2New) > Small)
        Y = Ro.Mag() + R.Mag() - (VarA * (1.0 - PsiOld * C3New) / sqrt(C2New));
      else
        Y = Ro.Mag() + R.Mag();
      /* ------- Check for negative values of y ------- */
      if ((VarA > 0.0) && (Y < 0.0))
      {
        YNegKtr = 1;
        while (1 == 1)
        {
          PsiNew = 0.8 * (1.0 / C3New) * 
                   (1.0 - (Ro.Mag() + R.Mag()) * sqrt(C2New) / VarA);

          /* ------ Find C2 and C3 functions ------ */
          FindC2C3(PsiNew, C2New, C3New);
          PsiOld = PsiNew;
          Lower  = PsiOld;
          if (fabs(C2New) > Small)
            Y = Ro.Mag() + R.Mag() - 
                (VarA * (1.0 - PsiOld * C3New) / sqrt(C2New));
          else
            Y = Ro.Mag() + R.Mag();
          if (Show == 'Y')
            if (FileOut != NULL)
              fprintf(FileOut, "%3d %10.5f %10.5f %10.5f %7.3f %9.5f %9.5f\n",
                      Loops, PsiOld, Y, XOld, dtNew, VarA, Upper, Lower);

          YNegKtr++;
          if ((Y >= 0.0) || (YNegKtr >= 10))
            break;
        }
      }

      if (YNegKtr < 10)
      {
        if (fabs(C2New) > Small)
          XOld = sqrt(Y / C2New);
        else
          XOld = 0.0;
        XOldCubed = XOld * XOld * XOld;
        dtNew     = XOldCubed * C3New + VarA * sqrt(Y);

        /* ----  Readjust upper and lower bounds ---- */
        if (dtNew < DtTU)
          Lower = PsiOld;
        if (dtNew > DtTU)
          Upper = PsiOld;
        PsiNew = (Upper + Lower) * 0.5;

        if (Show == 'Y')
          if (FileOut != NULL)
            fprintf(FileOut, "%3d %10.5f %10.5f %10.5f %7.3f %9.5f %9.5f\n",
                    Loops, PsiOld, Y, XOld, dtNew, VarA, Upper, Lower);

        /* -------------- Find c2 and c3 functions ---------- */
        FindC2C3(PsiNew, C2New, C3New);
        PsiOld = PsiNew;
        Loops++;

        /* ---- Make sure the first guess isn't too close --- */
        if ((fabs(dtNew - DtTU) < Small) && (Loops == 1))
          dtNew = DtTU - 1.0;
      }
      
      if ((fabs(dtNew - DtTU) < Small) || (Loops > NumIter) || (YNegKtr > 10))
        break;
    }

    if ((Loops >= NumIter) || (YNegKtr >= 10))
    {
      strcpy(Error, "GNotConv");
      if (YNegKtr >= 10)
        strcpy(Error, "Y negative");
    }
    else
    {
      /* ---- Use F and G series to find Velocity Vectors ----- */
      F    = 1.0 - Y / Ro.Mag();
      GDot = 1.0 - Y / R.Mag();
      G    = 1.0 / (VarA * sqrt(Y)); // 1 over G
      FDot = sqrt(Y) * (-R.Mag() - Ro.Mag() + Y) / (R.Mag() * Ro.Mag() * VarA);
      for (UINT i = 0; i <= 3; i++)
      {
        Vo.Set((R.Get(i) - F * Ro.Get(i)) * G, i);
        V.Set((GDot * R.Get(i) - Ro.Get(i)) * G, i);
      }
    }
  }
  else
    strcpy(Error, "impos 180ø");

/*
    ---- For Fig 6-14 dev with testgau.pas ---- 
    IF Error = 'ok' THEN Write( FileOut,PsiNew:14:7,DtTU*13.44685:14:7 )
     ELSE Write( FileOut,' 9999.0 ':14,DtTU*13.44685:14:7 );
*/
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE TARGET
|
|  This PROCEDURE accomplishes the targeting problem using KEPLER/PKEPLER and
|    Lambert.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInt        - Initial Position vector of Int ER
|    VInt        - Initial Velocity vector of Int ER/TU
|    RTgt        - Initial Position vector of Tgt ER
|    VTgt        - Initial Velocity vector of Tgt ER/TU
|    dm          - Direction of Motion for Gauss  'L','S'
|    Kind        - Type of propagator             'K','P'
|    DtTU        - Time of flight to the int      TU
|
|  Outputs       :
|    V1t         - Initial Transfer Velocity vec  ER/TU
|    V2t         - Final Transfer Velocity vec    ER/TU
|    DV1         - Initial Change Velocity vec    ER/TU
|    DV2         - Final Change Velocity vec      ER/TU
|    Error       - Error flag from Gauss          'ok', ...
|
|  Locals        :
|    TransNormal - CROSS product of trans orbit   ER
|    IntNormal   - CROSS product of int orbit     ER
|    R1Tgt       - Position vector after Dt, Tgt  ER
|    V1Tgt       - Velocity vector after Dt, Tgt  ER/TU
|    RIRT        - RInt[4] * R1Tgt[4]
|    CosDeltaNu  - Cosine of DeltaNu              rad
|    SinDeltaNu  - Sine of DeltaNu                rad
|    DeltaNu     - ANGLE between position vectors rad
|    i           - Index
|
|  Coupling      :
|    KEPLER      - Find R and V at future time
|    LAMBERTUNIV - Find velocity vectors at each END of transfer
|    LNCOM2      - Linear combination of two vectors and constants
|
|  References    :
|    Vallado       2001, 468-474, Alg 58
|
 -----------------------------------------------------------------------------*/
void Target
    (
      Vector RInt, Vector VInt, Vector RTgt, Vector VTgt, 
      char Dm, char Kind, Real DtTU,
      Vector& V1t, Vector& V2t, Vector& DV1, Vector& DV2, char* Error
    )
{
  Vector IntNormal(3), TransNormal(3), R1Tgt(3), V1Tgt(3);
  Real   Temp, RIRT, CosDeltaNu, SinDeltaNu, DeltaNu;

  /* ----------- Propagate TARGET forward in time ----------------- */
  switch (Kind)
  {
    case 'K':
    case 'k':
      Kepler(RTgt, VTgt, DtTU,  R1Tgt, V1Tgt, Error);
      break;
    case 'P':
    case 'p':
//      PKepler(RTgt, VTgt, DtTU,  R1Tgt, V1Tgt, Error);
      break;
    default:
      Kepler(RTgt, VTgt, DtTU,  R1Tgt, V1Tgt, Error);
      break;
  }

  /* ----------- Calculate transfer orbit between r's ------------- */
  if (strcmp(Error, "ok") == 0)
  {
    LambertUniv(RInt, R1Tgt, Dm, 'N', DtTU,  V1t, V2t, Error);

    if (strcmp(Error, "ok") == 0)
    {
      DV1 = V1t - VInt;
      DV2 = V1Tgt - V2t;
    }
    else
    {
      V1t.Clear();
      V2t.Clear();
      DV1.Clear();
      DV2.Clear();
    }
  }
}

/* Utility functions for LambertBattin, etc */
static Real k(Real v)
{
  Real d[21] = 
       {
            1.0 /     3.0,     4.0 /    27.0,
            8.0 /    27.0,     2.0 /     9.0,
           22.0 /    81.0,   208.0 /   891.0,
          340.0 /  1287.0,   418.0 /  1755.0,
          598.0 /  2295.0,   700.0 /  2907.0,
          928.0 /  3591.0,  1054.0 /  4347.0,
         1330.0 /  5175.0,  1480.0 /  6075.0,
         1804.0 /  7047.0,  1978.0 /  8091.0,
         2350.0 /  9207.0,  2548.0 / 10395.0,
         2968.0 / 11655.0,  3190.0 / 12987.0,
         3658.0 / 14391.0
       };
  Real del, delold, term, termold, temp, sum1;
  SINT i;

  /* ---- Process Forwards ---- */
  sum1    = d[0];
  delold  = 1.0;
  termold = d[0];
  i = 1;
  while ((i <= 20) && (fabs(termold) > 0.000001))
  {
    del  = 1.0 / ( 1.0 - d[i] * v * delold);
    term = termold * (del - 1.0);
    sum1 = sum1 + term;
    i++;
    delold  = del;
    termold = term;
  }
  return sum1;
}

static Real See(Real v)
{
  Real c[21] =
       {
           0.2,
           9.0 /   35.0,   16.0 /   63.0,
          25.0 /   99.0,   36.0 /  143.0,
          49.0 /  195.0,   64.0 /  255.0,
          81.0 /  323.0,  100.0 /  399.0,
         121.0 /  483.0,  144.0 /  575.0,
         169.0 /  675.0,  196.0 /  783.0,
         225.0 /  899.0,  256.0 / 1023.0,
         289.0 / 1155.0,  324.0 / 1295.0,
         361.0 / 1443.0,  400.0 / 1599.0,
         441.0 / 1763.0,  484.0 / 1935.0
       };
  Real term, termold, del, delold, sum1, temp, eta, SQRTOpv;
  SINT i;

  SQRTOpv = sqrt(1.0 + v);
  eta     = v / Power(1.0 + SQRTOpv, 2);

  /* ---- Process Forwards ---- */
  delold  = 1.0;
  termold = c[0];  // * eta
  sum1    = termold;
  i = 1;
  while ((i <= 20) && (fabs(termold) > 0.000001))
  {
    del  = 1.0 / (1.0 + c[i] * eta * delold);
    term = termold * (del - 1.0);
    sum1 = sum1 + term;
    i++;
    delold  = del;
    termold = term;
  }

  return ((1.0 / (8.0 * (1.0 + SQRTOpv))) * (3.0 + sum1 / (1.0 + eta * sum1)));
}
