/*-----------------------------------------------------------------

                               TESTASTF.PAS

|   This file contains sample subroutines to test each of the functions
|   in the astrodynamic libraries.
|
|                                                                            
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

        -----------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>

#include "ast2body.h"
#include "astiod.h"
#include "astmath.h"
#include "asttime.h"
#include "astutil.h"
#include "astreduc.h"
#include "constants.h"

FILE *InFile, *OutFile;
char *IFileName = "testastf.dat";
char *OFileName = "testastf.out";
char *Blank4 = "    ";
char  Title[131];
char  junque[10];
bool  keepitup = false;

/* ------------------------------------------------------------------------
|
|                       Individual PROCEDURES
|
|  Forward references for test subroutines.
|
|
  ------------------------------------------------------------------------*/

void TestADDVEC(void);
void TestADD3VEC(void);
void TestANGLE(void);
void TestARCCOSH(void);
void TestARCSINH(void);
void TestARCTANH(void);
void TestBINOMIAL(void);
void TestCROSS(void);
void TestCUBIC(void);
void TestDETERMINANT(void);
void TestDOT(void);
void TestFACTOR(void);
void TestFACTORIAL(void);
void TestFILEEXPPRINTMAT(void);
void TestFILEPRINTMAT(void);
void TestGETPART(void);
void TestGETPARTL(void);
void TestGetPartR(void);
void TestLNCOM1(void);
void TestLNCOM2(void);
void TestLNCOM3(void);
void TestMAG(void);
void TestMAKEMAT(void);
void TestMATADD(void);
void TestMATINVERSE(void);
void TestMATMULT(void);
void TestMATSCALE(void);
void TestMATSUB(void);
void TestMATTRANS(void);
void TestMAX(void);
void TestMIN(void);
void TestNORM(void);
void TestPLANE(void);
void TestPOLYFIT(void);
void TestPRINTMAT(void);
void TestQUADRATIC(void);
void TestQUARTIC(void);
void TestROT1(void);
void TestROT2(void);
void TestROT3(void);
void TestSINH(void);

/* -------------------- Time s ----------------------------------- */
void TestDAYS2MDHMS(void);
void TestDAYLIGHTST(void);
void TestDAYOFWEEK(void);
void TestDMS_RAD(void);
void TestFINDDAYS(void);
void TestGETINTDAY(void);
void TestGETINTMON(void);
void TestHMS_RAD(void);
void TestHMS_SEC(void);
void TestHMS_UT(void);
void TestINITTIME(void);
void TestINVJULIANDAY(void);
void TestJULIANDAY(void);
void TestJULIANDAYALL(void);
void TestLSTIME(void);
void TestMOONRISESET(void);
void TestSUNRISESET(void);
void TestUpCaseSt(void);

/* --------------------  2 body s -------------------------------- */
void TestCHECKHITEARTH(void);
void TestELORB(void);
void TestFINDC2C3(void);
void TestFINDTOF(void);
void TestGEOCGEOD(void);
void TestIJKTOLATLONA(void);
void TestKEPLER(void);
void TestLIGHT(void);
void TestMOON(void);
void TestNEWTONE(void);
void TestNEWTONM(void);
void TestNEWTONNU(void);
void TestPATH(void);
void TestRANDV(void);
void TestRNGAZ(void);
void TestSATFOV(void);
void TestSIGHT(void);
void TestSUN(void);

/* --------------------  IOD s ----------------------------------- */
void TestANGLESGAUSS(void);
void TestANGLESLAPLACE(void);
void TestGIBBS(void);
void TestHERRGIBBS(void);
void TestLAMBERTBATTIN(void);
void TestLAMBERTUNIV(void);
void TestRADEC_AZEL(void);
void TestRADEC_ELATLON(void);
void TestRV_ELATLON(void);
void TestRV_RADEC(void);
void TestRV_RAZEL(void);
void TestRV_TRADEC(void);
void TestRVSEZ_RAZEL(void);
void TestSITE(void);
void TestTARGET(void);

/* --------------------  REDUC's ---------------------------------- */
void TestConvTime(void);
void TestFK4(void);
void TestInitNutation(void);
void TestNutation(void);
void TestPolarM(void);
void TestPrecession(void);
void TestSidereal(void);
void TestTrueMean(void);

/* -------------------- MANV -------------------------------------- */
void TestHohmann(void);
void TestBiElliptic(void);
void TestOneTangent(void);
void TestIOnlyChg(void);
void TestNodeOnlyChg(void);
void TestIandNodeChg(void);
void TestMinCombinedPlaneChg(void);
void TestCombinedPlaneChg(void);
void TestRendezvous(void);
void TestNonCoplanarRendz(void);
void TestHillsR(void);
void TestHillsV(void);
void TestIJK_RSW(void);
void TestCow2Hill(void);
void TestPKEPLER(void);
void TestJ2DragPert(void);
void TestPredict(void);
void TestTestInitGravityField(void);
void TestLegPoly(void);
void TestDeriv(void);
void TestPertaccel(void);
void TestPderiv(void);
void TestRK4(void);
void TestRKF45(void);
void TestCowell(void);
void TestAtmos(void);
void TestNonlin(void);
void TestFindAtwaAtwb(void);
void TestLeastSquares(void);
void TestSequential(void);


/* ------------------------------------------------------------------------
|
|                       Individual PROCEDURES
|
|  These subroutines implement each function or procedure using a datafile
|    with the number of each routine to identify which subroutine to use.
|
|
  ------------------------------------------------------------------------*/
void TestADDVEC(void)
{
  float v1, v2, v3;
  Vector V1(3), V2(3), V3(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);

  V3 = V1 + V2;

  fprintf(OutFile, "  Results: %f %f %f\n", V3.Get(1), V3.Get(2), V3.Get(3));
  fgets(Title, 120, InFile);
}

void TestADD3VEC(void)
{
  float  v1, v2, v3;
  Vector V1(3), V2(3), V3(3), V4(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V3.Set((Real)v1, 1);
  V3.Set((Real)v2, 2);
  V3.Set((Real)v3, 3);

  V4 = V1 + V2 + V3;

  fprintf(OutFile, "  Results: %f %f %f\n", V4.Get(1), V4.Get(2), V4.Get(3));
  fgets(Title, 120, InFile);
}

void TestANGLE(void)
{
  float  v1, v2, v3;
  Vector V1(3), V2(3), V3(3), V4(3);
  Real   Theta;

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);

  Theta = V1.Angle(V2);

  fprintf(OutFile, "  Results: %f\n", Theta);
  fgets(Title, 120, InFile);
}

void TestARCCOSH(void)
{
  Real  AC;
  float XVal;

  fscanf(InFile, "%f", &XVal);

  AC = Arccosh(XVal);

  fprintf(OutFile, "  Results: %f\n", AC);
  fgets(Title, 120, InFile);
}

void TestARCSINH(void)
{
  float XVal;
  Real  AS;

  fscanf(InFile, "%f", &XVal);

  AS = Arctanh(XVal);

  fprintf(OutFile, "  Results: %f\n", AS);
  fgets(Title, 120, InFile);
}

void TestARCTANH(void)
{
  float XVal;
  Real  AT;

  fscanf(InFile, "%f", &XVal);

  AT = Arctanh(XVal);

  fprintf(OutFile, "  Results: %f\n", AT);
  fgets(Title, 120, InFile);
}

void TestBINOMIAL(void)
{
  UINT i, j;
  Real Bio;

  fscanf(InFile, "%d %d", &i, &j);

  Bio = Binomial(i, j);

  fprintf(OutFile, "  Results: %f\n", Bio);
  fgets(Title, 120, InFile);
}

void TestCROSS(void)
{
  float v1, v2, v3;
  Vector V1(3), V2(3), V3(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);

  V3 = V1.Cross(V2);

  fprintf(OutFile, "  Results: %f %f %f\n", V3.Get(1), V3.Get(2), V3.Get(3));
  fgets(Title, 120, InFile);
}

void TestCUBIC(void)
{
  float a, b, c, d;
  Real  R1r, R1i, R2r, R2i, R3r, R3i;

  fscanf(InFile, "%f %f %f %f", &a, &b, &c, &d);

  Cubic(a, b, c, d, R1r, R1i, R2r, R2i, R3r, R3i);

  fprintf(OutFile, "  Results: %f %f %f %f %f %f\n", 
                    R1r, R1i, R2r, R2i, R3r, R3i);
  fgets(Title, 120, InFile);
}

void TestDETERMINANT(void)
{
  SINT  Order;

  fscanf(InFile, "%d", &Order);
  fprintf(OutFile, "%d\n", Order);

  {
    Matrix Mat1(Order, Order), Mat2(Order, Order);
    float v;
    Real  D;

    for (UINT i = 1; i <= Order; i++)
      for (UINT j = 1; j <= Order; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }

    D = Mat1.Determinant();

    fprintf(OutFile, "  Results: %f\n", D);
  }
  fgets(Title, 120, InFile);
}

void TestDOT(void)
{
  Real Dt;
  float v1, v2, v3;
  Vector V1(3), V2(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);

  Dt = V1.Dot(V2);

  fprintf(OutFile, "  Results: %f\n", Dt);
  fgets(Title, 120, InFile);
}

void TestFACTOR(void)
{
  SINT NRootS;
  float v;

  fscanf(InFile, "%d", &NRootS);

  {
    Real Poly[NRootS], Roots[NRootS][NRootS];

    for (UINT i = 1; i <= NRootS; i++)
    {
      fscanf(InFile, "%f", &v);
      Poly[i-1] = (Real)v;
      fprintf(OutFile, "%f ", v);
    }
    fprintf(OutFile, "\n");

//    Factor(Poly, NRootS, (Real**)Roots);

    fprintf(OutFile, "  Results: \n");
    for (UINT i = 1; i <= NRootS; i++)
      fprintf(OutFile, "%f %f\n", Roots[i-1][0], Roots[i-1][1]);
  }
  fgets(Title, 120, InFile);
}

void TestFACTORIAL(void)
{
  UINT x;
  Real Fact;

  fscanf(InFile, "%d", &x);

  Fact = Factorial(x);

  fprintf(OutFile, "  Results: %f\n", Fact);
  fgets(Title, 120, InFile);
}

void TestFILEEXPPRINTMAT(void)
{
  SINT r, c;
  char title[15];

  fscanf(InFile, "%d %d", &r, &c);
  fprintf(OutFile, "%d %d\n", &r, &c);
  fgets(Title, 120, InFile);

  {
    Matrix m(r, c);
    float  v;

    for (UINT i = 1; i <= r; i++)
      for (UINT j = 1; j <= c; j++)
      {
        fscanf(InFile, "%f", &v);
        m.Set((Real)v, i, j);
      }

    FilePrint(m, "sample sale", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestFILEPRINTMAT(void)
{
  SINT r, c;
  char title[15];

  fscanf(InFile, "%d %d", &r, &c);
  fprintf(OutFile, "%d %d\n", &r, &c);
  fgets(Title, 120, InFile);

  {
    Matrix m(r, c);
    float  v;

    for (UINT i = 1; i <= r; i++)
      for (UINT j = 1; j <= c; j++)
      {
        fscanf(InFile, "%f", &v);
        m.Set((Real)v, i, j);
      }

    FilePrint(m, "sample sale", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestGETPART(void)
{
  SINT i, LocStart, Leng;
  char InStr[250];
  char *loc;

  fgets(InStr, 120, InFile);
  InStr[strlen(InStr)-2] = '\0';
  loc = RmSpcs(InStr);
  fprintf(OutFile, "%s\n", InStr);
  fscanf(InFile, "%d %d", &LocStart, &Leng);

  i = GetPart(loc, LocStart,Leng);

  fprintf(OutFile, "  Results: %d\n", i);
  fgets(Title, 120, InFile);
}

void TestGETPARTL(void)
{
  SINT LocStart, Leng;
  LINT i;
  char InStr[250];
  char *loc;

  fgets(InStr, 120, InFile);
  InStr[strlen(InStr)-2] = '\0';
  loc = RmSpcs(InStr);
  fprintf(OutFile, "%s\n", InStr);
  fscanf(InFile, "%d %d", &LocStart, &Leng);

  i = GetPartL(loc, LocStart,Leng);

  fprintf(OutFile, "  Results: %d\n", i);
  fgets(Title, 120, InFile);
}

void TestGetPartR(void)
{
  SINT LocStart, Leng;
  Real i;
  char InStr[250];
  char *loc;

  fgets(InStr, 120, InFile);
  InStr[strlen(InStr)-2] = '\0';
  loc = RmSpcs(InStr);
  fprintf(OutFile, "%s\n", InStr);
  fscanf(InFile, "%d %d", &LocStart, &Leng);

  i = GetPartR(loc, LocStart,Leng);

  fprintf(OutFile, "  Results: %f\n", i);
  fgets(Title, 120, InFile);
}

void TestLNCOM1(void)
{
  float m1, v1, v2, v3;
  Vector V1(3), V2(3);

  fscanf(InFile, "%f %f %f %f", &m1, &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);

  V2 = m1 * V1;

  fprintf(OutFile, "  Results: %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3));
  fgets(Title, 120, InFile);
}

void TestLNCOM2(void)
{
  float m1, m2, v1, v2, v3;
  Vector V1(3), V2(3), V3(3);

  fscanf(InFile, "%f %f", &m1, &m2);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);

  V3 = m1 * V1 + m2 * V2;

  fprintf(OutFile, "  Results: %f %f %f\n", V3.Get(1), V3.Get(2), V3.Get(3));
  fgets(Title, 120, InFile);
}

void TestLNCOM3(void)
{
  float m1, m2, m3, v1, v2, v3;
  Vector V1(3), V2(3), V3(3), V4(3);

  fscanf(InFile, "%f %f %f", &m1, &m2, &m3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V2.Set((Real)v1, 1);
  V2.Set((Real)v2, 2);
  V2.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V3.Set((Real)v1, 1);
  V3.Set((Real)v2, 2);
  V3.Set((Real)v3, 3);

  V4 = m1 * V1 + m2 * V2 + m3 * V3;

  fprintf(OutFile, "  Results: %f %f %f\n", V4.Get(1), V4.Get(2), V4.Get(3));
  fgets(Title, 120, InFile);
}

void TestMAG(void)
{
  float v1, v2, v3;
  Vector V(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V.Set((Real)v1, 1);
  V.Set((Real)v2, 2);
  V.Set((Real)v3, 3);

  fprintf(OutFile, "  Results: %f\n", V.Mag());
  fgets(Title, 120, InFile);
}

void TestMAKEMAT(void)
{
  float Angl;
  SINT  Numbr;

  fscanf(InFile, "%f %d", &Angl, &Numbr);

  fprintf(OutFile, "  Results: %f %d\n", Angl, Numbr);
  fgets(Title, 120, InFile);
}

void TestMATADD(void)
{
  SINT Mat1r, Mat1c;

  fscanf(InFile, "%d %d", &Mat1r, &Mat1c);
  fprintf(OutFile, "%d %d\n", Mat1r, Mat1c);

  {
    Matrix Mat1(Mat1r, Mat1c), Mat2(Mat1r, Mat1c), Mat3(Mat1r, Mat1c);
    float  v;

    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }
    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat2.Set((Real)v, i, j);
      }

    Mat3 = Mat1 + Mat2;

    FilePrint(Mat3, "  Results:", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestMATINVERSE(void)
{
  SINT  Order;
  float v;

  fscanf(InFile, "%d", &Order);
  fprintf(OutFile, "%d\n", Order);

  {
    Matrix Mat1(Order, Order), Mat2(Order, Order);

    for (UINT i = 1; i <= Order; i++)
      for (UINT j = 1; j <= Order; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }

    Mat2 = Mat1.Inverse();

    FilePrint(Mat2, "  Results:", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestMATMULT(void)
{
  SINT Mat1r, Mat1c, Mat2c;

  fscanf(InFile, "%d %d %d", &Mat1r, &Mat1c, &Mat2c);
  fprintf(OutFile, "%d %d %d\n", Mat1r, Mat1c, &Mat2c);

  {
    Matrix Mat1(Mat1r, Mat1c), Mat2(Mat1c, Mat2c), Mat3(Mat1r, Mat2c);
    float  v;

    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }
    for (UINT i = 1; i <= Mat1c; i++)
      for (UINT j = 1; j <= Mat2c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat2.Set((Real)v, i, j);
      }

    Mat3 = Mat1 * Mat2;

    FilePrint(Mat3, "  Results:", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestMATSCALE(void)
{
  SINT  Mat1r, Mat1c;
  float Scale;

  fscanf(InFile, "%d %d %f", &Mat1r, &Mat1c, &Scale);
  fprintf(OutFile, "%d %d %f\n", Mat1r, Mat1c, Scale);

  {
    Matrix Mat1(Mat1r, Mat1c), Mat2(Mat1r, Mat1c);
    float v;

    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }
    Mat2 = Mat1.Scale((Real)Scale);

    FilePrint(Mat2, "  Results:", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestMATSUB(void)
{
  SINT Mat1r, Mat1c;

  fscanf(InFile, "%d %d", &Mat1r, &Mat1c);
  fprintf(OutFile, "%d %d\n", Mat1r, Mat1c);

  {
    Matrix Mat1(Mat1r, Mat1c), Mat2(Mat1r, Mat1c), Mat3(Mat1r, Mat1c);
    float  v;

    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }
    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat2.Set((Real)v, i, j);
      }

    Mat3 = Mat1 - Mat2;

    FilePrint(Mat3, "  Results:", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestMATTRANS(void)
{
  SINT Mat1r, Mat1c;

  fscanf(InFile, "%d %d", &Mat1r, &Mat1c);
  fprintf(OutFile, "%d %d\n", Mat1r, Mat1c);

  {
    Matrix Mat1(Mat1r, Mat1c), Mat2(Mat1c, Mat1r);
    float  v;

    for (UINT i = 1; i <= Mat1r; i++)
      for (UINT j = 1; j <= Mat1c; j++)
      {
        fscanf(InFile, "%f", &v);
        Mat1.Set((Real)v, i, j);
      }

    Mat2 = Mat1.Transpose();

    FilePrint(Mat2, "  Results:", 4, OutFile);
  }
  fgets(Title, 120, InFile);
}

void TestMAX(void)
{
  Real Mx;
  SINT m, n;

  fscanf(InFile, "%d %d", &m, &n);

  Mx = Max(m, n);

  fprintf(OutFile, "  Results: %f\n", Mx);
  fgets(Title, 120, InFile);
}

void TestMIN(void)
{
  Real Mn;
  SINT m, n;

  fscanf(InFile, "%d %d", &m, &n);

  Mn = Min(m, n);

  fprintf(OutFile, "  Results: %f\n", Mn);
  fgets(Title, 120, InFile);
}

void TestNORM(void)
{
  float v1, v2, v3;
  Vector V1(3), V2(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);

  V2 = V1.Norm();

  fprintf(OutFile, "  Results: %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3));
  fgets(Title, 120, InFile);
}

void TestPLANE(void)
{
  float x1, y1, z1, x2, y2, z2, x3, y3, z3;
  Real  a, b, c, d;

  fscanf(InFile, "%f %f %f", &x1, &y1, &z1);
  fprintf(OutFile, "%f %f %f\n", x1, y1, z1);
  fscanf(InFile, "%f %f %f", &x2, &y2, &z2);
  fprintf(OutFile, "%f %f %f\n", x2, y2, z2);
  fscanf(InFile, "%f %f %f", &x3, &y3, &z3);
  fprintf(OutFile, "%f %f %f\n", x3, y3, z3);

  Plane((Real)x1, (Real)y1, (Real)z1, (Real)x2, (Real)y2, (Real)z2, (Real)x3,
        (Real)y3, (Real)z3, a,b,c,d);

  fprintf(OutFile, "  Results: %f %f %f %f\n", a, b, c, d);
  fgets(Title, 120, InFile);
}

void TestPOLYFIT(void)
{
  SINT   Degree, NumPts;
  Matrix DataPoints(10,2) ,Coeff(3,3);
  float  a, b;
  Real   MinX, MinY;

  fscanf(InFile, "%d %d", &Degree, &NumPts);
  fprintf(OutFile, "%d %d\n", Degree, NumPts);
  for (UINT i = 1; i <= NumPts; i++)
  {
    fscanf(InFile, "%f %f", &a, &b);
    DataPoints.Set(a, i, 1);
    DataPoints.Set(b, i, 2);
    fprintf(OutFile, "%f %f\n", DataPoints.Get(i, 1), DataPoints.Get(i, 2));
  }

  Polyfit(Degree, NumPts, DataPoints, Coeff, MinX, MinY);

  fprintf(OutFile, "  Results:\n");
    for (UINT i = 1; i <= NumPts; i++)
    fprintf(OutFile, "%f ", Coeff.Get(i, 1));
  fprintf(OutFile, "\n");
  fgets(Title, 120, InFile);
}

void TestPRINTMAT(void)
{
  SINT r, c;
  char title[15];

  fscanf(InFile, "%d %d", &r, &c);
  fprintf(OutFile, "%d %d\n", &r, &c);
  fgets(Title, 120, InFile);

  {
    Matrix m(r, c);
    float  v;

    for (UINT i = 1; i <= r; i++)
      for (UINT j = 1; j <= c; j++)
      {
        fscanf(InFile, "%f", &v);
        m.Set((Real)v, i, j);
      }

    Print(m, "sample", 4);
  }
  fgets(Title, 120, InFile);
}

void TestQUADRATIC(void)
{
  float a, b, c;
  Real  R1r, R1i, R2r, R2i;

  fscanf(InFile, "%f %f %f", &a, &b, &c);

  Quadratic(a, b, c, R1r, R1i, R2r, R2i);

  fprintf(OutFile, "  Results: %f %f %f %f\n", R1r, R1i, R2r, R2i);
  fgets(Title, 120, InFile);
}

void TestQUARTIC(void)
{
  float a, b, c, d, e;
  Real  R1r, R1i, R2r, R2i, R3r, R3i, R4r, R4i;

  fscanf(InFile, "%f %f %f %f %f", &a, &b, &c, &d, &e);

  Quartic(a, b, c, d, e, R1r, R1i, R2r, R2i, R3r, R3i, R4r, R4i);

  fprintf(OutFile, "  Results: %f %f %f %f\n", R1r, R1i, R2r, R2i);
  fgets(Title, 120, InFile);
}

void TestROT1(void)
{
  float r, v1, v2, v3;
  Vector V1(3), V2(3);

  fscanf(InFile, "%f %f %f %f", &r, &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);

  V2 = V1.Rot1(r);

  fprintf(OutFile, "  Results: %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3));
  fgets(Title, 120, InFile);
}

void TestROT2(void)
{
  float r, v1, v2, v3;
  Vector V1(3), V2(3);

  fscanf(InFile, "%f %f %f %f", &r, &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);

  V2 = V1.Rot2(r);

  fprintf(OutFile, "  Results: %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3));
  fgets(Title, 120, InFile);
}

void TestROT3(void)
{
  float r, v1, v2, v3;
  Vector V1(3), V2(3);

  fscanf(InFile, "%f %f %f %f", &r, &v1, &v2, &v3);
  V1.Set((Real)v1, 1);
  V1.Set((Real)v2, 2);
  V1.Set((Real)v3, 3);

  V2 = V1.Rot3(r);

  fprintf(OutFile, "  Results: %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3));
  fgets(Title, 120, InFile);
}

void TestSINH(void)
{
  float XVal;
  Real  S;

  fscanf(InFile, "%f", &XVal);

  S = sinh(XVal);

  fprintf(OutFile, "  Results: %f\n", S);
  fgets(Title, 120, InFile);

}


/* -------------------- Time s ----------------------------------- */
void TestDAYS2MDHMS(void)
{
  float Days;
  SINT  Year, Mon, Day, Hr, Min;
  Real Sec;

  fscanf(InFile, "%d %f", &Year, &Days);
  fprintf(OutFile, "%d %f\n", Year, Days);

  Days2MDHMS(Year, Days,  Mon, Day, Hr, Min, Sec);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%d %d %d %d %f\n", Mon, Day, Hr, Min, Sec);
  fgets(Title, 120, InFile);
}

void TestDAYLIGHTST(void)
{
}

void TestDAYOFWEEK(void)
{
}

void TestDMS_RAD(void)
{
  SINT  Deg, Min;
  Real  DMS;
  float Sec;
  char dir[5];

  fscanf(InFile, "%d %d %f %s", &Deg, &Min, &Sec, dir);
  fprintf(OutFile, "%d %d %f :%s:\n", Deg, Min, Sec, dir);

  if (strcmp(dir, "TOO") == 0)
    DMS_Rad(Deg,Min, (Real)Sec, TOO, DMS);
  else
    DMS_Rad(Deg,Min, (Real)Sec, FROM, DMS);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f\n", DMS);
  fgets(Title, 120, InFile);
}

void TestFINDDAYS(void)
{
  SINT  Year, Mon, Day, Hr, Min;
  Real  Days;
  float Sec;

  fscanf(InFile, "%d %d %d %d %d %f", &Year, &Mon, &Day, &Hr, &Min, &Sec);
  fprintf(OutFile, "%d %d %d %d %d %f\n", Year, Mon, Day, Hr, Min, Sec);

  FindDays(Year, Mon, Day, Hr, Min, Sec, Days);

  fprintf(OutFile, "  Results: %f\n", Days);
  fgets(Title, 120, InFile);
}

void TestGETINTDAY(void)
{
}

void TestGETINTMON(void)
{
}

void TestHMS_RAD(void)
{
  SINT  Hr, Min;
  Real  HMS;
  float Sec;
  char dir[5];

  fscanf(InFile, "%d %d %f %s", &Hr, &Min, &Sec, dir);
  fprintf(OutFile, "%d %d %f :%s:\n", Hr, Min, Sec, dir);

  if (strcmp(dir, "TOO") == 0)
    HMS_Rad(Hr,Min, (Real)Sec, TOO, HMS);
  else
    HMS_Rad(Hr,Min, (Real)Sec, FROM, HMS);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f\n", HMS);
  fgets(Title, 120, InFile);
}

void TestHMS_SEC(void)
{
  SINT  Hr, Min;
  Real  UTSec;
  float Sec;
  char dir[5];

  fscanf(InFile, "%d %d %f %s", &Hr, &Min, &Sec, dir);
  fprintf(OutFile, "%d %d %f :%s:\n", Hr, Min, Sec, dir);

  if (strcmp(dir, "TOO") == 0)
    HMS_Sec(Hr,Min, (Real)Sec, TOO, UTSec);
  else
    HMS_Sec(Hr,Min, (Real)Sec, FROM, UTSec);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f\n", UTSec);
  fgets(Title, 120, InFile);
}

void TestHMS_UT(void)
{
  SINT  Hr, Min;
  Real  UT;
  float Sec;
  char dir[5];

  fscanf(InFile, "%d %d %f %s", &Hr, &Min, &Sec, dir);
  fprintf(OutFile, "%d %d %f :%s:\n", Hr, Min, Sec, dir);

  if (strcmp(dir, "TOO") == 0)
    HMS_UT(Hr,Min, (Real)Sec, TOO, UT);
  else
    HMS_UT(Hr,Min, (Real)Sec, FROM, UT);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f\n", UT);
  fgets(Title, 120, InFile);
}

void TestINITTIME(void)
{
}

void TestINVJULIANDAY(void)
{
  SINT  Year, Mon, Day, Hr, Min;
  float JD;
  Real  Sec;

  fscanf(InFile, "%f", &JD);
  fprintf(OutFile, "%f\n", JD);

  InvJulianDay(JD, Year, Mon, Day, Hr, Min, Sec);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%d %d %d %d %d %f\n", Year, Mon, Day, Hr, Min, Sec);
  fgets(Title, 120, InFile);
}

void TestJULIANDAY(void)
{
  SINT  Year, Mon, Day, Hr, Min;
  float Sec;
  Real  JD;

  fscanf(InFile, "%d %d %d %d %d %f", &Year, &Mon, &Day, &Hr, &Min, &Sec);
  fprintf(OutFile, "%d %d %d %d %d %f", Year, Mon, Day, Hr, Min, Sec);

  JulianDay(Year,Mon,Day,Hr,Min, (Real)Sec, JD);

  fprintf(OutFile, "  Results: %f\n", JD);
  fgets(Title, 120, InFile);
}

void TestJULIANDAYALL(void)
{
  SINT  Year, Mon, Day, Hr, Min;
  float Sec;
  Real  JD;
  char  WhichType;

  fscanf(InFile, "%d %d %d %d %d %f %c",
                 &Year, &Mon, &Day, &Hr, &Min, &Sec, &WhichType);
  fprintf(OutFile, "%d %d %d %d %d %f", Year, Mon, Day, Hr, Min, Sec);

  JulianDayAll(Year,Mon,Day,Hr,Min, (Real)Sec, WhichType, JD);

  fprintf(OutFile, "  Results: %f\n", JD);
  fgets(Title, 120, InFile);
}

void TestLSTIME(void)
{
  float Lon, JD;
  Real  LST, GST;
  const Real Rad = 180.0 / PI;

  fscanf(InFile, "%f %f", &Lon, &JD);
  fprintf(OutFile, "%f %f\n", Lon, JD);

  LSTime(Lon, JD, LST, GST);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", LST * Rad, GST * Rad);
  fgets(Title, 120, InFile);
}

void TestMOONRISESET(void)
{
  float JD, Latgd, Lon;
  Real  UTMoonRise, UTMoonSet, MoonPhaseAng;
  const Real Rad = PI / 180.0;
  char  Error[12];

  fscanf(InFile, "%f %f %f", &Lon, &Latgd, &JD);
  fprintf(OutFile, "%f %f %f\n", Lon, Latgd, JD);
  Latgd = Latgd * Rad;
  Lon   = Lon * Rad;

  MoonRiseSet(JD, Latgd,Lon, UTMoonRise, UTMoonSet, MoonPhaseAng, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %s\n", UTMoonRise, UTMoonSet, MoonPhaseAng, Error);
  fgets(Title, 120, InFile);
}

void TestSUNRISESET(void)
{
  float JD, Latgd, Lon;
  Real  UTSunRise, UTSunSet;
  const Real Rad = PI / 180.0;
  char  Whichkind;
  char  Error[12];

  fscanf(InFile, "%f %f %c", &Lon, &Latgd, &JD, &Whichkind);
  fprintf(OutFile, "%f %f\n", Lon, Latgd, JD);
  Latgd = Latgd * Rad;
  Lon   = Lon * Rad;

  SunRiseSet(JD, Latgd,Lon, Whichkind, UTSunRise, UTSunSet, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %s\n", UTSunRise, UTSunSet, Error);
  fgets(Title, 120, InFile);
}

void TestUpCaseSt(void)
{
  char S[250], UpS[250];

  fscanf(InFile, "%s", S);
  fprintf(OutFile, "%s\n", S);

  UpCaseStr(S);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%s\n", S);
  fgets(Title, 120, InFile);
}


/* --------------------  2 body s -------------------------------- */
void TestCHECKHITEARTH(void)
{
  Vector RInt(3), V1t(3), RTgt(3), V2t(3);
  float  v1, v2, v3;
  char HitEarth;

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  RInt.Clear();
  RInt.Set(v1, 1);
  RInt.Set(v2, 2);
  RInt.Set(v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  V1t.Clear();
  V1t.Set(v1, 1);
  V1t.Set(v2, 2);
  V1t.Set(v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  RTgt.Clear();
  RTgt.Set(v1, 1);
  RTgt.Set(v2, 2);
  RTgt.Set(v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  V2t.Clear();
  V2t.Set(v1, 1);
  V2t.Set(v2, 2);
  V2t.Set(v3, 3);


  CheckHitEarth(RInt, V1t, RTgt, V2t, HitEarth);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%c\n", HitEarth);
  fgets(Title, 120, InFile);
}

void TestELORB(void)
{
  const Real Rad = 180.0 / PI;
  float v1, v2, v3;
  Vector R(3), V(3);
  Real P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat, TrueLon, LonPer;

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R.Set(v1, 1);
  R.Set(v2, 2);
  R.Set(v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  V.Set(v1, 1);
  V.Set(v2, 2);
  V.Set(v3, 3);

  ElOrb(R, V, P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat, TrueLon, LonPer);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%14.8f %14.8f %14.8f\n", P, A, Ecc);
  fprintf(OutFile, "%14.8f %14.8f %14.8f %14.8f %14.8f\n",
                    Incl * Rad, Omega * Rad, Argp * Rad, Nu * Rad, M * Rad);
  fprintf(OutFile, "%14.8f %14.8f %14.8f\n",
                    ArgLat * Rad, TrueLon * Rad, LonPer * Rad);
  fgets(Title, 120, InFile);
}

void TestFINDC2C3(void)
{
  float ZNew;
  Real  C2New, C3New;

  fscanf(InFile, "%f", &ZNew);
  fprintf(OutFile, "%f\n", ZNew);

  FindC2C3(ZNew, C2New, C3New);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", C2New, C3New);
  fgets(Title, 120, InFile);
}

void TestFINDTOF(void)
{
  Vector Ro(3), R(3);
  float  p, v1, v2, v3;
  Real  Tof;
  
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Ro.Set(v1, 1);
  Ro.Set(v2, 2);
  Ro.Set(v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3, &p);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3, p);
  R.Set(v1, 1);
  R.Set(v2, 2);
  R.Set(v3, 3);

  FindTOF(Ro, R, p, Tof);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f\n", Tof);
  fgets(Title, 120, InFile);
}

void TestGEOCGEOD(void)
{
  const Real Rad = 180.0 / PI;
  char dir[6];
  float Latgc;
  Real  Latgd;

  fscanf(InFile, "%f %s", &Latgc, dir);
  fprintf(OutFile, "%f %s\n", Latgc, dir);
  Latgc = Latgc / Rad;

  if (strcmp(dir, "TOO") == 0)
    GeocGeod((Real)Latgc, TOO, Latgd);
  else
    GeocGeod((Real)Latgc, FROM, Latgd);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f\n", Latgd * Rad);
  fgets(Title, 120, InFile);
}

void TestIJKTOLATLONA(void)
{
  const  Real Rad = 180.0 / PI;
  Vector R(3);
  float  JD, v1, v2, v3;
  Real   LatGc, LatGd, Lon, Hellp;

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R.Set(v1, 1);
  R.Set(v2, 2);
  R.Set(v3, 3);
  fscanf(InFile, "%f", &JD);
  fprintf(OutFile, "%f\n", JD);

  IJKtoLatLonA( R, JD, LatGc, LatGd, Lon, Hellp );

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f\n", LatGc * Rad, LatGd * Rad, Lon * Rad, Hellp);
  fgets(Title, 120, InFile);
}

void TestKEPLER(void)
{
  Vector Ro(3), Vo(3), R(3), V(3);
  float  DtTU, v1, v2, v3;
  char Error[12];

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Ro.Set(v1, 1);
  Ro.Set(v2, 2);
  Ro.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Vo.Set(v1, 1);
  Vo.Set(v2, 2);
  Vo.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f", &DtTU);
  fprintf(OutFile, "%f\n", DtTU);

  Kepler(Ro, Vo, DtTU, R, V, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f\n", R.Get(1), R.Get(2), R.Get(3));
  fprintf(OutFile, "%f %f %f\n", V.Get(1), V.Get(2), V.Get(3));
  fgets(Title, 120, InFile);
}

void TestLIGHT(void)
{
  float JD, v1, v2, v3;
  Vector R;
  char   WhichKind;
  char   Lit[4];

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  fscanf(InFile, "%f %c", &JD, &WhichKind);
  fprintf(OutFile, "%f %c", JD, WhichKind);

  Light(R, JD, WhichKind, Lit);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%s\n", Lit);
  fgets(Title, 120, InFile);
}

void TestMOON(void)
{
  float  JD;
  Real   RtAsc, Decl;
  Vector RMoon(3);

  fscanf(InFile, "%f", &JD);
  fprintf(OutFile, "%f\n", JD);

  Sun(JD, RMoon, RtAsc, Decl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n",
                    RMoon.Get(1),RMoon.Get(1),RMoon.Get(1), RMoon.Mag());
  fgets(Title, 120, InFile);
}

void TestNEWTONE(void)
{
  const Real Rad = 180.0 / PI;
  float Ecc, E0;
  Real  M, Nu;

  fscanf(InFile, "%f %f", &Ecc, &E0);
  fprintf(OutFile, "%f %f\n", Ecc, E0);
  E0 = E0 / Rad;

  NewtonE(Ecc, E0, M, Nu);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", M * Rad, Nu * Rad);
  fgets(Title, 120, InFile);
}

void TestNEWTONM(void)
{
  const Real Rad = 180.0 / PI;
  float Ecc, M;
  Real  E0, Nu;

  fscanf(InFile, "%f %f", &Ecc, &M);
  fprintf(OutFile, "%f %f\n", Ecc, M);
  M = M / Rad;

  NewtonE(Ecc, M, E0, Nu);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", E0 * Rad, Nu * Rad);
  fgets(Title, 120, InFile);
}

void TestNEWTONNU(void)
{
  const Real Rad = 180.0 / PI;
  float Ecc, Nu;
  Real  M, E0;

  fscanf(InFile, "%f %f", &Ecc, &Nu);
  fprintf(OutFile, "%f %f\n", Ecc, Nu);
  Nu = Nu / Rad;

  NewtonE(Ecc, Nu, E0, M);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", E0 * Rad, M * Rad);
  fgets(Title, 120, InFile);
}

void TestPATH(void)
{
  const Real Rad = 180.0 / PI;
  float LLat, LLon, Range, Az;
  Real  TLat, TLon;

  fscanf(InFile, "%f %f %f %f", &LLat, &LLon, &Range, &Az);
  fprintf(OutFile, "%f %f %f %f\n", LLat, LLon, Range, Az);
  LLat = LLat / Rad;
  LLon = LLon / Rad;
  Az   = Az / Rad;

  Path(LLat, LLon, Range, Az, TLat, TLon);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", TLat * Rad, TLon * Rad);
  fgets(Title, 120, InFile);
}

void TestRANDV(void)
{
  const Real Rad = PI / 180.0;
  float P, a, Ecc, Incl, Omega, Argp, Nu, ArgLat, TrueLon, LonPer;
  Vector R(3), V(3);

  fscanf(InFile, "%f %f %f %f %f %f %f %f %f %f",
         &P, &a, &Ecc, &Incl, &Omega, &Argp, &Nu, &ArgLat, &TrueLon, &LonPer);
  fprintf(OutFile, "%f %f %f %f %f %f %f %f %f %f\n",
          P, a, Ecc, Incl, &Omega, Argp, Nu, ArgLat, TrueLon, LonPer);
  Incl    = Incl * Rad;
  Omega   = Omega * Rad;
  Argp    = Argp * Rad;
  Nu      = Nu * Rad;
  ArgLat  = ArgLat * Rad;
  TrueLon = TrueLon * Rad;
  LonPer  = LonPer * Rad;

  RandV(P, Ecc, Incl, Omega, Argp, Nu, ArgLat, TrueLon, LonPer, R, V);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%14.8f %14.8f %14.8f %14.8f\n",
                    R.Get(1), R.Get(2), R.Get(3), R.Mag());
  fprintf(OutFile, "%14.8f %14.8f %14.8f %14.8f\n",
                    V.Get(1), V.Get(2), V.Get(3), V.Mag());
  fgets(Title, 120, InFile);
}

void TestRNGAZ(void)
{
  const Real Rad = 180.0 / PI;
  float LLat, LLon, TLat, TLon;
  Real  ToF, Range, Az;

  fscanf(InFile, "%f %f %f %f", &LLat, &LLon, &TLat, &TLon);
  fprintf(OutFile, "%f %f %f %f\n", LLat, LLon, TLat, TLon);
  LLat = LLat / Rad;
  LLon = LLon / Rad;
  TLat = TLat / Rad;
  TLon = TLon / Rad;

  RngAz(LLat, LLon, TLat, TLon, ToF, Range, Az);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f\n", ToF, Range, Az * Rad);
  fgets(Title, 120, InFile);
}

void TestSATFOV(void)
{
  const Real Rad = 180.0 / PI;
  float Incl, Az, SLatgd, SLon, SAlt, tFOV, EtaCtr, FOVMax;
  Real  TotalRng, RhoMax, RhoMin,TgtLat, TgtLon;

  fscanf(InFile, "%f %f %f %f %f %f %f %f",
                 &Incl, &Az, &SLatgd, &SLon, &SAlt, &tFOV, &EtaCtr, &FOVMax);
  fprintf(OutFile, "%f %f %f %f %f %f %f %f\n",
                    Incl, Az, SLatgd, SLon, SAlt, tFOV, EtaCtr, FOVMax);
  Incl   = Incl / Rad;
  Az     = Az / Rad;
  SLatgd = SLatgd / Rad;
  SLon   = SLon / Rad;

  SatFOV((Real)Incl, (Real)Az, (Real)SLatgd, (Real)SLon, (Real)SAlt,
         (Real)tFOV, (Real)EtaCtr, (Real)FOVMax, TotalRng,
         RhoMax, RhoMin, TgtLat, TgtLon);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f %f\n",
                    TotalRng, RhoMax, RhoMin, TgtLat, TgtLon);
  fgets(Title, 120, InFile);
}

void TestSIGHT(void)
{
  Vector R1(3), R2(3);
  float  v1, v2, v3;
  char   WhichKind;
  char   LOS[5];

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R1.Set(v1, 1);
  R1.Set(v2, 2);
  R1.Set(v3, 3);
  fscanf(InFile, "%f %f %f %c", &v1, &v2, &v3, &WhichKind);
  fprintf(OutFile, "%f %f %f %c\n", v1, v2, v3, WhichKind);
  R2.Set(v1, 1);
  R2.Set(v2, 2);
  R2.Set(v3, 3);

  Sight(R1,R2, WhichKind, LOS);
 
  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%s\n", LOS);
  fgets(Title, 120, InFile);
}

void TestSUN(void)
{
  float  JD;
  Real   RtAsc, Decl;
  Vector RSun(3);

  fscanf(InFile, "%f", &JD);
  fprintf(OutFile, "%f\n", JD);

  Sun(JD, RSun, RtAsc, Decl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", 
                    RSun.Get(1),RSun.Get(1),RSun.Get(1), RSun.Mag());
  fgets(Title, 120, InFile);
}


/* --------------------  IOD s ----------------------------------- */
void TestANGLESGAUSS(void)
{
  float Delta1, Delta2, Delta3, Alpha1, Alpha2, Alpha3, JD1, JD2, JD3;
  Vector RS1(3), RS2(3), RS3(3), r2(3), v2(3);

  fscanf(InFile, "%f %f %f %f %f %f %f %f %f",
       &Delta1, &Delta2, &Delta3, &Alpha1, &Alpha2, &Alpha3, &JD1, &JD2, &JD3);
  fprintf(OutFile, "%f %f %f %f %f %f %f %f %f\n",
          Delta1, Delta2, Delta3, Alpha1, Alpha2, Alpha3, JD1, JD2, JD3);

  AnglesGauss(Delta1, Delta2, Delta3, Alpha1, Alpha2, Alpha3, JD1, JD2, JD3,
                RS1, RS2, RS3, r2, v2);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", r2.Get(1), r2.Get(2), r2.Get(3), r2.Mag());
  fprintf(OutFile, "%f %f %f %f\n", v2.Get(1), v2.Get(2), v2.Get(3), v2.Mag());
  fgets(Title, 120, InFile);
}

void TestANGLESLAPLACE(void)
{
  float Delta1, Delta2, Delta3, Alpha1, Alpha2, Alpha3, JD1, JD2, JD3;
  Vector RS1(3), RS2(3), RS3(3), r2(3), v2(3);

  fscanf(InFile, "%f %f %f %f %f %f %f %f %f",
       &Delta1, &Delta2, &Delta3, &Alpha1, &Alpha2, &Alpha3, &JD1, &JD2, &JD3);
  fprintf(OutFile, "%f %f %f %f %f %f %f %f %f\n",
          Delta1, Delta2, Delta3, Alpha1, Alpha2, Alpha3, JD1, JD2, JD3);

  AnglesLaplace(Delta1, Delta2, Delta3, Alpha1, Alpha2, Alpha3, JD1, JD2, JD3,
                RS1, RS2, RS3, r2, v2);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", r2.Get(1), r2.Get(2), r2.Get(3), r2.Mag());
  fprintf(OutFile, "%f %f %f %f\n", v2.Get(1), v2.Get(2), v2.Get(3), v2.Mag());
  fgets(Title, 120, InFile);
}

void TestGIBBS(void)
{
  const Real Rad =  180.0 / PI;
  float v1, v2, v3;
  Real  Theta, Theta1, Copa;
  char  Error[15];
  Vector R1(3), R2(3), R3(3), V2(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R1.Set(v1, 1);
  R1.Set(v2, 2);
  R1.Set(v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R2.Set(v1, 1);
  R2.Set(v2, 2);
  R2.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R3.Set(v1, 1);
  R3.Set(v2, 2);
  R3.Set(v3, 3);

  Gibbs(R1, R2, R3, V2, Theta, Theta1, Copa, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3), V2.Mag());
  fprintf(OutFile, "angle %f %f %f %s\n",
                    Theta * Rad, Theta1 * Rad, Copa * Rad, Error);
  fgets(Title, 120, InFile);
}

void TestHERRGIBBS(void)
{
  const Real Rad =  180.0 / PI;
  float t1, t2, t3, v1, v2, v3;
  Real JD1, JD2, JD3, Theta, Theta1, Copa;
  char  Error[15];
  Vector R1(3), R2(3), R3(3), V2(3);

  fscanf(InFile, "%f %f %f %f", &v1, &v2, &v3, &t1);
  JD1 = 2448608.0L + t1 / 86400.0L;
  fprintf(OutFile, "%f %f %f %f\n", v1, v2, v3, JD1);
  R1.Set(v1, 1);
  R1.Set(v2, 2);
  R1.Set(v3, 3);
  fscanf(InFile, "%f %f %f %f", &v1, &v2, &v3, &t2);
  JD1 = 2448608.0L + t2 / 86400.0L;
  fprintf(OutFile, "%f %f %f %f\n", v1, v2, v3, JD2);
  R2.Set(v1, 1);
  R2.Set(v2, 2);
  R2.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %f %f %f", &v1, &v2, &v3, &t3);
  JD1 = 2448608.0L + t3 / 86400.0L;
  fprintf(OutFile, "%f %f %f %f\n", v1, v2, v3, JD3);
  R3.Set(v1, 1);
  R3.Set(v2, 2);
  R3.Set(v3, 3);

  HerrGibbs(R1,R2,R3, JD1, JD2, JD3, V2, Theta, Theta1, Copa, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", V2.Get(1), V2.Get(2), V2.Get(3), V2.Mag());
  fgets(Title, 120, InFile);
}

void TestLAMBERTBATTIN(void)
{
  char Dm, OverRev;
  char Error[15];
  float DtTU, v1, v2, v3;
  Vector Ro(3), R(3), Vo(3), V(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Ro.Set(v1, 1);
  Ro.Set(v2, 2);
  Ro.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R.Set(v1, 1);
  R.Set(v2, 2);
  R.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %c %c", &DtTU, &Dm, &OverRev);
  fprintf(OutFile, "%c %c %f\n", Dm, OverRev, DtTU);

  LambertBattin(Ro, R, Dm, OverRev, DtTU, Vo,V, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", Vo.Get(1), Vo.Get(2), Vo.Get(3), Vo.Mag());
  fprintf(OutFile, "%f %f %f %f\n", V.Get(1), V.Get(2), V.Get(3), V.Mag());
  fgets(Title, 120, InFile);
}

void TestLAMBERTUNIV(void)
{
  char Dm, OverRev;
  char Error[15];
  float DtTU, v1, v2, v3;
  Vector Ro(3), R(3), Vo(3), V(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Ro.Set(v1, 1);
  Ro.Set(v2, 2);
  Ro.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  R.Set(v1, 1);
  R.Set(v2, 2);
  R.Set(v3, 3);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%f %c %c", &DtTU, &Dm, &OverRev);
  printf("     %c %c %12.8f\n", Dm, OverRev, DtTU);
  fprintf(OutFile, "%c:%c %12.8f\n", Dm, OverRev, DtTU);

  LambertUniv(Ro, R, Dm, OverRev, DtTU, Vo,V, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", Vo.Get(1), Vo.Get(2), Vo.Get(3), Vo.Mag());
  fprintf(OutFile, "%f %f %f %f\n", V.Get(1), V.Get(2), V.Get(3), V.Mag());
  fgets(Title, 120, InFile);
}

void TestRADEC_AZEL(void)
{
  const Real Rad =  180.0 / PI;
  char  dir[5];
  float RtAsc, Decl, LST, LatGc;
  Real  Az, El;

  fscanf(InFile, "%f %f %f %f %s", &RtAsc, &Decl, &LST, &LatGc, dir);
  fprintf(OutFile, "%f %f %f %f %s\n", RtAsc, Decl, LST, LatGc, dir);

  if (strcmp(dir, "TOO") == 0)
    RaDec_AzEl((Real)RtAsc, (Real)Decl, (Real)LST, (Real)LatGc, TOO, Az, El);
  else
    RaDec_AzEl((Real)RtAsc, (Real)Decl, (Real)LST, (Real)LatGc, FROM, Az, El);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", Az * Rad, El * Rad);
  fgets(Title, 120, InFile);
}

void TestRADEC_ELATLON(void)
{
  const Real Rad =  180.0 / PI;
  char  dir[5];
  float  RtAsc, Decl;
  Real EclLat, EclLon;

  fscanf(InFile, "%f %f %s", &RtAsc, &Decl, dir);
  fprintf(OutFile, "%f %f %s\n", RtAsc, Decl, dir);
  RtAsc = RtAsc / Rad;
  Decl  = Decl  / Rad;

  if (strcmp(dir, "TOO") == 0)
    RaDec_ELatLon((Real)RtAsc, (Real)Decl, TOO, EclLat, EclLon);
  else
    RaDec_ELatLon((Real)RtAsc, (Real)Decl, FROM, EclLat, EclLon);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", EclLat * Rad, EclLon * Rad);

  RaDec_ELatLon((Real)RtAsc, (Real)Decl, FROM, EclLat, EclLon);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f\n", RtAsc * Rad, Decl * Rad);
  fgets(Title, 120, InFile);
}

void TestRV_ELATLON(void)
{
  const Real Rad           =  180.0 / PI;
  const Real RadiusEarthKm = 6378.1363;
  const Real VKmPerSec     =    7.905366149846;
  float v1, v2, v3;
  char  dir[5];
  Real  rr, EclLat, EclLon, DRr, DEclLat, DEclLon;
  Vector Rijk(3), Vijk(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Rijk.Set((Real)v1, 1);
  Rijk.Set((Real)v2, 2);
  Rijk.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f %s", &v1, &v2, &v3, dir);
  fprintf(OutFile, "%f %f %f %s\n", v1, v2, v3, dir);
  Vijk.Set((Real)v1, 1);
  Vijk.Set((Real)v2, 2);
  Vijk.Set((Real)v3, 3);

  Rijk = Rijk / RadiusEarthKm;
  Vijk = Vijk / VKmPerSec;

  if (strcmp(dir, "TOO") == 0)
    RV_ELatLon(Rijk, Vijk, TOO, rr, EclLat, EclLon, DRr, DEclLat, DEclLon);
  else
    RV_ELatLon(Rijk, Vijk, FROM, rr, EclLat, EclLon, DRr, DEclLat, DEclLon);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f %f %f\n",
                    rr, EclLat, EclLon, DRr, DEclLat, DEclLon);

  RV_ELatLon(Rijk, Vijk, FROM, rr, EclLat, EclLon, DRr, DEclLat, DEclLon);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n",
                    Rijk.Get(1), Rijk.Get(2), Rijk.Get(3), Rijk.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                    Vijk.Get(1), Vijk.Get(2), Vijk.Get(3), Vijk.Mag());
  fgets(Title, 120, InFile);
}

void TestRV_RADEC(void)
{
  const Real Rad           =  180.0 / PI;
  const Real RadiusEarthKm = 6378.1363;
  const Real VKmPerSec     =    7.905366149846;
  float v1, v2, v3;
  char  dir[5];
  Real  rr, RtAsc, Decl, DRr, DRtAsc, DDecl;
  Vector Rijk(3), Vijk(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  Rijk.Set((Real)v1, 1);
  Rijk.Set((Real)v2, 2);
  Rijk.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f %s", &v1, &v2, &v3, dir);
  fprintf(OutFile, "%f %f %f %s\n", v1, v2, v3, dir);
  Vijk.Set((Real)v1, 1);
  Vijk.Set((Real)v2, 2);
  Vijk.Set((Real)v3, 3);

  Rijk = Rijk / RadiusEarthKm;
  Vijk = Vijk / VKmPerSec;

  if (strcmp(dir, "TOO") == 0)
    RV_RaDec(Rijk, Vijk, TOO, rr, RtAsc, Decl, DRr, DRtAsc, DDecl);
  else
    RV_RaDec(Rijk, Vijk, FROM, rr, RtAsc, Decl, DRr, DRtAsc, DDecl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f %f %f\n",
                    rr, RtAsc * Rad, Decl * Rad, DRr, DRtAsc, DDecl);

  RV_RaDec(Rijk, Vijk, FROM, rr, RtAsc, Decl, DRr, DRtAsc, DDecl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n",
                    Rijk.Get(1), Rijk.Get(2), Rijk.Get(3), Rijk.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                    Vijk.Get(1), Vijk.Get(2), Vijk.Get(3), Vijk.Mag());
  fgets(Title, 120, InFile);
}

void TestRV_RAZEL(void)
{
  const Real Rad = 180.0 / PI;
  float v1, v2, v3, v4, v5;
  char  dir[5];
  Real  Lat, LST, Rho, Az, El, DRho, DAz, DEl; 
  Vector Rijk(3), Vijk(3), RS(3);

  fscanf(InFile, "%f %f %f %f", &v1, &v2, &v3, &v4);
  Rijk.Set((Real)v1, 1);
  Rijk.Set((Real)v2, 2);
  Rijk.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f %f", &v1, &v2, &v3, &v4);
  Vijk.Set((Real)v1, 1);
  Vijk.Set((Real)v2, 2);
  Vijk.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f %f %f %s", &v1, &v2, &v3, &v4, &v5, dir);
  fprintf(OutFile, "%f %f %f %s\n", v1, v2, v3, dir);
  RS.Set((Real)v1, 1);
  RS.Set((Real)v2, 2);
  RS.Set((Real)v3, 3);
  Lat = Lat / Rad;
  LST = LST / Rad;

    if (strcmp(dir, "TOO") == 0)
    RV_RAzEl(Rijk, Vijk, RS, Lat,LST, TOO, Rho, Az, El, DRho, DAz, DEl);
  else
    RV_RAzEl(Rijk, Vijk, RS, Lat,LST, FROM, Rho, Az, El, DRho, DAz, DEl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f %f %f\n", Rho, Az, El, DRho, DAz, DEl);

  RV_RAzEl(Rijk, Vijk, RS, Lat,LST, FROM, Rho, Az, El, DRho, DAz, DEl);

  fprintf(OutFile, "  Results:\n");
    fprintf(OutFile, "%f %f %f %f\n",
                    Rijk.Get(1), Rijk.Get(2), Rijk.Get(3), Rijk.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                    Vijk.Get(1), Vijk.Get(2), Vijk.Get(3), Vijk.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                    RS.Get(1), RS.Get(2), RS.Get(3), RS.Mag());
  fgets(Title, 120, InFile);
}

void TestRV_TRADEC(void)
{
  const Real Rad           = 180.0 / PI;
  const Real RadiusEarthKm = 6378.1363;
  const Real VKmPerSec     =    7.905366149846;
  float v1, v2, v3;
  char  dir[5];
  Real  Rho, TRtAsc, TDecl, DRho, DTRtAsc, DTDecl;
  Vector Rijk(3), Vijk(3), RS(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  Rijk.Set((Real)v1, 1);
  Rijk.Set((Real)v2, 2);
  Rijk.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  Vijk.Set((Real)v1, 1);
  Vijk.Set((Real)v2, 2);
  Vijk.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f %s", &v1, &v2, &v3, dir);
  fprintf(OutFile, "%f %f %f %s\n", v1, v2, v3, dir);
  RS.Set((Real)v1, 1);
  RS.Set((Real)v2, 2);
  RS.Set((Real)v3, 3);

  Rijk = Rijk / RadiusEarthKm;
  RS   = RS   / RadiusEarthKm;
  Vijk = Vijk / VKmPerSec;

  if (strcmp(dir, "TOO") == 0)
    RV_TRaDec(Rijk, Vijk, RS, TOO, Rho, TRtAsc, TDecl, DRho, DTRtAsc, DTDecl);
  else
    RV_TRaDec(Rijk, Vijk, RS, FROM, Rho, TRtAsc, TDecl, DRho, DTRtAsc, DTDecl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f %f %f\n",
                   Rho, TRtAsc * Rad, TDecl * Rad, DRho, DTRtAsc, DTDecl);

  RV_TRaDec(Rijk, Vijk, RS, FROM, Rho, TRtAsc, TDecl, DRho, DTRtAsc, DTDecl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n",
                    Rijk.Get(1), Rijk.Get(2), Rijk.Get(3), Rijk.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                    Vijk.Get(1), Vijk.Get(2), Vijk.Get(3), Vijk.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                    RS.Get(1), RS.Get(2), RS.Get(3), RS.Mag());
  fgets(Title, 120, InFile);
}

void TestRVSEZ_RAZEL(void)
{
  const Real Rad = 180.0 / PI;
  float v1, v2, v3;
  char  dir[5];
  Real  Rho, Az, El, DRho, DAz, DEl;
  Vector Rhosez(3), DRhosez(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  Rhosez.Set((Real)v1, 1);
  Rhosez.Set((Real)v2, 2);
  Rhosez.Set((Real)v3, 3);
  fscanf(InFile, "%f %f %f %s", &v1, &v2, &v3, dir);
  DRhosez.Set((Real)v1, 1);
  DRhosez.Set((Real)v2, 2);
  DRhosez.Set((Real)v3, 3);

  if (strcmp(dir, "TOO") == 0)
    RVSez_RAzEl(Rhosez, DRhosez, TOO, Rho, Az, El, DRho, DAz, DEl);
  else
    RVSez_RAzEl(Rhosez, DRhosez, FROM, Rho, Az, El, DRho, DAz, DEl);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f %f %f\n",
                   Rho, Az * Rad, El * Rad, DRho, DAz, DEl);

  RVSez_RAzEl(Rhosez, DRhosez, FROM, Rho, Az, El, DRho, DAz, DEl);

    fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n",
                 Rhosez.Get(1), Rhosez.Get(2), Rhosez.Get(3), Rhosez.Mag());
  fprintf(OutFile, "%f %f %f %f\n",
                 DRhosez.Get(1), DRhosez.Get(2), DRhosez.Get(3), DRhosez.Mag());
  fgets(Title, 120, InFile);
}

void TestSITE(void)
{
  const Real Rad = 180.0 / PI;
  const Real Reft = 20925644.0288713L;
  float Latgd, Alt, LST;
  Vector RS(3), VS(3);

  fscanf(InFile, "%f %f %f", &Latgd, &Alt, &LST);
  fprintf(OutFile, "%f %f %f\n", Latgd, Alt, LST);
  Latgd = Latgd / Rad;
  LST   = LST / Rad;
  Alt   = Alt / Reft;

  Site(Latgd, Alt, LST, RS, VS);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "%f %f %f %f\n", RS.Get(1), RS.Get(2), RS.Get(3), RS.Mag());
  fprintf(OutFile, "%f %f %f %f\n", VS.Get(1), VS.Get(2), VS.Get(3), VS.Mag());
  fgets(Title, 120, InFile);
}

void TestTARGET(void)
{
  float v1, v2, v3, v4, v5, v6;
  Real  DtTU;
  char Dm, Kind;
  char Error[15];
  Vector RInt(3), VInt(3), RTgt(3), VTgt(3), V1t(3), V2t(3), DV1(3), DV2(3);

  fscanf(InFile, "%f %f %f %f %f %f %c", &v1, &v2, &v3, &v4, &v5, &v6, &Kind);
  fprintf(OutFile, "%f %f %f %f %f %f %c\n", v1, v2, v3, v4, v5, v6, Kind);

//  Target(RInt, VInt, RTgt, VTgt, Dm, Kind, DtTU, V1t, V2t, DV1, DV2, Error);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "Not implemented\n");
  fgets(Title, 120, InFile);
}


/* --------------------  REDUC's ----------------------------------- */
void TestConvTime(void)
{
  SINT Year, Mon, Day, Hr, MIN, TimeZone;
  char TypeUTIn;
  char FileN1[64];
  char Error[15];
  Real DUT1, DAT, xp, yp, UT1, TUT1, JDUT1, UTC, TAI, TDT, TTDT,
       JDTDT, TDB, TTDB, JDTDB;
  float Sec;

  fscanf(InFile, "%s", FileN1);
  fprintf(OutFile, "%s\n", FileN1);
  fscanf(InFile, "%d %d %d %d %d %f %d %c",
                  &Year, &Mon, &Day, &Hr, &Min, &Sec, &TimeZone, &TypeUTIn);
  fprintf(OutFile, "%d %d %d %d %d %f %d %c",
                  Year,Mon,Day,Hr,Min,Sec,TimeZone,TypeUTIn);

//  ConvTime(FileN1, Year, Mon, Day, Hr, MIN, Sec, TimeZone, TypeUTIn,
//           DUT1, DAT, xp, yp, UT1, TUT1, JDUT1, UTC, TAI, TDT, TTDT,
//           JDTDT, TDB, TTDB, JDTDB, Error);

  fprintf(OutFile, "  Results:\n");
//  fprintf(OutFile, "%f %f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f\n",
//                    DUT1, DAT, xp, yp, UT1, TUT1, JDUT1, UTC, TAI, TDT,
//                    TTDT, JDTDT, TDB, TTDB, JDTDB);
  fprintf(OutFile, "Parameters differ between test and astreduc definition\n");
  fgets(Title, 120, InFile);
}

void TestFK4(void)
{
  char dir[5];
  float v1, v2, v3;
  Vector rJ2000(3), vJ2000(3), rFK4(3), vFK4(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "In r %f %f %f\n", v1, v2, v3);
  rJ2000.Set(v1, 1);
  rJ2000.Set(v2, 2);
  rJ2000.Set(v3, 3);
  fscanf(InFile, "%f %f %f%s", &v1, &v2, &v3,dir);
  fprintf(OutFile, "   v %f %f %f %s\n", v1, v2, v3, dir);
  vJ2000.Set(v1, 1);
  vJ2000.Set(v2, 2);
  vJ2000.Set(v3, 3);

  if (strcmp(dir, "TOO") == 0)
    FK4(rJ2000, vJ2000, TOO, rFK4, vFK4);
  else
    FK4(rJ2000, vJ2000, FROM, rFK4, vFK4);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "r %f %f %f\n", rFK4.Get(1), rFK4.Get(2), rFK4.Get(3));
  fprintf(OutFile, "v %f %f %f\n", vFK4.Get(1), vFK4.Get(2), vFK4.Get(3));
  fgets(Title, 120, InFile);
}

void TestInitNutation(void)
{
  char FileN1[65], FileN2[65];
  SINT IAr[5][106], iIAr[5][106], pIAr[10][112];
  Real RAr[4][106], iRAr[4][106], pRAr[4][112];

  fscanf(InFile, "%s", FileN1);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%s", FileN2);
  fprintf(OutFile, "%s %s\n", FileN1, FileN2);

//  InitNutation(FileN1, FileN2, IAr, RAr, iIAr, iRAr, pIAr, pRAr);

  fprintf(OutFile, "  Results:\n");
//  fprintf(OutFile, "%d %d\n", IAr[0][0], RAr[0][0] * 3600.0 / 0.0001);
//  fprintf(OutFile, "%d %d\n", IAr[1][0], RAr[1][0] * 3600.0 / 0.0001);
//  fprintf(OutFile, "%d %d\n", IAr[2][0], RAr[2][0] * 3600.0 / 0.0001);
  fgets(Title, 120, InFile);
}

void TestNutation(void)
{
  const Real Rad = 180.0 / PI;
  char  dir[5];
  Real  RAr[4][106], iRAr[4][4], pRAr[4][85];
  SINT  IAr[5][106], iIAr[5][4], pIAr[10][85];
  float DeltaPsi, TrueEps, TTDB, v1, v2, v3;
  char  FileN1[65], FileN2[65];
  Vector rTOD(3), vTOD(3), rMOD(3), vMOD(3);
  
  fscanf(InFile, "%s", FileN1);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%s", FileN2);
  fprintf(OutFile, "%s %s\n", FileN1, FileN2);
//  InitNutation(FileN1, FileN2, IAr, RAr, iIAr, iRAr, pIAr, pRAr);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  rMOD.Set(v1, 1);
  rMOD.Set(v2, 2);
  rMOD.Set(v3, 3);
  fscanf(InFile, "%f %f %f %f %f %f %s", &v1, &v2, &v3,
                  &DeltaPsi, &TrueEps, &TTDB, dir);
  fprintf(OutFile, "%f %f %f %s\n", v1, v2, v3, dir);
  vMOD.Set(v1, 1);
  vMOD.Set(v2, 2);
  vMOD.Set(v3, 3);


//  if (strcmp(dir, "TOO") == 0)
//    Nutation(rMOD, vMOD, TOO, rTOD, vTOD, DeltaPsi,TrueEps, TTDB, IAr,RAr,
//             iIAr,iRAr, pIAr,pRAr);
//  else
//    Nutation(rMOD, vMOD, FROM, rTOD, vTOD, DeltaPsi,TrueEps, TTDB, IAr,RAr,
//             iIAr,iRAr, pIAr,pRAr);

  fprintf(OutFile, "  Results:\n");
//  fprintf(OutFile, "r %f %f %f\n", rTOD.Get(1), rTOD.Get(2), rTOD.Get(3));
//  fprintf(OutFile, "v %f %f %f\n", vTOD.Get(1), vTOD.Get(2), vTOD.Get(3));
  fprintf(OutFile, "Parameters to InitNutation and Nutation wrong\n");
  fgets(Title, 120, InFile);
}

void TestPolarM(void)
{
  const Real Rad = 180.0 / PI;
  char dir[5];
  float xp, yp, v1, v2, v3;
  Vector rECEFwpm(3), vECEFwpm(3), rECEF(3), vECEF(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "In r %f %f %f\n", v1, v2, v3);
  rECEF.Set(v1, 1);
  rECEF.Set(v2, 2);
  rECEF.Set(v3, 3);
  fscanf(InFile, "%f %f %f %f %f %s", &v1, &v2, &v3, &xp, &yp, dir);
  fprintf(OutFile, "   v %f %f %f %s\n", v1, v2, v3, dir);
  vECEF.Set(v1, 1);
  vECEF.Set(v2, 2);
  vECEF.Set(v3, 3);
  xp = xp / (3600.0 * Rad);
  yp = yp / (3600.0 * Rad);

  if (strcmp(dir, "TOO") == 0)
    PolarM(rECEF, vECEF, TOO, rECEFwpm, vECEFwpm, xp, yp);
  else
    PolarM(rECEF, vECEF, FROM, rECEFwpm, vECEFwpm, xp, yp);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "r %f %f %f\n", 
                    rECEFwpm.Get(1), rECEFwpm.Get(2), rECEFwpm.Get(3));
  fprintf(OutFile, "v %f %f %f\n", 
                    vECEFwpm.Get(1), vECEFwpm.Get(2), vECEFwpm.Get(3));
  fgets(Title, 120, InFile);
}

void TestPrecession(void)
{
  char dir[5];
  float  TTDB, v1, v2, v3;
  Vector rJ2000(3), vJ2000(3), rMOD(3), vMOD(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "%f %f %f\n", v1, v2, v3);
  rJ2000.Set(v1, 1);
  rJ2000.Set(v2, 2);
  rJ2000.Set(v3, 3);
  fscanf(InFile, "%f %f %f %f %s", &v1, &v2, &v3, &TTDB, dir);
  fprintf(OutFile, "%f %f %f :%s: %f\n", v1, v2, v3, dir, TTDB);
  vJ2000.Set(v1, 1);
  vJ2000.Set(v2, 2);
  vJ2000.Set(v3, 3);

  rJ2000 = rJ2000 / 6378.1363;
  vJ2000 = vJ2000 / 7.905366149846;

  if (strcmp(dir, "TOO") == 0)
    Precession(rJ2000, vJ2000, TOO, rMOD, vMOD, TTDB);
  else
    Precession(rJ2000, vJ2000, FROM, rMOD, vMOD, TTDB);

  rMOD = rMOD * 6378.1363;
  vMOD = vMOD * 7.905366149846;

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "r %f %f %f\n", rMOD.Get(1), rMOD.Get(2), rMOD.Get(3));
  fprintf(OutFile, "v %f %f %f\n", vMOD.Get(1), vMOD.Get(2), vMOD.Get(3));
  fgets(Title, 120, InFile);
}

void TestSidereal(void)
{
  const Real Rad = 180.0 / PI;
  char dir[5];
  Real Omega, XLOD;
  float JDUT1, DeltaPsi, TrueEps, v1, v2, v3;
  Vector rTOD(3), vTOD(3), rECEF(3), vECEF(3);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "IN r %f %f %f\n", v1, v2, v3);
  rTOD.Set(v1, 1);
  rTOD.Set(v2, 2);
  rTOD.Set(v3, 3);
  fscanf(InFile, "%f %f %f %f %f %f %s",
                  &v1, &v2, &v3, &JDUT1, &DeltaPsi, &TrueEps, dir);
  fprintf(OutFile, "%f %f %f %f %f %f %s\n",
                    v1, v2, v3, JDUT1, DeltaPsi, TrueEps, dir);
  vTOD.Set(v1, 1);
  vTOD.Set(v2, 2);
  vTOD.Set(v3, 3);
  DeltaPsi = DeltaPsi / Rad;
  TrueEps  = TrueEps / Rad;

  rTOD = rTOD / 6378.1363;
  vTOD = vTOD / 7.905366149846;

  XLOD  = 0.0;
  Omega = 0.5;
  if (strcmp(dir, "TOO") == 0)
    Sidereal(rTOD, vTOD, TOO, rECEF, vECEF, JDUT1, Omega, 
             DeltaPsi, TrueEps, XLOD);
  else
    Sidereal(rTOD, vTOD, FROM, rECEF, vECEF, JDUT1, Omega, 
             DeltaPsi, TrueEps, XLOD);

  rECEF = rECEF * 6378.1363;
  vECEF = vECEF * 7.905366149846;

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "r %f %f %f\n", rECEF.Get(1), rECEF.Get(2), rECEF.Get(3));
  fprintf(OutFile, "v %f %f %f\n", vECEF.Get(1), vECEF.Get(2), vECEF.Get(3));
  fgets(Title, 120, InFile);
}

void TestTrueMean(void)
{
  const Real Rad = 180.0 / PI;
  char dir[5];
  char FileN1[65], FileN2[65];
  SINT IAr[5][106], iIAr[5][4], pIAr[10][85];
  float DeltaPsi, TTdb, TrueEps, v1, v2, v3;
  Real RAr[4][106], iRAr[4][4], pRAr[4][85];
  Vector rMOD(3), vMOD(3), rTM(3), vTM(3);
  
  fscanf(InFile, "%s", FileN1);
  fgets(Title, 120, InFile);
  fscanf(InFile, "%s", FileN2);
  fprintf(OutFile, "%s %s\n", FileN1, FileN2);

//  InitNutation(FileN1, FileN2, IAr, RAr, iIAr, iRAr, pIAr, pRAr);

  fscanf(InFile, "%f %f %f", &v1, &v2, &v3);
  fprintf(OutFile, "IN r %f %f %f\n", v1, v2, v3);
  rMOD.Set(v1, 1);
  rMOD.Set(v2, 2);
  rMOD.Set(v3, 3);
  fscanf(InFile, "%f %f %f %f %f %f %s",
                  &v1, &v2, &v3, &DeltaPsi, &TrueEps, &TTdb, dir);
  fprintf(OutFile, "%f %f %f %f %f %f %s\n",
                    v1, v2, v3, DeltaPsi, TrueEps, TTdb, dir);
  vMOD.Set(v1, 1);
  vMOD.Set(v2, 2);
  vMOD.Set(v3, 3);
  DeltaPsi = DeltaPsi / Rad;
  TrueEps  = TrueEps / Rad;

  if (strcmp(dir, "TOO") == 0)
    TrueMean(rMOD, vMOD, TOO, rTM, vTM, (Real)DeltaPsi, (Real)TrueEps,
             TTdb, IAr, RAr);
  else
    TrueMean(rMOD, vMOD, FROM, rTM, vTM, (Real)DeltaPsi, (Real)TrueEps,
             TTdb, IAr, RAr);

  fprintf(OutFile, "  Results:\n");
  fprintf(OutFile, "r %f %f %f\n", rTM.Get(1), rTM.Get(2), rTM.Get(3));
  fprintf(OutFile, "v %f %f %f\n", vTM.Get(1), vTM.Get(2), vTM.Get(3));
  fgets(Title, 120, InFile);
}


/* -------------------- MANV -------------------------------------- */
void TestHohmann(void)
{
}

void TestBiElliptic(void)
{
}

void TestOneTangent(void)
{
}

void TestIOnlyChg(void)
{
}

void TestNodeOnlyChg(void)
{
}

void TestIandNodeChg(void)
{
}

void TestMinCombinedPlaneChg(void)
{
}

void TestCombinedPlaneChg(void)
{
}

void TestRendezvous(void)
{
}

void TestNonCoplanarRendz(void)
{
}

void TestHillsR(void)
{
}

void TestHillsV(void)
{
}

void TestIJK_RSW(void)
{
}

void TestCow2Hill(void)
{
}

void TestPKEPLER(void)
{
}

void TestJ2DragPert(void)
{
}

void TestPredict(void)
{
}

void TestInitGravityField(void)
{
}

void TestLegPoly(void)
{
}

void TestDeriv(void)
{
}

void TestPertaccel(void)
{
}

void TestPderiv(void)
{
}

void TestRK4(void)
{
}

void TestRKF45(void)
{
}

void TestCowell(void)
{
}

void TestAtmos(void)
{
}

void TestNonlin(void)
{
}

void TestFindAtwaAtwb(void)
{
}

void TestLeastSquares(void)
{
}

void TestSequential(void)
{
}



/* --------------------------- M A I N ----------------------------- */
int main(void)
{
  SINT  i, incnt, NumChose;

  if ((InFile = fopen(IFileName, "r")) == NULL)
  {
    printf("Unable to open input test file %s\n", IFileName);
    exit (0);
  }

  if ((OutFile = fopen(OFileName, "w")) == NULL)
  {
    printf("Unable to open output test file %s\n", OFileName);
    exit (0);
  }

  InitTime();

  keepitup = true;
  while (keepitup == true)
  {
    fgets(Title, 120, InFile);
    Title[strlen(Title)-2] = '\0';
    sscanf(Title, "%d", &NumChose);
    if (NumChose != 0)
    {
      printf("------ xx %d %s Test Case --\n", NumChose, Title+7);
      fprintf(OutFile, 
              "-- Test Case -------- xx %d %s -----------------------\n",
              NumChose, Title+7);
  
      switch (NumChose)
      {
        /* -------------------- Tests for Math ---------------------- */
        case 2:
          TestGETPART();
          break;
        case 3:
          TestGETPARTL();
          break;
        case 4:
          TestGetPartR();
          break;
        case 5:
          TestFACTORIAL();
          break;
        case 6:
          TestBINOMIAL();
          break;
        case 7:
          TestMIN();
          break;
        case 8:
          TestMAX();
          break;
        case 9:
          TestPLANE();
          break;
        case 10:
          TestARCCOSH();
          break;
        case 11:
          TestSINH();
          break;
        case 12:
          TestARCSINH();
          break;
        case 13:
          TestARCTANH();
          break;
        case 14:
          TestDOT();
          break;
        case 15:
          TestCROSS();
          break;
        case 16:
          TestMAG();
          break;
        case 17:
          TestNORM();
          break;
        case 18:
          TestROT1();
          break;
        case 19:
          TestROT2();
          break;
        case 20:
          TestROT3();
          break;
        case 21:
          TestADDVEC();
          break;
        case 22:
          TestADD3VEC();
          break;
        case 23:
          TestLNCOM1();
          break;
        case 24:
          TestLNCOM2();
          break;
        case 25:
          TestLNCOM3();
          break;
        case 26:
          TestPOLYFIT();
          break;
        case 27:
          TestANGLE();
          break;
        case 29:
          TestFACTOR();
          break;
        case 30:
          TestQUADRATIC();
          break;
        case 31:
          TestCUBIC();
          break;
        case 32:
          TestQUARTIC();
          break;
        case 33:
          TestMATSCALE();
          break;
        case 34:
          TestMATMULT();
          break;
        case 35:
          TestMATADD();
          break;
        case 36:
          TestMATSUB();
          break;
        case 37:
          TestMATTRANS();
          break;
        case 38:
          TestMAKEMAT();
          break;
        case 41:
          TestMATINVERSE();
          break;
        case 42:
          TestPRINTMAT();
          break;
        case 43:
          TestFILEPRINTMAT();
          break;
        case 44:
          TestFILEEXPPRINTMAT();
          break;
        case 45:
          TestDETERMINANT();
          break;
        /* -------------------- Tests for Time ---------------------- */
        case 100:
          TestUpCaseSt();
          break;
        case 101:
          TestINITTIME();
          break;
        case 102:
          TestGETINTMON();
          break;
        case 103:
          TestGETINTDAY();
          break;
        case 104:
          TestDAYOFWEEK();
          break;
        case 105:
          TestDAYLIGHTST();
          break;
        case 106:
          TestJULIANDAY();
          break;
        case 107:
          TestJULIANDAYALL();
          break;
        case 108:
          TestDAYS2MDHMS();
          break;
        case 109:
          TestINVJULIANDAY();
          break;
        case 110:
          TestFINDDAYS();
          break;
        case 111:
          TestLSTIME();
          break;
        case 112:
          TestSUNRISESET();
          break;
        case 113:
          TestMOONRISESET();
          break;
        case 114:
          TestHMS_SEC();
          break;
        case 115:
          TestHMS_UT();
          break;
        case 116:
          TestHMS_RAD();
          break;
        case 117:
          TestDMS_RAD();
          break;
        /* -------------------- Tests for 2body --------------------- */
        case 150:
          TestELORB();
          break;
        case 151:
          TestRANDV();
          break;
        case 152:
          TestFINDC2C3();
          break;
        case 153:
          TestNEWTONE();
          break;
        case 154:
          TestNEWTONM();
          break;
        case 155:
          TestNEWTONNU();
          break;
        case 156:
          TestKEPLER();
          break;
        case 157:
          TestFINDTOF();
          break;
        case 158:
          TestIJKTOLATLONA();
          break;
        case 159:
          TestGEOCGEOD();
          break;
        case 160:
          TestSIGHT();
          break;
        case 161:
          TestSUN();
          break;
        case 162:
          TestMOON();
          break;
        case 163:
          TestLIGHT();
          break;
        case 164:
          TestCHECKHITEARTH();
          break;
        case 165:
          TestSATFOV();
          break;
        case 166:
          TestRNGAZ();
          break;
        case 167:
          TestPATH();
          break;
        /* -------------------- Tests for iod ----------------------- */
        case 200:
          TestSITE();
          break;
        case 201:
          TestANGLESLAPLACE();
          break;
        case 202:
          TestANGLESGAUSS();
          break;
        case 203:
          TestRV_RADEC();
          break;
        case 204:
          TestRV_TRADEC();
          break;
        case 205:
          TestRV_RAZEL();
          break;
        case 206:
          TestRV_ELATLON();
          break;
        case 207:
          TestRVSEZ_RAZEL();
          break;
        case 208:
          TestRADEC_ELATLON();
          break;
        case 209:
          TestRADEC_AZEL();
          break;
        case 210:
          TestGIBBS();
          break;
        case 211:
          TestHERRGIBBS();
          break;
        case 212:
          TestLAMBERTUNIV();
          break;
        case 213:
          TestLAMBERTBATTIN();
          break;
        case 214:
          TestTARGET();
          break;
        /* -------------------- Tests for reduc --------------------- */
        case 300:
          TestInitNutation();
          break;
        case 301:
          TestConvTime();
          break;
        case 302:
          TestPrecession();
          break;
        case 303:
          TestNutation();
          break;
        case 304:
          TestSidereal();
          break;
        case 305:
          TestPolarM();
          break;
        case 306:
          TestTrueMean();
          break;
        case 307:
          TestFK4();
          break;
        /* -------------------- Tests for manv --------------------- */
        case 401:
          TestHohmann();
          break;
        case 402:
          TestBiElliptic();
          break;
        case 403:
          TestOneTangent();
          break;
        case 404:
          TestIOnlyChg();
          break;
        case 405:
          TestNodeOnlyChg();
          break;
        case 406:
          TestIandNodeChg();
          break;
        case 407:
          TestMinCombinedPlaneChg();
          break;
        case 408:
          TestCombinedPlaneChg();
          break;
        case 409:
          TestRendezvous();
          break;
        case 410:
          TestNonCoplanarRendz();
          break;
        case 411:
          TestMinCombinedPlaneChg();
          break;
        case 412:
          TestHillsR();
          break;
        case 501:
          TestPKEPLER();
          break;
        case 502:
          TestJ2DragPert();
          break;
        case 503:
          TestPredict();
          break;
        case 504:
          TestInitGravityField();
          break;
        case 505:
          TestLegPoly();
          break;
        case 506:
          TestDeriv();
          break;
        case 507:
          TestPertaccel();
          break;
        case 508:
          TestPderiv();
          break;
        case 509:
          TestRK4();
          break;
        case 510:
          TestRKF45();
          break;
        case 511:
          TestCowell();
          break;
        case 512:
          TestAtmos();
          break;
        case 600:
          TestNonlin();
          break;
        case 601:
          TestFindAtwaAtwb();
          break;
        case 602:
          TestLeastSquares();
          break;
        case 603:
          TestSequential();
          break;
        default:
          break;
      }
    }

    fflush(stdout);
    if (NumChose >= 307)
      keepitup = false;
  }

  fclose(InFile);
  return 0;
}
