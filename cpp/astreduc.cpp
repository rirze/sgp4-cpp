/*     ----------------------------------------------------------------      
                                                                             

                               UNIT ASTREDUC;

                                                                             
    This file contains Astrodynamic procedures and functions to implement    
    reduction calculations. These routines are described in Ch3.             
                                                                             
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

#include "asttime.h"
#include "astreduc.h"

/* Do these as global so they can be accessed
   Each declaration is repeated, but commented out in the procedures */
Real DEpsFK5,DEps96, DPsi96P,DEps96P, DPsiFK5L, DEpsFK5L;
Matrix PrecMat(3, 3) ,NutMat(3, 3);

/******************************************************************************
  Create a file class for manipulating time records.
******************************************************************************/
typedef enum ctimefilestat {CTF_OK, CTF_CLOSE_ERR, CTF_OPEN_ERR} CTF_Stat;
typedef enum ctimefiletype {CTF_TEXT, CTF_BIN, CTR_WRITE, CTF_READ, CTF_APPEND}
             CTF_Type;
class CTimeFile
{
  private:
    char  filename[80];
    FILE *theFile;
    long  curpos;
    long  filesize;
    bool  infile;
    bool  outfile;
    bool  open;
  public:
    // Constructor/Destructor
    CTimeFile(void);
   ~CTimeFile(void);
    // User methods 
    CTF_Stat Close(void);
    long     FileSize(void);
    bool     isOpen(void);
    CTF_Stat Open(char *, CTF_Type);
};

// Constructor/Destructor
CTimeFile::CTimeFile(void)
{
  filename[0] = '\0';
  theFile     = NULL;
  curpos      = 0;
  filesize    = 0;
  infile      = false;
  outfile     = false;
  open        = false;
}

CTimeFile::~CTimeFile(void)
{
  if ((open == true) && (theFile != NULL))
    fclose(theFile);
}

// User methods 
CTF_Stat CTimeFile::Close(void)
{
  return CTF_OK;
}

long CTimeFile::FileSize(void)
{
  return 0;
}

bool CTimeFile::isOpen(void)
{
  if ((open == true) && (theFile != NULL))
    return true;
  else
    return false;
}

CTF_Stat CTimeFile::Open(char *fn, CTF_Type ct)
{
  return CTF_OK;
}

/******************************************************************************
  Implementation functions
******************************************************************************/
void ConvTime
    (
      char *FileN1, UINT Year, UINT Mon, UINT Day, UINT Hr, UINT Min, Real Sec,
      UINT TimeZone, char TypeUTIn, Real& DUT1, Real& DAT, Real& xp, Real& yp,
      Real& UT1, Real& TUT1, Real& JDUT1, Real& UTC, Real& JDUTC, Real& TAI,
      Real& TT, Real& TTT, Real& JDTT, Real& TDB, Real& TTDB, Real& JDTDB,
      Real& DDPsi, Real& DDEps, Real& XLOD, char *Error
    )
{
  const Real Deg2Rad = PI / 180.0;
  const Real Rad     = 180.0 / PI;

  FILE    *TimeFile;
  TimeRec  CurrTimeRec;
  Integer  LocalHr, HrTemp, MinTemp;
  long     curpos, endpos, filesize;
  SINT     temp;
  Real     Ratio, SecTemp, ME, Temp, JD;

  strcpy(Error, "ok");

  if ((TimeFile = fopen(FileN1, "r")) == NULL)
  {
    printf("Can't open time file %s\n", FileN1);
    exit(0);
  }
  JulianDay(Year, Mon, Day, 0, 0, 0.0, JD);

  Ratio = (Hr + Min / 60.0 + Sec / 3600.0) / 24.0;
  curpos = ftell(TimeFile);
  temp = fseek(TimeFile, 0, SEEK_END);
  endpos = ftell(TimeFile);
  filesize = (endpos - curpos + 1) / sizeof(TimeRec);
  temp = fseek(TimeFile, 0, SEEK_SET);
//  if ((SINT(JD - 2442917.5) > FileSize(TimeFile)) || (JD < 2442917.5))
  if ((SINT(JD - 2442917.5) > filesize) || (JD < 2442917.5))
    strcpy(Error, "DateOutRng");
  else
    temp = fseek(TimeFile, SINT(JD - 2442917.5) * sizeof(TimeRec), SEEK_SET);

  if (strcmp(Error, "ok") == 0)
  {
    fscanf(TimeFile, "%i %x %x %f %x %f %f %f %f %f",
                      &CurrTimeRec.Year, &CurrTimeRec.Mon, &CurrTimeRec.Day, 
                      &CurrTimeRec.DUT1, &CurrTimeRec.DAT, &CurrTimeRec.xp, 
                      &CurrTimeRec.yp, &CurrTimeRec.XLOD, &CurrTimeRec.DDPsi, 
                      &CurrTimeRec.DDEps);
    DUT1  = CurrTimeRec.DUT1; // sec
    DAT   = CurrTimeRec.DAT;
    xp    = CurrTimeRec.xp;   // "
    yp    = CurrTimeRec.yp;
    DDPsi = CurrTimeRec.DDPsi; // "
    DDEps = CurrTimeRec.DDEps;
    XLOD  = CurrTimeRec.XLOD;

    if (Show == 'I')
      if (FileOut != NULL)
      {
        fprintf(FileOut, "--------------- Before Interpolation\n");
        fprintf(FileOut, "Time is input in UT %c\n", TypeUTIn);
        fprintf(FileOut, "DUT1  %14.9f s\n", DUT1);
        fprintf(FileOut, "DAT   %14.9f\n", DAT);
        fprintf(FileOut, "xp    %14.9f\"\n", xp);
        fprintf(FileOut, "yp    %14.9f\n", yp);
        fprintf(FileOut, "DDPsi %14.9f\"\n", DDPsi);
        fprintf(FileOut, "DDEps %14.9f\n", DDEps);
        fprintf(FileOut, "XLOD  %14.9f\n", XLOD);
      }
 
      /* -- Do simple linear interpolation with next day to find values -- */
      if (FileOut != NULL)
        fprintf(FileOut, " NO INTERPOLATION\n");

    fscanf(TimeFile, "%i %x %x %f %x %f %f %f %f %f",
                      &CurrTimeRec.Year, &CurrTimeRec.Mon, &CurrTimeRec.Day, 
                      &CurrTimeRec.DUT1, &CurrTimeRec.DAT, &CurrTimeRec.xp, 
                      &CurrTimeRec.yp, &CurrTimeRec.XLOD, &CurrTimeRec.DDPsi, 
                      &CurrTimeRec.DDEps);
    DUT1  = DUT1 + Ratio * (CurrTimeRec.DUT1 - DUT1);
    DAT   = DAT  + Ratio * (CurrTimeRec.DAT  - DAT);
    xp    = xp   + Ratio * (CurrTimeRec.xp   - xp);
    yp    = yp   + Ratio * (CurrTimeRec.yp   - yp);
    DDPsi = DDPsi+ Ratio * (CurrTimeRec.DDPsi - DDPsi);
    DDEps = DDEps+ Ratio * (CurrTimeRec.DDEps - DDEps);
    XLOD  = XLOD + Ratio * (CurrTimeRec.XLOD - XLOD);

    xp    = xp / (3600.0 * Rad); // " to rad
    yp    = yp / (3600.0 * Rad);
  }
  else
  {
    /* ---- Set to zero, but still calculate ---- */
    DUT1  = 0.0;
    DAT   = 0.0;
    xp    = 0.0;
    yp    = 0.0;
    DDPsi = 0.0;
    DDEps = 0.0;
    XLOD  = 0.0;
    if (FileOut != NULL)
      fprintf(FileOut, "%s %18.2f too far\n", Error, JD);
  }
  fclose(TimeFile);

  /* -------------------- Start IF UT1 is known ------------------- */
  LocalHr = TimeZone + Hr;
  if (TypeUTIn == '1')
  {
    HMS_Sec(LocalHr, (Integer)Min, Sec, TOO, UT1);
    JulianDay(Year, Mon, Day, LocalHr, Min, Sec, JDUT1);
    TUT1 = (JDUT1 - 2451545.0) / 36525.0;

    UTC  = UT1 - DUT1;
    HMS_Sec(HrTemp, MinTemp, SecTemp, FROM, UTC);
    JulianDay(Year,Mon,Day, HrTemp, MinTemp, SecTemp, JDUTC);
  }
  else
  {
    /* ------------------ Start IF UTC is known ------------------- */
    HMS_Sec(LocalHr,(Integer) Min, Sec, TOO, UTC);
    JulianDay(Year, Mon, Day, LocalHr, Min, Sec, JDUTC);

    UT1 = UTC + DUT1;
    HMS_Sec(HrTemp, MinTemp, SecTemp, FROM, UT1);
    JulianDay(Year, Mon, Day, HrTemp, MinTemp, SecTemp, JDUT1);
    TUT1 = (JDUT1 - 2451545.0)/ 36525.0 ;
  }

  TAI = UTC + DAT;

  TT = TAI + 32.184; // SEC
  HMS_Sec(HrTemp, MinTemp, SecTemp, FROM, TT);
  JulianDay(Year,Mon,Day, HrTemp, MinTemp, SecTemp, JDTT);
  TTT = (JDTT - 2451545.0)/ 36525.0;

  ME = 357.5277233 + 35999.05034 * TTT;  // approx - should do with TTDB
  ME = Mod(ME, 360.0) * Deg2Rad;
  TDB = TT + 0.001658 * sin(ME) + 0.00001385 * sin(2.0 * ME);
  HMS_Sec(HrTemp, MinTemp, SecTemp, FROM, TDB);
  JulianDay(Year, Mon, Day, HrTemp, MinTemp, SecTemp, JDTDB);
  TTDB = (JDTDB - 2451545.0) / 36525.0;

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------ CONVTIME ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "Year    %4d\n", Year);
      fprintf(FileOut, "Month   %4d\n", Mon);
      fprintf(FileOut, "Day     %4d\n", Day);
      fprintf(FileOut, "Hr      %4d\n", Hr);
      fprintf(FileOut, "Min     %4d\n", Min);
      fprintf(FileOut, "Sec     %14.2f\n", Sec);
      fprintf(FileOut, "Time is input in UT%c\n", TypeUTIn);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "Output vars --------\n");
      fprintf(FileOut, "DUT1  %14.9fs\n", DUT1);
      fprintf(FileOut, "DAT   %14.9f\n", DAT);
      fprintf(FileOut, "xp    %14.9f\"\n", xp * 3600.0 * Rad);
      fprintf(FileOut, "yp    %14.9f\n", yp * 3600.0 * Rad);
      fprintf(FileOut, "DDPsi %14.9f\"\n", DDPsi);
      fprintf(FileOut, "DDEps %14.9f\n", DDEps);
      fprintf(FileOut, "XLOD  %14.9f\n", XLOD);
      fprintf(FileOut, "UT1   %14.9fs\n", UT1);
      fprintf(FileOut, "UTC   %14.9fs\n", UTC);
      fprintf(FileOut, "TAI   %14.9fs\n", TAI);
      fprintf(FileOut, "TT    %14.9fs\n", TT);
      fprintf(FileOut, "TDB   %14.9fs\n", TDB);
      fprintf(FileOut, "TUT1  %14.9fs\n", TUT1);
      fprintf(FileOut, "TTT   %14.9fs\n", TTT);
      fprintf(FileOut, "TTDB  %14.9fs\n", TTDB);
      fprintf(FileOut, "JDUT1 %14.9fs\n", JDUT1);
      fprintf(FileOut, "JDUTC %14.9fs\n", JDUTC);
      fprintf(FileOut, "JDTT  %14.9fs\n", JDTT);
      fprintf(FileOut, "JDTDB %14.9fs\n", JDTDB);
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE FK4
|
|  This PROCEDURE converts vectors between the B1950 and J2000 epochs.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rJ2000      - Initial J2000 Position vector       ER
|    vJ2000      - Initial J2000 Velocity vector      ER/TU
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    rECEF       - Position vector Earth Centered Earth Fixed    ER
|    vECEF       - Velocity vector Earth Centered Earth Fixed  ER/TU
|
|  Locals        :
|    r11,r12,r13 - Components of rotation matrix
|    r21,r22,r23 - Components of rotation matrix
|    r31,r32,r33 - Components of rotation matrix
|
|  Coupling      :
|    MAG         - Vector magnitude
|
|  References    :
|    Vallado       2001, 227-228
|
  ----------------------------------------------------------------------------*/
void FK4
    (
      Vector& rJ2000, Vector& vJ2000, Direction dir, Vector& rFK4, Vector& vFK4
    )
{
  Real r11, r12, r13, r21, r22, r23, r31, r32, r33;

  r11 =  0.9999256794956877;
  r12 = -0.0111814832204662;
  r13 = -0.0048590038153592;

  r21 =  0.0111814832391717;
  r22 =  0.9999374848933135;
  r23 = -0.0000271625947142;

  r31 =  0.0048590037723143;
  r32 = -0.0000271702937440;
  r33 =  0.9999881946043742;

  if (dir == TOO)
  {
    rFK4.Set(r11*rJ2000.Get(1) + r21*rJ2000.Get(2) + r31*rJ2000.Get(3), 1);
    rFK4.Set(r12*rJ2000.Get(1) + r22*rJ2000.Get(2) + r32*rJ2000.Get(3), 2);
    rFK4.Set(r13*rJ2000.Get(1) + r23*rJ2000.Get(2) + r33*rJ2000.Get(3), 3);
    vFK4.Set(r11*vJ2000.Get(1) + r21*vJ2000.Get(2) + r31*vJ2000.Get(3), 1);
    vFK4.Set(r12*vJ2000.Get(1) + r22*vJ2000.Get(2) + r32*vJ2000.Get(3), 2);
    vFK4.Set(r13*vJ2000.Get(1) + r23*vJ2000.Get(2) + r33*vJ2000.Get(3), 3);
  }
  else
  {
    rJ2000.Set(r11*rFK4.Get(1) + r12*rFK4.Get(2) + r13*rFK4(3), 1);
    rJ2000.Set(r21*rFK4.Get(1) + r22*rFK4.Get(2) + r23*rFK4(3), 2);
    rJ2000.Set(r31*rFK4.Get(1) + r32*rFK4.Get(2) + r33*rFK4(3), 3);
    vJ2000.Set(r11*vFK4.Get(1) + r12*vFK4.Get(2) + r13*vFK4(3), 1);
    vJ2000.Set(r21*vFK4.Get(1) + r22*vFK4.Get(2) + r23*vFK4(3), 2);
    vJ2000.Set(r31*vFK4.Get(1) + r32*vFK4.Get(2) + r33*vFK4(3), 3);
  }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE IAU2000GCRF
|
|  This PROCEDURE converts position and velocity vectors between the
|    Inertial Geocentric Celestial Reference Frame (GCRF), and the
|    Body Fixed, International Terrestrial Reference Frame (ITRF).
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rGCRF       - Initial GCRF ECI Position      ER
|    vGCRF       - Initial GCRF ECI Velocity      ER/TU
|    TTT         - Julian Centuries of TT         centuries
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    vITRF       - Initial ITRF Velocity          ER/TU
|    rITRF       - Initial ITRF Position          ER
|
|  Locals        :
|    TTT2        - TTT squared
|    TTT3        - TTT cubed
|    Zeta        - PRECESSION ANGLE               rad
|    z           - PRECESSION ANGLE               rad
|    Theta       - PRECESSION ANGLE               rad
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|
|  References    :
|    Vallado       2001,
|
  ----------------------------------------------------------------------------*/
void IAU2000GCRF
    (
      Vector& rGCRF, Vector& vGCRF, Direction dir, Vector& rITRF, Vector& vITRF,      Real JDUT1, Real TTT, Real xp, Real yp, IAr5x106 IAr00, RAr6x106 RAr00
    )
{
  const Real Deg2Rad       =    PI / 180.0;
  const Real OmegaEarth    =    0.05883359980154919;
  const Real RadiusEarthKM = 6378.1363;
  const Real VKmPerSec     =    7.905366149846;

  Vector v3a(3), r1(3), r2(3), r3(3), v1(3), v2(3), v3(3), Temp(3), Temp1(3);
  Real   zeta, z, theta, TTT2, TTT3;

  SINT deg, MIN;
  Real LonMer, LonVen, LonEar, LonMar, LonJup, LonSat, PrecRate, SA,
       DeltaPsi00, DeltaEps00, X, Y, s, a, aa, ac, AccDispl, deltapsi, deltaeps,
       Eps,TrueEps, SEC, TempVal, TTT4, rr, l, l1, F, D, Omega;
  Real PN[3][3];

  /* ------- Determine coefficients for IAU 2000(xx) NUTATION Theory -------- */
  TTT2 = TTT * TTT;
  TTT3 = TTT2 * TTT;
  TTT4 = TTT3 * TTT;
  // GetTime( th,tm,ts,ts1);
  // Tem1Time = th*3600 + tm*60 + ts + ts1/100;

  Eps = 23.439291 - 0.0130042 * TTT - 0.000000164 * TTT2 + 0.000000504 * TTT3;
  Eps = Mod(Eps, 360.0);
  Eps = Eps * Deg2Rad;

  /* -------------- Evaluate Brown's fundamental arguments ---------------- */
  rr    = 360.0; // deg
  l     = 134.96340251 + (1325 * rr  + 198.8675605) * TTT + 0.0088553 * TTT2 +
          0.000014343 * TTT3 - 0.00000006797 * TTT4;
  l1    =  357.52910918 + (99 * rr  + 359.0502911) * TTT - 0.0001537 * TTT2 -
          0.000000038 * TTT3 - 0.00000000319 * TTT4;
  F     =   93.27209062 + (1342 * rr  +  82.0174577) * TTT - 0.0035420 * TTT2 +
          0.000000288 * TTT3 + 0.00000000116 * TTT4;
  D     =  297.85019547 + (1236 * rr  + 307.1114469) * TTT - 0.0017696 * TTT2 +
          0.000001831 * TTT3 - 0.00000000880 * TTT4;
  Omega =  125.04455501 - (5 * rr  + 134.1361851) * TTT + 0.0020756 * TTT2 +
          0.000002139 * TTT3 - 0.00000001650 * TTT4;

  l     = Mod(l, 360.0)     * Deg2Rad;
  l1    = Mod(l1, 360.0)    * Deg2Rad;
  F     = Mod(F, 360.0)     * Deg2Rad;
  D     = Mod(D, 360.0)     * Deg2Rad;
  Omega = Mod(Omega, 360.0) * Deg2Rad;

  DeltaPsi00 = 0.0;
  DeltaEps00 = 0.0;

  for (UINT i = 106; i > 0; i--)
  {
    TempVal = IAr00[0][i-1] * l + IAr00[1][i-1] * l1 + IAr00[2][i-1] * F +
              IAr00[3][i-1] * D + IAr00[4][i-1] * Omega;
    DeltaPsi00 = DeltaPsi00 + (RAr00[0][i-1] + 
                 RAr00[1][i-1] * TTT) * sin(TempVal) + 
                 RAr00[4][i-1] * TTT * cos(TempVal);
    DeltaEps00 = DeltaEps00 + (RAr00[2][i-1] + 
                 RAr00[3][i-1] * TTT) * cos(TempVal) +
                 RAr00[5][i-1] * TTT * sin(TempVal);
  }
  TrueEps = Eps + DeltaEps00;

  /* ------------------------ Find Stellar Angle ------------------- */
  SA = 2.0 * PI * (0.779057273264 + 1.00273781191135448 * (JDUT1 - 2451545.0));
  TempVal = SA;
  SA = Mod(SA, 2.0*PI);  // Rad

  X = 2004.3109 * TTT - 0.42665 * TTT2 - 0.198656 * TTT3 + 0.0000140 * TTT4 +
      0.00006 * TTT2 * cos(Omega) + 0.00204 * TTT2 * sin(Omega) + 
      0.00016 * TTT2 * cos( 2*(F - D + Omega));
  X = X * Deg2Rad / 3600.0 + sin(TrueEps) * DeltaPsi00; // " to rad
  Y = -0.00013 - 22.40992 * TTT2 + 0.001836 * TTT3 + 0.0011130 * TTT4 -
      0.00231 * TTT2 * cos(Omega) - 0.00014 * TTT2 * cos(2 * (F - D + Omega));
  Y = Y * Deg2Rad / 3600.0 + DeltaEps00; // " to rad
  
  a = 0.5 + 0.125 * (X * X + Y * Y);  // "
  s = -X * Y * 0.5 + 0.00385 * TTT - 0.07259 * TTT3 - 0.00264 * sin(Omega) -
      0.00006 * sin(2 * Omega) + 0.00074 * TTT2 * sin(Omega) +
      0.0006 * TTT2 * sin(2 * (F - D + Omega));

  a = a * Deg2Rad / 3600.0; // " to rad
  s = s * Deg2Rad / 3600.0; // " to rad

  aa = 0.12; // "
  ac = 0.26;
  AccDispl = 0.0015*(ac*ac / 1.2 + aa*aa) * TTT * Deg2Rad / 3600.0;  // " to rad

  PN[0][0] = 1.0 - a * X * X;
  PN[0][1] = -a * X * Y;
  PN[0][2] = X;
  PN[1][0] = PN[0][1];
  PN[1][1] = 1.0 - a * Y * Y;
  PN[1][2] = Y;
  PN[2][0] = -X;
  PN[2][1] = -Y;
  PN[2][2] = 1.0 - a * (X * X + Y * Y);

  if (dir == TOO)
  {
    r1.Set(PN[0][0] * rGCRF.Get(1) + PN[0][1] * rGCRF.Get(2) + 
           PN[0][2] * rGCRF.Get(3), 1);
    r1.Set(PN[1][0] * rGCRF.Get(1) + PN[1][1] * rGCRF.Get(2) + 
           PN[1][2] * rGCRF.Get(3), 2);
    r1.Set(PN[2][0] * rGCRF.Get(1) + PN[2][1] * rGCRF.Get(2) + 
           PN[2][2] * rGCRF.Get(3), 3);
    r2    = r1.Rot3(s);
    r3    = r2.Rot3(SA);
    Temp  = r3.Rot3(AccDispl);
    Temp1 = Temp.Rot1(-yp);
    rITRF = Temp1.Rot2(-xp);

    v1.Set(PN[0][0] * vGCRF.Get(1) + PN[0][1] * vGCRF.Get(2) + 
           PN[0][2] * vGCRF.Get(3), 1);
    v1.Set(PN[1][0] * vGCRF.Get(1) + PN[1][1] * vGCRF.Get(2) + 
           PN[1][2] * vGCRF.Get(3), 2);
    v1.Set(PN[2][0] * vGCRF.Get(1) + PN[2][1] * vGCRF.Get(2) + 
           PN[2][2] * vGCRF.Get(3), 3);
    v2    = v1.Rot3(s);
    v3a   = v2.Rot3(SA);
    Temp  = v3.Rot3(AccDispl);
    Temp1 = Temp.Rot1(-yp);
    vITRF = Temp1.Rot2(-xp);
  }
  else
  {
    Temp  = rITRF.Rot2(xp);
    Temp1 = Temp.Rot1(yp);
    r3    = Temp1.Rot3(-AccDispl);
    r2    = r3.Rot3(-SA);
    r1    = r2.Rot3(-s);
    rGCRF.Set(PN[0][0] * r1.Get(1) + PN[1][0] * r1.Get(2) + 
              PN[2][0] * r1.Get(3), 1);
    rGCRF.Set(PN[0][1] * r1.Get(1) + PN[1][1] * r1.Get(2) + 
              PN[2][1] * r1.Get(3), 2);
    rGCRF.Set(PN[0][2] * r1.Get(1) + PN[1][2] * r1.Get(2) + 
              PN[2][2] * r1.Get(3), 3);

    Temp  = vITRF.Rot2(xp);
    Temp1 = Temp.Rot1(yp);
    v3    = Temp1.Rot3(-AccDispl);
    v3a   = v3;
    v3a.Set(v3.Get(1) - r3.Get(2) * OmegaEarth, 1);
    v3a.Set(v3.Get(2) - r3.Get(1) * OmegaEarth, 2);
    v2 = v3a.Rot3(-SA);
    v1 = v2.Rot3(-s);
    vGCRF.Set(PN[0][0] * v1.Get(1) + PN[1][0] * v1.Get(2) + 
              PN[2][0] * v1.Get(3), 1);
    vGCRF.Set(PN[0][1] * v1.Get(1) + PN[1][1] * v1.Get(2) + 
              PN[2][1] * v1.Get(3), 2);
    vGCRF.Set(PN[0][2] * v1.Get(1) + PN[1][2] * v1.Get(2) + 
              PN[2][2] * v1.Get(3), 3);

  }
 
  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  IAU2000 ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "JDUT1 %18.12fs\n", JDUT1);
      fprintf(FileOut, "TTT   %14.8s\n", TTT);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "l       %18.12fø %18.12f\n", l / Deg2Rad, l);
      fprintf(FileOut, "l1      %18.12fø %18.12f\n", l1 / Deg2Rad, l1);
      fprintf(FileOut, "f       %18.12fø %18.12f\n", F / Deg2Rad, F);
      fprintf(FileOut, "d       %18.12fø %18.12f\n", D / Deg2Rad, D);
      fprintf(FileOut, "omega   %18.12fø %18.12f\n", Omega / Deg2Rad, Omega);
      DMS_Rad(deg, MIN, SEC, FROM, Eps);
      fprintf(FileOut, "DeltaPsi00 %18.12fø %18.12f\n", 
                        DeltaPsi00 / Deg2Rad, DeltaPsi00);
      fprintf(FileOut, "Eps     %18.12fø %3d%3d%9.5f%18.12f\n", 
                        Eps / Deg2Rad, deg, MIN, SEC, Eps);
      fprintf(FileOut, "DeltaEps00 %18.12fø%18.12f\n", 
                        DeltaEps00 / Deg2Rad, DeltaEps00);
      fprintf(FileOut, "TrueEps %18.12fø\n", TrueEps / Deg2Rad);
      fprintf(FileOut, "x       %18.12fø %18.12f\n", X / Deg2Rad, X);
      fprintf(FileOut, "y       %18.12fø %18.12f\n", Y / Deg2Rad, Y);
      fprintf(FileOut, "a       %18.12fø %18.12f\n", a / Deg2Rad, a);
      fprintf(FileOut, "s       %18.12fø %18.12f\n", s / Deg2Rad, s);
      fprintf(FileOut, "r2 aftSA%18.12f%18.12f%18.12f\n", 
                        r2.Get(1) * RadiusEarthKM, r2.Get(2) * RadiusEarthKM,
                        r2.Get(3) * RadiusEarthKM);
      fprintf(FileOut, "V2      %18.12f%18.12f%18.12f\n",
                       v2.Get(1) * VKmPerSec, v2.Get(2) * VKmPerSec, 
                       v2.Get(3) * VKmPerSec);
      fprintf(FileOut, "SA      %18.12f%18.12f\n", SA / Deg2Rad, SA);
      fprintf(FileOut, "SA      %18.12f%18.12f\n", TempVal / Deg2Rad, TempVal);
      fprintf(FileOut, "R3 aftPM%18.12f%18.12f%18.12f\n",
                        r3.Get(1) * RadiusEarthKM, r3.Get(2) * RadiusEarthKM, 
                        r3.Get(3) * RadiusEarthKM);
      fprintf(FileOut, "V3      %18.12f%18.12f%18.12f\n",
                        v3.Get(1) * RadiusEarthKM, v3.Get(2) * RadiusEarthKM, 
                        v3.Get(3) * RadiusEarthKM);
      fprintf(FileOut, "sp      %18.12fø %18.12f\n", 
                        AccDispl / Deg2Rad, AccDispl);
      fprintf(FileOut, "xp %12.9f yp %12.9f rad\n", xp, yp);
      fprintf(FileOut, "xp %12.9f yp %12.9f rad\n", 
                        xp * 180.0 / PI, yp * 180.0 / PI);
      fprintf(FileOut, "Output vars --------\n");
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE INITNUTATION
|
|  This procedure initializes the nutation matricies needed for reduction
|    calculations. The routine needs the filename of the files as input.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|                                               
|  Inputs          Description                    Range / Units
|    FileN1      - Filename for nutation terms
|    FileN2      - Filename for planetary nutation terms
|
|  Outputs       :
|                - Arrays of data for each technique
|
|  Locals        :
|    Convrt      - Changes milli-arcseconds to degrees
|
|  Coupling      :
|    None        -
|
 -----------------------------------------------------------------------------*/
void InitNutation
    (
      char *FileN1, char *FileN2, IAr5x106 IAr50, RAr4x106 RAr50,
      IAr5x106 IAr80, RAr4x106 RAr80, IAr5x263 IAr96, RAr6x263 RAr96,
      IAr10x112 pIAr96, RAr4x112 pRAr96
    )
{
  Real  Convrt;
  FILE *InFile;
  char  dummy[128];

  if ((InFile = fopen(FileN1, "r")) == NULL)
  {
    printf("Unable to open input file %s\n", FileN1);
    exit(0);
  }

//  /* --------- Read in 2000 IAU Theory of Nutation Terms ---------- */
//  Convrt = 0.0001 / 3600.0; // 0.0001" to deg
//  fscanf(InFile, "%s\n", dummy);  // Header
//  fscanf(InFile, "%s\n", dummy);  // Header
//  for (UINT i = 0; i < 106; i++)
//  {
//    fscanf(InFile, "%f %f %f %f %f %f %f %f %f %f %f", 
//                    IAr00[0][i],IAr00[1][i],IAr00[2][i],IAr00[3][i],IAr00[4][i],
//      RAr00[0][i],RAr00[1][i],RAr00[2][i],RAr00[3][i],RAr00[4][i],RAr00[5][i]);
//    RAr00[0][i] = Convrt * RAr00[0][i];
//    RAr00[1][i] = Convrt * RAr00[1][i];
//    RAr00[2][i] = Convrt * RAr00[2][i];
//    RAr00[3][i] = Convrt * RAr00[3][i];
//    RAr00[4][i] = Convrt * RAr00[4][i];
//    RAr00[5][i] = Convrt * RAr00[5][i];
//  }

  /* --------- Read in 1950 IAU Theory of Nutation Terms ---------- */
  Convrt = 0.0001 / 3600.0; // 0.0001" to deg
  fscanf(InFile, "%s\n", dummy);  // Header
  fscanf(InFile, "%s\n", dummy);  // Header
  for (UINT i = 0; i < 106; i++)
  {
    fscanf(InFile, "%f %f %f %f %f %f %f %f %f %f %f",
                    IAr50[0][i],IAr50[1][i],IAr50[2][i],IAr50[3][i],IAr50[4][i],
                    RAr50[0][i],RAr50[1][i],RAr50[2][i],RAr50[3][i]);
    RAr50[0][i] = Convrt * RAr50[0][i];
    RAr50[1][i] = Convrt * RAr50[1][i];
    RAr50[2][i] = Convrt * RAr50[2][i];
    RAr50[3][i] = Convrt * RAr50[3][i];
  }

  /* --------- Read in 1980 IAU Theory of Nutation Terms ---------- */
  Convrt = 0.0001 / 3600.0; // 0.0001" to deg
  fscanf(InFile, "%s\n", dummy);  // Header
  fscanf(InFile, "%s\n", dummy);  // Header
  for (UINT i = 0; i < 106; i++)
  {
    fscanf(InFile, "%f %f %f %f %f %f %f %f %f %f %f",
                    IAr80[0][i],IAr80[1][i],IAr80[2][i],IAr80[3][i],IAr80[4][i],
                    RAr80[0][i],RAr80[1][i],RAr80[2][i],RAr80[3][i]);
    RAr80[0][i] = Convrt * RAr80[0][i];
    RAr80[1][i] = Convrt * RAr80[1][i];
    RAr80[2][i] = Convrt * RAr80[2][i];
    RAr80[3][i] = Convrt * RAr80[3][i];
  }

  /* --------- Read in 1996 IAU Theory of Nutation Terms ---------- */
  Convrt = 0.000001 / 3600.0; // 0.000001" to deg
  fscanf(InFile, "%s\n", dummy);  // Header
  fscanf(InFile, "%s\n", dummy);  // Header
  for (UINT i = 0; i < 263; i++)
  {
    fscanf(InFile, "%f %f %f %f %f %f %f %f %f %f %f",
                    IAr96[0][i],IAr96[1][i],IAr96[2][i],IAr96[3][i],IAr96[4][i],
      RAr96[0][i],RAr96[1][i],RAr96[2][i],RAr96[3][i],RAr96[4][i],RAr96[5][i]);
    RAr96[0][i] = Convrt * RAr96[0][i];
    RAr96[1][i] = Convrt * RAr96[1][i];
    RAr96[2][i] = Convrt * RAr96[2][i];
    RAr96[3][i] = Convrt * RAr96[3][i];
    RAr96[4][i] = Convrt * RAr96[4][i];
    RAr96[5][i] = Convrt * RAr96[5][i];
  }
  fclose(InFile);

/* ----------- Read in 1996 Planetary Nutation Terms ------------ */
  Convrt = 0.000001 / 3600.0; // 0.000001" to deg
  if ((InFile = fopen(FileN2, "r")) == NULL)
  {
    printf("Unable to open input file %s\n", FileN1);
    exit(0);
  }
  fscanf(InFile, "%s\n", dummy);  // Header
  fscanf(InFile, "%s\n", dummy);  // Header
  for (UINT i = 0; i < 112; i++)
  {
    fscanf(InFile, "%f %f %f %f %f %f %f %f %f %f %f",
           pIAr96[0][i],pIAr96[1][i],pIAr96[2][i],pIAr96[3][i],pIAr96[4][i],
           pRAr96[0][i],pRAr96[1][i],pRAr96[2][i],pRAr96[3][i]);
    pRAr96[0][i] = Convrt * pRAr96[0][i];
    pRAr96[1][i] = Convrt * pRAr96[1][i];
    pRAr96[2][i] = Convrt * pRAr96[2][i];
    pRAr96[3][i] = Convrt * pRAr96[3][i];
  }
  fclose(InFile);

}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE NUT80FK5
|
|  This PROCEDURE converts position and velocity vectors between the
|    Mean Equator Mean Equinox of date (MOD) and the True equator True equinox
|    of date (TOD).  The results take into account the effects of NUTATION using
|    the 1980 IAU Theory of Nutation.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rMOD        - Position vector of date
|                    Mean Equator, Mean Equinox   ER
|    vMOD        - Velocity vector of date
|                    Mean Equator, Mean Equinox   ER/TU
|    TTT         - Julian Centuries of TDB
|    iAr80       - Array of data coefficients
|    rAr80       - Array of data coefficients
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    rTOD        - Position vector of date
|                    True Equator, True Equinox   ER
|    vTOD        - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|    DPsiFK5     - NUTATION ANGLE                 rad
|    TrueEpsFK5  - True obliquity of the ecliptic rad
|    Omega       -                                rad
|
|  Locals        :
|    TTT2        - TTT squared
|    TTT3        - TTT cubed
|    TTT4        - TTT to the fourth
|    Eps         - Mean obliquity of the ecliptic rad
|    l           -                                rad
|    ll          -                                rad
|    F           -                                rad
|    D           -                                rad
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|    REALMOD     -
|
|  References    :
|    Vallado       2001, 216-219, Eq 3-654
|
  ----------------------------------------------------------------------------*/
void Nut80FK5
    (
      Vector& rMOD, Vector& vMOD, Direction dir, Vector& rTOD, Vector& vTOD, 
      Real& DPsiFK5, Real& TrueEpsFK5, Real& Omega, Real TTT, 
      IAr5x106 IAr80, RAr4x106 RAr80
    )
{
  const Real Deg2Rad = PI / 180.0;

  Vector Temp(3), Temp1(3);
  SINT   Deg, Min;
  Real   Dp, De, TTT4, Sec, TempVal, TTT2, TTT3, Eps, rr, l, l1, F, D;

  /* ------- Determine coefficients for IAU 1980 NUTATION Theory ---------- */
  TTT2 = TTT  * TTT;
  TTT3 = TTT2 * TTT;
  TTT4 = TTT2 * TTT2;
  Eps  = -46.8150 * TTT - 0.00059 * TTT2 + 0.001813 * TTT3 + 84381.448; // " 
  Eps  = Mod( Eps / 3600.0, 360.0);  // deg
  Eps  = Eps * Deg2Rad;  // rad

  l     =  134.96340251 + (1717915923.2178*TTT + 31.8792*TTT2 + 
        0.051635 * TTT3 - 0.00024470 * TTT4) / 3600.0; // deg and ", deg result
  l1    =  357.52910918 + (129596581.0481 * TTT - 0.5532 * TTT2 - 
                           0.000136 * TTT3 - 0.00001149*  TTT4) / 3600.0;
  F     =   93.27209062 + (1739527262.8478 * TTT - 12.7512 * TTT2 +
          0.001037 * TTT3 + 0.00000417 * TTT4) / 3600.0;
  D     =  297.85019547 + (1602961601.2090 * TTT - 6.3706 * TTT2 +
          0.006593 * TTT3 - 0.00003169 * TTT4) / 3600.0;
  Omega =  125.04455501 + (-6962890.2665 * TTT + 7.4722 * TTT2 +
          0.007702 * TTT3 - 0.00005939 * TTT4) / 3600.0;

  l     = Mod(l,    360.0) * Deg2Rad; // rad
  l1    = Mod(l1,   360.0) * Deg2Rad;
  F     = Mod(F,    360.0) * Deg2Rad;
  D     = Mod(D,    360.0) * Deg2Rad;
  Omega = Mod(Omega,360.0) * Deg2Rad;

  DPsiFK5  = 0.0;
  DPsiFK5L = 0.0;
  DEpsFK5  = 0.0;
  DEpsFK5L = 0.0;

  for (UINT i = 106; i > 0; i--)
  {
    TempVal = IAr80[0][i-1] * l + IAr80[1][i-1] * l1 + IAr80[2][i-1] * F + 
              IAr80[3][i-1] * D + IAr80[4][i-1] * Omega;
    DPsiFK5 = DPsiFK5 + (RAr80[0][i-1] + RAr80[1][i-1] * TTT) * sin(TempVal);
    DEpsFK5 = DEpsFK5 + (RAr80[2][i-1] + RAr80[3][i-1] * TTT) * cos(TempVal);

    if (i <= 4)
    {
      DPsiFK5L = DPsiFK5L + (RAr80[0][i-1] + RAr80[1][i-1]*TTT) * sin(TempVal);
      DEpsFK5L = DEpsFK5L + (RAr80[2][i-1] + RAr80[3][i-1]*TTT) * cos(TempVal);
    }
  }

  /* --------------- Find NUTATION Parameters --------------------- */
  Dp         = DPsiFK5;
  De         = DEpsFK5;
  DPsiFK5    = Mod(DPsiFK5, 360.0) * Deg2Rad; // rad
  DEpsFK5    = Mod(DEpsFK5, 360.0) * Deg2Rad;
  TrueEpsFK5 = Eps + DEpsFK5;

  if (dir == TOO)
  {
    Temp  = rMOD.Rot1(Eps);
    Temp1 = Temp.Rot3(-DPsiFK5);
    rTOD  = Temp1.Rot1(-TrueEpsFK5);
    Temp  = vMOD.Rot1(Eps);
    Temp1 = Temp.Rot3(-DPsiFK5);
    vTOD  = Temp1.Rot1(-TrueEpsFK5);
  }
  else
  {
    Temp  = rTOD.Rot1(TrueEpsFK5);
    Temp1 = Temp.Rot3(DPsiFK5);
    rMOD  = Temp1.Rot1(-Eps);
    Temp  = vTOD.Rot1(TrueEpsFK5);
    Temp1 = Temp.Rot3(DPsiFK5);
    vMOD  = Temp1.Rot1(-Eps);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  80 FK5 NUTATION ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "TTT     %18.12f\n", TTT);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "l       %18.12fø%18.12f", l  / Deg2Rad, l);
      fprintf(FileOut, "l1      %18.12fø%18.12f", l1 / Deg2Rad, l1);
      fprintf(FileOut, "f       %18.12fø%18.12f", F  / Deg2Rad, F);
      fprintf(FileOut, "d       %18.12fø%18.12f", D  / Deg2Rad, D);
      DMS_Rad(Deg, Min, Sec, FROM, Eps);
      fprintf(FileOut, "Eps     %18.12fø%3d%3d%9.5f%18.12f", 
                        Eps / Deg2Rad, Deg, Min, Sec, Eps);
      fprintf(FileOut, "DEpsFK5 %18.12fø%18.12f", DEpsFK5 / Deg2Rad, DEpsFK5);
      fprintf(FileOut, "Output vars --------\n");
      fprintf(FileOut, "DPsiFK5 %18.12fø%18.12f", DPsiFK5 / Deg2Rad, DPsiFK5);
      DMS_Rad(Deg, Min, Sec, FROM, TrueEpsFK5);
      fprintf(FileOut, "TrueEpsFK5 %18.12fø%3d%3d%9.5f",
                        TrueEpsFK5 / Deg2Rad, Deg, Min, Sec);
      fprintf(FileOut, "omega   %18.12fø%18.12f", Omega / Deg2Rad, Omega);
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE NUT80GCRF
|
|  This PROCEDURE converts position and velocity vectors between the GCRF
|    Mean Equator Mean Equinox of date (MOD) and the True equator True equinox
|    of date (TOD).  The results take into account the effects of NUTATION using
|    the 1980 IAU Theory of Nutation.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rMOD        - Position vector of date
|                    Mean Equator, Mean Equinox   ER
|    vMOD        - Velocity vector of date
|                    Mean Equator, Mean Equinox   ER/TU
|    TTT         - Julian Centuries of TDB
|    DdPsi       - Delta NUTATION ANGLE            "
|    DdEps       - Delta obliquity of the ecliptic "
|    iAr80       - Array of data coefficients
|    rAr80       - Array of data coefficients
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    rTOD        - Position vector of date
|                    True Equator, True Equinox   ER
|    vTOD        - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|    DPsiFK5     - NUTATION ANGLE                 rad
|    TrueEpsFK5  - True obliquity of the ecliptic rad
|    Omega       -                                rad
|
|  Locals        :
|    TTT2        - TTT squared
|    TTT3        - TTT cubed
|    TTT4        - TTT to the fourth
|    Eps         - Mean obliquity of the ecliptic rad
|    l           -                                rad
|    ll          -                                rad
|    F           -                                rad
|    D           -                                rad
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|    REALMOD     -
|
|  References    :
|    Vallado       2001,
|
  ----------------------------------------------------------------------------*/
void Nut80GCRF
    (
      Vector& rMOD, Vector& vMOD, Direction dir, Vector& rTOD, Vector& vTOD, 
      Real& DPsiFK5, Real& TrueEpsFK5, Real& Omega, Real TTT, Real DDPsi, 
      Real DDEps, IAr5x106 IAr80, RAr4x106 RAr80
    )
{
  const Real Deg2Rad = PI / 180.0;

  Vector Temp(3), Temp1(3);
  SINT   Deg, Min;
  Real   Dp, De, TTT4, Sec, TempVal, TTT2, TTT3, Eps, rr, l, l1, F, D;

  /* ------- Determine coefficients for IAU 1980 NUTATION Theory ---------- */
  TTT2 = TTT  * TTT;
  TTT3 = TTT2 * TTT;
  TTT4 = TTT2 * TTT2;
  Eps  = -46.8150 * TTT - 0.00059 * TTT2 + 0.001813 * TTT3 + 84381.448; // " 
  Eps  = Mod( Eps / 3600.0, 360.0);  // deg
  Eps  = Eps * Deg2Rad;  // rad

  l     =  134.96340251 + (1717915923.2178*TTT + 31.8792*TTT2 + 
        0.051635 * TTT3 - 0.00024470 * TTT4) / 3600.0; // deg and ", deg result
  l1    =  357.52910918 + (129596581.0481 * TTT - 0.5532 * TTT2 - 
                           0.000136 * TTT3 - 0.00001149*  TTT4) / 3600.0;
  F     =   93.27209062 + (1739527262.8478 * TTT - 12.7512 * TTT2 +
          0.001037 * TTT3 + 0.00000417 * TTT4) / 3600.0;
  D     =  297.85019547 + (1602961601.2090 * TTT - 6.3706 * TTT2 +
          0.006593 * TTT3 - 0.00003169 * TTT4) / 3600.0;
  Omega =  125.04455501 + (-6962890.2665 * TTT + 7.4722 * TTT2 +
          0.007702 * TTT3 - 0.00005939 * TTT4) / 3600.0;

  l     = Mod(l,    360.0) * Deg2Rad; // rad
  l1    = Mod(l1,   360.0) * Deg2Rad;
  F     = Mod(F,    360.0) * Deg2Rad;
  D     = Mod(D,    360.0) * Deg2Rad;
  Omega = Mod(Omega,360.0) * Deg2Rad;

  DPsiFK5  = 0.0;
  DEpsFK5  = 0.0;

  for (UINT i = 106; i > 0; i--)
  {
    TempVal = IAr80[0][i-1] * l + IAr80[1][i-1] * l1 + IAr80[2][i-1] * F + 
              IAr80[3][i-1] * D + IAr80[4][i-1] * Omega;
    DPsiFK5 = DPsiFK5 + (RAr80[0][i-1] + RAr80[1][i-1] * TTT) * sin(TempVal);
    DEpsFK5 = DEpsFK5 + (RAr80[2][i-1] + RAr80[3][i-1] * TTT) * cos(TempVal);
  }

  /* --------------- Find NUTATION Parameters --------------------- */
  Dp         = DPsiFK5;
  De         = DEpsFK5;
  DPsiFK5    = DPsiFK5 + DDPsi / 3600.0; // deg
  DEpsFK5    = DEpsFK5 + DDEps / 3600.0;
  DPsiFK5    = Mod(DPsiFK5, 360.0) * Deg2Rad; // rad
  DEpsFK5    = Mod(DEpsFK5, 360.0) * Deg2Rad;
  TrueEpsFK5 = Eps + DEpsFK5;

  if (dir == TOO)
  {
    Temp  = rMOD.Rot1(Eps);
    Temp1 = Temp.Rot3(-DPsiFK5);
    rTOD  = Temp1.Rot1(-TrueEpsFK5);
    Temp  = vMOD.Rot1(Eps);
    Temp1 = Temp.Rot3(-DPsiFK5);
    vTOD  = Temp1.Rot1(-TrueEpsFK5);
  }
  else
  {
    Temp  = rTOD.Rot1(TrueEpsFK5);
    Temp1 = Temp.Rot3(DPsiFK5);
    rMOD  = Temp1.Rot1(-Eps);
    Temp  = vTOD.Rot1(TrueEpsFK5);
    Temp1 = Temp.Rot3(DPsiFK5);
    vMOD  = Temp1.Rot1(-Eps);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  80 GCRF NUTATION ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "TTT     %18.12f\n", TTT);
      fprintf(FileOut, "DDPsi   %18.12f\n", DDPsi);
      fprintf(FileOut, "DDEps   %18.12f\n", DDEps);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "l       %18.12fø%18.12f", l  / Deg2Rad, l);
      fprintf(FileOut, "l1      %18.12fø%18.12f", l1 / Deg2Rad, l1);
      fprintf(FileOut, "f       %18.12fø%18.12f", F  / Deg2Rad, F);
      fprintf(FileOut, "d       %18.12fø%18.12f", D  / Deg2Rad, D);
      DMS_Rad(Deg, Min, Sec, FROM, Eps);
      fprintf(FileOut, "Eps     %18.12fø%3d%3d%9.5f%18.12f", 
                        Eps / Deg2Rad, Deg, Min, Sec, Eps);
      fprintf(FileOut, "DPsi106 %18.12fø%18.12f", Dp, Dp * Deg2Rad);
      fprintf(FileOut, "DEps106 %18.12fø%18.12f", De, De * Deg2Rad);
      fprintf(FileOut, "DEpsFK5 %18.12fø%18.12f", DEpsFK5 / Deg2Rad, DEpsFK5);
      fprintf(FileOut, "Output vars --------\n");
      fprintf(FileOut, "DPsiFK5 %18.12fø%18.12f", DPsiFK5 / Deg2Rad, DPsiFK5);
      DMS_Rad(Deg, Min, Sec, FROM, TrueEpsFK5);
      fprintf(FileOut, "TrueEpsFK5 %18.12fø%3d%3d%9.5f",
                        TrueEpsFK5 / Deg2Rad, Deg, Min, Sec);
      fprintf(FileOut, "omega   %18.12fø%18.12f", Omega / Deg2Rad, Omega);
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE NUT96FK5
|
|  This PROCEDURE converts position and velocity vectors between the
|    Mean Equator Mean Equinox of date (MOD) and the True equator True equinox
|    of date (TOD).  The results take into account the effects of NUTATION using
|    the 1996 IAU Theory of Nutation.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rMOD        - Position vector of date
|                    Mean Equator, Mean Equinox   ER
|    vMOD        - Velocity vector of date
|                    Mean Equator, Mean Equinox   ER/TU
|    TTT         - Julian Centuries of TDB
|    iAr96       - Array of data coefficients
|    rAr96       - Array of data coefficients
|    piAr96      - Array of data coefficients planetary
|    prAr96      - Array of data coefficients
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    rTOD        - Position vector of date
|                    True Equator, True Equinox   ER
|    vTOD        - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|    DPsi96      - NUTATION ANGLE                 rad
|    TrueEps96   - True obliquity of the ecliptic rad
|    Omega       -                                rad
|
|  Locals        :
|    TTT2        - TTT squared
|    TTT3        - TTT cubed
|    TTT4        - TTT to the fourth
|    Eps         - Mean obliquity of the ecliptic rad
|    l           -                                rad
|    ll          -                                rad
|    F           -                                rad
|    D           -                                rad
|    DEps96      - Change in obliquity            rad
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|    REALMOD     -
|
|  References    :
|    Vallado       2001,
|
  ----------------------------------------------------------------------------*/
void Nut96FK5
    (
      Vector& rMOD, Vector& vMOD, Direction dir, Vector& rTOD, Vector& vTOD, 
      Real& DPsi96, Real& TrueEps96, Real& Omega, Real TTT, 
      IAr5x263 IAr96, RAr6x263 RAr96, IAr10x112 pIAr96, RAr4x112 pRAr96
    )
{
  const Real Deg2Rad = PI / 180.0;

  Vector Temp(3), Temp1(3);
  SINT   Deg, Min;
  Real   LonMer, LonVen, LonEar, LonMar, LonJup, LonSat, PrecRate,
         Dp, De, Sec, TempVal, TTT2, TTT3, TTT4, Eps, rr, l, l1, F, D;

  /* ------- Determine coefficients for IAU 1996 NUTATION Theory ---------- */
  TTT2 = TTT  * TTT;
  TTT3 = TTT2 * TTT;
  TTT4 = TTT2 * TTT2;
  Eps = 23.439291 - 0.0130042 * TTT - 0.000000164 * TTT2 + 0.000000504 * TTT3;
  Eps = Mod(Eps, 360.0);
  Eps = Eps * Deg2Rad;

  rr    = 360.0; // deg
  l     =  134.96340251 + (1325 * rr  + 198.8675605) * TTT + 0.0088553  * TTT2 +
                          0.000014343 * TTT3 - 0.00000006797  *  TTT4;
  l1    =  357.52910918 + (99 * rr  + 359.0502911) * TTT - 0.0001537 * TTT2 -
                          0.000000038 * TTT3 - 0.00000000319 * TTT4;
  F     =   93.27209062 + (1342 * rr  +  82.0174577) * TTT - 0.0035420 * TTT2 +
                          0.000000288 * TTT3 + 0.00000000116 * TTT4;
  D     =  297.85019547 + (1236 * rr  + 307.1114469) * TTT - 0.0017696 * TTT2 +
                          0.000001831 * TTT3 - 0.00000000880 * TTT4;
  Omega =  125.04455501 - (5 * rr  + 134.1361851) * TTT + 0.0020756 * TTT2 +
                          0.000002139 * TTT3 - 0.00000001650 * TTT4;

  l     = Mod(l,    360.0) * Deg2Rad; // rad
  l1    = Mod(l1,   360.0) * Deg2Rad;
  F     = Mod(F,    360.0) * Deg2Rad;
  D     = Mod(D,    360.0) * Deg2Rad;
  Omega = Mod(Omega,360.0) * Deg2Rad;

  DPsi96 = 0.0;
  DEps96 = 0.0;

  for (UINT i = 263; i > 0; i--)
  {
    TempVal = IAr96[0][i-1] * l + IAr96[1][i-1] * l1 + IAr96[2][i-1] * F + 
              IAr96[3][i-1] * D + IAr96[4][i-1] * Omega;
    DPsi96 = DPsi96 + (RAr96[0][i-1] + RAr96[1][i-1] * TTT) * sin(TempVal);
    DEps96 = DEps96 + (RAr96[2][i-1] + RAr96[3][i-1] * TTT) * cos(TempVal);
  }

  /* --- Determine coefficients for 1996 Planetary NUTATION Theory ----- */
  LonVen   = 181.979800853  +    58517.8156748  * TTT;   // deg
  LonEar   = 100.466448494  +    35999.3728521  * TTT;
  LonMar   = 355.433274605  +    19140.299314   * TTT;
  LonJup   =  34.351483900  +     3034.90567464 * TTT;
  LonSat   =  50.0774713998 +     1222.11379404 * TTT;
  PrecRate =   1.39697137214 * TTT + 0.0003086  * TTT2;

  LonVen   = Mod(LonVen,   360.0) * Deg2Rad;  // rad
  LonEar   = Mod(LonEar,   360.0) * Deg2Rad;
  LonMar   = Mod(LonMar,   360.0) * Deg2Rad;
  LonJup   = Mod(LonJup,   360.0) * Deg2Rad;
  LonSat   = Mod(LonSat,   360.0) * Deg2Rad;
  PrecRate = Mod(PrecRate, 360.0) * Deg2Rad;

  DPsi96P  = 0.0;
  DEps96P  = 0.0;
  for (UINT i = 112; i > 0; i--)
  {
    TempVal = pIAr96[0][i-1] * LonVen + pIAr96[1][i-1] * LonEar + 
              pIAr96[2][i-1] * LonMar + pIAr96[3][i-1] * LonJup + 
              pIAr96[4][i-1] * LonSat + pIAr96[5][i-1] * PrecRate +
              pIAr96[6][i-1] * D      + pIAr96[7][i-1] * F +
              pIAr96[8][i-1] * l      + pIAr96[9][i-1] * Omega;
    DPsi96P = DPsi96P + (pRAr96[0][i-1] + pRAr96[1][i-1] * TTT) * sin(TempVal);
    DEps96P = DEps96P + (pRAr96[2][i-1] + pRAr96[3][i-1] * TTT) * cos(TempVal);
  }

  /* --------------- Find NUTATION Parameters --------------------- */
  /* --- Add in Offset and Rate terms since J2000 is referenced --- */
  Dp        = DPsi96;
  De        = DEps96;
  DPsi96    = DPsi96 - 0.0431 / 3600.0 - 0.2957 / 3600.0 * TTT + DPsi96P;
  DEps96    = DEps96 - 0.0051 / 3600.0 - 0.0227 / 3600.0 * TTT + DEps96P;
  DPsi96    = Mod(DPsi96, 360.0) * Deg2Rad; // rad
  DEps96    = Mod(DEps96, 360.0) * Deg2Rad;
  TrueEps96 = Eps + DEps96;

  if (dir == TOO)
  {
    Temp  = rMOD.Rot1(Eps);
    Temp1 = Temp.Rot3(-DPsi96);
    rTOD  = Temp1.Rot1(-TrueEps96);
    Temp  = vMOD.Rot1(Eps);
    Temp1 = Temp.Rot3(-DPsi96);
    vTOD  = Temp1.Rot1(-TrueEps96);
  }
  else
  {
    Temp  = rTOD.Rot1(TrueEps96);
    Temp1 = Temp.Rot3(DPsi96);
    rMOD  = Temp1.Rot1(-Eps);
    Temp  = vTOD.Rot1(TrueEps96);
    Temp1 = Temp.Rot3(DPsi96);
    vMOD  = Temp1.Rot1(-Eps);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  96 GCRF NUTATION ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "l       %18.12fø%18.12f", l  / Deg2Rad, l);
      fprintf(FileOut, "l1      %18.12fø%18.12f", l1 / Deg2Rad, l1);
      fprintf(FileOut, "f       %18.12fø%18.12f", F  / Deg2Rad, F);
      fprintf(FileOut, "d       %18.12fø%18.12f", D  / Deg2Rad, D);
      DMS_Rad(Deg, Min, Sec, FROM, Eps);
      fprintf(FileOut, "Eps     %18.12fø%3d%3d%9.5f%18.12f",
                        Eps / Deg2Rad, Deg, Min, Sec, Eps);
      fprintf(FileOut, "DPsi263 %18.12fø%18.14f", Dp, Dp * Deg2Rad);
      fprintf(FileOut, "DEps263 %18.12fø%18.14f", De, De * Deg2Rad);
      fprintf(FileOut, "DPsi96p %18.12fø", DPsi96P);
      fprintf(FileOut, "DEps96p %18.12fø", DEps96P);
      fprintf(FileOut, "DEps96 %18.12fø%18.12f", DEps96 / Deg2Rad, DEps96);
      fprintf(FileOut, "Output vars --------\n");
      fprintf(FileOut, "DPsi96 %18.12fø%18.12f", DPsi96 / Deg2Rad, DPsi96);
      DMS_Rad(Deg, Min, Sec, FROM, TrueEps96);
      fprintf(FileOut, "TrueEps96%18.12fø%3d%3d%9.5f",
                        TrueEps96 / Deg2Rad, Deg, Min, Sec);
      fprintf(FileOut, "omega   %18.12fø%18.12f", Omega / Deg2Rad, Omega);
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE POLARM
|
|  This PROCEDURE converts position and velocity vectors  between the
|    True equator True equinox of date (ECEF) and the True equator True equinox
|    of date (ITRF).  The results take into account the effects of Polar
|    Motion.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rECEF       - Position vector of date
|                    True Equator, True Equinox   ER
|    vECEF       - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|    TTT         - Julian Centuries of TDB        centuries
|    Direction   - Which set of vars to output    FROM  TOO
|    xp          - Polar motion coefficient       rad
|    yp          - Polar motion coefficient       rad
|
|  Outputs       :
|    rITRF       - Position vector of date
|                    True Equator, True Equinox   ER
|    vITRF       - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|
|  Locals        :
|    None.
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|
|  References    :
|    Vallado       2001, 219-220, Eq 3-68
|
  ----------------------------------------------------------------------------*/
void PolarM
    (
      Vector& rECEF, Vector& vECEF, Direction dir, 
      Vector& rITRF, Vector& vITRF, Real xp, Real yp
    )
{
  if (dir == TOO)
  {
    rITRF.Set(rECEF.Get(1) + xp * rECEF.Get(3), 1);
    rITRF.Set(rECEF.Get(2) - yp * rECEF.Get(3), 2);
    rITRF.Set(rECEF.Get(3) - xp * rECEF.Get(1) + yp * rECEF.Get(2), 3);
    vITRF.Set(vECEF.Get(1) + xp * vECEF.Get(3), 1);
    vITRF.Set(vECEF.Get(2) - yp * vECEF.Get(3), 2);
    vITRF.Set(vECEF.Get(3) - xp * vECEF.Get(1) + yp * vECEF.Get(2), 3);
  }
  else
  {
    rECEF.Set(rITRF.Get(1) - xp * rITRF.Get(3), 1);
    rECEF.Set(rITRF.Get(2) + yp * rITRF.Get(3), 2);
    rECEF.Set(rITRF.Get(3) + xp * rITRF.Get(1) - yp * rITRF.Get(2), 3);
    vECEF.Set(vITRF.Get(1) - xp * vITRF.Get(3), 1);
    vECEF.Set(vITRF.Get(2) + yp * vITRF.Get(3), 2);
    vECEF.Set(vITRF.Get(3) + xp * vITRF.Get(1) - yp * vITRF.Get(2), 3);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  POLAR MOTION ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "xp %12.9f yp %12.9f rad\n", xp, yp);
      fprintf(FileOut, "xp %12.9f yp %12.9f rad\n",
                        xp * 3600.0 * 180.0 * PI, yp * 3600.0 * 180.0 * PI);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "Output vars --------\n");
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE PRECESSION
|
|  This PROCEDURE converts position and velocity vectors between the
|    J2000 epoch (FK5 ECI) frame (Mean Equator Mean Equinox of J2000) and the
|    Mean equator Mean equinox of date (MOD).
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rJ2000      - Initial J2000 ECI Position     ER
|    vJ2000      - Initial J2000 ECI Velocity     ER/TU
|    TTT         - Julian Centuries of TT         centuries
|    Direction   - Which set of vars to output    FROM  TOO
|
|  Outputs       :
|    rMOD        - Position vector of date
|                    Mean Equator, Mean Equinox   ER
|    vMOD        - Velocity vector of date
|                    Mean Equator, Mean Equinox   ER/TU
|
|  Locals        :
|    TTT2        - TTT squared
|    TTT3        - TTT cubed
|    Zeta        - PRECESSION ANGLE               rad
|    z           - PRECESSION ANGLE               rad
|    Theta       - PRECESSION ANGLE               rad
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|
|  References    :
|    Vallado       2001, 214-215, Eq 3-57
|
  ----------------------------------------------------------------------------*/
void Precession
    (
      Vector& rJ2000, Vector& vJ2000, Direction dir, 
      Vector& rMOD, Vector& vMOD, Real TTT
    )
{
  const Real Deg2Rad = PI / 180.0;

  Vector Temp(3), Temp1(3);
  Real   Zeta, Z, Theta, TTT2, TTT3;

  /* ---------------------- 1976 IAU PRECESSION angles --------------------- */
  TTT2 = TTT  * TTT;
  TTT3 = TTT2 * TTT;
  Zeta  = 2306.2181 * TTT + 0.30188 * TTT2 + 0.017998 * TTT3; // "
  Theta = 2004.3109 * TTT - 0.42665 * TTT2 - 0.041833 * TTT3;
  Z     = 2306.2181 * TTT + 1.09468 * TTT2 + 0.018203 * TTT3;

  Zeta  = Zeta  * Deg2Rad / 3600.0;  // rad
  Theta = Theta * Deg2Rad / 3600.0;
  Z     = Z     * Deg2Rad / 3600.0;

  /* ---------------------- Perform rotations --------------------- */
  if (dir == TOO)
  {
    Temp1 = rJ2000.Rot3(-Zeta);
    Temp  = Temp1.Rot2(Theta);
    rMOD  = Temp.Rot3(-Z);
    Temp1 = vJ2000.Rot3(-Zeta);
    Temp  = Temp1.Rot2(Theta);
    vMOD  = Temp.Rot3(-Z);
  }
  else
  {
    Temp   = rMOD.Rot3(Z);
    Temp1  = Temp.Rot2(-Theta);
    rJ2000 = Temp1.Rot3(Zeta);
    Temp   = vMOD.Rot3(Z);
    Temp1  = Temp.Rot2(-Theta);
    vJ2000 = Temp1.Rot3(Zeta);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------ PRECESSION Rotation Data ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "TTT     %18.12f\n", TTT);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "Output vars --------\n");
      fprintf(FileOut, "Zeta    %18.12fø%18.12f\n", Zeta / Deg2Rad, Zeta);
      fprintf(FileOut, "Z       %18.12fø%18.12f\n", Z / Deg2Rad, Z);
      fprintf(FileOut, "Theta   %18.12fø%18.12f\n", Theta / Deg2Rad, Theta);
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE SIDEREAL
|
|  This PROCEDURE converts position and velocity vectors between the
|    True equator True equinox of date (ECI) and the True equator True equinox
|    of date (ECEF).  The results take into account the effects of SIDEREAL
|    Time.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    rTOD        - Position vector of date
|                    True Equator, True Equinox   ER
|    vTOD        - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|    Direction   - Which set of vars to output    FROM  TOO
|    JDUT1       - Julian Date of UT1             days from 4713 BC
|    Omega       -                                rad
|    DPsi        - NUTATION ANGLE                 rad
|    TrueEps     - True obliquity of the ecliptic rad
|
|  Outputs       :
|    rECEF       - Position vector of date
|                    True Equator, True Equinox   ER
|    vECEF       - Velocity vector of date
|                    True Equator, True Equinox   ER/TU
|
|  Locals        :
|    GMST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
|    AST         - Apparent ST                    0 to 2Pi rad
|    Hr          - hour                           hr
|    MIN         - Minutes                        MIN
|    SEC         - seconds                        SEC
|    Temp        - Temporary vector
|    TempVal     - Temporary variable
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|    MAG         - Magnitude of a vector
|
|  References    :
|    Vallado       2001, 219, Eq 3-65 to 3-66
|
  ----------------------------------------------------------------------------*/
void Sidereal
    (
      Vector& rTOD, Vector& vTOD, Direction dir, Vector& rECEF, Vector& vECEF,
      Real JDUT1, Real Omega, Real DPsi, Real TrueEps, Real XLOD
    )
{
  const Real TwoPi      = 2.0 * PI;
  const Real Deg2Rad    = PI / 180.0;

  Vector Temp(3);
  SINT   Hr, Min;
  Real   TUT1, TUMin, ThetaSA, EqE, GMST, AST, TempVal, Sec, OmegaEarth;

  /* ------------------------ Find Mean GMST ----------------------- */
  GMST = GSTime(JDUT1);

  /* ------------------------ Find Mean AST ----------------------- */
  if (JDUT1 < 2450450.5)  // 1 Jan 1997
    EqE = DPsi * cos(TrueEps); // rad
  else
    EqE = DPsi * cos(TrueEps) + 0.00264 * Deg2Rad / 3600.0 * sin(Omega) +
          0.000063 * Deg2Rad / 3600.0 * sin(2.0 * Omega);
  AST = GMST + EqE;

  TUT1    = (JDUT1 - 2451545.0) / 36525.0;
  ThetaSA = 7.29211514670698E-05 * (1.0 - XLOD / 86400.0); // rad

  TUMin =       13.446849855110;
  TUMin = sqrt(Power(6378.1363, 3) / 398600.4418) / 60.0;

  OmegaEarth = ThetaSA * TUMin * 60.0;

  if (dir == TOO)
  {
    rECEF = rTOD.Rot3(AST);
    Temp  = vTOD.Rot3(AST);
    vECEF.Set(Temp.Get(1) + rECEF.Get(2) * OmegaEarth, 1);  // RatioSidSol
    vECEF.Set(Temp.Get(2) + rECEF.Get(1) * OmegaEarth, 2);  // RatioSidSol
    vECEF.Set(Temp.Get(3), 3);
  }
  else
  {
    Temp.Set(vECEF.Get(1) - rECEF.Get(2) * OmegaEarth, 1);
    Temp.Set(vECEF.Get(2) - rECEF.Get(1) * OmegaEarth, 2);
    Temp.Set(vECEF.Get(3), 3);
    vTOD = Temp.Rot3(-AST);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  SIDEREAL TIME ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "JDUT1   %18.12f\n", JDUT1);
      fprintf(FileOut, "omega   %18.12fø%18.12f\n", Omega / Deg2Rad, Omega);
      fprintf(FileOut, "DPsi    %18.12fø%18.12f\n", DPsi  / Deg2Rad, DPsi);
      fprintf(FileOut, "TrueEps %18.12fø\n", TrueEps  / Deg2Rad);
      fprintf(FileOut, "XLOD    %14.9\n", XLOD);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "ThetaSA %20.16f rad/tu %20.16f\n", 
                        ThetaSA, ThetaSA / (TUMin / 1440.0)/86400.0);
      fprintf(FileOut, "GMST    %18.16f%18.16f\n", GMST * 180.0 / PI, GMST);
      fprintf(FileOut, "AST     %18.16f%18.16f\n",  AST * 180.0 / PI, GMST);
      fprintf(FileOut, "OmegaEarth%18.16fxx\n", OmegaEarth);
      HMS_Rad(Hr, Min, Sec, FROM, GMST);
      fprintf(FileOut, "GMST (hms) %3d%3d%9.5f\n", Hr, Min, Sec);
      HMS_Rad(Hr, Min, Sec, FROM, AST);
      fprintf(FileOut, "AST (hms) %3d%3d%9.5f\n", Hr, Min, Sec);
      TempVal = AST - GMST;
      HMS_Rad(Hr, Min, Sec, FROM, TempVal);
      fprintf(FileOut, " EQEQ %14.8fsec\n", Sec);
      fprintf(FileOut, "Eqe    %20.14f %f\n", EqE * 180.0 / PI, EqE);
      fprintf(FileOut, "Output vars --------\n");
    }
}

/*----------------------------------------------------------------------------
|
|                           PROCEDURE TRUEMEAN
|
|  This PROCEDURE converts position and velocity vectors between the
|    NORAD True Equator Mean Equinox of date and the Mean equator Mean equinox
|    of date (ECI).  The results approximate the effects of NUTATION.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units        Ex Value
|    rMOD        - Position vector of date
|                    Mean Equator, Mean Equinox   ER
|    vMOD        - Velocity vector of date
|                    Mean Equator, Mean Equinox   ER/TU
|    TTT         - Julian Centuries of TDB
|    Direction   - Which set of vars to output    FROM  TOO            TOO
|
|  Outputs       :
|    rTm         - Position vector of date
|                    Mean Equator, True Equinox   ER
|    vTm         - Velocity vector of date
|                    Mean Equator, True Equinox   ER/TU
|    DPsiSp  - NUTATION ANGLE                 rad
|    TrueEpsSp   - True obliquity of the ecliptic rad
|
|  Locals        :
|    TTT2        - TTT squared
|    TTT3        - TTT cubed
|    Eps         - Mean obliquity of the ecliptic rad
|    l           -                                rad
|    ll          -                                rad
|    F           -                                rad
|    D           -                                rad
|    Omega       -                                rad
|    DEps    - Change in obliquity            rad
|
|  Coupling      :
|    ROT2        - Rotation about the second axis
|    ROT3        - Rotation about the third axis
|    REALMOD     -
|
|  References    :
|    Vallado       2001, 227
|
  ----------------------------------------------------------------------------*/
void TrueMean
    (
      Vector& rMod, Vector& vMod, Direction dir, Vector& rTM, Vector& vTM, 
      Real& DPsiTM, Real& TrueEpsTM, Real TTT, IAr5x106 IAr80, RAr4x106 RAr80
    )
{
  const Real Deg2Rad = PI / 180.0;

  SINT i;
  Real DEpsTM,TempVal, TTT2, TTT3, Eps, rr, l, l1, F, D, Omega;

  /* ---- Determine coefficients for IAU 1980 NUTATION Theory ----- */
  TTT2 = TTT  * TTT;
  TTT3 = TTT2 * TTT;
  Eps  = 23.439291 - 0.0130042 * TTT - 0.000000164 * TTT2 + 0.000000504 * TTT3;
  Eps  = Mod(Eps, 360.0);
  Eps  = Eps * Deg2Rad;

  rr    = 360.0;  // deg
  l     =  134.9629814 + (1325 * rr + 198.8673981) * TTT + 0.0086972 * TTT2 + 
             0.00001778 * TTT3;
  l1    = 357.5277233 + (99 * rr + 359.05034) * TTT - 0.00016028 * TTT2 - 
            0.00000333 * TTT3;
  F     =   93.2719103 + (1342 * rr +  82.0175381) * TTT - 0.0036825 * TTT2 + 
            0.00000306 * TTT3;
  D     =  297.8503631 + (1236 * rr + 307.111480) * TTT - 0.00191417 * TTT2 + 
            0.00000528 * TTT3;
  Omega =  125.0445222 - (5 * rr + 134.1362608) * TTT + 0.0020708 * TTT2 + 
             0.00000222 * TTT3;
  l     = Mod(l,     360.0)  * Deg2Rad;
  l1    = Mod(l1,    360.0)  * Deg2Rad;
  F     = Mod(F,     360.0 ) * Deg2Rad;
  D     = Mod(D,     360.0 ) * Deg2Rad;
  Omega = Mod(Omega, 360.0 ) * Deg2Rad;

  DPsiTM = 0.0;
  DEpsTM = 0.0;

  for (UINT ii = 1; i <=4; i++)
  {
    switch (ii)
    {
      case 1:
      case 2:
        i = ii;
        break;
      case 3:
        i = 9;
        break;
      case 4:
        i = 31;
        break;
    }

    i = ii;

    TempVal = IAr80[0][i-1] * l + IAr80[1][i-1] * l1 + IAr80[2][i-1] * F + 
              IAr80[3][i-1] * D + IAr80[4][i-1] * Omega;
    DPsiTM  = DPsiTM + (RAr80[0][i-1] + RAr80[1][i-1] * TTT) * sin(TempVal);
    DEpsTM  = DEpsTM + (RAr80[2][i-1] + RAr80[3][i-1] * TTT) * cos(TempVal);
  }

  /* --------------- Find Approx Nutation Parameters ---------------- */
  DPsiTM    = Mod( DPsiTM, 360.0 ) * Deg2Rad;
  DEpsTM    = Mod( DEpsTM, 360.0 ) * Deg2Rad;
  TrueEpsTM = Eps + DEpsTM;

  if (dir == TOO)
  {
    rTM.Set(rMod.Get(1) - DPsiTM * sin(Eps) * rMod.Get(3), 1);
    rTM.Set(rMod.Get(2) - DEpsTM * rMod.Get(3), 2);
    rTM.Set(rMod.Get(3) + DPsiTM * sin(Eps) * rMod.Get(1) + 
                          DEpsTM * rMod.Get(2), 3);
    vTM.Set(vMod.Get(1) - DPsiTM * sin(Eps) * vMod.Get(3), 1);
    vTM.Set(vMod.Get(2) - DEpsTM * vMod.Get(3), 2);
    vTM.Set(vMod.Get(3) + DPsiTM * sin(Eps) * vMod.Get(1) + 
                          DEpsTM * vMod.Get(2), 3);
  }
  else
  {
    rMod.Set(rTM.Get(1) + DPsiTM * sin(Eps) * rTM.Get(3), 1);
    rMod.Set(rTM.Get(2) + DEpsTM * rTM.Get(3), 2);
    rMod.Set(rTM.Get(3) - DPsiTM * sin(Eps) * rTM.Get(1) - 
                          DEpsTM * rTM.Get(2), 3);
    vMod.Set(vTM.Get(1) + DPsiTM * sin(Eps) * vTM.Get(3), 1);
    vMod.Set(vTM.Get(2) + DEpsTM * vTM.Get(3), 2);
    vMod.Set(vTM.Get(3) - DPsiTM * sin(Eps) * vTM.Get(1) - 
                          DEpsTM * vTM.Get(2), 3);
  }

  if (Show == 'I')
    if (FileOut != NULL)
    {
      fprintf(FileOut, "------  TRUE MEAN ----------\n");
      fprintf(FileOut, "Input vars ---------\n");
      fprintf(FileOut, "TTT     %18.12f\n", TTT);
      fprintf(FileOut, "Intermediate vars --\n");
      fprintf(FileOut, "Output vars --------\n");
      fprintf(FileOut, "DelPsi  %18.12fø%18.12f\"\n",
                        DPsiTM / Deg2Rad, DPsiTM * 3600.0 / Deg2Rad);
      fprintf(FileOut, "TrueEps %18.12fø\n", TrueEpsTM / Deg2Rad);
      fprintf(FileOut, "DelEps  %18.12fø%18.12f\"\n",
                        DEpsTM / Deg2Rad, DEpsTM * 3600.0 / Deg2Rad);
    }
}