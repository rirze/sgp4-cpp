/*----------------------------------------------------------------      


                               UNIT ASTMANV


    This file contains fundamental Astrodynamic procedures and functions     
    relating to orbit transfer calculations. Ch 6 describes each of these    
    routines.                                                                

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

#include "asttime.h"  // for show and fileout
#include "astmanv.h"

/*------------------------------------------------------------------------------
|
|                           PROCEDURE ONETANGENT
|
|  This PROCEDURE calculates the delta v's for a One Tangent transfer for either
|    circle to circle, or ellipse to ellipse.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInit       - Initial position magnitude     ER
|    RFinal      - Final position magnitude       ER
|    eInit       - Eccentricity of first orbit
|    eFinal      - Eccentricity of final orbit
|    NuInit      - True Anomaly of first orbit    0 or Pi rad
|    Nu2         - True Anomaly of second orbit   Same quad as NuInit, rad
|    NuFinal     - True Anomaly of final orbit    0 or Pi rad
|
|  OutPuts       :
|    DeltaVa     - Change in velocity at point a  ER / TU
|    DeltaVb     - Change in velocity at point b  ER / TU
|    DtTU        - Time of Flight for the transf  TU
|
|  Locals        :
|    SME1        - Mech Energy of first orbit     ER2 / TU
|    SME2        - Mech Energy of transfer orbit  ER2 / TU
|    SME3        - Mech Energy of final orbit     ER2 / TU
|    VInit       - Velocity of first orbit at a   ER / TU
|    vTransa     - Velocity of trans orbit at a   ER / TU
|    vTransb     - Velocity of trans orbit at b   ER / TU
|    VFinal      - Velocity of final orbit at b   ER / TU
|    aInit       - Semimajor axis of first orbit  ER
|    aTrans      - Semimajor axis of Trans orbit  ER
|    aFinal      - Semimajor axis of final orbit  ER
|    E           - Ecc anomaly of trans at b      rad
|    Ratio       - Ratio of initial to final
|                    orbit radii
|
|  Coupling      :
|    ATAN2       - Arc tangent rountine that solves quadrant ambiguities
|
|  References    :
|    Vallado       2001, 317-321, Alg 38, Ex 6-3
|
 -----------------------------------------------------------------------------*/
void OneTangent
    (
      Real rInit, Real rFinal, Real eInit, Real eFinal, Real NuInit, 
      Real NuTran, Real& DeltaVa, Real& DeltaVb, Real& DtTU
    )
{
  Real EAInit, VInit, VTrana, VTranb, VFinal, eTran, aInit, aTran, aFinal, 
       fpaTranb, fpaFinal, E, Sinv, Cosv, Ratio;

  /* --------------------  Initialize values   ------------------- */
  DeltaVa = 0.0;
  DeltaVb = 0.0;
  DtTU    = 0.0;
  Ratio   = rInit / rFinal;
  if (fabs(NuInit) < 0.01)   // check 0 or 180
  {
    eTran  = (Ratio - 1.0) / (cos(NuTran) - Ratio);  // init at perigee
    EAInit = 0.0;
  }
  else
  {
    eTran  = (Ratio - 1.0) / (cos(NuTran) + Ratio);  // init at apogee
    EAInit = PI;
  }
  if (eTran >= 0.0)
  {
    aInit  = (rInit * (1.0 + eInit * cos(NuInit))) / (1.0 - eInit * eInit);
    aFinal = (rFinal * (1.0 + eFinal * cos(NuTran))) / (1.0 - eFinal * eFinal);
              // nutran is used since it = nufinal!!
aInit  = rInit;
aFinal = rFinal;
    if (fabs(eTran - 1.0) > 0.000001)
      if (fabs(NuInit) < 0.01)     // check 0 or 180
        aTran = (rInit * (1.0 + eTran * cos(NuInit))) / 
                (1.0 - eTran * eTran); // per
      else
        aTran = (rInit * (1.0 + eTran * cos(NuInit))) /
                (1.0 + eTran * eTran); // apo
    else
      aTran = 999999.9;  // Infinite for Parabolic orbit

    /* -----------------  Find Delta v at point a  ----------------- */
    VInit   = sqrt(2.0 / rInit - (1.0 / aInit));
    VTrana  = sqrt(2.0 / rInit - (1.0 / aTran));
    DeltaVa = fabs(VTrana - VInit);

    /* -----------------  Find Delta v at point b  ----------------- */
    VFinal   = sqrt(2.0 / rFinal - (1.0 / aFinal));
    VTranb   = sqrt(2.0 / rFinal - (1.0 / aTran));
    fpaTranb = atan((eTran * sin(NuTran)) / (1.0 + eTran * cos(NuTran)));

    /* -----------------  Find Delta v at point b  ----------------- */
    VFinal   = sqrt(2.0 / rFinal - (1.0 / aFinal));
    VTranb   = sqrt(2.0 / rFinal - (1.0 / aTran));
    fpaTranb = atan((eTran * sin(NuTran)) / (1.0 + eTran * cos(NuTran)));

    DeltaVb = sqrt(VTranb * VTranb + VFinal * VFinal - 
                   2.0 * VTranb * VFinal * cos(fpaTranb - fpaFinal));

    /* ----------------  Find Transfer Time of Flight  ------------- */
    if (eTran < 0.99999)
    {
      Sinv = (sqrt(1.0 - eTran * eTran) * sin(NuTran)) / 
             (1.0 + eTran * cos(NuTran));
      Cosv = (eTran + cos(NuTran)) / (1.0 + eTran * cos(NuTran));
    }
    else
      if (fabs(eTran - 1.0) < 0.000001)
      {
        // Parabolic DtTU
      }
      else
      {
        // Hyperbolic DtTU
      }

    if ((Show == 'Y') || (Show == 'S'))
      if (FileOut != NULL)
      {
        fprintf(FileOut, " R-1 %10.6f\n", Ratio);
        fprintf(FileOut, "aInit  %10.6f ER\n", aInit);
        fprintf(FileOut, "aTran  %10.6f ER\n", aTran);
        fprintf(FileOut, "aFinal %10.6f ER\n", aFinal);
        fprintf(FileOut, "VInit  %10.6f  VTrana %10.6f  DVa %10.6f\n", 
                         VInit, VTrana, DeltaVa);
        fprintf(FileOut, " fpaTtran %10.6f  fpaFinal = %10.6f\n",
                         fpaTranb * 57.2955, fpaFinal * 57.2955);
        fprintf(FileOut, "Vtranb %10.6f  VFinal  %10.6f  DVb %10.6f\n",
                         VTranb, VFinal, DeltaVb);
        fprintf(FileOut, "Total DV  %10.6f %12.7f ER/TU\n",
                         (DeltaVa+DeltaVb), (DeltaVa+DeltaVb) * 7.90536599);
        fprintf(FileOut, "Total DtTU %11.6f %13.6f  min %11.6f hrs\n",
                     DtTU, DtTU * 13.44685115881, DtTU * 13.44685115881 / 60);
        fprintf(FileOut, " eTran   %10.6f  E = %10.6f  aTran %10.6f\n",
                         eTran, E * 57.2955, aTran);
      }
  }
  else
    printf("The one tangent burn is not possible for this case \n");
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE BIELLIPTIC
|
|  This PROCEDURE calculates the delta v's for a Bi-elliptic transfer for either
|    circle to circle, or ellipse to ellipse.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInit       - Initial position magnitude     ER
|    R2          - Interim orbit magnitude        ER
|    RFinal      - Final position magnitude       ER
|    eInit       - Eccentricity of first orbit
|    eFinal      - Eccentricity of final orbit
|    NuInit      - True Anomaly of first orbit    0 or Pi rad
|    NuFinal     - True Anomaly of final orbit    0 or Pi rad, Opp of NuInit
|
|  OutPuts       :
|    DeltaVa     - Change in velocity at point a  ER / TU
|    DeltaVb     - Change in velocity at point b  ER / TU
|    DtTU        - Time of Flight for the trans   TU
|
|  Locals        :
|    SME1        - Mech Energy of first orbit     ER2 / TU
|    SME2        - Mech Energy of transfer orbit  ER2 / TU
|    SME3        - Mech Energy of final orbit     ER2 / TU
|    VInit       - Velocity of first orbit at a   ER / TU
|    vTransa     - Velocity of trans orbit at a   ER / TU
|    vTransb     - Velocity of trans orbit at b   ER / TU
|    VFinal      - Velocity of final orbit at b   ER / TU
|    aInit       - Semimajor axis of first orbit  ER
|    aTrans      - Semimajor axis of Trans orbit  ER
|    aFinal      - Semimajor axis of final orbit  ER
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 311-314, Alg 37, Ex 6-2
|
 -----------------------------------------------------------------------------*/
void BiElliptic
    (
      Real rInit, Real Rb, Real rFinal, Real eInit, Real eFinal, Real NuInit,
      Real NuFinal, Real& DeltaVa, Real& DeltaVb, Real& DeltaVc, Real& DtTU
    )
{
  Real VInit, VTran1a, VTran1b, VTran2b, VTran2c, VFinal, 
       aInit, aTran1, aTran2, aFinal;

  /* --------------------  Initialize values   ------------------- */
  aInit  = (rInit * (1.0 + eInit * cos(NuInit))) / (1.0 - eInit * eInit);
  aTran1 = (rInit + Rb) * 0.5;
  aTran2 = (Rb + rFinal) * 0.5;
  aFinal = (rFinal * (1.0+eFinal * cos(NuFinal))) / (1.0 - eFinal * eFinal);

  DeltaVa = 0.0;
  DeltaVb = 0.0;
  DeltaVc = 0.0;
  DtTU    = 0.0;

  if ((eInit < 1.0) && (eFinal < 1.0)) 
  {
    /* -----------------  Find Delta v at point a  ----------------- */
    VInit   = sqrt(2.0 / rInit - (1.0 / aInit));
    VTran1a = sqrt(2.0 / rInit - (1.0 / aTran1));
    DeltaVa = fabs(VTran1a - VInit);

    /* -----------------  Find Delta v at point b  ----------------- */
    VTran1b = sqrt(2.0 / Rb - (1.0 / aTran1));
    VTran2b = sqrt(2.0 / Rb - (1.0 / aTran2));
    DeltaVb = fabs(VTran1b - VTran2b );

    /* -----------------  Find Delta v at point c  ----------------- */
    VTran2c = sqrt(2.0 / rFinal - (1.0 / aTran2));
    VFinal  = sqrt(2.0 / rFinal - (1.0 / aFinal));
    DeltaVc = fabs(VFinal - VTran2c );

    /* ----------------  Find Transfer Time of Flight  ------------- */
    DtTU = PI * sqrt(aTran1 * aTran1 * aTran1) +
           PI * sqrt(aTran2 * aTran2 * aTran2);

    if ((Show == 'Y') || (Show == 'S'))
      if (FileOut != NULL)
      {
        fprintf(FileOut, "aInit  %10.6f ER\n", aInit);
        fprintf(FileOut, "aTran1 %10.6f ER\n", aTran1);
        fprintf(FileOut, "aTran2 %10.6f ER\n", aTran2);
        fprintf(FileOut, "aFinal %10.6f ER\n", aFinal);
        fprintf(FileOut, "VInit   %10.6f  VTran1a %10.6f  DVa %10.6f\n",
                         VInit, VTran1a, DeltaVa);
        fprintf(FileOut, "VTran1b %10.6f  VTran2b %10.6f  DVb %10.6f\n",
                         VTran1b, VTran2b, DeltaVb);
        fprintf(FileOut, "VTran2c %10.6f  VFinal  %10.6f  DVc %10.6f\n",
                         VTran2c, VFinal, DeltaVc);
        fprintf(FileOut, "Total DV  %10.6f %12.7f ER/TU\n", 
                         (DeltaVa + DeltaVb + DeltaVc), 
                         (DeltaVa + DeltaVb + DeltaVc) * 7.90536599);
        fprintf(FileOut, "Total DtTU %11.6f %13.6f min %11.6f hrs\n",
                      DtTU, DtTU * 13.44685115881, DtTU * 13.44685115881/60.0);
      }
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE COMBINEDPLANECHG
|
|  This PROCEDURE calculates the delta v's and the change in inclination
|    necessary for the change in velocity when traveling between two
|    non-coplanar orbits.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInit       - Initial position magnitude     ER
|    RFinal      - Final position magnitude       ER
|    eInit       - Ecc of first orbit
|    e2          - Ecc of trans orbit
|    eFinal      - Ecc of final orbit
|    NuInit      - True Anomaly of first orbit    0 or Pi rad
|    NuFinal     - True Anomaly of final orbit    0 or Pi rad
|    iInit       - Incl of the first orbit        rad
|    iFinal      - Incl of the second orbit       rad
|
|  OutPuts       :
|    Deltai1     - Amount of incl chg req at a    rad
|    DeltaVa     - Change in velocity at point a  ER / TU
|    DeltaVb     - Change in velocity at point b  ER / TU
|    DtTU        - Time of Flight for the trans   TU
|    NumIter     - Number of iterations
|
|  Locals        :
|    SME1        - Mech Energy of first orbit     ER2 / TU
|    SME2        - Mech Energy of transfer orbit  ER2 / TU
|    SME3        - Mech Energy of final orbit     ER2 / TU
|    VInit       - Velocity of first orbit at a   ER / TU
|    vTransa     - Velocity of trans orbit at a   ER / TU
|    vTransb     - Velocity of trans orbit at b   ER / TU
|    VFinal      - Velocity of final orbit at b   ER / TU
|    aInit       - Semimajor axis of first orbit  ER
|    aTrans      - Semimajor axis of Trans orbit  ER
|    aFinal      - Semimajor axis of final orbit  ER
|    e2          - Eccentricity of second orbit
|
|  Coupling      :
|    POWER       - Raise a base to a power
|    ARCSIN      - Arc sine routine
|
|  References    :
|    Vallado       2001, 336-339, Alg 42, Table 6-3
|
 -----------------------------------------------------------------------------*/
void CombinedPlaneChg
    (
      Real RInit, Real RFinal, Real eInit, Real e2, Real eFinal, Real NuInit, 
      Real Nu2a, Real Nu2b, Real NuFinal, Real Deltai,
      Real& DeltaVa, Real& DeltaVb, Real& DtTU
    )
{
  Real SME1, SME2, SME3, VInit, vTransa, vTransb, VFinal, E, Eo, Sinv, Cosv,
       a1, a2, a3, fpa1, fpa2a, fpa2b, fpa3;

  /* --------------------  Initialize values   ------------------- */
  a1 = (RInit * (1.0 + eInit * cos(NuInit))) / (1.0 - eInit * eInit);
  if (fabs(e2 - 1.0) > 0.000001)
  {
    a2   = (RInit * (1.0+ e2 * cos(Nu2a))) / (1.0 - e2*e2);
    SME2 = -1.0 / (2.0 * a2);
  }
  else
  {
    a2   = 999999.9;  // Undefined for Parabolic orbit
    SME2 = 0.0;
  }
  a3   = (RFinal * (1.0 + eFinal * cos(NuFinal))) / (1.0 - eFinal * eFinal);
  SME1 = -1.0 / (2.0 * a1);
  SME3 = -1.0 / (2.0 * a3);

  /* -----------------  Find Delta v at point a  ----------------- */
  VInit   = sqrt(2.0 * ((1.0 / RInit) + SME1));
  vTransa = sqrt(2.0 * ((1.0 / RInit) + SME2));
  fpa2a   = atan((e2 * sin(Nu2a)) / (1.0 + e2 * cos(Nu2a)));
  fpa1    = atan((eInit * sin(NuInit)) / (1.0 + eInit * cos(NuInit)));
  DeltaVa = sqrt(vTransa * vTransa + VInit * VInit - 2.0 * vTransa * VInit*
              (sin(fpa2a) * sin(fpa1) + cos(fpa2a) * cos(fpa1) * cos(Deltai)));

  /* -----------------  Find Delta v at point b  ----------------- */
  VFinal  = sqrt(2.0 * ((1.0 / RFinal) + SME3));
  vTransb = sqrt(2.0 * ((1.0 / RFinal) + SME2));
  fpa2b   = atan((e2 * sin(Nu2b)) / ( 1.0 + e2 * cos(Nu2b)));
  fpa3    = atan((eFinal * sin(NuFinal)) / (1.0 + eFinal * cos(NuFinal)));
  DeltaVb = sqrt(vTransb * vTransb + VFinal * VFinal - 2.0 * vTransb * VFinal *
               (sin(fpa2b) * sin(fpa3) + cos(fpa2b) * cos(fpa3) * cos(Deltai)));

  /* ----------------  Find Transfer Time of Flight  ------------- */
  Sinv = (sqrt(1.0 - e2 * e2) * sin(Nu2b)) / (1.0 + e2 * cos(Nu2b));
  Cosv = (e2 + cos(Nu2b)) / (1.0 + e2 * cos(Nu2b));
  E    = atan2(Sinv, Cosv);
  Sinv = (sqrt(1.0-e2*e2) * sin(Nu2a)) / (1.0 + e2 * cos(Nu2a));
  Cosv = (e2 + cos(Nu2a)) / (1.0 + e2 * cos(Nu2a));
  Eo   = atan2(Sinv, Cosv);
  DtTU = sqrt(a2 * a2 * a2) * ((E - e2 * sin(E)) - (Eo - e2 * sin(Eo)));
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE COW2Hill
|
|  This PROCEDURE finds the equivalent relative motion vector given a geocentric
|    vector.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rtgt        - Position vector of tgt         ER
|    Vtgt        - Velocity Vector of tgt         ER / TU
|    Rint        - Position vector of int         ER
|    Vint        - Velocity Vector of int         ER / TU
|
|  Outputs       :
|    RHill       - Position vector of int rel to
|                  target                         ER
|    VHill       - Velocity Vector of int rel to
|                  target                         ER / TU
|
|  Locals        :
|    None.
|
|
|
|
|
|
|  Coupling      :
|    CROSS       - Cross product of two vectors
|    NORM        - Unit vector
|    IJK_RSW     -
|    LNCOM2      - Linear combination of two scalars and two vectors
|    MAG         - Magnitude of a vector
|    ROT3        - Rotation about the 3rd axis
|
|  References    :
|    Vallado       2001, 395-397
|
 -----------------------------------------------------------------------------*/
void Cow2Hill
    (
      Vector rtgt, Vector vtgt, Vector rint, Vector vint, 
      Vector& RHill, Vector& VHill
    )
{
  Vector R(3), S(3), W(3), RV(3), Rtgt(3), Vtgt(3),
         Rtemp(3), Vtemp(3), Rt(3), Vt(3), Ri(3), Vi(3);
  Real   Halfpi, angly, anglz, temp, LST, Lat, angl;

  /* ---- Form RSW unit vectors for transformation ---- */
  R  = rtgt.Norm();
  RV = rtgt.Cross(vtgt);
  W  = RV.Norm();
  S  = W.Cross(R);

  Rtgt = rtgt;
  Vtgt = vtgt;
  IJK_RSW(Rtgt, Vtgt, R, S, W, TOO, Rt, Vt);  // IJK to RSW
  IJK_RSW(rint, vint, R, S, W, TOO, Ri, Vi);  // IJK to RSW

  /* ---- Determine z offset to correct vector ---- */
  if (fabs(Ri.Get(3)) > 0.0000001)
  {
    anglz = atan(Ri.Get(3) / Rt.Mag());  // coord sys based at tgt
    Ri    = Ri.Rot2(-anglz);
    Vi    = Vi.Rot2(-anglz);  // should be ROT2(a), but opp
  }
  else
    anglz = 0.0;

  /* ---- Determine y offset to correct vector ---- */
  if (fabs(Ri.Get(2)) > 0.0000001)
  {
    angly = atan(Ri.Get(2) / Rt.Mag()); 
    Ri = Ri.Rot3(-angly);
    Vi = Vi.Rot3(-angly);  // should be ROT3(-a), but opp, but sign for later
  }
  else
    angly = 0.0;

  /* ---- Do all 3 here ---- */
  /* RHill = Rt - Ri; */
  VHill = Vi - Vt;

  /* ---- Now add in corrections ---- */
  RHill.Set(Ri.Get(1) - Rt.Mag(), 1);
  RHill.Set(angly * Ri.Mag(), 2);
  RHill.Set(anglz * Ri.Mag(), 3);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HILL2COW
|
|  This PROCEDURE finds the equivalent geocentric vector given the target and
|    relative motion vectors.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rtgt        - Position vector of tgt         ER
|    Vtgt        - Velocity Vector of tgt         ER / TU
|    RHill       - Position vector of int rel to
|                  target                         ER
|    VHill       - Velocity Vector of int rel to
|                  target                         ER / TU
|
|  Outputs       :
|    Rint        - Position vector of int         ER
|    Vint        - Velocity Vector of int         ER / TU
|
|  Locals        :
|    None.
|
|
|
|
|
|
|  Coupling      :
|    IJKRSW      - Form translation matrix given position and velocity vectors
|    ROT2        - Rotation about the 2nd axis
|    ROT3        - Rotation about the 3rd axis
|    MAG         - Magnitude of a vector
|    NORM        - Unit vector
|    CROSS       - Cross product of two vectors
|
|  References    :
|    Vallado       2001, 395-397
|
 -----------------------------------------------------------------------------*/
void Hill2Cow
    (
      Vector rTgtijk, Vector vTgtijk, Vector RHill, Vector VHill,
      Vector& rIntijk, Vector& vIntijk
    )
{
  Real   angly, anglz;
  Vector RTem(3), VTem(3), RTgtrsw(3), VTgtrsw(3), RV(3), R(3), S(3), W(3);
  Vector RTgtijk(3), VTgtijk(3);

  /* ---- Form RSW unit vectors for transformation ---- */
  R  = rTgtijk.Norm();
  RV = rTgtijk.Cross(vTgtijk);
  W  = RV.Norm();
  S  = W.Cross(R);

  RTgtijk = rTgtijk;
  VTgtijk = vTgtijk;
  IJK_RSW(RTgtijk, VTgtijk, R, S, W, TOO, RTgtrsw, VTgtrsw);  // IJK to RSW

  RTem.Set(RTgtrsw.Get(1) + RHill.Get(1), 1);  // in RSW
  RTem.Set(RTgtrsw.Get(2), 2);
  RTem.Set(RTgtrsw.Get(3), 3);
  VTem.Set(VTgtijk.Get(1) + VHill.Get(1), 1);
  VTem.Set(VTgtijk.Get(2) + VHill.Get(2), 2);
  VTem.Set(VTgtijk.Get(3) + VHill.Get(3), 3);

  /* ---- Now perform rotation to fix y ---- */
  if (fabs(RHill.Get(2)) > 0.0000001)
  {
    // rtgt, but IF x non-zero, needs extra
    angly = atan(RHill.Get(2) / RTem.Mag());
    RTem = RTem.Rot3(-angly);  // should be ROT3(a) but opp
    VTem = VTem.Rot3(-angly);
  }
  else
    angly = 0.0;

  /* ---- Now perform rotation to fix z ---- */
  if (fabs(RHill.Get(3)) > 0.0000001)
  {
    anglz = atan(RHill.Get(3) / RTem.Mag());
    RTem = RTem.Rot2(anglz);  // should be ROT4(-a) but opp
    VTem = VTem.Rot2(anglz);
  }
  else
    anglz = 0.0;

  /* ---- Now back to IJK via MatMult!! ---- */
  IJK_RSW(rIntijk, vIntijk, R, S, W, FROM, RTem, VTem);  // RSW to IJK
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HILLSR
|
|  This PROCEDURE calculates various position information for Hills equations.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - Initial Position vector of INT ER
|    V           - Initial Velocity Vector of INT ER / TU
|    Alt         - Altitude of TGT satellite      ER
|    DtTU        - Desired Time                   TU
|
|  Outputs       :
|    RInit       - Final Position vector of INT   ER
|    VInit       - Final Velocity Vector of INT   ER / TU
|
|  Locals        :
|    nt          - Angular velocity times time    rad
|    Omega       -
|    Sinnt       - Sine of nt
|    Cosnt       - Cosine of nt
|    Radius      - Magnitude of range vector      ER
|
|  Coupling      :
|    MAG         - Magnitude of a vector
|
|  References    :
|    Vallado       2001, 377-381, Alg 47, Ex 6-14
|
 -----------------------------------------------------------------------------*/
void HillsR
    (
      Vector R, Vector V, Real Alt, Real DtTU,
      Vector& RInit, Vector& VInit
    )
{
  Real SinNt, CosNt, Omega, nt, Radius;

  /* ---------------- Initialize the orbit elements -------------- */
  Radius = 1.0 + Alt;
  Omega  = sqrt(1.0 / (Radius * Radius * Radius));
  nt     = Omega * DtTU;
  CosNt  = cos(nt);
  SinNt  = sin(nt);

  /* ---------------- Determine new positions  ------------------- */
  RInit.Set((V.Get(1) / Omega) * SinNt -
             ((2.0 * V.Get(2) / Omega) + 3.0 * R.Get(1)) * CosNt +
             ((2.0 * V.Get(2) / Omega) + 4.0 * R.Get(1)), 1);
  RInit.Set((2.0 * V.Get(1) / Omega ) * CosNt +
             ((4.0 * V.Get(2) / Omega) + 6.0 * R.Get(1)) * SinNt +
             (R.Get(2) - (2.0 * V.Get(1) / Omega)) -
             (3.0 * V.Get(2) + 6.0 * Omega * R.Get(1)) * DtTU, 2);
  RInit.Set(R.Get(3) * CosNt + (V.Get(3) / Omega) * SinNt, 3);

  /* ---------------- Determine new velocities  ------------------ */
  VInit.Set(V.Get(1) * CosNt + 
            (2.0 * V.Get(2) + 3.0 * Omega * R.Get(1)) * SinNt, 1);
  VInit.Set(-2.0 * V.Get(1) * SinNt + 
            (4.0 * V.Get(2) + 6.0 * Omega * R.Get(1)) * CosNt -
            (3.0 * V.Get(2) + 6.0 * Omega * R.Get(1)), 2);
  VInit.Set(-R.Get(3) * Omega * SinNt + V.Get(3) * CosNt, 3);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HILLSV
|
|  This PROCEDURE calculates initial velocity for Hills equations.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    R           - Initial Position vector of INT ER
|    Alt         - Altitude of TGT satellite      ER
|    DtTU        - Desired Time                   TU
|
|  Outputs       :
|    V           - Initial Velocity Vector of INT ER / TU
|
|  Locals        :
|    Numer       -
|    Denom       -
|    nt          - Angular velocity times time    rad
|    Omega       -
|    Sinnt       - Sine of nt
|    Cosnt       - Cosine of nt
|    Radius      - Magnitude of range vector      ER
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 390-395, Eq 6-56, Ex 6-15
|
 -----------------------------------------------------------------------------*/
void HillsV(Vector R, Real Alt, Real DtTU, Vector& V)
{
  Real Numer, Denom, SinNt, CosNt, Omega, nt, Radius;

  /* --------------- Initialize the orbit elements --------------- */
  Radius = 1.0 + Alt;
  Omega  = sqrt(1.0 / (Radius * Radius * Radius));
  nt     = Omega * DtTU;
  CosNt  = cos(nt);
  SinNt  = sin(nt);

  /* ---------------- Determine initial Velocity ----------------- */
  Numer = ((6.0 * R.Get(1) * (nt - SinNt) - R.Get(2)) * Omega * SinNt - 
          2.0 * Omega * R.Get(1) * (4.0 - 3.0 * CosNt) * (1.0 - CosNt));
  Denom = (4.0 * SinNt - 3.0 * nt) * SinNt + 
          4.0 * (1.0 - CosNt) * (1.0 - CosNt);

  if (fabs(Denom) > 0.000001)
    V.Set(Numer / Denom, 2);
  else
    V.Set(0.0, 2);
  if (fabs(SinNt) > 0.000001)
    V.Set(-(Omega * R.Get(1) * (4.0 - 3.0 * CosNt) + 
          2.0 * (1.0 - CosNt) * V.Get(2)) / (SinNt), 1);
  else
    V.Set(0.0, 1);
  V.Set(-R.Get(3) * Omega * Cot(nt), 3);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE HOHMANN
|
|  This PROCEDURE calculates the delta v's for a Hohmann transfer for either
|    circle to circle, or ellipse to ellipse.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInit       - Initial position magnitude     ER
|    RFinal      - Final position magnitude       ER
|    eInit       - Eccentricity of first orbit
|    eFinal      - Eccentricity of final orbit
|    NuInit      - True Anomaly of first orbit    0 or Pi rad
|    NuFinal     - True Anomaly of final orbit    0 or Pi rad
|
|  OutPuts       :
|    DeltaVa     - Change in velocity at point a  ER / TU
|    DeltaVb     - Change in velocity at point b  ER / TU
|    DtTU        - Time of Flight for the trans   TU
|
|  Locals        :
|    SME1        - Mech Energy of first orbit     ER2 / TU
|    SME2        - Mech Energy of transfer orbit  ER2 / TU
|    SME3        - Mech Energy of final orbit     ER2 / TU
|    VInit       - Velocity of first orbit at a   ER / TU
|    vTransa     - Velocity of trans orbit at a   ER / TU
|    vTransb     - Velocity of trans orbit at b   ER / TU
|    VFinal      - Velocity of final orbit at b   ER / TU
|    aInit       - Semimajor axis of first orbit  ER
|    aTrans      - Semimajor axis of Trans orbit  ER
|    aFinal      - Semimajor axis of final orbit  ER
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 308-313, Alg 36, Ex 6-1
|
 -----------------------------------------------------------------------------*/
void Hohmann
    (
      Real rInit, Real rFinal, Real eInit, Real eFinal, Real NuInit,
      Real NuFinal, Real& DeltaVa, Real& DeltaVb, Real& DtTU
    )
{
  Real VInit, VTrana, VTranb, VFinal, aInit, aTran, aFinal;

  /* --------------------  Initialize values   ------------------- */
  aInit  = (rInit * (1.0 + eInit * cos(NuInit))) / (1.0 - eInit * eInit);
  aTran  = (rInit + rFinal) / 2.0;
  aFinal = (rFinal * (1.0 + eFinal * cos(NuFinal))) / (1.0 - eFinal * eFinal);
  DeltaVa = 0.0;
  DeltaVb = 0.0;
  DtTU    = 0.0;

  if ((eInit < 1.0) || (eFinal < 1.0))
  {
    /* -----------------  Find Delta v at point a  -------------- */
    VInit   = sqrt(2.0 / rInit - (1.0 / aInit));
    VTrana  = sqrt(2.0 / rInit - (1.0 / aTran));
    DeltaVa = fabs(VTrana - VInit);

    /* -----------------  Find Delta v at point b  -------------- */
    VFinal  = sqrt(2.0 / rFinal - (1.0 / aFinal));
    VTranb  = sqrt(2.0 / rFinal - (1.0 / aTran));
    DeltaVb = fabs(VFinal - VTranb );

    /* ----------------  Find Transfer Time of Flight  ---------- */
    DtTU = PI * sqrt(aTran * aTran * aTran);

    if ((Show == 'Y') || (Show == 'S'))
      if (FileOut != NULL)
      {
        fprintf(FileOut, "aInit  %10.6f ER\n", aInit);
        fprintf(FileOut, "aTran  %10.6f ER\n", aTran);
        fprintf(FileOut, "aFinal  %10.6f ER\n", aFinal);
        fprintf(FileOut, "Init  %10.6f VTrana %10.6 DVa %10.6f\n",
                         VInit, VTrana, DeltaVa);
        fprintf(FileOut, "Vtranb %10.6f VFinal  %10.6f DVb %10.6f\n",
                         VTranb, VFinal, DeltaVb);
        fprintf(FileOut, "Total DV  %10.6f %12.7f ER/TU\n",
                         (DeltaVa+DeltaVb), (DeltaVa+DeltaVb) * 7.90536599);
        fprintf(FileOut, "Total DtTU %11.6f %13.6f min %11.6f hrs\n",
                    DtTU, DtTU * 13.44685115881, DtTU * 13.44685115881 / 60.0);
      }
  }

}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE IandNodeChg
|
|  This PROCEDURE calculates the delta v's for a change in inclination and
|    longitude of ascending node.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    VInit       - Initial velocity vector        ER/TU
|    iInit       - Initial inclination            rad
|    fpa         - Flight path angle              rad
|    DeltaOmega  - Change in Node                 Rad
|    DeltaI      - Change in inclination          Rad
|    RFinal      - Final position magnitude       ER
|
|  OutPuts       :
|    iFinal      - Final inclination              rad
|    Deltav      - Change in velocity             ER/TU
|
|  Locals        :
|    ArgLat      - Argument of latitude           rad
|    ArgLat1     - Final Argument of latitude     rad
|    Theta       -
|
|  Coupling      :
|    ARCCOS      - Arc cosine FUNCTION
|
|  References    :
|    Vallado       2001, 334-335, Alg 41, Ex 6-6
|
 -----------------------------------------------------------------------------*/
void IandNodeChg
    (
      Real iInit, Real DeltaOmega, Real Deltai, Real VInit, Real fpa,
      Real& DeltaV, Real& iFinal
    )
{
  const Real Rad = 180.0 / PI;

  Real ArgLat, ArgLat1, Theta;

  iFinal = iInit - Deltai;
  Theta  = acos(cos(iInit) * cos(iFinal) +
                sin(iInit) * sin(iFinal) * cos(DeltaOmega));
  DeltaV = 2.0 * VInit * cos(fpa) * sin(0.5 * Theta);

  ArgLat = acos((sin(iFinal) * cos(DeltaOmega) -
                cos(Theta) * sin(iInit)) / (sin(Theta) * cos(iInit)));
  ArgLat1 = acos((cos(iInit) * sin(iFinal) -
                 sin(iInit) * cos(iFinal) * cos(DeltaOmega)) / sin(Theta));

  if (Show == 'Y')
    if (FileOut != NULL)
      fprintf(FileOut, "proce  dv %11.7f %11.7f km/s Theta %11.7f %11.6f ArgLat1 = %11.6f IF %11.6f\n",
                       DeltaV, DeltaV * 7.9053659986, Theta * Rad, 
                       ArgLat * Rad, ArgLat1 * Rad);
}

void IJK_RSW
    (
      Vector& Rijk, Vector& Vijk, Vector R, Vector S, Vector W,
      Direction dir, Vector& Rrsw, Vector& Vrsw
    )
{
  Matrix RotMat(3, 3);

  if (dir == FROM)
  {
    /* ------ Form Rotation Matrix ------ */
    for (UINT i = 1; i <= 3; i++)
    {
      RotMat.Set(R.Get(i), i, 1);
      RotMat.Set(S.Get(i), i, 2);
      RotMat.Set(W.Get(i), i, 3);
    }

    /* ------- Do multiplication -------- */
    Rijk = RotMat * Rrsw;
    Vijk = RotMat * Vrsw;
  }
  else
  {
    /* ------ Form Rotation Matrix ------ */
    for (UINT i = 1; i <= 3; i++)
    {
      RotMat.Set(R.Get(i), i, 1);
      RotMat.Set(S.Get(i), i, 2);
      RotMat.Set(W.Get(i), i, 3);
    }

    /* ------- Do multiplication -------- */
    Rrsw = RotMat * Rijk;
    Vrsw = RotMat * Vijk;
  }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE IONLYCHG
|
|  This PROCEDURE calculates the delta v's for a change in inclination only.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    DeltaI      - Change in inclination          rad
|    VInit       - Initial velocity vector        ER/TU
|    fpa         - Flight path angle              rad
|
|  OutPuts       :
|    DeltaVionly - answer
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
|  References    :
|    Vallado       2001, 328-331, Alg 39, Ex 6-4
|
 -----------------------------------------------------------------------------*/
void IOnlyChg(Real Deltai, Real VInit, Real fpa, Real& DeltaViOnly)
{
  DeltaViOnly = 2.0 * VInit * cos(fpa) * sin(0.5 * Deltai);
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE MINCOMBINEDPLANECHG
|
|  This PROCEDURE calculates the delta v's and the change in inclination
|    necessary for the minimum change in velocity when traveling between two
|    non-coplanar orbits.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    RInit       - Initial position magnitude     ER
|    RFinal      - Final position magnitude       ER
|    eInit       - Ecc of first orbit
|    e2          - Ecc of trans orbit
|    eFinal      - Ecc of final orbit
|    NuInit      - True Anomaly of first orbit    0 or Pi rad
|    NuFinal     - True Anomaly of final orbit    0 or Pi rad
|    iInit       - Incl of the first orbit        rad
|    iFinal      - Incl of the second orbit       rad
|
|  OutPuts       :
|    Deltai1     - Amount of incl chg req at a    rad
|    DeltaVa     - Change in velocity at point a  ER / TU
|    DeltaVb     - Change in velocity at point b  ER / TU
|    DtTU        - Time of Flight for the trans   TU
|    NumIter     - Number of iterations
|
|  Locals        :
|    SME1        - Mech Energy of first orbit     ER2 / TU
|    SME2        - Mech Energy of transfer orbit  ER2 / TU
|    SME3        - Mech Energy of final orbit     ER2 / TU
|    VInit       - Velocity of first orbit at a   ER / TU
|    vTransa     - Velocity of trans orbit at a   ER / TU
|    vTransb     - Velocity of trans orbit at b   ER / TU
|    VFinal      - Velocity of final orbit at b   ER / TU
|    aInit       - Semimajor axis of first orbit  ER
|    aTrans      - Semimajor axis of Trans orbit  ER
|    aFinal      - Semimajor axis of final orbit  ER
|    e2          - Eccentricity of second orbit
|
|  Coupling      :
|    POWER       - Raise a base to a power
|    ARCSIN      - Arc sine routine
|
|  References    :
|    Vallado       2001, 336-339, Alg 42, Table 6-3
|
 -----------------------------------------------------------------------------*/
void MinCombinedPlaneChg
    (
      Real RInit, Real RFinal, Real eInit, Real eFinal, Real NuInit,
      Real NuFinal, Real iInit, Real iFinal,
      Real& Deltai, Real& Deltai1, Real& DeltaVa, Real& DeltaVb, Real& DtTU
    )
{
  const Real Rad = 180.0 / PI;
  Integer    NumIter;
  Real       DeltaiNew, DVold, Temp, TDi, SME1, SME2, SME3, 
             VInit, V1t, V3t, VFinal, a1, a2, a3;

  /* --------------------  Initialize values   -------------------- */
  a1   = (RInit * (1.0 + eInit * cos(NuInit))) / (1.0 - eInit* eInit);
  a2   = 0.5 * (RInit + RFinal);
  a3   = (RFinal * (1.0 + eFinal * cos(NuFinal))) / (1.0 - eFinal * eFinal);
  SME1 = -1.0 / (2.0 * a1);
  SME2 = -1.0 / (2.0 * a2);
  SME3 = -1.0 / (2.0 * a3);

  /* ----------- Find velocities -------- */
  VInit  = sqrt(2.0 * ((1.0 / RInit) + SME1));
  V1t    = sqrt(2.0 * ((1.0 / RInit) + SME2));

  VFinal = sqrt(2.0 * ((1.0 / RFinal) + SME3));
  V3t    = sqrt(2.0 * ((1.0 / RFinal) + SME2));

  /* ----------- Find the optimum change of inclination ----------- */
  TDi = iFinal - iInit;

/*
  Temp = (1.0 / TDi) * atan((Power(RFinal / RInit, 1.5) - cos(TDi)) / sin(TDi));
*/

  Temp = (1.0 / TDi) * atan(sin(TDi) / (Power(RFinal / RInit, 1.5) + cos(TDi)));

  DeltaVa = sqrt(V1t * V1t + VInit * VInit - 
                 2.0 * V1t * VInit * cos(Temp * TDi));
  DeltaVb = sqrt(V3t * V3t + VFinal * VFinal - 
                 2.0 * V3t * VFinal * cos(TDi*(1.0-Temp)));

  Deltai  = Temp * TDi;
  Deltai1 = TDi * (1.0 - Temp);

  /* ----------------  Find Transfer Time of Flight  -------------- */
  DtTU = PI * sqrt(a2 * a2 * a2);

  if (Show == 'Y')
    if (FileOut != NULL)
    {
      DVold = fabs(V1t - VInit) + 
              sqrt(V3t * V3t + VFinal * VFinal - 2.0 * V3t * VFinal * cos(TDi));
      fprintf(FileOut, "s = %11.7f  this uses Di in rad\n", Temp);
      fprintf(FileOut, "RInit %14.7f %14.7f  RFinal %14.7f %14.7f\n",
                        RInit, RInit * 6378.137, RFinal, RFinal * 6378.137);
      fprintf(FileOut, "Deltai1 %13.7f %13.7f\n", Deltai * Rad, Deltai1 * Rad);
      fprintf(FileOut, "DeltaVa %13.7f DeltaVb %13.f ER/TU\n",
                        DeltaVa, DeltaVb);
      fprintf(FileOut, "DeltaVa %13.7f DeltaVb %13.7f km/s\n",
                        DeltaVa * 7.905365998, DeltaVb * 7.905365998);
      fprintf(FileOut, "%13.7f  m/s\n", 
                        1000 * (DeltaVa + DeltaVb) * 7.905365998);
      fprintf(FileOut, "Dv old way %13.7f m/s\n", 1000 * DVold * 7.905365998);
      fprintf(FileOut, "DtTU     %13.7f min\n", DtTU * 13.446851158);
    }

  /* ----- Iterate to find the optimum change of inclination ----- */
  DeltaiNew = Deltai;         // 1st guess, 0.01 to 0.025 seems good
  Deltai1   = 100.0;          // IF going to smaller orbit, should be
  NumIter   = 0;              // 1.0 - 0.025

  while (fabs(DeltaiNew - Deltai1) > 0.000001)
  {
    Deltai1 = DeltaiNew;
    DeltaVa = sqrt(V1t*V1t + VInit*VInit - 2.0 * V1t * VInit * cos(Deltai1));

    DeltaVb = sqrt(V3t*V3t + VFinal*VFinal - 
                   2.0 * V3t * VFinal * cos(TDi-Deltai1));
    DeltaiNew = asin((DeltaVa * VFinal * V3t * sin(TDi - Deltai1)) / 
                     (VInit * V1t * DeltaVb));
    NumIter++;
  }
  if (FileOut != NULL)
    fprintf(FileOut, "Iter Di %14.6f ø %3d %13.7f\n",
                      Deltai1 * Rad, (DeltaVa + DeltaVb) * 7905.365998);
}

void NodeOnlyChg
    (
      Real iInit, Real ecc, Real DeltaOmega, Real VInit, Real fpa, Real incl, 
      Real& iFinal, Real& DeltaV
    )
{
}

void NonCoplanarRendz
    (
      Real PhaseNew, Real Deltai, Real Delta2Node, Real LonTrue, Real RInt, 
      Real RTgt, Integer kTgt, Integer kInt,
      Real& TTrans, Real& TPhase, Real& DVPhase, Real& DVTrans1, Real& DVTrans2
    )
{
  const Real TwoPi = 2.0 * PI;
  const Real Rad   = 180.0 / PI;

  Real AngVelInt, AngVelTgt, ATrans, APhase, Lead, LeadNew, TNode, LonTrueNew,
       VInt, VTgt, VPhase, VTrans1, VTrans2;

  AngVelInt = sqrt(1.0 / (RInt * RInt * RInt));
  AngVelTgt = sqrt(1.0 / (RTgt * RTgt * RTgt));
  ATrans    = (RInt + RTgt) * 0.5;
  TTrans    = PI * sqrt(ATrans * ATrans * ATrans);

  Lead = AngVelTgt * TTrans;

  TNode = Delta2Node / AngVelInt;

  LonTrueNew = LonTrue + AngVelTgt * TNode; // fix     PhaseNew:= 13.5/Rad;
  LeadNew = PI + PhaseNew;

  TPhase = (LeadNew - Lead + TwoPi * kTgt) / AngVelTgt;

  APhase = Power((TPhase/(TwoPi * kInt)), 2.0 / 3.0);

  /* -----------------  Find Deltav's  ----------------- */
  VInt    = sqrt(1.0 / RInt);
  VPhase  = sqrt(2.0 / RInt - 1.0/APhase);
  DVPhase = VPhase - VInt;

  VTrans1  = sqrt(2.0 / RInt - 1.0 / ATrans);
  DVTrans1 = VTrans1 - VPhase;

  VTrans2  = sqrt(2.0 / RTgt - 1.0 / ATrans);
  VTgt     = sqrt(1.0 / RTgt);
  DVTrans2 = sqrt(VTgt * VTgt + VTrans2 * VTrans2 - 
                  2.0 * VTgt * VTrans2 * cos(Deltai));

  if ((Show == 'Y') || (Show == 'S'))
    if (FileOut != NULL)
    {
      fprintf(FileOut, "    A trans  = %12.8f ER\n", ATrans);
      fprintf(FileOut, "    Lead Angle N  = %12.8f ø\n", LeadNew * Rad);
      fprintf(FileOut, "    Time to Node  = %12.8f TU\n", TNode);
      fprintf(FileOut, "    Angle Init    = %12.8f ø\n", Lead * Rad);
      fprintf(FileOut, "    Phasenew      = %12.8f ø\n", PhaseNew * Rad);
      fprintf(FileOut, "    Lead Angle N  = %12.8f ø\n", LeadNew * Rad);
      fprintf(FileOut, "    DtTU Phaseing = %12.8f TU\n", TPhase);
      fprintf(FileOut, "    A phase       = %12.8f ER\n", APhase);
      fprintf(FileOut, "    vint tgt = %11.7f %11.7f\n", VInt, VTgt);
      fprintf(FileOut, "    vtrans1,2 = %11.7f %11.7f vp\n",
                        VTrans1, VTrans2);
      fprintf(FileOut, "    DtTU Transfer =  %12.8f TU\n", TTrans);
      fprintf(FileOut, "    AngVelTgt     =  %12.8f rad/TU kTgt %3d\n",
                        AngVelTgt, kTgt);
      fprintf(FileOut, "    LonTrueNew    = %12.8f ø\n", LonTrueNew * Rad);
      fprintf(FileOut, "    DeltaV  ph    = %12.8f %12.8f\n",
                        DVPhase, DVPhase * 7.905365);
      fprintf(FileOut, "    DeltaV  t1    = %12.8f %12.8f\n",
                        DVTrans1, DVTrans1 * 7.905365);
      fprintf(FileOut, "    DeltaV  t2    = %12.8f %12.8f\n",
                        DVTrans2, DVTrans2 * 7.905365);
    }
}

/*------------------------------------------------------------------------------
|
|                           PROCEDURE RENDEZVOUS
|
|  This PROCEDURE calculates parameters for a Hohmann transfer rendezvous.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs          Description                    Range / Units
|    Rcs1        - Radius of circular orbit int   ER
|    Rcs2        - Radius of circular orbit tgt   ER
|    eInit       - Ecc of first orbit
|    eFinal      - Ecc of final orbit
|    NuInit      - True Anomaly of first orbit    0 or Pi rad
|    NuFinal     - True Anomaly of final orbit    0 or Pi rad
|    PhaseI      - Initial phase angle (Tgt-Int)  +(ahead) or -(behind) rad
|    NumRevs     - Number of revs to wait
|    kTgt        -
|    kInt        -
|
|  OutPuts       :
|    PhaseF      - Final Phase Angle              rad
|    WaitTime    - Wait before next intercept opp TU
|    Deltav      - Change in velocity             ER/TU
|
|  Locals        :
|    DtTUTrans   - Time of flight of trans orbit  TU
|    ATrans      - Semimajor axis of trans orbit  ER
|    AngVelTgt   - Angular velocity of target     rad / TU
|    AngVelInt   - Angular velocity of int        rad / TU
|    LeadAng     - Lead Angle                     rad
|
|  Coupling      :
|    POWER       - Raise a base to a power
|
|  References    :
|    Vallado       2001, 343-350, Alg 44, Alg 45, Ex 6-8, Ex 6-9
|
 -----------------------------------------------------------------------------*/
void Rendezvous
    (
      Real Rcs1, Real Rcs3, Real PhaseI, Real eInit, Real eFinal, Real NuInit, 
      Real NuFinal, Integer kTgt, Integer kInt,
      Real& PhaseF, Real& WaitTime, Real& DeltaV
    )
{
  const Real TwoPi = 2.0 * PI;
  const Real Rad   = 180.0 / PI;

  Real PeriodTrans, Rp, DtTUTrans, LeadAng, aTrans, AngVelTgt, AngVelInt,
       a1, a2, a3, VInit, vTransa, VFinal, vTransb, SME1, SME2, SME3, DeltaVa,
       DeltaVb, VInt, VTrans;

  aTrans    = (Rcs1 + Rcs3) / 2.0;
  DtTUTrans = PI * sqrt(aTrans * aTrans * aTrans);
  AngVelInt = 1.0 / (sqrt(Rcs1 * Rcs1 * Rcs1));
  AngVelTgt = 1.0 / (sqrt(Rcs3 * Rcs3 * Rcs3));
  VInt      = sqrt(1.0 / Rcs1);

  /* ---------- Check for satellites in the same orbits ----------- */
  if (fabs(AngVelInt - AngVelTgt) < 0.000001)
  {
    PeriodTrans = (kTgt * TwoPi + PhaseI) / AngVelTgt;
    aTrans      = Power((PeriodTrans / (TwoPi * kInt)), 2.0 / 3.0 );
    Rp          = 2.0 * aTrans - Rcs1;
    if (Rp < 1.0)
      printf("Error - The transfer orbit intersects the Earth\n");
    VTrans   = sqrt((2.0 / Rcs1) - (1.0 / aTrans));
    DeltaV   = 2.0 * (VTrans - VInt);
    WaitTime = 0.0;
    if (FileOut != NULL)
      fprintf(FileOut, "tpha %11.7f  vint %11.f  phi %11.7f\n",
                        PeriodTrans, VInt, PhaseI * Rad);
    PhaseF   = PhaseI;
    WaitTime = PeriodTrans;
    LeadAng  = 0.0;
  }
  else
  {
    LeadAng  = AngVelTgt * DtTUTrans;
    PhaseF   = LeadAng - PI;
    WaitTime = (PhaseF - PhaseI + 2.0 * PI * kTgt) / (AngVelInt - AngVelTgt);

    a1   = (Rcs1 * (1.0 + eInit * cos(NuInit))) / (1.0 - eInit * eInit);
    a2   = ( Rcs1 + Rcs3) / 2.0;
    a3   = (Rcs3 * (1.0 + eFinal * cos(NuFinal))) / (1.0 - eFinal * eFinal);
    SME1 = -1.0 / (2.0 * a1);
    SME2 = -1.0 / (2.0 * a2);
    SME3 = -1.0 / (2.0 * a3);
    /* -----------------  Find Delta v at point a  ------------------ */
    VInit   = sqrt(2.0 * ((1.0 / Rcs1) + SME1));
    vTransa = sqrt(2.0 * ((1.0 / Rcs1) + SME2));
    DeltaVa = fabs(vTransa - VInit);

    /* -----------------  Find Delta v at point b  ------------------ */
    VFinal  = sqrt(2.0 * ((1.0 / Rcs3) + SME3));
    vTransb = sqrt(2.0 * ((1.0 / Rcs3) + SME2));
    DeltaVb = fabs(VFinal - vTransb);
    DeltaV  = DeltaVa + DeltaVb;

    if ((Show == 'Y') || (Show == 'S'))
      if (FileOut != NULL)
      {
        fprintf(FileOut, "    A transfer  = %12.8f ER\n", aTrans);
        fprintf(FileOut, "    DtTU Transfer= %12.8f TU %11.7f\n",
                          DtTUTrans, DtTUTrans * 13.44685115);
        fprintf(FileOut, "    AngVelTgt   = %12.8f rad/TU kTgt %3d %3d\n",
                          AngVelTgt, kTgt, kInt);
        fprintf(FileOut, "    AngVelInt   = %12.8f rad/TU\n", AngVelInt);
        fprintf(FileOut, "     Angle  = %12.8f ø\n", LeadAng * Rad);
      }
  }
}