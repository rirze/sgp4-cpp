/*     ---------------------------------------------------------------- 

                             UNIT CONSTS;

  This file contains all of the constants and conversions used for         
  mathematical procedures.                                                 

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
     ---------------------------------------------------------------- */
#include <math.h>
#include "constants.h"

/*
 * Global Parameter References
*/
Real RadiusEarthNM, RadiusEarthFt, TUMin,         TUSec,     TUDay,
     OmegaEarthr,   OmegaEarth,    RadiusEarthKm, VFtPerSec, VKmPerSec,
     EESqrd,        DegPerSec,     Flat,          Mu,        MuKms,
     RadPerDay,     Rad,           TwoPi,         Pi;

/*-----------------------------------------------------------------------------
|
|                           PROCEDURE ConvertRNum
|
|  This procedure converts a number using one of the global conversion
|    factors.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs        :
|    RNum        - Number to convert
|    Numbr       - Code Number of the conversion                 1 to 11
|    WhichWay    - Which way to convert, divide or multiply        D M
|
|  OutPuts       :
|    RNum        - The converted number
|
|  Locals        :
|    None.
|
|  Coupling      :
|    GetConvNBR    Gets the conversion factor
|
 -----------------------------------------------------------------------------*/
void ConvertRNum(Real& RNum, Byte Numbr, Byte WhichWay)
{
  if (Numbr != 0)
    if (WhichWay == 'D')
      RNum = RNum / GetConvNbr(Numbr);
    else
      RNum = RNum * GetConvNbr(Numbr);
}

/*-----------------------------------------------------------------------------
|
|                           FUNCTION GETCONVNBR
|
|  This function gets the conversion value from the number which is passed in.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs        :
|    Numbr       - Code Number for the particular conversion 1 to 11
|
|  OutPuts       :
|                  The conversion number
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
Real GetConvNbr(Byte n)
{
  Real res;

  switch (n)
  {
    case 0:
      res = 1.0;  // Leave Unchanged
      break;
    case 1:
      res = RadiusEarthNM;
      break;
    case 2:
      res = RadiusEarthFt;
      break;
    case 3:
      res = RadiusEarthKm;
      break;
    case 4:
      res = TUMin;
      break;
    case 5:
      res = TUSec;
      break;
    case 6:
      res = Rad;
      break;
    case 7:
      res = OmegaEarth;
      break;
    case 8:
      res = VFtPerSec;
      break;
    case 9:
      res = VKmPerSec;
      break;
    case 10:
      res = EESqrd;
      break;
    case 11:
      res = DegPerSec;
      break;
    default:
      res = 1.0;  // Leave Unchanged
      break;
  }
  return res;
}

/* -----------------------------------------------------------------------------
|
|                           PROCEDURE WGS
|
|  This procedure initializes all the constants and conversions for a
|    program.
|
|  Author        : David Vallado                  303-344-6037    1 Mar 2001
|
|  Inputs        :
|    WhichOne    - Which WGS to use       2, 96, 84, 72
|
|  OutPuts       :
|    None.
|
|  Locals        :
|    None.
|
|  Coupling      :
|    None.
|
 -----------------------------------------------------------------------------*/
void WGS(WGS_Type WhichOne)
{
  switch (WhichOne)
  {
    case 2:
      // -- JGM-2, JGM-3, and EGM-96 values  --
      RadiusEarthKm = 6378.136300000;   // km
      Flat          = 1.0 / 298.257;
      OmegaEarthr   = 7.2921158553E-5;  // rad/sec
      MuKms         = 398600.4415;      // km3/s2
      break;
    case 72:
      // -- WGS 72 values --
      RadiusEarthKm = 6378.135000000;   // km
      Flat          = 1.0 / 298.26;
      OmegaEarthr   = 7.2921151470-5;   // rad/sec
      MuKms         = 398600.5;         // km3/s2
      break;
    case 84:
      // -- WGS 84 values --
      RadiusEarthKm = 6378.137000000;   // km
      Flat          = 1.0 / 298.2572235630;
      OmegaEarthr   = 7.2921158553E-5;  // rad/sec
      MuKms         = 398600.47;        // km3/s2
      break;
    case 96:
      // -- WGS-84/EGM-96 values --
      RadiusEarthKm = 6378.137000000;   // km
      Flat          = 1.0 / 298.2572235630;
      OmegaEarthr   = 7.2921158553E-5;  // rad/sec
      MuKms         = 398600.4418;      // km3/s2
      break;
    default:
      // -- JGM-2, JGM-3, and EGM-96 values  --
      RadiusEarthKm = 6378.136300000;   // km
      Flat          = 1.0 / 298.257;
      OmegaEarthr   = 7.2921158553E-5;  // rad/sec
      MuKms         = 398600.4415;      // km3/s2
      break;
  }

  TwoPi = 2.0 * PI;
  Rad   = 180.0 / PI;

  RadiusEarthNM = RadiusEarthKm / 1.852;
  RadiusEarthFt = RadiusEarthKm * 1000.0 / 0.3048;
  TUSec         = sqrt(RadiusEarthKm * RadiusEarthKm * RadiusEarthKm / MuKms);

  TUMin         = TUSec / 60.0;
  TUDay         = TUSec / 86400.0;
  OmegaEarth    = OmegaEarthr * TUSec;
  VKmPerSec     = sqrt(MuKms / RadiusEarthKm);
  VFtPerSec     = VKmPerSec * 1000.0 / 0.3048;
  EESqrd        = 2.0 * Flat - Flat * Flat;

  DegPerSec     = Rad / TUSec;      // This is actually BACKWARDS
  /* The conversion between deg/sec and rad/tu is opposite, but
     since the menu ONLY / on entry, and * on display, the constant
     is flipped.  Defined as TUSec/Rad, the correct way, the
     NORMAL conversion is 0.5 deg/sec / degpersec = rad/tu */
  RadPerDay     = 1.002737909350795 * TwoPi;  //OmegaEarth / TUDay
  Mu            = 1.0;     // in DU3/TU2
}
