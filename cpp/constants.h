#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_
/* ---------------------------------------------------------------- 


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
                                INTERFACE
   ---------------------------------------------------------------- 
*/

#include "astmath.h"

typedef enum {WGS2 = 2, WGS72 = 72, WGS84 = 84, WGS96 = 96} WGS_Type;

extern Real RadiusEarthNM, RadiusEarthFt, TUMin,         TUSec,     TUDay,
       OmegaEarthr,   OmegaEarth,    RadiusEarthKm, VFtPerSec, VKmPerSec,
       EESqrd,        DegPerSec,     Flat,          Mu,        MuKms,
       RadPerDay,     Rad,           TwoPi,         Pi;
/*
 * Forward References
*/
void ConvertRNum(Real&, Byte, Byte);
Real GetConvNbr(Byte);
void WGS(WGS_Type);

#endif
