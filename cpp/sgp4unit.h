#ifndef _SGP4UNIT_
#define _SGP4UNIT_
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

                                  INTERFACE

       ----------------------------------------------------------------      */
#include <math.h>
#include <stdio.h>

#include "astmath.h"
#include "constants.h"

typedef struct NearEarthRecord
{
  SINT Isimp, Method;
  Real Aycof,   CON41,  Cc1,      Cc4,   Cc5,   D2, D3, D4,    Delmo, Eta,
       ArgpDot, Omgcof, Sinmao,   T,     T2cof, T3cof,  T4cof, T5cof, X1mth2,
       X7thm1,  MDot,   OmegaDot, Xlcof, Xmcof, OmegaCF;
} NearEarthType;

typedef struct DeepSpaceRecord
{
  SINT IRez;
  Real D2201, D2211, D3210, D3222, D4410, D4422, D5220, D5232, D5421, D5433,
       Dedt,  Del1,  Del2,  Del3,  Didt,  Dmdt,  Dnodt, Domdt, E3,    Ee2,
       Peo,   Pgho,  Pho,   Pinco, Plo,   Se2,   Se3,   Sgh2,  Sgh3,  Sgh4,
       Sh2,   Sh3,   Si2,   Si3,   Sl2,   Sl3,   Sl4,   GSTo,  Xfact, Xgh2,
       Xgh3,  Xgh4,  Xh2,   Xh3,   Xi2,   Xi3,   Xl2,   Xl3,   Xl4,   Xlamo,
       Zmol,  Zmos,  Atime, Xli,   Xni;
} DeepSpaceType;

typedef struct ElSetRecord
{
  LINT SatNum;
  char Class;
  Byte EpochYr;
  NearEarthType NEValues;
  DeepSpaceType DSValues;
  Real a, Altp, Alta, EpochDays, JDSatEpoch, NDDot, NDot, BStar, RCSe,
       Inclo, Omegao, Ecco, Argpo, Mo, No;
  LINT EpochTyNumRev;
} ElSetRec;

void DPPer
    (
      Real, Real, Real, Real, Real, Real,
      Real, Real, Real, Real, Real, Real,
      Real, Real, Real, Real, Real, Real,
      Real, Real, Real, Real, Real, Real,
      Real, Real, Real, Real, Real, Real,
      Real, Real, 
      SINT,
      Real&, Real&, Real&, Real&, Real&
    );
void DSInit
    (
      Real, Real, Real, Real, Real, Real, Real, 
      Real, Real, Real, Real, Real, Real, Real, 
      Real, Real, Real, Real, Real, Real, Real, 
      Real, Real, Real, Real, Real, Real, Real, 
      Real, Real, Real, Real, Real, Real, Real, 
      Real, Real, Real, Real, Real, Real,
      SINT,
      Real&, Real&, Real&, Real&, Real&, 
      Real&, Real&, Real&, Real&, Real&, 
      Real&, Real&, Real&, Real&, Real&, 
      Real&, Real&, Real&, Real&, Real&, 
      Real&, Real&, Real&, Real&
    );
void DSCom
    (
      Real,  Real,  Real,  Real,  Real,  Real, Real,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&
    );
void DSpace
    (
      SINT, 
      Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
      Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
      Real, Real, 
      Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&

    );
void InitL
    (
      LINT, Real, Real, Real, Real&, SINT&,
      Real&, Real&, Real&, Real&, Real&, Real&, Real&, 
      Real&, Real&, Real&, Real&, Real&, Real&, Real&
    );
void SGP4(ElSetRec, Vector&, Vector&, SINT&);

void SGP4Init
    (
      LINT, SINT, Real, Real&, Real&, Real&, Real&, Real&, Real&, Real&,
      SINT&, NearEarthType&, DeepSpaceType&
    );

void TwoLine2RV
    (
      char*, char*, char, ElSetRec, Vector&, Vector&, FILE*, 
      char*, char*, char*, char*
    );

void TwoLine2RVSGP4
    (
      char *, char, ElSetRec, Vector&, Vector&, FILE*, char*, char*, char*
    );

#endif
