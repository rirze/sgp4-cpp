#ifndef _ASTREDUC_
#define _ASTREDUC_
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

                                  INTERFACE

       ----------------------------------------------------------------      */
#include "astmath.h"

typedef struct timerec
{
  SINT Year;
  Byte Mon;
  Byte Day;
  Real DUT1;
  Byte DAT;
  Real xp;
  Real yp;
  Real XLOD;
  Real DDPsi;
  Real DDEps;
} TimeRec;

typedef struct datarec
{
  Real   JD;
  Real   DPsiFK5;
  Real   DEpsFk5;
  Real   DPsi96;
  Real   DEps96;
  Matrix PNmI;
} DataRec;

typedef SINT IAr5x106[5][106];
typedef Real RAr4x106[5][106];
typedef SINT IAr5x263[5][263];
typedef Real RAr6x263[6][263];
typedef Real RAr6x106[5][106];
typedef SINT IAr10x112[10][112];
typedef Real RAr4x112[4][112];

void ConvTime
    (
      char *, UINT, UINT, UINT, UINT, UINT, Real, UINT, char,
      Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&, Real&, char *
    );
void FK4
    (
      Vector&, Vector&, Direction, Vector&, Vector&
    );
void IAU2000GCRF
    (
      Vector, Vector, Direction, Vector&, Vector&, Real, Real, Real, Real,
      IAr5x106, RAr6x106
    );
void InitNutation
    (
      char *, char *, IAr5x106, RAr4x106, IAr5x106, RAr4x106,
      IAr5x263, RAr6x263, IAr10x112, RAr4x112
    );
void Nut80FK5
    (
      Vector&, Vector&, Direction, Vector&, Vector&, Real&, Real&, Real&,
      Real, IAr5x106, RAr4x106
    );
void Nut80GCRF
    (
      Vector&, Vector&, Direction, Vector&, Vector&, Real&, Real&, Real&,
      Real, Real, Real, IAr5x106, RAr4x106
    );
void Nut96FK5
    (
      Vector&, Vector&, Direction, Vector&, Vector&, Real&, Real&, Real&,
      Real, IAr5x263, RAr6x263, IAr10x112, RAr4x112
    );
void PolarM
    (
      Vector&, Vector&, Direction, Vector&, Vector&,Real, Real
    );
void Precession
    (
      Vector&, Vector&, Direction, Vector&, Vector&, Real
    );
void Sidereal
    (
      Vector&, Vector&, Direction, Vector&, Vector&,
      Real, Real, Real, Real, Real
    );
void TrueMean
    (
      Vector&, Vector&, Direction, Vector&, Vector&, Real&, Real&,
      Real, IAr5x106, RAr4x106
    );

#endif