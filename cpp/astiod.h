#ifndef _ASTIOD_H_
#define _ASTIOD_H_
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
                                                                             
    Current :                                                                
              14 May 01  David Vallado                                       
                           2nd edition baseline                              
    Changes :                                                                
              23 Nov 87  David Vallado                                       
                           Original baseline                                 
                                                                             
       ----------------------------------------------------------------      



       ----------------------------------------------------------------     

                                  INTERFACE

       ----------------------------------------------------------------      */
#include "astmath.h"
#include "asttime.h"

/* --------- only for testgau test ----------- */
void Site(Real, Real, Real, Vector&, Vector&);

/* ------------------- Angles-only techniques --------------------- */
void AnglesGauss
    (
      Real, Real, Real, Real, Real, Real, Real, Real, Real,
      Vector, Vector, Vector, Vector&, Vector&
    );

void AnglesLaplace
    (
      Real, Real, Real, Real, Real, Real, Real, Real, Real,
      Vector, Vector, Vector, Vector&, Vector&
    );

/* -------------------- Conversion techniques --------------------- */
void RaDec_AzEl(Real&, Real&, Real&, Real&, Direction, Real&, Real&);

void RaDec_ELatLon(Real&, Real&, Direction, Real&, Real&);

void RV_ELatLon
    (
      Vector&, Vector&, Direction, 
      Real&, Real&, Real&, Real&, Real&, Real&
    );

void RV_RaDec
    (
      Vector&, Vector&, Direction, 
      Real&, Real&, Real&, Real&, Real&, Real&
    );

void RV_RAzEl
    (
      Vector&, Vector&, Vector&, Real, Real, Direction, 
      Real&, Real&, Real&, Real&, Real&, Real&
    );

void RV_TRaDec
    (
      Vector&, Vector&, Vector&, Direction, 
      Real&, Real&, Real&, Real&, Real&, Real&
    );

void RVSez_RAzEl
    (
      Vector&, Vector&, Direction, 
      Real&, Real&, Real&, Real&, Real&, Real&
    );

/* -------------------- Three vector techniques ------------------- */
void Gibbs(Vector, Vector, Vector, Vector&, Real&, Real&, Real&, char*);

void HerrGibbs
    (
      Vector, Vector, Vector, Real, Real, Real, 
      Vector&, Real&, Real&, Real&, char*
    );

/* ----------------------- Lambert techniques -------------------- */
void LambertBattin(Vector, Vector, char, char, Real, Vector&, Vector&, char*);

void LambertUniv(Vector, Vector, char, char, Real, Vector&, Vector&, char*);

void Target
    (
      Vector, Vector, Vector, Vector, char, char, Real, 
      Vector&, Vector&, Vector&, Vector&, char*
    );

#endif
