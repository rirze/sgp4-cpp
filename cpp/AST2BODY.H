#ifndef _AST2BODY_H_
#define _AST2BODY_H_
/*     ----------------------------------------------------------------      
                                                                           

                               UNIT AST2BODY;

                                                                           
    This file contains fundamental Astrodynamic procedures and functions     
    using 2-body dynamics. The routines span a wide range of material, and   
    they come from chapters 2, 3, 5, and 11.                                 
                                                                             
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

       ----------------------------------------------------------------      */
#include "astmath.h"
#include "asttime.h"

void CheckHitEarth(Vector, Vector, Vector, Vector, char&);

void ElOrb
    (
      Vector, Vector,
      Real&, Real&, Real&, Real&, Real&, Real&,
      Real&, Real&, Real&, Real&, Real&
    );

void FindC2C3(Real, Real&, Real&);

void FindTOF(Vector, Vector, Real, Real&);

void GeocGeod(Real&, Direction, Real&);

void IJKtoLatLonA(Vector, Real, Real&, Real&, Real&, Real&);

void IJKtoLatLonB(Vector, Real, Real&, Real&, Real&, Real&);

void IJKtoLatLonE(Vector, Real, Real&, Real&, Real&, Real&);

void Light(Vector, Real, char, char *);

void Kepler(Vector, Vector, Real, Vector&, Vector&, char *);

void Moon(Real, Vector&, Real&, Real&);

void MoonIll(Real, Real, Real&);

void NewtonE(Real, Real, Real&, Real&);

void NewtonM(Real, Real, Real&, Real&);

void NewtonU(Real, Real, Real&, Real&);

void Path(Real, Real, Real, Real, Real&, Real&);

void RandV
    (
      Real, Real, Real, Real, Real, Real, Real, Real, Real,
      Vector&, Vector&
    );

void RngAz(Real, Real, Real, Real, Real, Real&, Real&);

void Sight(Vector, Vector, char, char *);

void Sun(Real, Vector&, Real&, Real&);

void SunIll(Real, Real, Real, Real&, Real&, Real&);

void SatFOV
    (
      Real, Real, Real, Real, Real, Real, Real,
      Real&, Real&, Real&, Real&, Real&, Real&
    );

#endif