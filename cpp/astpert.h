#ifndef _ASTPERT_
#define _ASTPERT_
/*     ----------------------------------------------------------------      
                                                                             

                               UNIT ASTPERT;

                                                                             
    This file contains fundamental Astrodynamic procedures and functions     
    allowing the calcualtions of perturbations. Most of these routines are   
    discussed in Ch 8.                                                       
                                                                             
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

void Atmos(Vector, Real&);
void Cowell
    (
      Vector, Vector, Real, Real, Real, char*, SINT, Real, Vector&, Vector&
    );
void Deriv(Matrix, Matrix&);
void FullGeop(Vector, Vector, Real, SINT, SINT, Real, Matrix, Matrix, Vector&);
void InitGravityField(SINT, Matrix&, Matrix&);
void J2DragPert(Real, Real, Real, Real, Real&, Real&, Real&);
void LegPoly(Real, Matrix&);
void PDeriv(Real, Matrix, char*, SINT, Real, Matrix&);
void PertAccel(Vector, Vector, Real, SINT, SINT, Real, Vector&);
void PKepler(Vector, Vector, Real, Real, Real, Vector&, Vector&);
void Predict
    (
      Real, Real, Real, Vector, Vector, Vector, char,
      Real&, Real&, Real&, Real&, Real&, char*
    );
void RK4(Real, Real, Matrix, char*, SINT, Real, Matrix&);
void RK45(Real, Real&, Matrix&, char*, SINT, Real, Matrix&);

#endif