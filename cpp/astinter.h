/*----------------------------------------------------------------      
                                                                             

                               UNIT ASTINTER;

                                                                             
    This file contains fundamental Astrodynamic procedures and functions     
    relating to interplanetary calculations.                                 
                                                                             
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

void Interplanetary
    (
      Real, Real, Real, Real, Real, Real, Real, Real&, Real&, Real&, Real&
    );
void PlanetRV(Integer, char *, char *, Real, Vector&, Vector&);
