#ifndef _ASTMANV_
#define _ASTMANV_
/*----------------------------------------------------------------      
                                                                             

                               UNIT ASTMANV;

                                                                             
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

                                  INTERFACE

       ----------------------------------------------------------------      */
#include "astmath.h"

/* ---------- Routines for Orbit Transfer calculations --------- */
void OneTangent(Real, Real, Real, Real, Real, Real, Real&, Real&, Real&);
void BiElliptic
    (
      Real, Real, Real, Real, Real, Real, Real, Real&, Real&, Real&, Real&
    );
void CombinedPlaneChg
    (
      Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,
      Real&, Real&, Real&
    );
void Cow2Hill(Vector, Vector, Vector, Vector, Vector&, Vector&);
void Hill2Cow(Vector, Vector, Vector, Vector, Vector&, Vector&);
void HillsR(Vector, Vector, Real, Real, Real, Vector&, Vector&);
void HillsV(Vector, Real, Real, Vector&);
void Hohmann(Real, Real, Real, Real, Real, Real, Real&, Real&, Real&);
void IandNodeChg(Real, Real, Real, Real, Real, Real&, Real&);
void IJK_RSW
    (
      Vector&, Vector&, Vector, Vector, Vector, Direction, Vector&, Vector&
    );
void IOnlyChg(Real, Real, Real, Real&);
void MinCombinedPlaneChg
    (
      Real, Real, Real, Real, Real, Real, Real, Real,
      Real&, Real&, Real&, Real&, Real&
    );
void NodeOnlyChg(Real, Real, Real, Real, Real, Real, Real& , Real&);
void NonCoplanarRendz
    (
      Real, Real, Real, Real, Real, Real,
      Integer, Integer,
      Real&, Real&, Real&, Real&, Real&
    );
void Rendezvous
    (
      Real, Real, Real, Real, Real, Real, Real,
      Integer, Integer,
      Real&, Real&, Real&
    );

#endif