#ifndef _ASTTIME_H_
#define _ASTTIME_H_

/*     ----------------------------------------------------------------      
                                                                             

                               UNIT ASTTIME;

                                                                             
    This file contains fundamental Astrodynamic procedures and functions     
    relating to the time functions. These routines are discussed in Ch       
    3 and Ch 5.                                                              
                                                                             
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
#include "astutil.h"

/*
 * Object Definitions
*/
typedef char * Str3;
typedef char * Str10;
typedef char * Str11;
typedef char * Str12;

typedef enum {FROM, TOO} Direction;

/*
 * Global Interfaces
*/
void    DayLightSt(Integer, Real, Integer&, Integer&, Real&, Real&);
Integer DayOfWeek(Real);
void    Days2MDHMS(Integer, Real, 
                   Integer&, Integer&, Integer&, Integer&, Real&);
void    DMS_Rad(Integer&, Integer&, Real&, Direction, Real&);
void    FindDays(Integer, Integer, Integer, Integer, Integer, Real, Real&);
Byte    GetIntDay(Str3);
Byte    GetIntMon(Str3);
Real    GSTime(Real);
Real    GSTim0(Integer);
void    HMS_Rad(Integer&, Integer&, Real&, Direction, Real&);
void    HMS_Sec(Integer&, Integer&, Real&, Direction, Real&);
void    HMS_UT(Integer&, Integer&, Real&, Direction, Real&);
void    InitTime(void);
void    InvJulianDay(Real, Integer&, Integer&, Integer&, Integer&, Integer&,
                     Real&);
void    JulianDay(Integer, Integer, Integer, Integer, Integer, Real, Real&);
void    JulianDayAll(Integer, Integer, Integer, Integer, Integer, Real, char, 
                     Real&);
void    LSTime(Real, Real, Real&, Real&);
void    MoonRiseSet(Real, Real, Real, Real&, Real&, Real&, char*);
void    SunRiseSet(Real, Real, Real, char, Real&, Real&, char*);

#endif
