#ifndef _eopspw_
#define _eopspw_
/*     ----------------------------------------------------------------
*
*                              eopspw.h
*
*  this file contains common routines for the obtaining earth orientation and
*  space weather data.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2007
*                            by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com
*
*     *****************************************************************
*
*  current :
*            21 mar 08  david vallado
*                           misc fixes
*  changes :
*            14 dec 05  david vallado
*                           misc fixes
*            21 oct 05  david vallado
*                           original baseline
*     *****************************************************************     */

#include <math.h>
#include <io.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "astmath.h"
#include "ast2body.h"
#include "asttime.h"
#include "msiscom.h"
#include "msis00.h"

/*    *****************************************************************
*     type definitions
*     *****************************************************************     */

// note the deltapsi/deltaeps are for iau76/fk5 from eop
// x/y/s/a and deltapsi/deltaeps are for computational speed in nutation calcs

typedef struct eopdata
  {
       double xp,   yp,  dut1, lod,  ddpsi,    ddeps,    dx,   dy;
       int    year, mon, day,  mjd,  dat;
       double x,    y,   s,    deltapsi, deltaeps;
  } eopdata;

typedef struct spwdata
  {
       double adjf10,  adjctrf81, adjlstf81, obsf10,   obsctrf81, obslstf81, cp;
       int    year,    mon,       day,       bsrn,     nd,        avgap,     c9,
              isn,     q,         aparr[8],  kparr[8], sumkp;
  } spwdata;

const int eopsize = 25000; // 25000 if from 62
const int spwsize = 2200; // 25000 if from 62



/*    *****************************************************************
*     routines
*     *****************************************************************    */

void initspw
     (
       spwdata spwarr[spwsize],
       double& jdspwstart
     );

void initeop
     (
       eopdata eoparr[eopsize],
       double& jdeopstart
     );

void findeopparam
     (
       double  jd,       double mfme,     char interp,
       eopdata eoparr[eopsize],           double jdeopstart,
       double& dut1,     int& dat,
       double& lod,      double& xp,      double& yp,
       double& ddpsi,    double& ddeps,   double& dx,   double& dy,
       double& x,        double& y,       double& s,
       double& deltapsi, double& deltaeps
     );

void findatmosparam
     (
       double jd, double mfme, char interp, char fluxtype, char f81type, char inputtype,
       spwdata spwarr[spwsize], double jdspwstart,
       double& f107, double& f107bar,
       double& ap, double& avgap, double aparr[8],
       double& kp, double& sumkp, double kparr[8]
     );

double kp2ap
       (
         double kpin
       );

double ap2kp
       (
         double apin
       );

// this routine will probably move elsewhere, but fits here ok for now
void interfaceatmos
     (
       double jde, double mfme, double recef[3],
       char interp, char fluxtype, char f81type, char inputtype,
       msistype& msis00r,
       spwdata spwarr[spwsize], double jdspwstart
     );

#endif

