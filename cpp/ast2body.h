#ifndef _ast2body_h_
#define _ast2body_h_
/* --------------------------------------------------------------------
*
*                                ast2body.h
*
*    this file contains miscallaneous two-body motion functions.
*
*    current :
*               4 may 09  david vallado
*                           misc updates
*    changes :
*              23 feb 07  david vallado
*                           3rd edition baseline
*               6 dec 05  david vallado
*                           add ijk2ll
*              20 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
  ----------------------------------------------------------------------      */

#include <string.h>

#include "asttime.h"
#include "astmath.h"
    

void rv2coe
     (
       double r[3], double v[3],
       double& p, double& a, double& ecc, double& incl, double& omega, double& argp,
       double& nu, double& m, double& arglat, double& truelon, double& lonper
     );

void coe2rv
     (
       double p, double ecc, double incl, double omega, double argp, double nu,
       double arglat, double truelon, double lonper,
       double r[3],  double v[3]
     );

void findc2c3
     (
       double znew,
       double& c2new, double& c3new
     );

void kepler
     (
       double ro[3], double vo[3], double dtseco, double r[3], double v[3]
     );

void rv2rsw
     (
       double r[3], double v[3],
       double rrsw[3], double vrsw[3], std::vector< std::vector<double> > &transmat
     );

void  rv2ntw
      (
        double r[3], double v[3],
        double rntw[3], double vntw[3], std::vector< std::vector<double> > &transmat
      );

void newtonm
     (
       double ecc, double m, double& e0, double& nu
     );

void newtonnu
     (
       double ecc, double nu,
       double& e0, double& m
     );

void gc_gd
     (
       double&    latgc ,
       edirection direct,
       double&    latgd
     );

void ijk2ll
     (
       double recef[3], double jdut1,
       double& latgc, double& latgd, double& lon, double& hellp
     );

void sun
     (
       double jd,
       double rmoon[3], double& rtasc, double& decl
     );

void moon
     (
       double jd,
       double rmoon[3], double& rtasc, double& decl
     );


#endif

