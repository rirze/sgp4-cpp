#ifndef _astiod_h_
#define _astiod_h_
/* -----------------------------------------------------------------------------
*
*                              astIOD.h
*
*   this file contains fundamental astrodynamic procedures and functions
*   relating to the initial orbit determination techniques. see ch 7 for
*   a complete discussion of these routines.
*
*                           companion code for
*              fundamentals of astrodynamics and applications
*                                    2013
*                             by david vallado
*
*      (w) 719-573-2600, email dvallado@agi.com
*
*      *****************************************************************
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              21 mar 08  david vallado
*                           misc updates
*              14 feb 08  david vallado
*                           fix razel conversions
*              31 may 07  david vallado
*                           3rd edition baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
--------------------------------------------------------------------------  */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "D:/Codes/LIBRARY/CPP/Libraries/astMath/astMath/astMath.h"  // pi, infinite, undefined
#include "D:/Codes/LIBRARY/CPP/Libraries/astTime/astTime/astTime.h"  // pi, twopi, edirection
#include "D:/Codes/LIBRARY/CPP/Libraries/ast2Body/ast2Body/ast2Body.h"


namespace astIOD 
{
	void site
		(
		double latgd, double lon, double alt,
		double rsecef[3], double vsecef[3]
		);

	/* ------------------------ angles-only techniques -------------------------- */
	void anglesgauss
		(
		double tdecl1, double tdecl2, double tdecl3,
		double trtasc1, double trtasc2, double trtasc3,
		double jd1, double jd2, double jd3,
		double rs1[3], double rs2[3], double rs3[3], double r2[3], double v2[3]
		);

	//void angleslaplace
	//    (
	//      double, double, double, double, double, double, double, double, double,
	//      double[], double[], double[], double[], double[]
	//    );

	/* ------------------------- conversion techniques -------------------------- */
	void radec_azel
		(
		double& rtasc, double& decl, double& lst, double& latgd,
		edirection direct,
		double& az, double& el
		);

	void radec_elatlon
		(
		double& rtasc, double& decl,
		edirection direct,
		double& ecllat, double& ecllon
		);

	void rv_elatlon
		(
		double rijk[3], double vijk[3],
		edirection direct,
		double& rr, double& ecllat, double& ecllon,
		double& drr, double& decllat, double& decllon
		);

	void rv_radec
		(
		double rijk[3], double vijk[3], edirection direct,
		double& rr, double& rtasc, double& decl, double& drr, double& drtasc, double& ddecl
		);

	void rv_razel
		(
		double recef[3], double vecef[3], double rsecef[3], double latgd, double lon,
		edirection direct,
		double& rho, double& az, double& el, double& drho, double& daz, double& del
		);

	void rv_tradec
		(
		double rijk[3], double vijk[3], double rsijk[3],
		edirection direct,
		double& rho, double& trtasc, double& tdecl,
		double& drho, double& dtrtasc, double& dtdecl
		);

	void rvsez_razel
		(
		double rhosez[3], double drhosez[3],
		edirection direct,
		double& rho, double& az, double& el, double& drho, double& daz, double& del
		);

	void gibbs
		(
		double r1[3], double r2[3], double r3[3],
		double v2[3], double& theta, double& theta1, double& copa, char error[12]
		);

	void herrgibbs
		(
		double r1[3], double r2[3], double r3[3], double jd1, double jd2, double jd3,
		double v2[3], double& theta, double& theta1, double& copa, char error[12]
		);

	static double kbattin(double v);

	static double seebattin(double v);

	void lambertuniv
		(
		double ro[3], double r[3], char dm, char df, int nrev, double dtsec,
		double vo[3], double v[3], int& error, FILE *outfile
		);

	void lambertbattin
		(
		double r1[3], double r2[3], double v1[3], char dm, char df, int nrev, double dtsec,
		double v1t[3], double v2t[3], int& error, FILE *outfile
		);

	void target
		(
		double rint[3], double vint[3], double rtgt[3], double vtgt[3],
		char dm, char df, char kind, double dtsec,
		double v1t[3], double v2t[3], double dv1[3], double dv2[3], char error[12]
		);


}  // namespace 

#endif


