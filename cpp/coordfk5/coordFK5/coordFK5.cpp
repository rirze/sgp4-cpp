/*       ----------------------------------------------------------------
*
*                              coordFK5.cpp
*
*  this file contains routines for iau-76/fk5 transformations.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2013
*                            by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              22 may 09  david vallado
*                           add all transformation
*              21 jan 08  david vallado
*                           fix matrix operations
*              30 may 07  david vallado
*                           3rd edition baseline
*              12 dec 05  david vallado
*                           misc updates
*               1 dec 05  david vallado
*                           misc to match final circular, proper compute order, etc.
*               5 may 05  david vallado
*                           misc fixes and clean up
*              24 mar 05  david vallado
*                           conversion to c++
*              21 feb 05  david vallado
*                           original baseline
*     *****************************************************************       */

#include "coordFK5.h"

namespace coordFK5
{

// simple function to store help levels that can be set on/off during execution
	void sethelp
		(
		char& iauhelp, char iauopt
		)
	{
		static char iaustore;

		if (iauopt != ' ')
		{
			iauhelp = iauopt;
			iaustore = iauopt;
		}
		else
			iauhelp = iaustore;
	}


/* ----------------------------------------------------------------------------
*
*                           function ddpsiddeps_dxdy
*
*  this function transforms the offset corrections between the ddpsi/ddeps and
*  the dx/dy.
*
*  author        : david vallado                  719-573-2600   11 nov 2005
*
*  revisions
*
*  inputs        : description                    range / units
*    ttt         - julian centuries of tt         centuries
*    ddpsi, ddeps- offsets to iau2000a            rad
*    direct      - direction of transfer          eFrom, eTo
*
*  outputs       :
*    dx, dy      - offsets to iau2000a            rad
*
*  locals        :
*    psia        - cannonical precession angle    rad
*    wa          - cannonical precession angle    rad
*    epsa        - cannonical precession angle    rad
*    chia        - cannonical precession angle    rad
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2007, 217
*    lieske et al. 1977, A&A 58, pp. 1-16
*    IERS Conventions 2000, chap. 5, herring et al. 2002, JGR 107, B4
* --------------------------------------------------------------------------- */

	void ddpsiddeps_dxdy
		(
		double ttt, double& ddpsi, double& ddeps,
		edirection direct,
		double& dx, double& dy
		)
	{
		double psia, chia, epsa, epsa1;
		const double a2r = 4.84813681109535993e-6;

		//     current values of precession
		psia = ((((-0.0000000951 * ttt + 0.000132851) * ttt - 0.00114045) * ttt - 1.0790069) * ttt + 5038.481507) * ttt; // "

		// value from website code - intermediate value
		epsa1 = 84381.448 - 46.8402 * ttt - 0.00059 * ttt * ttt + 0.001813 * ttt * ttt * ttt;
		epsa = ((((-0.0000000434 * ttt - 0.000000576) * ttt + 0.00200340) * ttt - 0.0001831) * ttt - 46.836769) * ttt + 84381.406;
		// use new value of eo or old one?

		chia = ((((-0.0000000560 * ttt + 0.000170663) * ttt - 0.00121197) * ttt - 2.3814292) * ttt + 10.556403) * ttt;

		psia = psia * a2r;
		chia = chia * a2r;
		epsa = epsa * a2r;

		if (direct == eTo)
		{
			//	Chapter 5, equation (23)
			dx = ddpsi * sin(epsa) + ddeps * (psia * cos(epsa) - chia);
			dy = ddeps - ddpsi * sin(epsa) * (psia * cos(epsa) - chia);
		}
		else
		{
			//	Chapter 5, equation (23) reversed
			ddpsi = dx / sin(epsa) - dy * (psia * cos(epsa) - chia) / sin(epsa);
			ddeps = dy + dx  *  (psia  *  cos(epsa) - chia);
		}
	}   // procedure ddpsiddeps_dxdy



/* -----------------------------------------------------------------------------
*
*                           function iau80in
*
*  this function initializes the nutation matricies needed for reduction
*    calculations. the routine needs the filename of the files as input.
*
*  author        : david vallado                  719-573-2600   27 may 2002
*
*  revisions
*    vallado     - conversion to c++                             21 feb 2005
*
*  inputs          description                    range / units
*
*  outputs       :
*    iau80rec    - record containing the iau80 constants rad
*
*  locals        :
*    convrt      - conversion factor milli arcsec to radians
*    i,j         - index
*
*  coupling      :
*    none        -
*
*  references    :
* --------------------------------------------------------------------------- */

	void iau80in
		(
		iau80data& iau80rec,
		char EopLoc[85]
		)
	{
		FILE *infile;
		double convrt;
		int i, j, ret;

		// ------------------------  implementation   -------------------
		convrt = 0.0001 * pi / (180 * 3600.0); 	// 0.0001" to rad

#ifdef _MSC_VER
		fopen_s(&infile, EopLoc, "r");
#else
		infile = fopen(EopLoc, "r");
#endif

		for (i = 1; i <= 106; i++)
		{
#ifdef _MSC_VER
			ret = fscanf_s(infile, "%d %d %d %d %d %lf %lf %lf %lf %d \n ",
				&iau80rec.iar80[i][1], &iau80rec.iar80[i][2], &iau80rec.iar80[i][3],
				&iau80rec.iar80[i][4], &iau80rec.iar80[i][5],
				&iau80rec.rar80[i][1], &iau80rec.rar80[i][2], &iau80rec.rar80[i][3],
				&iau80rec.rar80[i][4], &j);
#else
			ret = fscanf(infile, "%d %d %d %d %d %lf %lf %lf %lf %d \n ",
				&iau80rec.iar80[i][1], &iau80rec.iar80[i][2], &iau80rec.iar80[i][3],
				&iau80rec.iar80[i][4], &iau80rec.iar80[i][5],
				&iau80rec.rar80[i][1], &iau80rec.rar80[i][2], &iau80rec.rar80[i][3],
				&iau80rec.rar80[i][4], &j);
#endif
			if (ret == EOF)
			{
				break;      /* get out of loop reading lines – found end of file prematurely */
			}

			for (j = 1; j <= 4; j++)
				iau80rec.rar80[i][j] = iau80rec.rar80[i][j] * convrt;
		}
		fclose(infile);
	}  // procedure iau80in

/* -----------------------------------------------------------------------------
*
*                           function fundarg
*
*  this function calulates the delauany variables and planetary values for
*  several theories.
*
*  author        : david vallado                  719-573-2600   16 jul 2004
*
*  revisions
*    vallado     - conversion to c++                             23 nov 2005
*
*  inputs          description                    range / units
*    ttt         - julian centuries of tt
*    opt         - method option                  e00a, e00b, e96, e80
*
*  outputs       :
*    l           - delaunay element               rad
*    ll          - delaunay element               rad
*    f           - delaunay element               rad
*    d           - delaunay element               rad
*    omega       - delaunay element               rad
*    planetary longitudes                         rad
*
*  locals        :
*
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2013, 225, Eq 3-82
* --------------------------------------------------------------------------- */

	void fundarg
		(
		double ttt, eOpt opt,
		double& l, double& l1, double& f, double& d, double& omega,
		double& lonmer, double& lonven, double& lonear, double& lonmar,
		double& lonjup, double& lonsat, double& lonurn, double& lonnep,
		double& precrate
		)
	{
		double deg2rad, oo3600;
		char iauhelp;

		sethelp(iauhelp, ' ');
		deg2rad = pi / 180.0;
		oo3600 = 1.0 / 3600.0;

		// ---- determine coefficients for icrs nutation theories ----
		// ---- iau-2000a theory
		if (opt == e00a)
		{
			// ------ form the delaunay fundamental arguments in deg
			l = ((((-0.00024470 * ttt + 0.051635) * ttt + 31.8792) * ttt + 1717915923.2178) * ttt + 485868.249036) * oo3600;
			l1 = ((((-0.00001149 * ttt + 0.000136) * ttt - 0.5532) * ttt + 129596581.0481) * ttt + 1287104.793048) * oo3600;
			f = ((((+0.00000417 * ttt - 0.001037) * ttt - 12.7512) * ttt + 1739527262.8478) * ttt + 335779.526232) * oo3600;
			d = ((((-0.00003169 * ttt + 0.006593) * ttt - 6.3706) * ttt + 1602961601.2090) * ttt + 1072260.703692) * oo3600;
			omega = ((((-0.00005939 * ttt + 0.007702) * ttt + 7.4722) * ttt - 6962890.5431) * ttt + 450160.398036) * oo3600;

			// ------ form the planetary arguments in deg
			lonmer = (908103.259872 + 538101628.688982  * ttt) * oo3600;
			lonven = (655127.283060 + 210664136.433548  * ttt) * oo3600;
			lonear = (361679.244588 + 129597742.283429  * ttt) * oo3600;
			lonmar = (1279558.798488 + 68905077.493988  * ttt) * oo3600;
			lonjup = (123665.467464 + 10925660.377991  * ttt) * oo3600;
			lonsat = (180278.799480 + 4399609.855732  * ttt) * oo3600;
			lonurn = (1130598.018396 + 1542481.193933  * ttt) * oo3600;
			lonnep = (1095655.195728 + 786550.320744  * ttt) * oo3600;
			precrate = ((1.112022 * ttt + 5028.8200) * ttt) * oo3600;
		}

		// ---- iau-2000b theory
		if (opt == e00b)
		{
			// ------ form the delaunay fundamental arguments in deg
			l = (1717915923.2178  * ttt + 485868.249036) * oo3600;
			l1 = (129596581.0481  * ttt + 1287104.79305) * oo3600;
			f = (1739527262.8478  * ttt + 335779.526232) * oo3600;
			d = (1602961601.2090  * ttt + 1072260.70369) * oo3600;
			omega = (-6962890.5431  * ttt + 450160.398036) * oo3600;

			// ------ form the planetary arguments in deg
			lonmer = 0.0;
			lonven = 0.0;
			lonear = 0.0;
			lonmar = 0.0;
			lonjup = 0.0;
			lonsat = 0.0;
			lonurn = 0.0;
			lonnep = 0.0;
			precrate = 0.0;
		}

		// ---- iau-1996 theory
		if (opt == e96)
		{
			// ------ form the delaunay fundamental arguments in deg
			l = ((((-0.00024470 * ttt + 0.051635) * ttt + 31.8792) * ttt + 1717915923.2178) * ttt) * oo3600 + 134.96340251;
			l1 = ((((-0.00001149 * ttt - 0.000136) * ttt - 0.5532) * ttt + 129596581.0481) * ttt) * oo3600 + 357.52910918;
			f = ((((+0.00000417 * ttt + 0.001037) * ttt - 12.7512) * ttt + 1739527262.8478) * ttt) * oo3600 + 93.27209062;
			d = ((((-0.00003169 * ttt + 0.006593) * ttt - 6.3706) * ttt + 1602961601.2090) * ttt) * oo3600 + 297.85019547;
			omega = ((((-0.00005939 * ttt + 0.007702) * ttt + 7.4722) * ttt - 6962890.2665) * ttt) * oo3600 + 125.04455501;
			// ------ form the planetary arguments in deg
			lonmer = 0.0;
			lonven = 181.979800853 + 58517.8156748   * ttt;
			lonear = 100.466448494 + 35999.3728521   * ttt;
			lonmar = 355.433274605 + 19140.299314    * ttt;
			lonjup = 34.351483900 + 3034.90567464  * ttt;
			lonsat = 50.0774713998 + 1222.11379404  * ttt;
			lonurn = 0.0;
			lonnep = 0.0;
			precrate = (0.0003086 * ttt + 1.39697137214) * ttt;
		}

		// ---- iau-1980 theory
		if (opt == e80)
		{
			// ------ form the delaunay fundamental arguments in deg
			l = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) * oo3600 + 134.96298139;
			l1 = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) * oo3600 + 357.52772333;
			f = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) * oo3600 + 93.27191028;
			d = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) * oo3600 + 297.85036306;
			omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) * oo3600 + 125.04452222;
			// ------ form the planetary arguments in deg
			// iers tn13 shows no planetary
			// seidelmann shows these equations
			// circ 163 shows no planetary
			// ???????
			lonmer = 252.3 + 149472.0  * ttt;
			lonven = 179.9 + 58517.8  * ttt;
			lonear = 98.4 + 35999.4  * ttt;
			lonmar = 353.3 + 19140.3  * ttt;
			lonjup = 32.3 + 3034.9  * ttt;
			lonsat = 48.0 + 1222.1  * ttt;
			lonurn = 0.0;
			lonnep = 0.0;
			precrate = 0.0;
		}

		// ---- convert units from deg to rad
		l = fmod(l, 360.0)      *  deg2rad;
		l1 = fmod(l1, 360.0)     *  deg2rad;
		f = fmod(f, 360.0)      *  deg2rad;
		d = fmod(d, 360.0)      *  deg2rad;
		omega = fmod(omega, 360.0)  *  deg2rad;

		lonmer = fmod(lonmer, 360.0) * deg2rad;
		lonven = fmod(lonven, 360.0) * deg2rad;
		lonear = fmod(lonear, 360.0) * deg2rad;
		lonmar = fmod(lonmar, 360.0) * deg2rad;
		lonjup = fmod(lonjup, 360.0) * deg2rad;
		lonsat = fmod(lonsat, 360.0) * deg2rad;
		lonurn = fmod(lonurn, 360.0) * deg2rad;
		lonnep = fmod(lonnep, 360.0) * deg2rad;
		precrate = fmod(precrate, 360.0) * deg2rad;

		if (iauhelp == 'y')
		{
			printf("fa %11.7f  %11.7f  %11.7f  %11.7f  %11.7f deg \n", l * 180 / pi, l1 * 180 / pi, f * 180 / pi, d * 180 / pi, omega * 180 / pi);
			printf("fa %11.7f  %11.7f  %11.7f  %11.7f deg \n", lonmer * 180 / pi, lonven * 180 / pi, lonear * 180 / pi, lonmar * 180 / pi);
			printf("fa %11.7f  %11.7f  %11.7f  %11.7f deg \n", lonjup * 180 / pi, lonsat * 180 / pi, lonurn * 180 / pi, lonnep * 180 / pi);
			printf("fa %11.7f  \n", precrate * 180 / pi);
		}
	}  // procedure fundarg

/* -----------------------------------------------------------------------------
*
*                           function precess
*
*  this function calulates the transformation matrix that accounts for the effects
*    of precession. both the 1980 and 2000 theories are handled. note that the
*    required parameters differ a little.
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  revisions
*    vallado     - conversion to c++                             21 feb 2005
*    vallado     - misc updates, nomenclature, etc               23 nov 2005
*
*  inputs          description                    range / units
*    ttt         - julian centuries of tt
*    opt         - method option                  e00a, e00b, e96, e80
*
*  outputs       :
*    prec        - transformation matrix for mod - j2000 (80 only)
*    psia        - cannonical precession angle    rad    (00 only)
*    wa          - cannonical precession angle    rad    (00 only)
*    epsa        - cannonical precession angle    rad    (00 only)
*    chia        - cannonical precession angle    rad    (00 only)
*    prec        - matrix converting from "mod" to gcrf
*
*  locals        :
*    convrt      - conversion factor arcsec to radians
*    zeta        - precession angle               rad
*    z           - precession angle               rad
*    theta       - precession angle               rad
*    oblo        - obliquity value at j2000 epoch "//
*
*  coupling      :
*    none        -
*
*  references    :
*    vallado       2013, 213, 226
* --------------------------------------------------------------------------- */

	void precess
		(
		double ttt, eOpt opt,
		double& psia, double& wa, double& epsa, double& chia,
		std::vector< std::vector<double> > &prec
		)
	{
		prec.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = prec.begin(); it != prec.end(); ++it)
			it->resize(3);
		std::vector< std::vector<double> > p1, p2, p3, p4, tr1, tr2;
		double convrt, zeta, theta, z, coszeta, sinzeta, costheta, sintheta,
			cosz, sinz, oblo;
		char iauhelp;
		sethelp(iauhelp, ' ');

		convrt = pi / (180.0 * 3600.0);

		// ------------------- iau 77 precession angles --------------------
		if ((opt == e80) | (opt == e96))
		{
			oblo = 84381.448; // "
			psia = ((-0.001147 * ttt - 1.07259) * ttt + 5038.7784) * ttt; // "
			wa = ((-0.007726 * ttt + 0.05127) * ttt) + oblo;
			epsa = ((0.001813 * ttt - 0.00059) * ttt - 46.8150) * ttt + oblo;
			chia = ((-0.001125 * ttt - 2.38064) * ttt + 10.5526) * ttt;

			zeta = ((0.017998 * ttt + 0.30188) * ttt + 2306.2181) * ttt; // "
			theta = ((-0.041833 * ttt - 0.42665) * ttt + 2004.3109) * ttt;
			z = ((0.018203 * ttt + 1.09468) * ttt + 2306.2181) * ttt;
		}
		// ------------------ iau 00 precession angles -------------------
		else
		{
			oblo = 84381.406; // "
			psia = ((((-0.0000000951 * ttt + 0.000132851) * ttt - 0.00114045) * ttt - 1.0790069) * ttt + 5038.481507) * ttt; // "
			wa = ((((0.0000003337 * ttt - 0.000000467) * ttt - 0.00772503) * ttt + 0.0512623) * ttt - 0.025754) * ttt + oblo;
			epsa = ((((-0.0000000434 * ttt - 0.000000576) * ttt + 0.00200340) * ttt - 0.0001831) * ttt - 46.836769) * ttt + oblo;
			chia = ((((-0.0000000560 * ttt + 0.000170663) * ttt - 0.00121197) * ttt - 2.3814292) * ttt + 10.556403) * ttt;

			zeta = ((((-0.0000003173 * ttt - 0.000005971) * ttt + 0.01801828) * ttt + 0.2988499) * ttt + 2306.083227) * ttt + 2.650545; // "
			theta = ((((-0.0000001274 * ttt - 0.000007089) * ttt - 0.04182264) * ttt - 0.4294934) * ttt + 2004.191903) * ttt;
			z = ((((0.0000002904 * ttt - 0.000028596) * ttt + 0.01826837) * ttt + 1.0927348) * ttt + 2306.077181) * ttt - 2.650545;
		}

		// convert units to rad
		psia = psia  * convrt;
		wa = wa    * convrt;
		oblo = oblo  * convrt;
		epsa = epsa  * convrt;
		chia = chia  * convrt;

		zeta = zeta  * convrt;
		theta = theta * convrt;
		z = z     * convrt;

		if ((opt == e80) | (opt == e96))
		{
			coszeta = cos(zeta);
			sinzeta = sin(zeta);
			costheta = cos(theta);
			sintheta = sin(theta);
			cosz = cos(z);
			sinz = sin(z);

			// ----------------- form matrix  mod to gcrf ------------------
			prec[0][0] = coszeta * costheta * cosz - sinzeta * sinz;
			prec[0][1] = coszeta * costheta * sinz + sinzeta * cosz;
			prec[0][2] = coszeta * sintheta;
			prec[1][0] = -sinzeta * costheta * cosz - coszeta * sinz;
			prec[1][1] = -sinzeta * costheta * sinz + coszeta * cosz;
			prec[1][2] = -sinzeta * sintheta;
			prec[2][0] = -sintheta * cosz;
			prec[2][1] = -sintheta * sinz;
			prec[2][2] = costheta;

			// alternate approach
			// p1 = rot3mat( z );
			// p2 = rot2mat( -theta );
			// p3 = rot3mat( zeta );
			// prec = p3 * p2*p1;
		}
		else
		{
			astMath::rot3mat(-chia, p1);
			astMath::rot1mat(wa, p2);
			astMath::rot3mat(psia, p3);
			astMath::rot1mat(-oblo, p4);
			astMath::matmult(p4, p3, tr1, 3, 3, 3);
			astMath::matmult(tr1, p2, tr2, 3, 3, 3);
			astMath::matmult(tr2, p1, prec, 3, 3, 3);
		}

		if (iauhelp == 'y')
		{
			printf("pr %11.7f  %11.7f  %11.7f %11.7fdeg \n", psia * 180 / pi, wa * 180 / pi, epsa * 180 / pi, chia * 180 / pi);
			printf("pr %11.7f  %11.7f  %11.7fdeg \n", zeta * 180 / pi, theta * 180 / pi, z * 180 / pi);
		}
	} // procedure precess

/* -----------------------------------------------------------------------------
*
*                           function nutation
*
*  this function calulates the transformation matrix that accounts for the
*    effects of nutation.
*
*  author        : david vallado                  719-573-2600   27 jun 2002
*
*  revisions
*    vallado     - consolidate with iau 2000                     14 feb 2005
*    vallado     - conversion to c++                             21 feb 2005
*
*  inputs          description                    range / units
*    ttt         - julian centuries of tt
*    ddpsi       - delta psi correction to gcrf   rad
*    ddeps       - delta eps correction to gcrf   rad
*    iau80rec    - record containing the iau80 constants rad
*
*  outputs       :
*    deltapsi    - nutation in longitude angle   rad
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    nut         - transform matrix for tod -     mod
*
*  locals        :
*    iar80       - integers for fk5 1980
*    rar80       - reals for fk5 1980
*    l           -                                rad
*    ll          -                                rad
*    f           -                                rad
*    d           -                                rad
*    deltaeps    - change in obliquity            rad
*
*  coupling      :
*    fundarg     - find fundamental arguments
*    fmod      - modulus division
*
*  references    :
*    vallado       2013, 213, 224
* --------------------------------------------------------------------------- */

	void nutation
		(
		double ttt, double ddpsi, double ddeps,
		const iau80data &iau80rec,
		double& deltapsi, double& deltaeps, double& trueeps, double& meaneps, double& omega,
		std::vector< std::vector<double> > &nut
		)
	{
		nut.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = nut.begin(); it != nut.end(); ++it)
			it->resize(3);
		// locals
		double deg2rad, l, l1, f, d,
			lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate,
			cospsi, sinpsi, coseps, sineps, costrueeps, sintrueeps;
		int  i;
		double tempval;

		char iauhelp;
		sethelp(iauhelp, ' ');

		deg2rad = pi / 180.0;

		// ---- determine coefficients for iau 1980 nutation theory ----
		meaneps = ((0.001813  * ttt - 0.00059) * ttt - 46.8150) * ttt + 84381.448;
		meaneps = fmod(meaneps / 3600.0, 360.0);
		meaneps = meaneps  *  deg2rad;

		fundarg(ttt, e80, l, l1, f, d, omega,
			lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);

		deltapsi = 0.0;
		deltaeps = 0.0;
		for (i = 106; i >= 1; i--)
		{
			tempval = iau80rec.iar80[i][1] * l + iau80rec.iar80[i][2] * l1 + iau80rec.iar80[i][3] * f +
				iau80rec.iar80[i][4] * d + iau80rec.iar80[i][5] * omega;
			deltapsi = deltapsi + (iau80rec.rar80[i][1] + iau80rec.rar80[i][2] * ttt)  *  sin(tempval);
			deltaeps = deltaeps + (iau80rec.rar80[i][3] + iau80rec.rar80[i][4] * ttt) * cos(tempval);
		}

		// --------------- find nutation parameters --------------------
		deltapsi = fmod(deltapsi + ddpsi, twopi);
		deltaeps = fmod(deltaeps + ddeps, twopi);

		trueeps = meaneps + deltaeps;

		cospsi = cos(deltapsi);
		sinpsi = sin(deltapsi);
		coseps = cos(meaneps);
		sineps = sin(meaneps);
		costrueeps = cos(trueeps);
		sintrueeps = sin(trueeps);

		nut[0][0] = cospsi;
		nut[0][1] = costrueeps * sinpsi;
		nut[0][2] = sintrueeps * sinpsi;
		nut[1][0] = -coseps * sinpsi;
		nut[1][1] = costrueeps * coseps * cospsi + sintrueeps * sineps;
		nut[1][2] = sintrueeps * coseps * cospsi - sineps * costrueeps;
		nut[2][0] = -sineps * sinpsi;
		nut[2][1] = costrueeps * sineps * cospsi - sintrueeps * coseps;
		nut[2][2] = sintrueeps * sineps * cospsi + costrueeps * coseps;

		// alternate approach
		// n1 = rot1mat( trueeps );
		// n2 = rot3mat( deltapsi );
		// n3 = rot1mat( -meaneps );
		// nut = n3 * n2 * n1;

		if (iauhelp == 'y')
			printf("meaneps %11.7f dp  %11.7f de  %11.7f te  %11.7f  \n", meaneps * 180 / pi, deltapsi * 180 / pi,
		         	deltaeps * 180 / pi, trueeps * 180 / pi);
	}  // procedure nutation

/* -----------------------------------------------------------------------------
*
*                           function sidereal
*
*  this function calulates the transformation matrix that accounts for the
*    effects of sidereal time. Notice that deltaspi should not be moded to a
*    positive number because it is multiplied rather than used in a
*    trigonometric argument.
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  revisions
*    vallado     - fix units on kinematic terms                   5 sep 2002
*    vallado     - add terms                                     30 sep 2002
*    vallado     - consolidate with iau 2000                     14 feb 2005
*    vallado     - conversion to c++                             21 feb 2005
*
*  inputs          description                    range / units
*    jdut1       - julian centuries of ut1        days
*    deltapsi    - nutation angle                 rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       - long of asc node of moon       rad
*    lod         - length of day                  sec
*    eqeterms    - terms for ast calculation      0,2
*
*  outputs       :
*    st          - transformation matrix for pef - tod
*    stdot       - transformation matrix for pef - tod rate
*
*  locals        :
*    gmst         - mean greenwich sidereal time   0 to 2pi rad
*    ast         - apparent gmst                   0 to 2pi rad
*    hr          - hour                           hr
*    min         - minutes                        min
*    sec         - seconds                        sec
*    temp        - temporary vector
*    tempval     - temporary variable
*
*  coupling      :
*
*  references    :
*    vallado       2013, 212, 223
* --------------------------------------------------------------------------- */

	void sidereal
		(
		double jdut1, double deltapsi, double meaneps, double omega,
		double lod, int eqeterms,
		std::vector< std::vector<double> > &st,
		std::vector< std::vector<double> > &stdot
		)
	{
		st.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = st.begin(); it != st.end(); ++it)
			it->resize(3);
		stdot.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = stdot.begin(); it != stdot.end(); ++it)
			it->resize(3);
		// locals
		double gmst, ast, thetasa, omegaearth, sinast, cosast;
		char iauhelp;

		sethelp(iauhelp, ' ');

		// ------------------------ find gmst --------------------------
		gmst = astTime::gstime(jdut1);

		// ------------------------ find mean ast ----------------------
		if ((jdut1 > 2450449.5) && (eqeterms > 0))
		{
			ast = gmst + deltapsi *  cos(meaneps) + 0.00264 * pi / (3600 * 180) * sin(omega)
				+ 0.000063 * pi / (3600 * 180) * sin(2.0  * omega);
		}
		else
			ast = gmst + deltapsi *  cos(meaneps);

		ast = fmod(ast, 2.0 * pi);

		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth = thetasa;

		sinast = sin(ast);
		cosast = cos(ast);

		st[0][0] = cosast;
		st[0][1] = -sinast;
		st[0][2] = 0.0;
		st[1][0] = sinast;
		st[1][1] = cosast;
		st[1][2] = 0.0;
		st[2][0] = 0.0;
		st[2][1] = 0.0;
		st[2][2] = 1.0;

		// compute sidereal time rate matrix
		stdot[0][0] = -omegaearth * sinast;
		stdot[0][1] = -omegaearth * cosast;
		stdot[0][2] = 0.0;
		stdot[1][0] = omegaearth * cosast;
		stdot[1][1] = -omegaearth * sinast;
		stdot[1][2] = 0.0;
		stdot[2][0] = 0.0;
		stdot[2][1] = 0.0;
		stdot[2][2] = 0.0;

		if (iauhelp == 'y')
			printf("st gmst %11.7f ast %11.7f ome  %11.7f \n", gmst * 180 / pi, ast * 180 / pi, omegaearth * 180 / pi);

	}  // procedure sidereal

/* -----------------------------------------------------------------------------
*
*                           function polarm
*
*  this function calulates the transformation matrix that accounts for polar
*    motion. both the 1980 and 2000 theories are handled. note that the rotation
*    order is different between 1980 and 2000 .
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  revisions
*    vallado     - conversion to c++                             23 nov 2005
*
*  inputs          description                    range / units
*    xp          - polar motion coefficient       rad
*    yp          - polar motion coefficient       rad
*    ttt         - julian centuries of tt (00 theory only)
*    opt         - method option                  e00a, e00b, e96, e80
*
*  outputs       :
*    pm          - transformation matrix for itrf - pef
*
*  locals        :
*    convrt      - conversion from arcsec to rad
*    sp          - s prime value
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2013, 212, 223
* --------------------------------------------------------------------------- */

	void polarm
		(
		double xp, double yp, double ttt, eOpt opt, std::vector< std::vector<double> > &pm
		)
	{
		pm.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = pm.begin(); it != pm.end(); ++it)
			it->resize(3);
		double convrt, cosxp, cosyp, sinxp, sinyp, sp, cossp, sinsp;

		//convrt = pi / (180.0 * 3600.0);
		convrt = 1.0;  // do when read data in one time

		cosxp = cos(xp * convrt);
		sinxp = sin(xp * convrt);
		cosyp = cos(yp * convrt);
		sinyp = sin(yp * convrt);

		if ((opt == e80) | (opt == e96))
		{
			pm[0][0] = cosxp;
			pm[0][1] = 0.0;
			pm[0][2] = -sinxp;
			pm[1][0] = sinxp  *  sinyp;
			pm[1][1] = cosyp;
			pm[1][2] = cosxp  *  sinyp;
			pm[2][0] = sinxp  *  cosyp;
			pm[2][1] = -sinyp;
			pm[2][2] = cosxp  *  cosyp;

			// alternate approach
			// a1 = rot2mat(xp);
			// a2 = rot1mat(yp);
			// pm = a2 * a1;
		}
		else
		{
			// approximate sp value in rad
			sp = -47.0e-6 * ttt * pi / (180.0 * 3600.0);
			cossp = cos(sp);
			sinsp = sin(sp);

			// form the matrix
			pm[0][0] = cosxp * cossp;
			pm[0][1] = -cosyp * sinsp + sinyp * sinxp * cossp;
			pm[0][2] = -sinyp * sinsp - cosyp * sinxp * cossp;
			pm[1][0] = cosxp * sinsp;
			pm[1][1] = cosyp * cossp + sinyp * sinxp * sinsp;
			pm[1][2] = sinyp * cossp - cosyp * sinxp * sinsp;
			pm[2][0] = sinxp;
			pm[2][1] = -sinyp * cosxp;
			pm[2][2] = cosyp * cosxp;

			// alternate approach
			// a1 = rot1mat(yp);
			// a2 = rot2mat(xp);
			// a3 = rot3mat(-sp);
			// pm = a3 * a2 * a1;
		}
	}  // procedure polarm

/* -----------------------------------------------------------------------------
*
*                           function framebias
*
*  this function calulates the transformation matrix that accounts for frame
*    bias.
*
*  author        : david vallado                  719-573-2600   19ep05
*
*  revisions
*
*  inputs          description                    range / units
*    opt         - frame bias method option       'j' j2000, 'f' fk5
*
*  outputs       :
*    term1       - alpha delta o                  rad
*    term2       - psi deltab sin(eps deltao)     rad
*    term3       - eps delta b                    rad
*    fb          - frame bias matrix              rad
*
*  locals        :
*    convrt      - conversion from arcsec to rad
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2007, 217
* --------------------------------------------------------------------------- */

	void framebias
		(
		char opt,
		double& term1, double& term2, double& term3, std::vector< std::vector<double> > &fb
		)
	{
		fb.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = fb.begin(); it != fb.end(); ++it)
			it->resize(3);
		double convrt;
		convrt = pi / (3600.0 * 180.0);

		// j2000 version referred to iau76/80 theory
		if (opt == 'j')
		{
			term1 = -0.0146 * convrt;
			term2 = -0.016617 * convrt;
			term3 = -0.0068192 * convrt;
			fb[0][0] = 0.99999999999999;
			fb[0][1] = 0.00000007078279;
			fb[0][2] = -0.00000008056149;
			fb[1][0] = -0.00000007078280;
			fb[1][1] = 1.0;
			fb[1][2] = -0.00000003306041;
			fb[2][0] = 0.00000008056149;
			fb[2][1] = 0.00000003306041;
			fb[2][2] = 1.0;
		}

		// fk5 version - catalog origin
		if (opt == 'f')
		{
			term1 = -0.0229 * convrt;
			term2 = 0.0091 * convrt;
			term3 = -0.0199 * convrt;
			fb[0][0] = 0.99999999999999;
			fb[0][1] = 0.00000011102234;
			fb[0][2] = 0.00000004411803;
			fb[1][0] = -0.00000011102233;
			fb[1][1] = 0.99999999999999;
			fb[1][2] = -0.00000009647793;
			fb[2][0] = -0.00000004411804;
			fb[2][1] = 0.00000009647792;
			fb[2][2] = 0.99999999999999;
		}

	}  // procedure framebias




/* ----------------------------------------------------------------------------
*
*                           function itrf_gcrf
*
*  this function transforms a vector between the earth fixed (itrf) frame, and
*    the gcrf mean equator mean equinox. this is the preferrred method to
*    accomplish the new iau 2000 resolutions and uses the eop corrections
*
*  author        : david vallado                  719-573-2600   23 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*    ritrf       - position vector earth fixed    km
*    vitrf       - velocity vector earth fixed    km/s
*    aitrf       - acceleration vector earth fixedkm/s2
*    direct      - direction of transfer          eFrom, eTo
*    iau80rec    - record containing the iau80 constants rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       rad
*    yp          - polar motion coefficient       rad
*    eqeterms    - terms for ast calculation      0,2
*    ddpsi       - delta psi correction to gcrf   rad
*    ddeps       - delta eps correction to gcrf   rad
*    deltapsi    - nutation angle                 rad
*    deltaeps    - nutation angle                 rad
*
*  outputs       :
*    rgcrf       - position vector gcrf            km
*    vgcrf       - velocity vector gcrf            km/s
*    agcrf       - acceleration vector gcrf        km/s2
*    trans       - matrix for pef - gcrf
*
*  locals        :
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    prec        - matrix for mod - gcrf
*    nut         - matrix for tod - mod
*    st          - matrix for pef - tod
*    stdot       - matrix for pef - tod rate
*    pm          - matrix for itrf - pef
*
*  coupling      :
*   precess      - rotation for precession
*   nutation     - rotation for nutation
*   sidereal     - rotation for sidereal time
*   polarm       - rotation for polar motion
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void itrf_gcrf
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		const iau80data &iau80rec,
		double ttt, double jdut1, double lod, double xp,
		double yp, int eqeterms, double ddpsi, double ddeps,
		double deltapsi, double deltaeps,
		std::vector< std::vector<double> > &trans
		)
	{
		trans.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
			it->resize(3);
        // locals
		double psia, wa, epsa, chia,
			trueeps, meaneps, omega, thetasa, omegaearth[3], omgxr[3], omgxomgxr[3],
			rpef[3], vpef[3], apef[3], omgxv[3], tempvec1[3], tempvec[3];
		double  deg2rad;
		//std::vector< std::vector<double> > prec, st, stdot, pm, pmp, stp, nutp, precp, temp;
		std::vector< std::vector<double> > prec = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > nut = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > st = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stdot = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pm = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > precp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > nutp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pmp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(3, std::vector<double>(3));

		deg2rad = pi / 180.0;

		// ---- find matrices
		coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
		coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);
		coordFK5::sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms, st, stdot);
		coordFK5::polarm(xp, yp, ttt, e80, pm);

		// ---- perform transformations
		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth[0] = 0.0;
		omegaearth[1] = 0.0;
		omegaearth[2] = thetasa;

		if (direct == eTo)
		{
			astMath::matvecmult(pm, ritrf, rpef);
			astMath::matmult(prec, nut, temp, 3, 3, 3);
			astMath::matmult(temp, st, trans, 3, 3, 3);
			astMath::matvecmult(trans, rpef, rgcrf);

			astMath::matvecmult(pm, vitrf, vpef);
			astMath::cross(omegaearth, rpef, omgxr);
			astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
			astMath::matvecmult(trans, tempvec1, vgcrf);

			astMath::matvecmult(pm, aitrf, apef);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
			astMath::matvecmult(trans, tempvec1, agcrf);
		}
		else
		{
			astMath::mattrans(pm, pmp, 3, 3);
			astMath::mattrans(st, stp, 3, 3);
			astMath::mattrans(nut, nutp, 3, 3);
			astMath::mattrans(prec, precp, 3, 3);

			astMath::matmult(stp, nutp, temp, 3, 3, 3);
			astMath::matmult(temp, precp, trans, 3, 3, 3);
			astMath::matvecmult(trans, rgcrf, rpef);
			astMath::matvecmult(pmp, rpef, ritrf);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::matvecmult(trans, vgcrf, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::matvecmult(pmp, vpef, vitrf);

			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::matvecmult(trans, agcrf, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
			astMath::matvecmult(pmp, apef, aitrf);
		}
	}  // procedure itrf_gcrf


/* ----------------------------------------------------------------------------
*
*                           function itrf_j2k
*
*  this function transforms a vector between the earth fixed (itrf) frame, and
*    the j2k mean equator mean equinox (iau-76, fk5). note that the delta
*    corrections are not included for this approach.
*
*  author        : david vallado                  719-573-2600   30 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*    ritrf       - position vector earth fixed    km
*    vitrf       - velocity vector earth fixed    km/s
*    aitrf       - acceleration vector earth fixedkm/s2
*    direct      - direction of transfer          eFrom, eTo
*    iau80rec    - record containing the iau80 constants rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       rad
*    yp          - polar motion coefficient       rad
*    eqeterms    - terms for ast calculation      0,2
*
*  outputs       :
*    rj2k        - position vector j2k            km
*    vj2k        - velocity vector j2k            km/s
*    aj2k        - acceleration vector j2k        km/s2
*
*  locals        :
*    deltapsi    - nutation angle                 rad
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    prec        - matrix for mod - j2k
*    nut         - matrix for tod - mod
*    st          - matrix for pef - tod
*    stdot       - matrix for pef - tod rate
*    pm          - matrix for itrf - pef
*
*  coupling      :
*   precess      - rotation for precession
*   nutation     - rotation for nutation
*   sidereal     - rotation for sidereal time
*   polarm       - rotation for polar motion
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void itrf_j2k
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rj2k[3], double vj2k[3], double aj2k[3],
		const iau80data &iau80rec,
		double ttt, double jdut1, double lod,
		double xp, double yp, int eqeterms,
		std::vector< std::vector<double> > &trans
		)
	{
		trans.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
			it->resize(3);
		// locals
		//std::vector< std::vector<double> > prec, nut, st, stdot, pm, temp,pmp, stp, nutp, precp;
		std::vector< std::vector<double> > prec = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > nut = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > st = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stdot = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pm = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > precp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > nutp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pmp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		double psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps,
			omega, thetasa, omegaearth[3], rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3],
			omgxv[3], tempvec1[3], tempvec[3];

		// ---- find matrices
		coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
		coordFK5::nutation(ttt, 0.0, 0.0, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);
		coordFK5::sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms, st, stdot);
		coordFK5::polarm(xp, yp, ttt, e80, pm);

		// ---- perform transformations
		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth[0] = 0.0;
		omegaearth[1] = 0.0;
		omegaearth[2] = thetasa;

		if (direct == eTo)
		{
			astMath::matvecmult(pm, ritrf, rpef);
			astMath::matmult(prec, nut, temp, 3, 3, 3);
			astMath::matmult(temp, st, trans, 3, 3, 3);
			astMath::matvecmult(trans, rpef, rj2k);

			astMath::matvecmult(pm, vitrf, vpef);
			astMath::cross(omegaearth, rpef, omgxr);
			astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
			astMath::matvecmult(trans, tempvec1, vj2k);

			astMath::matvecmult(pm, aitrf, apef);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
			astMath::matvecmult(trans, tempvec1, aj2k);
		}
		else
		{
			astMath::mattrans(pm, pmp, 3, 3);
			astMath::mattrans(st, stp, 3, 3);
			astMath::mattrans(nut, nutp, 3, 3);
			astMath::mattrans(prec, precp, 3, 3);

			astMath::matmult(stp, nutp, temp, 3, 3, 3);
			astMath::matmult(temp, precp, trans, 3, 3, 3);
			astMath::matvecmult(trans, rj2k, rpef);
			astMath::matvecmult(pmp, rpef, ritrf);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::matvecmult(trans, vj2k, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::matvecmult(pmp, vpef, vitrf);

			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::matvecmult(trans, aj2k, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
			astMath::matvecmult(pmp, apef, aitrf);
		}
	} // procedure itrf_j2k

/* ----------------------------------------------------------------------------
*
*                           function itrf_mod
*
*  this function transforms a vector between the earth fixed (itrf) frame, and
*    the gcrf mean equator mean equinox (j2000).
*
*  author        : david vallado                  719-573-2600   22 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*    ritrf       - position vector earth fixed    km
*    vitrf       - velocity vector earth fixed    km/s
*    aitrf       - acceleration vector earth fixedkm/s2
*    direct      - direction of transfer          eFrom, eTo
*    iau80rec    - record containing the iau80 constants rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       rad
*    yp          - polar motion coefficient       rad
*    eqeterms    - terms for ast calculation      0,2
*    ddpsi       - delta psi correction to gcrf   rad
*    ddeps       - delta eps correction to gcrf   rad
*
*  outputs       :
*    rmod        - position vector mod            km
*    vmod        - velocity vector mod            km/s
*    amod        - acceleration vector mod        km/s2
*
*  locals        :
*    deltapsi    - nutation angle                 rad
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    nut         - matrix for tod - mod
*    st          - matrix for pef - tod
*    stdot       - matrix for pef - tod rate
*    pm          - matrix for itrf - pef
*
*  coupling      :
*   nutation     - rotation for nutation
*   sidereal     - rotation for sidereal time
*   polarm       - rotation for polar motion
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void itrf_mod
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rmod[3], double vmod[3], double amod[3],
		const iau80data &iau80rec,
		double ttt, double jdut1, double lod, double xp,
		double yp, int eqeterms, double ddpsi, double ddeps,
		std::vector< std::vector<double> > &trans
		)
	{
		trans.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
			it->resize(3);
		//std::vector< std::vector<double> > nut, st, stdot, pm, temp[3][3], tempmat,	pmp, stp, nutp, precp;
		std::vector< std::vector<double> > nut = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > st = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stdot = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pm = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > precp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > nutp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pmp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > tempmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
		double deltapsi, deltaeps, trueeps, meaneps, omega, thetasa, omegaearth[3],
			rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3],
			omgxv[3], tempvec1[3], tempvec[3];

		// ---- find matrices
		coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);
		coordFK5::sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms, st, stdot);
		coordFK5::polarm(xp, yp, ttt, e80, pm);

		// ---- perform transformations
		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth[0] = 0.0;
		omegaearth[1] = 0.0;
		omegaearth[2] = thetasa;

		if (direct == eTo)
		{
			astMath::matvecmult(pm, ritrf, rpef);
			astMath::matmult(nut, st, tempmat, 3, 3, 3);
			astMath::matvecmult(tempmat, rpef, rmod);

			astMath::matvecmult(pm, vitrf, vpef);
			astMath::cross(omegaearth, rpef, omgxr);
			astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
			astMath::matvecmult(tempmat, tempvec1, vmod);

			astMath::matvecmult(pm, aitrf, apef);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
			astMath::matvecmult(tempmat, tempvec1, amod);
		}
		else
		{
			astMath::mattrans(pm, pmp, 3, 3);
			astMath::mattrans(st, stp, 3, 3);
			astMath::mattrans(nut, nutp, 3, 3);

			astMath::matmult(stp, nutp, tempmat, 3, 3, 3);
			astMath::matvecmult(tempmat, rmod, rpef);
			astMath::matvecmult(pmp, rpef, ritrf);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::matvecmult(tempmat, vmod, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::matvecmult(pmp, vpef, vitrf);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::matvecmult(tempmat, amod, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
			astMath::matvecmult(pmp, apef, aitrf);
		}
	}  // procedure itrf_mod



/* ----------------------------------------------------------------------------
*
*                           function itrf_tod
*
*  this function transforms a vector between the earth fixed (itrf) frame, and
*    the true of date (tod).
*
*  author        : david vallado                  719-573-2600   22 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*    ritrf       - position vector earth fixed    km
*    vitrf       - velocity vector earth fixed    km/s
*    aitrf       - acceleration vector earth fixedkm/s2
*    direct      - direction of transfer          eFrom, eTo
*    iau80rec    - record containing the iau80 constants rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    xp          - polar motion coefficient       rad
*    yp          - polar motion coefficient       rad
*    eqeterms    - terms for ast calculation      0,2
*    ddpsi       - delta psi correction to gcrf   rad
*    ddeps       - delta eps correction to gcrf   rad
*
*  outputs       :
*    rtod        - position vector tod            km
*    vtod        - velocity vector tod            km/s
*    atod        - acceleration vector tod        km/s2
*
*  locals        :
*    deltapsi    - nutation angle                 rad
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    nut         - matrix for tod - mod
*    st          - matrix for pef - tod
*    stdot       - matrix for pef - tod rate
*    pm          - matrix for itrf - pef
*
*  coupling      :
*   nutation     - rotation for nutation
*   sidereal     - rotation for sidereal time
*   polarm       - rotation for polar motion
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void itrf_tod
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rtod[3], double vtod[3], double atod[3],
		const iau80data &iau80rec,
		double ttt, double jdut1, double lod, double xp,
		double yp, int eqeterms, double ddpsi, double ddeps,
		std::vector< std::vector<double> > &trans
		)
	{
		trans.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
			it->resize(3);
		// locals
		//std::vector< std::vector<double> > nut, st, stdot, pm, temp[3][3], tempmat,	pmp, stp, nutp, precp;
		std::vector< std::vector<double> > nut = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > st = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stdot = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pm = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > precp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > nutp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pmp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > tempmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
		double deltapsi, deltaeps, trueeps, meaneps, omega, thetasa, omegaearth[3],
			rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3],
			omgxv[3], tempvec1[3], tempvec[3];

		// ---- find matrices
		coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);
		coordFK5::sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms, st, stdot);
		coordFK5::polarm(xp, yp, ttt, e80, pm);

		// ---- perform transformations
		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth[0] = 0.0;
		omegaearth[1] = 0.0;
		omegaearth[2] = thetasa;

		if (direct == eTo)
		{
			astMath::matvecmult(pm, ritrf, rpef);
			astMath::matvecmult(st, rpef, rtod);

			astMath::matvecmult(pm, vitrf, vpef);
			astMath::cross(omegaearth, rpef, omgxr);
			astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
			astMath::matvecmult(st, tempvec1, vtod);

			astMath::matvecmult(pm, aitrf, apef);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
			astMath::matvecmult(st, tempvec1, atod);
		}
		else
		{
			astMath::mattrans(pm, pmp, 3, 3);
			astMath::mattrans(st, stp, 3, 3);

			astMath::matvecmult(stp, rtod, rpef);
			astMath::matvecmult(pmp, rpef, ritrf);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::matvecmult(stp, vtod, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::matvecmult(pmp, vpef, vitrf);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::matvecmult(stp, atod, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
			astMath::matvecmult(pmp, apef, aitrf);
		}
	}  // procedure itrf_tod


/* ----------------------------------------------------------------------------
*
*                           function itrf_pef
*
*  this function transforms a vector between the earth fixed (itrf) frame, and
*    the psuedo earth fixed (pef).
*
*  author        : david vallado                  719-573-2600   22 nov 2005
*
*  revisions
*
*  inputs          description                    range / units
*    ritrf       - position vector earth fixed    km
*    vitrf       - velocity vector earth fixed    km/s
*    aitrf       - acceleration vector earth fixedkm/s2
*    direct      - direction of transfer          eFrom, eTo
*    ttt         - julian centuries of tt         centuries
*    xp          - polar motion coefficient       rad
*    yp          - polar motion coefficient       rad
*
*  outputs       :
*    rpef        - position vector pef            km
*    vpef        - velocity vector pef            km/s
*    apef        - acceleration vector pef        km/s2
*
*  locals        :
*    pm          - matrix for itrf - pef
*
*  coupling      :
*   polarm       - rotation for polar motion
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void itrf_pef
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rpef[3], double vpef[3], double apef[3],
		double ttt, double xp, double yp,
		std::vector< std::vector<double> > &trans
		)
	{
		trans.resize(3);  // rows
		for (std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
			it->resize(3);
        // locals
		//std::vector< std::vector<double> > pm, temp[3][3], pmp;
		std::vector< std::vector<double> > pm = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pmp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(3, std::vector<double>(3));

		// ---- find matrices
		coordFK5::polarm(xp, yp, ttt, e80, pm);

		// ---- perform transformations
		if (direct == eTo)
		{
			astMath::matvecmult(pm, ritrf, rpef);

			astMath::matvecmult(pm, vitrf, vpef);

			astMath::matvecmult(pm, aitrf, apef);
		}
		else
		{
			astMath::mattrans(pm, pmp, 3, 3);

			astMath::matvecmult(pmp, rpef, ritrf);

			astMath::matvecmult(pmp, vpef, vitrf);

			astMath::matvecmult(pmp, apef, aitrf);
		}
	} // procedure itrf_pef


/* ----------------------------------------------------------------------------
*
*                           function teme_ecef
*
*  this function transforms a vector between the earth fixed(ITRF) frame and the
*  true equator mean equniox frame(teme).the results take into account
*    the effects of sidereal time, and polar motion.
*
*  author        : david vallado                  719 - 573 - 2600   30 oct 2017
*
*  revisions
*
*  inputs          description                    range / units
*    rteme        - position vector teme            km
*    vteme        - velocity vector teme            km / s
*    ateme        - acceleration vector teme        km / s2
*    direct       - direction of transfer           eFrom, 'TOO '
*    ttt          - julian centuries of tt          centuries
*    jdut1        - julian date of ut1              days from 4713 bc
*    lod          - excess length of day            sec
*    xp           - polar motion coefficient        rad
*    yp           - polar motion coefficient        rad
*    eqeterms     - use extra two terms(kinematic) after 1997  0, 2
*
*  outputs       :
*    recef        - position vector earth fixed     km
*    vecef        - velocity vector earth fixed     km / s
*    aecef        - acceleration vector earth fixed km / s2
*
*  locals :
*    st           - matrix for pef - tod
*    pm           - matrix for ecef - pef
*
*  coupling :
*   gstime        - greenwich mean sidereal time    rad
*   polarm        - rotation for polar motion       pef - ecef
*
*  references :
*    vallado       2013, 231 - 233
* ----------------------------------------------------------------------------*/

	void teme_ecef
		(
		double rteme[3], double vteme[3], double ateme[3],
		edirection direct,
		double recef[3], double vecef[3], double aecef[3],
		double ttt, double jdut1, double lod, double xp, double yp, int eqeterms
		)
	{
		double deg2rad, omega, gmstg, thetasa;
		//trans.resize(3);  // rows
		//for (std::vector< std::vector<double> >::iterator it = trans.begin(); it != trans.end(); ++it)
		//	it->resize(3);
		std::vector< std::vector<double> > st = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > stdot = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > temp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > tempmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > pm, pmp, stp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		double omegaearth[3], rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3],
			omgxv[3], tempvec1[3], tempvec[3], gmst;

		deg2rad = pi / 180.0;

		// find omeage from nutation theory
		omega = 125.04452222 + (-6962890.5390 *ttt + 7.455 *ttt*ttt + 0.008 *ttt*ttt*ttt) / 3600.0;
		omega = fmod(omega, 360.0) * deg2rad;

		// ------------------------find gmst--------------------------
		gmst = astTime::gstime(jdut1);

		// teme does not include the geometric terms here
		// after 1997, kinematic terms apply
		if ((jdut1 > 2450449.5) && (eqeterms > 0))
		{
			gmstg = gmst
				+ 0.00264*pi / (3600.0 * 180.0)*sin(omega)
				+ 0.000063*pi / (3600.0 * 180.0)*sin(2.0 *omega);
		}
		else
			gmstg = gmst;
		gmstg = fmod(gmstg, 2.0*pi);

		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth[0] = 0.0;
		omegaearth[1] = 0.0;
		omegaearth[2] = thetasa;
		
		st[0][0] = cos(gmstg);
		st[0][1] = -sin(gmstg);
		st[0][2] = 0.0;
		st[1][0] = sin(gmstg);
		st[1][1] = cos(gmstg);
		st[1][2] = 0.0;
		st[2][0] = 0.0;
		st[2][1] = 0.0;
		st[2][2] = 1.0;

		// compute sidereal time rate matrix
		stdot[0][0] = -omegaearth[2] * sin(gmstg);
		stdot[0][1] = -omegaearth[2] * cos(gmstg);
		stdot[0][2] = 0.0;
		stdot[1][0] = omegaearth[2] * cos(gmstg);
		stdot[1][1] = -omegaearth[2] * sin(gmstg);
		stdot[1][2] = 0.0;
		stdot[2][0] = 0.0;
		stdot[2][1] = 0.0;
		stdot[2][2] = 0.0;

		coordFK5::polarm(xp, yp, ttt, e80, pm);

		if (direct == eTo)
		{
			astMath::mattrans(pm, pmp, 3, 3);
			astMath::mattrans(st, stp, 3, 3);

			astMath::matvecmult(stp, rteme, rpef);
			astMath::matvecmult(pmp, rpef, recef);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::matvecmult(stp, vteme, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::matvecmult(pmp, vpef, vecef);

			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::matvecmult(stp, ateme, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
			astMath::matvecmult(pmp, apef, aecef);
		}
		else
		{
			astMath::matvecmult(pm, recef, rpef);
			astMath::matvecmult(st, rpef, rteme);

			astMath::matvecmult(pm, vecef, vpef);
			astMath::cross(omegaearth, rpef, omgxr);
			astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
			astMath::matvecmult(st, tempvec1, vteme);

			astMath::matvecmult(pm, aecef, apef);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
			astMath::matvecmult(st, tempvec1, ateme);
		}

		//fprintf(1, 'st gmst %11.8f ast %11.8f ome  %11.8f \n', gmst * 180 / pi, ast * 180 / pi, omegaearth * 180 / pi);

	}  // procedure teme_ecef


/* ----------------------------------------------------------------------------
*
*                           function pef_gcrf
*
*  this function transforms a vector between the pseudo earth fixed frame (pef),
*    and the mean equator mean equinox (j2000) frame.
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  revisions
*    vallado     - add terms for ast calculation                 30 sep 2002
*    vallado     - consolidate with iau 2000                     14 feb 2005
*    vallado     - conversion to c++                             21 feb 2005
*
*  inputs          description                    range / units
*    rpef        - position pseudo earth fixed    km
*    vpef        - velocity pseudo earth fixed    km/s
*    apef        - acceleration pseudo earth fixedkm/s2
*    direct      - direction of transfer          eFrom, 'TOO '
*    iau80rec    - record containing the iau80 constants rad
*    ttt         - julian centuries of tt         centuries
*    jdut1       - julian date of ut1             days from 4713 bc
*    lod         - excess length of day           sec
*    terms       - number of terms for ast calculation 0,2
*
*  outputs       :
*    rgcrf        - position vector gcrf            km
*    vgcrf        - velocity vector gcrf            km/s
*    agcrf        - acceleration vector gcrf        km/s2
*
*  locals        :
*    prec        - matrix for gcrf - mod
*    deltapsi    - nutation angle                 rad
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    nut         - matrix for mod - tod
*    st          - matrix for tod - pef
*    stdot       - matrix for tod - pef rate
*
*  coupling      :
*   precess      - rotation for precession        mod - gcrf
*   nutation     - rotation for nutation          tod - mod
*   sidereal     - rotation for sidereal time     pef - tod
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void pef_gcrf
		(
		double rpef[3], double vpef[3], double apef[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		const iau80data &iau80rec,
		double ttt, double jdut1, double lod, int eqeterms,
		double ddpsi, double ddeps
		)
	{
		std::vector< std::vector<double> > prec, nut, st, stdot, temp, tempmat, stp, nutp, precp = 
			std::vector< std::vector<double> >(3, std::vector<double>(3));
		double psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps,
			omega, thetasa, omegaearth[3], omgxr[3], omgxomgxr[3],
			omgxv[3], tempvec1[3], tempvec[3];

		coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
		coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);
		coordFK5::sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms, st, stdot);

		thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		omegaearth[0] = 0.0;
		omegaearth[1] = 0.0;
		omegaearth[2] = thetasa;

		if (direct == eTo)
		{
			astMath::matmult(prec, nut, temp, 3, 3, 3);
			astMath::matmult(temp, st, tempmat, 3, 3, 3);
			astMath::matvecmult(tempmat, rpef, rgcrf);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
			astMath::matvecmult(tempmat, tempvec1, vgcrf);

			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
			astMath::matvecmult(tempmat, tempvec1, agcrf);
		}
		else
		{
			astMath::mattrans(st, stp, 3, 3);
			astMath::mattrans(nut, nutp, 3, 3);
			astMath::mattrans(prec, precp, 3, 3);

			astMath::matmult(stp, nutp, temp, 3, 3, 3);
			astMath::matmult(temp, precp, tempmat, 3, 3, 3);
			astMath::matvecmult(tempmat, rgcrf, rpef);

			astMath::cross(omegaearth, rpef, omgxr);
			astMath::matvecmult(tempmat, vgcrf, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);

			astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
			astMath::cross(omegaearth, vpef, omgxv);
			astMath::cross(omegaearth, omgxr, omgxomgxr);
			astMath::matvecmult(tempmat, agcrf, tempvec1);
			astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
			astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
		}
	}  // procedure pef_gcrf

/* -----------------------------------------------------------------------------
*
*                           function tod_gcrf
*
*  this function transforms a vector between the true equator true equinox frame
*    of date (tod), and the mean equator mean equinox (j2000) frame.
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  revisions
*    vallado     - consolidate with iau 2000                     14 feb 2005
*    vallado     - conversion to c++                             21 feb 2005
*
*  inputs          description                    range / units
*    rtod        - position vector of date
*                    true equator, true equinox   km
*    vtod        - velocity vector of date
*                    true equator, true equinox   km/s
*    atod        - acceleration vector of date
*                    true equator, true equinox   km/s2
*    direct      - direction of transfer          eFrom, 'TOO '
*    iau80rec    - record containing the iau80 constants rad
*    ttt         - julian centuries of tt         centuries
*
*  outputs       :
*    rgcrf        - position vector gcrf            km
*    vgcrf        - velocity vector gcrf            km/s
*    agcrf        - acceleration vector gcrf        km/s2
*
*  locals        :
*    deltapsi    - nutation angle                 rad
*    trueeps     - true obliquity of the ecliptic rad
*    meaneps     - mean obliquity of the ecliptic rad
*    omega       -                                rad
*    nut         - matrix for mod - tod
*
*  coupling      :
*   precess      - rotation for precession        mod - gcrf
*   nutation     - rotation for nutation          tod - mod
*
*  references    :
*    vallado       2013, 228, Alg 24
* ----------------------------------------------------------------------------*/

	void tod_gcrf
		(
		double rtod[3], double vtod[3], double atod[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		const iau80data &iau80rec,
		double ttt, double ddpsi, double ddeps
		)
	{
		std::vector< std::vector<double> > prec, nut, nutp, precp = std::vector< std::vector<double> >(3, std::vector<double>(3));
		std::vector< std::vector<double> > tempmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
		double psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps, omega;

		coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
		coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);

		if (direct == eTo)
		{
			astMath::matmult(prec, nut, tempmat, 3, 3, 3);
			astMath::matvecmult(tempmat, rtod, rgcrf);
			astMath::matvecmult(tempmat, vtod, vgcrf);
			astMath::matvecmult(tempmat, atod, agcrf);
		}
		else
		{
			astMath::mattrans(nut, nutp, 3, 3);
			astMath::mattrans(prec, precp, 3, 3);
			astMath::matmult(nutp, precp, tempmat, 3, 3, 3);

			astMath::matvecmult(tempmat, rgcrf, rtod);
			astMath::matvecmult(tempmat, vgcrf, vtod);
			astMath::matvecmult(tempmat, agcrf, atod);
		}
	}  // procedure tod_gcrf

/* -----------------------------------------------------------------------------
*
*                           function mod_gcrf
*
*  this function transforms a vector between the mean equator mean equinox of
*    date (mod) and the mean equator mean equinox (j2000) frame.
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  revisions
*    vallado     - consolidate with iau 2000                     14 feb 2005
*    vallado     - conversion to c++                             21 feb 2005
*
*  inputs          description                    range / units
*    rmod        - position vector of date
*                    mean equator, mean equinox   km
*    vmod        - velocity vector of date
*                    mean equator, mean equinox   km/s
*    amod        - acceleration vector of date
*                    mean equator, mean equinox   km/s2
*    direct      - direction of transfer          eFrom, 'TOO '
*    ttt         - julian centuries of tt         centuries
*
*  outputs       :
*    rgcrf        - position vector gcrf            km
*    vgcrf        - velocity vector gcrf            km/s
*    agcrf        - acceleration vector gcrf        km/s2
*
*  locals        :
*    none.
*
*  coupling      :
*   precess      - rotation for precession        mod - gcrf
*
*  references    :
*    vallado       2013, 228, Alg 24
* --------------------------------------------------------------------------- */

	void mod_gcrf
		(
		double rmod[3], double vmod[3], double amod[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		double ttt
		)
	{
		std::vector< std::vector<double> > prec, precp;
		double psia, wa, epsa, chia;

		coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);

		if (direct == eTo)
		{
			astMath::matvecmult(prec, rmod, rgcrf);
			astMath::matvecmult(prec, vmod, vgcrf);
			astMath::matvecmult(prec, amod, agcrf);
		}
		else
		{
			astMath::mattrans(prec, precp, 3, 3);

			astMath::matvecmult(precp, rgcrf, rmod);
			astMath::matvecmult(precp, vgcrf, vmod);
			astMath::matvecmult(precp, agcrf, amod);
		}
	}  // procedure mod_gcrf



/* ----------------------------------------------------------------------------
*
*                           function teme_eci
*
*  this function transforms a vector from the true equator mean equinox system,
*  (teme) to the mean equator mean equinox (j2000) system.
*
*  author        : david vallado                  719 - 573 - 2600   30 oct 2017
*
*  inputs        description                    range / units
*    rteme       - position vector of date
*                  true equator, mean equinox      km
*    vteme       - velocity vector of date
*                  true equator, mean equinox      km / s
*    ateme       - acceleration vector of date
*                  true equator, mean equinox      km / s2
*    ttt         - julian centuries of tt          centuries
*    ddpsi       - delta psi correction to gcrf    rad
*    ddeps       - delta eps correction to gcrf    rad
*
*  outputs       :
*    reci        - position vector eci             km
*    veci        - velocity vector eci             km / s
*    aeci        - acceleration vector eci         km / s2
*
*  locals :
*    prec        - matrix for eci - mod
*    nutteme     - matrix for mod - teme - an approximation for nutation
*    eqeg        - rotation for equation of equinoxes(geometric terms only)
*    tm          - combined matrix for teme2eci
*
*  coupling      :
*   precess      - rotation for precession         eci - mod
*   nutation     - rotation for nutation           eci - tod
*
*  references :
*    vallado       2013, 231 - 233
* ----------------------------------------------------------------------------*/

	void teme_eci
		(
		double rteme[3], double vteme[3], double ateme[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		const iau80data &iau80rec,
		double ttt, double ddpsi, double ddeps
		)
	{
		std::vector< std::vector<double> > prec, nut, temp, tempmat, nutp, precp, eqep;
		std::vector< std::vector<double> > eqe = std::vector< std::vector<double> >(3, std::vector<double>(3));
		double psia, wa, epsa, chia, deltapsi, deltaeps, trueeps, meaneps, omega, eqeg;

		coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
		coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);

		// ------------------------find eqeg----------------------
		// rotate teme through just geometric terms
		eqeg = deltapsi* cos(meaneps);
		eqeg = fmod(eqeg, twopi);

		eqe[0][0] = cos(eqeg);
		eqe[0][1] = sin(eqeg);
		eqe[0][2] = 0.0;
		eqe[1][0] = -sin(eqeg);
		eqe[1][1] = cos(eqeg);
		eqe[1][2] = 0.0;
		eqe[2][0] = 0.0;
		eqe[2][1] = 0.0;
		eqe[2][2] = 1.0;

		if (direct == eTo)
		{
			astMath::mattrans(eqe, eqep, 3, 3);

			astMath::matmult(nut, eqep, temp, 3, 3, 3);
			astMath::matmult(prec, temp, tempmat, 3, 3, 3);

			astMath::matvecmult(tempmat, rteme, rgcrf);
			astMath::matvecmult(tempmat, vteme, vgcrf);
			astMath::matvecmult(tempmat, ateme, agcrf);
		}
		else
		{
			astMath::mattrans(nut, nutp, 3, 3);
			astMath::mattrans(prec, precp, 3, 3);

			astMath::matmult(nutp, precp, temp, 3, 3, 3);
			astMath::matmult(eqe, temp, tempmat, 3, 3, 3);

			astMath::matvecmult(tempmat, rgcrf, rteme);
			astMath::matvecmult(tempmat, vgcrf, vteme);
			astMath::matvecmult(tempmat, agcrf, ateme);
		}

	}  //  teme_eci


}  // namespace

