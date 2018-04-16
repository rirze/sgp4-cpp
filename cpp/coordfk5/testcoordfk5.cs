// test the coordfk5 rountines

	// ------------------------------------------  test coords  ------------------------------------------
	double recef[3] = { -1033.4793830,  7901.2952754,  6380.3565958 };
	double vecef[3] = { -3.225636520, -2.872451450, 5.531924446 };
	double aecef[3] = { 0.001, 0.002, 0.003 };
	std::vector< std::vector<double> > st00, stdot00, st, stdot;
	std::vector< std::vector<double> > nut, nut00;
	std::vector< std::vector<double> > prec;
	std::vector< std::vector<double> > pm, trans;
	double receft[3], veceft[3], aeceft[3],
		rmod[3], vmod[3], amod[3], rmod20[3], vmod20[3], amod20[3], 
		rtod[3], vtod[3], atod[3], rtod20[3], vtod20[3], atod20[3], 
	    reci20[3], veci20[3], aeci20[3], reci00[3], veci00[3], aeci00[3],
		rpef[3], vpef[3], apef[3], rpef20[3], vpef20[3], apef20[3], rpef00[3], vpef00[3], apef00[3],
		reci[3], veci[3], aeci[3], recit[3], vecit[3], aecit[3], recit20[3], vecit20[3], aecit20[3],
		rteme[3], vteme[3], ateme[3];
    iau80data iau80rec;
	char EopLoc[85] = "../nut80.dat";  //D:/Codes/LIBRARY/CPP/TestSGP4DC/
	coordFK5::iau80in(iau80rec, EopLoc);

	double conv, dat, xp, yp, ddeps, ddy, dut1, lod, ddx, ddpsi, ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdttf, tcg, tdb, ttdb, jdtdb, jdtdbf, tcb, timezone, sec;
	double deltapsi, deltaeps, trueeps, meaneps, omega, psia, wa, epsa, chia, deltapsi00, deltaeps00, trueeps00, meaneps00, omega00;
	int year, mon, day, hr, minute, eqeterms;
	conv = pi / (180.0*3600.0);
	eqeterms = 2;

	year = 2004;
	mon = 4;
	day = 6;
	hr = 7;
	minute = 51;
	sec = 28.386009;
	dut1 = -0.4399619;  // sec
	dat = 32;         // sec
	xp = -0.140682 * conv;  // " to rad
	yp = 0.333309 * conv;
	lod = 0.0015563;
	ddpsi = -0.052195 * conv;  // " to rad
	ddeps = -0.003875 * conv;
	ddx = -0.000205 * conv;  // " to rad
	ddy = -0.000136 * conv;
	timezone = 0.0;
	
	astTime::convtime(year, mon, day, hr, minute, sec, timezone, dut1, dat, ut1, tut1, jdut1, jdut1f, utc, tai, tt, ttt, jdtt, jdttf, tcg, tdb, ttdb, jdtdb, jdtdbf, tcb);
	printf("ut1 %8.6f tut1 %16.12f jdut1 %18.11f ", ut1, tut1, jdut1 + jdut1f);
	astTime::hms_sec(hr, minute, sec, eFrom, ut1);

	printf("hms %3i %3i %8.6f \n", hr, minute, sec);
	printf("utc %8.6f ", utc);
	astTime::hms_sec(hr, minute, sec, eFrom, utc);
	printf("hms %3i %3i %8.6f \n", hr, minute, sec);
	printf("tai %8.6f ", tai);
	astTime::hms_sec(hr, minute, sec, eFrom, tai);
	printf("hms %3i %3i %8.6f \n", hr, minute, sec);
	printf("tt  %8.6f ttt  %16.12f jdtt  %18.11f ", tt, ttt, jdtt + jdttf);
	astTime::hms_sec(hr, minute, sec, eFrom, tt);
	printf("hms %3i %3i %8.6f \n", hr, minute, sec);
	printf("tdb %8.6f ttdb %16.12f jdtdb %18.11f\n", tdb, ttdb, jdtdb + jdtdbf);

	// --------precess - transformation matrix for precession
	printf("precession matrix \n");
	coordFK5::precess(ttt, e80, psia, wa, epsa, chia, prec);
	printf("%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n\n", prec[0][0], prec[0][1], prec[0][2], prec[1][0], prec[1][1], prec[1][2], prec[2][0], prec[2][1], prec[2][2]);
	
	// --------nutation - transformation matrix for nutation
	printf("nutation matrix \n");
	coordFK5::nutation(ttt, 0.0, 0.0, iau80rec, deltapsi00, deltaeps00, trueeps00, meaneps00, omega00, nut00);
	printf("%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n\n", nut00[0][0], nut00[0][1], nut00[0][2], nut00[1][0], nut00[1][1], nut00[1][2], nut00[2][0], nut00[2][1], nut00[2][2]);
	coordFK5::nutation(ttt, ddpsi, ddeps, iau80rec, deltapsi, deltaeps, trueeps, meaneps, omega, nut);
	printf("%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n\n", nut[0][0], nut[0][1], nut[0][2], nut[1][0], nut[1][1], nut[1][2], nut[2][0], nut[2][1], nut[2][2]);

	// --------sidereal - transformation matrix for sidereal time
	printf("sidereal time matrix \n");
	coordFK5::sidereal(jdut1 + jdut1f, deltapsi, meaneps, omega, lod, 0, st00, stdot00);
	printf("%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n\n", st00[0][0], st00[0][1], st00[0][2], st00[1][0], st00[1][1], st00[1][2], st00[2][0], st00[2][1], st00[2][2]);
	coordFK5::sidereal(jdut1 + jdut1f, deltapsi, meaneps, omega, lod, eqeterms, st, stdot);
	printf("%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n\n", st[0][0], st[0][1], st[0][2], st[1][0], st[1][1], st[1][2], st[2][0], st[2][1], st[2][2]);
	// note you could also have two calculations with the wo corr nutation values
	
	// --------polarm - transformation matrix for polar motion
	printf("polar motion matrix \n");
	coordFK5::polarm(xp, yp, ttt, e80, pm);
	printf("%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n%16.11f %16.11f %16.11f\n\n", pm[0][0], pm[0][1], pm[0][2], pm[1][0], pm[1][1], pm[1][2], pm[2][0], pm[2][1], pm[2][2]);

	//------------------ -
	// sample runs for book
//	[recigg, vecigg, aecig] = iau00f2iS(recef, vecef, aecef, ttt, jdut1 + jdut1frac, lod, xp, yp, "c", ddx, ddy);
//	printf("GCRF          IAU-2006 CIO %14.7f %14.7f %14.7f", recigg);
//	printf(" v %14.9f %14.9f %14.9f", vecigg);
//	printf(" a %14.9f %14.9f %14.9f\n", aecig);

	printf("\n\n ============== convert various coordinate systems from ecef =================== \n");
	printf("ITRF          IAU-76/FK5   %14.7f %14.7f %14.7f", recef[0], recef[1], recef[2]);
	printf(" v %14.9f %14.9f %14.9f", vecef[0], vecef[1], vecef[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aecef[0], aecef[1], aecef[2]);

	// --------pef transformations
	coordFK5::itrf_pef(recef, vecef, aecef, eTo, rpef, vpef, apef, ttt, xp, yp, trans);
	printf("PEF                        %14.7f %14.7f %14.7f", rpef[0], rpef[1], rpef[2]);
	printf(" v %14.9f %14.9f %14.9f", vpef[0], vpef[1], vpef[2]);
	printf(" a %14.9f %14.9f %14.9f\n", apef[0], apef[1], apef[2]);

	// --------tod transformations
	coordFK5::itrf_tod(recef, vecef, aecef, eTo, rtod, vtod, atod, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, ddpsi, ddeps, trans);
	printf("TOD 2 w corr  IAU-76/FK5   %14.7f %14.7f %14.7f", rtod[0], rtod[1], rtod[2]);
	printf(" v %14.9f %14.9f %14.9f", vtod[0], vtod[1], vtod[2]);
	printf(" a %14.9f %14.9f %14.9f\n", atod[0], atod[1], atod[2]);

	coordFK5::itrf_tod(recef, vecef, aecef, eTo, rtod20, vtod20, atod20, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, 0.0, 0.0, trans);
	printf("TOD 2 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f", rtod20[0], rtod20[1], rtod20[2]);
	printf(" v %14.9f %14.9f %14.9f", vtod20[0], vtod20[1], vtod20[2]);
	printf(" a %14.9f %14.9f %14.9f\n", atod20[0], atod20[1], atod20[2]);

	// --------mod transformations
	coordFK5::itrf_mod(recef, vecef, aecef, eTo, rmod, vmod, amod, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, ddpsi, ddeps, trans);
	printf("MOD 2 w corr  IAU-76/FK5   %14.7f %14.7f %14.7f", rmod[0], rmod[1], rmod[2]);
	printf(" v %14.9f %14.9f %14.9f", vmod[0], vmod[1], vmod[2]);
	printf(" a %14.9f %14.9f %14.9f\n", amod[0], amod[1], amod[2]);

	coordFK5::itrf_mod(recef, vecef, aecef, eTo, rmod20, vmod20, amod20, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, 0.0, 0.0, trans);
	printf("MOD 2 wo corr IAU-76/FK5   %14.7f %14.7f %14.7f", rmod20[0], rmod20[1], rmod20[2]);
	printf(" v %14.9f %14.9f %14.9f", vmod20[0], vmod20[1], vmod20[2]);
	printf(" a %14.9f %14.9f %14.9f\n", amod20[0], amod20[1], amod20[2]);

	// --------gcrf transformations
	coordFK5::itrf_gcrf(recef, vecef, aecef, eTo, reci, veci, aeci, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, ddpsi, ddeps, deltapsi, deltaeps, trans);
	printf("GCRF 2 w corr IAU-76/FK5   %14.7f %14.7f %14.7f", reci[0], reci[1], reci[2]);
	printf(" v %14.9f %14.9f %14.9f", veci[0], veci[1], veci[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeci[0], aeci[1], aeci[2]);

	// --------teme transformations
	coordFK5::teme_ecef(rteme, vteme, ateme, eFrom, recef, vecef, aecef, ttt, jdut1 + jdut1f, lod, xp, yp, eqeterms);
	printf("TEME                       %14.7f %14.7f %14.7f", rteme[0], rteme[1], rteme[2]);
	printf(" v %14.9f %14.9f %14.9f", vteme[0], vteme[1], vteme[2]);
	printf(" a %14.9f %14.9f %14.9f\n", ateme[0], ateme[1], ateme[2]);


	// ----------------------- ind to gcrf transformations ---------------------------
	// --------pef transformations
	coordFK5::pef_gcrf(rpef, vpef, apef, eTo, recit, vecit, aecit, iau80rec, ttt, jdut1 + jdut1f, lod, eqeterms, ddpsi, ddeps);
	printf("pef-eci                    %14.7f %14.7f %14.7f", recit[0], recit[1], recit[2]);
	printf(" v %14.9f %14.9f %14.9f", vecit[0], vecit[1], vecit[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aecit[0], aecit[1], aecit[2]);

	// --------tod transformations
	coordFK5::tod_gcrf(rtod, vtod, atod, eTo, recit, vecit, aecit, iau80rec, ttt, ddpsi, ddeps);
	printf("tod-eci                    %14.7f %14.7f %14.7f", recit[0], recit[1], recit[2]);
	printf(" v %14.9f %14.9f %14.9f", vecit[0], vecit[1], vecit[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aecit[0], aecit[1], aecit[2]);

	coordFK5::tod_gcrf(rtod20, vtod20, atod20, eTo, recit20, vecit20, aecit20, iau80rec, ttt, 0.0, 0.0);
	printf("tod20-j2keci                %14.7f %14.7f %14.7f", recit20[0], recit20[1], recit20[2]);
	printf(" v %14.9f %14.9f %14.9f", vecit20[0], vecit20[1], vecit20[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aecit20[0], aecit20[1], aecit20[2]);

	// --------mod transformations
	coordFK5::mod_gcrf(rmod, vmod, amod, eTo, recit, vecit, aecit, ttt);
	printf("mod-eci                  %14.7f %14.7f %14.7f", recit[0], recit[1], recit[2]);
	printf(" v %14.9f %14.9f %14.9f", vecit[0], vecit[1], vecit[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aecit[0], aecit[1], aecit[2]);

	// --------teme transformations
	coordFK5::teme_eci(rteme, vteme, ateme, eTo, recit, vecit, aecit, iau80rec, ttt, 0.0, 0.0);
	printf("teme-j2keci                 %14.7f %14.7f %14.7f", recit[0], recit[1], recit[2]);
	printf(" v %14.9f %14.9f %14.9f", vecit[0], vecit[1], vecit[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aecit[0], aecit[1], aecit[2]);

	// ----------------------- reverse transformations ---------------------------
	// --------pef transformations
	coordFK5::itrf_pef(receft, veceft, aeceft, eFrom, rpef, vpef, apef, ttt, xp, yp, trans);
	printf("pef - ecef                 %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	// --------tod transformations
	coordFK5::itrf_tod(receft, veceft, aeceft, eFrom, rtod, vtod, atod, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, ddpsi, ddeps, trans);
	printf("tod - ecef                 %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	coordFK5::itrf_tod(receft, veceft, aeceft, eFrom, rtod20, vtod20, atod20, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, 0.0, 0.0, trans);
	printf("tod20 - ecef               %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	// --------mod transformations
	coordFK5::itrf_mod(receft, veceft, aeceft, eFrom, rmod, vmod, amod, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, ddpsi, ddeps, trans);
	printf("mod - ecef                 %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	coordFK5::itrf_mod(receft, veceft, aeceft, eFrom, rmod20, vmod20, amod20, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, 0.0, 0.0, trans);
	printf("mod20 - ecef               %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	// --------gcrf transformations
	coordFK5::itrf_gcrf(receft, veceft, aeceft, eFrom, reci, veci, aeci, iau80rec, ttt, jdut1 + jdut1f, lod, xp, yp, 2, ddpsi, ddeps, deltapsi, deltaeps, trans);
	printf("eci - ecef                 %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	// --------teme transformations
	coordFK5::teme_ecef(rteme, vteme, ateme, eTo, receft, veceft, aeceft, ttt, jdut1 + jdut1f, lod, xp, yp, eqeterms);
	printf("teme-ecef                  %14.7f %14.7f %14.7f", receft[0], receft[1], receft[2]);
	printf(" v %14.9f %14.9f %14.9f", veceft[0], veceft[1], veceft[2]);
	printf(" a %14.9f %14.9f %14.9f\n", aeceft[0], aeceft[1], aeceft[2]);

	// reverse eci transformations
	// --------pef transformations
	coordFK5::pef_gcrf(rpef, vpef, apef, eFrom, reci, veci, aeci, iau80rec, ttt, jdut1 + jdut1f, lod, eqeterms, ddpsi, ddeps);
	printf("eci - pef                %14.7f %14.7f %14.7f", rpef[0], rpef[1], rpef[2]);
	printf(" v %14.9f %14.9f %14.9f", vpef[0], vpef[1], vpef[2]);
	printf(" a %14.9f %14.9f %14.9f\n", apef[0], apef[1], apef[2]);

	// --------tod transformations
	coordFK5::tod_gcrf(rtod, vtod, atod, eFrom, reci, veci, aeci, iau80rec, ttt, ddpsi, ddeps);
	printf("eci - tod                %14.7f %14.7f %14.7f", rtod[0], rtod[1], rtod[2]);
	printf(" v %14.9f %14.9f %14.9f", vtod[0], vtod[1], vtod[2]);
	printf(" a %14.9f %14.9f %14.9f\n", atod[0], atod[1], atod[2]);

	coordFK5::tod_gcrf(rtod20, vtod20, atod20, eFrom, reci, veci, aeci, iau80rec, ttt, 0.0, 0.0);
	printf("eci - tod20              %14.7f %14.7f %14.7f", rtod20[0], rtod20[1], rtod20[2]);
	printf(" v %14.9f %14.9f %14.9f", vtod20[0], vtod20[1], vtod20[2]);
	printf(" a %14.9f %14.9f %14.9f\n", atod20[0], atod20[1], atod20[2]);

	// --------mod transformations
	coordFK5::mod_gcrf(rmod, vmod, amod, eFrom, reci, veci, aeci, ttt);
	printf("eci - mod                %14.7f %14.7f %14.7f", rmod[0], rmod[1], rmod[2]);
	printf(" v %14.9f %14.9f %14.9f", vmod[0], vmod[1], vmod[2]);
	printf(" a %14.9f %14.9f %14.9f\n", amod[0], amod[1], amod[2]);

	coordFK5::mod_gcrf(rmod20, vmod20, amod20, eFrom, reci, veci, aeci, ttt);
	printf("eci - mod20              %14.7f %14.7f %14.7f", rmod20[0], rmod20[1], rmod20[2]);
	printf(" v %14.9f %14.9f %14.9f", vmod20[0], vmod20[1], vmod20[2]);
	printf(" a %14.9f %14.9f %14.9f\n", amod20[0], amod20[1], amod20[2]);

	// --------teme transformations
	printf("ecit20               %14.7f %14.7f %14.7f", recit20[0], recit20[1], recit20[2]);
	printf(" v %14.9f %14.9f %14.9f\n\n", vecit20[0], vecit20[1], vecit20[2]);

	coordFK5::teme_eci(rteme, vteme, ateme, eFrom, recit20, vecit20, aecit20, iau80rec, ttt, 0.0, 0.0);
	printf("j2keci - teme               %14.7f %14.7f %14.7f", rteme[0], rteme[1], rteme[2]);
	printf(" v %14.9f %14.9f %14.9f", vteme[0], vteme[1], vteme[2]);
	printf(" a %14.9f %14.9f %14.9f\n", ateme[0], ateme[1], ateme[2]);




	//tm = fk4;
	//r1950 = tm" * reci20;
	//v1950 = tm" * veci20;
	//printf("eci - FK4                  %14.7f %14.7f %14.7f", r1950);
	//printf(" v %14.9f %14.9f %14.9f \n", v1950);

	//reci20i = tm * r1950;
	//veci20i = tm * v1950;
	//printf("FK4 - eci                  %14.7f %14.7f %14.7f", reci20i);
	//printf(" v %14.9f %14.9f %14.9f \n", veci20i);
