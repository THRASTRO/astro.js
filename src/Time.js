/*############################################################################*/
// AstroJS Time Module (c) 2016 by Marcel Greter
// https://www.github.com/mgreter/astrojs/LICENSE
/*############################################################################*/
// Implementation is mostly one to one from libnova
// http://libnova.sourceforge.net/group__sidereal.html
// http://libnova.sourceforge.net/group__dynamical.html
// http://libnova.sourceforge.net/group__heliocentric.html
// http://libnova.sourceforge.net/group__nutation.html
/*############################################################################*/
'use strict';

(function (exports)
{

	/*############################################################################*/
	// constant factors for calculations
	/*############################################################################*/

	// dynamical time offset in seconds for
	// every second year from 1620 to 1992
	var delta_t = [
		124.0, 115.0, 106.0, 98.0, 91.0,
		85.0, 79.0, 74.0, 70.0, 65.0,
		62.0, 58.0, 55.0, 53.0, 50.0,
		48.0, 46.0, 44.0, 42.0, 40.0,
		37.0, 35.0, 33.0, 31.0, 28.0,
		26.0, 24.0, 22.0, 20.0, 18.0,
		16.0, 14.0, 13.0, 12.0, 11.0,
		10.0, 9.0, 9.0, 9.0, 9.0,
		9.0, 9.0, 9.0, 9.0, 10.0,
		10.0, 10.0, 10.0, 10.0, 11.0,
		11.0, 11.0, 11.0, 11.0, 11.0,
		11.0, 12.0, 12.0, 12.0, 12.0,
		12.0, 12.0, 13.0, 13.0, 13.0,
		13.0, 14.0, 14.0, 14.0, 15.0,
		15.0, 15.0, 15.0, 16.0, 16.0,
		16.0, 16.0, 16.0, 17.0, 17.0,
		17.0, 17.0, 17.0, 17.0, 17.0,
		17.0, 16.0, 16.0, 15.0, 14.0,
		13.7, 13.1, 12.7, 12.5, 12.5,
		12.5, 12.5, 12.5, 12.5, 12.3,
		12.0, 11.4, 10.6, 9.6, 8.6,
		7.5, 6.6, 6.0, 5.7, 5.6,
		5.7, 5.9, 6.2, 6.5, 6.8,
		7.1, 7.3, 7.5, 7.7, 7.8,
		7.9, 7.5, 6.4, 5.4, 2.9,
		1.6, -1.0, -2.7, -3.6, -4.7,
		-5.4, -5.2, -5.5, -5.6, -5.8,
		-5.9, -6.2, -6.4, -6.1, -4.7,
		-2.7, 0.0, 2.6, 5.4, 7.7,
		10.5, 13.4, 16.0, 18.2, 20.2,
		21.2, 22.4, 23.5, 23.9, 24.3,
		24.0, 23.9, 23.9, 23.7, 24.0,
		24.3, 25.3, 26.2, 27.3, 28.2,
		29.1, 30.0, 30.7, 31.4, 32.2,
		33.1, 34.0, 35.0, 36.5, 38.3,
		40.2, 42.2, 44.5, 46.5, 48.5,
		50.5, 52.2, 53.8, 54.9, 55.8,
		56.9, 58.3
	];

	// arguments and coefficients taken from table 21A on page 133
	var nut_args = [
		[ 0.0,  0.0,  0.0,  0.0,  1.0], [-2.0,  0.0,  0.0,  2.0,  2.0], [ 0.0,  0.0,  0.0,  2.0,  2.0],
		[ 0.0,  0.0,  0.0,  0.0,  2.0], [ 0.0,  1.0,  0.0,  0.0,  0.0], [ 0.0,  0.0,  1.0,  0.0,  0.0],
		[-2.0,  1.0,  0.0,  2.0,  2.0], [ 0.0,  0.0,  0.0,  2.0,  1.0], [ 0.0,  0.0,  1.0,  2.0,  2.0],
		[-2.0, -1.0,  0.0,  2.0,  2.0], [-2.0,  0.0,  1.0,  0.0,  0.0], [-2.0,  0.0,  0.0,  2.0,  1.0],
		[ 0.0,  0.0, -1.0,  2.0,  2.0], [ 2.0,  0.0,  0.0,  0.0,  0.0], [ 0.0,  0.0,  1.0,  0.0,  1.0],
		[ 2.0,  0.0, -1.0,  2.0,  2.0], [ 0.0,  0.0, -1.0,  0.0,  1.0], [ 0.0,  0.0,  1.0,  2.0,  1.0],
		[-2.0,  0.0,  2.0,  0.0,  0.0], [ 0.0,  0.0, -2.0,  2.0,  1.0], [ 2.0,  0.0,  0.0,  2.0,  2.0],
		[ 0.0,  0.0,  2.0,  2.0,  2.0], [ 0.0,  0.0,  2.0,  0.0,  0.0], [-2.0,  0.0,  1.0,  2.0,  2.0],
		[ 0.0,  0.0,  0.0,  2.0,  0.0], [-2.0,  0.0,  0.0,  2.0,  0.0], [ 0.0,  0.0, -1.0,  2.0,  1.0],
		[ 0.0,  2.0,  0.0,  0.0,  0.0], [ 2.0,  0.0, -1.0,  0.0,  1.0], [-2.0,  2.0,  0.0,  2.0,  2.0],
		[ 0.0,  1.0,  0.0,  0.0,  1.0], [-2.0,  0.0,  1.0,  0.0,  1.0], [ 0.0, -1.0,  0.0,  0.0,  1.0],
		[ 0.0,  0.0,  2.0, -2.0,  0.0], [ 2.0,  0.0, -1.0,  2.0,  1.0], [ 2.0,  0.0,  1.0,  2.0,  2.0],
		[ 0.0,  1.0,  0.0,  2.0,  2.0], [-2.0,  1.0,  1.0,  0.0,  0.0], [ 0.0, -1.0,  0.0,  2.0,  2.0],
		[ 2.0,  0.0,  0.0,  2.0,  1.0], [ 2.0,  0.0,  1.0,  0.0,  0.0], [-2.0,  0.0,  2.0,  2.0,  2.0],
		[-2.0,  0.0,  1.0,  2.0,  1.0], [ 2.0,  0.0, -2.0,  0.0,  1.0], [ 2.0,  0.0,  0.0,  0.0,  1.0],
		[ 0.0, -1.0,  1.0,  0.0,  0.0], [-2.0, -1.0,  0.0,  2.0,  1.0], [-2.0,  0.0,  0.0,  0.0,  1.0],
		[ 0.0,  0.0,  2.0,  2.0,  1.0], [-2.0,  0.0,  2.0,  0.0,  1.0], [-2.0,  1.0,  0.0,  2.0,  1.0],
		[ 0.0,  0.0,  1.0, -2.0,  0.0], [-1.0,  0.0,  1.0,  0.0,  0.0], [-2.0,  1.0,  0.0,  0.0,  0.0],
		[ 1.0,  0.0,  0.0,  0.0,  0.0], [ 0.0,  0.0,  1.0,  2.0,  0.0], [ 0.0,  0.0, -2.0,  2.0,  2.0],
		[-1.0, -1.0,  1.0,  0.0,  0.0], [ 0.0,  1.0,  1.0,  0.0,  0.0], [ 0.0, -1.0,  1.0,  2.0,  2.0],
		[ 2.0, -1.0, -1.0,  2.0,  2.0], [ 0.0,  0.0,  3.0,  2.0,  2.0], [ 2.0, -1.0,  0.0,  2.0,  2.0]
	];
	var nut_coeffs = [
		[-171996.0, -174.2, 92025.0,8.9], [-13187.0, -1.6, 5736.0, -3.1], [-2274.0, -0.2, 977.0, -0.5],
		[2062.0,   0.2,-895.0,  0.5], [1426.0,  -3.4,  54.0,  -0.1], [ 712.0,   0.1,  -7.0,   0.0],
		[-517.0,   1.2, 224.0, -0.6], [-386.0,  -0.4, 200.0,   0.0], [-301.0,   0.0, 129.0,  -0.1],
		[ 217.0, - 0.5, -95.0,  0.3], [-158.0,   0.0,   0.0,   0.0], [ 129.0,   0.1, -70.0,   0.0],
		[ 123.0,   0.0, -53.0,  0.0], [  63.0,   0.0,   0.0,   0.0], [  63.0,   0.1, -33.0,   0.0],
		[ -59.0,   0.0,  26.0,  0.0], [ -58.0,  -0.1,  32.0,   0.0], [ -51.0,   0.0,  27.0,   0.0],
		[  48.0,   0.0,   0.0,  0.0], [  46.0,   0.0, -24.0,   0.0], [ -38.0,   0.0,  16.0,   0.0],
		[ -31.0,   0.0,  13.0,  0.0], [  29.0,   0.0,   0.0,   0.0], [  29.0,   0.0, -12.0,   0.0],
		[  26.0,   0.0,   0.0,  0.0], [ -22.0,   0.0,   0.0,   0.0], [  21.0,   0.0, -10.0,   0.0],
		[  17.0, - 0.1,   0.0,  0.0], [  16.0,   0.0,  -8.0,   0.0], [ -16.0,   0.1,   7.0,   0.0],
		[ -15.0,   0.0,   9.0,  0.0], [ -13.0,   0.0,   7.0,   0.0], [ -12.0,   0.0,   6.0,   0.0],
		[  11.0,   0.0,   0.0,  0.0], [ -10.0,   0.0,   5.0,   0.0], [  -8.0,   0.0,   3.0,   0.0],
		[   7.0,   0.0,  -3.0,  0.0], [  -7.0,   0.0,   0.0,   0.0], [  -7.0,   0.0,   3.0,   0.0],
		[  -7.0,   0.0,   3.0,  0.0], [   6.0,   0.0,   0.0,   0.0], [   6.0,   0.0,  -3.0,   0.0],
		[   6.0,   0.0,  -3.0,  0.0], [  -6.0,   0.0,   3.0,   0.0], [  -6.0,   0.0,   3.0,   0.0],
		[   5.0,   0.0,   0.0,  0.0], [  -5.0,   0.0,   3.0,   0.0], [  -5.0,   0.0,   3.0,   0.0],
		[  -5.0,   0.0,   3.0,  0.0], [   4.0,   0.0,   0.0,   0.0], [   4.0,   0.0,   0.0,   0.0],
		[   4.0,   0.0,   0.0,  0.0], [  -4.0,   0.0,   0.0,   0.0], [  -4.0,   0.0,   0.0,   0.0],
		[  -4.0,   0.0,   0.0,  0.0], [   3.0,   0.0,   0.0,   0.0], [  -3.0,   0.0,   0.0,   0.0],
		[  -3.0,   0.0,   0.0,  0.0], [  -3.0,   0.0,   0.0,   0.0], [  -3.0,   0.0,   0.0,   0.0],
		[  -3.0,   0.0,   0.0,  0.0], [  -3.0,   0.0,   0.0,   0.0], [  -3.0,   0.0,   0.0,   0.0]
	];

	/*############################################################################*/
	// return nutation parameters on given julian days
	/*############################################################################*/

	var FD0 = 297.85036 * DEG2RAD, FD1 = 445267.111480 * DEG2RAD, FD2 = - 0.0019142 * DEG2RAD, FD3 = DEG2RAD / 189474.0,
	    FM0 = 357.52772 * DEG2RAD, FM1 = 35999.050340 * DEG2RAD, FM2 = - 0.0086972 * DEG2RAD, FM3 = - DEG2RAD / 300000.0,
	    FMM0 = 134.96298 * DEG2RAD, FMM1 = 477198.867398 * DEG2RAD, FMM2 = 0.0086972 * DEG2RAD, FMM3 = DEG2RAD / 56250.0,
	    FF0 = 93.2719100 * DEG2RAD, FF1 = 483202.017538 * DEG2RAD, FF2 = - 0.0036825 * DEG2RAD, FF3 = DEG2RAD / 327270.0,
	    FO0 = 125.04452 * DEG2RAD, FO1 = - 1934.136261 * DEG2RAD, FO2 = 0.0020708 * DEG2RAD, FO3 = DEG2RAD / 450000.0;

	function getNutation(JD)
	{

		// should we bother recalculating nutation (caching)
		// if (fabs(JD - c_JD) > EPOCH_THRESHOLD) {
		// set the new epoch //
		// c_JD = JD;

		// get julian ephemeris day
		var JDE = JD2JDE(JD);

		// calc T and dependants
		var T = JD2J2K(JDE) / 100,
		    T2 = T * T, T3 = T2 * T;

		// calculate D,M,M',F and Omega (convert to radians)
		var D = FD0 + FD1 * T + FD2 * T2 + FD3 * T3;
		var M = FM0 + FM1 * T + FM2 * T2 + FM3 * T3;
		var MM = FMM0 + FMM1 * T + FMM2 * T2 + FMM3 * T3;
		var F = FF0 + FF1 * T + FF2 * T2 + FF3 * T3;
		var O = FO0 + FO1 * T + FO2 * T2 + FO3 * T3;

		// init variables for coeff sums
		// get accumulated in for loop
		var c_longitude = 0, c_obliquity = 0;
		var coeff_sine, coeff_cos, arg;

		// calc sum of terms in table 21A
		for (var i = 0; i < nut_coeffs.length; i++) {
			// calc coefficients of sine and cosine
			coeff_sine = (nut_coeffs[i][0] +
				(nut_coeffs[i][1] * T));
			coeff_cos = (nut_coeffs[i][2] +
				(nut_coeffs[i][3] * T));
			// calculate main argument
			arg = nut_args[i][0] * D
			    + nut_args[i][1] * M
			    + nut_args[i][2] * MM
			    + nut_args[i][3] * F
			    + nut_args[i][4] * O;
			// sum up oscillating elements
			c_longitude += coeff_sine * Math.sin(arg);
			c_obliquity += coeff_cos * Math.cos(arg);
		}

		// c_longitude /= 100000.0;
		// c_obliquity /= 100000.0;

		// calculate mean ecliptic - Meeus 2nd edition, eq. 22.2
		var c_ecliptic = 23.0 + 26.0 / 60.0 + 21.448 / 3600.0
		               - 46.81500 / 3600.0 * T
		               - 0.000590 / 3600.0 * T2
		               + 0.001813 / 3600.0 * T3;

		// change to radiants
		c_ecliptic *= DEG2RAD;
		// use special conversion factor
		// change time (hours?) to to arcsecs
		c_longitude /= 3600 * 100 * 100 / DEG2RAD;
		c_obliquity /= 3600 * 100 * 100 / DEG2RAD;

		// Uncomment this if function should return
		// true obliquity rather than mean obliquity
		/* c_ecliptic += c_obliquity; */

		// } // caching

		// return results
		return {
			longitude: TURN(c_longitude), // verified
			obliquity: TURN(c_obliquity), // verified
			ecliptic: CYCLE(c_ecliptic) // unverified
		}

	}
	// EO nutation

	// Stephenson and Houlden for years prior to 948 A.D.
	function getDynamicalTimeDiffSH1(JD)
	{
		// number of centuries from 948
		var E = (JD - 2067314.5) / 36525.0;
		return 1830.0 - 405.0 * E + 46.5 * E * E;
	}

	// Stephenson and Houlden for years between 948 A.D. and 1600 A.D.
	function getDynamicalTimeDiffSH2(JD)
	{
		/* number of centuries from 1850 */
		var t = (JD - 2396758.5) / 36525.0;
		return 22.5 * t * t;
	}

	// Table 9.a pg 72 for years 1600..1992
	// uses interpolation formula 3.3 on pg 25
	function getDynamicalTimeDiffTable(JD)
	{
		// get no days since 1620 and divide by 2 years
		var i = Math.floor((JD - 2312752.5) / 730.5);
		// get the base interpolation factor in the table
		if (i > (delta_t.length - 2)) i = delta_t.length - 2;
		// calculate a, b, c, n
		var a = delta_t[i+1] - delta_t[i], b = delta_t[i+2] - delta_t[i+1],
		    c = a - b, n = ((JD - (2312752.5 + (730.5 * i))) / 730.5);
		return delta_t[i + 1] + n / 2 * (a + b + n * c);
	}

	// get the dynamical time diff in the near past / future 1992 .. 2010
	/// uses interpolation formula 3.3 on pg 25
	function getDynamicalTimeDiffNear(JD)
	{
		// TD for 1990, 2000, 2010
		var delta_T = [ 56.86, 63.83, 70.0 ];
		// calculate TD by interpolating value
		var a = delta_T[1] - delta_T[0],
		    b = delta_T[2] - delta_T[1],
		    c = b - a;
		// get number of days since 2000 and divide by 10 years
		var n = (JD - 2451544.5) / 3652.5;
		return delta_T[1] + (n / 2) * (a + b + n * c);
	}

	// uses equation 9.1 pg 73 to calc JDE for other JD values
	function getDynamicalTimeDiffElse(JD)
	{
		var a = (JD - 2382148.0);
		return -15.0 + a * a / 41048480.0;
	}

	// Calculates the dynamical time (TD) difference in seconds
	// (delta T) from universal time. Equation 9.1 on pg 73.
	// check when JD is, and use corresponding formula
	function getDynamicalTimeDiff(JD)
	{
		// check for date < 948 A.D.
		if (JD < 2067314.5) {
			// Stephenson and Houlden
			return getDynamicalTimeDiffSH1(JD);
		}
		// check for date 948..1600 A.D.
		else if (JD >= 2067314.5 && JD < 2305447.5) {
			// Stephenson and Houlden
			return getDynamicalTimeDiffSH2(JD);
		}
		// check for value in table 1620..1992
		else if (JD >= 2312752.5 && JD < 2448622.5) {
			// interpolation of table
			return getDynamicalTimeDiffTable(JD);
		}
		// check for near future 1992..2010
		else if (JD >= 2448622.5 && JD <= 2455197.5) {
			// use value interpolation
			return getDynamicalTimeDiffNear(JD);
		}
		// other time period (outside)
		return getDynamicalTimeDiffElse(JD);
	}

	// Calculates the Julian Ephemeris Day(JDE)
	function JD2JDE(JD)
	{
		// from the given Julian Day (JD)
		return JD + getDynamicalTimeDiff(JD) / JD2SEC;
	}

	// Calculates the Julian Day(JDE)
	function JDE2JD(JDE)
	{
		// from the given Julian Ephemeris Day (JD)
		return JDE - getDynamicalTimeDiff(JDE) / JD2SEC;
	}

	/*############################################################################*/
	// Sidereal time calculations
	/*############################################################################*/

	// factors for sidereal time offset
	var FST1 = 280.46061837 * DEG2RAD,
	    FST2 = 360.98564736629 * DEG2RAD,
	    FST3 = 0.000387933 * DEG2RAD,
	    FST4 = - DEG2RAD / 38710000.0,
	    // ratio longitude nutation
	    FLON = 1; // 15.0 * DEG2HMS;

	// returns rad (you want HMS?)
	function getMeanSiderealTime(JD)
	{
		var T = JD2J2K(JD) / 100;
		// calculate mean angle (in radians)
		var sidereal = FST1 + FST2 *(JD - 2451545.0)
			           + FST3 * T * T + FST4 * T * T * T;
		// clamp to valid range
		return CYCLE(sidereal);
	}

	// returns rad (you want HMS?)
	function getApparentSiderealTime(JD)
	{
		// get nutation parameters
		var nutation = getNutation(JD);
		// get the mean sidereal time and add
		return getMeanSiderealTime(JD) +
			// add corrections for nutation in longitude
			// and for the true obliquity o the ecliptic
			nutation.longitude / FLON *
			Math.cos(nutation.obliquity);
	}

	/*############################################################################*/
	// END OF AstroJS Time Module
	/*############################################################################*/

	// date converter
	exports.JD2JDE = JD2JDE;
	exports.JDE2JD = JDE2JD;

	// sidereal time calculations
	exports.getMeanSiderealTime = getMeanSiderealTime;
	exports.getApparentSiderealTime = getApparentSiderealTime;

	// calculate earth nutations
	exports.getNutation = getNutation;

})(this);
