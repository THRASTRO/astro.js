/*############################################################################*/
// AstroJS Coordinate Module (c) 2016 by Marcel Greter
// https://www.github.com/mgreter/astrojs/LICENSE
/*############################################################################*/
// https://en.wikipedia.org/wiki/Cartesian_coordinate_system (xyz)
// https://en.wikipedia.org/wiki/Spherical_coordinate_system (lbr)
// https://en.wikipedia.org/wiki/Cylindrical_coordinate_system (zlp)
/*############################################################################*/
'use strict';

(function (exports) {

	/*############################################################################*/
	// private converter and helper functions
	/*############################################################################*/

	// getter functions for calculated coordinates
	function cyl() { return { p: this.p, l: this.l, z: this.z }; }
	function cart() { return { x: this.x, y: this.y, z: this.z }; }
	function sph() { return { l: this.l, b: this.b, i: this.i, r: this.r }; }

	/*############################################################################*/
	// These functions act as auto-loaders. They will convert the
	// coordinates once and then replace the getter function.
	/*############################################################################*/

	// Cartesian to Spherical
	function cart2sph()
	{
		this.r = Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
		this.i = Math.acos(this.z / this.r), this.b = PI/2 - this.i;
		this.l = Math.atan2(this.y, this.x);
		return (this.sph = sph).call(this);
	}

	// Spherical to Cartesian
	function sph2cart()
	{
		var rsb = this.r * Math.sin(this.i);
		this.x = rsb * Math.cos(this.l);
		this.y = rsb * Math.sin(this.l);
		this.z = this.r * Math.cos(this.i);
		return (this.cart = cart).call(this);
	}

	// Cartesian to Cylindrical
	function cart2cyl()
	{
		this.p = Math.sqrt(this.x*this.x + this.y*this.y);
		this.l = Math.atan2(this.y, this.x);
		return (this.cyl = cyl).call(this);
	}

	// Cylindrical to Cartesian
	function cyl2cart()
	{
		this.x = this.p * Math.cos(this.l);
		this.y = this.p * Math.sin(this.l);
		return (this.cart = cart).call(this);
	}

	// Cylindrical to Spherical
	function cyl2sph()
	{
		this.i = Math.atan2(this.p, this.z), this.b = PI/2 - this.i;
		this.r = Math.sqrt(this.p*this.p + this.z*this.z);
		return (this.sph = sph).call(this);
	}

	// Spherical to Cylindrical
	function sph2cyl()
	{
		this.p = this.r * Math.cos(this.b);
		this.z = this.r * Math.sin(this.b);
		return (this.cyl = cyl).call(this);
	}

	/*############################################################################*/
	// Class Constructor for Coordinate
	/*############################################################################*/

	/******************************************************************************/
	// Coordinate supports conversion between Cartesian, Spherical
	// and Cylindrical Coordinates. Use either for the constructor.
	/******************************************************************************/
	function Coord(state)
	{
		// from Cartesian coordinates
		if ('x' in state && 'y' in state && 'z' in state) {
			this.x = state.x, this.y = state.y, this.z = state.z;
			this.cart = cart, this.sph = cart2sph, this.cyl = cart2cyl;
		}
		// from Cylindrical coordinates
		else if ('p' in state && 'l' in state && 'z' in state) {
			this.p = state.p, this.l = state.l, this.z = state.z;
			this.cart = cyl2cart, this.sph = cyl2sph, this.cyl = cyl;
		}
		// from Spherical coordinates
		else if ('l' in state && 'r' in state) {
			// store given parameters, then ...
			this.l = state.l, this.r = state.r;
			// ... inclination or elevation
			if ('i' in state || 'b' in state) {
				this.cart = sph2cart, this.sph = sph, this.cyl = sph2cyl;
				if ('b' in state) this.b = state.b, this.i = PI/2 - this.b;
				else if ('i' in state) this.i = state.i, this.b = PI/2 - this.i;
			}
		}
		// ToDo: error out on invalid input?
		// Silently ignores ra/dec without dist!
	}

	/*############################################################################*/
	// conversion between various formats
	/*############################################################################*/
	// Equatorial - declination (d) and hour angle (h)
	//            - declination (d) and right ascension (a)
	// Ecliptic - ecliptic latitude (ß) and ecliptic longitude (?)
	// Galactic - galactic latitude (b) and galactic longitude (l)
	// Observer - altitude (a) or elevation and azimuth (A)
	// Note: VSOP87 coordinates are heliocentric ecliptic
	/*############################################################################*/
	// ecliptic to equatorial (compensate axial tilt)
	// equatorial to galactic (compensate precession)
	/*############################################################################*/
	// Mean equinox of date: Equator rotated by precession to its position
	// at "date", but free from the small periodic oscillations of nutation
	// True equinox of date: This is the actual intersection of the two
	// planes at any particular moment, with all motions accounted for.
	/*############################################################################*/

	// ecliptic to equatorial coordinates
	// pass obliquity of ecliptic (axial tilt)
	Coord.prototype.ecl2equ = function ecl2equ(tilt)
	{
		// convert via Cartesian
		var cart = this.cart(),
		    sin_e = Math.sin(tilt),
		    cos_e = Math.cos(tilt);
		// return new coordinate
		return new Coord({
			x: cart.x,
			y: cart.y * cos_e - cart.z * sin_e,
			z: cart.z * cos_e + cart.y * sin_e,
		});
	}

	// equatorial to ecliptic coordinates
	// pass obliquity of ecliptic (axial tilt)
	Coord.prototype.equ2ecl = function equ2ecl(tilt)
	{
		// convert via Cartesian
		var cart = this.cart(),
		    sin_e = Math.sin(tilt),
		    cos_e = Math.cos(tilt);
		// return new coordinate
		return new Coord({
			x: cart.x,
			y: cart.y * cos_e + cart.z * sin_e,
			z: cart.z * cos_e - cart.y * sin_e,
		});
	}

	var posangle = 32.932 * DEG2RAD;
	var pole_ra = 192.859508 * DEG2RAD;
	var pole_dec = 27.128336 * DEG2RAD;

	// North galactic pole (B1950)
	// var pole_ra = 192.25 * DEG2RAD;
	// var pole_dec = 27.4 * DEG2RAD;
	// var posangle = 33.0 * DEG2RAD;

	// galactic to equatorial coordinates
	// equatorial coordinate has only l/b
	Coord.prototype.gal2equ = function gal2equ()
	{
		// pre-calculations
		var sph = this.sph(),
		    l = sph.l, b = sph.b,
		    sin_b = Math.sin(b),
		    cos_b = Math.cos(b),
		    sin_pdc = Math.sin(pole_dec),
		    cos_pdc = Math.cos(pole_dec),
		    sin_pos = Math.sin(l - posangle),
		    cos_pos = Math.cos(l - posangle),
		    sincos_bp = cos_b * sin_pos;
		// ToDo: compensate for precision if epoch differs
		// to galactic coordinates after precession compensated
		var B = Math.asin( sin_b * sin_pdc + sincos_bp * cos_pdc );
		var L = Math.atan2(cos_b * cos_pos, sin_b * cos_pdc - sincos_bp * sin_pdc);
		// return object in valid range
		// ToDo: should `b` be half turn?
		return { l: CYCLE(L + pole_ra), b: TURN(B) };
	}

	// equatorial to galactic coordinates
	// gallactic coordinate has only l/b
	Coord.prototype.equ2gal = function equ2gal()
	{
		// pre-calculations
		var sph = this.sph(),
		    l = sph.l, b = sph.b,
		    sin_b = Math.sin(b),
		    cos_b = Math.cos(b),
		    tan_b = sin_b / cos_b,
		    sin_pdc = Math.sin(pole_dec),
		    cos_pdc = Math.cos(pole_dec),
		    sin_pos = Math.sin(pole_ra - l),
		    cos_pos = Math.cos(pole_ra - l),
		    sincos_bp = cos_b * sin_pos;
		// ToDo: compensate for precision if epoch differs
		// to equatorial coordinates after precession compensated
		var B = Math.asin( sin_b * sin_pdc + cos_b * cos_pdc * cos_pos );
		var L = Math.atan2( sin_pos, cos_pos * sin_pdc - tan_b * cos_pdc);
		// return object in valid range
		// ToDo: should `b` be half turn?
		return { l: CYCLE(posangle - L - PI/2), b: TURN(B) };
	}

	/*############################################################################*/
	// Precession to convert from the international terrestrial reference
	// system (ITRF) to the international celestial reference frame (ICRF)
	/*############################################################################*/

	// precession factors for J2000 frame
	var PFA = 2306.2181*DEG2RAD, PFB = 1.39656*DEG2RAD, PFC = 0.000139*DEG2RAD,
	    PFD = 0.30188*DEG2RAD, PFE = 0.000344*DEG2RAD, PFF = 0.017998*DEG2RAD,
	    PFG = 1.09468*DEG2RAD, PFH = 0.000066*DEG2RAD, PFI = 0.018203*DEG2RAD,
	    PFJ = 2004.3109*DEG2RAD, PFK = 0.85330*DEG2RAD, PFL = 0.000217*DEG2RAD,
	    PFM = 0.42665*DEG2RAD, PFN = 0.000217*DEG2RAD, PFP = 0.041833*DEG2RAD;

	Coord.prototype.precess = function precess(JD, epoch)
	{
		// Change original ra and dec to radians
		var mean_ra = this.l, mean_dec = this.b;
		// calc t, zeta, eta and theta for J2000.0 Equ 20.3
		var eta, zeta, theta, T, T2,
		    t = (JD - epoch) / 360000.0,
		    t2 = t * t, t3 = t2 *t;
		// optimize for epoch == 0
		if (t) {
			T = epoch / 360000.0, T2 = T * T;
			eta = (PFA + PFB * T - PFC * T2) * t
			     + (PFG + PFH * T) * t2 + PFI * t3;
			zeta = (PFA + PFB * T - PFC * T2) * t
			     + (PFD - PFE * T) * t2 + PFF * t3;
			theta = (PFJ - PFK * T - PFN * T2) * t
				    - (PFM + PFN * T) * t2 - PFP * t3;
		}
		else {
			eta = PFA * t + PFG * t2 + PFP * t3;
			zeta = PFA * t + PFD * t2 + PFF * t3;
			theta = PFJ * t - PFN * t2 - PFP * t3;
		}

		// calc A,B,C equ 20.4
		var sin_md = Math.sin(mean_dec), cos_md = Math.cos(mean_dec),
		    mean_term = cos_md * Math.cos(mean_ra + zeta),
		    A = cos_md * Math.sin(mean_ra + zeta),
		    B = Math.cos(theta) * mean_term - Math.sin(theta) * Math.sin(mean_dec),
		    C = Math.sin (theta) * mean_term + Math.cos(theta) * Math.sin(mean_dec);

		// calculate right ascension
		var ra = Math.atan2(A, B) + eta, dec;

		// check for object near celestial pole
		if (mean_dec > (0.4 * PI) || mean_dec < (-0.4 * PI)) {
			// close to pole (better precision)
			dec = Math.acos(Math.sqrt(A * A + B * B));
			// check 0 <= acos() <= PI
			if (mean_dec < 0.) { dec *= -1;  }
		}
		// not close to pole
		else { dec = Math.asin(C); }

		// return object in valid range
		// ToDo: should `b` be half turn?
		var coord = { l: CYCLE(ra), b: TURN(dec) };
		if ('r' in this) coord.r = this.r;
		return coord;
	}

	/*############################################################################*/
	// END OF AstroJS Coordinate Module
	/*############################################################################*/

	exports.Coord = Coord;

})(this);
