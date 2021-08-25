/*############################################################################*/
// AstroJS Orbit Module (c) 2016 by Marcel Greter
// https://www.github.com/mgreter/astrojs/LICENSE
/*############################################################################*/
// http://www.lns.cornell.edu/~seb/celestia/orbital-parameters.html
// Similar module: https://github.com/jordanstephens/kepler.js
/*############################################################################*/
'use strict';

(function (exports)
{

	// Alias for math functions
	// Better compressibility
	// Saves about 0.4kb
	var abs = Math.abs;
	var pow = Math.pow;
	var sin = Math.sin;
	var cos = Math.cos;
	var sqrt = Math.sqrt;
	var cbrt = Math.cbrt;
	var atan2 = Math.atan2;
	var asin = Math.asin;
	var acos = Math.acos;

	/******************************************************************************/
	// Orbits can be created from orbital elements (6 independent parameters)
	// or from state vectors (position and velocity). In both cases we need 6
	// arguments to fully define an orbit. Most related parameters will be
	// calculated on construction. But certain parameters as the true anomaly
	// are a bit more complex to calculate from orbital elements, if only the
	// mean anomaly is given as initial parameters. They are lazy calculated
	// on demand, either by using the getter or invoking updateElements(true).
	/******************************************************************************/
	// To get state vectors for different times we just add a linear factor
	// to the mean anomaly (related to mean motion and period). To get the
	// actual state vectors from the mean anomaly we need to compute the
	// solution for Kepler's Equation via Newton-Raphson solver loop.
	/******************************************************************************/
	// Keep in mind that all angular parameters are in rads and that all
	// other units are linked to each other. This means that you can change
	// the units by adjusting the gravitational parameter. We use the solar
	// gravitational parameter by default with the units: au^3/(solm*day^2).
	// Normally in physics this parameter is measured in m^3/(kg*s^2). All
	// distance and time units must therefore be in the same units as `G`.
	/******************************************************************************/
	// ToDo: finish and add tests for circular and hyperbola orbits.
	/******************************************************************************/

	function Orbit(arg)
	{

		var orbit = this;

		// internal use
		if (arg.empty) {
			return orbit;
		}

		// the gravitational parameter is very important and must
		// match all others in terms of the units. In physics the
		// units are by default m^3/(kg*s^2), but for astronomical
		// uses we most often want to use au^3/(solm*day^2).
		orbit._G = arg.G || Orbit.GMP.sun,
			orbit._t = arg.epoch || arg.t || 0;
		orbit.translate = arg.translate || null;

		// create state vectors from input parameters
		if ('x' in arg && 'y' in arg && 'z' in arg) {
			orbit._r = new Vector3(arg.x, arg.y, arg.z);
		}
		else if ('rx' in arg && 'ry' in arg && 'rz' in arg) {
			orbit._r = new Vector3(arg.rx, arg.ry, arg.rz);
		}
		if ('X' in arg && 'Y' in arg && 'Z' in arg) {
			orbit._v = new Vector3(arg.X, arg.Y, arg.Z);
		}
		else if ('vx' in arg && 'vy' in arg && 'vz' in arg) {
			orbit._v = new Vector3(arg.vx, arg.vy, arg.vz);
		}

		if ('r' in arg && 'v' in arg) {
			var r = arg.r, v = arg.v;
			// create new vectors (want a clone anyway)
			orbit._r = new Vector3(r.x, r.y, r.z);
			orbit._v = new Vector3(v.x, v.y, v.z);
		}

		// set status for state vectors
		orbit.vectors = !!(orbit._r && orbit._v);

		// get other parameters from arg
		// unknown parameters are ignored
		var vars = 'aeiMnPLcCwWOkhqpG'.split('');
		for (var i in vars) {
			if (arg.hasOwnProperty(vars[i])) {
				orbit['_' + vars[i]] = arg[vars[i]];
			}
		}

		// orbit may has name
		if ('name' in arg) {
			orbit.name = arg.name;
		}

		// resolve if not lazy
		if (!orbit.lazy) {
			orbit.updateElements();
			orbit.requireElements();
		}

	}
	// EO Orbit ctor

	// use mass system
	// has julian days
	Orbit.GMP = GMJD;

	// optimized for minifier
	var Klass = Orbit.prototype;

	/******************************************************************************/
	// clone the existing orbital
	/******************************************************************************/

	// copy existing orbit and return clone
	Klass.clone = function clone(dt)
	{
		var orbit = this;
		// create an empty object (internal use)
		var clone = new Orbit({ empty: true });
		// avoid unnecessary cloning of objects
		if (dt) { delete orbit._v; delete orbit._r; }
		// process all properties
		for (var key in orbit) {
			if (orbit.hasOwnProperty(key)) {
				// invoke nested clone functions (vectors)
				if (orbit[key] && typeof orbit[key].clone == 'function')
				{ clone[key] = orbit[key].clone(); }
				else { clone[key] = orbit[key]; }
			}
		}
		// optionally call update
		if (dt) clone.update(dt);
		// return copy
		return clone;
	}

	/******************************************************************************/
	// check current object for valid resolved state
	/******************************************************************************/
	Klass.checkElements = function checkElements(full)
	{
		var orbit = this;
		// most basic check for initial pre-calculations
		if (!orbit.elements) throw ('Orbitals not calculated');
		// do some basic checks for needed orbital elements
		if (!('_i' in orbit)) throw ('Orbit is missing inclination (i)');
		if (!('_e' in orbit)) throw ('Orbit is missing eccentricity (e)');
		if (!('_a' in orbit)) throw ('Orbit is missing semi-major axis (a)');
		if (!('_b' in orbit)) throw ('Orbit is missing semi-minor axis (b)');
		if (!('_l' in orbit)) throw ('Orbit is missing semi-latus rectum (l)');
		if (!('_c' in orbit)) throw ('Orbit is missing periapsis (c)');
		if (!('_C' in orbit)) throw ('Orbit is missing apoapsis (C)');
		if (!('_L' in orbit)) throw ('Orbit is missing mean longitude (L)');
		if (!('_O' in orbit)) throw ('Orbit is missing right ascending node (O)');
		if (!('_T' in orbit)) throw ('Orbit is missing time of periapsis (w)');
		if (!('_t' in orbit)) throw ('Orbit is missing time of reference epoch (t)');
		if (!('_w' in orbit)) throw ('Orbit is missing argument of periapsis (w)');
		if (!('_W' in orbit)) throw ('Orbit is missing longitude of the periapsis (W)');
		if (!('_n' in orbit)) throw ('Orbit is missing mean motion (n)');
		if (!('_M' in orbit)) throw ('Orbit is missing mean anomaly (M)');
		if (!('_P' in orbit)) throw ('Orbit is missing orbital period (P)');
		if (!('_A' in orbit)) throw ('Orbit is missing angular momentum (A)');
		if (!('_c' in orbit)) throw ('Orbit is missing pericenter (c)');
		if (!('_C' in orbit)) throw ('Orbit is missing apocenter (C)');
		// check if eccentricity is inside valid boundaries (0-1)
		if (orbit._e < 0) throw ('Negative eccentricity is invalid');
		if (orbit._e > 1) throw ('Eccentricity must not be hyperbolic (> 1)');
		if (orbit._e == 1) throw ('Eccentricity must not be parabolic (== 1)');
		// state vectors (r and v) are not tested here ...
		if (full) if (!('_m' in orbit)) throw ('Orbit is missing true anomaly (m)');
		if (full) if (!('_B' in orbit)) throw ('Orbit is missing radial velocity (B)');
		if (full) if (!('_E' in orbit)) throw ('Orbit is missing eccentric anomaly (E)');
	}

	/******************************************************************************/
	// ensure certain orbital elements exist: aeiWLMOw
	/******************************************************************************/
	Klass.requireElements = function requireElements(full)
	{
		// invoke resolver (fully)
		this.updateElements(full);
		this.checkElements(full);
		// chain-able
		return this;
	}

	/******************************************************************************/
	// resolve orbital elements as far as possible
	/******************************************************************************/
	Klass.updateElements = function updateElements(full)
	{

		var orbit = this;
		// only resolve elements once
		if (orbit.elements) return orbit;
		// will calculate now
		orbit.elements = true;

		/********************************************************/
		// resolve from state vectors (r, v)
		/********************************************************/

		// check main state vectors
		if ('_r' in orbit && '_v' in orbit) {

			var r = orbit._r, v = orbit._v,
				G = orbit._G, rl = r.length();

			// specific relative angular momentum
			orbit._A3 = r.clone().cross(v);
			orbit._A2 = orbit._A3.lengthSq();
			orbit._A = sqrt(orbit._A2);
			// calculate radial velocity
			orbit._B = r.dot(v) / rl;

			// calculate eccentricity vector
			var e3 = orbit._e3 = r.clone().multiplyScalar(v.lengthSq() - (G / rl))
				.sub(v.clone().multiplyScalar(rl * orbit._B)).multiplyScalar(1 / G);
			// get eccentricity value from vector
			var e2 = e3.lengthSq(), e = orbit._e = sqrt(e2);
			// get inclination (i) via orbital momentum vector
			orbit._i = acos(orbit._A3.z / orbit._A);

			// calculate semilatus rectum (ℓ)
			orbit._l = orbit._A2 / G;
			// and periapsis (c) and apoapsis (C)
			orbit._c = orbit._l / (1 + e);
			orbit._C = orbit._l / (1 - e);
			// and finally semi-major axis (a)
			orbit._a = (orbit._C + orbit._c) / 2;

			// pre-calculate node line
			var nx = - orbit._A3.y,
				ny = + orbit._A3.x;
			var nl = sqrt(nx * nx + ny * ny);

			// calculate ascending node (O)
			var omega = nl == 0 ? 0 : acos(nx / nl);
			orbit._O = ny < 0 ? (TAU - omega) : omega;

			// calculate argument of periapsis (ω)
			var nedot = nx * e3.x +
				ny * e3.y;
			if (nl === 0 || e === 0) { orbit._w = 0; }
			else { orbit._w = acos(nedot / nl / e); }
			if (e3.z < 0) { orbit._w *= -1; }

			// calculate true anomaly
			var u; // argument of latitude
			// case for circular orbit
			// and without inclination
			if (e === 0 && nl === 0) {
				// orbit needs a test case
				u = acos(r.x / rl);
			}
			// circular orbit
			// with inclination
			else if (e === 0) {
				// orbit needs a test case
				var nrdot = nx * r.x + ny * r.y;
				u = acos(nrdot / (nl * rl));
			}
			// elliptic orbit
			else {
				var redot = e3.x * r.x + e3.y * r.y + e3.z * r.z;
				u = acos(redot / e / rl);
			}
			// bring into correct range via simple check
			var m = orbit._m = orbit._B < 0 ? (TAU - u) : u;

			// for elliptic orbits
			if (e < 1) {
				// calculate eccentric anomaly
				var E = orbit._E = CYCLE(atan2(
					sqrt(1 - e * e) * sin(m),
					e + cos(m)
				));
				// calculate mean anomaly
				orbit._M = E - e * sin(E);
			}
			// for hyperbolic orbits
			else if (e > 1) {
				// calculate eccentric anomaly
				var E = orbit._E = CYCLE(atan2(
					sqrt(1 - e * e) * sin(m),
					e + cos(m)
				));
				// calculate mean anomaly
				orbit._M = E - e * sin(E);
			}

			// calculate mean anomaly
			orbit._M = E - e * sin(E);

		}
		// EO state vectors

		/********************************************************/
		// VSOP87 uses some exotic orbital elements
		// not sure if they have a technical name
		/********************************************************/

		// VSOP arguments (k/h -> W/e)
		// Called pi (W) in VSOP87
		// k = e*cos(W) [rad]
		// h = e*sin(W) [rad]
		if ('_k' in orbit && '_h' in orbit) {
			// periapsis longitude directly from p and q
			orbit._W = atan2(orbit._h, orbit._k);
			// orbit._e = orbit._k / cos(orbit._W);
			orbit._e = orbit._h / sin(orbit._W);
		}

		// VSOP arguments (q/p -> O/i)
		// Called omega (O) in VSOP87
		// q = sin(i/2)*cos(O) [rad]
		// p = sin(i/2)*sin(O) [rad]
		if ('_q' in orbit && '_p' in orbit) {
			// ascending node directly from p and q
			orbit._O = atan2(orbit._p, orbit._q);
			// values for inclination
			var d = orbit._p - orbit._q,
				// using the faster but equivalent form for
				// dt = sin(orbit._O) - cos(orbit._O);
				dt = - sqrt(2) * sin(PI / 4 - orbit._O);
			// now calculate inclination
			orbit._i = 2 * asin(d / dt);
		}

		/********************************************************/
		// directly related parameters for a
		// more stuff is calculated later on
		/********************************************************/

		if ('_n' in orbit && !('_a' in orbit)) {
			// mean motion is translated directly to size via
			orbit._a = cbrt(orbit._G / orbit._n / orbit._n);
		}
		if ('_P' in orbit && !('_a' in orbit)) {
			// period is translated directly to size via G
			var PTAU = orbit._P / TAU; // reuse for square
			orbit._a = cbrt(orbit._G * PTAU * PTAU);
		}

		/********************************************************/
		// semi-major axis and eccentricity
		// from apocenter and pericenter
		// C = a * (1 + e), c = a * (1 - e)
		/********************************************************/
		if (!('_a' in orbit)) { // semi-major axis
			if ('_c' in orbit && '_C' in orbit) { orbit._a = (orbit._C + orbit._c) / 2; }
			else if ('_e' in orbit && '_C' in orbit) { orbit._a = orbit._C / (1 + orbit._e); }
			else if ('_c' in orbit && '_e' in orbit) { orbit._a = orbit._c / (1 - orbit._e); }
		}
		if (!('_e' in orbit)) { // eccentricity
			if ('_c' in orbit && '_C' in orbit) { orbit._e = 1 - orbit._c / orbit._a; }
			else if ('_a' in orbit && '_C' in orbit) { orbit._e = orbit._C / orbit._a - 1; }
			else if ('_c' in orbit && '_a' in orbit) { orbit._e = 1 - orbit._c / orbit._a; }
		}

		/********************************************************/
		// orbital size parameters
		// a * l = b * b
		// l = a * (1 - e*e)
		// b = a * sqrt(1 - e*e)
		// need 2 independent arguments
		// there are 6 valid combination
		/********************************************************/

		// 3 valid options with e
		if ('_e' in orbit) {
			var e2term = 1 - orbit._e * orbit._e;
			if ('_a' in orbit) {
				orbit._l = orbit._a * e2term;
				orbit._b = sqrt(orbit._a * orbit._l);
			}
			else if ('_b' in orbit) {
				orbit._a = orbit._b / sqrt(e2term);
				orbit._l = orbit._a * e2term;
			}
			else if ('_l' in orbit) {
				orbit._a = orbit._l / e2term;
				orbit._b = sqrt(orbit._a * orbit._l);
			}
		}
		// 2 valid options with a
		else if ('_a' in orbit) {
			if ('_b' in orbit) {
				orbit._l = orbit._b * orbit._b / orbit._a;
				orbit._e = sqrt(1 - orbit._l / orbit._a);
			}
			else if ('_l' in orbit) {
				orbit._e = sqrt(1 - orbit._l / orbit._a);
				orbit._b = sqrt(orbit._a * orbit._l);
			}
		}
		// only one valid options left
		else if ('_b' in orbit && '_l' in orbit) {
			orbit.a = orbit._b * orbit._b / orbit._l;
			orbit._e = sqrt(1 - orbit._l / orbit._a);
		}

		// calculate orbital period (P)
		if ('_a' in orbit && !('_P' in orbit)) {
			// calculate dependants via gravitational parameter
			orbit._P = (TAU / sqrt(orbit._G)) * pow(orbit._a, 1.5);
		}
		// calculate mean motion (n)
		if ('_P' in orbit && !('_n' in orbit)) {
			orbit._n = TAU / orbit._P;
		}
		// specific relative angular momentum (A)
		if ('_P' in orbit && !('_A' in orbit)) {
			if ('_a' in orbit && '_b' in orbit) {
				orbit._A = TAU * orbit._a * orbit._b / orbit._P;
			}
		}

		// apsis from eccentricity and semi-major axis
		if ('_a' in orbit && '_e' in orbit) {
			if (!('_c' in orbit)) orbit._c = orbit._a * (1 - orbit._e);
			if (!('_C' in orbit)) orbit._C = orbit._a * (1 + orbit._e);
		}

		/********************************************************/
		// orientation parameters
		// W = L - M = O + w
		// need 3 independent arguments
		// there are 8 valid combination
		/********************************************************/

		// 5 valid options with L
		if ('_L' in orbit) {
			if ('_M' in orbit) {
				// calculate the main dependants
				orbit._W = orbit._L - orbit._M;
				// calculate additional dependants
				if ('_O' in orbit) orbit._w = orbit._W - orbit._O;
				else if ('_w' in orbit) orbit._O = orbit._W - orbit._w;
			}
			else if ('_W' in orbit) {
				// calculate the main dependants
				orbit._M = orbit._L - orbit._W;
				// calculate additional dependants
				if ('_O' in orbit) orbit._w = orbit._W - orbit._O;
				else if ('_w' in orbit) orbit._O = orbit._W - orbit._w;
			}
			else if ('_O' in orbit && '_w' in orbit) {
				// calculate the main dependants
				orbit._W = orbit._O + orbit._w;
				orbit._M = orbit._L - orbit._W;
			}
			else {
				throw ('Orbit incomplete')
			}
		}
		// 2 valid options with O
		else if ('_O' in orbit) {
			if ('_w' in orbit) {
				// calculate the main dependants
				orbit._W = orbit._O + orbit._w;
				// calculate additional dependants
				if ('_M' in orbit) orbit._L = orbit._M + orbit._W;
				// the L case is already handled in the very first if
				// else if ('_L' in orbit) orbit._M = orbit._L - orbit._W;
			}
			else if ('_W' in orbit) {
				// calculate the main dependants
				orbit._w = orbit._W - orbit._O;
				// calculate additional dependants
				if ('_M' in orbit) orbit._L = orbit._M + orbit._W;
				// the L case is already handled in the very first if
				// else if ('_L' in orbit) orbit._M = orbit._L - orbit._W;
			}
			else {
				throw ('Orbit incomplete')
			}
		}
		// only one valid options left
		else if ('_M' in orbit && '_w' in orbit && '_W' in orbit) {
			orbit._L = orbit._W + orbit._M;
			orbit._O = orbit._W - orbit._w;
		}

		/********************************************************/
		// optionally do expensive calculations
		/********************************************************/

		if (full) orbit.m();

		/********************************************************/
		// do range corrections
		/********************************************************/

		// bring parameters into correct ranges
		if ('_W' in orbit) orbit._W = TURN(orbit._W); // verified
		if ('_L' in orbit) orbit._L = CYCLE(orbit._L); // verified
		if ('_M' in orbit) orbit._M = CYCLE(orbit._M); // verified
		if ('_O' in orbit) orbit._O = CYCLE(orbit._O); // verified
		if ('_w' in orbit) orbit._w = TURN(orbit._w); // verified

		/********************************************************/
		// TODO: better time epoch support
		/********************************************************/
		if ('_M' in orbit && '_n') {
			orbit._T = - orbit._M / orbit._n;
		}

		// chain-able
		return orbit;

	}

	/*############################################################################*/
	// Next functions are needed to calculate eccentricity anomaly and
	// related parameter (i.e. from mean anomaly to true anomaly).
	/*############################################################################*/

	/******************************************************************************/
	// In orbital mechanics, eccentric anomaly (E) is an angular parameter that
	// defines the position of a body that is moving along an elliptic Kepler
	// orbit. The eccentric anomaly (E) is one of three angular parameters
	// ("anomalies") that define a position along an orbit, the other two
	// being the true anomaly (m) and the mean anomaly (M).
	/******************************************************************************/
	Klass.E = function eccentricAnomaly(dt)
	{

		var orbit = this;
		var epoch = orbit._t || 0;

		// return cached
		if ('_E' in orbit && !dt) {
			return orbit._E;
		}

		// basic parameters
		var e = orbit._e, // must
			m = orbit._m, // either
			M = orbit._M; // ... or

		// from true anomaly (m)
		// much easier calculation
		if (!dt && e != null && m != null) {
			return orbit._E = CYCLE(atan2(
				sqrt(1 - e * e) * sin(m),
				e + cos(m)
			));
		}

		// from mean anomaly (M)
		// M = L - W (W = O + w)
		// more expensive solver loop
		if (e != null && M != null) {
			// advance mean anomaly for new time offset
			if (dt) M = CYCLE(orbit._n * (dt - orbit._T - epoch));
			// prepare for solution solver
			var E = e < 0.8 ? M : PI, F;
			var F = E - e * sin(M) - M;
			// Newton-Raphson method to solve
			// f(E) = M - E + e * sin(E) = 0
			var f, dfdE, dE = 1;
			for (var it = 0; abs(dE) > EPSILON && it < MAXLOOP; ++it) {
				f = M - E + e * sin(E);
				dfdE = e * cos(E) - 1.0;
				dE = f / dfdE;
				E -= dE; // next iteration
			}
			// clamp range
			E = CYCLE(E);
			// cache for zero time
			if (!dt) orbit._E = E;
			// return result
			return E;
		}

		// remove in release version
		throw ('Invalid orbital state');

		// fail gracefully
		return null;
	}

	/******************************************************************************/
	// In celestial mechanics, true anomaly is an angular parameter that defines the
	// position of a body moving along a Keplerian orbit. It is the angle between the
	// direction of periapsis and the current position of the body, as seen from the
	// main focus of the ellipse (the point around which the object orbits).
	/******************************************************************************/
	Klass.m = function trueAnomaly(E)
	{
		var orbit = this;
		if (arguments.length === 0) {
			// return cached
			if ('_m' in orbit) {
				return orbit._m;
			}
			// from eccentric anomaly
			// most expensive step
			var hE = orbit.E() / 2;
			// calculate the true anomaly
			return orbit._m = CYCLE(2 * atan2(
				sqrt(1 + orbit._e) * sin(hE),
				sqrt(1 - orbit._e) * cos(hE)
			));
		}
		else {
			// calculate the true anomaly
			return CYCLE(2 * atan2(
				sqrt(1 + orbit._e) * sin(E / 2),
				sqrt(1 - orbit._e) * cos(E / 2)
			));
		}
	}

	/*############################################################################*/
	// additional helper functions (expensive)
	// only use them if you really need too
	// no direct way from orbital elements yet
	/*############################################################################*/

	/******************************************************************************/
	// The radial velocity of an object with respect to a given point is the rate of
	// change of the distance between the object and the point. That is, the radial
	// velocity is the component of the object's velocity that points in the direction
	// of the radius connecting the object and the point. In astronomy, the point is
	// usually taken to be the observer on Earth, so the radial velocity then denotes
	// the speed with which the object moves away from or approaches the Earth.
	/******************************************************************************/
	// Klass.vt = function tangentialVelocity() {}
	/******************************************************************************/
	Klass.B = function radialVelocity()
	{
		var orbit = this;
		// return cached value
		if ('_B' in orbit) return orbit._B;
		// calculate state vectors
		var r = orbit.r(), v = orbit.v();
		// get radial velocity from state vectors
		return orbit._B = r.dot(v) / r.length();
	};

	/******************************************************************************/
	// In celestial mechanics, the eccentricity vector of a Kepler orbit is the
	// dimensionless vector with direction pointing from apoapsis to periapsis and
	// with magnitude equal to the orbit's scalar eccentricity. For Kepler orbits
	// the eccentricity vector is a constant of motion. Its main use is in the
	// analysis of almost circular orbits, as perturbing (non-Keplerian) forces
	// on an actual orbit will cause the osculating eccentricity vector to change
	// continuously. For the eccentricity and argument of periapsis parameters,
	// eccentricity zero (circular orbit) corresponds to a singularity.
	/******************************************************************************/
	Klass.e3 = function eccentricity3()
	{
		var orbit = this;
		// return cached value
		if ('_e3' in orbit) return orbit._e3;
		// force state vector calculation
		var r = orbit.r(), v = orbit.v();
		if (r !== null && v !== null) {
			return orbit._e3 = r.clone().multiplyScalar(
				v.lengthSq() - (orbit._G / r.length())
			).sub(
				v.clone().multiplyScalar(r.length() * orbit.B())
				).multiplyScalar(1 / orbit._G);
		}
	}

	/*############################################################################*/
	// getter functions for elements that are available after resolving
	/*############################################################################*/

	// Orbit inclination is the minimum angle between a reference plane and the
	// orbital plane or axis of direction of an object in orbit around another
	// object. The inclination is one of the six orbital parameters describing
	// the shape and orientation of a celestial orbit. It is the angular distance
	// of the orbital plane from the plane of reference (usually the primary's
	// equator or the ecliptic), normally stated in degrees. In the Solar System,
	// orbital inclination is usually stated with respect to Earth's orbit.
	Klass.i = function inclination() { return this._i; }

	// The orbital eccentricity of an astronomical object is a parameter that
	// determines the amount by which its orbit around another body deviates
	// from a perfect circle. A value of 0 is a circular orbit, values between
	// 0 and 1 form an elliptical orbit, 1 is a parabolic escape orbit, and
	// greater than 1 is a hyperbola. The term derives its name from the
	// parameters of conic sections, as every Kepler orbit is a conic section.
	// It is normally used for the isolated two-body problem, but extensions
	// exist for objects following a rosette orbit through the galaxy.
	Klass.e = function eccentricity() { return this._e; }

	// In geometry, the major axis of an ellipse is its longest diameter:
	// a line segment that runs through the center and both foci, with ends
	// at the widest points of the perimeter. The semi-major axis is on half
	// of the major axis, and thus runs from the centre, through a focus, and
	// to the perimeter. Essentially, it is the radius of an orbit at the
	// orbit's two most distant points. For the special case of a circle, the
	// semi-major axis is the radius. One can think of the semi-major axis as
	// an ellipse's long radius. The length of the semi-major axis a of an
	// ellipse is related to the semi-minor axis's length b through the
	// eccentricity e and the semilatus rectum ℓ.
	Klass.a = function semiMajorAxis() { return this._a; }
	Klass.b = function semiMinorAxis() { return this._b; }
	Klass.l = function semilatusRectum() { return this._l; }

	// The point of a body's elliptical orbit about the system's centre of mass where
	// the distance between the body and the centre of mass is at its maximum.
	// For some celestial bodies, specialised terms are used (aphelion/apogee).
	Klass.C = function apocenter() { return this._C; }
	// The point of a body's elliptical orbit about the system's centre of mass where
	// the distance between the body and the centre of mass is at its minimum.
	// For some celestial bodies, specialised terms are used (perihelion/perigee).
	Klass.c = function pericenter() { return this._c; }

	// Mean longitude is the ecliptic longitude at which an orbiting body
	// could be found if its orbit were circular and free of perturbations.
	// While nominally a simple longitude, in practice the mean longitude
	// is a hybrid angle.
	Klass.L = function meanLongitude() { return this._L; }

	// The longitude of the ascending node (☊ or Ω) is one of the orbital
	// elements used to specify the orbit of an object in space. It is the
	// angle from a reference direction, called the origin of longitude,
	// to the direction of the ascending node, measured in a reference plane.
	// The ascending node is the point where the orbit of the object passes
	// through the plane of reference.
	Klass.O = function ascendingNode() { return this._O; }

	// The argument of periapsis (also called argument of perifocus or
	// argument of pericenter), is one of the orbital elements (ω) of
	// an orbiting body. Parametrically, ω is the angle from the body's
	// ascending node to its periapsis, measured in the direction of motion.
	Klass.w = function argOfPeriapsis() { return this._w; }

	// The time of pericenter passage when the orbiting body passes through
	// the pericenter, closest to the central body. The mean anomaly is 0.0
	// when the orbiting body is at pericenter, so defining the orbital
	// elements at the epoch T (the time of the pericenter passage)
	// eliminates the need to determine the mean anomaly.
	Klass.T = function timeOfPericenter() { return this._T; }

	// In celestial mechanics, the longitude of the periapsis (ϖ) of
	// an orbiting body is the longitude (measured from the point of
	// the vernal equinox) at which the periapsis (closest approach
	// to the central body) would occur if the body's inclination were
	// zero. For motion of a planet around the Sun, this position could
	// be called longitude of perihelion. The longitude of periapsis is
	// a compound angle, with part of it being measured in the plane of
	// reference and the rest being measured in the plane of the orbit.
	// Likewise, any angle derived from the longitude of periapsis (e.g.
	// mean longitude and true longitude) will also be compound.
	Klass.W = function longitudeOfPeriapsis() { return this._W; }

	// In orbital mechanics, mean motion is the angular speed required
	// for a body to complete one orbit, assuming constant speed in a
	// circular orbit which completes in the same time as the variable
	// speed, elliptical orbit of the actual body.
	Klass.n = function meanMotion() { return this._n; }

	// In a closed system, no torque can be exerted on any matter
	// without the exertion on some other matter of an equal and
	// opposite torque. Hence, angular momentum can be exchanged
	// between objects in a closed system, but total angular momentum
	// before and after an exchange remains constant (is conserved)
	Klass.A = function angularMomentum() { return this._A; }

	// In celestial mechanics, the mean anomaly is an angle used in
	// calculating the position of a body in an elliptical orbit in the
	// classical two-body problem. It is the angular distance from the
	// pericenter which a fictitious body would have if it moved in a
	// circular orbit, with constant speed, in the same orbital period
	// as the actual body in its elliptical orbit. It is also the
	// product of mean motion and time since pericenter passage.
	Klass.M = function meanAnomaly() { return this._M; }

	// The orbital period is the time taken for a given object to make one
	// complete orbit around another object. When mentioned without further
	// qualification in astronomy this refers to the sidereal period of an
	// astronomical object, which is calculated with respect to the stars.
	Klass.P = function orbitalPeriod() { return this._P; }

	/*############################################################################*/
	// convert to vsop parameters (never seen them anywhere else)
	// those are calculated on demand, so you need to call me first
	/*############################################################################*/

	// q: sin(i/2)*cos(O) (vsop)
	Klass.q = function vsopPeriapsis()
	{
		var orbit = this;
		// return cached calculation
		if ('_q' in orbit) return orbit._q;
		// calculation given in vsop87 example.f
		return orbit._q = cos(orbit._O)
			* sin(orbit._i / 2);
	};

	// p: sin(i/2)*sin(O) (vsop)
	Klass.p = function vsopApoapsis()
	{
		var orbit = this;
		// return cached calculation
		if ('_p' in orbit) return orbit._p;
		// calculation given in vsop87 example.f
		return orbit._p = sin(orbit._O)
			* sin(orbit._i / 2);
	};

	// k: e*cos(W) (vsop)
	Klass.k = function vsopK()
	{
		var orbit = this;
		// return cached calculation
		if ('_k' in orbit) return orbit._k;
		// calculation given in vsop87 example.f
		return orbit._k = orbit._e * cos(orbit._W);
	};

	// h: e*sin(W) (vsop)
	Klass.h = function vsopH()
	{
		var orbit = this;
		// return cached calculation
		if ('_h' in orbit) return orbit._h;
		// calculation given in vsop87 example.f
		return orbit._h = orbit._e * sin(orbit._W);
	};

	/*############################################################################*/
	// calculate state vectors from orbital elements
	// TODO: implement for different time offsets (TBD)
	/*############################################################################*/

	// update orbital state for new epoch
	// advance position and reset some states
	// for dt = P, mean anomaly does not change
	Klass.update = function update(dt)
	{
		var orbit = this;
		// check if elements are already resolved
		if (!orbit.elements) orbit.resolveElements(true);
		// advance mean anomaly for new time
		orbit._M = CYCLE(orbit._n * (dt - orbit._T));
		// adjust mean longitude for new time
		orbit._L = CYCLE(orbit._M + orbit._W);
		// reset dependent parameters
		delete orbit._E; delete orbit._m;
		// invalidate state vectors
		delete orbit._r; delete orbit._v;
		// chainable
		return orbit;
	}

	/*
	// orbital elements to spherical position (lon/lat/r)
	Klass.sph = function spherical (dt)
	{
		// check if elements are already resolved
		if (!this.elements) this.resolveElements(true);
		dt = dt || 0; // time offset
	}
	*/

	// orbital elements to rectangular position (x/y/z)
	Klass.state = function state(dt)
	{

		var orbit = this, mat;
		var epoch = orbit._t || 0;

		// check if elements are already resolved
		if (!orbit.elements) orbit.resolveElements(true);

		dt = dt || 0; // time offset
		// return cached result for our epoch
		// ToDo: also cache last epoch offsets?
		if (!dt && '_r' in orbit && '_v' in orbit)
		{
			// state result
			return {
				time: dt,
				r: orbit._r,
				v: orbit._v,
				epoch: epoch,
				orbit: orbit
			};
		}

		var e = orbit._e, a = orbit._a, i = orbit._i,
			O = orbit._O, w = orbit._w, M = orbit._M,
			E = orbit.E(dt), // eccentric anomaly
			m = orbit.m(E); // true anomaly

		// Distance to true anomaly position
		var r = a * (1.0 - e * cos(E));
		var vf = sqrt(orbit._G * a) / r;

		// Perifocal reference plane
		var rx = r * cos(m),
			ry = r * sin(m),
			vx = vf * - sin(E),
			vy = vf * sqrt(1.0 - e * e) * cos(E);

		// Pre-calculate elements for rotation matrix
		var sinO = sin(O), cosO = cos(O),
			sinI = sin(i), cosI = cos(i),
			sinW = sin(w), cosW = cos(w),
			sinWcosO = sinW * cosO, sinWsinO = sinW * sinO,
			cosWcosO = cosW * cosO, cosWsinO = cosW * sinO,
			sinWcosI = sinW * cosI,
			FxX = cosWcosO - sinWcosI * sinO,
			FyX = cosWsinO + sinWcosI * cosO,
			FxY = cosWsinO * cosI + sinWcosO,
			FyY = cosWcosO * cosI - sinWsinO,
			FzX = sinW * sinI,
			FzY = cosW * sinI;

		// Equatorial position
		var r = new Vector3(
			rx * FxX - ry * FxY,
			rx * FyX + ry * FyY,
			rx * FzX + ry * FzY
		);
		// Equatorial velocity
		var v = new Vector3(
			vx * FxX - vy * FxY,
			vx * FyX + vy * FyY,
			vx * FzX + vy * FzY
		);
		// cache results if no time offset
		if (!dt) orbit._r = r, orbit._v = v;
		// return result object with reference
		// ToDo: maybe add position object
		// To calculate ra/dec and more stuff
		return {
			r: r,
			v: v,
			time: dt,
			epoch: epoch,
			orbit: orbit
		};

	}
	// EO state

	/*############################################################################*/
	// getters for state vectors
	/*############################################################################*/

	// r: position state vector
	Klass.r = function position3(dt)
	{
		// calculate at epoch time
		// state will cache results
		var state = this.state(dt);
		// return vector
		return state.r;
	}

	// v: velocity state vector
	Klass.v = function velocity3(dt)
	{
		// calculate at epoch time
		// state will cache results
		var state = this.state(dt);
		// return vector
		return state.v;
	}

	/*############################################################################*/
	// add fully named getter functions if preferred
	/*############################################################################*/

	// orbital elements that are always calculated
	// these can be used as orbital input arguments
	Klass.inclination = Klass.i;
	Klass.eccentricity = Klass.e;
	Klass.semiMajorAxis = Klass.a;
	Klass.semiMinorAxis = Klass.b;
	Klass.semilatusRectum = Klass.l;
	Klass.meanLongitude = Klass.L;
	Klass.ascendingNode = Klass.O;
	Klass.argOfPericenter = Klass.w;
	Klass.argOfPeriapsis = Klass.w;
	Klass.timeOfPericenter = Klass.T;
	Klass.longitudeOfPericenter = Klass.W;
	Klass.meanMotion = Klass.n;
	Klass.meanAnomaly = Klass.M;
	Klass.orbitalPeriod = Klass.P;
	Klass.angularMomentum = Klass.A;

	// additional orbital elements only on demand
	// these must not be used as orbital input arguments
	Klass.trueAnomaly = Klass.m;
	Klass.radialVelocity = Klass.B;
	Klass.eccentricAnomaly = Klass.E;
	// additional parameters are calculated on demand
	Klass.apoapsis = Klass.C;
	Klass.aphelion = Klass.C;
	Klass.apocenter = Klass.C;
	Klass.periapsis = Klass.c;
	Klass.perhelion = Klass.c;
	Klass.pericenter = Klass.c;

	// state vectors at epoch time
	Klass.position = Klass.r;
	Klass.velocity = Klass.v;

	/*############################################################################*/
	// END OF AstroJS Orbit Module
	/*############################################################################*/

	exports.Orbit = exports.Orbital = Orbit;

})(this);
