/*############################################################################*/
// AstroJS Orbit Module (c) 2016 by Marcel Greter
// https://www.github.com/mgreter/astrojs/LICENSE
/*############################################################################*/
// http://www.lns.cornell.edu/~seb/celestia/orbital-parameters.html
// Similar module: https://github.com/jordanstephens/kepler.js
/*############################################################################*/
'use strict';

(function (exports) {

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

		// internal use
		if (arg.empty) {
			return this;
		}

		// the gravitational parameter is very important and must
		// match all others in terms of the units. In physics the
		// units are by default m^3/(kg*s^2), but for astronomical
		// uses we most often want to use au^3/(solm*day^2).
		this._G = GMP.sun;

		// create state vectors from input parameters
		if ('x' in arg && 'y' in arg && 'z' in arg) {
			this._r = new Vector3(arg.x, arg.y, arg.z);
		}
		else if ('rx' in arg && 'ry' in arg && 'rz' in arg) {
			this._r = new Vector3(arg.rx, arg.ry, arg.rz);
		}
		if ('X' in arg && 'Y' in arg && 'Z' in arg) {
			this._v = new Vector3(arg.X, arg.Y, arg.Z);
		}
		else if ('vx' in arg && 'vy' in arg && 'vz' in arg) {
			this._v = new Vector3(arg.vx, arg.vy, arg.vz);
		}

		if ('r' in arg && 'v' in arg) {
			var r = arg.r, v = arg.v;
			// create new vectors (want a clone anyway)
			this._r = new Vector3(r.x, r.y, r.z);
			this._v = new Vector3(v.x, v.y, v.z);
		}

		// set status for state vectors
		this.vectors = !!(this._r && this._v);

		// get other parameters from arg
		// unknown parameters are ignored
		var vars = 'aeiMnPLcCwWOkhqpG'.split('');
		for (var i in vars) {
			if (arg.hasOwnProperty(vars[i])) {
				this['_' + vars[i]] = arg[vars[i]];
			}
		}

		// orbit may has name
		if ('name' in arg) {
			this.name = arg.name;
		}

		// resolve if not lazy
		if (!this.lazy) {
			this.updateElements();
			this.requireElements();
		}

	}
	// EO Orbit ctor

	// optimized for minifier
	var Klass = Orbit.prototype;

	/******************************************************************************/
	// clone the existing orbital
	/******************************************************************************/

	// copy existing orbit and return clone
	Klass.clone = function clone(dt)
	{
		// create an empty object (internal use)
		var clone = new Orbit({ empty: true });
		// avoid unnecessary cloning of objects
		if (dt) { delete this._v; delete this._r; }
		// process all properties
		for (var key in this) {
			if (this.hasOwnProperty(key)) {
				// invoke nested clone functions (vectors)
				if (typeof this[key].clone == 'function')
				{ clone[key] = this[key].clone(); }
				else { clone[key] = this[key]; }
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
		// most basic check for initial pre-calculations
		if (!this.elements) throw('Orbitals not calculated');
		// do some basic checks for needed orbital elements
		if (!('_i' in this)) throw('Orbit is missing inclination (i)');
		if (!('_e' in this)) throw('Orbit is missing eccentricity (e)');
		if (!('_a' in this)) throw('Orbit is missing semi-major axis (a)');
		if (!('_b' in this)) throw('Orbit is missing semi-minor axis (b)');
		if (!('_l' in this)) throw('Orbit is missing semi-latus rectum (l)');
		if (!('_c' in this)) throw('Orbit is missing periapsis (c)');
		if (!('_C' in this)) throw('Orbit is missing apoapsis (C)');
		if (!('_L' in this)) throw('Orbit is missing mean longitude (L)');
		if (!('_O' in this)) throw('Orbit is missing right ascending node (O)');
		if (!('_T' in this)) throw('Orbit is missing time of periapsis (w)');
		if (!('_w' in this)) throw('Orbit is missing argument of periapsis (w)');
		if (!('_W' in this)) throw('Orbit is missing longitude of the periapsis (W)');
		if (!('_n' in this)) throw('Orbit is missing mean motion (n)');
		if (!('_M' in this)) throw('Orbit is missing mean anomaly (M)');
		if (!('_P' in this)) throw('Orbit is missing orbital period (P)');
		if (!('_A' in this)) throw('Orbit is missing angular momentum (A)');
		if (!('_c' in this)) throw('Orbit is missing pericenter (c)');
		if (!('_C' in this)) throw('Orbit is missing apocenter (C)');
		// check if eccentricity is inside valid boundaries (0-1)
		if (this._e < 0) throw('Negative eccentricity is invalid');
		if (this._e > 1) throw('Eccentricity must not be hyperbolic (> 1)');
		if (this._e == 1) throw('Eccentricity must not be parabolic (== 1)');
		// state vectors (r and v) are not tested here ...
		if (full) if (!('_m' in this)) throw('Orbit is missing true anomaly (m)');
		if (full) if (!('_B' in this)) throw('Orbit is missing radial velocity (B)');
		if (full) if (!('_E' in this)) throw('Orbit is missing eccentric anomaly (E)');
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

		// only resolve elements once
		if (this.elements) return this;
		// will calculate now
		this.elements = true;

		/********************************************************/
		// resolve from state vectors (r, v)
		/********************************************************/

		// check main state vectors
		if ('_r' in this && '_v' in this) {

			var r = this._r, v = this._v,
			    G = this._G, rl = r.length();

			// specific relative angular momentum
			this._A3 = r.clone().cross(v);
			this._A2 = this._A3.lengthSq();
			this._A = Math.sqrt(this._A2);
			// calculate radial velocity
			this._B = r.dot(v) / rl;

			// calculate eccentricity vector
			var e3 = this._e3 = r.clone().multiplyScalar(v.lengthSq() - (G / rl))
				.sub(v.clone().multiplyScalar(rl * this._B)).multiplyScalar(1 / G);
			// get eccentricity value from vector
			var e2 = e3.lengthSq(), e = this._e = Math.sqrt(e2);
			// get inclination (i) via orbital momentum vector
			this._i = Math.acos(this._A3.z / this._A);

			// calculate semilatus rectum (ℓ)
			this._l = this._A2 / G;
			// and periapsis (c) and apoapsis (C)
			this._c = this._l / (1 + e);
			this._C = this._l / (1 - e);
			// and finally semi-major axis (a)
			this._a = (this._C + this._c) / 2;

			// pre-calculate node line
			var nx = - this._A3.y,
			    ny = + this._A3.x;
			var nl = Math.sqrt(nx*nx + ny*ny);

			// calculate ascending node (O)
			var omega = nl == 0 ? 0 : Math.acos(nx / nl);
			this._O = ny < 0 ? (TAU - omega) : omega;

			// calculate argument of periapsis (ω)
			var nedot = nx * e3.x +
			            ny * e3.y;
			if (nl === 0 || e === 0) { this.w = 0; }
			else { this._w = Math.acos(nedot / nl / e); }
			if (e3.z < 0) { this._w *= -1; }

			// calculate true anomaly
			var u; // argument of latitude
			// case for circular orbit
			// and without inclination
			if (e === 0 && nl === 0) {
				// this needs a test case
				u = Math.acos(r.x / rl);
			}
			// circular orbit
			// with inclination
			else if (e === 0) {
				// this needs a test case
				var nrdot = nx * r.x + ny * r.y;
				u = Math.acos(nrdot / (nl * rl));
			}
			// elliptic orbit
			else {
				var redot = e3.x * r.x + e3.y * r.y + e3.z * r.z;
				u = Math.acos(redot / e / rl);
			}
			// bring into correct range via simple check
			var m = this._m = this._B < 0 ? (TAU - u) : u;

			// calculate eccentric anomaly
			var E = this._E = CYCLE(Math.atan2(
				Math.sqrt(1 - e*e) * Math.sin(m),
				e + Math.cos(m)
			));

			// calculate mean anomaly
			this._M = E - e * Math.sin(E);

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
		if ('_k' in this && '_h' in this) {
			// periapsis longitude directly from p and q
			this._W = Math.atan2(this._h, this._k);
			// this._e = this._k / Math.cos(this._W);
			this._e = this._h / Math.sin(this._W);
		}

		// VSOP arguments (q/p -> O/i)
		// Called omega (O) in VSOP87
		// q = sin(i/2)*cos(O) [rad]
		// p = sin(i/2)*sin(O) [rad]
		if ('_q' in this && '_p' in this) {
			// ascending node directly from p and q
			this._O = Math.atan2(this._p, this._q);
			// values for inclination
			var d = this._p - this._q,
				// using the faster but equivalent form for
				// dt = Math.sin(this._O) - Math.cos(this._O);
				dt = - Math.sqrt(2) * Math.sin(PI / 4 - this._O);
			// now calculate inclination
			this._i = 2 * Math.asin(d / dt);
		}

		/********************************************************/
		// directly related parameters for a
		// more stuff is calculated later on
		/********************************************************/

		if ('_n' in this && !('_a' in this)) {
			// mean motion is translated directly to size via
			this._a = Math.cbrt(this._G / this._n / this._n);
		}
		if ('_P' in this && !('_a' in this)) {
			// period is translated directly to size via G
			var PTAU = this._P / TAU; // reuse for square
			this._a = Math.cbrt(this._G * PTAU * PTAU);
		}

		/********************************************************/
		// semi-major axis and eccentricity
		// from apocenter and pericenter
		// C = a * (1 + e), c = a * (1 - e)
		/********************************************************/
		if (!('_a' in this)) { // semi-major axis
			if ('_c' in this && '_C' in this) { this._a = (this._C + this._c) / 2; }
			else if ('_e' in this && '_C' in this) { this._a = this._C / (1 + this._e); }
			else if ('_c' in this && '_e' in this) { this._a = this._c / (1 - this._e); }
		}
		if (!('_e' in this)) { // eccentricity
			if ('_c' in this && '_C' in this) { this._e = 1 - this._c / this._a; }
			else if ('_a' in this && '_C' in this) { this._e = this._C / this._a - 1; }
			else if ('_c' in this && '_a' in this) { this._e = 1 - this._c / this._a; }
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
		if ('_e' in this) {
			var e2term = 1 - this._e * this._e;
			if ('_a' in this) {
				this._l = this._a * e2term;
				this._b = Math.sqrt(this._a * this._l);
			}
			else if ('_b' in this) {
				this._a = this._b / Math.sqrt(e2term);
				this._l = this._a * e2term;
			}
			else if ('_l' in this) {
				this._a = this._l / e2term;
				this._b = Math.sqrt(this._a * this._l);
			}
		}
		// 2 valid options with a
		else if ('_a' in this) {
			if ('_b' in this) {
				this._l = this._b * this._b / this._a;
				this._e = Math.sqrt(1 - this._l / this._a);
			}
			else if ('_l' in this) {
				this._e = Math.sqrt(1 - this._l / this._a);
				this._b = Math.sqrt(this._a * this._l);
			}
		}
		// only one valid options left
		else if ('_b' in this && '_l' in this) {
			this.a = this._b * this._b / this._l;
			this._e = Math.sqrt(1 - this._l / this._a);
		}

		// calculate orbital period (P)
		if ('_a' in this && !('_P' in this)) {
			// calculate dependants via gravitational parameter
			this._P = (TAU / Math.sqrt(this._G)) * Math.pow(this._a, 1.5);
		}
		// calculate mean motion (n)
		if ('_P' in this && !('_n' in this)) {
			this._n = TAU / this._P;
		}
		// specific relative angular momentum (A)
		if ('_P' in this && !('_A' in this)) {
			if ('_a' in this && '_b' in this) {
				this._A = TAU * this._a * this._b / this._P;
			}
		}

		// apsis from eccentricity and semi-major axis
		if ('_a' in this && '_e' in this) {
			if (!('_c' in this)) this._c = this._a * ( 1 - this._e );
			if (!('_C' in this)) this._C = this._a * ( 1 + this._e );
		}

		/********************************************************/
		// orientation parameters
		// W = L - M = O + w
		// need 3 independent arguments
		// there are 8 valid combination
		/********************************************************/

		// 5 valid options with L
		if ('_L' in this) {
			if ('_M' in this) {
				// calculate the main dependants
				this._W = this._L - this._M;
				// calculate additional dependants
				if ('_O' in this) this._w = this._W - this._O;
				else if ('_w' in this) this._O = this._W - this._w;
			}
			else if ('_W' in this) {
				// calculate the main dependants
				this._M = this._L - this._W;
				// calculate additional dependants
				if ('_O' in this) this._w = this._W - this._O;
				else if ('_w' in this) this._O = this._W - this._w;
			}
			else if ('_O' in this && '_w' in this) {
				// calculate the main dependants
				this._W = this._O + this._w;
				this._M = this._L - this._W;
			}
			else {
				throw('Orbit incomplete')
			}
		}
		// 2 valid options with O
		else if ('_O' in this) {
			if ('_w' in this) {
				// calculate the main dependants
				this._W = this._O + this._w;
				// calculate additional dependants
				if ('_M' in this) this._L = this._M + this._W;
				// the L case is already handled in the very first if
				// else if ('_L' in this) this._M = this._L - this._W;
			}
			else if ('_W' in this) {
				// calculate the main dependants
				this._w = this._W - this._O;
				// calculate additional dependants
				if ('_M' in this) this._L = this._M + this._W;
				// the L case is already handled in the very first if
				// else if ('_L' in this) this._M = this._L - this._W;
			}
			else {
				throw('Orbit incomplete')
			}
		}
		// only one valid options left
		else if ('_M' in this && '_w' in this && '_W' in this) {
			this._L = this._W + this._M;
			this._O = this._W - this._w;
		}

		/********************************************************/
		// optionally do expensive calculations
		/********************************************************/

		if (full) this.m();

		/********************************************************/
		// do range corrections
		/********************************************************/

		// bring parameters into correct ranges
		if ('_W' in this) this._W = TURN(this._W); // verified
		if ('_L' in this) this._L = CYCLE(this._L); // verified
		if ('_M' in this) this._M = CYCLE(this._M); // verified
		if ('_O' in this) this._O = CYCLE(this._O); // verified
		if ('_w' in this) this._w = TURN(this._w); // verified

		/********************************************************/
		// TODO: better time epoch support
		/********************************************************/
		if ('_M' in this && '_n') {
			this._T = - this._M / this._n;
		}

		// chain-able
		return this;

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

		// return cached
		if ('_E' in this && !dt) {
			return this._E;
		}

		// basic parameters
		var e = this._e, // must
		    m = this._m, // either
		    M = this._M; // ... or

		// from true anomaly (m)
		// much easier calculation
		if (!dt && e != null && m != null) {
			return this._E = CYCLE(Math.atan2(
				Math.sqrt(1 - e*e) * Math.sin(m),
				e + Math.cos(m)
			));
		}

		// from mean anomaly (M)
		// M = L - W (W = O + w)
		// more expensive solver loop
		if (e != null && M != null) {
			// advance mean anomaly for new time offset
			if (dt) M = CYCLE(this._n * (dt - this._T));
			// prepare for solution solver
			var E = e < 0.8 ? M : PI, F;
			var F = E - e * Math.sin(M) - M;
			// Newton-Raphson method to solve
			// f(E) = M - E + e * sin(E) = 0
			var f, dfdE, dE = 1;
			for (var it = 0; Math.abs(dE) > EPSILON && it < MAXLOOP; ++it) {
				f = M - E + e * Math.sin(E);
				dfdE = e * Math.cos(E) - 1.0;
				dE = f / dfdE;
				E -= dE; // next iteration
			}
			// clamp range
			E = CYCLE(E);
			// cache for zero time
			if (!dt) this._E = E;
			// return result
			return E;
		}

		// remove in release version
		throw('Invalid orbital state');

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
		if (arguments.length === 0) {
			// return cached
			if ('_m' in this) {
				return this._m;
			}
			// from eccentric anomaly
			// most expensive step
			var hE = this.E() / 2;
			// calculate the true anomaly
			return this._m = CYCLE(2 * Math.atan2(
				Math.sqrt(1+this._e) * Math.sin(hE),
				Math.sqrt(1-this._e) * Math.cos(hE)
			));
		}
		else {
			// calculate the true anomaly
			return CYCLE(2 * Math.atan2(
				Math.sqrt(1+this._e) * Math.sin(E/2),
				Math.sqrt(1-this._e) * Math.cos(E/2)
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
		// return cached value
		if ('_B' in this) return this._B;
		// calculate state vectors
		var r = this.r(), v = this.v();
		// get radial velocity from state vectors
		return this._B = r.dot(v) / r.length();
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
		// return cached value
		if ('_e3' in this) return this._e3;
		// force state vector calculation
		var r = this.r(), v = this.v();
		if (r !== null && v !== null) {
			return this._e3 = r.clone().multiplyScalar(
				v.lengthSq() - (this._G / r.length())
			).sub(
				v.clone().multiplyScalar(r.length() * this.B())
			).multiplyScalar(1 / this._G);
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
		// return cached calculation
		if ('_q' in this) return this._q;
		// calculation given in vsop87 example.f
		return this._q = Math.cos(this._O)
			* Math.sin(this._i / 2);
	};

	// p: sin(i/2)*sin(O) (vsop)
	Klass.p = function vsopApoapsis()
	{
		// return cached calculation
		if ('_p' in this) return this._p;
		// calculation given in vsop87 example.f
		return this._p = Math.sin(this._O)
			* Math.sin(this._i / 2);
	};

	// k: e*cos(W) (vsop)
	Klass.k = function vsopK()
	{
		// return cached calculation
		if ('_k' in this) return this._k;
		// calculation given in vsop87 example.f
		return this._k = this._e * Math.cos(this._W);
	};

	// h: e*sin(W) (vsop)
	Klass.h = function vsopH()
	{
		// return cached calculation
		if ('_h' in this) return this._h;
		// calculation given in vsop87 example.f
		return this._h = this._e * Math.sin(this._W);
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
		// check if elements are already resolved
		if (!this.elements) this.resolveElements(true);
		// advance mean anomaly for new time
		this._M = CYCLE(this._n * (dt - this._T));
		// adjust mean longitude for new time
		this._L = CYCLE(this._M + this._W);
		// reset dependent parameters
		delete this._E; delete this._m;
		// invalidate state vectors
		delete this._r; delete this._v;
		// chainable
		return this;
	}

	// orbital elements to spherical position (lon/lat/d)
	Klass.state = function spherical (dt)
	{
	}

	// orbital elements to rectangular position (x/y/z)
	Klass.state = function state (dt)
	{

		// check if elements are already resolved
		if (!this.elements) this.resolveElements(true);

		var e = this._e, a = this._a, i = this._i,
		    O = this._O, w = this._w, M = this._M,
		    E = this.E(), // eccentric anomaly
		    m = this.m(), // true anomaly
		    t = dt || 0; // time offset

		// Distance to true anomaly position
		var r = a * (1.0 - e * Math.cos(E));
		var vf = Math.sqrt(this._G * a) / r;

		// Perifocal reference plane
		var rx = r * Math.cos(m),
		    ry = r * Math.sin(m),
		    vx = vf * - Math.sin(E),
		    vy = vf * Math.sqrt(1.0 - e*e) * Math.cos(E);

		// Pre-calculate elements for rotation matrix
		var sinO = Math.sin(O), cosO = Math.cos(O),
		    sinI = Math.sin(i), cosI = Math.cos(i),
		    sinW = Math.sin(w), cosW = Math.cos(w),
		    sinWcosO = sinW*cosO, sinWsinO = sinW*sinO,
		    cosWcosO = cosW*cosO, cosWsinO = cosW*sinO,
		    FxX = (cosW*cosO - sinW*sinO*cosI),
		    FyX = (cosW*sinO + sinW*cosO*cosI),
		    FxY = (cosW*sinO*cosI + sinW*cosO),
		    FyY = (cosW*cosO*cosI - sinW*sinO),
		    FzX = (sinW*sinI),
		    FzY = (cosW*sinI);

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
		if (!t) this._r = r, this._v = v;
		// return result object with reference
		// ToDo: maybe add position object
		// To calculate ra/dec and more stuff
		return { r: r, v: v, time: t, orbit: this };

	}
	// EO state

	/*############################################################################*/
	// getters for state vectors
	/*############################################################################*/

	// r: position state vector
	Klass.r = function position3()
	{
		// return cached
		if ('_r' in this) {
			return this._r;
		}
		// calculate at epoch time
		var state = this.state(0);
		// cache the result
		this._r = state.r;
		this._v = state.v;
		// return vector
		return state.r;
	}

	// v: velocity state vector
	Klass.v = function velocity3()
	{
		// return cached
		if ('_v' in this) {
			return this._v;
		}
		// calculate at epoch time
		var state = this.state(0);
		// cache the result
		this._r = state.r;
		this._v = state.v;
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
