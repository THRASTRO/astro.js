# AstroJS

JavaScript Library for Astronomical Calculations.
Related to https://github.com/mgreter/ephem.js

# Orbit.js

Orbits can be created from orbital elements (6 independent parameters)
or from state vectors (position and velocity). In both cases we need 6
arguments to fully define an orbit. Most related parameters will be
calculated on construction. But certain parameters as the true anomaly
are a bit more complex to calculate from orbital elements, if only the
mean anomaly is given as initial parameters. They are lazy calculated
on demand, either by using the getter or invoking updateElements(true).

To get state vectors for different times we just add a linear factor
to the mean anomaly (related to mean motion and period). To get the
actual state vectors from the mean anomaly we need to compute the
solution for Kepler's Equation via Newton-Raphson solver loop.

Keep in mind that all angular parameters are in rads and that all
other units are linked to each other. This means that you can change
the units by adjusting the gravitational parameter. We use the solar
gravitational parameter by default with the units: `au^3/(solm*day^2)`.
Normally in physics this parameter is measured in `m^3/(kg*s^2)`. All
distance and time units must therefore be in the same units as `G`.

## Orbit Constructors

```js
// from standard orbital elements
// https://en.wikipedia.org/wiki/Orbital_elements
orbits.push({
	name: 'Planet',
	e: 0.09331680132087833, // unitless
	a: 1.5236773948307722, // au
	i: 0.032283817314135946, // rad
	O: 0.8649518974751991, // rad
	W: -0.4171558510055842, // rad
	L: 6.203874882756421 // rad
	// default G is the sun with au and d
	// G: 2.9591220836841438269e-04
});

// satellite orbiting earth
var orbit = new Orbit({
	x: 5052.4587, // km
	y: 1056.2713, // km
	z: 5011.6366, // km
	X: 3.8589872, // km/s
	Y: 4.2763114, // km/s
	Z: -4.8070493, // km/s
	// this G is for the earth in km and s
	G: 398600.44
});

// from orbital state vectors
// https://en.wikipedia.org/wiki/Orbital_state_vectors
new Orbit({
	name: 'State',
	rx: 1.390712873654179, // au
	ry: -0.013411960276849655, // au
	rz: -0.034462732955625755, // au
	vx: 0.0006714786862244502, // au/d
	vy: 0.015187272121968753, // au/d
	vz: 0.0003016547612144848, // au/d
});

// from VSOP parameters
// https://en.wikipedia.org/wiki/VSOP_(planets)
new Orbit({
	name: 'VSOP',
	a: 1.5236773948307722, // au
	L: 6.203874882756421, // rad
	h: -0.037808407518263365,
	k: 0.08531441689241746,
	q: 0.01047042574,
	p: 0.01228449307
});
```

## Query Orbital Elements

After construction, most parameters are already calculated. You may access them
directly, but it's recommended call their getters instead, to ensure  the
requested element was really calculated. If you want the most performance
you may access them directly though (Just be careful!).

Pretty much every orbital element is stored as a one char parameter. Thus it
is important that you make sure you access the correct variable, it is easy
to confuse certain variables. Consult the following list if needed:

## Variable definitions

### orbital elements that are always calculated
These can be used as orbital input arguments
```js
i: inclination
e: eccentricity
a: semiMajorAxis
b: semiMinorAxis
l: semilatusRectum
L: meanLongitude
O: ascendingNode
w: argOfPeriapsis
T: timeOfPericenter
W: longitudeOfPericenter
n: meanMotion
M: meanAnomaly
P: orbitalPeriod
A: angularMomentum
```

### additional orbital elements only on demand
These must not be used as orbital input arguments
```js
m: trueAnomaly
B: radialVelocity
E: eccentricAnomaly
```

### additional parameters are calculated on demand
Only access them after calling their method once
```js
C: apoapsis, aphelion, apocenter
c: periapsis, perhelion, pericenter
```

### state vectors at epoch time
Available after calling `state` method.
```js
r: position
v: velocity
```

### `E = function eccentricAnomaly(dt)`
In orbital mechanics, eccentric anomaly (E) is an angular parameter that
defines the position of a body that is moving along an elliptic Kepler
orbit. The eccentric anomaly (E) is one of three angular parameters
("anomalies") that define a position along an orbit, the other two
being the true anomaly (m) and the mean anomaly (M).

### `m = function trueAnomaly(E)`
In celestial mechanics, true anomaly is an angular parameter that defines the
position of a body moving along a Keplerian orbit. It is the angle between the
direction of periapsis and the current position of the body, as seen from the
main focus of the ellipse (the point around which the object orbits).

### `B = function radialVelocity()`
The radial velocity of an object with respect to a given point is the rate of
change of the distance between the object and the point. That is, the radial
velocity is the component of the object's velocity that points in the direction
of the radius connecting the object and the point. In astronomy, the point is
usually taken to be the observer on Earth, so the radial velocity then denotes
the speed with which the object moves away from or approaches the Earth.

### `e3 = function eccentricity3()`
In celestial mechanics, the eccentricity vector of a Kepler orbit is the
dimensionless vector with direction pointing from apoapsis to periapsis and
with magnitude equal to the orbit's scalar eccentricity. For Kepler orbits
the eccentricity vector is a constant of motion. Its main use is in the
analysis of almost circular orbits, as perturbing (non-Keplerian) forces
on an actual orbit will cause the osculating eccentricity vector to change
continuously. For the eccentricity and argument of periapsis parameters,
eccentricity zero (circular orbit) corresponds to a singularity.

### `i = function inclination()`
Orbit inclination is the minimum angle between a reference plane and the
orbital plane or axis of direction of an object in orbit around another
object. The inclination is one of the six orbital parameters describing
the shape and orientation of a celestial orbit. It is the angular distance
of the orbital plane from the plane of reference (usually the primary's
equator or the ecliptic), normally stated in degrees. In the Solar System,
orbital inclination is usually stated with respect to Earth's orbit.

### `e = function eccentricity()`
The orbital eccentricity of an astronomical object is a parameter that
determines the amount by which its orbit around another body deviates
from a perfect circle. A value of 0 is a circular orbit, values between
0 and 1 form an elliptical orbit, 1 is a parabolic escape orbit, and
greater than 1 is a hyperbola. The term derives its name from the
parameters of conic sections, as every Kepler orbit is a conic section.
It is normally used for the isolated two-body problem, but extensions
exist for objects following a rosette orbit through the galaxy.

### `a = function semiMajorAxis()`
In geometry, the major axis of an ellipse is its longest diameter:
a line segment that runs through the center and both foci, with ends
at the widest points of the perimeter. The semi-major axis is on half
of the major axis, and thus runs from the centre, through a focus, and
to the perimeter. Essentially, it is the radius of an orbit at the
orbit's two most distant points. For the special case of a circle, the
semi-major axis is the radius. One can think of the semi-major axis as
an ellipse's long radius. The length of the semi-major axis a of an
ellipse is related to the semi-minor axis's length b through the
eccentricity e and the semilatus rectum ℓ.

### `b = function semiMinorAxis()`
The semi-minor axis b is related to the semi-major axis a through the
eccentricity e and the semi-latus rectum ℓ.

### `l = function semilatusRectum()`
The semi-latus rectum ℓ is related to the semi-major axis a through the
eccentricity e and the semi-minor axis b.

### `C = function apocenter()`
The point of a body's elliptical orbit about the system's centre of mass where
the distance between the body and the centre of mass is at its maximum.
For some celestial bodies, specialised terms are used (aphelion/apogee).

### `c = function pericenter()`
The point of a body's elliptical orbit about the system's centre of mass where
the distance between the body and the centre of mass is at its minimum.
For some celestial bodies, specialised terms are used (perihelion/perigee).

### `L = function meanLongitude()`
Mean longitude is the ecliptic longitude at which an orbiting body
could be found if its orbit were circular and free of perturbations.
While nominally a simple longitude, in practice the mean longitude
is a hybrid angle.

### `O = function ascendingNode()`
The longitude of the ascending node (☊ or Ω) is one of the orbital
elements used to specify the orbit of an object in space. It is the
angle from a reference direction, called the origin of longitude,
to the direction of the ascending node, measured in a reference plane.
The ascending node is the point where the orbit of the object passes
through the plane of reference.

### `w = function argOfPeriapsis()`
The argument of periapsis (also called argument of perifocus or
argument of pericenter), is one of the orbital elements (ω) of
an orbiting body. Parametrically, ω is the angle from the body's
ascending node to its periapsis, measured in the direction of motion.

### `T = function timeOfPericenter()`
The time of pericenter passage when the orbiting body passes through
the pericenter, closest to the central body. The mean anomaly is 0.0
when the orbiting body is at pericenter, so defining the orbital
elements at the epoch T (the time of the pericenter passage)
eliminates the need to determine the mean anomaly.

### `W = function longitudeOfPeriapsis()`
In celestial mechanics, the longitude of the periapsis (ϖ) of
an orbiting body is the longitude (measured from the point of
the vernal equinox) at which the periapsis (closest approach
to the central body) would occur if the body's inclination were
zero. For motion of a planet around the Sun, this position could
be called longitude of perihelion. The longitude of periapsis is
a compound angle, with part of it being measured in the plane of
reference and the rest being measured in the plane of the orbit.
Likewise, any angle derived from the longitude of periapsis (e.g.
mean longitude and true longitude) will also be compound.

### `n = function meanMotion()`
In orbital mechanics, mean motion is the angular speed required
for a body to complete one orbit, assuming constant speed in a
circular orbit which completes in the same time as the variable
speed, elliptical orbit of the actual body.

### `A = function angularMomentum()`
In a closed system, no torque can be exerted on any matter
without the exertion on some other matter of an equal and
opposite torque. Hence, angular momentum can be exchanged
between objects in a closed system, but total angular momentum
before and after an exchange remains constant (is conserved)

### `M = function meanAnomaly()`
In celestial mechanics, the mean anomaly is an angle used in
calculating the position of a body in an elliptical orbit in the
classical two-body problem. It is the angular distance from the
pericenter which a fictitious body would have if it moved in a
circular orbit, with constant speed, in the same orbital period
as the actual body in its elliptical orbit. It is also the
product of mean motion and time since pericenter passage.

### `P = function orbitalPeriod()`
The orbital period is the time taken for a given object to make one
complete orbit around another object. When mentioned without further
qualification in astronomy this refers to the sidereal period of an
astronomical object, which is calculated with respect to the stars.

### `epoch = function epoch()`
In chronology and periodization, an epoch or reference epoch is an
instant in time chosen as the origin of a particular calendar era.
The "epoch" serves as a reference point from which time is measured. 

### VSOP parameters
```
// q: sin(i/2)*cos(O) (vsop)
q = function vsopPeriapsis()
// p: sin(i/2)*sin(O) (vsop)
p = function vsopApoapsis()
// k: e*cos(W) (vsop)
k = function vsopK()
// h: e*sin(W) (vsop)
h = function vsopH()
```

Descriptions taken from https://en.wikipedia.org

# Time.js

### `getMeanSiderealTime(JD)`

Calculate the apparent sidereal time at the meridian of Greenwich of a given date.

### `getApparentSiderealTime(JD)`

Calculate the mean sidereal time at the meridian of Greenwich of a given date.

### `getNutation(JD)`

Calculate nutation of longitude and obliquity in degrees from Julian Ephemeris Day

One to one implementations from [libnova][1]

[1]: http://libnova.sourceforge.net/

# Coord.js

- Cartesian (x, y, z)
- Spherical (r, θ, φ)
- Cylindrical (ρ, φ, z)

Some variables are the same in multiple representations (`z` and `φ`). For the
conversion between cylindrical and Cartesian coordinates, it is convenient to
assume that the reference plane of the former is the Cartesian x–y plane (with
equation z = 0), and the cylindrical axis is the Cartesian z axis. Then the z
coordinate is the same in both systems, and the correspondence between
cylindrical (ρ,φ) and Cartesian (x,y) are the same as for polar coordinates.

## Constructor

```js
// cartesian (preferred)
new Coord({
	x: 4.001168777115359,
	y: 2.9385853724645656,
	z: -0.10178505133708564
});

// spherical with elevation
// this is the VSOP87 format
new Coord({
	r: 4.9653867444245865, // radius (r)
	l: 36.29468623590961 * DEG2RAD, // phi (φ)
	b: -1.1745809849948283 * DEG2RAD, // theta (θ)
});

// spherical with inclination
new Coord({
	r: 4.9653797207164, // radius (r)
	l: 36.29468623590961 * DEG2RAD, // phi (φ)
	i: 91.17458536812916 * DEG2RAD, // theta (θ)
});

// cylindrical
new Coord({
	p: 4.9643363679575, // radial dist (ρ)
	l: 36.29474187359682 * DEG2RAD, // azimuth (φ)
	z: -5.831853859137424 * DEG2RAD // height (z)
});
```

## Variables

- `p`: radial distance ρ is the Euclidean distance from the z axis to the point P.
- `l`: The azimuth φ is the angle between the reference direction on the chosen plane and the line from the origin to the projection of P on the plane.
- `z`: The height z is the signed distance from the chosen plane to the point P. Same as the cartesian z value.
- `x`: Cartesian x value.
- `y`: Cartesian y value.
- `b`: Spherical θ as elevation.
- `i`: Spherical θ as inclination.
- `r`: Spherical radius.

## Methods

### `cart()`

Return [Cartesian Coordinates][2]

[2]: https://en.wikipedia.org/wiki/Cartesian_coordinate_system (xyz)

### `sph()`

Return [Spherical Coordinates][3]

[3]: https://en.wikipedia.org/wiki/Spherical_coordinate_system (lbr)

### `cyl()`

Return [Cylindrical Coordinates][4]

[4]: https://en.wikipedia.org/wiki/Cylindrical_coordinate_system (zlp)

### `ecl2equ(tilt)`

Ecliptic to Equatorial Coordinates; Pass obliquity of ecliptic (axial tilt)

### `equ2ecl(tilt)`

Equatorial to Ecliptic Coordinates; Pass obliquity of ecliptic (axial tilt)

### `gal2equ()`

Galactic to Equatorial Coordinates

### `equ2gal()`

Equatorial to Galactic Coordinates