(function () {

	// most values verified via wikipedia
	// initial values are from vsop87_mer(0)

	var orbit = {
		// state vectors
		rx: 1.390712873654179, // au
		ry: -0.013411960276849655, // au
		rz: -0.034462732955625755, // au
		vx: 0.0006714786862244502, // au
		vy: 0.015187272121968753, // au
		vz: 0.0003016547612144848, // au
		// vsop orbital elements
		a: 1.5236773948307722, // au
		L: 6.203874882756421, // rad
		h: -0.037808407518263365, // ?
		k: 0.08531441689241746, // ?
		q: 0.01047042574, // ?
		p: 0.01228449307, // ?
		// main orbital elements
		e: 0.09331680132087833, // unitless
		i: 0.032283817314135946, // rad
		O: 0.8649518974751991, // rad
		w: -1.2821077484807832, // rad
		W: -0.4171558510055842, // rad
		n: 0.00914622550420718, // rad
		m: 0.40725574492289923, // rad
		M: 0.3378454265824189, // rad
		// additional parameters
		c: 1.3814926941002357, // au
		C: 1.6658620955613088, // au
		P: 686.9703031364554 // d
	};

	var orbits = [];
	var orbitals = [];

	orbits.push({
		name: 'VSOP',
		a: orbit.a, // au
		L: orbit.L, // rad
		h: orbit.h, // rad
		k: orbit.k, // rad
		q: orbit.q, // rad
		p: orbit.p  // rad
	});

	orbits.push({
		name: 'State',
		rx: orbit.rx, // au
		ry: orbit.ry, // au
		rz: orbit.rz, // au
		vx: orbit.vx, // au
		vy: orbit.vy, // au
		vz: orbit.vz  // au
	});

	orbits.push({
		name: 'Planet V1',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		L: orbit.L  // rad
	});

	orbits.push({
		name: 'Planet V2',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		w: orbit.w, // rad
		W: orbit.W, // rad
		L: orbit.L  // rad
	});

	orbits.push({
		name: 'Planet V3',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'Planet V4',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		w: orbit.w, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'Planet V5',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		M: orbit.M, // rad
		L: orbit.L  // rad
	});

	orbits.push({
		name: 'Planet V6',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		w: orbit.w, // rad
		M: orbit.M, // rad
		L: orbit.L  // rad
	});

	orbits.push({
		name: 'Planet V7',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		w: orbit.w, // rad
		L: orbit.L  // rad
	});

	orbits.push({
		name: 'Planet V8',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		w: orbit.w, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'ISS V1',
		c: orbit.c, // au
		C: orbit.C, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'ISS V2',
		c: orbit.c, // au
		e: orbit.e, // rad
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'ISS V3',
		c: orbit.c, // au
		a: orbit.a, // rad
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'ISS V4',
		C: orbit.C, // au
		e: orbit.e, // rad
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'ISS V5',
		C: orbit.C, // au
		a: orbit.a, // rad
		i: orbit.i, // rad
		O: orbit.O, // rad
		W: orbit.W, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'Asteroid',
		e: orbit.e, // unitless
		a: orbit.a, // au
		i: orbit.i, // rad
		O: orbit.O, // rad
		L: orbit.L, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'TLE V1',
		e: orbit.e, // unitless
		i: orbit.i, // rad
		O: orbit.O, // rad
		w: orbit.w, // rad
		n: orbit.n, // rad
		M: orbit.M  // rad
	});

	orbits.push({
		name: 'TLE V2',
		e: orbit.e, // unitless
		i: orbit.i, // rad
		O: orbit.O, // rad
		w: orbit.w, // rad
		P: orbit.P, // days
		M: orbit.M  // rad
	});


	QUnit.module( 'Orbitals JD', {
			beforeEach: function () {

				// change mass system
				Orbit.GMP = GMJD;

				// reset array
				orbitals.length = 0;

				// create orbits of all variations
				for (var i = 0; i < orbits.length; i ++) {
					orbitals.push(new Orbit(orbits[i]));
				}

				// create exact copy of orbital
				var clone = orbitals[0].clone();
				clone.name = 'Cloned';
				orbitals.push(clone);

				// rotate orbit by one full cycle
				var cycled = clone.clone(orbit.P);
				cycled.name = 'Cycled';
				orbitals.push(cycled);

			}
		});

		QUnit.test( 'Inclination (i)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].i() * RAD2DEG, 1.8497264786713628, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Eccentricity (e)', function( assert )
		{
			var EPSILON = 1e-7; // local epsilon
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].e(), 0.09331680132087901, EPSILON, orbitals[i].name);
				assert.close(orbitals[i].e3().x, 0.0852789, EPSILON, orbitals[i].name + ' (x)');
				assert.close(orbitals[i].e3().y, -0.0377781, EPSILON, orbitals[i].name + ' (y)');
				assert.close(orbitals[i].e3().z, -0.0028874, EPSILON, orbitals[i].name + ' (z)');
			}
		});

		QUnit.test( 'Semi-Major Axis (a)', function( assert )
		{
			var EPSILON = 1e-5; // local epsilon
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].a(), 1.523679, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Semi-Minor Axis (b)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].b(), 1.5170287783679177, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Semi-Latus Rectum (l)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].l(), 1.5104091733618323, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Mean Longitude (L)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].L() * RAD2DEG, 355.4558474091613, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Longitude of the ascending node (O)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].O() * RAD2DEG, 49.558093207161185, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Argument of Perihelion (w)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].w() * RAD2DEG, -73.45936286896935, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Longitude of Periapsis (W)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].W() * RAD2DEG, -23.901269661808165, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Mean Motion (n)', function( assert )
		{
			var EPSILON = 1e-6; // local epsilon
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].n(), 0.009146944846617423, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Mean Anomaly (M)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].M() * RAD2DEG, 19.357117070969462, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Orbit Period (P)', function( assert )
		{
			var EPSILON = 1e-10; // local epsilon
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].P(), 686.9703031364554, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Angular Momentum (A)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].A(), 0.021141156875379624, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Periapsis (c)', function( assert )
		{
			var EPSILON = 1e-4; // local epsilon
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].c(), 1.3814, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Apoapsis (C)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].C(), 1.66586209556131, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Position state (r)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].r().x, orbit.rx, EPSILON, orbitals[i].name + ' position (x)');
				assert.close(orbitals[i].r().y, orbit.ry, EPSILON, orbitals[i].name + ' position (y)');
				assert.close(orbitals[i].r().z, orbit.rz, EPSILON, orbitals[i].name + ' position (z)');
			}
		});

		QUnit.test( 'Velocity state (v)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].v().x, orbit.vx, EPSILON, orbitals[i].name + ' velocity (X)');
				assert.close(orbitals[i].v().y, orbit.vy, EPSILON, orbitals[i].name + ' velocity (Y)');
				assert.close(orbitals[i].v().z, orbit.vz, EPSILON, orbitals[i].name + ' velocity (Z)');
			}
		});

		QUnit.test( 'True Anomaly (m)', function( assert )
		{
			var EPSILON = 1e-9; // local epsilon
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].m() * RAD2DEG, 23.334035367007363, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Radial Velocity (B)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].B(), 0.0005173553813483017, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'Eccentric Anomaly (E)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].E(), 0.371741701957453, EPSILON, orbitals[i].name);
			}
		});

		QUnit.test( 'VSOP parameters (khqp)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].k(), orbit.k, EPSILON, orbitals[i].name + ' VSOP k');
				assert.close(orbitals[i].h(), orbit.h, EPSILON, orbitals[i].name + ' VSOP h');
				assert.close(orbitals[i].q(), orbit.q, EPSILON, orbitals[i].name + ' VSOP q');
				assert.close(orbitals[i].p(), orbit.p, EPSILON, orbitals[i].name + ' VSOP p');
			}
		});

})();
