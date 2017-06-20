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


	QUnit.module( 'Epoch', {
			beforeEach: function () {

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

		QUnit.test( 'Position state +P (r)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].r(orbit.P).x, orbit.rx, EPSILON, orbitals[i].name + ' position (x)');
				assert.close(orbitals[i].r(orbit.P).y, orbit.ry, EPSILON, orbitals[i].name + ' position (y)');
				assert.close(orbitals[i].r(orbit.P).z, orbit.rz, EPSILON, orbitals[i].name + ' position (z)');
			}
		});

		QUnit.test( 'Velocity state +P (v)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].v(orbit.P).x, orbit.vx, EPSILON, orbitals[i].name + ' velocity (X)');
				assert.close(orbitals[i].v(orbit.P).y, orbit.vy, EPSILON, orbitals[i].name + ' velocity (Y)');
				assert.close(orbitals[i].v(orbit.P).z, orbit.vz, EPSILON, orbitals[i].name + ' velocity (Z)');
			}
		});

		QUnit.test( 'Position state +42JD (r)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].r(42).x, 1.2878811331585023, EPSILON, orbitals[i].name + ' position (x)');
				assert.close(orbitals[i].r(42).y, 0.6058984859445274, EPSILON, orbitals[i].name + ' position (y)');
				assert.close(orbitals[i].r(42).z, -0.01896131454009689, EPSILON, orbitals[i].name + ' position (z)');
			}
		});

		QUnit.test( 'Velocity state +42JD (v)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].v(42).x, -0.005421667544623961, EPSILON, orbitals[i].name + ' velocity (X)');
				assert.close(orbitals[i].v(42).y, 0.013856217072479677, EPSILON, orbitals[i].name + ' velocity (Y)');
				assert.close(orbitals[i].v(42).z, 0.0004235313959266824, EPSILON, orbitals[i].name + ' velocity (Z)');
			}
		});

		QUnit.test( 'Position state +42JD -P (r)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].r(42-orbit.P).x, 1.2878811331585023, EPSILON, orbitals[i].name + ' position (x)');
				assert.close(orbitals[i].r(42-orbit.P).y, 0.6058984859445274, EPSILON, orbitals[i].name + ' position (y)');
				assert.close(orbitals[i].r(42-orbit.P).z, -0.01896131454009689, EPSILON, orbitals[i].name + ' position (z)');
			}
		});

		QUnit.test( 'Velocity state +42JD -P (v)', function( assert )
		{
			for (var i = 0; i < orbitals.length; i++) {
				assert.close(orbitals[i].v(42-orbit.P).x, -0.005421667544623961, EPSILON, orbitals[i].name + ' velocity (X)');
				assert.close(orbitals[i].v(42-orbit.P).y, 0.013856217072479677, EPSILON, orbitals[i].name + ' velocity (Y)');
				assert.close(orbitals[i].v(42-orbit.P).z, 0.0004235313959266824, EPSILON, orbitals[i].name + ' velocity (Z)');
			}
		});


})();
