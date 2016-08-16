(function () {

	// checked via web converter
	// and via vsop87[ab]_jup(0)

	// cartesian
	// preferred
	var ctrlcart = {
		x: 4.001168777115359,
		y: 2.9385853724645656,
		z: -0.10178505133708564
	};

	// spherical with elevation
	// this is the VSOP87 format
	var ctrlsphe = {
		r: 4.9653867444245865, // radius (r)
		l: 36.29468623590961 * DEG2RAD, // phi (φ)
		b: -1.1745809849948283 * DEG2RAD, // theta (θ)
	};

	// spherical with inclination
	var ctrlsphi = {
		r: 4.9653797207164, // radius (r)
		l: 36.29468623590961 * DEG2RAD, // phi (φ)
		i: 91.17458536812916 * DEG2RAD, // theta (θ)
	};

	// cylindrical
	var ctrlcyl = {
		p: 4.9643363679575, // radial dist (ρ)
		l: 36.29474187359682 * DEG2RAD, // azimuth (φ)
		z: -5.831853859137424 * DEG2RAD // height (z)
	};

	// local epsilon
	var EPSILON = 1e-5;

	// test modules
	var tsts = [
		['Cartesian', ctrlcart ],
		['Cylindrical', ctrlcyl ],
		['Spherical (Elevation)', ctrlsphe ],
		['Spherical (Inclination)', ctrlsphi ],
	];

	// loop all test modules
	for (var i = 0; i < tsts.length; i ++) {
		// create new variable scope
		(function (i)
		{

			// local test variables
			var sph, cyl, cart;

			// define a new test module
			QUnit.module( tsts[i][0], {
				beforeEach: function () {
					// convert between coordinates
					sph = new Coord( tsts[i][1] ).sph();
					cyl = new Coord( tsts[i][1] ).cyl();
					cart = new Coord( tsts[i][1] ).cart();
				}
			}, function () {

				QUnit.test( 'To Spherical', function( assert )
				{

					assert.close(sph.r, ctrlsphi.r, EPSILON, 'Radius r (r)');
					assert.close(sph.b, ctrlsphe.b, EPSILON, 'Elevation θ (b)');
					assert.close(sph.i, ctrlsphi.i, EPSILON, 'Inclination θ (l)');
					assert.close(sph.l, ctrlsphe.l, EPSILON, 'Azimuth φ (l)');
				});

				QUnit.test( 'To Cartesian', function( assert )
				{
					assert.close(cart.x, ctrlcart.x, EPSILON, 'Position x (x)');
					assert.close(cart.y, ctrlcart.y, EPSILON, 'Position y (y)');
					assert.close(cart.z, ctrlcart.z, EPSILON, 'Position z (z)');
				});

				QUnit.test( 'To Cylindrical', function( assert )
				{
					assert.close(cyl.p, ctrlcyl.p, EPSILON, 'Radial Distance ρ (p)');
					assert.close(cyl.l, ctrlcyl.l, EPSILON, 'Azimuth φ (l)');
					assert.close(cyl.z, ctrlcyl.z, EPSILON, 'Elevation z (z)');
				});

			});

		})(i);

	}

})();