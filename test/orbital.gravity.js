(function () {

	QUnit.module( "Gravity", function () {

		// http://www.castor2.ca/05_OD/01_Gauss/14_Kepler/index.html

		var orbit = new Orbit({
			x: 5052.4587,
			y: 1056.2713,
			z: 5011.6366,
			X: 3.8589872,
			Y: 4.2763114,
			Z: -4.8070493,
			G: 398600.44
		});

		var cloned = orbit.clone(10);

		QUnit.test( 'Satellite Earth', function( assert )
		{

			assert.close(
				orbit.i() * RAD2DEG, 71.048202, 1e-6,
				"State vector to inclination (i)"
			);

			assert.close(
				orbit.e3().x, 0.011926, 1e-6,
				"State vector to eccentricity (x)"
			);
			assert.close(
				orbit.e3().y, 0.0031623, 1e-6,
				"State vector to eccentricity (y)"
			);
			assert.close(
				orbit.e3().z, 0.0101645, 1e-6,
				"State vector to eccentricity (z)"
			);
			assert.close(
				orbit.e(), 0.0159858, 1e-6,
				"State vector to eccentricity (e)"
			);

			assert.close(
				orbit.a(), 7310.8163295, 1e-6,
				"State vector to semi-major axis (a)"
			);
			assert.close(
				orbit.b(), 7309.88213535005, 1e-6,
				"State vector to semi-minor axis (b)"
			);
			assert.close(
				orbit.l(), 7308.948060569815, 1e-6,
				"State vector to semi-latus rectum (l)"
			);

			assert.close(
				orbit.L() * RAD2DEG, 344.01120082626187, 1e-2,
				"State vector to mean longitude (L)"
			);

			assert.close(
				orbit.O() * RAD2DEG, 211.2837710, 1e-6,
				"State vector to right ascending node (O)"
			);

			assert.close(
				orbit.w() * RAD2DEG, 137.75619, 1e-4,
				"State vector to argument of perigee (w)"
			);

			assert.close(
				orbit.T(), -6134.095906842509, 1e-4,
				"State vector to time of pericenter passage (T)"
			);

			assert.close(
				orbit.W() * RAD2DEG, -10.960123956443512, 1e-2,
				"State vector to longitude of pericenter (W)"
			);

			assert.close(
				orbit.n() * RAD2DEG, 0.05786856452419324, 1e-2,
				"State vector to mean motion (n)"
			);

			assert.close(
				orbit.M() * RAD2DEG, 354.9723980687702, 1e-2,
				"State vector to mean anomaly (M)"
			);

			assert.close(
				orbit.P(), 6220.9941262, 1e-6,
				"State vector to period (P)"
			);

			assert.close(
				orbit.A(), 53975.4565787, 1e-6,
				"State vector to orbital momentum (A)"
			);

			assert.close(
				orbit.m() * RAD2DEG, 354.80860, 1e-2,
				"State vector to true anomaly (m)"
			);

			assert.close(
				orbit.B(), -0.0106840, 1e-6,
				"State vector to radial velocity (B)"
			);

			assert.close(
				orbit.E() * RAD2DEG , 354.8908317970553, 1e-2,
				"State vector to eccentric anomaly (E)"
			);

			assert.close(
				orbit.C(), 7427.686218831075, 1e-2,
				"State vector to distance of pericenter (C)"
			);

			assert.close(
				orbit.C(), 7427.686218831075, 1e-2,
				"State vector to distance of apocenter (c)"
			);

			assert.close(
				orbit.r().x, 5052.4587, 1e-2,
				"State vector recalc position (x)"
			);
			assert.close(
				orbit.r().y, 1056.2713, 1e-2,
				"State vector recalc position (y)"
			);
			assert.close(
				orbit.r().z, 5011.6366, 1e-2,
				"State vector recalc position (z)"
			);

			assert.close(
				orbit.v().x, 3.8589872, 1e-2,
				"State vector recalc velocity (x)"
			);
			assert.close(
				orbit.v().y, 4.2763114, 1e-2,
				"State vector recalc velocity (y)"
			);
			assert.close(
				orbit.v().z, -4.8070493, 1e-2,
				"State vector recalc velocity (z)"
			);

		});


		QUnit.test( 'Satellite Cloned', function( assert )
		{

			assert.close(
				cloned.i() * RAD2DEG, 71.048202, 1e-6,
				"State vector to inclination (i)"
			);

			assert.close(
				cloned.e3().x, 0.011926, 1e-6,
				"State vector to eccentricity (x)"
			);
			assert.close(
				cloned.e3().y, 0.0031623, 1e-6,
				"State vector to eccentricity (y)"
			);
			assert.close(
				cloned.e3().z, 0.0101645, 1e-6,
				"State vector to eccentricity (z)"
			);
			assert.close(
				cloned.e(), 0.0159858, 1e-6,
				"State vector to eccentricity (e)"
			);

			assert.close(
				cloned.a(), 7310.8163295, 1e-6,
				"State vector to semi-major axis (a)"
			);
			assert.close(
				cloned.b(), 7309.88213535005, 1e-6,
				"State vector to semi-minor axis (b)"
			);
			assert.close(
				cloned.l(), 7308.948060569815, 1e-6,
				"State vector to semi-latus rectum (l)"
			);

			assert.close(
				cloned.L() * RAD2DEG, 344.58988647150386, 1e-2,
				"State vector to mean longitude (L)"
			);

			assert.close(
				cloned.O() * RAD2DEG, 211.2837710, 1e-6,
				"State vector to right ascending node (O)"
			);

			assert.close(
				cloned.w() * RAD2DEG, 137.75619, 1e-4,
				"State vector to argument of perigee (w)"
			);

			assert.close(
				cloned.T(), -6134.095906842509, 1e-4,
				"State vector to time of pericenter passage (T)"
			);

			assert.close(
				cloned.W() * RAD2DEG, -10.960123956443512, 1e-2,
				"State vector to longitude of pericenter (W)"
			);

			assert.close(
				cloned.n() * RAD2DEG, 0.05786856452419324, 1e-2,
				"State vector to mean motion (n)"
			);

			assert.close(
				cloned.M() * RAD2DEG, 355.55001042794737, 1e-2,
				"State vector to mean anomaly (M)"
			);

			assert.close(
				cloned.P(), 6220.9941262, 1e-6,
				"State vector to period (P)"
			);

			assert.close(
				cloned.A(), 53975.4565787, 1e-6,
				"State vector to orbital momentum (A)"
			);

			assert.close(
				cloned.m() * RAD2DEG, 355.4049928137755, 1e-2,
				"State vector to true anomaly (m)"
			);

			assert.close(
				cloned.B(), -0.0106840, 1e-6,
				"State vector to radial velocity (B)"
			);

			assert.close(
				cloned.E() * RAD2DEG , 355.47779398375076, 1e-2,
				"State vector to eccentric anomaly (E)"
			);

			assert.close(
				cloned.C(), 7427.686218831075, 1e-2,
				"State vector to distance of pericenter (C)"
			);

			assert.close(
				cloned.C(), 7427.686218831075, 1e-2,
				"State vector to distance of apocenter (c)"
			);

			assert.close(
				cloned.r().x, 5090.7774707576045, 1e-2,
				"State vector recalc position (x)"
			);
			assert.close(
				cloned.r().y, 1098.9771183585046, 1e-2,
				"State vector recalc position (y)"
			);
			assert.close(
				cloned.r().z, 4963.29873666537, 1e-2,
				"State vector recalc position (z)"
			);

			assert.close(
				cloned.v().x, 3.8046961084148414, 1e-2,
				"State vector recalc velocity (x)"
			);
			assert.close(
				cloned.v().y, 4.264773647122109, 1e-2,
				"State vector recalc velocity (y)"
			);
			assert.close(
				cloned.v().z, -4.86043482005546, 1e-2,
				"State vector recalc velocity (z)"
			);

		});

	});

})();

