(function () {
	'use strict';

	// local epsilon
	var EPSILON = 1e-5;

	var equ = new Coord({ l: 2.65511899, b: 0.20805518, r: 1 });
	var gal = new Coord({ l: 3.95341300, b: 0.85419900, r: 1 });

	var galA = new Coord({ l: 0.0000000000 * DEG2RAD, b: 0.0000000000 * DEG2RAD, r: 1 });
	var equA = new Coord({ l: 266.40506655 * DEG2RAD, b: -28.93616241 * DEG2RAD, r: 1 });
	var equB = new Coord({ l: 266.41700000 * DEG2RAD, b: -29.00800000 * DEG2RAD, r: 1 });
	var galB = new Coord({ l: 359.94430560 * DEG2RAD, b: -0.046194444 * DEG2RAD, r: 1 });


	var galequ, galequA, galequB;
	QUnit.module( 'Galactic', {

		beforeEach: function () {
			galequ = gal.gal2equ();
			galequA = galA.gal2equ();
			galequB = galB.gal2equ();
		}
	}, function () {

		QUnit.test( 'To Equatorial', function( assert )
		{
			assert.close(galequ.l, equ.l, EPSILON, 'Right Ascension α (l)');
			assert.close(galequ.b, equ.b, EPSILON, 'Declination δ (b)');
			assert.close(galequA.l, equA.l, EPSILON, 'Right Ascension α (l)');
			assert.close(galequA.b, equA.b, EPSILON, 'Declination δ (b)');
			assert.close(galequB.l, equB.l, EPSILON, 'Right Ascension α (l)');
			assert.close(galequB.b, equB.b, EPSILON, 'Declination δ (b)');
		});

	});

	var equgal, equgalA, equgalB;
	QUnit.module( 'Equatorial', {

		beforeEach: function () {
			equgal = equ.equ2gal();
			equgalA = equA.equ2gal();
			equgalB = equB.equ2gal();
		}
	}, function () {

		QUnit.test( 'To Galactic', function( assert )
		{
			assert.close(equgal.l, gal.l, EPSILON, 'Right Ascension α (l)');
			assert.close(equgal.b, gal.b, EPSILON, 'Declination δ (b)');
			assert.close(equgalA.l, galA.l, EPSILON, 'Right Ascension α (l)');
			assert.close(equgalA.b, galA.b, EPSILON, 'Declination δ (b)');
			assert.close(equgalB.l, galB.l, EPSILON, 'Right Ascension α (l)');
			assert.close(equgalB.b, galB.b, EPSILON, 'Declination δ (b)');
		});

	});

})();