/*    
 * Copyright (c) 2014 K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch.hpp>

#include <SAM/twoBodyMethods.hpp>

namespace sam
{
namespace tests
{

typedef double Real;

TEST_CASE( "Convert semi-major axis to mean motion", "[semi-major-axis-to-mean-motion]" )
{
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.

    // Set satellite mass [kg].
    const Real satelliteMass = 1.0e3;

    // Set gravitational parameter of Earth [m^3 s^-2].
    const Real earthGravitationalParameter = 6.67259e-11 * 5.9736e24;

    // Set distance between Earth center and satellite [m].
    const Real distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Set expected mean motion [rad/s].
    const Real expectedMeanMotion = 7.2921e-5;

    // Compute mean motion.
    const Real meanMotion = computeKeplerMeanMotion(
        distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Check if computed mean motion matches expected value.
    REQUIRE( meanMotion == Approx( expectedMeanMotion ).epsilon( 1.0e-7 ) );
}

TEST_CASE( "Convert semi-major axis to orbital period", "[semi-major-axis-to-orbital-period]" )
{
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.

    // Set satellite mass [kg].
    const Real satelliteMass = 1.0e3;

    // Set gravitational parameter of Earth [m^3 s^-2].
    const Real earthGravitationalParameter = 6.67259e-11 * 5.9736e24;

    // Set distance between Earth center and satellite [m].
    const Real distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Set expected orbital period [s].
    const Real expectedOrbitalPeriod = 86164.09054;

    // Compute orbital period of satellite.
    const Real orbitalPeriod = computeKeplerOrbitalPeriod(
        distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Check if computed orbital period matches expected value.
    REQUIRE( orbitalPeriod == Approx( expectedOrbitalPeriod ).epsilon( 1.0e-5 ) );
}

TEST_CASE( "Compute circular velocity", "[circular-velocity]" )
{
    SECTION( "Test zero-semi-major-axis error" )
    {
        REQUIRE_THROWS( computeCircularVelocity( 0.0, 0.0 ) );
    }
    
    SECTION( "Test orbits around the Earth" )
    {
        // Reference data obtained from Wertz (2001).

        // Set Earth equatorial radius [m].
        const Real earthRadius = 6378136.0;

        // Set Earth's gravitational parameter [m^3 s^-2].
        const Real earthGravitationalParameter = 3.98600441e14;

        // Set altitudes [km].
        const Real altitudes[ 5 ] = { 0.0, 200.0, 500.0, 1000.0, 35786.0 };

        // Set expected circular velocities [km/s].
        const Real expectedCircularVelocities[ 5 ] = { 7.905, 7.784, 7.613, 7.350, 3.075 };

        for ( unsigned int i = 0; i < 5; i++ )
        {
            // Compute circular velocity [m/s].
            const Real computedCircularVelocity 
                = computeCircularVelocity( 
                    earthRadius + altitudes[ i ] * 1.0e3, earthGravitationalParameter ); 

            //! Check if computed circular velocity matches expected value.
            REQUIRE( ( computedCircularVelocity / 1.0e3 )
                     == Approx( expectedCircularVelocities[ i ] ).epsilon( 1.0e-4 ) );
        }
    }

    SECTION( "Test orbit around Mars" )
    {
        // Reference data obtained from Wikipedia (2014a, 2014b).

        // Set Mars mean radius [km].
        const Real marsRadius = 3389.5;

        // Set Mars' gravitational parameter [km^3 s^-2].
        const Real marsGravitationalParameter = 42828.0;

        // Set radius of orbit [km].
        const Real radius = marsRadius + 200.0;

        // Set expected circular velocity [km/s].
        // This value was computed "by hand" using Julia.
        const Real expectedCircularVelocity = 3.454195532696839;

        const Real computedCircularVelocity 
            = computeCircularVelocity( radius, marsGravitationalParameter ); 

        //! Check if computed circular velocity matches expected value.
        REQUIRE( computedCircularVelocity == expectedCircularVelocity );
    }
}

} // namespace tests
} // namespace sam

/*!
 * References
 *  Vallado, D. A., McClain, W. D. Fundamentals of astrodynamics and applications, 2nd Edition,
 *   Kluwer Academic Publishers, The Netherlands, 2004. 
 *  Wertz, J.R. Mission Geometry: Orbit and Constellation Design and Management, Mircocosm Press,
 *   El Segundo, CA, 2001.
 *  Wikipedia. Mars, http://en.wikipedia.org/wiki/Mars, last updated: 2 Oct 2014, last accessed:
 *   7 Oct 2014a.
 *  Wikipedia. Standard gravitational parameter, 
 *   http://en.wikipedia.org/wiki/Standard_gravitational_parameter, last updated: 2 Oct 2014,
 *   last accessed: 7 Oct 2014b.
 */
