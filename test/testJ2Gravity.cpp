/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Abhishek Agrawal, Delft University of Technology
 *                                           (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <limits>
#include <vector>

#include <catch.hpp>

#include "Astro/j2Gravity.hpp"
#include "Astro/twoBodyMethods.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector;

TEST_CASE( "First order, orbit-averaged Keplerian element rate of change due to J2 perturbation",
           "[first-order-orbit-averaged-j2]" )
{
    const Real pi = 3.14159265358979323846;
    // Earth's equatorial radius in [km].
    const Real earthEquatorialRadius = 6378.13649;
    // Earth's gravitational parameter [km^3 s^-2].
    const Real earthGravitationalParameter = 398600.4418;

    SECTION( "Test for Shuttle (LEO) orbit" )
    {
        Vector keplerianElements( 6 );
        keplerianElements[ 0 ] = 6700.0;
        keplerianElements[ 1 ] = 0.0;
        keplerianElements[ 2 ] = 28.0 * pi / 180.0 ;
        keplerianElements[ 3 ] = 0.0;
        keplerianElements[ 4 ] = 0.0;
        keplerianElements[ 5 ] = 0.0;

        // Mean motion or orbiting object in radians per second.
        const Real semiMajorAxis = keplerianElements[ 0 ];
        const Real massOrbitingBody = 0.0;
        const Real meanMotion = astro::computeKeplerMeanMotion( semiMajorAxis,
                                                                earthGravitationalParameter,
                                                                massOrbitingBody );
        // Mean motion in degrees per day.
        const Real meanMotionDegreesPerDay = ( meanMotion * 86400.0 ) * ( 180.0 / pi );

        Real longitudeAscendingNodeDot = 0.0;
        Real argumentOfPeriapsisDot = 0.0;

        const Real expectedLongitudeAscendingNodeDot = -7.35;
        const Real expectedArgumentOfPeriapsisDot = 12.05;

        computeFirstOrderAveragedEffectJ2Perturbation(
            keplerianElements,
            meanMotionDegreesPerDay,
            earthEquatorialRadius,
            longitudeAscendingNodeDot,
            argumentOfPeriapsisDot );

        REQUIRE( expectedLongitudeAscendingNodeDot
                    == Approx( longitudeAscendingNodeDot ).epsilon( 1.0e-2 ) );
        REQUIRE( expectedArgumentOfPeriapsisDot
                    == Approx( argumentOfPeriapsisDot ).epsilon( 1.0e-2 ) );
    }

    SECTION( "Test for GPS (HEO) orbit" )
    {
        Vector keplerianElements( 6 );
        keplerianElements[ 0 ] = 26600.0;
        keplerianElements[ 1 ] = 0.0;
        keplerianElements[ 2 ] = 60.0 * pi / 180.0 ;
        keplerianElements[ 3 ] = 0.0;
        keplerianElements[ 4 ] = 0.0;
        keplerianElements[ 5 ] = 0.0;

        // Mean motion or orbiting object in radians per second.
        const Real semiMajorAxis = keplerianElements[ 0 ];
        const Real massOrbitingBody = 0.0;
        const Real meanMotion = astro::computeKeplerMeanMotion( semiMajorAxis,
                                                                earthGravitationalParameter,
                                                                massOrbitingBody );
        // Mean motion in degrees per day.
        const Real meanMotionDegreesPerDay = ( meanMotion * 86400.0 ) * ( 180.0 / pi );

        Real longitudeAscendingNodeDot = 0.0;
        Real argumentOfPeriapsisDot = 0.0;

        const Real expectedLongitudeAscendingNodeDot = -0.033;
        const Real expectedArgumentOfPeriapsisDot = 0.008;

        computeFirstOrderAveragedEffectJ2Perturbation(
            keplerianElements,
            meanMotionDegreesPerDay,
            earthEquatorialRadius,
            longitudeAscendingNodeDot,
            argumentOfPeriapsisDot );

        REQUIRE( expectedLongitudeAscendingNodeDot
                    == Approx( longitudeAscendingNodeDot ).epsilon( 1.0e-3 ) );
        REQUIRE( expectedArgumentOfPeriapsisDot
                    == Approx( argumentOfPeriapsisDot ).epsilon( 1.0e-3 ) );
    }

    SECTION( "Test for Geo-Stationary (GEO) orbit" )
    {
        Vector keplerianElements( 6 );
        keplerianElements[ 0 ] = 42160.0;
        keplerianElements[ 1 ] = 0.0;
        keplerianElements[ 2 ] = 0.0;
        keplerianElements[ 3 ] = 0.0;
        keplerianElements[ 4 ] = 0.0;
        keplerianElements[ 5 ] = 0.0;

        // Mean motion or orbiting object in radians per second.
        const Real semiMajorAxis = keplerianElements[ 0 ];
        const Real massOrbitingBody = 0.0;
        const Real meanMotion = astro::computeKeplerMeanMotion( semiMajorAxis,
                                                                earthGravitationalParameter,
                                                                massOrbitingBody );
        // Mean motion in degrees per day.
        const Real meanMotionDegreesPerDay = ( meanMotion * 86400.0 ) * ( 180.0 / pi );

        Real longitudeAscendingNodeDot = 0.0;
        Real argumentOfPeriapsisDot = 0.0;

        const Real expectedLongitudeAscendingNodeDot = -0.013;
        const Real expectedArgumentOfPeriapsisDot = 0.026;

        computeFirstOrderAveragedEffectJ2Perturbation(
            keplerianElements,
            meanMotionDegreesPerDay,
            earthEquatorialRadius,
            longitudeAscendingNodeDot,
            argumentOfPeriapsisDot );

        REQUIRE( expectedLongitudeAscendingNodeDot
                    == Approx( longitudeAscendingNodeDot ).epsilon( 1.0e-3 ) );
        REQUIRE( expectedArgumentOfPeriapsisDot
                    == Approx( argumentOfPeriapsisDot ).epsilon( 1.0e-3 ) );
    }
}

} // namespace tests
} // namespace astro

/*!
 * References
 * Wertz, J.R., et al. Space mission analysis and design, Third Edition,
 * ISBN 1-881883-10-8, Space Technology Library.
 */
