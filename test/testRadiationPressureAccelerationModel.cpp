/*
 * Copyright (c) 2014-2018 Kartik Kumar (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2016 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2016 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <iostream>

#include <cmath>
#include <vector>

#include <catch2/catch.hpp>

#include "astro/radiationPressureAccelerationModel.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector;

TEST_CASE( "Compute radiation pressure for complete absorption at 1 AU", "[radiation_pressure]" )
{
    // Test data obtained from Montenbruck & Gill (2000).

    // Set expected radiation pressure for complete absorption using energy flux at the Earth (1 AU)
    // [N m^-2].
    const Real expectedRadiationPressure = 4.560e-6;

    // Set tolerance for error between expected value and computed value.
    const Real tolerance = 1.0e-4;

    // Set energy flux at 1 AU (McCarthy, 1996).
    const Real energyFlux = 1367.0;

    // Check if computed radiation pressure matches expected value.
    REQUIRE( computeAbsorptionRadiationPressure( energyFlux )
             == Approx( expectedRadiationPressure ).epsilon( tolerance ) );
}

TEST_CASE( "Compute radiation pressure at Mercury by using 1 AU as reference",
           "[radiation_pressure]")
{
    // Test data obtained from Wikipedia (2018).

    // Set Mercury distance.
    const Real distance = 0.2;

    // Set expected radiation pressure for perfect reflectance at 0.2 AU [N m^-2].
    const Real expectedRadiationPressure = 227.0e-6;

    // Set tolerance for error between expected value and computed value.
    const Real tolerance = 1.0e-15;

    // Set reference distance to 1 AU.
    const Real referenceDistance = 1.0;

    // Set reference radiation pressure at 1 AU.
    const Real referenceRadiationPressure = 9.08e-6;

    // Check if computed radiation pressure matches expected value.
    REQUIRE( computeRadiationPressure( referenceRadiationPressure, referenceDistance, distance )
             ==  Approx( expectedRadiationPressure ).epsilon( tolerance ) );
}

TEST_CASE( "Compute radiation pressure acceleration for arbitrary case",
           "[radiation_pressure, acceleration, models]" )
{
    // // Set expected radiation pressure acceleration vector [m/s^2].
    // Vector expectedAcceleration( 3 );
    // expectedAcceleration[ 0 ] = -2.964e-06;
    // expectedAcceleration[ 1 ] = 0.0;
    // expectedAcceleration[ 2 ] = 0.0;

    // // Set 1 AU in metres [m].
    // const Real astronomicalUnitInMeters = 1.49598e11;

//     // Set tolerance for error between expected value and computed value.
//     const Real tolerance = 1.0e-15;

//     // Solar Radiation Pressure at 1 AU [N m^-2].
//     const Real solarRadiationPressure = 4.56e-6;

//     // Set radiation pressure coefficient.
//     const Real radiationPressureCoefficient = 1.0 + 0.3;

//     // Set absorbing area [m^2].
//     const Real area = 2.0;

//     // Set mass [kg].
//     const Real mass = 4.0;

//     // Set unit vector pointing from S/C to the Sun.
//     Vector vectorToSource( 3 );
//     vectorToSource[ 0 ] = astronomicalUnitInMeters;
//     vectorToSource[ 1 ] = 0.0;
//     vectorToSource[ 2 ] = 0.0;

//     // Compute the unit vector to the Sun.
//     const Real normVectorToSource = std::sqrt( vectorToSource[ 0 ] * vectorToSource[ 0 ]
//                                                + vectorToSource[ 1 ] * vectorToSource[ 1 ]
//                                                + vectorToSource[ 2 ] * vectorToSource[ 2 ] );

//     const Real squaredNormVectorToSource = normVectorToSource * normVectorToSource;

//     Vector unitVectorToSource( 3 );
//     unitVectorToSource[ 0 ] = vectorToSource[ 0 ] / normVectorToSource;
//     unitVectorToSource[ 1 ] = vectorToSource[ 1 ] / normVectorToSource;
//     unitVectorToSource[ 2 ] = vectorToSource[ 2 ] / normVectorToSource;

//     // Set radiation pressure at target [N/m^2].
//     const double radiationPressureAtTarget =  solarRadiationPressure
//                                               * astronomicalUnitInMeters * astronomicalUnitInMeters
//                                               / squaredNormVectorToSource;

//      // Compute the radiation pressure [m s^-2].
//     const Vector computedAcceleration
//         = computeSolarRadiationPressureAcceleration(
//             radiationPressureAtTarget,      // Solar radiation pressure                    [N m^-2]
//             radiationPressureCoefficient,   // Radiation pressure coefficient              [adim]
//             unitVectorToSource,             // Unit vector pointing from S/C to sun (3x1)  [adim]
//             area,                           // Absorbing area of S/C                       [m^2]
//             mass );                         // Mass of the S/C                             [kg]

//     // Check if computed acceleration matches expected values.
//     REQUIRE( computedAcceleration[ 0 ]
//              == Approx( expectedAcceleration[ 0 ] ).epsilon( tolerance ) );
//     REQUIRE( computedAcceleration[ 1 ]
//              == Approx( expectedAcceleration[ 1 ] ).epsilon( tolerance ) );
//     REQUIRE( computedAcceleration[ 2 ]
//              == Approx( expectedAcceleration[ 2 ] ).epsilon( tolerance ) );
}

} // namespace tests
} // namespace astro

/*!
 * References
 *  Montenbruck, O, and Gill, E. Satellite orbits: models, methods and applications. Springer
 *      Science & Business Media, 2012.
 *  Wikipedia. Radiation pressure. https://en.wikipedia.org/wiki/Radiation_pressure, 2018.
 */
