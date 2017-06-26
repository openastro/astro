/*
 * Copyright (c) 2014-2015 Kartik Kumar, Dinamica Srl
 * Copyright (c) 2014-2015 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2015 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2015 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cmath>
#include <vector>

#include <catch.hpp>

#include "Astro/centralBodyAccelerationModel.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector;

TEST_CASE( "Compute acceleration vector for a GEO spacecraft using a custom made MATLAB script",
           "[central_gravity, acceleration, models]" )
{
    // Benchmark values for test case obtained using custom MATLAB script: https://goo.gl/8BQztN

    // Set expected acceleration vector [km s^-2].
    Vector expectedAcceleration( 3 );
    expectedAcceleration[ 0 ] = -2.242096133923724e-4;
    expectedAcceleration[ 1 ] =  0.0;
    expectedAcceleration[ 2 ] =  0.0;

     // Set tolerance = error between expected value and computed value.
    const Real tolerance = 1.0e-15;

    // Set value of gravitational parameter of central body [km^3 s^-2].
    // In this case the used value is that of the planet Earth (Wertz, 1999; page 132).
    const Real gravitationalParameter = 3.986005e5;

    // Set position vector of the spacecraft relative to the origin of the reference frame [km].
    // In the following example the Earth-fixed reference frame is used ((Wertz, 1999; page 96). 
    // Moreover, the spacecraft is assumed to be in positioned in Geostationary Orbit 
    // (GEO) and on the Greenwich Meridian. 
    Vector positionVector( 3 );
    positionVector[ 0 ] = 4.2164e4; 
    positionVector[ 1 ] = 0.0;
    positionVector[ 2 ] = 0.0;

    // Compute the acceleration vector [km s^-2].
    const Vector computedAcceleration = computeCentralBodyAcceleration( gravitationalParameter,
                                                                        positionVector );

    // Check if computed acceleration matches expected values.
    REQUIRE( computedAcceleration[ 0 ] 
             == Approx( expectedAcceleration[ 0 ] ).epsilon( tolerance ) );
    REQUIRE( computedAcceleration[ 1 ]
             == Approx( expectedAcceleration[ 1 ] ).epsilon( tolerance ) );
    REQUIRE( computedAcceleration[ 2 ]
             == Approx( expectedAcceleration[ 2 ] ).epsilon( tolerance ) );
}

TEST_CASE( "Compute acceleration vector for arbitrary case using Tudat library values",
           "[central_gravity, acceleration, models]" )
{
    // Test case computed using the data of the Tudat library calculated using the planet Mercury as 
    // a central body. More info can be found at the following link: http://goo.gl/MyXBJE.

    // Set expected acceleration vector [m s^-2].
    Vector expectedAcceleration( 3 );
    expectedAcceleration[ 0 ] = -6.174552714649318e-2;
    expectedAcceleration[ 1 ] =  3.024510782481964e-1;
    expectedAcceleration[ 2 ] = -1.228994266291893e-1;

     // Set tolerance = error between expected value and computed value.
    const Real tolerance = 1.0e-15;

    // Set value of gravitational parameter of central body [m^3 s^-2].
    const Real gravitationalParameter = 2.2032e13;

    // Set position vector of the body relative to the origin of the reference frame [m].
    Vector positionVector( 3 );
    positionVector[ 0 ] =  1513.3e3;
    positionVector[ 1 ] = -7412.67e3;
    positionVector[ 2 ] =  3012.1e3;

    // Compute the acceleration vector [m s^-2].
    const Vector computedAcceleration = computeCentralBodyAcceleration( gravitationalParameter,   
                                                                        positionVector);

    // Check if computed acceleration matches expected values.
    REQUIRE( computedAcceleration[ 0 ]
             == Approx( expectedAcceleration[ 0 ] ).epsilon( tolerance ) );
    REQUIRE( computedAcceleration[ 1 ]
             == Approx( expectedAcceleration[ 1 ] ).epsilon( tolerance ) );
    REQUIRE( computedAcceleration[ 2 ]
             == Approx( expectedAcceleration[ 2 ] ).epsilon( tolerance ) );
}

} // namespace tests
} // namespace astro

/*
 * Wertz, J.R. & Larson, W.J., Space mission analysis and design, Springer Netherlands
 */
