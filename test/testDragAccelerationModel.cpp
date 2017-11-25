/*
 * Copyright (c) 2015 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cmath>
#include <vector>

#include <catch.hpp>

#include "Astro/dragAccelerationModel.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector;

TEST_CASE( "Obtain drag acceleration: test 1", "[drag, acceleration, models]" )
{
    // Set expected drag acceleration vector [m/s].
    Vector expectedDragAcceleration( 3 );
    expectedDragAcceleration[ 0 ] = 0.107800109999944e-4;
    expectedDragAcceleration[ 1 ] = 0.0;
    expectedDragAcceleration[ 2 ] = 0.000154000157143e-4;

    // Set epsilon = error between expected value and computed value.
    const Real epsilon = 1.0e-10;

    // Set drag coefficient.
    const Real dragCoefficient = 2.2;

    // Set atmospheric density [kg/m^3].
    const Real atmosphericDensity = 2.0e-11;

    // Velocity vector in [m/s].
    Vector velocity( 3 );
    velocity[ 0 ] = 7000.0;
    velocity[ 1 ] = 0.0;
    velocity[ 2 ] = 10.0;

    // Set drag area [m^2].
    const Real dragArea = 5.0;

    // Set mass [kg].
    const Real mass = 500.0;

     //! Compute drag acceleration.
    const Vector dragAcceleration = computeDragAcceleration( dragCoefficient, 
                                                             atmosphericDensity, 
                                                             velocity,
                                                             dragArea,
                                                             mass ); 

    // Check if computed mean motion matches expected value.
    REQUIRE( std::fabs( dragAcceleration[ 0 ] - expectedDragAcceleration[ 0 ] ) <= epsilon );
    REQUIRE( std::fabs( dragAcceleration[ 1 ] - expectedDragAcceleration[ 1 ] ) <= epsilon );
    REQUIRE( std::fabs( dragAcceleration[ 2 ] - expectedDragAcceleration[ 2 ] ) <= epsilon );
}

TEST_CASE( "Obtain drag acceleration: test 2", "[obtain-drag-acceleration-2]" )
{
    // Reference: http://tudat.tudelft.nl/

    // Set expected drag acceleration vector [m/s].
    Vector expectedDragAcceleration( 3 );
    expectedDragAcceleration[ 0 ] = 0.0;
    expectedDragAcceleration[ 1 ] = 0.0;
    expectedDragAcceleration[ 2 ] =  267.4211815284975;

    // Set epsilon = error between expected value and computed value.
    const Real epsilon = 1.0e-10;

    // Set drag coefficient.
    const Real dragCoefficient = 1.1;

    // Set atmospheric density [kg/m3].
    const Real atmosphericDensity = 3.5e-5;

    // Velocity vector in [m/s].
    Vector velocity( 3 );
    velocity[ 0 ] = 0.0;
    velocity[ 1 ] = 0.0;
    velocity[ 2 ] = 3491.0;

    // Set drag area [m^2]
    const Real dragArea = 2.2;

    // Set mass [kg]
    const Real mass = 1.93;

     //! Compute drag acceleration.
    const Vector dragAcceleration = computeDragAcceleration( dragCoefficient, 
                                                             atmosphericDensity, 
                                                             velocity,
                                                             dragArea,
                                                             mass ); 

    // Check if computed mean motion matches expected value.
    REQUIRE( std::fabs(dragAcceleration[ 0 ] - expectedDragAcceleration[ 0 ]) <= epsilon );
    REQUIRE( std::fabs(dragAcceleration[ 1 ] - expectedDragAcceleration[ 1 ]) <= epsilon );
    REQUIRE( std::fabs(dragAcceleration[ 2 ] - expectedDragAcceleration[ 2 ]) <= epsilon );
}

} // namespace tests
} // namespace astro
