/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <catch.hpp>

#include <Eigen/Core>

#include "Astro/relativeMotion.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef Eigen::Matrix< Real, Eigen::Dynamic, 1 > Vector;
static const Real tolerance = 100.0 * std::numeric_limits< Real >::epsilon( );

TEST_CASE( "Test Clohessy-Wiltshire solution", "[clohessy-wiltshire,relative-motion]" )
{
    // Scenarios taken from chapter 3, Fehse (2003).

    const Real earthGravitationalParameter = 398600.14e9;       // m^3 s^-2
    const Real targetSemiMajorAxis = 7200.0e3;                  // m
    const Real targetMeanMotion
        = std::sqrt( earthGravitationalParameter
            / ( targetSemiMajorAxis * targetSemiMajorAxis * targetSemiMajorAxis ) );
    const Real finalTime = 100.0;                               // s

    SECTION( "Test free motion" )
    {
        Vector thrustAcceleration( 3 );
        thrustAcceleration[ 0 ] = 0.0;
        thrustAcceleration[ 1 ] = 0.0;
        thrustAcceleration[ 2 ] = 0.0;

        SECTION( "Test different altitude" )
        {
            Vector initialState( 6 );
            initialState[ 0 ] = 0.0;
            initialState[ 1 ] = 0.0;
            initialState[ 2 ] = 10.0;
            initialState[ 3 ] = 1.5 * targetMeanMotion * 10.0;
            initialState[ 4 ] = 0.0;
            initialState[ 5 ] = 0.0;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = 1.5 * targetMeanMotion * 10.0 * finalTime;
            expectedFinalState[ 1 ] = 0.0;
            expectedFinalState[ 2 ] = 10.0;
            expectedFinalState[ 3 ] = 1.5 * targetMeanMotion * 10.0;
            expectedFinalState[ 4 ] = 0.0;
            expectedFinalState[ 5 ] = 0.0;

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }

        SECTION( "Test different altitude with target's velocity" )
        {
            Vector initialState( 6 );
            initialState[ 0 ] = 0.0;
            initialState[ 1 ] = 0.0;
            initialState[ 2 ] = 10.0;
            initialState[ 3 ] = 0.0;
            initialState[ 4 ] = 0.0;
            initialState[ 5 ] = 0.0;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = 6.0 * 10.0 * ( targetMeanMotion * finalTime
                - std::sin( targetMeanMotion * finalTime ) );
            expectedFinalState[ 1 ] = 0.0;
            expectedFinalState[ 2 ]
                = 10.0 * ( 4.0 - 3.0 * std::cos( targetMeanMotion * finalTime ) );
            expectedFinalState[ 3 ] = 6.0 * 10.0 * targetMeanMotion
                * ( 1.0 - std::cos( targetMeanMotion * finalTime ) );
            expectedFinalState[ 4 ] = 0.0;
            expectedFinalState[ 5 ]
                = 3.0 * 10.0 * targetMeanMotion * std::sin( targetMeanMotion * finalTime );

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }

        SECTION( "Test out-of-plane free drift" )
        {
            Vector initialState( 6 );
            initialState[ 0 ] = 0.0;
            initialState[ 1 ] = 20.0;
            initialState[ 2 ] = 0.0;
            initialState[ 3 ] = 0.0;
            initialState[ 4 ] = 0.0;
            initialState[ 5 ] = 0.0;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = 0.0;
            expectedFinalState[ 1 ] = 20.0 * std::cos( targetMeanMotion * finalTime );
            expectedFinalState[ 2 ] = 0.0;
            expectedFinalState[ 3 ] = 0.0;
            expectedFinalState[ 4 ]
                = -20.0 * targetMeanMotion * std::sin( targetMeanMotion * finalTime );
            expectedFinalState[ 5 ] = 0.0;

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }

        SECTION( "Test along-track free drift" )
        {
            Vector initialState( 6 );
            initialState[ 0 ] = 5.0;
            initialState[ 1 ] = 0.0;
            initialState[ 2 ] = 0.0;
            initialState[ 3 ] = 0.0;
            initialState[ 4 ] = 0.0;
            initialState[ 5 ] = 0.0;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = 5.0;
            expectedFinalState[ 1 ] = 0.0;
            expectedFinalState[ 2 ] = 0.0;
            expectedFinalState[ 3 ] = 0.0;
            expectedFinalState[ 4 ] = 0.0;
            expectedFinalState[ 5 ] = 0.0;

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }
    }

    SECTION( "Test constant-force motion applied at the target" )
    {
        Vector initialState( 6 );
        initialState[ 0 ] = 0.0;
        initialState[ 1 ] = 0.0;
        initialState[ 2 ] = 0.0;
        initialState[ 3 ] = 0.0;
        initialState[ 4 ] = 0.0;
        initialState[ 5 ] = 0.0;

        SECTION( "Test along-track force" )
        {
            Vector thrustAcceleration( 3 );
            thrustAcceleration[ 0 ] = 2.161;
            thrustAcceleration[ 1 ] = 0.0;
            thrustAcceleration[ 2 ] = 0.0;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ]
                = thrustAcceleration[ 0 ] / ( targetMeanMotion * targetMeanMotion )
                  * ( 4.0 * ( 1.0 - std::cos( targetMeanMotion * finalTime ) )
                      - 1.5 * targetMeanMotion * targetMeanMotion * finalTime * finalTime );
            expectedFinalState[ 1 ] = 0.0;
            expectedFinalState[ 2 ]
                = 2.0 * thrustAcceleration[ 0 ] / ( targetMeanMotion * targetMeanMotion )
                  * ( std::sin( targetMeanMotion * finalTime ) - targetMeanMotion * finalTime );
            expectedFinalState[ 3 ]
                = thrustAcceleration[ 0 ] / ( targetMeanMotion * targetMeanMotion )
                  * ( 4.0 * targetMeanMotion * std::sin( targetMeanMotion * finalTime )
                      - 3.0 * targetMeanMotion * targetMeanMotion * finalTime );
            expectedFinalState[ 4 ] = 0.0;
            expectedFinalState[ 5 ]
                = 2.0 * thrustAcceleration[ 0 ] / targetMeanMotion
                  * ( std::cos( targetMeanMotion * finalTime ) - 1.0 );

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }

        SECTION( "Test cross-track force" )
        {
            Vector thrustAcceleration( 3 );
            thrustAcceleration[ 0 ] = 0.0;
            thrustAcceleration[ 1 ] = -4.915;
            thrustAcceleration[ 2 ] = 0.0;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = 0.0;
            expectedFinalState[ 1 ]
                = thrustAcceleration[ 1 ] / ( targetMeanMotion * targetMeanMotion )
                  * ( 1.0 - std::cos( targetMeanMotion * finalTime ) );
            expectedFinalState[ 2 ] = 0.0;
            expectedFinalState[ 3 ] = 0.0;
            expectedFinalState[ 4 ]
                = thrustAcceleration[ 1 ] / targetMeanMotion
                  * std::sin( targetMeanMotion * finalTime );
            expectedFinalState[ 5 ] = 0.0;

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }

        SECTION( "Test radial force" )
        {
            Vector thrustAcceleration( 3 );
            thrustAcceleration[ 0 ] = 0.0;
            thrustAcceleration[ 1 ] = 0.0;
            thrustAcceleration[ 2 ] = 5.842;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ]
                = 2.0 * thrustAcceleration[ 2 ] / ( targetMeanMotion * targetMeanMotion )
                  * ( targetMeanMotion * finalTime - std::sin( targetMeanMotion * finalTime ) );
            expectedFinalState[ 1 ] = 0.0;
            expectedFinalState[ 2 ]
                = thrustAcceleration [ 2 ] / ( targetMeanMotion * targetMeanMotion )
                  * ( 1.0 - std::cos( targetMeanMotion * finalTime ) );
            expectedFinalState[ 3 ]
                = 2.0 * thrustAcceleration [ 2 ] / targetMeanMotion
                  * ( 1.0 - std::cos( targetMeanMotion * finalTime ) );
            expectedFinalState[ 4 ] = 0.0;
            expectedFinalState[ 5 ]
                = thrustAcceleration[ 2 ] / targetMeanMotion
                  * std::sin( targetMeanMotion * finalTime );

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }
    }

    SECTION( "Test internal snapshot" )
    {
        // The tests presented here are only to confirm internal consistency by testing against a
        // snapshot of the results obtained at the time of writing the test.

        SECTION( "Test arbitrary free motion" )
        {
            Vector thrustAcceleration( 3 );
            thrustAcceleration[ 0 ] = 0.0;
            thrustAcceleration[ 1 ] = 0.0;
            thrustAcceleration[ 2 ] = 0.0;

            Vector initialState( 6 );
            initialState[ 0 ] = 15.613;
            initialState[ 1 ] = -1.6136;
            initialState[ 2 ] = 43.123;
            initialState[ 3 ] = -1.35;
            initialState[ 4 ] = 0.612;
            initialState[ 5 ] = -5.699;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = -177.220096793708;
            expectedFinalState[ 1 ] = 59.4861383364061;
            expectedFinalState[ 2 ] = -511.134488585295;
            expectedFinalState[ 3 ] = -2.49554339092689;
            expectedFinalState[ 4 ] = 0.608907076143407;
            expectedFinalState[ 5 ] = -5.3762829430049;

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }

        SECTION( "Test arbitrary constant-force motion" )
        {
            Vector thrustAcceleration( 3 );
            thrustAcceleration[ 0 ] = 3.415;
            thrustAcceleration[ 1 ] = -1.556;
            thrustAcceleration[ 2 ] = 8.821;

            Vector initialState( 6 );
            initialState[ 0 ] = 15.613;
            initialState[ 1 ] = -1.6136;
            initialState[ 2 ] = 43.123;
            initialState[ 3 ] = -1.35;
            initialState[ 4 ] = 0.612;
            initialState[ 5 ] = -5.699;

            Vector expectedFinalState( 6 );
            expectedFinalState[ 0 ] = 19873.9479718838;
            expectedFinalState[ 1 ] = -7713.59262479076;
            expectedFinalState[ 2 ] = 42378.8990414098;
            expectedFinalState[ 3 ] = 427.649888486457;
            expectedFinalState[ 4 ] = -154.714292723343;
            expectedFinalState[ 5 ] = 839.895191983715;

            Vector computedFinalState = propagateClohessyWiltshireSolution( initialState,
                                                                            finalTime,
                                                                            targetMeanMotion,
                                                                            thrustAcceleration );

            for ( unsigned int i = 0; i < 6; ++i )
            {
                REQUIRE( computedFinalState[ i ]
                            == Approx( expectedFinalState[ i ] ).epsilon( tolerance ) );
            }
        }
    }
}

} // namespace tests
} // namespace astro

/*!
 * References
 *  Fehse, W. (2003) Automated Rendezvous and Docking of Spacecraft, Cambridge Aerospace Series 16,
 *    Cambridge University Press, Cambridge, United Kingdom.
 */
