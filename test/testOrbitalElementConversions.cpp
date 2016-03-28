/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <limits>
#include <vector>

#include <catch.hpp>

#include <SML/sml.hpp>

#include "Astro/orbitalElementConversions.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector;

TEST_CASE( "Convert Cartesian elements to Keplerian elements",
           "[cartesian-to-keplerian-elements]" )
{
    SECTION( "Test elliptical orbit around the Earth" )
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set Earth gravitational parameter [m^3 s^-2] .
        const Real earthGravitationalParameter = 3.986004415e14;

        // Set Cartesian elements.
        Vector cartesianElements( 6 );
        cartesianElements[ 0 ] = 3.75e6;
        cartesianElements[ 1 ] = 4.24e6;
        cartesianElements[ 2 ] = -1.39e6;
        cartesianElements[ 3 ] = -4.65e3;
        cartesianElements[ 4 ] = -2.21e3;
        cartesianElements[ 5 ] = 1.66e3;

        // Set expected Keplerian elements.
        Vector keplerianElements( 6 );
        keplerianElements[ 0 ] = 3.707478199246163e6;
        keplerianElements[ 1 ] = 0.949175203660321;
        keplerianElements[ 2 ] = 0.334622356632438;
        keplerianElements[ 3 ] = 2.168430616511167;
        keplerianElements[ 4 ] = 1.630852596545341;
        keplerianElements[ 5 ] = 3.302032232567084;

        // Compute Keplerian elements.
        Vector result( 6 );
        result = convertCartesianToKeplerianElements< Real, Vector >(
            cartesianElements, earthGravitationalParameter );

        // Loop through vectors and require that each element of the result matches the expected
        // value to the given tolerance.
        for ( unsigned int i = 0; i < result.size( ); i++ )
        {
            REQUIRE( keplerianElements[ i ] == Approx( result[ i ] ).epsilon( 1.0e-14 ) );
        }
    }
}

TEST_CASE( "Convert Keplerian elements to Cartesian elements",
           "[keplerian-to-cartesian-elements]" )
{
    // Set PI.
    const Real pi = 3.14159265358979323846;

    SECTION( "Test elliptical orbit around the Earth" )
    {
        SECTION( "Test using ODTBX")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            // Set Earth gravitational parameter [m^3/s^2] .
            const Real earthGravitationalParameter = 3.986004415e14;

            // Set Keplerian elements [m,-,rad,rad,rad,rad].
            Vector keplerianElements( 6 );
            keplerianElements[ semiMajorAxisIndex ]            = 8000.0 * 1000.0;
            keplerianElements[ eccentricityIndex ]             = 0.23;
            keplerianElements[ inclinationIndex ]              = 20.6 / 180.0 * pi;
            keplerianElements[ argumentOfPeriapsisIndex ]      = 274.78 / 180.0 * pi;
            keplerianElements[ longitudeOfAscendingNodeIndex ] = 108.77 / 180.0 * pi;
            keplerianElements[ trueAnomalyIndex ]              = 46.11 / 180.0 * pi;

            // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
            Vector cartesianElements( 6 );
            cartesianElements[ xPositionIndex ] = 2.021874804243437e6;
            cartesianElements[ yPositionIndex ] = 6.042523817035284e6;
            cartesianElements[ zPositionIndex ] = -1.450371183512575e6;
            cartesianElements[ xVelocityIndex ] = -7.118283509842652e3;
            cartesianElements[ yVelocityIndex ] = 4.169050171542199e3;
            cartesianElements[ zVelocityIndex ] = 2.029066072016241e3;

            // Compute Cartesian elements.
            Vector result = convertKeplerianToCartesianElements< Real, Vector >(
                keplerianElements, earthGravitationalParameter );

            // Loop through vectors and require that each element of the result matches the expected
            // value to the given tolerance.
            for ( unsigned int i = 0; i < result.size( ); i++ )
            {
                REQUIRE( cartesianElements[ i ] == Approx( result[ i ] ).epsilon( 1.0e-14 ) );
            }
        }

        SECTION( "Test using data from course on Mission Geometry and Orbit Design at TU Delft" )
        {
            // The benchmark data is obtained from the lecture notes (Noomen, 2014).

            SECTION( "Test using example 1" )
            {
                // Set Earth gravitational parameter [m^3/s^2] .
                // const Real earthGravitationalParameter = 3.986004415e14;
                const Real earthGravitationalParameter = 3.98600441e14;

                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                Vector keplerianElements( 6 );
                keplerianElements[ semiMajorAxisIndex ]            = 6787746.891;
                keplerianElements[ eccentricityIndex ]             = 0.000731104;
                keplerianElements[ inclinationIndex ]              = 51.68714486 / 180.0 * pi;
                keplerianElements[ argumentOfPeriapsisIndex ]      = 74.21987137 / 180.0 * pi;
                keplerianElements[ longitudeOfAscendingNodeIndex ] = 127.5486706 / 180.0 * pi;
                keplerianElements[ trueAnomalyIndex ]              = 24.10027677 / 180.0 * pi;

                // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
                Vector cartesianElements( 6 );
                cartesianElements[ xPositionIndex ]                = -2700816.14;
                cartesianElements[ yPositionIndex ]                = -3314092.80;
                cartesianElements[ zPositionIndex ]                = 5266346.42;
                cartesianElements[ xVelocityIndex ]                = 5168.606550;
                cartesianElements[ yVelocityIndex ]                = -5597.546618;
                cartesianElements[ zVelocityIndex ]                = -868.878445;

                // Compute Cartesian elements.
                Vector result = convertKeplerianToCartesianElements< Real, Vector >(
                    keplerianElements, earthGravitationalParameter );

                // Loop through vectors and require that each element of the result matches the expected
                // value to the given tolerance.
                for ( unsigned int i = 0; i < result.size( ); i++ )
                {
                    REQUIRE( cartesianElements[ i ] == Approx( result[ i ] ).epsilon( 1.0e-9 ) );
                }
            }

            SECTION( "Test using example 2" )
            {
                // Set Earth gravitational parameter [m^3/s^2].
                const Real earthGravitationalParameter = 3.98600441e14;

                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                Vector keplerianElements( 6 );
                keplerianElements[ semiMajorAxisIndex ]            = 7096137.00;
                keplerianElements[ eccentricityIndex ]             = 0.0011219;
                keplerianElements[ inclinationIndex ]              = 92.0316 / 180.0 * pi;
                keplerianElements[ argumentOfPeriapsisIndex ]      = 120.6878 / 180.0 * pi;
                keplerianElements[ longitudeOfAscendingNodeIndex ] = 296.1384 / 180.0 * pi;
                keplerianElements[ trueAnomalyIndex ]              = 239.5437 / 180.0 * pi;

                // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
                Vector cartesianElements( 6 );
                cartesianElements[ xPositionIndex ]                = 3126974.99;
                cartesianElements[ yPositionIndex ]                = -6374445.74;
                cartesianElements[ zPositionIndex ]                = 28673.59;
                cartesianElements[ xVelocityIndex ]                = -254.91197;
                cartesianElements[ yVelocityIndex ]                = -83.30107;
                cartesianElements[ zVelocityIndex ]                = 7485.70674;

                // Compute Cartesian elements.
                Vector result = convertKeplerianToCartesianElements< Real, Vector >(
                    keplerianElements, earthGravitationalParameter );

                // Loop through vectors and require that each element of the result matches the expected
                // value to the given tolerance.
                for ( unsigned int i = 0; i < result.size( ); i++ )
                {
                    REQUIRE( cartesianElements[ i ] == Approx( result[ i ] ).epsilon( 1.0e-3 ) );
                }
            }

        }
    }
}

TEST_CASE( "Convert true anomaly to eccentric anomaly" , "[true-to-eccentric-anomaly]" )
{
    SECTION( "Test run-time errors" )
    {
        SECTION( "Test negative eccentricities" )
        {
            REQUIRE_THROWS( convertTrueAnomalyToEllipticalEccentricAnomaly( 1.234, -0.152 ) );
            REQUIRE_THROWS( convertTrueAnomalyToHyperbolicEccentricAnomaly( 1.234, -0.152 ) );
            REQUIRE_THROWS( convertTrueAnomalyToEccentricAnomaly( 1.234, -0.152 ) );
        }

        SECTION( "Test parabolic eccentricity" )
        {
            REQUIRE_THROWS( convertTrueAnomalyToEllipticalEccentricAnomaly( 1.234, 1.0 ) );
            REQUIRE_THROWS( convertTrueAnomalyToHyperbolicEccentricAnomaly( 1.234, 1.0 ) );
            REQUIRE_THROWS( convertTrueAnomalyToEccentricAnomaly( 1.234, 1.0 ) );
        }

        SECTION( "Test hyperbolic eccentricity for elliptical conversion" )
        {
            REQUIRE_THROWS( convertTrueAnomalyToEllipticalEccentricAnomaly( 1.234, 2.345 ) );

        }

        SECTION( "Test elliptical eccentricity for hyperbolic conversion" )
        {
            REQUIRE_THROWS( convertTrueAnomalyToHyperbolicEccentricAnomaly( 1.234, 0.152 ) );
        }
    }

    SECTION( "Test elliptical orbits" )
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set eccentricities [-].
        const Real eccentricities[ 3 ] = { 0.146, 0.0, 0.0 };

        // Set true anomalies [rad].
        const Real trueAnomalies[ 3 ]
            = { 82.16 / 180.0 * sml::SML_PI,
                160.43 / 180.0 * sml::SML_PI,
                0.0 };

        // Set expected eccentric anomalies [rad].
        const Real expectedEccentricAnomalies[ 3 ]
            = { 1.290237398010989,
                2.800031718974503,
                0.0 };

        // Loop through all the cases, compute the eccentric anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 3; i++ )
        {
            // Compute eccentric anomaly [rad].
            const Real computedEllipticalEccentricAnomaly
                = convertTrueAnomalyToEllipticalEccentricAnomaly(
                    trueAnomalies[ i ], eccentricities[ i ] );

            const Real computedEccentricAnomaly
                = convertTrueAnomalyToEccentricAnomaly(
                    trueAnomalies[ i ], eccentricities[ i ] );

            // Check if computed eccentric anomaly matches the expected value.
            REQUIRE( computedEllipticalEccentricAnomaly
                     == Approx( expectedEccentricAnomalies[ i ] ).epsilon(
                                    std::numeric_limits< Real >::epsilon( ) ) );

            REQUIRE( computedEccentricAnomaly
                     == Approx( expectedEccentricAnomalies[ i ] ).epsilon(
                                    std::numeric_limits< Real >::epsilon( ) ) );
        }
    }

    SECTION( "Test general hyperbolic orbits" )
    {
        // The benchmark data is obtained from (Fortescue, 2003).

        // Set eccentricities [-].
        const Real eccentricities[ 1 ] = { 3.0 };

        // Set true anomalies [rad].
        const Real trueAnomalies[ 1 ] = { 0.5291 };

        // Set expected eccentric anomalies [rad].
        const Real expectedEccentricAnomalies[ 1 ] = { 0.3879 };

        // Loop through all the cases, compute the eccentric anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 1; i++ )
        {
            // Compute eccentric anomaly [rad].
            const Real computedHyperbolicEccentricAnomaly
                = convertTrueAnomalyToHyperbolicEccentricAnomaly(
                    trueAnomalies[ i ], eccentricities[ i ] );

            const Real computedEccentricAnomaly
                = convertTrueAnomalyToEccentricAnomaly(
                    trueAnomalies[ i ], eccentricities[ i ] );

            // Check if computed eccentric anomaly matches the expected value.
            REQUIRE( computedHyperbolicEccentricAnomaly
                     == Approx( expectedEccentricAnomalies[ i ] ).epsilon( 1.0e-5 ) );

            REQUIRE( computedEccentricAnomaly
                     == Approx( expectedEccentricAnomalies[ i ] ).epsilon( 1.0e-5 ) );
        }
    }
}

TEST_CASE( "Convert eccentric anomaly to mean anomaly" , "[eccentric-to-mean-anomaly]" )
{
    SECTION( "Test run-time errors" )
    {
        SECTION( "Test negative eccentricities" )
        {
            REQUIRE_THROWS( convertEllipticalEccentricAnomalyToMeanAnomaly( 1.234, -0.152 ) );
            REQUIRE_THROWS( convertHyperbolicEccentricAnomalyToMeanAnomaly( 1.234, -0.152 ) );
            REQUIRE_THROWS( convertEccentricAnomalyToMeanAnomaly( 1.234, -0.152 ) );
        }

        SECTION( "Test parabolic eccentricity" )
        {
            REQUIRE_THROWS( convertEllipticalEccentricAnomalyToMeanAnomaly( 1.234, 1.0 ) );
            REQUIRE_THROWS( convertHyperbolicEccentricAnomalyToMeanAnomaly( 1.234, 1.0 ) );
            REQUIRE_THROWS( convertEccentricAnomalyToMeanAnomaly( 1.234, 1.0 ) );
        }

        SECTION( "Test hyperbolic eccentricity for elliptical conversion" )
        {
            REQUIRE_THROWS( convertEllipticalEccentricAnomalyToMeanAnomaly( 1.234, 2.345 ) );

        }

        SECTION( "Test elliptical eccentricity for hyperbolic conversion" )
        {
            REQUIRE_THROWS( convertHyperbolicEccentricAnomalyToMeanAnomaly( 1.234, 0.152 ) );
        }
    }

    SECTION( "Test elliptical orbits" )
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set eccentricities [-].
        const Real eccentricities[ 3 ] = { 0.541, 0.0, 0.0 };

        // Set elliptical eccentric anomalies [rad].
        const Real ellipticalEccentricAnomalies[ 3 ]
            = { 176.09 / 180.0 * sml::SML_PI,
                320.12 / 180.0 * sml::SML_PI,
                0.0 };

        // Set expected mean anomalies [rad].
        const Real expectedMeanAnomalies[ 3 ]
            = { 3.036459804491048,
                5.587148001484247,
                0.0 };

        // Loop through all the cases, compute the mean anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 3; i++ )
        {
            // Compute mean anomaly [rad].
            const Real computedMeanAnomaly
                = convertEllipticalEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomalies[ i ], eccentricities[ i ] );

            const Real computedMeanAnomalyWrapper
                = convertEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomalies[ i ], eccentricities[ i ] );

            // Check if computed mean anomaly matches the expected value.
            REQUIRE( computedMeanAnomaly
                     == Approx( expectedMeanAnomalies[ i ] ).epsilon(
                                    std::numeric_limits< Real >::epsilon( ) ) );

            REQUIRE( computedMeanAnomalyWrapper
                     == Approx( expectedMeanAnomalies[ i ] ).epsilon(
                                    std::numeric_limits< Real >::epsilon( ) ) );
        }
    }

    SECTION( "Test hyperbolic orbits" )
    {
        // The benchmark data is obtained from (Vallado, 2004).

        // Set eccentricities [-].
        const Real eccentricities[ 1 ] = { 2.4 };

        // Set hyperbolic eccentric anomalies [rad].
        const Real hyperbolicEccentricAnomalies[ 1 ] = { 1.6013761449 };

        // Set expected mean anomalies [rad].
        const Real expectedMeanAnomalies[ 1 ] = { 235.4 / 180.0 * sml::SML_PI };

        // Loop through all the cases, compute the mean anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 1; i++ )
        {
            // Compute mean anomaly [rad].
            const Real computedMeanAnomaly
                = convertHyperbolicEccentricAnomalyToMeanAnomaly(
                    hyperbolicEccentricAnomalies[ i ], eccentricities[ i ] );

            const Real computedMeanAnomalyWrapper
                = convertEccentricAnomalyToMeanAnomaly(
                    hyperbolicEccentricAnomalies[ i ], eccentricities[ i ] );

            // Check if computed mean anomaly matches the expected value.
            REQUIRE( computedMeanAnomaly
                     == Approx( expectedMeanAnomalies[ i ] ).epsilon( 1.0e-10 ) );

            REQUIRE( computedMeanAnomalyWrapper
                     == Approx( expectedMeanAnomalies[ i ] ).epsilon( 1.0e-10 ) );
        }
    }
}

} // namespace tests
} // namespace astro

/*!
 * References
 *  Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition, Wiley, England, 2003.
 *  NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *   Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *   31st January, 2012.
 *  Noomen, R. Space Mission Design: Basics, ae4-878 (Mission Geometry & Orbit Design), Delft
 *   University of Technology, The Netherlands, 2014.
 */
