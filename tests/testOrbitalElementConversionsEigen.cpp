/*    
 * Copyright (c) 2014 K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <limits>

#include <catch.hpp>

#include <Eigen/Core>

#include <SML/sml.hpp>

#include <SAM/orbitalElementConversions.hpp>

namespace sam
{
namespace tests
{

typedef double Real;
typedef Eigen::Matrix< Real, Eigen::Dynamic, 1 > Vector;

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
        cartesianElements[ xPositionIndex ] = 3.75e6;
        cartesianElements[ yPositionIndex ] = 4.24e6;
        cartesianElements[ zPositionIndex ] = -1.39e6;
        cartesianElements[ xVelocityIndex ] = -4.65e3;
        cartesianElements[ yVelocityIndex ] = -2.21e3;
        cartesianElements[ zVelocityIndex ] = 1.66e3;

        // Set expected Keplerian elements.
        Vector keplerianElements( 6 );
        keplerianElements[ semiMajorAxisIndex ]            = 3.707478199246163e6;
        keplerianElements[ eccentricityIndex ]             = 0.949175203660321;
        keplerianElements[ inclinationIndex ]              = 0.334622356632438;
        keplerianElements[ argumentOfPeriapsisIndex ]      = 2.168430616511167;
        keplerianElements[ longitudeOfAscendingNodeIndex ] = 1.630852596545341;
        keplerianElements[ trueAnomalyIndex ]              = 3.302032232567084;

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
} // namespace sam

/*!
 * References
 *  NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *   Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *   31st January, 2012.
 *  Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition, Wiley, England, 2003. 
 */
