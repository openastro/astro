/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#define REAL double

#include <limits>
#include <vector>

#include <catch.hpp>

#include <SML/sml.hpp>

#include <SAM/constants.hpp>
#include <SAM/orbitalElementConversions.hpp>

namespace sam
{
namespace unit_tests
{

typedef std::vector< REAL > Vector6;
typedef std::vector< REAL > Vector3;

TEST_CASE( "Convert Cartesian elements to Keplerian elements", "[cartesian-to-keplerian-elements]")
{
    SECTION( "Test elliptical orbit around the Earth" )
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set Earth gravitational parameter [m^3 s^-2] .
        const REAL earthGravitationalParameter = 3.986004415e14;

        // Set Cartesian elements.
        Vector6 cartesianElements( 6 );
        cartesianElements[ xPositionIndex ] = 3.75e6;
        cartesianElements[ yPositionIndex ] = 4.24e6;
        cartesianElements[ zPositionIndex ] = -1.39e6;
        cartesianElements[ xVelocityIndex ] = -4.65e3;
        cartesianElements[ yVelocityIndex ] = -2.21e3;
        cartesianElements[ zVelocityIndex ] = 1.66e3;

        // Set expected Keplerian elements.
        Vector6 keplerianElements( 6 );
        keplerianElements[ semiMajorAxisIndex ] = 3.707478199246163e6;
        keplerianElements[ eccentricityIndex ] = 0.949175203660321;
        keplerianElements[ inclinationIndex ] = 0.334622356632438;
        keplerianElements[ argumentOfPeriapsisIndex ] = 2.168430616511167;
        keplerianElements[ longitudeOfAscendingNodeIndex ] = 1.630852596545341;
        keplerianElements[ trueAnomalyIndex ] = 3.302032232567084;

        // Compute Keplerian elements.
        Vector6 result( 6 );
        result = convertCartesianToKeplerianElements< REAL, Vector6, Vector3 >( 
            cartesianElements, earthGravitationalParameter );

        // Loop through vectors and require that each element of the result matches the expected 
        // value to the given tolerance.
        for ( unsigned int i = 0; i < result.size( ); i++ )
        {
            REQUIRE( keplerianElements[ i ] == Approx( result[ i ] ).epsilon( 1.0e-14 ) );        
        }        
    } 
}

TEST_CASE( "Convert true anomaly to eccentric anomaly" , "[true-to-eccentric-anomaly]")
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
        const REAL eccentricities[ 3 ] = { 0.146, 0.0, 0.0 };

        // Set true anomalies [rad].
        const REAL trueAnomalies[ 3 ] 
            = { 82.16 / 180.0 * sml::SML_PI, 
                160.43 / 180.0 * sml::SML_PI,
                0.0 };

        // Set expected eccentric anomalies [rad].
        const REAL expectedEccentricAnomalies[ 3 ] 
            = { 1.290237398010989, 
                2.800031718974503,
                0.0 };

        // Loop through all the cases, compute the eccentric anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 3; i++ )
        {
            // Compute eccentric anomaly [rad].
            const REAL computedEllipticalEccentricAnomaly 
                = convertTrueAnomalyToEllipticalEccentricAnomaly( 
                    trueAnomalies[ i ], eccentricities[ i ] );

            const REAL computedEccentricAnomaly 
                = convertTrueAnomalyToEccentricAnomaly( 
                    trueAnomalies[ i ], eccentricities[ i ] );

            // Check if computed eccentric anomaly matches the expected value.
            REQUIRE( computedEllipticalEccentricAnomaly 
                     == Approx( expectedEccentricAnomalies[ i ] ).epsilon( 
                                    std::numeric_limits< REAL >::epsilon( ) ) );        

            REQUIRE( computedEccentricAnomaly 
                     == Approx( expectedEccentricAnomalies[ i ] ).epsilon( 
                                    std::numeric_limits< REAL >::epsilon( ) ) );              
        }
    }

    SECTION( "Test general hyperbolic orbits" )
    {
        // The benchmark data is obtained from (Fortescue, 2003).

        // Set eccentricities [-].
        const REAL eccentricities[ 1 ] = { 3.0 };        

        // Set true anomalies [rad].
        const REAL trueAnomalies[ 1 ] = { 0.5291 };        

        // Set expected eccentric anomalies [rad].
        const REAL expectedEccentricAnomalies[ 1 ] = { 0.3879 };        

        // Loop through all the cases, compute the eccentric anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 1; i++ )
        {
            // Compute eccentric anomaly [rad].
            const REAL computedHyperbolicEccentricAnomaly 
                = convertTrueAnomalyToHyperbolicEccentricAnomaly( 
                    trueAnomalies[ i ], eccentricities[ i ] );

            const REAL computedEccentricAnomaly 
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

TEST_CASE( "Convert eccentric anomaly to mean anomaly" , "[eccentric-to-mean-anomaly]")
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
        const REAL eccentricities[ 3 ] = { 0.541, 0.0, 0.0 };

        // Set elliptical eccentric anomalies [rad].
        const REAL ellipticalEccentricAnomalies[ 3 ] 
            = { 176.09 / 180.0 * sml::SML_PI,
                320.12 / 180.0 * sml::SML_PI,
                0.0 };

        // Set expected mean anomalies [rad].
        const REAL expectedMeanAnomalies[ 3 ] 
            = { 3.036459804491048,
                5.587148001484247,
                0.0 };

        // Loop through all the cases, compute the mean anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 3; i++ )
        {
            // Compute mean anomaly [rad].
            const REAL computedMeanAnomaly 
                = convertEllipticalEccentricAnomalyToMeanAnomaly( 
                    ellipticalEccentricAnomalies[ i ], eccentricities[ i ] );

            const REAL computedMeanAnomalyWrapper 
                = convertEccentricAnomalyToMeanAnomaly( 
                    ellipticalEccentricAnomalies[ i ], eccentricities[ i ] );

            // Check if computed mean anomaly matches the expected value.
            REQUIRE( computedMeanAnomaly 
                     == Approx( expectedMeanAnomalies[ i ] ).epsilon( 
                                    std::numeric_limits< REAL >::epsilon( ) ) );        

            REQUIRE( computedMeanAnomalyWrapper 
                     == Approx( expectedMeanAnomalies[ i ] ).epsilon( 
                                    std::numeric_limits< REAL >::epsilon( ) ) );              
        }
    }    

    SECTION( "Test hyperbolic orbits" )
    {
        // The benchmark data is obtained from (Vallado, 2004).

        // Set eccentricities [-].
        const REAL eccentricities[ 1 ] = { 2.4 };

        // Set hyperbolic eccentric anomalies [rad].
        const REAL hyperbolicEccentricAnomalies[ 1 ] = { 1.6013761449 };

        // Set expected mean anomalies [rad].
        const REAL expectedMeanAnomalies[ 1 ] = { 235.4 / 180.0 * sml::SML_PI };    

        // Loop through all the cases, compute the mean anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for ( unsigned int i = 0; i < 1; i++ )
        {
            // Compute mean anomaly [rad].
            const REAL computedMeanAnomaly 
                = convertHyperbolicEccentricAnomalyToMeanAnomaly( 
                    hyperbolicEccentricAnomalies[ i ], eccentricities[ i ] );

            const REAL computedMeanAnomalyWrapper 
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

TEST_CASE( "Convert semi-major axis to mean motion", "[semi-major-axis-to-mean-motion]")
{
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit.

    // Set satellite mass [kg].
    const REAL satelliteMass = 1.0e3;

    // Set gravitational parameter of Earth [m^3 s^-2].
    const REAL earthGravitationalParameter = SAM_GRAVITATIONAL_CONSTANT * 5.9736e24;

    // Set distance between Earth center and satellite [m].
    const REAL distanceBetweenSatelliteAndEarth = 4.2164e7;

    // Set expected mean motion [rad/s].
    const REAL expectedMeanMotion = 7.2921e-5;

    // Compute mean motion.
    const REAL meanMotion = computeKeplerMeanMotion(
        distanceBetweenSatelliteAndEarth, earthGravitationalParameter, satelliteMass );

    // Check if computed mean motion matches expected value.
    REQUIRE( meanMotion == Approx( expectedMeanMotion ).epsilon( 1.0e-7 ) );
}

} // namespace unit_tests
} // namespace sam

/*!
 * References
 *  NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *   Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *   31st January, 2012.
 *  Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition, Wiley, England, 2003. 
 *  Vallado, D. A., McClain, W. D. Fundamentals of astrodynamics and applications, 2nd Edition,
 *   Kluwer Academic Publishers, The Netherlands, 2004. 
 */
