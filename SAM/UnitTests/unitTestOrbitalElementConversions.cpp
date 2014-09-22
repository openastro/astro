/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#define REAL double

#include <vector>

#include <catch.hpp>

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

        // Loop through vectors and require that each element of the result matches the expected value
        // to the given tolerance.
        for ( unsigned int i = 0; i < result.size( ); i++ )
        {
            REQUIRE( keplerianElements[ i ] == Approx( result[ i ] ).epsilon( 1.0e-14 ) );        
        }        
    } 
}

} // namespace unit_tests
} // namespace sam

/*!
 * References
 *      NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *          Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *          31st January, 2012.
 */
