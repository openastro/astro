/*
 * Copyright (c) 2015 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cmath>
#include <vector>

#include <catch.hpp>

#include "Astro/constants.hpp"
#include "Astro/solarRadiationPressureAccelerationModel.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector;

TEST_CASE( "Obtain the solar radiation pressure acceleration: test 1", "[solar_radiation_Pressure, acceleration, models]" )
{
    // Set expected solar radiation pressure acceleration vector [m/s^2].
    Vector expectedAcceleration( 3 );
    expectedAcceleration[ 0 ] = -2.964e-06;
    expectedAcceleration[ 1 ] = 0.0;
    expectedAcceleration[ 2 ] = 0.0;

    // Set epsilon = error between expected value and computed value.
    const Real epsilon = 1.0e-10;

    // Solar Radiation Pressure at 1 AU [N m^-2]
    const Real radPres = ;

    // Set Radiation Pressure Coef.
    const Real radPresCoef = 1.0 + 0.3;

    // Set Absorbing Area [m^2]
    const Real area = 2.0;    

    // Set mass [kg].
    const Real mass = 4.0;

    // Velocity vector in [m/s].
    Vector velocity( 3 );
    velocity[ 0 ] = ASTRO_AU_IN_KM;
    velocity[ 1 ] = 0.0;
    velocity[ 2 ] = 0.0;

     // Compute the Solar radiation pressure [m^-2]
    const Vector computedAcceleration = computeSolarRadiationPressureAcceleration(radPres,          // Solar Radiation Presure                     [N m^-2]
                                                                                  radPresCoef,      // Radiation Pressure Coefficiente             [adim]
                                                                                  vectorToSource,   // Unit vector pointing from S/C to sun (3x1)  [adim]                                            
                                                                                  area,             // Absorbing Area of S/C                       [m^2]
                                                                                  mass );           // Mass of the S/C                             [kg]    


    // Check if computed mean motion matches expected value.
    REQUIRE( std::fabs( computedAcceleration[ 0 ] - expectedAcceleration[ 0 ] ) <= epsilon );
    REQUIRE( std::fabs( computedAcceleration[ 1 ] - expectedAcceleration[ 1 ] ) <= epsilon );
    REQUIRE( std::fabs( computedAcceleration[ 2 ] - expectedAcceleration[ 2 ] ) <= epsilon );
}

} // namespace tests
} // namespace astro