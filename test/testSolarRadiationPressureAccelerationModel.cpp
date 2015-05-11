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

TEST_CASE( "Obtain the Solar Pres. Acc.: test 1", "[solar_radiation_pressure, acceleration, models]" )
{
    // Set expected solar radiation pressure acceleration vector [m/s^2].
    Vector expectedAcceleration( 3 );
    expectedAcceleration[ 0 ] = -2.964e-06;
    expectedAcceleration[ 1 ] = 0.0;
    expectedAcceleration[ 2 ] = 0.0;

    // Set 1 AU in metres [m].
    const double astronomicalUnitInMeters = 1.49598e11;

    // Set epsilon = error between expected value and computed value.
    const Real epsilon = 1.0e-15;

    // Solar Radiation Pressure at 1 AU [N m^-2]
    const Real radPres = 4.56e-6;

    // Set Radiation Pressure Coef.
    const Real radPresCoef = 1.0 + 0.3;

    // Set Absorbing Area [m^2]
    const Real area = 2.0;    

    // Set mass [kg].
    const Real mass = 4.0;

    // Unit vector poitn from S/C to sun 
    Vector vectorToSource( 3 );
    vectorToSource[ 0 ] = astronomicalUnitInMeters;
    vectorToSource[ 1 ] = 0.0;
    vectorToSource[ 2 ] = 0.0;

    // Compute the 
    const Real normVectorToSource =  sqrt(   vectorToSource[ 0 ] * vectorToSource[ 0 ]
                                           + vectorToSource[ 1 ] * vectorToSource[ 1 ]
                                           + vectorToSource[ 2 ] * vectorToSource[ 2 ] );

    const Real squaredNormVectorToSource = normVectorToSource * normVectorToSource;

    Vector unitVectorToSource( 3 );
    unitVectorToSource[ 0 ] = vectorToSource[ 0 ] / normVectorToSource;
    unitVectorToSource[ 1 ] = vectorToSource[ 1 ] / normVectorToSource;
    unitVectorToSource[ 2 ] = vectorToSource[ 2 ] / normVectorToSource;
    

   // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget =  radPres
                                             * astronomicalUnitInMeters * astronomicalUnitInMeters
                                             / squaredNormVectorToSource;

     // Compute the Solar radiation pressure [m^-2]
    const Vector computedAcceleration = computeSolarRadiationPressureAcceleration(radiationPressureAtTarget,  // Solar Radiation Presure                     [N m^-2]
                                                                                  radPresCoef,      // Radiation Pressure Coefficiente             [adim]
                                                                                  unitVectorToSource, // Unit vector pointing from S/C to sun (3x1)  [adim]                                            
                                                                                  area,             // Absorbing Area of S/C                       [m^2]
                                                                                  mass );           // Mass of the S/C                             [kg]    

    // Check if computed mean motion matches expected value.
    REQUIRE( std::fabs( computedAcceleration[ 0 ] - expectedAcceleration[ 0 ] ) <= epsilon );
    REQUIRE( std::fabs( computedAcceleration[ 1 ] - expectedAcceleration[ 1 ] ) <= epsilon );
    REQUIRE( std::fabs( computedAcceleration[ 2 ] - expectedAcceleration[ 2 ] ) <= epsilon );
}

} // namespace tests
} // namespace astro