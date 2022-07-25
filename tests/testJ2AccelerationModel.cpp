/*
 * Copyright (c) 2014-2022 Kartik Kumar (me@kartikkumar.com)
 * Copyright (c) 2014-2015 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2015 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2015 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <vector>

#include "astro/centralBodyAccelerationModel.hpp"
#include "astro/j2AccelerationModel.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef std::vector<Real> Vector;

TEST_CASE("Compute central+J2 acceleration vector for spacecraft around Mercury",
           "[central_gravity, j2_gravity, acceleration, models]")
{
    // Benchmark values for test case obtained using the gravityzonal() function in MATLAB.

    // Set expected acceleration vector (central+J2) [km s^-2].
    Vector expectedAcceleration(3);
    expectedAcceleration[0] = -6.174568462599339e-02;
    expectedAcceleration[1] = 3.024518496375884e-01;
    expectedAcceleration[2] = -1.229017246366501e-01;

     // Set tolerance = error between expected value and computed value.
    const Real tolerance = 1.0e-15;

    // Set value of gravitational parameter of central body [km^3 s^-2].
    const Real gravitationalParameter = 2.2032e13;

    // Se the J2 coefficient of the spherical harmonics expansion and the corresponding equatorial
    // radius [km].
    const Real equatorialRadius       = 2439.0e3;
    const Real j2Coefficient          = 0.00006;

    // Set position vector of the orbiting body relative to the origin of the reference frame [km].
    Vector position(3);
    position[0] = 1513.3e3;
    position[1] = -7412.67e3;
    position[2] = 3012.1e3;

    // Compute the acceleration vector due to the central body [km s^-2].
    const Vector computedCentralBodyAcceleration
        = computeCentralBodyAcceleration(gravitationalParameter, position);

    // Compute the acceleration vector due to J2 [km s^-2].
    const Vector computedJ2Acceleration = computeJ2Acceleration(gravitationalParameter,
                                                                position,
                                                                equatorialRadius,
                                                                j2Coefficient);

    // Check if computed acceleration matches expected values.
    REQUIRE((computedCentralBodyAcceleration[0] + computedJ2Acceleration[0])
                == Catch::Approx(expectedAcceleration[0]).epsilon(tolerance));
    REQUIRE((computedCentralBodyAcceleration[1] + computedJ2Acceleration[1])
                == Catch::Approx(expectedAcceleration[1]).epsilon(tolerance));
    REQUIRE((computedCentralBodyAcceleration[2] + computedJ2Acceleration[2])
                == Catch::Approx(expectedAcceleration[2]).epsilon(tolerance));
}

} // namespace tests
} // namespace astro
