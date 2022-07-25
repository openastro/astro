/*
 * Copyright (c) 2014-2022 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <limits>
#include <vector>

#include "astro/orbitalElementConversions.hpp"

namespace astro
{
namespace tests
{

typedef double Real;
typedef int Integer;
typedef std::vector<Real> Vector;

TEST_CASE("Convert Cartesian elements to Keplerian elements",
           "[cartesian-to-keplerian-elements]")
{
    SECTION("Test elliptical orbit around the Earth")
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set Earth gravitational parameter [m^3 s^-2].
        const Real gravitationalParameter = 3.986004415e14;

        // Set Cartesian elements [m].
        Vector cartesianElements(6);
        cartesianElements[0] = 3.75e6;
        cartesianElements[1] = 4.24e6;
        cartesianElements[2] = -1.39e6;
        cartesianElements[3] = -4.65e3;
        cartesianElements[4] = -2.21e3;
        cartesianElements[5] = 1.66e3;

        // Set expected Keplerian elements [m, -, rad, rad, rad, rad].
        Vector keplerianElements(6);
        keplerianElements[0] = 3.707478199246163e6;
        keplerianElements[1] = 0.949175203660321;
        keplerianElements[2] = 0.334622356632438;
        keplerianElements[3] = 2.168430616511167;
        keplerianElements[4] = 1.630852596545341;
        keplerianElements[5] = 3.302032232567084;

        const Vector result = convertCartesianToKeplerianElements<Real, Vector>(
            cartesianElements, gravitationalParameter);

        REQUIRE(keplerianElements[0] == Catch::Approx(result[0]).epsilon(1.0e-14));
        REQUIRE(keplerianElements[1] == Catch::Approx(result[1]).epsilon(1.0e-14));
        REQUIRE(keplerianElements[2] == Catch::Approx(result[2]).epsilon(1.0e-14));
        REQUIRE(keplerianElements[3] == Catch::Approx(result[3]).epsilon(1.0e-14));
        REQUIRE(keplerianElements[4] == Catch::Approx(result[4]).epsilon(1.0e-14));
        REQUIRE(keplerianElements[5] == Catch::Approx(result[5]).epsilon(1.0e-14));
    }

    SECTION("Test equatorial, circular orbit around Venus")
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set Venus gravitational parameter [m^3 s^-2].
        const Real gravitationalParameter = 3.2485504415e14;

        // Set Cartesian elements [m].
        Vector cartesianElements(6);
        cartesianElements[0] = 5.580537430785387e6;
        cartesianElements[1] = 2.816487703435473e6;
        cartesianElements[2] = 0.0;
        cartesianElements[3] = -3.248092722413634e3;
        cartesianElements[4] = 6.435711753323540e3;
        cartesianElements[5] = 0.0;

        // Set expected Keplerian elements [m, -, rad, rad, rad, rad].
        Vector keplerianElements(6);
        keplerianElements[0] = 6.251e6;
        keplerianElements[1] = 0.0;
        keplerianElements[2] = 0.0;
        keplerianElements[3] = 0.0;
        keplerianElements[4] = 0.0;
        keplerianElements[5] = 26.78 / 180.0 * 3.14159265358979323846;

        const Vector result = convertCartesianToKeplerianElements<Real, Vector>(
            cartesianElements, gravitationalParameter);

        REQUIRE(keplerianElements[0] == Catch::Approx(result[0]).epsilon(1.0e-15));
        REQUIRE(std::fabs(result[1]) < std::numeric_limits<Real>::epsilon());
        REQUIRE(std::fabs(result[2]) < std::numeric_limits<Real>::epsilon());
        REQUIRE(std::isnan(result[3]));
        REQUIRE(std::isnan(result[4]));
        REQUIRE(keplerianElements[5] == Catch::Approx(result[5]).epsilon(1.0e-15));
    }

    // @TODO: find a way to test assert statement in convertCartesianToKeplerianElements()

}

TEST_CASE("Convert Keplerian elements to Cartesian elements",
           "[keplerian-to-cartesian-elements]")
{
    const Real pi = 3.14159265358979323846;

    SECTION("Test elliptical orbit around the Earth")
    {
        SECTION("Test using ODTBX")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            // Set Earth gravitational parameter [m^3/s^2] .
            const Real earthGravitationalParameter = 3.986004415e14;

            // Set Keplerian elements [m,-,rad,rad,rad,rad].
            Vector keplerianElements(6);
            keplerianElements[semiMajorAxisIndex]            = 8000.0 * 1000.0;
            keplerianElements[eccentricityIndex]             = 0.23;
            keplerianElements[inclinationIndex]              = 20.6 / 180.0 * pi;
            keplerianElements[argumentOfPeriapsisIndex]      = 274.78 / 180.0 * pi;
            keplerianElements[longitudeOfAscendingNodeIndex] = 108.77 / 180.0 * pi;
            keplerianElements[trueAnomalyIndex]              = 46.11 / 180.0 * pi;

            // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
            Vector cartesianElements(6);
            cartesianElements[xPositionIndex] = 2.021874804243437e6;
            cartesianElements[yPositionIndex] = 6.042523817035284e6;
            cartesianElements[zPositionIndex] = -1.450371183512575e6;
            cartesianElements[xVelocityIndex] = -7.118283509842652e3;
            cartesianElements[yVelocityIndex] = 4.169050171542199e3;
            cartesianElements[zVelocityIndex] = 2.029066072016241e3;

            // Compute Cartesian elements.
            Vector result = convertKeplerianToCartesianElements<Real, Vector>(
                keplerianElements, earthGravitationalParameter);

            // Loop through vectors and require that each element of the result matches the expected
            // value to the given tolerance.
            for (unsigned int i = 0; i < result.size(); i++)
            {
                REQUIRE(cartesianElements[i] == Catch::Approx(result[i]).epsilon(1.0e-14));
            }
        }

        SECTION("Test using data from course on Mission Geometry and Orbit Design at TU Delft")
        {
            // The benchmark data is obtained from the lecture notes (Noomen, 2014).

            SECTION("Test using example 1")
            {
                // Set Earth gravitational parameter [m^3/s^2] .
                const Real earthGravitationalParameter = 3.98600441e14;

                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                Vector keplerianElements(6);
                keplerianElements[semiMajorAxisIndex]            = 6787746.891;
                keplerianElements[eccentricityIndex]             = 0.000731104;
                keplerianElements[inclinationIndex]              = 51.68714486 / 180.0 * pi;
                keplerianElements[argumentOfPeriapsisIndex]      = 74.21987137 / 180.0 * pi;
                keplerianElements[longitudeOfAscendingNodeIndex] = 127.5486706 / 180.0 * pi;
                keplerianElements[trueAnomalyIndex]              = 24.10027677 / 180.0 * pi;

                // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
                Vector cartesianElements(6);
                cartesianElements[xPositionIndex]                = -2700816.14;
                cartesianElements[yPositionIndex]                = -3314092.80;
                cartesianElements[zPositionIndex]                = 5266346.42;
                cartesianElements[xVelocityIndex]                = 5168.606550;
                cartesianElements[yVelocityIndex]                = -5597.546618;
                cartesianElements[zVelocityIndex]                = -868.878445;

                // Compute Cartesian elements.
                Vector result = convertKeplerianToCartesianElements<Real, Vector>(
                    keplerianElements, earthGravitationalParameter);

                // Loop through vectors and require that each element of the result matches the
                // expected value to the given tolerance.
                for (unsigned int i = 0; i < result.size(); i++)
                {
                    REQUIRE(cartesianElements[i] == Catch::Approx(result[i]).epsilon(1.0e-9));
                }
            }

            SECTION("Test using example 2")
            {
                // Set Earth gravitational parameter [m^3/s^2].
                const Real earthGravitationalParameter = 3.98600441e14;

                // Set Keplerian elements [m,-,rad,rad,rad,rad].
                Vector keplerianElements(6);
                keplerianElements[semiMajorAxisIndex]            = 7096137.00;
                keplerianElements[eccentricityIndex]             = 0.0011219;
                keplerianElements[inclinationIndex]              = 92.0316 / 180.0 * pi;
                keplerianElements[argumentOfPeriapsisIndex]      = 120.6878 / 180.0 * pi;
                keplerianElements[longitudeOfAscendingNodeIndex] = 296.1384 / 180.0 * pi;
                keplerianElements[trueAnomalyIndex]              = 239.5437 / 180.0 * pi;

                // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
                Vector cartesianElements(6);
                cartesianElements[xPositionIndex]                = 3126974.99;
                cartesianElements[yPositionIndex]                = -6374445.74;
                cartesianElements[zPositionIndex]                = 28673.59;
                cartesianElements[xVelocityIndex]                = -254.91197;
                cartesianElements[yVelocityIndex]                = -83.30107;
                cartesianElements[zVelocityIndex]                = 7485.70674;

                // Compute Cartesian elements.
                Vector result = convertKeplerianToCartesianElements<Real, Vector>(
                    keplerianElements, earthGravitationalParameter);

                // Loop through vectors and require that each element of the result matches the
                // expected value to the given tolerance.
                for (unsigned int i = 0; i < result.size(); i++)
                {
                    REQUIRE(cartesianElements[i] == Catch::Approx(result[i]).epsilon(1.0e-3));
                }
            }

        }
    }
}

TEST_CASE("Convert true anomaly to eccentric anomaly" , "[true-to-eccentric-anomaly]")
{
    SECTION("Test elliptical orbits")
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set eccentricities [-].
        const Real eccentricities[3] = {0.146, 0.0, 0.0};

        // Set true anomalies [rad].
        const Real trueAnomalies[3]
            = {82.16 / 180.0 * 3.14159265358979323846,
               160.43 / 180.0 * 3.14159265358979323846,
               0.0};

        // Set expected eccentric anomalies [rad].
        const Real expectedEccentricAnomalies[3]
            = {1.290237398010989,
               2.800031718974503,
               0.0};

        // Loop through all the cases, compute the eccentric anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for (unsigned int i = 0; i < 3; i++)
        {
            // Compute eccentric anomaly [rad].
            const Real computedEllipticalEccentricAnomaly
                = convertTrueAnomalyToEllipticalEccentricAnomaly(
                    trueAnomalies[i], eccentricities[i]);

            const Real computedEccentricAnomaly
                = convertTrueAnomalyToEccentricAnomaly(
                    trueAnomalies[i], eccentricities[i]);

            // Check if computed eccentric anomaly matches the expected value.
            REQUIRE(computedEllipticalEccentricAnomaly
                        == Catch::Approx(expectedEccentricAnomalies[i]).epsilon(
                            std::numeric_limits<Real>::epsilon()));

            REQUIRE(computedEccentricAnomaly
                        == Catch::Approx(expectedEccentricAnomalies[i]).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }
    }

    SECTION("Test general hyperbolic orbits")
    {
        // The benchmark data is obtained from (Fortescue, 2003).

        // Set eccentricities [-].
        const Real eccentricities[1] = {3.0};

        // Set true anomalies [rad].
        const Real trueAnomalies[1] = {0.5291};

        // Set expected eccentric anomalies [rad].
        const Real expectedEccentricAnomalies[1] = {0.3879};

        // Loop through all the cases, compute the eccentric anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for (unsigned int i = 0; i < 1; i++)
        {
            // Compute eccentric anomaly [rad].
            const Real computedHyperbolicEccentricAnomaly
                = convertTrueAnomalyToHyperbolicEccentricAnomaly(
                    trueAnomalies[i], eccentricities[i]);

            const Real computedEccentricAnomaly
                = convertTrueAnomalyToEccentricAnomaly(
                    trueAnomalies[i], eccentricities[i]);

            // Check if computed eccentric anomaly matches the expected value.
            REQUIRE(computedHyperbolicEccentricAnomaly
                        == Catch::Approx(expectedEccentricAnomalies[i]).epsilon(1.0e-5));

            REQUIRE(computedEccentricAnomaly
                        == Catch::Approx(expectedEccentricAnomalies[i]).epsilon(1.0e-5));
        }
    }

    // @TODO: find a way to test assert statements

}

TEST_CASE("Convert eccentric anomaly to mean anomaly" , "[eccentric-to-mean-anomaly]")
{
    SECTION("Test elliptical orbits")
    {
        // The benchmark data is obtained by running ODTBX (NASA, 2012).

        // Set eccentricities [-].
        const Real eccentricities[3] = {0.541, 0.0, 0.0};

        // Set elliptical eccentric anomalies [rad].
        const Real ellipticalEccentricAnomalies[3]
            = {176.09 / 180.0 * 3.14159265358979323846,
               320.12 / 180.0 * 3.14159265358979323846,
               0.0};

        // Set expected mean anomalies [rad].
        const Real expectedMeanAnomalies[3]
            = {3.036459804491048,
               5.587148001484247,
               0.0};

        // Loop through all the cases, compute the mean anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for (unsigned int i = 0; i < 3; i++)
        {
            // Compute mean anomaly [rad].
            const Real computedMeanAnomaly
                = convertEllipticalEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomalies[i], eccentricities[i]);

            const Real computedMeanAnomalyWrapper
                = convertEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomalies[i], eccentricities[i]);

            // Check if computed mean anomaly matches the expected value.
            REQUIRE(computedMeanAnomaly
                        == Catch::Approx(expectedMeanAnomalies[i]).epsilon(
                            std::numeric_limits<Real>::epsilon()));

            REQUIRE(computedMeanAnomalyWrapper
                        == Catch::Approx(expectedMeanAnomalies[i]).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }
    }

    SECTION("Test hyperbolic orbits")
    {
        // The benchmark data is obtained from (Vallado, 2004).

        // Set eccentricities [-].
        const Real eccentricities[1] = {2.4};

        // Set hyperbolic eccentric anomalies [rad].
        const Real hyperbolicEccentricAnomalies[1] = {1.6013761449};

        // Set expected mean anomalies [rad].
        const Real expectedMeanAnomalies[1] = {235.4 / 180.0 * 3.14159265358979323846};

        // Loop through all the cases, compute the mean anomaly and check that the result
        // matches the expected value. Both the direct and wrapper functions are tested.
        for (unsigned int i = 0; i < 1; i++)
        {
            // Compute mean anomaly [rad].
            const Real computedMeanAnomaly
                = convertHyperbolicEccentricAnomalyToMeanAnomaly(
                    hyperbolicEccentricAnomalies[i], eccentricities[i]);

            const Real computedMeanAnomalyWrapper
                = convertEccentricAnomalyToMeanAnomaly(
                    hyperbolicEccentricAnomalies[i], eccentricities[i]);

            // Check if computed mean anomaly matches the expected value.
            REQUIRE(computedMeanAnomaly
                        == Catch::Approx(expectedMeanAnomalies[i]).epsilon(1.0e-10));

            REQUIRE(computedMeanAnomalyWrapper
                        == Catch::Approx(expectedMeanAnomalies[i]).epsilon(1.0e-10));
        }
    }

    // @TODO: find a way to test assert statements

}

TEST_CASE("Convert eccentric anomaly to true anomaly" , "[eccentric-to-true-anomaly]")
{
    const Real pi = 3.14159265358979323846;

    SECTION("Test elliptical orbits")
    {
        SECTION("Test general case")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            const Real eccentricity = 0.639;
            const Real ellipticEccentricAnomaly = 239.45 / 180.0 * pi;
            const Real expectedTrueAnomaly = 3.665218735816221;

            // Compute true anomaly, modulo 2*pi.
            const Real computedTrueAnomaly
                    = convertEllipticalEccentricAnomalyToTrueAnomaly(
                        ellipticEccentricAnomaly, eccentricity) + 2.0 * pi;

            REQUIRE(computedTrueAnomaly
                    == Catch::Approx(expectedTrueAnomaly).epsilon(
                        std::numeric_limits<Real>::epsilon()));
        }

        SECTION("Test wrapper function")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            const Real eccentricity = 0.639;
            const Real ellipticEccentricAnomaly = 239.45 / 180.0 * pi;
            const Real expectedTrueAnomaly = 3.665218735816221;

            // Compute true anomaly, modulo 2*pi.
            const Real computedTrueAnomaly
                    = convertEccentricAnomalyToTrueAnomaly(
                        ellipticEccentricAnomaly, eccentricity) + 2.0 * pi;

            REQUIRE(computedTrueAnomaly
                        == Catch::Approx(expectedTrueAnomaly).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }
    }

    SECTION("Test circular orbits")
    {
        SECTION("Test general case")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            const Real eccentricity = 0.0;
            const Real ellipticEccentricAnomaly = -99.54 / 180.0 * pi;
            const Real expectedTrueAnomaly = 4.545884569744431;

            // Compute true anomaly, modulo 2*pi.
            const Real computedTrueAnomaly
                = convertEllipticalEccentricAnomalyToTrueAnomaly(
                    ellipticEccentricAnomaly, eccentricity) + 2.0 * pi;

            REQUIRE(computedTrueAnomaly
                        == Catch::Approx(expectedTrueAnomaly).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }

        SECTION("Test wrapper function")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            const Real eccentricity = 0.0;
            const Real ellipticEccentricAnomaly = -99.54 / 180.0 * pi;
            const Real expectedTrueAnomaly = 4.545884569744431;

            // Compute true anomaly, modulo 2*pi.
            const Real computedTrueAnomaly
                    = convertEccentricAnomalyToTrueAnomaly(
                        ellipticEccentricAnomaly, eccentricity) + 2.0 * pi;

            REQUIRE(computedTrueAnomaly
                        == Catch::Approx(expectedTrueAnomaly).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }

        SECTION("Test at periapsis (E=0)")
        {
            const Real eccentricity = 0.0;
            const Real ellipticEccentricAnomaly = 0.0;
            const Real expectedTrueAnomaly = 0.0;

            // Compute true anomaly, modulo 2*pi.
            const Real computedTrueAnomaly
                = convertEllipticalEccentricAnomalyToTrueAnomaly(ellipticEccentricAnomaly,
                                                                 eccentricity);

            REQUIRE(computedTrueAnomaly
                        == Catch::Approx(expectedTrueAnomaly).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }

        SECTION("Test wrapper function")
        {
            // The benchmark data is obtained by running ODTBX (NASA, 2012).

            const Real eccentricity = 0.0;
            const Real ellipticEccentricAnomaly = 0.0;
            const Real expectedTrueAnomaly = 0.0;

            // Compute true anomaly, modulo 2*pi.
            const Real computedTrueAnomaly
                    = convertEccentricAnomalyToTrueAnomaly(
                        ellipticEccentricAnomaly, eccentricity);

            REQUIRE(computedTrueAnomaly
                        == Catch::Approx(expectedTrueAnomaly).epsilon(
                            std::numeric_limits<Real>::epsilon()));
        }
    }

    SECTION("Test hyperbolic orbits")
    {
        SECTION("Test general case")
        {
            // The benchmark data is obtained from (Fortescue, 2003).

            const double eccentricity = 3.0;
            const double hyperbolicEccentricAnomaly = 0.3879;
            const double expectedTrueAnomaly = 0.5291;

            // Compute true anomaly.
            const Real computedTrueAnomaly
                = convertHyperbolicEccentricAnomalyToTrueAnomaly(hyperbolicEccentricAnomaly,
                                                                 eccentricity);

            REQUIRE(computedTrueAnomaly == Catch::Approx(expectedTrueAnomaly).epsilon(1.0e-5));
        }

        SECTION("Test wrapper function")
        {
            // The benchmark data is obtained from (Fortescue, 2003).

            const Real eccentricity = 3.0;
            const Real hyperbolicEccentricAnomaly = 0.3879;
            const Real expectedTrueAnomaly = 0.5291;

            // Compute true anomaly.
            const Real computedTrueAnomaly
                = convertEccentricAnomalyToTrueAnomaly(hyperbolicEccentricAnomaly, eccentricity);

            REQUIRE(computedTrueAnomaly == Catch::Approx(expectedTrueAnomaly).epsilon(1.0e-5));
        }
    }
}

TEST_CASE("Compute Kepler function for elliptical orbits" , "[kepler-function]")
{
    // Set elliptical eccentric anomalies [rad].
    const Real ellipticalEccentricAnomalies[3]
            = { 176.09 / 180.0 * 3.14159265358979323846,
                320.12 / 180.0 * 3.14159265358979323846,
                0.0 };

    // Set eccentricities [-].
    const Real eccentricities[3] = {0.541, 0.0, 0.0};

    // Set mean anomalies [rad].
    const Real expectedMeanAnomalies[3]
        = {3.036459804491048,
           5.587148001484247,
           0.0};

    for (unsigned int i = 0; i < 3; ++i)
    {
        REQUIRE(computeEllipticalKeplerFunction(ellipticalEccentricAnomalies[i],
                                                eccentricities[i],
                                                expectedMeanAnomalies[i])
                 < 1.0e-15);
    }
}

TEST_CASE("Compute 1st derivative of Kepler function for elliptical orbits" , "[kepler-function]")
{
    // Note that the test values in this section do not come from literature. They were generated
    // manually.

    // Set elliptical eccentric anomalies [rad].
    const Real ellipticalEccentricAnomalies[3]
            = {0.0,
               2.89735,
               -1.7274};

    // Set eccentricities [-].
    const Real eccentricities[3] = {0.0, 0.3782, 0.79442};

    // Set expected function values.
    const Real expectedFunctionValues[3] = {1.0,
                                            1.3669753060972498,
                                            1.1239011971120707};

    for (unsigned int i = 0; i < 3; ++i)
    {
        REQUIRE(computeFirstDerivativeEllipticalKeplerFunction(ellipticalEccentricAnomalies[i],
                                                               eccentricities[i])
                    == Catch::Approx(expectedFunctionValues[i]).epsilon(
                        std::numeric_limits<Real>::epsilon()));
    }
}

TEST_CASE("Convert elliptical mean anomaly to eccentric anomaly for elliptical orbits" ,
           "[mean-to-eccentric-anomaly]")
{
    const Real pi = 3.14159265358979323846;

    SECTION("Test circular orbits")
    {
        const Real eccentricity = 0.0;
        const Real ellipticalMeanAnomaly = 1.0472;
        const Real expectedEccentricAnomaly = 1.0472;

        // Compute eccentric anomaly, modulo 2*pi.
        const Real computedEccentricAnomaly
                = convertEllipticalMeanAnomalyToEccentricAnomaly<Real, Integer>(
                    eccentricity, ellipticalMeanAnomaly);

        REQUIRE(computedEccentricAnomaly
                    == Catch::Approx(expectedEccentricAnomaly).epsilon(
                        std::numeric_limits<Real>::epsilon()));
    }

    SECTION("Test arbitrary elliptical orbits using GTOP data.")
    {
        // The benchmark data is obtained from GTOP.
        const Real eccentricities[3] = {0.01671, 0.43582, 0.78514};
        const Real ellipticalMeanAnomalies[3] = {60.0 / 180.0 * pi,
                                                 90.0 / 180.0 * pi,
                                                 120.0 / 180.0 * pi};
        const Real expectedEccentricAnomalies[3] = {1.06178920406832,
                                                    1.97200731113253,
                                                    2.5392410896466};

        for (unsigned int i = 0; i < 3; ++i)
        {
            const Real computedEccentricAnomaly
                = convertEllipticalMeanAnomalyToEccentricAnomaly<Real, Integer>(
                    eccentricities[i], ellipticalMeanAnomalies[i]);

            Real computedEccentricAnomalyShifted = std::fmod(computedEccentricAnomaly, 2.0 * pi);
            if (computedEccentricAnomalyShifted < 0.0)
            {
                computedEccentricAnomalyShifted += 2.0 * pi;
            }

            REQUIRE(computedEccentricAnomalyShifted
                        == Catch::Approx(expectedEccentricAnomalies[i]).epsilon(
                            10.0 * std::numeric_limits<Real>::epsilon()));
        }
    }

    SECTION("Test arbitrary elliptical orbits using PyKEP data.")
    {
        // The benchmark data is obtained from PyKEP.
        const Real eccentricities[4] = {0.5132, 0.0, 0.223, 0.991};
        const Real ellipticalMeanAnomalies[4] = {2.5746, 4.47712, -3.39915, 0.5571};
        const Real expectedEccentricAnomalies[4] = {2.76387035891018,
                                                       4.47712,
                                                       -3.35247173243822 + 2.0 * pi,
                                                       1.54783886054501};

        for (unsigned int i = 0; i < 4; ++i)
        {
            const Real computedEccentricAnomaly
                = convertEllipticalMeanAnomalyToEccentricAnomaly<Real, Integer>(
                    eccentricities[i], ellipticalMeanAnomalies[i]);

            Real computedEccentricAnomalyShifted = std::fmod(computedEccentricAnomaly, 2.0 * pi);
            if (computedEccentricAnomalyShifted <0.0)
            {
                computedEccentricAnomalyShifted += 2.0 * pi;
            }

            REQUIRE(computedEccentricAnomalyShifted
                        == Catch::Approx(expectedEccentricAnomalies[i]).epsilon(
                            10.0 * std::numeric_limits<Real>::epsilon()));
        }
    }
}

} // namespace tests
} // namespace astro

/*!
 * References
 *  Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition, Wiley, England, 2003.
 *  Advanced Concepts Team, ESA. Global Trajectory Optimization Problems Database (GTOP),
 *    https://www.esa.int/gsp/ACT/projects/gtop/gtop.html, last accessed: 2nd January, 2019.
 *  NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *    Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *    31st January, 2012.
 *  Noomen, R. Space Mission Design: Basics, ae4-878 (Mission Geometry & Orbit Design), Delft
 *    University of Technology, The Netherlands, 2014.
 *  pykep Development Team, ESA. pykep, https://esa.github.io/pykep, last accessed: 2nd January,
 *    2019.
 */
