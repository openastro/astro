/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Abhishek Agrawal, Delft University of Technology
 *                                          (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef J2_GRAVITY_HPP
#define J2_GRAVITY_HPP

#include <cmath>

#include <Astro/stateVectorIndices.hpp>

namespace astro
{

//! Compute orbit-averaged rate of change in the orbital elements due to J2 perturbation.
/*!
 * Computes the first order, orbit-averaged rate of change in longitude of ascending node and
 * argument of periapsis due to the J2 perturbation.
 *
 * @tparam      Real                Real type
 * @tparam      Vector6             Vector of 6 elements
 *
 * @param[in]   keplerianElements   Vector containing Keplerian elements
 *                                  N.B.: Order of elements and units is important!
 *                                  keplerianElements( 0 ) = semiMajorAxis                [km]  <br>
 *                                  keplerianElements( 1 ) = eccentricity                 [-]   <br>
 *                                  keplerianElements( 2 ) = inclination                  [rad] <br>
 *                                  keplerianElements( 3 ) = argument of periapsis        [rad] <br>
 *                                  keplerianElements( 4 ) = longitude of ascending node  [rad] <br>
 *                                  keplerianElements( 5 ) = true anomaly                 [rad]
 * @param[in]  meanMotion           meanMotion of the orbiting object [deg/day]
 * @param[in]  earthMeanRadius      Earth mean radius [km]
 *
 * @param[out] longitudeAscendingNodeDot Stores rate of change in longitude of ascending node
 *                                       due to J2 perturbation [deg/day]
 * @param[out] argumentOfPeriapsisDot    Stores rate of change in argument of periapsis
 *                                       due to J2 perturbation [deg/day]
 */
 template< typename Real, typename Vector6 >
 void computeFirstOrderAveragedEffectJ2Perturbation(
    const Vector6& keplerianElements,
    const Real     meanMotion,
    const Real     earthMeanRadius,
    const Real&    longitudeAscendingNodeDot,
    const Real&    argumentOfPeriapsisDot )
 {
    const Real J2 = 0.00108263;
    const Real semiMajorAxis = keplerianElements[ astro::semiMajorAxisIndex ];
    const Real inclination   = keplerianElements[ astro::inclinationIndex ];
    const Real eccentricity  = keplerianElements[ astro::eccentricityIndex ];

    longitudeAscendingNodeDot
        = -1.5 * meanMotion * J2
          * ( earthMeanRadius / semiMajorAxis ) * ( earthMeanRadius / semiMajorAxis )
          * std::cos( inclination )
          / ( ( 1 - eccentricity * eccentricity ) * ( 1 - eccentricity * eccentricity ) );

    argumentOfPeriapsisDot
        = 0.75 * meanMotion * J2
          * ( earthMeanRadius / semiMajorAxis ) * ( earthMeanRadius / semiMajorAxis )
          * ( 4 - 5 * std::sin( inclination ) * std::sin( inclination ) )
          / ( ( 1 - eccentricity * eccentricity ) * ( 1 - eccentricity * eccentricity ) );
 }

} // namespace_astro

#endif // J2_GRAVITY_HPP
