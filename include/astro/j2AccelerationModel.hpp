/*
 * Copyright (c) 2014-2025 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

#include <cmath>

namespace astro
{

//! Compute gravitational acceleration due to J2.
/*!
 * Compute gravitational acceleration at given position vector subject to an irregular gravity
 * field. The acceleration due to the J2-coefficient is given by (Melman, 2012):
 *
 * \f{eqnarray*}{
 *      {a}_{gravity,x} &=& -\mu_{2}*\frac{x_{1}-x_{2}}{r^{3}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (1-5*\hat{z}^2)  \\
 *      {a}_{gravity,y} &=& -\mu_{2}*\frac{y_{1}-y_{2}}{r^{3}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (1-5*\hat{z}^2)  \\
 *      {a}_{gravity,z} &=& \frac{-\mu_{2}}{r^{2}} * \frac{3}{2}*J_{2}
 *                          * \left(\frac{R}{r}\right)^{2} * (3-5*\hat{z}^2) * \hat{z}  \\
 * \f}
 *
 * where \f$\mu_{2}\f$ is the gravitational parameter of the central body,
 * \f$\hat{z} = \frac{z}{r}\f$, \f$x\f$, \f$y\f$ and \f$z\f$ are the Cartesian position components,
 * \f$r\f$ is the radial position, and \f$J_{2}\f$ is the second zonal coefficient of the gravity
 * field.
 *
 * The positions and accelerations are given with respect to an inertial (barycentric) reference
 * frame.
 *
 * @tparam Real                       Real type
 * @tparam Vector                     Vector type
 * @param[in] position                Position vector of body subject to J2-acceleration  [m]
 * @param[in] gravitationalParameter  Gravitational parameter of central body             [m^3 s^-2]
 * @param[in] equatorialRadius        Equatorial radius of central body, in formulation
 *                                    of spherical harmonics expansion                    [m]
 * @param[in] j2Coefficient           Unnormalized J2-coefficient of spherical harmonics
 *                                    expansion                                           [-]
 * @return                            J2 gravitational acceleration                       [m s^-2]
 */
template <typename Real, typename Vector3>
Vector3 computeJ2Acceleration(const Real     gravitationalParameter,
                              const Vector3& position,
                              const Real     equatorialRadius,
                              const Real     j2Coefficient)
{
    Vector3 acceleration = position;

    const Real positionNormSquared = position[0] * position[0]
                                     + position[1] * position[1]
                                     + position[2] * position[2];
    const Real positionNorm = std::sqrt(positionNormSquared);

    const Real scaledZSquared = position[2] * position[2] / positionNormSquared;

    const Real preMultiplier = -gravitationalParameter
                                / (positionNormSquared * positionNormSquared * positionNorm)
                                * 1.5 * j2Coefficient * equatorialRadius * equatorialRadius;

    acceleration[0] = preMultiplier * position[0] * (1.0 - 5.0 * scaledZSquared);
    acceleration[1] = preMultiplier * position[1] * (1.0 - 5.0 * scaledZSquared);
    acceleration[2] = preMultiplier * position[2] * (3.0 - 5.0 * scaledZSquared);

    return acceleration;
}

} // namespace astro
