/*
 * Copyright (c) 2014-2018 Kartik Kumar (me@kartikkumar.com)
 * Copyright (c) 2014-2015 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2015 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2015 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_CENTRAL_BODY_ACCELERATION_MODEL_HPP
#define ASTRO_CENTRAL_BODY_ACCELERATION_MODEL_HPP

#include <cmath>

namespace astro
{

//! Compute the acceleration of a point mass body orbiting a uniform central body.
/*!
 * Computes the acceleration of a point mass body orbiting a uniform central body based on Newton's
 * second law and his law of gravitation.
 *
 * The expression of this acceleration is based on the two-body equation of motion, which represents
 * the relative equation of motion of a body as it orbits the central body.
 * \f[
 *      \vec{a}_{gravity}=-\frac{\mu}{r^{3}}\vec{r}
 * \f]
 *
 * where \f$\mu\f$ is a gravitational parameter (e.g., \f$\mu=GM=398,600.5 km^{3}s^-{2}\f$
 * in the case of the Earth), \f$r\f$ is the distance from the origin of the reference frame to the
 * body (i.e., the magnitude of the position vector \f$\vec{r}\f$), and \f$\vec{r}\f$ is the
 * position vector of the body relative to the origin of the reference frame.
 *
 * @tparam     Real                           Real type
 * @tparam     Vector3                        3-vector type
 * @param[in]  gravitationalParameter         Gravitational parameter of central body  [km^3 s^-2]
 * @param[in]  position                       Position vector of the orbiting body     [km]
 * @return                                    Acceleration vector                      [km s^-2]
 */
template< typename Real, typename Vector3 >
Vector3 computeCentralBodyAcceleration( const Real     gravitationalParameter,
                                        const Vector3& position )
{
    Vector3 acceleration = position;

    const Real positionNorm = std::sqrt( position[ 0 ] * position[ 0 ]
                                         + position[ 1 ] * position[ 1 ]
                                         + position[ 2 ] * position[ 2 ] );
    const Real preMultiplier = -gravitationalParameter
                                 / ( positionNorm * positionNorm * positionNorm );

    acceleration[ 0 ] = preMultiplier * position[ 0 ];
    acceleration[ 1 ] = preMultiplier * position[ 1 ];
    acceleration[ 2 ] = preMultiplier * position[ 2 ];

    return acceleration;
}

} // namespace astro

#endif // ASTRO_CENTRAL_BODY_ACCELERATION_MODEL_HPP
