/*
 * Copyright (c) 2014-2015 Kartik Kumar, Dinamica Srl
 * Copyright (c) 2014-2015 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2015 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2015 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_CENTRAL_BODY_ACCELERATION_MODEL_HPP
#define ASTRO_CENTRAL_BODY_ACCELERATION_MODEL_HPP

#include <cmath>

#include <sml/sml.hpp>

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
 * @param[in]  positionVector                 Position vector of the body              [km]
 * @return                                    Acceleration vector                      [km s^-2]
 */
template< typename Real, typename Vector3 >
Vector3 computeCentralBodyAcceleration( const Real     gravitationalParameter,
                                        const Vector3& positionVector )
{
    Vector3 bodyAcceleration = positionVector;

    // Compute the magnitude of the position vector.
    const double positionNorm = sml::norm< double >( positionVector );

    const double preMultiplier = -gravitationalParameter
                                 / ( positionNorm * positionNorm * positionNorm );

    // Compute the acceleration.
    bodyAcceleration[ 0 ] = preMultiplier * positionVector[ 0 ];
    bodyAcceleration[ 1 ] = preMultiplier * positionVector[ 1 ];
    bodyAcceleration[ 2 ] = preMultiplier * positionVector[ 2 ];

    return bodyAcceleration;
}

} // namespace astro

#endif // ASTRO_CENTRAL_BODY_ACCELERATION_MODEL_HPP
