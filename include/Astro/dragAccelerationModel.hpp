/*
 * Copyright (c) 2014-2015 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_DRAG_ACCELERATION_MODEL_HPP
#define ASTRO_DRAG_ACCELERATION_MODEL_HPP

#include <SML/linearAlgebra.hpp>

namespace astro
{

//! Compute drag acceleration on a cannonball.
/*!
 * Computes drag acceleration using a cannonball model. The model for the acceleration is given by
 *
 * \f[
 *      a_{drag} = \frac{1}{2} \frac{C_{d}}{m} \rho S V \vec{V}
 * \f]
 * 
 * where \f$a_{drag}\f$ is the drag acceleration, \f$C_{d}\f$ is the drag coefficient, \f$m\f$ is
 * the mass, \f$\rho\f$ is the atmospheric density, \f$\vec{V}\f$ is the relative velocity with 
 * respect to the body-fixed frame and \f$S\f$ is the drag area, that is, the projected area
 * of the object perpendicular to \f$\vec{V}\f$.
 *
 * @param[in]  dragCoefficient      Drag coefficient               [adim]
 * @param[in]  atmosphericDensity   Atmospheric density            [kg m^-3]
 * @param[in]  velocity             Velocity vector (3x1)          [m s^-1]
 * @param[in]  dragArea             Drag area                      [m^2]
 * @param[in]  mass                 Mass                           [kg]
 * @return                          Drag acceleration vector (3x1) [m s^-2]                                                      
 */
template< typename Real, typename Vector3 >
Vector3 computeDragAcceleration( const Real     dragCoefficient, 
                                 const Real     atmosphericDensity, 
                                 const Vector3& velocity,
                                 const Real     dragArea,
                                 const Real     mass ); 

//! Compute drag acceleration on a cannonball.
template< typename Real, typename Vector3 >
Vector3 computeDragAcceleration( const Real     dragCoefficient, 
                                 const Real     atmosphericDensity, 
                                 const Vector3& velocity,
                                 const Real     dragArea,
                                 const Real     mass )
{
    Vector3 dragAcceleration = velocity;

    // Compute the squared norm of the velocity.
    const Real normVelocity = sml::norm< Real >( velocity );

    // Compute a premultiplier so that it does not have to be written several times.
    const Real preMultiplier = 0.5 * dragCoefficient * atmosphericDensity 
                               * dragArea * normVelocity / mass;

    // Compute the drag acceleration vector.
    dragAcceleration[ 0 ] = preMultiplier * velocity[ 0 ];
    dragAcceleration[ 1 ] = preMultiplier * velocity[ 1 ];
    dragAcceleration[ 2 ] = preMultiplier * velocity[ 2 ];

    return dragAcceleration;
}

} // namespace astro

#endif // ASTRO_DRAG_ACCELERATION_MODEL_HPP
