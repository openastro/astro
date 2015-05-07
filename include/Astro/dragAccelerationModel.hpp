/*
 * Copyright (c) 2014-2015 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <SML/linearAlgebra.hpp>

namespace astro
{

//! Compute drag acceleration on a cannonball.
/*!
 * Computes drag acceleration using a cannonball model. The model for the acceleration is given by
 *
 * \f[
 *      a_{drag} = \frac{1}{2} \frac{C_{d}}{m} \rho V S \vec{V}
 * \f]
 * 
 * where \f$a_{drag}\f$ is the drag acceleration, \f$C_{d}\f$ is the drag coefficient, \f$m\f$ is
 * the mass, \f$\rho\f$ is the atmospheric density, \$f\vec{V}$f is the relative velocity with 
 * respect to the body-fixed frame  and \f$S\f$ is the drag area, that is, the projected area
 * of the object perpendicular to \$f\vec{V}$f.
 *
 * @param  dragCoefficient                     Drag coefficient        [adim]
 * @param  atmosphericDensity                  Atmospheric Density     [kg * m^-3]
 * @param  velocity                            Velocity vector 3x3     [m s^-1]
 * @param  dragArea                            Drag area               [m^2]
 * @param  mass                                Mass                    [kg]
 * @return                                        
 */
 template< typename Real, typename Vector3 >
 Vector3 computeDragAcceleration( Real    dragCoefficient, 
                                  Real    atmosphericDensity, 
                                  Vector3 velocity,
                                  Real    dragArea,
                                  Real    mass ); 


 //! Compute drag acceleration on a cannonball.
 template< typename Real, typename Vector3 >
 Vector3 computeDragAcceleration( Real    dragCoefficient, 
                                  Real    atmosphericDensity, 
                                  Vector3 velocity,
                                  Real    dragArea,
                                  Real    mass )
 {
    Vector3 dragAcceleration( 3 );

    // Compute the squared norm of the velocity
    Real normVelocity = sml::norm< Real >( velocity );

    // Compute a premultiplier so that it does not have to be written several times
    Real preMultiplier = 0.5 * dragCoefficient * atmosphericDensity * dragArea * normVelocity/mass;

    // Compute the drag acceleration vector
    dragAcceleration[ 0 ] = preMultiplier * velocity[ 0 ];
    dragAcceleration[ 1 ] = preMultiplier * velocity[ 1 ];
    dragAcceleration[ 2 ] = preMultiplier * velocity[ 2 ];

    return dragAcceleration;
 }

} // namespace astro