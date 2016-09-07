/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2016 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2016 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_SOLAR_PRESSURE_ACCELERATION_MODEL_HPP
#define ASTRO_SOLAR_PRESSURE_ACCELERATION_MODEL_HPP

#include <sml/linearAlgebra.hpp>

namespace astro
{

//! Compute solar radiation pressure acceleration on a cannonball.
/*!
 * Compute solar radiation pressure acceleration on a cannonball model. The model for the solar
 * radiation pressure acceleration is given by:
 *
 * \f[
 *      a_{solar pressure} =  - P  C_{R}  \frac{A_{S/C}}{m} \vec{u}
 * \f]
 *
 * where \f$P\f$ is the solar radiotion pressure, \f$C_{R}\f$ is the radiation pressure
 * coefficient, \f$ A_{S/C} \f$ is the absorbing area of teh satellite, which is \f$ \pi R^2 \f$
 * for the cannonball, \f$ m \f$ is the mass of the satellite, and \f$ \vec{u} \f$ is the unit
 * vector pointing from the satellite toward the Sun.
 *
 * @tparam    Real                            Floating-point type
 * @tparam    Vector3                         3-vector type
 * @param[in] radiationPressure               Solar radiation pressure                    [N m^-2]
 * @param[in] radiationPressureCoefficient    Radiation pressure coefficient              [adim]
 * @param[in] vectorToSource                  Unit vector pointing from S/C to sun (3x1)  [adim]
 * @param[in] area                            Absorbing area of S/C                       [m^2]
 * @param[in] mass                            Mass of S/C                                 [kg]
 * @return                                    Solar radiation pressure acceleration       [m s^-2]
 */
template< typename Real, typename Vector3 >
Vector3 computeSolarRadiationPressureAcceleration( const Real     radiationPressure,
                                                   const Real     radiationPressureCoefficient,
                                                   const Vector3& vectorToSource,
                                                   const Real     area,
                                                   const Real     mass )
{

    // Initialize the output acceleration.
    //@todo: [JRMM] Change the initialization
    Vector3 radiationPressureAcceleration = vectorToSource;

    // Compute the solar pressure radiation magnitude.
    const Real radiationPressureMagnitude = -radiationPressure * radiationPressureCoefficient
                                            * area / mass;

    // Compute the radiation pressure acceleration vector.
    radiationPressureAcceleration[ 0 ] = radiationPressureMagnitude * vectorToSource[ 0 ];
    radiationPressureAcceleration[ 1 ] = radiationPressureMagnitude * vectorToSource[ 1 ];
    radiationPressureAcceleration[ 2 ] = radiationPressureMagnitude * vectorToSource[ 2 ];

    return radiationPressureAcceleration;
};

} // namespace astro

#endif // ASTRO_SOLAR_PRESSURE_ACCELERATION_MODEL_HPP
