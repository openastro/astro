/*
 * Copyright (c) 2014-2015 Juan Manuel Romero Martin (romero.martin.jm@gmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_SOLAR_PRESSURE_ACCELERATION_MODEL_HPP
#define ASTRO_SOLAR_PRESSURE_ACCELERATION_MODEL_HPP

#include <SML/linearAlgebra.hpp>

namespace astro
{

//! Compute solar radiation pressure acceleration on a cannonball.
/*!
 * Compute solar radiation pressure acceleration on a cannonball model. The model for the solar 
 * radiation pressure acceleration is given by
 *
 * \f[
 *      a_{solar pressure} =  - P  C_{R}  \frac{A_{S/C}}{m} \vec{u}
 * \f]
 * 
 * where \f$P\f$ is the solar radiotion pressure, \f$C_{R}\f$ is the radiation pressure coefficient, 
 * \f$ A_{S/C} \f$ is the absorbing area of teh satellite, which is \f$ \pi R^2 \f$ for the 
 * connonball, \f$ m \f$ is the mass of the satellite, and \f$ \vec{u} \f$ is the unit vector 
 * pointing from the satellite toward the sun.
 *
 * @param[in] radiationPressure        Solar Radiation Presure                     [N m^-2]
 * @param[in] radiationPressureCoef    Radiation Pressure Coefficiente             [adim]
 * @param[in] vectorToSource           Unit vector pointing from S/C to sun (3x1)  [adim]                     
 * @param[in] area                     Absorbing Area of S/C                       [m^2]
 * @param[in] mass                     Mass of the S/C                             [kg]    
 * return The Solar radiation pressure                                             [m^-2]
*/
template< typename Real, typename Vector3 >
Vector3 computeSolarRadiationPressureAcceleration(  const Real     radiationPressure, 
                                                    const Real     radiationPressureCoef, 
                                                    const Vector3& vectorToSource,                                
                                                    const Real     area,
                                                    const Real     mass );

//! Compute solar radiation pressure acceleration on a cannonball.
template< typename Real, typename Vector3 >
Vector3 computeSolarRadiationPressureAcceleration(  const Real     radiationPressure, 
                                                    const Real     radiationPressureCoefficient, 
                                                    const Vector3& vectorToSource,                                
                                                    const Real     area,
                                                    const Real     mass )
{
  
    // Instance the output accelation 
    Vector3 radiationPressureAcceleration = vectorToSource; //@todo: [JRMM] Change the initialization

    // Compute the solar pressure radiation magnitude
    const Real radPressureMagnitude = -radiationPressure * radiationPressureCoefficient * area / mass;

    // Compute the radiation pressure acceleration vector.
    radiationPressureAcceleration[ 0 ] = radPressureMagnitude * vectorToSource[ 0 ];
    radiationPressureAcceleration[ 1 ] = radPressureMagnitude * vectorToSource[ 1 ];
    radiationPressureAcceleration[ 2 ] = radPressureMagnitude * vectorToSource[ 2 ];

    return radiationPressureAcceleration;
};

} // namespace astro

#endif // ASTRO_SOLAR_PRESSURE_ACCELERATION_MODEL_HPP
