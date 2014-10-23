/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef SAM_TWO_BODY_METHODS_HPP
#define SAM_TWO_BODY_METHODS_HPP

#include <cmath>
#include <stdexcept>
 
#include "SAM/constants.hpp"

namespace sam
{

//! Compute mean motion.
/*!
 * Computes the two-body mean motion of an orbiting body that follows a conic section
 * (Kepler orbit). The mass of the orbiting body is set to that of a test particle by default.
 *
 * The two-body mean motion \f$n\f$ is computed as follows:
 * \f[
 *      n = \sqrt{\frac{\mu}{a^{3}}}
 * \f]
 * where \f$\mu\f$ is the sum of the gravitational parameters of the two bodies (for the default
 * case of an orbiting test particle, this is the gravitational parameter of the central body) and 
 * \f$a\f$ is the semi-major axis of the Kepler orbit. 
 *
 * @tparam Real                                Real type  
 * @param  semiMajorAxis                       Semi-major axis of Kepler orbit         [m]
 * @param  gravitationalParameterOfCentralBody Gravitational parameter of central body [m^3 s^-2]
 * @param  massOfOrbitingBody                  Mass of orbiting body                   [kg]
 * @return                                     Two-body mean motion                    [rad/s]
 */
template< typename Real > 
Real computeKeplerMeanMotion( const Real semiMajorAxis,
                              const Real gravitationalParameterOfCentralBody,
                              const Real massOfOrbitingBody = 0.0 );

//! Compute circular velocity.
/*!
 * Computes circular orbital velocity \f$V_{c}\f$ for a (massless) body in a Kepler orbit as 
 * follows:
 * 
 * \f[
 *      V_{c} = \sqrt{\frac{2\mu}{a}}
 * \f]
 * where \f$\mu\f$ is the gravitational parameter of the central body) and \f$a\f$ is the 
 * semi-major axis of the Kepler orbit.
 *
 * If the semi-major axis is zero (to machine precision), an exception is thrown. 
 *
 * @tparam Real                                Real type  
 * @param  semiMajorAxis                       Semi-major axis of Kepler orbit         [m]
 * @param  gravitationalParameterOfCentralBody Gravitational parameter of central body [m^3 s^-2] 
 * @return                                     Circular velocity
 */
template< typename Real > 
Real computeCircularVelocity( const Real semiMajorAxis, 
                              const Real gravitationalParameterOfCentralBody );
//! Compute mean motion.
template< typename Real > 
Real computeKeplerMeanMotion( const Real semiMajorAxis,
                              const Real gravitationalParameterOfCentralBody,
                              const Real massOfOrbitingBody )
{
    return std::sqrt( ( ( SAM_GRAVITATIONAL_CONSTANT * massOfOrbitingBody )
                        + gravitationalParameterOfCentralBody ) 
                      / ( semiMajorAxis * semiMajorAxis * semiMajorAxis ) );
}

//! Compute circular velocity.
template< typename Real > 
Real computeCircularVelocity( const Real semiMajorAxis, 
                              const Real gravitationalParameterOfCentralBody )
{
    if ( std::fabs( semiMajorAxis ) < std::numeric_limits< Real >::epsilon( ) )
    {
        throw std::runtime_error( "ERROR: Semi-major axis is zero!" );
    }

    return std::sqrt( gravitationalParameterOfCentralBody / semiMajorAxis );
}

} // namespace sam

#endif // SAM_TWO_BODY_METHODS_HPP
