/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP
#define SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP

#include <cmath>
 
#include <SML/sml.hpp>

#include "SAM/constants.hpp"

namespace sam
{

//! Cartesian element array indices.
enum CartesianElementIndices
{
    xPositionIndex = 0,
    yPositionIndex = 1,
    zPositionIndex = 2,
    xVelocityIndex = 3,
    yVelocityIndex = 4,
    zVelocityIndex = 5
};

//! Cartesian element array indices.
enum KeplerianElementIndices
{
    semiMajorAxisIndex = 0,
    semiLatusRectumIndex = 0,    
    eccentricityIndex = 1,
    inclinationIndex = 2,
    argumentOfPeriapsisIndex = 3,
    longitudeOfAscendingNodeIndex = 4,
    trueAnomalyIndex = 5,
    meanAnomalyIndex = 5
};

//! Convert Cartesian elements to Keplerian elements.
/*!
 * Converts a given set of Cartesian elements (position, velocity) to classical (osculating) 
 * Keplerian elements. See Chobotov (2006) for a full derivation of the conversion.
 *
 * The tolerance is given a default value. It should not be changed unless required for specific 
 * scenarios. Below this tolerance value for eccentricity and inclination, the orbit is considered 
 * to be a limit case. Essentially, special solutions are then used for parabolic, circular 
 * inclined, non-circular equatorial, and circular equatorial orbits. These special solutions are 
 * required because of singularities in the classical Keplerian elements. If  high precision is 
 * required near these singularities, users are encouraged to consider using other elements, such 
 * as Modified Equinoctial Elements (MEE). It should be noted that MEE also suffer from 
 * singularities, but not for zero eccentricity and inclination.
 * 
 * WARNING: If eccentricity is 1.0 within 1.0e-15, keplerianElements( 0 ) = semi-latus rectum,
 *          since the orbit is parabolic.
 * WARNING: If eccentricity is 0.0 within 1.0e-15, argument of periapsis is set to 0.0, since the
 *          orbit is circular.
 * WARNING: If inclination is 0.0 within 1.0e-15, longitude of ascending node is set to 0.0, since 
 *          the orbit is equatorial.
 *
 * @tparam  Real              Real type
 * @tparam  Vector6           6-Vector type
 * @tparam  Vector3           3-Vector type
 * @param   cartesianElements Vector containing Cartesian elements. 
 *                            N.B.: Order of elements is important!
 *                            cartesianElements( 0 ) = x-position coordinate [m]
 *                            cartesianElements( 1 ) = y-position coordinate [m]
 *                            cartesianElements( 2 ) = z-position coordinate [m]
 *                            cartesianElements( 3 ) = x-velocity coordinate [m/s]
 *                            cartesianElements( 4 ) = y-velocity coordinate [m/s]
 *                            cartesianElements( 5 ) = z-velocity coordinate [m/s]
 * @param   gravitationalParameter Gravitational parameter of central body [m^3 s^-2]
 * @param   tolerance         Tolerance used to check for limit cases 
 *                            (zero eccentricity, inclination)
 * @return  Converted vector of Keplerian elements.
 *          N.B.: Order of elements is important! 
 *          keplerianElements( 0 ) = semiMajorAxis                [m]
 *          keplerianElements( 1 ) = eccentricity                 [-]
 *          keplerianElements( 2 ) = inclination                  [rad]
 *          keplerianElements( 3 ) = argument of periapsis        [rad]
 *          keplerianElements( 4 ) = longitude of ascending node  [rad]
 *          keplerianElements( 5 ) = true anomaly                 [rad]
 */
template< typename Real, typename Vector6, typename Vector3 >
Vector6 convertCartesianToKeplerianElements(
    const Vector6& cartesianElements, const Real gravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( ) )
{
    Vector6 keplerianElements( 6 );

    // Set position and velocity vectors.
    Vector3 position( 3 );
    position[ 0 ] = cartesianElements[ 0 ];
    position[ 1 ] = cartesianElements[ 1 ];
    position[ 2 ] = cartesianElements[ 2 ];

    Vector3 velocity( 3 );
    velocity[ 0 ] = cartesianElements[ 3 ];
    velocity[ 1 ] = cartesianElements[ 4 ];
    velocity[ 2 ] = cartesianElements[ 5 ];     

    // Compute orbital angular momentum vector.
    const Vector3 angularMomentum( sml::cross( position, velocity ) );

    // Compute semi-latus rectum.
    const Real semiLatusRectum 
        = sml::squaredNorm< Real >( angularMomentum ) / gravitationalParameter;

    // Compute unit vector to ascending node.
    Vector3 ascendingNodeUnitVector 
        = sml::normalize< Real >( 
            sml::cross( sml::getZUnitVector< Vector3 >( ), 
                        sml::normalize< Real >( angularMomentum ) ) );        

    // Compute eccentricity vector.
    Vector3 eccentricityVector 
        = sml::add( sml::multiply( sml::cross( velocity, angularMomentum ), 
                                   1.0 / gravitationalParameter ),
                    sml::multiply( sml::normalize< Real >( position ), -1.0 ) );

    // Store eccentricity.
    keplerianElements[ eccentricityIndex ] = sml::norm< Real >( eccentricityVector );        

    // Compute and store semi-major axis.
    // Check if orbit is parabolic. If it is, store the semi-latus rectum instead of the
    // semi-major axis.
    if ( std::fabs( keplerianElements[ eccentricityIndex ] - 1.0 ) < tolerance )
    {
        keplerianElements[ semiLatusRectumIndex ] = semiLatusRectum;
    }

    // Else the orbit is either elliptical or hyperbolic, so store the semi-major axis.
    else
    {
        keplerianElements[ semiMajorAxisIndex ] 
            = semiLatusRectum / ( 1.0 - keplerianElements[ eccentricityIndex ]
                                  * keplerianElements[ eccentricityIndex ] );
    }

    // Compute and store inclination.
    keplerianElements[ inclinationIndex ]
        = std::acos( angularMomentum[ zPositionIndex ] 
                     / sml::norm< Real >( angularMomentum ) );    

    // Compute and store longitude of ascending node.
    // Define the quadrant condition for the argument of perigee.
    Real argumentOfPeriapsisQuandrantCondition = eccentricityVector[ zPositionIndex ];

    // Check if the orbit is equatorial. If it is, set the vector to the line of nodes to the
    // x-axis.
    if ( std::fabs( keplerianElements[ inclinationIndex ] ) < tolerance )
    {
        ascendingNodeUnitVector = sml::getXUnitVector< Vector3 >( );

        // If the orbit is equatorial, eccentricityVector_z is zero, therefore the quadrant
        // condition is taken to be the y-component, eccentricityVector_y.
        argumentOfPeriapsisQuandrantCondition = eccentricityVector[ yPositionIndex ];
    }

    // Compute and store the resulting longitude of ascending node.
    keplerianElements[ longitudeOfAscendingNodeIndex ] 
        = std::acos( ascendingNodeUnitVector[ xPositionIndex] );

    // Check if the quandrant is correct.
    if ( ascendingNodeUnitVector[ yPositionIndex ] < 0.0 )
    {
        keplerianElements[ longitudeOfAscendingNodeIndex ] 
            = 2.0 * sml::SML_PI - keplerianElements[ longitudeOfAscendingNodeIndex ];
    }

    // Compute and store argument of periapsis.
    // Define the quadrant condition for the true anomaly.
    Real trueAnomalyQuandrantCondition = sml::dot< Real >( position, velocity );

    // Check if the orbit is circular. If it is, set the eccentricity vector to unit vector
    // pointing to the ascending node, i.e. set the argument of periapsis to zero.
    if ( std::fabs( keplerianElements[ eccentricityIndex ] ) < tolerance )
    {
        eccentricityVector = ascendingNodeUnitVector;

        keplerianElements[ argumentOfPeriapsisIndex ] = 0.0;

        // Check if orbit is also equatorial and set true anomaly quandrant check condition
        // accordingly.
        if ( ascendingNodeUnitVector == sml::getXUnitVector< Vector3 >( ) )
        {
            // If the orbit is circular, dot( position, velocity ) = 0, therefore this value
            // cannot be used as a quadrant condition. Moreover, if the orbit is equatorial,
            // position_z is also zero and therefore the quadrant condition is taken to be the
            // y-component, position_y.
            trueAnomalyQuandrantCondition = position[ yPositionIndex ];
        }

        else
        {
            // If the orbit is circular, dot( position, velocity ) = 0, therefore the quadrant
            // condition is taken to be the z-component of the position, position_z.
            trueAnomalyQuandrantCondition = position[ zPositionIndex ];
        }
    }

    // Else, compute the argument of periapsis as the angle between the eccentricity vector and
    // the unit vector to the ascending node.
    else
    {
        keplerianElements[ argumentOfPeriapsisIndex ]
            = std::acos( sml::dot< Real >( sml::normalize< Real >( eccentricityVector ), ascendingNodeUnitVector ) );

        // Check if the quadrant is correct.
        if ( argumentOfPeriapsisQuandrantCondition < 0.0 )
        {
            keplerianElements[ argumentOfPeriapsisIndex ] 
                = 2.0 * sml::SML_PI - keplerianElements[ argumentOfPeriapsisIndex ];
        }
    }    

    // Compute dot-product of position and eccentricity vectors.
    Real dotProductPositionAndEccentricityVectors
        = sml::dot< Real >( sml::normalize< Real >( position ), 
                            sml::normalize< Real >(eccentricityVector ) );

    // Check if the dot-product is one of the limiting cases: 0.0 or 1.0
    // (within prescribed tolerance).
    if ( std::fabs( 1.0 - dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors = 1.0;
    }

    if ( std::fabs( dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors = 0.0;
    }

    // Compute and store true anomaly.
    keplerianElements[ trueAnomalyIndex ] 
        = std::acos( dotProductPositionAndEccentricityVectors );    

    // Check if the quandrant is correct.
    if ( trueAnomalyQuandrantCondition < 0.0 )
    {
        keplerianElements[ trueAnomalyIndex ]
            = 2.0 * sml::SML_PI - keplerianElements[ trueAnomalyIndex ];
    }

    return keplerianElements;
}

//! Convert true anomaly to elliptical eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical orbits (0 <= eccentricity < 1.0).
 * 
 * This implementation performs a check on the eccentricity and throws an error if it is 
 * non-elliptical, i.e., eccentricity < 0.0 and eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real         Real type
 * @param  trueAnomaly  True anomaly                 [rad]
 * @param  eccentricity Eccentricity                 [-]
 * @return              Elliptical eccentric anomaly [rad]
 */
template< typename Real >
Real convertTrueAnomalyToEllipticalEccentricAnomaly( const Real trueAnomaly, 
                                                     const Real eccentricity )
{
    if ( eccentricity >= 1.0 || eccentricity < 0.0 )
    {
        throw std::runtime_error( "ERROR: Eccentricity is non-elliptical!" );
    }

    // Compute sine and cosine of eccentric anomaly.
    const Real sineOfEccentricAnomaly 
        = std::sqrt( 1.0 - std::pow( eccentricity, 2.0 ) )
          * std::sin( trueAnomaly ) / ( 1.0 + eccentricity * std::cos( trueAnomaly ) );
    const Real cosineOfEccentricAnomaly 
        = ( eccentricity + std::cos( trueAnomaly ) ) 
          / ( 1.0 + eccentricity * std::cos( trueAnomaly ) );

    // Return elliptical eccentric anomaly.
    return std::atan2( sineOfEccentricAnomaly, cosineOfEccentricAnomaly );
}

//! Convert true anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts true anomaly to hyperbolic eccentric anomaly for hyperbolic orbits
 * (eccentricity > 1.0). 
 *
 * This implementation performs a check on the eccentricity and throws an error if it is 
 * non-hyperbolic, i.e., eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real         Real type 
 * @param  trueAnomaly  True anomaly                 [rad]
 * @param  eccentricity Eccentricity                 [-]
 * @return              Hyperbolic eccentric anomaly [rad]
 */
template< typename Real > 
Real convertTrueAnomalyToHyperbolicEccentricAnomaly( const Real trueAnomaly,
                                                     const Real eccentricity )
{
    if ( eccentricity <= 1.0 )
    {
        throw std::runtime_error( "ERROR: Eccentricity is non-hyperbolic!" );
    }

    // Compute hyperbolic sine and hyperbolic cosine of hyperbolic eccentric anomaly.
    const Real hyperbolicSineOfHyperbolicEccentricAnomaly
        = std::sqrt( std::pow( eccentricity, 2.0 ) - 1.0 )
          * std::sin( trueAnomaly ) / ( 1.0 + std::cos( trueAnomaly ) );

    const Real hyperbolicCosineOfHyperbolicEccentricAnomaly
        = ( std::cos( trueAnomaly ) + eccentricity ) / ( 1.0 + std::cos( trueAnomaly ) );

    // Return hyperbolic eccentric anomaly.
    // The inverse hyperbolic tangent is computed here manually, since the atanh() function is not
    // available in older C++ compilers.
    const Real angle
        = hyperbolicSineOfHyperbolicEccentricAnomaly 
          / hyperbolicCosineOfHyperbolicEccentricAnomaly;
    return 0.5 * ( std::log( 1.0 + angle ) - std::log( 1.0 - angle ) );
}

//! Convert true anomaly to eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical and hyperbolic orbits
 * (eccentricity < 1.0 && eccentricity > 1.0). This function is essentially a wrapper for
 * functions that treat each case. It should be used in cases where the eccentricity of the orbit
 * is not known a priori. 
 *
 * This implementation performs a check on the eccentricity and throws an error for 
 * eccentricity < 0.0 and parabolic orbits, which have not been implemented.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @sa convertTrueAnomalyToEllipticalEccentricAnomaly, 
 *     convertTrueAnomalyToHyperbolicEccentricAnomaly
 * @tparam Real         Real type 
 * @param  trueAnomaly  True anomaly       [rad]
 * @param  eccentricity Eccentricity       [-]
 * @return              Eccentric anomaly  [rad]
 */
template< typename Real >   
Real convertTrueAnomalyToEccentricAnomaly( const Real trueAnomaly, const Real eccentricity )
{
    Real eccentricAnomaly = 0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < 0.0 )
    {
        throw std::runtime_error( "ERROR: Eccentricity is negative!" );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity - 1.0 ) < std::numeric_limits< Real >::epsilon( ) )
    {
        throw std::runtime_error( "ERROR: Parabolic orbits have not been implemented!" );
    }

    // Check if orbit is elliptical and compute eccentric anomaly.
    else if ( eccentricity >= 0.0 && eccentricity < 1.0 )
    {
        eccentricAnomaly 
            = convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly, eccentricity );
    }

    // Check if orbit is hyperbolic and compute eccentric anomaly.
    else if ( eccentricity > 1.0 )
    {
        eccentricAnomaly 
            = convertTrueAnomalyToHyperbolicEccentricAnomaly( trueAnomaly, eccentricity );
    }

    return eccentricAnomaly;
}

//! Convert elliptical eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical orbits (0 <= eccentricity < 1.0).
 *
 * This implementation performs a check on the eccentricity and throws an error if it is 
 * non-elliptical, i.e., eccentricity < 0.0 and eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real                       Real type  
 * @param  ellipticalEccentricAnomaly Elliptical eccentric anomaly [rad] 
 * @param  eccentricity               Eccentricity                 [-]
 * @return                            Mean anomaly                 [rad]
 */
template< typename Real >
Real convertEllipticalEccentricAnomalyToMeanAnomaly( const Real ellipticalEccentricAnomaly,
                                                     const Real eccentricity )
{
    if ( eccentricity >= 1.0 || eccentricity < 0.0 )
    {
        throw std::runtime_error( "ERROR: Eccentricity is non-elliptical!" );
    }

    return ellipticalEccentricAnomaly - eccentricity * std::sin( ellipticalEccentricAnomaly );
}

//! Convert hyperbolic eccentric anomaly to mean anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to mean anomaly for hyperbolic orbits
 * (eccentricity > 1.0). 
 *
 * This implementation performs a check on the eccentricity and throws an error if it is 
 * non-hyperbolic, i.e., eccentricity >= 1.0.
 *
 * The equations used can be found in (Chobotov, 2002).
 *
 * @tparam Real                       Real type  
 * @param  hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly [rad]
 * @param  eccentricity               Eccentricity                 [-]
 * @return                            Mean anomaly                 [rad]
 */
template< typename Real >
Real convertHyperbolicEccentricAnomalyToMeanAnomaly(
    const Real hyperbolicEccentricAnomaly, const Real eccentricity )
{
    if ( eccentricity <= 1.0 )
    {
        throw std::runtime_error( "ERROR: Eccentricity is non-hyperbolic!" );
    }

    return eccentricity * std::sinh( hyperbolicEccentricAnomaly ) - hyperbolicEccentricAnomaly;
}

//! Convert eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical and hyperbolic orbits
 * (eccentricity < 1.0 && eccentricity > 1.0). This function is essentially a wrapper for
 * functions that treat each case. It should be used in cases where the eccentricity of the orbit
 * is not known a priori. 
 *
 * Currently, this implementation performs a check on the eccentricity and throws an error for 
 * eccentricity < 0.0 and parabolic orbits, which have not been implemented. 
 *
 * The equations used can be found in (Chobotov, 2002).
 * 
 * @sa convertEllipticalEccentricAnomalyToMeanAnomaly, 
 *     convertHyperbolicEccentricAnomalyToMeanAnomaly
 * @tparam Real             Real type  
 * @param  eccentricity     Eccentricity      [-]
 * @param  eccentricAnomaly Eccentric anomaly [rad]
 * @return                  Mean anomaly      [rad]
 */
template< typename Real >
Real convertEccentricAnomalyToMeanAnomaly( const Real eccentricAnomaly, const Real eccentricity )
{
    Real meanAnomaly = 0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < 0.0 )
    {
        throw std::runtime_error( "ERROR: Eccentricity is negative!" );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity - 1.0 ) < std::numeric_limits< Real >::epsilon( ) )
    {
        throw std::runtime_error( "ERROR: Parabolic orbits have not been implemented!" );
    }

    // Check if orbit is elliptical and compute mean anomaly.
    else if ( eccentricity >= 0.0 && eccentricity < 1.0 )
    {
        meanAnomaly 
            = convertEllipticalEccentricAnomalyToMeanAnomaly( eccentricAnomaly, eccentricity );
    }

    // Check if orbit is hyperbolic and compute mean anomaly.
    else if ( eccentricity > 1.0 )
    {
        meanAnomaly 
            = convertHyperbolicEccentricAnomalyToMeanAnomaly( eccentricAnomaly, eccentricity );
    }

    return meanAnomaly;
}

//! Compute two-body mean motion.
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
                              const Real massOfOrbitingBody = 0.0 )
{
    return std::sqrt( ( ( SAM_GRAVITATIONAL_CONSTANT * massOfOrbitingBody )
                        + gravitationalParameterOfCentralBody ) 
                      / ( semiMajorAxis * semiMajorAxis * semiMajorAxis ) );
}

} // namespace sam

#endif // SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP

/*!
 * References
 *  Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002. 
 */
