/*
 * Copyright (c) 2014-2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2016 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2016 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_ORBITAL_ELEMENT_CONVERSIONS_HPP
#define ASTRO_ORBITAL_ELEMENT_CONVERSIONS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <SML/sml.hpp>

#include "Astro/stateVectorIndices.hpp"

namespace astro
{

//! Convert Cartesian elements to Keplerian elements.
/*!
 * Converts a given set of Cartesian elements (position, velocity) to classical (osculating)
 * Keplerian elements. See Chobotov (2006) for a full derivation of the conversion.
 *
 * The tolerance is given a default value. It should not be changed unless required for specific
 * scenarios. Below this tolerance value for eccentricity and inclination, the orbit is considered
 * to be a limit case. Essentially, special solutions are then used for parabolic, circular
 * inclined, non-circular equatorial, and circular equatorial orbits. These special solutions are
 * required because of singularities in the classical Keplerian elements. If high precision is
 * required near these singularities, users are encouraged to consider using other elements, such
 * as Modified Equinoctial Elements (MEE). It should be noted that MEE also suffer from
 * singularities, but not for zero eccentricity and inclination.
 *
 * WARNING: If eccentricity is 1.0 within tolerance, keplerianElements( 0 ) = semi-latus rectum,
 *          since the orbit is parabolic.
 * WARNING: If eccentricity is 0.0 within tolerance, argument of periapsis is set to 0.0, since the
 *          orbit is circular.
 * WARNING: If inclination is 0.0 within tolerance, longitude of ascending node is set to 0.0, since
 *          since the orbit is equatorial.
 *
 * @tparam  Real                   Real type
 * @tparam  Vector                 Vector type
 * @param   cartesianElements      Vector containing Cartesian elements                     <br>
 *                                 N.B.: Order of elements is important!                    <br>
 *                                 cartesianElements( 0 ) = x-position coordinate       [m] <br>
 *                                 cartesianElements( 1 ) = y-position coordinate       [m] <br>
 *                                 cartesianElements( 2 ) = z-position coordinate       [m] <br>
 *                                 cartesianElements( 3 ) = x-velocity coordinate     [m/s] <br>
 *                                 cartesianElements( 4 ) = y-velocity coordinate     [m/s] <br>
 *                                 cartesianElements( 5 ) = z-velocity coordinate     [m/s]
 * @param   gravitationalParameter Gravitational parameter of central body       [m^3 s^-2]
 * @param   tolerance              Tolerance used to check for limit cases
 *                                 (eccentricity, inclination)
 * @return  Converted vector of Keplerian elements                                          <br>
 *          N.B.: Order of elements is important!                                           <br>
 *          keplerianElements( 0 ) = semiMajorAxis                                      [m] <br>
 *          keplerianElements( 1 ) = eccentricity                                       [-] <br>
 *          keplerianElements( 2 ) = inclination                                      [rad] <br>
 *          keplerianElements( 3 ) = argument of periapsis                            [rad] <br>
 *          keplerianElements( 4 ) = longitude of ascending node                      [rad] <br>
 *          keplerianElements( 5 ) = true anomaly                                     [rad]
 */
template< typename Real, typename Vector6 >
Vector6 convertCartesianToKeplerianElements(
    const Vector6& cartesianElements,
    const Real gravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( ) )
{
    // Check that the Cartesian elements vector contains exactly 6 elemenets and otherwise throw
    // and error.
    if ( cartesianElements.size( ) != 6 )
    {
        throw std::runtime_error(
            "ERROR: Cartesian elements vector has more or less than 6 elements!" );
    }

    Vector6 keplerianElements = cartesianElements;

    // Set position and velocity vectors.
    std::vector< Real > position( 3 );
    position[ xPositionIndex ] = cartesianElements[ xPositionIndex ];
    position[ yPositionIndex ] = cartesianElements[ yPositionIndex ];
    position[ zPositionIndex ] = cartesianElements[ zPositionIndex ];

    std::vector< Real > velocity( 3 );
    velocity[ 0 ] = cartesianElements[ xVelocityIndex ];
    velocity[ 1 ] = cartesianElements[ yVelocityIndex ];
    velocity[ 2 ] = cartesianElements[ zVelocityIndex ];

    // Compute orbital angular momentum vector.
    const std::vector< Real > angularMomentum( sml::cross( position, velocity ) );

    // Compute semi-latus rectum.
    const Real semiLatusRectum
        = sml::squaredNorm< Real >( angularMomentum ) / gravitationalParameter;

    // Compute unit vector to ascending node.
    std::vector< Real > ascendingNodeUnitVector
        = sml::normalize< Real >(
            sml::cross( sml::getZUnitVector< std::vector< Real > >( ),
                        sml::normalize< Real >( angularMomentum ) ) );

    // Compute eccentricity vector.
    std::vector< Real > eccentricityVector
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
        = std::acos( angularMomentum[ zPositionIndex ] / sml::norm< Real >( angularMomentum ) );

    // Compute and store longitude of ascending node.
    // Define the quadrant condition for the argument of perigee.
    Real argumentOfPeriapsisQuandrantCondition = eccentricityVector[ zPositionIndex ];

    // Check if the orbit is equatorial. If it is, set the vector to the line of nodes to the
    // x-axis.
    if ( std::fabs( keplerianElements[ inclinationIndex ] ) < tolerance )
    {
        ascendingNodeUnitVector = sml::getXUnitVector< std::vector< Real > >( );

        // If the orbit is equatorial, eccentricityVector_z is zero, therefore the quadrant
        // condition is taken to be the y-component, eccentricityVector_y.
        argumentOfPeriapsisQuandrantCondition = eccentricityVector[ yPositionIndex ];
    }

    // Compute and store the resulting longitude of ascending node.
    keplerianElements[ longitudeOfAscendingNodeIndex ]
        = std::acos( ascendingNodeUnitVector[ xPositionIndex ] );

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
        if ( ascendingNodeUnitVector == sml::getXUnitVector< std::vector< Real > >( ) )
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
            = std::acos( sml::dot< Real >( sml::normalize< Real >( eccentricityVector ),
                         ascendingNodeUnitVector ) );

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
                            sml::normalize< Real >( eccentricityVector ) );

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
    keplerianElements[ trueAnomalyIndex ] = std::acos( dotProductPositionAndEccentricityVectors );

    // Check if the quandrant is correct.
    if ( trueAnomalyQuandrantCondition < 0.0 )
    {
        keplerianElements[ trueAnomalyIndex ]
            = 2.0 * sml::SML_PI - keplerianElements[ trueAnomalyIndex ];
    }

    return keplerianElements;
}

//! Convert Keplerian elements to Cartesian elements.
/*!
 * Converts a given set of Keplerian (osculating) elements to Cartesian elements (position,
 * velocity). See Chobotov (2006) for a full derivation of the conversion.
 *
 * WARNING: If eccentricity is 1.0 within tolerance, the user should provide
 *          keplerianElements( 0 ) = semi-latus rectum, since the orbit is parabolic.
 *
 * @tparam  Real                   Real type
 * @tparam  Vector                 Vector type
 * @param   keplerianElements      Vector containing Keplerian elemenets                    <br>
 *                                 N.B.: Order of elements is important!                    <br>
 *                                 keplerianElements( 0 ) = semiMajorAxis               [m] <br>
 *                                 keplerianElements( 1 ) = eccentricity                [-] <br>
 *                                 keplerianElements( 2 ) = inclination               [rad] <br>
 *                                 keplerianElements( 3 ) = argument of periapsis     [rad] <br>
 *                                 keplerianElements( 4 ) = longitude of              [rad]
 *                                                          ascending node            [rad] <br>
 *                                 keplerianElements( 5 ) = true anomaly              [rad]
 * @param   gravitationalParameter Gravitational parameter of central body       [m^3 s^-2]
 * @param   tolerance              Tolerance used to check for limit case of eccentricity
 * @return Converted vector of Cartesian elements                                           <br>
 *                                 N.B.: Order of elements is important!                    <br>
 *                                 cartesianElements( 0 ) = x-position coordinate       [m] <br>
 *                                 cartesianElements( 1 ) = y-position coordinate       [m] <br>
 *                                 cartesianElements( 2 ) = z-position coordinate       [m] <br>
 *                                 cartesianElements( 3 ) = x-velocity coordinate     [m/s] <br>
 *                                 cartesianElements( 4 ) = y-velocity coordinate     [m/s] <br>
 *                                 cartesianElements( 5 ) = z-velocity coordinate     [m/s]
 */
template< typename Real, typename Vector6 >
Vector6 convertKeplerianToCartesianElements(
    const Vector6& keplerianElements, const Real gravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( ) )
{
    Vector6 cartesianElements = keplerianElements;

    const Real semiMajorAxis                    = keplerianElements[ semiMajorAxisIndex ];
    const Real eccentricity                     = keplerianElements[ eccentricityIndex ];
    const Real inclination                      = keplerianElements[ inclinationIndex ];
    const Real argumentOfPeriapsis              = keplerianElements[ argumentOfPeriapsisIndex ];
    const Real longitudeOfAscendingNode         = keplerianElements[ longitudeOfAscendingNodeIndex ];
    const Real trueAnomaly                      = keplerianElements[ trueAnomalyIndex ];

    // Pre-compute sines and cosines of angles for efficient computation.
    const Real cosineOfInclination              = std::cos( inclination );
    const Real sineOfInclination                = std::sin( inclination );
    const Real cosineOfArgumentOfPeriapsis      = std::cos( argumentOfPeriapsis );
    const Real sineOfArgumentOfPeriapsis        = std::sin( argumentOfPeriapsis );
    const Real cosineOfLongitudeOfAscendingNode = std::cos( longitudeOfAscendingNode );
    const Real sineOfLongitudeOfAscendingNode   = std::sin( longitudeOfAscendingNode );
    const Real cosineOfTrueAnomaly              = std::cos( trueAnomaly );
    const Real sineOfTrueAnomaly                = std::sin( trueAnomaly );

    // Compute semi-latus rectum in the case the orbit is not a parabola.
    Real semiLatusRectum = 0.0;
    if ( std::fabs( eccentricity - 1.0 ) > tolerance  )
    {
        semiLatusRectum = semiMajorAxis * ( 1.0 - eccentricity * eccentricity );
    }

    // Else set the semi-latus rectum as the first element in the vector of Keplerian elements.
    else
    {
        semiLatusRectum = keplerianElements[ 0 ];
    }

    // Compute the magnitude of the orbital radius, measured from the focal point.
    const Real radiusMagnitude = semiLatusRectum / ( 1.0 + eccentricity * cosineOfTrueAnomaly );

    // Define position and velocity in the perifocal coordinate system.
    const Real xPositionPerifocal = radiusMagnitude * cosineOfTrueAnomaly;
    const Real yPositionPerifocal = radiusMagnitude * sineOfTrueAnomaly;
    const Real xVelocityPerifocal
        = -std::sqrt( gravitationalParameter / semiLatusRectum ) * sineOfTrueAnomaly;
    const Real yVelocityPerifocal
        = std::sqrt( gravitationalParameter / semiLatusRectum )
          * ( eccentricity + cosineOfTrueAnomaly );

    // Compute scalar components of rotation matrix to rotate from periforcal to Earth-Centered
    // Inertial (ECI) frame.
    const Real rotationMatrixComponent11
        = ( cosineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis
            - sineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis * cosineOfInclination );
    const Real rotationMatrixComponent12
        = ( -cosineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis
              -sineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis * cosineOfInclination );

    const Real rotationMatrixComponent21
        = ( sineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis
            + cosineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis * cosineOfInclination );
    const Real rotationMatrixComponent22
        = ( -sineOfLongitudeOfAscendingNode * sineOfArgumentOfPeriapsis
              + cosineOfLongitudeOfAscendingNode * cosineOfArgumentOfPeriapsis * cosineOfInclination );

    const Real rotationMatrixComponent31 = ( sineOfArgumentOfPeriapsis * sineOfInclination );
    const Real rotationMatrixComponent32 = ( cosineOfArgumentOfPeriapsis * sineOfInclination );

    // Compute Cartesian position and velocities.
    cartesianElements[ xPositionIndex ] = rotationMatrixComponent11 * xPositionPerifocal
                                          + rotationMatrixComponent12 * yPositionPerifocal;

    cartesianElements[ yPositionIndex ] = rotationMatrixComponent21 * xPositionPerifocal
                                          + rotationMatrixComponent22 * yPositionPerifocal;

    cartesianElements[ zPositionIndex ] = rotationMatrixComponent31 * xPositionPerifocal
                                          + rotationMatrixComponent32 * yPositionPerifocal;

    cartesianElements[ xVelocityIndex ] = rotationMatrixComponent11 * xVelocityPerifocal
                                          + rotationMatrixComponent12 * yVelocityPerifocal;

    cartesianElements[ yVelocityIndex ] = rotationMatrixComponent21 * xVelocityPerifocal
                                          + rotationMatrixComponent22 * yVelocityPerifocal;

    cartesianElements[ zVelocityIndex ] = rotationMatrixComponent31 * xVelocityPerifocal
                                          + rotationMatrixComponent32 * yVelocityPerifocal;

    return cartesianElements;
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

} // namespace astro

#endif // ASTRO_ORBITAL_ELEMENT_CONVERSIONS_HPP

/*!
 * References
 *  Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 */
