/*
 * Copyright (c) 2016 Kartik Kumar, Dinamica Srl (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_RELATIVE_MOTION_HPP
#define ASTRO_RELATIVE_MOTION_HPP

#include <cmath>

#include "Astro/stateVectorIndices.hpp"

namespace astro
{

//! Propagate constant-force Clohessy-Wiltshire solution.
/*!
 * Propagates the constant-force Clohessy-Wiltshire solution to the Hill equations. The
 * Clohessy-Wiltshire solution is a closed-form solution to the equations of motion describing the
 * relative motion between a target and a chaser. The solution is predicated on the assumption that
 * the target is on a circular orbit. The solution is expressed in terms of the Hill frame. The
 * reader is referred to Fehse (2003) for a concise statement of the problem and an in-depth
 * explanation of the derivation of this solution.
 *
 * The Clohessy-Wiltshire solution can be evaluated at a specified final time (\f$t_{2}\f$), given
 * an initial state (\f$\bar{x}_{0}\f$) in the Hill frame, to obtain the chaser's final state. A
 * constant thrust acceleration can be specified, i.e., a constant acceleration is added to the
 * right-hand side of the equations of motion.
 *
 * It should be noted that the definition of the Hill frame conforms to the norm outlined by Fehse
 * (2003). The user should be careful to adhere to this when executing this function.
 *
 * @sa computeUnperturbedHillEquations
 * @tparam  Real                   Real type
 * @tparam  Vector6                6-Vector type
 * @tparam  Vector3                3-Vector type
 * @param   initialState           Cartesian state at initial epoch in Hill frame
 * @param   finalTime              Final epoch
 * @param   targetMeanMotion       Mean motion of target orbit (circular)
 * @param   thrustAcceleration     Thrust acceleration to be applied to chaser, in Hill frame
 * @return                         Final state in Hill frame
 */
template< typename Real, typename Vector6, typename Vector3 >
Vector6 propagateClohessyWiltshireSolution( const Vector6& initialState,
                                            const Real finalTime,
                                            const Real targetMeanMotion,
                                            const Vector3& thrustAcceleration )
{
    Vector6 finalState = initialState;

    finalState[ astro::xPositionIndex ]
        = // unperturbed terms
          ( 4.0 / targetMeanMotion * initialState[ astro::xVelocityIndex ]
            - 6.0 * initialState[ astro::zPositionIndex ] )
            * std::sin( targetMeanMotion * finalTime )
          - 2.0 * initialState[ astro::zVelocityIndex ] / targetMeanMotion
            * std::cos( targetMeanMotion * finalTime )
          + ( 6.0 * targetMeanMotion * initialState[ astro::zPositionIndex ]
              - 3.0 * initialState[ astro::xVelocityIndex ] ) * finalTime
          + ( initialState[ astro::xPositionIndex ]
              + 2.0 * initialState[ astro::zVelocityIndex ] / targetMeanMotion )
          // constant force terms
          + 2.0 / ( targetMeanMotion * targetMeanMotion )
            * thrustAcceleration[ astro::zPositionIndex ]
            * ( targetMeanMotion * finalTime - std::sin( targetMeanMotion * finalTime ) )
          + thrustAcceleration[ astro::xPositionIndex ]
            * ( 4.0 / ( targetMeanMotion * targetMeanMotion )
                * ( 1.0 - std::cos( targetMeanMotion * finalTime ) )
                - 1.5 * finalTime * finalTime );

    finalState[ astro::yPositionIndex ]
        = // unperturbed terms
          initialState[ astro::yPositionIndex ] * std::cos( targetMeanMotion * finalTime )
          + initialState[ astro::yVelocityIndex ] / targetMeanMotion
            * std::sin( targetMeanMotion * finalTime )
          // constant force terms
          + thrustAcceleration[ astro::yPositionIndex ] / ( targetMeanMotion * targetMeanMotion )
            * ( 1.0 - std::cos( targetMeanMotion * finalTime ) );

    finalState[ astro::zPositionIndex ]
        = // unperturbed terms
          ( 2.0 * initialState[ astro::xVelocityIndex ] / targetMeanMotion
            - 3.0 * initialState[ astro::zPositionIndex ] )
            * std::cos( targetMeanMotion * finalTime )
          + initialState[ astro::zVelocityIndex ] / targetMeanMotion
          * std::sin( targetMeanMotion * finalTime )
          + ( 4.0 * initialState[ astro::zPositionIndex ]
            - 2.0 * initialState[ astro::xVelocityIndex ] / targetMeanMotion )
          // constant force terms
          + 2.0 / ( targetMeanMotion * targetMeanMotion )
            * thrustAcceleration[ astro::xPositionIndex ]
            * ( std::sin( targetMeanMotion * finalTime ) - targetMeanMotion * finalTime )
          + thrustAcceleration[ astro::zPositionIndex ] / ( targetMeanMotion * targetMeanMotion )
            * ( 1.0 - std::cos( targetMeanMotion * finalTime ) );

    finalState[ astro::xVelocityIndex ]
        = // unperturbed terms
          targetMeanMotion * ( 4.0 / targetMeanMotion * initialState[ astro::xVelocityIndex ]
            - 6.0 * initialState[ astro::zPositionIndex ] )
            * std::cos( targetMeanMotion * finalTime )
          + 2.0 * initialState[ astro::zVelocityIndex ]
            * std::sin( targetMeanMotion * finalTime )
          + ( 6.0 * targetMeanMotion * initialState[ astro::zPositionIndex]
              - 3.0 * initialState[ astro::xVelocityIndex ] )
          // constant force terms
          + thrustAcceleration[ astro::zPositionIndex ]
            * 2.0 / ( targetMeanMotion * targetMeanMotion )
            * ( targetMeanMotion - targetMeanMotion * std::cos( targetMeanMotion * finalTime ) )
          + thrustAcceleration[ astro::xPositionIndex ]
            * ( 4.0 / targetMeanMotion * std::sin( targetMeanMotion * finalTime )
                - 3.0 * finalTime );

    finalState[ astro::yVelocityIndex ]
        = // unperturbed terms
          -targetMeanMotion * initialState[ astro::yPositionIndex ]
            * std::sin( targetMeanMotion * finalTime )
          + initialState[ astro::yVelocityIndex ] * std::cos( targetMeanMotion * finalTime )
          // constant force terms
          + thrustAcceleration[ astro::yPositionIndex ] / targetMeanMotion
            * std::sin( targetMeanMotion * finalTime );

    finalState[ astro::zVelocityIndex ]
        = // unperturbed terms
          -targetMeanMotion * ( 2.0 * initialState[ astro::xVelocityIndex ] / targetMeanMotion
            - 3.0 * initialState[ astro::zPositionIndex ] )
            * std::sin( targetMeanMotion * finalTime )
          + initialState[ astro::zVelocityIndex ] * std::cos( targetMeanMotion * finalTime )
          // constant force terms
          + thrustAcceleration[ astro::xPositionIndex ]
            * 2.0 / ( targetMeanMotion * targetMeanMotion )
            * ( targetMeanMotion * std::cos( targetMeanMotion * finalTime ) - targetMeanMotion )
          + thrustAcceleration[ astro::zPositionIndex ] / targetMeanMotion
            * std::sin( targetMeanMotion * finalTime );

    return finalState;
}

} // namespace astro

#endif // ASTRO_RELATIVE_MOTION_HPP

/*!
 * References
 *  Fehse, W. (2003) Automated Rendezvous and Docking of Spacecraft, Cambridge Aerospace Series 16,
 *    Cambridge University Press, Cambridge, United Kingdom.
 */
