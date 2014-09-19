/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#ifndef SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP
#define SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP

namespace sam
{

template< typename Real, typename Vector6, typename Vector3 >
Vector6 convertCartesianToKeplerianElements(
    const Vector6& cartesianElements, const Real centralBodyGravitationalParameter,
    const Real tolerance = 10.0 * std::numeric_limits< Real >::epsilon( ) )
{
    
}

// convert true anomaly to eccentric anomaly

// convert eccentric anomaly to mean anomaly

} // namespace sam

#endif // SAM_ORBITAL_ELEMENT_CONVERSIONS_HPP
