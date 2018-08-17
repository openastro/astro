/*
 * Copyright (c) 2014-2018 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch2/catch.hpp>

#include "astro/stateVectorIndices.hpp"

namespace astro
{
namespace tests
{

TEST_CASE( "Test definition of Cartesian element indices", "[constants]")
{
    REQUIRE( xPositionIndex == 0 );
    REQUIRE( yPositionIndex == 1 );
    REQUIRE( zPositionIndex == 2 );
    REQUIRE( xVelocityIndex == 3 );
    REQUIRE( yVelocityIndex == 4 );
    REQUIRE( zVelocityIndex == 5 );
}

TEST_CASE( "Test definition of Keplerian element indices", "[constants]")
{
    REQUIRE( semiMajorAxisIndex            == 0 );
    REQUIRE( semiLatusRectumIndex          == 0 );
    REQUIRE( eccentricityIndex             == 1 );
    REQUIRE( inclinationIndex              == 2 );
    REQUIRE( argumentOfPeriapsisIndex      == 3 );
    REQUIRE( longitudeOfAscendingNodeIndex == 4 );
    REQUIRE( trueAnomalyIndex              == 5 );
    REQUIRE( meanAnomalyIndex              == 5 );
}

} // namespace tests
} // namespace astro
