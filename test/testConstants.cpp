/*
 * Copyright (c) 2014-2015 Kartik Kumar, Dinamica Srl
 * Copyright (c) 2014-2015 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2015 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2015 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch.hpp>

#include <Astro/constants.hpp>

namespace astro
{
namespace tests
{

TEST_CASE( "Test definition of constants", "[constants]")
{
    REQUIRE( ASTRO_GRAVITATIONAL_CONSTANT == 6.67259e-11 );
    REQUIRE( ASTRO_JULIAN_DAY_IN_SECONDS  == 86400.0     );
    REQUIRE( ASTRO_JULIAN_YEAR_IN_DAYS    == 365.25      );
    REQUIRE( ASTRO_JULIAN_YEAR_IN_SECONDS == 3.15576e7   );
    REQUIRE( ASTRO_AU_IN_KM               == 149597870.7 );
}

} // namespace tests
} // namespace astro
