/*    
 * Copyright (c) 2014, K. Kumar (me@kartikkumar.com)
 * All rights reserved.
 */

#include <catch.hpp>

#include <SAM/constants.hpp>

namespace sam
{
namespace unit_tests
{

TEST_CASE( "Test definition of constants", "[constants]")
{
    REQUIRE( SAM_GRAVITATIONAL_CONSTANT == 6.67259e-11 );
    REQUIRE( SAM_JULIAN_DAY_IN_SECONDS == 86400.0 );
    REQUIRE( SAM_JULIAN_YEAR_IN_DAYS == 365.25 );
    REQUIRE( SAM_JULIAN_YEAR_IN_SECONDS == 3.15576e7 );
}

} // namespace unit_tests
} // namespace sam
