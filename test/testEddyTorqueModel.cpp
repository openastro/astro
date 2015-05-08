/*
 * Copyright (c) 
 */
#include <cmath>

#include <vector>

#include <catch.hpp>

#include <Astro/eddyCurrentModel.hpp>

using namespace std;

namespace astro
{

namespace tests
{

typedef double Real;
typedef std::vector< Real > Vector3;

TEST_CASE( "Obtain eddy current torque: test 1", "[obtain-eddy-torque-1]" )
{
    // Set expected eddy current torque vector [N * m].
    Vector3 expectedEddyTorque( 3 );
    expectedEddyTorque[ 0 ] =  0.095000000000000;
    expectedEddyTorque[ 1 ] =  0.065000000000000;
    expectedEddyTorque[ 2 ] = -0.149000000000000;

    // Set magnetic moment vector [A * m^2].
    Vector3 magneticMoment( 3 );
    magneticMoment[ 0 ] = 100.0;
    magneticMoment[ 1 ] = 1000.0;
    magneticMoment[ 2 ] = 500.0;

    // Set magnetic field vector [T].
    Vector3 magneticField( 3 );
    magneticField[ 0 ] = 150e-6;
    magneticField[ 1 ] = 10e-6;
    magneticField[ 2 ] = 100e-6;

    // Set epsilon = error between expected value and computed value.
    const Real epsilon = 1.0e-10;

     //! Compute eddy current torque.
    const Vector3 eddyTorque = computeEddyTorque( magneticMoment, 
                                                  magneticField ); 

    // Check if computed torque matches expected value.
    REQUIRE( std::fabs(eddyTorque[ 0 ] - expectedEddyTorque[ 0 ]) <= epsilon );
    REQUIRE( std::fabs(eddyTorque[ 1 ] - expectedEddyTorque[ 1 ]) <= epsilon );
    REQUIRE( std::fabs(eddyTorque[ 2 ] - expectedEddyTorque[ 2 ]) <= epsilon );
}

TEST_CASE( "Obtain eddy current torque: test 2", "[obtain-eddy-torque-2]" )
{
    // Set expected eddy current torque vector [N * m].
    Vector3 expectedEddyTorque( 3 );
    expectedEddyTorque[ 0 ] = 0.0;
    expectedEddyTorque[ 1 ] = 0.0;
    expectedEddyTorque[ 2 ] = 0.0;

    // Set magnetic moment vector [A * m^2].
    Vector3 magneticMoment( 3 );
    magneticMoment[ 0 ] = 0.0;
    magneticMoment[ 1 ] = 0.0;
    magneticMoment[ 2 ] = 1150.0;

    // Set magnetic field vector [T].
    Vector3 magneticField( 3 );
    magneticField[ 0 ] = 0.0;
    magneticField[ 1 ] = 0.0;
    magneticField[ 2 ] = 127e-6;

    // Set epsilon = error between expected value and computed value.
    const Real epsilon = 1.0e-10;

     //! Compute eddy current torque.
    const Vector3 eddyTorque = computeEddyTorque( magneticMoment, 
                                                  magneticField ); 

    // Check if computed torque matches expected value.
    REQUIRE( std::fabs(eddyTorque[ 0 ] - expectedEddyTorque[ 0 ]) <= epsilon );
    REQUIRE( std::fabs(eddyTorque[ 1 ] - expectedEddyTorque[ 1 ]) <= epsilon );
    REQUIRE( std::fabs(eddyTorque[ 2 ] - expectedEddyTorque[ 2 ]) <= epsilon );
}


} // namespace tests
} // namespace astro
