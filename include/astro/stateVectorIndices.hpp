/*
 * Copyright (c) 2014-2022 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

namespace astro
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

//! Keplerian element array indices.
enum KeplerianElementIndices
{
    semiMajorAxisIndex            = 0,
    semiLatusRectumIndex          = 0,
    eccentricityIndex             = 1,
    inclinationIndex              = 2,
    argumentOfPeriapsisIndex      = 3,
    longitudeOfAscendingNodeIndex = 4,
    trueAnomalyIndex              = 5,
    meanAnomalyIndex              = 5
};

} // namespace astro
