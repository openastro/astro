/*
 * Copyright (c) 2014-2015 Kartik Kumar, Dinamica Srl
 * Copyright (c) 2014-2015 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2015 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2015 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTRO_CONSTANTS_HPP
#define ASTRO_CONSTANTS_HPP

namespace astro
{

//! Gravitational constant [m^3 s^-2] (Standish, 1995).
const static double ASTRO_GRAVITATIONAL_CONSTANT = 6.67259e-11;

//! Julian day in seconds (NASA, 2012).
const static double ASTRO_JULIAN_DAY_IN_SECONDS = 86400.0;

//! Julian year in days (NASA, 2012).
const static double ASTRO_JULIAN_YEAR_IN_DAYS = 365.25;

//! Julian year in seconds.
const static double ASTRO_JULIAN_YEAR_IN_SECONDS = 3.15576e7;

//! Astronautical Unit in km
const static double ASTRO_AU_IN_KM = 149597870.7;

} // namespace astro

#endif // ASTRO_CONSTANTS_HPP

/*!
 * References
 *  Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards", in Highlights
 *   of Astronomy (I. Appenzeller, ed.), Table 1, Kluwer Academic Publishers, Dordrecht.
 */
