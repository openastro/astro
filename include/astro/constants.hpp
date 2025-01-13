/*
 * Copyright (c) 2014-2025 Kartik Kumar (me@kartikkumar.com)
 * Copyright (c) 2014-2016 Marko Jankovic, DFKI GmbH
 * Copyright (c) 2014-2016 Natalia Ortiz, University of Southampton
 * Copyright (c) 2014-2016 Juan Romero, University of Strathclyde
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

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

//! Astronautical Unit in km (NASA, 2012).
const static double ASTRO_AU_IN_KM = 149597870.7;

//! Start of Gregorian epoch in Julian days (Ramsey, 2016).
const static double ASTRO_GREGORIAN_EPOCH_IN_JULIAN_DAYS = 1721425.5;

//! Speed of light [m s^-1] (Wikipedia, 2018).
const static double ASTRO_SPEED_OF_LIGHT = 299792458.0;

} // namespace astro

/*!
 * References
 *  Standish, E.M. (1995) "Report of the IAU WGAS Sub-Group on Numerical Standards", in Highlights
 *   of Astronomy (I. Appenzeller, ed.), Table 1, Kluwer Academic Publishers, Dordrecht.
 *  NASA (2012) Jet Propulsion Laboratory, Solar System Dynamics Astrodynamic Constants
 *   accessed 15 March 2016 on: http://ssd.jpl.nasa.gov/?constants
 *  Ramsey, C.B. (2016), Oxford Radiocarbon Accelerator Unit, Research Lab for
 *	 Archaeology, Dyson Perrins Building, South Parks Road, Oxford, OX1 3QY United Kingdom.
 *   accessed 15 March 2016 on: https://c14.arch.ox.ac.uk/oxcalhelp/hlp_analysis_calend.html
 */
