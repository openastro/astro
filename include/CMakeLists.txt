# Copyright (c) 2014-2025 Kartik Kumar (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE or copy at http://opensource.org/licenses/MIT

# The CMake setup for this project is based off of the following source:
# - https://cliutils.gitlab.io/modern-cmake

# Add interface library since this is a header-only library
add_library(astro_lib INTERFACE)
target_include_directories(astro_lib INTERFACE .)
