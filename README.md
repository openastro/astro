Astro
===

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT) [![Build Status](https://travis-ci.org/kartikkumar/astro.svg?branch=master)](https://travis-ci.org/kartikkumar/astro) [![Coverity Scan Build Status](https://scan.coverity.com/projects/3697/badge.svg)](https://scan.coverity.com/projects/3697) [![Coverage Status](https://coveralls.io/repos/kartikkumar/astro/badge.png)](https://coveralls.io/r/kartikkumar/astro)

Astro is a C++ template library that provides some basic astrodynamics functionality. It's intended to be lightweight and simple to use in other projects. A CMake module is available to make it easy to include SML in other CMake-based projects: [FindAstro.cmake](https://github.com/kartikkumar/cmake-modules/Modules/FindAstro.cmake).

Features
------

  - Header-only
  - Orbital element conversions
  - Methods related to the 2-Body Problem
  - Useful physical constants
  - Full suite of tests

Requirements
------

To install this project, please ensure that you have installed the following (install guides are provided on the respective websites):

  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/), [MinGW](http://www.mingw.org/)
  - [CMake](http://www.cmake.org)
  - [Doxygen](http://www.doxygen.org "Doxygen homepage") (optional)
  - [Gcov](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html) (optional)
  - [LCOV](http://ltp.sourceforge.net/coverage/lcov.php) (optional)

In addition, Astro depends on the following libraries:

  - [SML](https://www.github.com/kartikkumar/sml) (math library)
  - [CATCH](https://www.github.com/philsquared/Catch) (unit testing library necessary for `BUILD_TESTS` option)
  - [Eigen](http://eigen.tuxfamily.org/) (linear algebra library necessary for `BUILD_TESTS_WITH_EIGEN` option)

These dependencies will be downloaded and configured automagically if not already present locally (requires an internet connection).

Installation
------

Run the following commands to download, build, and install this project.

    git clone https://www.github.com/kartikkumar/sml
    cd sml
    git submodule init && git submodule update
    mkdir build && cd build
    cmake .. && cmake --build .

To install the header files, run the following from within the `build` directory:

    make install

Note that dependencies are installed by fetching them online, in case they cannot be detected on your local system. If the build process fails, check the error log given. Typically, building fails due to timeout. Simply run the `cmake --build .` command once more.

Build options
-------------

You can pass the following, general command-line options when running CMake:

  - `-DCMAKE_INSTALL_PREFIX[=$install_dir]`: set path prefix for install script (`make install`); if not set, defaults to usual locations
  - `-DBUILD_DOCS[=ON|OFF (default)]`: build the [Doxygen](http://www.doxygen.org "Doxygen homepage") documentation ([LaTeX](http://www.latex-project.org/) must be installed with `amsmath` package)
  - `-DBUILD_TESTS[=ON|OFF (default)]`: build tests (execute tests from build-directory using `ctest -V`)
  - `-DBUILD_DEPENDENCIES[=ON|OFF (default)]`: force local build of dependencies, instead of first searching system-wide using `find_package()`

The following commands are conditional and can only be set if `BUILD_TESTS = ON`:

 - `-DBUILD_TESTS_WITH_EIGEN[=ON|OFF (default)]`: build tests using [Eigen](http://eigen.tuxfamily.org/) (execute tests from build-directory using `ctest -V`)
 - `-DBUILD_COVERAGE[=ON|OFF (default)]`: build code coverage using [Gcov](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html) and [LCOV](http://ltp.sourceforge.net/coverage/lcov.php) (both must be installed; requires [GCC](https://gcc.gnu.org/) compiler; execute coverage analysis from build-directory using `make coverage`)

Pass these options either directly to the `cmake ..` build command or run `ccmake ..` instead to bring up the interface that can be used to toggle options.

Contributing
------------

Once you've made your great commits:

1. [Fork](https://github.com/kartikkumar/astro/fork) Astro
2. Create a topic branch - `git checkout -b my_branch`
3. Push to your branch - `git push origin my_branch`
4. Create a [Pull Request](http://help.github.com/pull-requests/) from your branch
5. That's it!

Disclaimer
------

The copyright holders are not liable for any damage(s) incurred due to improper use of Astro.

TODO
------

  - Add other basic orbital element conversions
  - Extend test suite
  - Figure out better way (avoiding code duplication) to build `STL`-based and `Eigen`-based tests in the same build tree
  - Add version detection in CMake module so that find_package respects minimum version required.
  - Find a way to provide an option to clean installation.
