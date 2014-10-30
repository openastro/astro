SAM
===

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT) [![Build Status](https://travis-ci.org/kartikkumar/sam.svg?branch=master)](https://travis-ci.org/kartikkumar/sam)

`Simple Astrodynamics Library (SAM)` is a C++ template library that provides some basic astrodynamics functionality.

Features
------

  - Header-only
  - Orbital element conversions
  - Methods related to the 2-Body Problem
  - Useful physical constants 
  - Full suite of tests

Requirements
------

This project requires the following:

  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/), [MinGW](http://www.mingw.org/)
  - [CMake](http://www.cmake.org)
 
Installation
------

Run the following commands to download, build, and install this project.

    git clone https://www.github.com/kartikkumar/sam
    cd sam
    git submodule init && git submodule update
    mkdir build
    cd build
    cmake ..
    cmake --build .

To install the header files and libraries, run the following from within the `build` directory:

    make install

Note that dependencies are installed by fetching them online, in case they cannot be detected on your local system. If the build process fails, check the error log given. Typically, building fails due to timeout. Simply run the `cmake --build .` command once more.

Build options
-------------

You can pass the follow command-line options when running `CMake`:

  - `-DBUILD_DOCS=on`: build the [Doxygen](http://www.doxygen.org "Doxygen homepage") documentation ([LaTeX](http://www.latex-project.org/) must be installed)
  - `-DBUILD_TESTS=on`: build tests (execute tests from build-directory using `make test`)
  - `-DBUILD_WITH_EIGEN=on`: build tests using [Eigen](http://eigen.tuxfamily.org/)  
  - `-DBUILD_SHARED_LIBS=on`: build shared libraries instead of static
  - `-DCMAKE_INSTALL_PREFIX`: set path prefix for install script (`make install`); if not set, defaults to usual locations.

Contributing
------------

Once you've made your great commits:

1. [Fork](https://github.com/kartikkumar/sam/fork) SAM
2. Create a topic branch - `git checkout -b my_branch`
3. Push to your branch - `git push origin my_branch`
4. Create a [Pull Request](http://help.github.com/pull-requests/) from your branch
5. That's it!

License
------

See `LICENSE.md`.

Disclaimer
------

The copyright holders are not liable for any damage(s) incurred due to improper use of `SAM`.

Contact
------

Shoot an [email](mailto:me@kartikkumar.com?subject=SAM) if you have any questions.

TODO
------

 - Add other basic orbital element conversions
 - Extend test suite
 - Figure out how to build tests with standard library and Eigen in the same build tree
