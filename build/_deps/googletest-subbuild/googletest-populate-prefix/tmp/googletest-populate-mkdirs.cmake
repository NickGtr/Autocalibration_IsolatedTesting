# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-src"
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-build"
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix"
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Data/gautt/bundle_adjustment/autocalibration_project/build/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
