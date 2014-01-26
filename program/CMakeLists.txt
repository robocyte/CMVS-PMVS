# CMake powered build system for CMVS
# Create by Pierre Moulon
# Modified by David Ullmann
# First draft the 4 june 2010
# Second draft 24 sept 2010
# Third draft (linux) 12 oct 2010 

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(CMVS-PMVS2 C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


# Options
if(MSVC)
  option(MSVC_USE_STATIC_CRT "MS Visual Studio: Link to static C-Run Time Library." ON)
endif(MSVC)
option(USE_OPENMP "Use OpenMP for parallelization" ON)


# Packages
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  message(STATUS "Found Eigen version: " ${EIGEN3_VERSION} " in " ${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

find_package(OpenMP QUIET) 
if (OPENMP_FOUND) 
  message(STATUS "Found OpenMP")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
  option(USE_OPENMP "Use OpenMP for parallelization" OFF)
endif(OPENMP_FOUND)
 
 
# Includes
if(WIN32)
  set(ADDITIONAL_WIN_INCLUDES 
	  thirdParty/miniBoost
	  thirdParty/jpeg)
endif(WIN32)

include_directories(program/base/image
					program/base/numeric
					program/base/pmvs
					thirdParty/
					thirdParty/graclus1.2/metisLib
					${ADDITIONAL_WIN_INCLUDES}
					${EIGEN3_INCLUDE_DIR})

					
# Config
if(MSVC)
  # Disable warning about the insecurity of using "std::copy"
  add_definitions("/wd4996")

  list(APPEND C_CXX_FLAGS
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(FLAG_VAR ${C_CXX_FLAGS})
    if(MSVC_USE_STATIC_CRT)
      if(${FLAG_VAR} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${FLAG_VAR} "${${FLAG_VAR}}")
      endif(${FLAG_VAR} MATCHES "/MD")
    else(MSVC_USE_STATIC_CRT)
      if(${FLAG_VAR} MATCHES "/MT")
        string(REGEX REPLACE "/MT" "/MD" ${FLAG_VAR} "${${FLAG_VAR}}")
      endif(${FLAG_VAR} MATCHES "/MT")
    endif(MSVC_USE_STATIC_CRT)
  ENDFOREACH()
endif(MSVC)


add_subdirectory(base)
add_subdirectory(main)
add_subdirectory(thirdParty)
