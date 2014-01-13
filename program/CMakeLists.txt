# CMake powered build system for CMVS
# Create by Pierre Moulon
# First draft the 4 june 2010
# Second draft 24 sept 2010
# Third draft (linux) 12 oct 2010 

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CMVS-PMVS2 C CXX)

# Options
# Multithreading using OpenMP
OPTION(OPENMP "Enable OpenMP)" ON)

IF (MSVC)
  OPTION(MSVC_USE_STATIC_CRT
    "MS Visual Studio: Link to static C-Run Time Library." ON)
ENDIF (MSVC)

IF (WIN32)
  SET(ADDITIONAL_WIN_INCLUDES 
        ./thirdParty/miniBoost
        ./thirdParty/jpeg)
ENDIF (WIN32)

INCLUDE_DIRECTORIES(
  ./program/base/image
  ./program/base/numeric
  ./program/base/pmvs

  ./thirdParty/
  ./thirdParty/graclus1.2/metisLib
  ${ADDITIONAL_WIN_INCLUDES}
  ${EIGEN_INCLUDE_DIRS}
)

IF (WIN32)

ELSE (WIN32)
  FIND_PACKAGE(Eigen REQUIRED)
  FIND_PACKAGE(JPEG REQUIRED)
  FIND_PACKAGE(Boost REQUIRED)
  SET(PRECOMPILED_LIBS
    ${JPEG_LIBRARY})
ENDIF (WIN32)

IF (OPENMP)
    FIND_PACKAGE(OpenMP QUIET)
    IF (OPENMP_FOUND)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ENDIF (OPENMP_FOUND)
ENDIF (OPENMP)

IF (MSVC)
  # Disable warning about the insecurity of using "std::copy"
  ADD_DEFINITIONS("/wd4996")

  LIST(APPEND C_CXX_FLAGS
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  FOREACH(FLAG_VAR ${C_CXX_FLAGS})
    IF (MSVC_USE_STATIC_CRT)
      IF (${FLAG_VAR} MATCHES "/MD")
        STRING(REGEX REPLACE "/MD" "/MT" ${FLAG_VAR} "${${FLAG_VAR}}")
      ENDIF (${FLAG_VAR} MATCHES "/MD")
    ELSE (MSVC_USE_STATIC_CRT)
      IF (${FLAG_VAR} MATCHES "/MT")
        STRING(REGEX REPLACE "/MT" "/MD" ${FLAG_VAR} "${${FLAG_VAR}}")
      ENDIF (${FLAG_VAR} MATCHES "/MT")
    ENDIF (MSVC_USE_STATIC_CRT)
  ENDFOREACH()
ENDIF (MSVC)

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(thirdParty)

