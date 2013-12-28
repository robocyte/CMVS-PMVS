# CMake powered build system for CMVS
# Create by Pierre Moulon
# First draft the 4 june 2010
# Second draft 24 sept 2010
# Third draft (linux) 12 oct 2010 

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CMVS-PMVS2 C CXX)

IF(WIN32)

  SET(ADDITIONAL_WIN_INCLUDES 
        ./thirdParty/miniBoost
        ./thirdParty/jpeg)

ENDIF(WIN32)

INCLUDE_DIRECTORIES(
  ./program/base/image
  ./program/base/numeric
  ./program/base/pmvs

  ./thirdParty/
  ./thirdParty/gsl-1.13/include
  # http://eris.liralab.it/iCub/downloads/packages/windows/msvc9/
  ./thirdParty/graclus1.2/metisLib
  ./thirdParty/lmfit-3.2/lib
  ${ADDITIONAL_WIN_INCLUDES}
)


IF(WIN32)


ELSE(WIN32)

  FIND_LIBRARY(GSL_LIB NAMES gsl)
  FIND_LIBRARY(GSLCBLAS_LIB NAMES gslcblas)
  FIND_PACKAGE(JPEG REQUIRED)
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(Threads REQUIRED)
  FIND_PACKAGE(Boost REQUIRED)
  SET(PRECOMPILED_LIBS
    ${JPEG_LIBRARY}
    ${GSL_LIB}
    ${GSLCBLAS_LIB}
    ${BLAS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})

ENDIF (WIN32)

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(thirdParty)

