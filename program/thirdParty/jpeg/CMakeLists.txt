PROJECT(jpeg)

SET(LIBJPEG_SYSDEP jmemnobs.c)

SET(LIBJPEG_COMMON jaricom.c jcomapi.c jutils.c jerror.c jmemmgr.c)

SET(LIBJPEG_CLIB
  jcapimin.c jcapistd.c jcarith.c  jctrans.c  jcparam.c 
  jdatadst.c jcinit.c   jcmaster.c jcmarker.c jcmainct.c 
  jcprepct.c jccoefct.c jccolor.c  jcsample.c jchuff.c 
  jcdctmgr.c jfdctfst.c jfdctflt.c jfdctint.c)

SET(LIBJPEG_DLIB
  jdapimin.c jdapistd.c jdarith.c  jdtrans.c  jdatasrc.c 
  jdmaster.c jdinput.c  jdmarker.c jdhuff.c   jdmainct.c 
  jdcoefct.c jdpostct.c jddctmgr.c jidctfst.c jidctflt.c 
  jidctint.c jdsample.c jdcolor.c  jquant1.c  jquant2.c 
  jdmerge.c)

IF(WIN32)
  ADD_DEFINITIONS( -DJPEGSTATIC -D_CRT_SECURE_NO_WARNINGS)
  # TODO(keir): CMake ignores this, but it should disable the following:
  # C:\...\basetsd.h(78): warning C4142: benign redefinition of type
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4142")

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
ENDIF(WIN32)

FILE(GLOB LIBJPEG_H *.h )

ADD_LIBRARY(jpeg ${LIBJPEG_SYSDEP} ${LIBJPEG_COMMON} ${LIBJPEG_CLIB} ${LIBJPEG_DLIB} ${LIBJPEG_H} )

