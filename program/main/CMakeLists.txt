add_executable(pmvs2 pmvs2.cc)

if(WIN32)
target_link_libraries(pmvs2 pmvs_lib image_lib ${PRECOMPILED_LIBS} jpeg)
else(WIN32)
target_link_libraries(pmvs2 pmvs_lib image_lib ${PRECOMPILED_LIBS})
endif(WIN32)


add_executable(genOption genOption.cc)

add_executable(cmvs cmvs.cc)
if(WIN32)
target_link_libraries(cmvs graclus_lib cmvs_lib image_lib ${PRECOMPILED_LIBS} jpeg)
else(WIN32)
target_link_libraries(cmvs cmvs_lib graclus_lib image_lib ${PRECOMPILED_LIBS})
endif(WIN32)
