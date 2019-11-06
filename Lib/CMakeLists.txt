if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH ZIRAN_RELATIVE_PATH ${ZIRAN_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${ZIRAN_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_subdirectory(Ziran)
add_subdirectory(MPM)
add_subdirectory(SPGrid)
