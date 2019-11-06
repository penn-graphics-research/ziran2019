#Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH ZIRAN_RELATIVE_PATH ${ZIRAN_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${ZIRAN_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

if(OPENVDB_FOUND)
  find_package(Threads)
  file(GLOB MPM_DEV_INSTANTIATIONS "instantiations/*.cpp")
  file(GLOB MPM_FORCE_HELPERS "Force/*.cpp")
  add_library(mpmdev STATIC ${MPM_DEV_INSTANTIATIONS} ${MPM_FORCE_HELPERS})
  target_compile_options(mpmdev PUBLIC "-Wno-unused-but-set-parameter")
  target_compile_options(mpmdev PUBLIC "-Wno-unused-but-set-variable")
  target_compile_options(mpmdev PUBLIC "-Wno-unused-variable")
  target_link_libraries(mpmdev spgrid partio ziran openvdb half tbb metis ${CMAKE_THREAD_LIBS_INIT})

endif()

