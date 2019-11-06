#Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH ZIRAN_RELATIVE_PATH ${ZIRAN_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${ZIRAN_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

if(OPENVDB_FOUND)

  add_ziran_executable(mpm mpm.cpp)
  target_compile_options(mpm PRIVATE ${ZIRAN_CXXFLAGS_LTO})
  target_link_libraries(mpm ziran mpmdev tbbmalloc tbbmalloc_proxy ${ZIRAN_CXXFLAGS_LTO})

  if(EXISTS "mpm_tests.cpp")
      add_ziran_executable(mpm_tests mpm_tests.cpp)
      target_link_libraries(mpm_tests gtest ziran mpmdev)
      add_memcheck_test(NAME mpm_tests COMMAND $<TARGET_FILE:mpm_tests>)
  endif()
  
endif()

