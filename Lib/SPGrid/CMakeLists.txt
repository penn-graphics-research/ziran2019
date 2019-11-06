#Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
    file(RELATIVE_PATH ZIRAN_RELATIVE_PATH ${ZIRAN_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    configure_file(${ZIRAN_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_library(spgrid STATIC
        Core/SPGrid_Geometry.cpp
        Core/SPGrid_Utilities.cpp)

target_compile_options(spgrid PUBLIC "-Wno-unused-local-typedefs")

target_include_directories(spgrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)