find_package(TBB QUIET)
if(TBB_FOUND)
  add_library(tbb SHARED IMPORTED GLOBAL)
  set_target_properties(tbb PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS}
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS}
    IMPORTED_LOCATION ${TBB_LIBRARIES}
    )
  add_library(tbbmalloc SHARED IMPORTED GLOBAL)
  set_target_properties(tbbmalloc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS}
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS}
    IMPORTED_LOCATION ${TBB_MALLOC}
    )
  add_library(tbbmalloc_proxy SHARED IMPORTED GLOBAL)
  set_target_properties(tbbmalloc_proxy PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS}
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${TBB_INCLUDE_DIRS}
    IMPORTED_LOCATION ${TBB_MALLOC_PROXY}
    )
endif()

find_package(Half QUIET)
if(HALF_FOUND)
  add_library(half SHARED IMPORTED GLOBAL)
  set_target_properties(half PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${HALF_INCLUDE_DIRS}
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${HALF_INCLUDE_DIRS}
    IMPORTED_LOCATION ${HALF_LIBRARIES}
    )
endif()

find_package(OpenVDB QUIET)
if(OPENVDB_FOUND AND TBB_FOUND AND HALF_FOUND)
  add_library(openvdb SHARED IMPORTED GLOBAL)
  set_target_properties(openvdb PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${OPENVDB_INCLUDE_DIRS}
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${OPENVDB_INCLUDE_DIRS}
    IMPORTED_LOCATION ${OPENVDB_LIBRARIES}
    INTERFACE_COMPILE_DEFINITIONS "OPENVDB_USE_BLOSC;ZIRAN_WITH_VDB"
    )
endif()

find_package(CUDA QUIET)
if(CUDA_FOUND)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};${GENCODE};-lineinfo;--expt-extended-lambda;--expt-relaxed-constexpr;--default-stream per-thread PARENT_SCOPE)
set(CUDA_SEPARABLE_COMPILATION ON PARENT_SCOPE)
endif()

find_package(Git)
# Enable DownloadProject CMake module
include(DownloadProject)
include(ExternalProject)

ExternalProject_Add(voro_external
        URL http://math.lbl.gov/voro++/download/dir/voro++-0.4.6.tar.gz
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/voro_external/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make
        INSTALL_COMMAND ""
        )
MESSAGE( STATUS "variable is: " ${CMAKE_CURRENT_BINARY_DIR} )
SET(VORO_EXTERNAL_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/voro_external/src" PARENT_SCOPE)


download_project(
  PROJ kdtree
  GIT_REPOSITORY https://github.com/nvmd/libkdtree.git
  GIT_TAG 7bb7e830d6899214e9e896f920483ddb39c43f7b
  UPDATE_DISCONNECTED 1
  )
add_library(kdtree INTERFACE)
target_include_directories(kdtree SYSTEM INTERFACE ${kdtree_SOURCE_DIR})

download_project(
  PROJ partio
  GIT_REPOSITORY https://github.com/wdas/partio.git
  GIT_TAG v1.1.0
  UPDATE_DISCONNECTED 1
  PATCH_COMMAND ${GIT_EXECUTABLE} apply ${CMAKE_CURRENT_SOURCE_DIR}/partio.patch
  )

add_subdirectory(${partio_SOURCE_DIR} ${partio_BINARY_DIR})
set_target_properties(partio PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${partio_SOURCE_DIR}/src/lib"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${partio_SOURCE_DIR}/src/lib"
    INTERFACE_LINK_LIBRARIES "z"
)

download_project(
  PROJ Tick
  GIT_REPOSITORY https://github.com/pfultz2/Tick.git
  GIT_TAG b82af5461e95ddc6b652ac2456976cc912ab2f01
  UPDATE_DISCONNECTED 1
  )

add_library(Tick INTERFACE)
target_include_directories(Tick SYSTEM INTERFACE ${Tick_SOURCE_DIR})

find_package(Metis)
if(METIS_FOUND)
  add_library(metis SHARED IMPORTED GLOBAL)
  if(METIS_INCLUDE_DIRS)
  set_target_properties(metis PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIRS}
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIRS}
    IMPORTED_LOCATION ${METIS_LIBRARIES}
    )
  else()
  set_target_properties(metis PROPERTIES
    IMPORTED_LOCATION ${METIS_LIBRARIES}
    )
  endif()
endif()

# allow visibility in more files
set(TBB_FOUND ${TBB_FOUND} PARENT_SCOPE)
set(HALF_FOUND ${HALF_FOUND} PARENT_SCOPE)
set(OPENVDB_FOUND ${OPENVDB_FOUND} PARENT_SCOPE)
set(METIS_FOUND ${METIS_FOUND} PARENT_SCOPE)
set(CUDA_FOUND ${CUDA_FOUND} PARENT_SCOPE)
