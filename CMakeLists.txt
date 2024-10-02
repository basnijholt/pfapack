cmake_minimum_required(VERSION 3.10)
project(pfapack Fortran C)

# RPATH settings
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Find required packages
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Add Fortran source files
file(GLOB FORTRAN_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/external/fortran/*.f"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/fortran/*.f90"
)

set(C_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/skbpf10.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/skbpfa.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/skbtrd.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/skpf10.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/skpfa.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/skpfa_batched.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/sktrd.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/sktrf.c
)

# Create Fortran library
add_library(pfapack_fortran SHARED ${FORTRAN_SOURCES})
target_link_libraries(pfapack_fortran PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
set_target_properties(pfapack_fortran PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    INSTALL_RPATH "@loader_path"
    BUILD_WITH_INSTALL_RPATH TRUE
)

# Create C library
add_library(pfapack_c SHARED ${C_SOURCES})
target_link_libraries(pfapack_c PUBLIC pfapack_fortran)
set_target_properties(pfapack_c PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    INSTALL_RPATH "@loader_path"
    BUILD_WITH_INSTALL_RPATH TRUE
)

# Set include directories for C library
target_include_directories(pfapack_c PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/fortran
    ${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface
)

# Install libraries
install(TARGETS pfapack_fortran pfapack_c
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)

# Copy header files
file(GLOB HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/external/c_interface/*.h")
install(FILES ${HEADER_FILES} DESTINATION include)
