
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

# Using the vcpkg submodule
set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/external/vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "Vcpkg toolchain file")

project(SimulationSkeleton LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Default to Debug
if(NOT CMAKE_BUILD_TYPE)
    message("Build type not specified, using Debug")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type not specified, using Debug" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(imgui CONFIG REQUIRED)
find_package(SFML CONFIG REQUIRED COMPONENTS graphics )
find_package(ImGui-SFML CONFIG REQUIRED)

add_executable(ParticleLife
    main.cpp
)

target_link_libraries(ParticleLife
    imgui::imgui
    sfml-graphics
    ImGui-SFML::ImGui-SFML
)

# Copy imgui.ini
add_custom_command(
    TARGET ParticleLife POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/imgui.ini
        ${CMAKE_CURRENT_BINARY_DIR})
