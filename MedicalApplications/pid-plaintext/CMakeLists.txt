cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 1.0 LANGUAGES CXX)

# add_executable(pid pid.cpp)

# add_executable(pid-controller pid-controller.cpp)

# add_executable(bergman bergman.cpp)

add_executable(pid-matlab-plaintext pid-matlab-plaintext.cpp)

find_package(Matlab REQUIRED)

# Link Microsoft SEAL
# target_link_libraries(pid SEAL::seal)

# target_link_libraries(pid-controller SEAL::seal)

# target_link_libraries(bergman SEAL::seal)

target_link_libraries(pid-serialize-test SEAL::seal)

# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)
project ("CMakeProject1" C CXX)
# Find Matlab
FIND_PACKAGE(Matlab REQUIRED)
# Include sub-projects.
INCLUDE_DIRECTORIES("${Matlab_ROOT_DIR}/extern/include")
add_subdirectory ("src")  #Main Project