cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 1.0 LANGUAGES CXX)

# add_executable(pid.o pid.cpp)

add_executable(pid-bergman.o pid-bergman.cpp)

find_package(SEAL 3.6 REQUIRED)

# find_package(Boost REQUIRED)

# Link Microsoft SEAL
# target_link_libraries(pid.o SEAL::seal)

target_link_libraries(pid-bergman.o SEAL::seal)

