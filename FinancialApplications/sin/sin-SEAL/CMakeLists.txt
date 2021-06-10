cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 1.0 LANGUAGES CXX)

add_executable(sin.o sin.cpp)

find_package(SEAL 3.6 REQUIRED)

target_link_libraries(sin.o SEAL::seal)

