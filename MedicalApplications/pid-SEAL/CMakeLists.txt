cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 3.4.4 LANGUAGES CXX)

add_executable(pid pid)

find_package(SEAL 3.4.4 REQUIRED)

# Link Microsoft SEAL
target_link_libraries(pid SEAL::seal)