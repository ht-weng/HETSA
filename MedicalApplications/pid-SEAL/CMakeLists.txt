cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 1.0 LANGUAGES CXX)

add_executable(pid pid.cpp)

add_executable(pid_client pid_client.cpp)

add_executable(pid_server pid_server.cpp)

find_package(SEAL 3.6 REQUIRED)

# Link Microsoft SEAL
target_link_libraries(pid SEAL::seal)

target_link_libraries(pid_client SEAL::seal)

target_link_libraries(pid_server SEAL::seal)
