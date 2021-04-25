cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 1.0 LANGUAGES CXX)

add_executable(client client.cpp)

add_executable(server server.cpp)

find_package(SEAL 3.6 REQUIRED)

find_package(Boost REQUIRED)

# Link Microsoft SEAL
target_link_libraries(client SEAL::seal)

target_link_libraries(server SEAL::seal)


