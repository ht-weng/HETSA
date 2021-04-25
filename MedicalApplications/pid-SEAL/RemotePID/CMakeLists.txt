cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 1.0 LANGUAGES CXX)

# add_executable(pid-controller.o pid-controller.cpp)

# add_executable(bergman.o bergman.cpp)

add_executable(pid-serialize-test pid-serialize-test.cpp)

find_package(SEAL 3.6 REQUIRED)

find_package(Boost REQUIRED)

# Link Microsoft SEAL

# target_link_libraries(pid-controller.o SEAL::seal)

# target_link_libraries(bergman.o SEAL::seal)

target_link_libraries(pid-serialize-test SEAL::seal)

