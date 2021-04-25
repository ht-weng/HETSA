# Minimum CMake required
cmake_minimum_required(VERSION 2.8.12)

# Project
project(protobuf-examples)

set(example /home/justin/Public/HE/HETSA-inside/MedicalApplications/protobuf)

# Find required protobuf package
find_package(protobuf CONFIG REQUIRED)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)
set(${example}_SRCS ${example}.cc)
set(${example}_PROTOS addressbook.proto)


#Executable setup
set(executable_name ${example}_cpp)
add_executable(${executable_name} ${${example}_SRCS} ${${example}_PROTOS})
target_link_libraries(${executable_name} protobuf::libprotobuf)
protobuf_generate(TARGET ${executable_name})