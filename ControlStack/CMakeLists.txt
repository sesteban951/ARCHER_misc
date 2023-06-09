cmake_minimum_required(VERSION 3.10)

# set the project name
project(hopper)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -pthread")

SET(MUJOCO_DIR /home/amber/repos/mujoco-2.2.1 CACHE PATH "Path to Mujoco dir")
SET(PINOCCHIO_DIR /home/amber/repos/pinocchio_install CACHE PATH "Path to pinocchio dir")
SET(MANIF_DIR /home/amber/repos/manif_install CACHE PATH "Path to manif dir")
list(APPEND CMAKE_PREFIX_PATH ${PINOCCHIO_DIR})

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(glfw3 REQUIRED)
find_package(pinocchio REQUIRED HINTS ${PINOCCHIO_DIR}/lib)
find_package(manif REQUIRED HINTS ${MANIF_DIR})
find_package(roscpp REQUIRED)
find_package(geometry_msgs REQUIRED)
include_directories(${PINOCCHIO_DIR}/include)

include_directories(${MUJOCO_DIR}/include)
include_directories(${MUJOCO_DIR}/include/mujoco)
find_library(MUJOCO_LIB mujoco REQUIRED HINTS ${MUJOCO_DIR}/lib)

add_library(mujoco SHARED IMPORTED GLOBAL)
set_property(TARGET mujoco PROPERTY IMPORTED_LOCATION ${CMAKE_INCLUDE_DIR}/lib/libmujoco.so)

add_library(MPC SHARED
    src/Hopper.cpp
    src/MPC.cpp
)


add_library(QP SHARED
    src/Hopper.cpp
    src/QP.cpp
    ) ###

target_link_libraries(MPC Eigen3::Eigen)
target_link_libraries(MPC osqp::osqp)
target_link_libraries(MPC OsqpEigen::OsqpEigen)
target_include_directories(MPC PUBLIC ${manif_INCLUDE_DIRS})
target_link_libraries(MPC ${PINOCCHIO_LIBRARIES})
target_link_directories(MPC PUBLIC ${PINOCCHIO_INCLUDE_DIRS})


target_link_libraries(QP Eigen3::Eigen) ###
target_link_libraries(QP osqp::osqp) ###
target_link_libraries(QP OsqpEigen::OsqpEigen) ###
target_include_directories(QP PUBLIC ${manif_INCLUDE_DIRS}) ###
#target_link_libraries(QP ${PINOCCHIO_LIBRARIES}) ###
#target_link_directories(QP PUBLIC ${PINOCCHIO_INCLUDE_DIRS}) ###

# add the executable
add_executable(hopper_ctrl src/Controller.cpp)
target_link_libraries(hopper_ctrl MPC)
target_link_libraries(hopper_ctrl QP) ###
target_link_libraries(hopper_ctrl yaml-cpp)

add_executable(hopper_ctrl_hardware src/Controller_Hardware.cpp)
target_link_libraries(hopper_ctrl_hardware MPC)
target_link_libraries(hopper_ctrl_hardware QP) ###
target_link_libraries(hopper_ctrl_hardware yaml-cpp)
target_include_directories(hopper_ctrl_hardware PUBLIC ${roscpp_INCLUDE_DIRS})
target_link_libraries(hopper_ctrl_hardware ${roscpp_LIBRARIES})

add_executable(hopper_sim src/Simulator.cpp)
target_link_libraries(hopper_sim PUBLIC Eigen3::Eigen)
target_link_libraries(hopper_sim PUBLIC glfw)
target_link_libraries(hopper_sim PUBLIC yaml-cpp)
target_link_libraries(hopper_sim PUBLIC ${MUJOCO_LIB})         ###############################
#target_link_libraries(hopper_sim PUBLIC ${MUJOCO_DIR})
target_link_directories(hopper_sim PUBLIC ${MUJOCO_PATH})





