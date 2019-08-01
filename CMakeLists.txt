cmake_minimum_required(VERSION 2.8.3)
project(attitude_controller)

find_package(catkin_simple REQUIRED
      mavros_msgs
      roscpp
      geometry_msgs
      )

catkin_simple()

# activate c++ 11
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ELSE()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF()

cs_add_library(mpc_solver
  qpoases/SRC/Bounds.cpp
  qpoases/SRC/Constraints.cpp
  qpoases/SRC/CyclingManager.cpp
  qpoases/SRC/Indexlist.cpp
  qpoases/SRC/MessageHandling.cpp
  qpoases/SRC/QProblem.cpp
  qpoases/SRC/QProblemB.cpp
  qpoases/SRC/SubjectTo.cpp
  qpoases/SRC/Utils.cpp
  qpoases/SRC/EXTRAS/SolutionAnalysis.cpp
  model/nonlinear_mpc_control_export/acado_qpoases_interface.cpp
  model/nonlinear_mpc_control_export/acado_integrator.c
  model/nonlinear_mpc_control_export/acado_solver.c
  model/nonlinear_mpc_control_export/acado_auxiliary_functions.c)

target_include_directories(mpc_solver PUBLIC
  model/nonlinear_mpc_control_export/
  qpoases
  qpoases/INCLUDE
  qpoases/SRC
  )

cs_install()
cs_export()

cs_add_executable(attitude_controller_node
  src/attitude_controller_node.cpp)

target_link_libraries(attitude_controller_node
	mpc_solver)

cs_install()

cs_export()