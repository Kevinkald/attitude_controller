# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kevin/catkin_ws/src/attitude_controller/model

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kevin/catkin_ws/src/attitude_controller/model/build

# Include any dependencies generated for this target.
include CMakeFiles/quadNMPC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/quadNMPC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/quadNMPC.dir/flags.make

CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o: CMakeFiles/quadNMPC.dir/flags.make
CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o: ../quad_nmpc_sim.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kevin/catkin_ws/src/attitude_controller/model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o -c /home/kevin/catkin_ws/src/attitude_controller/model/quad_nmpc_sim.cpp

CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kevin/catkin_ws/src/attitude_controller/model/quad_nmpc_sim.cpp > CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.i

CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kevin/catkin_ws/src/attitude_controller/model/quad_nmpc_sim.cpp -o CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.s

CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.requires:

.PHONY : CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.requires

CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.provides: CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.requires
	$(MAKE) -f CMakeFiles/quadNMPC.dir/build.make CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.provides.build
.PHONY : CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.provides

CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.provides.build: CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o


# Object files for target quadNMPC
quadNMPC_OBJECTS = \
"CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o"

# External object files for target quadNMPC
quadNMPC_EXTERNAL_OBJECTS =

../quadNMPC: CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o
../quadNMPC: CMakeFiles/quadNMPC.dir/build.make
../quadNMPC: /home/kevin/ACADOtoolkit/build/lib/libacado_toolkit_s.so
../quadNMPC: CMakeFiles/quadNMPC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/kevin/catkin_ws/src/attitude_controller/model/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../quadNMPC"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quadNMPC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/quadNMPC.dir/build: ../quadNMPC

.PHONY : CMakeFiles/quadNMPC.dir/build

CMakeFiles/quadNMPC.dir/requires: CMakeFiles/quadNMPC.dir/quad_nmpc_sim.cpp.o.requires

.PHONY : CMakeFiles/quadNMPC.dir/requires

CMakeFiles/quadNMPC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/quadNMPC.dir/cmake_clean.cmake
.PHONY : CMakeFiles/quadNMPC.dir/clean

CMakeFiles/quadNMPC.dir/depend:
	cd /home/kevin/catkin_ws/src/attitude_controller/model/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kevin/catkin_ws/src/attitude_controller/model /home/kevin/catkin_ws/src/attitude_controller/model /home/kevin/catkin_ws/src/attitude_controller/model/build /home/kevin/catkin_ws/src/attitude_controller/model/build /home/kevin/catkin_ws/src/attitude_controller/model/build/CMakeFiles/quadNMPC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/quadNMPC.dir/depend
