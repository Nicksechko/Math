# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /opt/clion-2019.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2019.2/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dragon-fly/Math

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dragon-fly/Math/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Math.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Math.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Math.dir/flags.make

CMakeFiles/Math.dir/main.cpp.o: CMakeFiles/Math.dir/flags.make
CMakeFiles/Math.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dragon-fly/Math/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Math.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Math.dir/main.cpp.o -c /home/dragon-fly/Math/main.cpp

CMakeFiles/Math.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Math.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dragon-fly/Math/main.cpp > CMakeFiles/Math.dir/main.cpp.i

CMakeFiles/Math.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Math.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dragon-fly/Math/main.cpp -o CMakeFiles/Math.dir/main.cpp.s

CMakeFiles/Math.dir/fraction.cpp.o: CMakeFiles/Math.dir/flags.make
CMakeFiles/Math.dir/fraction.cpp.o: ../fraction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dragon-fly/Math/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Math.dir/fraction.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Math.dir/fraction.cpp.o -c /home/dragon-fly/Math/fraction.cpp

CMakeFiles/Math.dir/fraction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Math.dir/fraction.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dragon-fly/Math/fraction.cpp > CMakeFiles/Math.dir/fraction.cpp.i

CMakeFiles/Math.dir/fraction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Math.dir/fraction.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dragon-fly/Math/fraction.cpp -o CMakeFiles/Math.dir/fraction.cpp.s

# Object files for target Math
Math_OBJECTS = \
"CMakeFiles/Math.dir/main.cpp.o" \
"CMakeFiles/Math.dir/fraction.cpp.o"

# External object files for target Math
Math_EXTERNAL_OBJECTS =

Math: CMakeFiles/Math.dir/main.cpp.o
Math: CMakeFiles/Math.dir/fraction.cpp.o
Math: CMakeFiles/Math.dir/build.make
Math: CMakeFiles/Math.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dragon-fly/Math/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Math"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Math.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Math.dir/build: Math

.PHONY : CMakeFiles/Math.dir/build

CMakeFiles/Math.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Math.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Math.dir/clean

CMakeFiles/Math.dir/depend:
	cd /home/dragon-fly/Math/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dragon-fly/Math /home/dragon-fly/Math /home/dragon-fly/Math/cmake-build-debug /home/dragon-fly/Math/cmake-build-debug /home/dragon-fly/Math/cmake-build-debug/CMakeFiles/Math.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Math.dir/depend

