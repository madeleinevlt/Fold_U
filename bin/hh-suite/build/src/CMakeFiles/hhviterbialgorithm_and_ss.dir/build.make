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
CMAKE_SOURCE_DIR = /home/mikov/Documents/hh-suite

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mikov/Documents/hh-suite/build

# Include any dependencies generated for this target.
include src/CMakeFiles/hhviterbialgorithm_and_ss.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/hhviterbialgorithm_and_ss.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/hhviterbialgorithm_and_ss.dir/flags.make

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/flags.make
src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o: ../src/hhviterbialgorithm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o -c /home/mikov/Documents/hh-suite/src/hhviterbialgorithm.cpp

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/hhviterbialgorithm.cpp > CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.i

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/hhviterbialgorithm.cpp -o CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.s

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.requires:

.PHONY : src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.requires

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.provides: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/hhviterbialgorithm_and_ss.dir/build.make src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.provides.build
.PHONY : src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.provides

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.provides.build: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o


hhviterbialgorithm_and_ss: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o
hhviterbialgorithm_and_ss: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/build.make

.PHONY : hhviterbialgorithm_and_ss

# Rule to build all files generated by this target.
src/CMakeFiles/hhviterbialgorithm_and_ss.dir/build: hhviterbialgorithm_and_ss

.PHONY : src/CMakeFiles/hhviterbialgorithm_and_ss.dir/build

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/requires: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o.requires

.PHONY : src/CMakeFiles/hhviterbialgorithm_and_ss.dir/requires

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/clean:
	cd /home/mikov/Documents/hh-suite/build/src && $(CMAKE_COMMAND) -P CMakeFiles/hhviterbialgorithm_and_ss.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/hhviterbialgorithm_and_ss.dir/clean

src/CMakeFiles/hhviterbialgorithm_and_ss.dir/depend:
	cd /home/mikov/Documents/hh-suite/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mikov/Documents/hh-suite /home/mikov/Documents/hh-suite/src /home/mikov/Documents/hh-suite/build /home/mikov/Documents/hh-suite/build/src /home/mikov/Documents/hh-suite/build/src/CMakeFiles/hhviterbialgorithm_and_ss.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/hhviterbialgorithm_and_ss.dir/depend

