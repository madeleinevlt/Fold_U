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
include src/CMakeFiles/CS_OBJECTS.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/CS_OBJECTS.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/CS_OBJECTS.dir/flags.make

src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o: ../src/cs/aa.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o -c /home/mikov/Documents/hh-suite/src/cs/aa.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/aa.cc > CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/aa.cc -o CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o


src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o: ../src/cs/as.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o -c /home/mikov/Documents/hh-suite/src/cs/as.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/as.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/as.cc > CMakeFiles/CS_OBJECTS.dir/cs/as.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/as.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/as.cc -o CMakeFiles/CS_OBJECTS.dir/cs/as.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o


src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o: ../src/cs/assert_helpers.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o -c /home/mikov/Documents/hh-suite/src/cs/assert_helpers.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/assert_helpers.cc > CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/assert_helpers.cc -o CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o


src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o: ../src/cs/blosum_matrix.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o -c /home/mikov/Documents/hh-suite/src/cs/blosum_matrix.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/blosum_matrix.cc > CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/blosum_matrix.cc -o CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o


src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o: ../src/cs/getopt_pp.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o -c /home/mikov/Documents/hh-suite/src/cs/getopt_pp.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/getopt_pp.cc > CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/getopt_pp.cc -o CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o


src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o: ../src/cs/log.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o -c /home/mikov/Documents/hh-suite/src/cs/log.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/log.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/log.cc > CMakeFiles/CS_OBJECTS.dir/cs/log.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/log.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/log.cc -o CMakeFiles/CS_OBJECTS.dir/cs/log.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o


src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o: src/CMakeFiles/CS_OBJECTS.dir/flags.make
src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o: ../src/cs/application.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o -c /home/mikov/Documents/hh-suite/src/cs/application.cc

src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CS_OBJECTS.dir/cs/application.cc.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/cs/application.cc > CMakeFiles/CS_OBJECTS.dir/cs/application.cc.i

src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CS_OBJECTS.dir/cs/application.cc.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/cs/application.cc -o CMakeFiles/CS_OBJECTS.dir/cs/application.cc.s

src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.requires:

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.requires

src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.provides: src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.requires
	$(MAKE) -f src/CMakeFiles/CS_OBJECTS.dir/build.make src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.provides.build
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.provides

src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.provides.build: src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o


CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o
CS_OBJECTS: src/CMakeFiles/CS_OBJECTS.dir/build.make

.PHONY : CS_OBJECTS

# Rule to build all files generated by this target.
src/CMakeFiles/CS_OBJECTS.dir/build: CS_OBJECTS

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/build

src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o.requires
src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o.requires
src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o.requires
src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o.requires
src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o.requires
src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o.requires
src/CMakeFiles/CS_OBJECTS.dir/requires: src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o.requires

.PHONY : src/CMakeFiles/CS_OBJECTS.dir/requires

src/CMakeFiles/CS_OBJECTS.dir/clean:
	cd /home/mikov/Documents/hh-suite/build/src && $(CMAKE_COMMAND) -P CMakeFiles/CS_OBJECTS.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/clean

src/CMakeFiles/CS_OBJECTS.dir/depend:
	cd /home/mikov/Documents/hh-suite/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mikov/Documents/hh-suite /home/mikov/Documents/hh-suite/src /home/mikov/Documents/hh-suite/build /home/mikov/Documents/hh-suite/build/src /home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/CS_OBJECTS.dir/depend

