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
include lib/ffindex/src/CMakeFiles/ffindex.dir/depend.make

# Include the progress variables for this target.
include lib/ffindex/src/CMakeFiles/ffindex.dir/progress.make

# Include the compile flags for this target's objects.
include lib/ffindex/src/CMakeFiles/ffindex.dir/flags.make

lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o: lib/ffindex/src/CMakeFiles/ffindex.dir/flags.make
lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o: ../lib/ffindex/src/ffindex.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/ffindex.dir/ffindex.c.o   -c /home/mikov/Documents/hh-suite/lib/ffindex/src/ffindex.c

lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ffindex.dir/ffindex.c.i"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/mikov/Documents/hh-suite/lib/ffindex/src/ffindex.c > CMakeFiles/ffindex.dir/ffindex.c.i

lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ffindex.dir/ffindex.c.s"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/mikov/Documents/hh-suite/lib/ffindex/src/ffindex.c -o CMakeFiles/ffindex.dir/ffindex.c.s

lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.requires:

.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.requires

lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.provides: lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.requires
	$(MAKE) -f lib/ffindex/src/CMakeFiles/ffindex.dir/build.make lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.provides.build
.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.provides

lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.provides.build: lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o


lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o: lib/ffindex/src/CMakeFiles/ffindex.dir/flags.make
lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o: ../lib/ffindex/src/ffutil.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/ffindex.dir/ffutil.c.o   -c /home/mikov/Documents/hh-suite/lib/ffindex/src/ffutil.c

lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/ffindex.dir/ffutil.c.i"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/mikov/Documents/hh-suite/lib/ffindex/src/ffutil.c > CMakeFiles/ffindex.dir/ffutil.c.i

lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/ffindex.dir/ffutil.c.s"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/mikov/Documents/hh-suite/lib/ffindex/src/ffutil.c -o CMakeFiles/ffindex.dir/ffutil.c.s

lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.requires:

.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.requires

lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.provides: lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.requires
	$(MAKE) -f lib/ffindex/src/CMakeFiles/ffindex.dir/build.make lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.provides.build
.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.provides

lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.provides.build: lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o


# Object files for target ffindex
ffindex_OBJECTS = \
"CMakeFiles/ffindex.dir/ffindex.c.o" \
"CMakeFiles/ffindex.dir/ffutil.c.o"

# External object files for target ffindex
ffindex_EXTERNAL_OBJECTS =

lib/libffindex.a: lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o
lib/libffindex.a: lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o
lib/libffindex.a: lib/ffindex/src/CMakeFiles/ffindex.dir/build.make
lib/libffindex.a: lib/ffindex/src/CMakeFiles/ffindex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library ../../libffindex.a"
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && $(CMAKE_COMMAND) -P CMakeFiles/ffindex.dir/cmake_clean_target.cmake
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ffindex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/ffindex/src/CMakeFiles/ffindex.dir/build: lib/libffindex.a

.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/build

lib/ffindex/src/CMakeFiles/ffindex.dir/requires: lib/ffindex/src/CMakeFiles/ffindex.dir/ffindex.c.o.requires
lib/ffindex/src/CMakeFiles/ffindex.dir/requires: lib/ffindex/src/CMakeFiles/ffindex.dir/ffutil.c.o.requires

.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/requires

lib/ffindex/src/CMakeFiles/ffindex.dir/clean:
	cd /home/mikov/Documents/hh-suite/build/lib/ffindex/src && $(CMAKE_COMMAND) -P CMakeFiles/ffindex.dir/cmake_clean.cmake
.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/clean

lib/ffindex/src/CMakeFiles/ffindex.dir/depend:
	cd /home/mikov/Documents/hh-suite/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mikov/Documents/hh-suite /home/mikov/Documents/hh-suite/lib/ffindex/src /home/mikov/Documents/hh-suite/build /home/mikov/Documents/hh-suite/build/lib/ffindex/src /home/mikov/Documents/hh-suite/build/lib/ffindex/src/CMakeFiles/ffindex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/ffindex/src/CMakeFiles/ffindex.dir/depend

