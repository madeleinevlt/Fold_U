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
include src/CMakeFiles/hhsearch.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/hhsearch.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/hhsearch.dir/flags.make

src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o: src/CMakeFiles/hhsearch.dir/flags.make
src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o: ../src/hhsearch_app.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o -c /home/mikov/Documents/hh-suite/src/hhsearch_app.cpp

src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hhsearch.dir/hhsearch_app.cpp.i"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mikov/Documents/hh-suite/src/hhsearch_app.cpp > CMakeFiles/hhsearch.dir/hhsearch_app.cpp.i

src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hhsearch.dir/hhsearch_app.cpp.s"
	cd /home/mikov/Documents/hh-suite/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mikov/Documents/hh-suite/src/hhsearch_app.cpp -o CMakeFiles/hhsearch.dir/hhsearch_app.cpp.s

src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.requires:

.PHONY : src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.requires

src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.provides: src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/hhsearch.dir/build.make src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.provides.build
.PHONY : src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.provides

src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.provides.build: src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o


# Object files for target hhsearch
hhsearch_OBJECTS = \
"CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o"

# External object files for target hhsearch
hhsearch_EXTERNAL_OBJECTS = \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhblits.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhdecl.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhhit.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhmatrices.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhsearch.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhalign.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhhitlist.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhposteriordecoder.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhutil.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/util.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhalignment.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhforwardalgorithm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhhmm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhposteriordecoderrunner.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhviterbialgorithm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhfullalignment.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhhmmsimd.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhposteriormatrix.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhviterbi.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhbacktracemac.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhmacalgorithm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhprefilter.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhviterbimatrix.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhbackwardalgorithm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhdatabase.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhhalfalignment.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhviterbirunner.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/HH_OBJECTS.dir/hhfunc.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/hhviterbialgorithm_with_celloff.dir/hhviterbialgorithm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o" \
"/home/mikov/Documents/hh-suite/build/src/CMakeFiles/hhviterbialgorithm_with_celloff_and_ss.dir/hhviterbialgorithm.cpp.o"

bin/hhsearch: src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/aa.cc.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/as.cc.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/assert_helpers.cc.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/blosum_matrix.cc.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/getopt_pp.cc.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/log.cc.o
bin/hhsearch: src/CMakeFiles/CS_OBJECTS.dir/cs/application.cc.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhblits.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhdecl.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhhit.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhmatrices.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhsearch.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhalign.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhhitlist.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhposteriordecoder.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhutil.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/util.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhalignment.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhforwardalgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhhmm.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhposteriordecoderrunner.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhviterbialgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhfullalignment.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhhmmsimd.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhposteriormatrix.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhviterbi.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhbacktracemac.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhmacalgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhprefilter.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhviterbimatrix.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhbackwardalgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhdatabase.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhhalfalignment.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhviterbirunner.cpp.o
bin/hhsearch: src/CMakeFiles/HH_OBJECTS.dir/hhfunc.cpp.o
bin/hhsearch: src/CMakeFiles/hhviterbialgorithm_with_celloff.dir/hhviterbialgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/hhviterbialgorithm_and_ss.dir/hhviterbialgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/hhviterbialgorithm_with_celloff_and_ss.dir/hhviterbialgorithm.cpp.o
bin/hhsearch: src/CMakeFiles/hhsearch.dir/build.make
bin/hhsearch: lib/libffindex.a
bin/hhsearch: src/CMakeFiles/hhsearch.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mikov/Documents/hh-suite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/hhsearch"
	cd /home/mikov/Documents/hh-suite/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hhsearch.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/hhsearch.dir/build: bin/hhsearch

.PHONY : src/CMakeFiles/hhsearch.dir/build

src/CMakeFiles/hhsearch.dir/requires: src/CMakeFiles/hhsearch.dir/hhsearch_app.cpp.o.requires

.PHONY : src/CMakeFiles/hhsearch.dir/requires

src/CMakeFiles/hhsearch.dir/clean:
	cd /home/mikov/Documents/hh-suite/build/src && $(CMAKE_COMMAND) -P CMakeFiles/hhsearch.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/hhsearch.dir/clean

src/CMakeFiles/hhsearch.dir/depend:
	cd /home/mikov/Documents/hh-suite/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mikov/Documents/hh-suite /home/mikov/Documents/hh-suite/src /home/mikov/Documents/hh-suite/build /home/mikov/Documents/hh-suite/build/src /home/mikov/Documents/hh-suite/build/src/CMakeFiles/hhsearch.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/hhsearch.dir/depend

