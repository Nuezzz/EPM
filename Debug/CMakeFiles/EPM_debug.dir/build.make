# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/ze/Tool/EPM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/ze/Tool/EPM/Debug

# Include any dependencies generated for this target.
include CMakeFiles/EPM_debug.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/EPM_debug.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/EPM_debug.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/EPM_debug.dir/flags.make

CMakeFiles/EPM_debug.dir/Main/src/main.o: CMakeFiles/EPM_debug.dir/flags.make
CMakeFiles/EPM_debug.dir/Main/src/main.o: ../Main/src/main.c
CMakeFiles/EPM_debug.dir/Main/src/main.o: CMakeFiles/EPM_debug.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/ze/Tool/EPM/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/EPM_debug.dir/Main/src/main.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/EPM_debug.dir/Main/src/main.o -MF CMakeFiles/EPM_debug.dir/Main/src/main.o.d -o CMakeFiles/EPM_debug.dir/Main/src/main.o -c /media/ze/Tool/EPM/Main/src/main.c

CMakeFiles/EPM_debug.dir/Main/src/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/EPM_debug.dir/Main/src/main.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /media/ze/Tool/EPM/Main/src/main.c > CMakeFiles/EPM_debug.dir/Main/src/main.i

CMakeFiles/EPM_debug.dir/Main/src/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/EPM_debug.dir/Main/src/main.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /media/ze/Tool/EPM/Main/src/main.c -o CMakeFiles/EPM_debug.dir/Main/src/main.s

# Object files for target EPM_debug
EPM_debug_OBJECTS = \
"CMakeFiles/EPM_debug.dir/Main/src/main.o"

# External object files for target EPM_debug
EPM_debug_EXTERNAL_OBJECTS =

bin/EPM_debug: CMakeFiles/EPM_debug.dir/Main/src/main.o
bin/EPM_debug: CMakeFiles/EPM_debug.dir/build.make
bin/EPM_debug: Control/libcontrol.a
bin/EPM_debug: Lattice/liblattice.a
bin/EPM_debug: Reader/libreader.a
bin/EPM_debug: List/liblist.a
bin/EPM_debug: Control/libcontrol.a
bin/EPM_debug: CMakeFiles/EPM_debug.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/ze/Tool/EPM/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable bin/EPM_debug"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/EPM_debug.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/EPM_debug.dir/build: bin/EPM_debug
.PHONY : CMakeFiles/EPM_debug.dir/build

CMakeFiles/EPM_debug.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/EPM_debug.dir/cmake_clean.cmake
.PHONY : CMakeFiles/EPM_debug.dir/clean

CMakeFiles/EPM_debug.dir/depend:
	cd /media/ze/Tool/EPM/Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/ze/Tool/EPM /media/ze/Tool/EPM /media/ze/Tool/EPM/Debug /media/ze/Tool/EPM/Debug /media/ze/Tool/EPM/Debug/CMakeFiles/EPM_debug.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/EPM_debug.dir/depend
