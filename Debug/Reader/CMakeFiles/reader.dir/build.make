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
include Reader/CMakeFiles/reader.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Reader/CMakeFiles/reader.dir/compiler_depend.make

# Include the progress variables for this target.
include Reader/CMakeFiles/reader.dir/progress.make

# Include the compile flags for this target's objects.
include Reader/CMakeFiles/reader.dir/flags.make

Reader/CMakeFiles/reader.dir/src/reader.o: Reader/CMakeFiles/reader.dir/flags.make
Reader/CMakeFiles/reader.dir/src/reader.o: ../Reader/src/reader.c
Reader/CMakeFiles/reader.dir/src/reader.o: Reader/CMakeFiles/reader.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/ze/Tool/EPM/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object Reader/CMakeFiles/reader.dir/src/reader.o"
	cd /media/ze/Tool/EPM/Debug/Reader && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT Reader/CMakeFiles/reader.dir/src/reader.o -MF CMakeFiles/reader.dir/src/reader.o.d -o CMakeFiles/reader.dir/src/reader.o -c /media/ze/Tool/EPM/Reader/src/reader.c

Reader/CMakeFiles/reader.dir/src/reader.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/reader.dir/src/reader.i"
	cd /media/ze/Tool/EPM/Debug/Reader && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /media/ze/Tool/EPM/Reader/src/reader.c > CMakeFiles/reader.dir/src/reader.i

Reader/CMakeFiles/reader.dir/src/reader.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/reader.dir/src/reader.s"
	cd /media/ze/Tool/EPM/Debug/Reader && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /media/ze/Tool/EPM/Reader/src/reader.c -o CMakeFiles/reader.dir/src/reader.s

Reader/CMakeFiles/reader.dir/src/reader_utils.o: Reader/CMakeFiles/reader.dir/flags.make
Reader/CMakeFiles/reader.dir/src/reader_utils.o: ../Reader/src/reader_utils.c
Reader/CMakeFiles/reader.dir/src/reader_utils.o: Reader/CMakeFiles/reader.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/ze/Tool/EPM/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object Reader/CMakeFiles/reader.dir/src/reader_utils.o"
	cd /media/ze/Tool/EPM/Debug/Reader && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT Reader/CMakeFiles/reader.dir/src/reader_utils.o -MF CMakeFiles/reader.dir/src/reader_utils.o.d -o CMakeFiles/reader.dir/src/reader_utils.o -c /media/ze/Tool/EPM/Reader/src/reader_utils.c

Reader/CMakeFiles/reader.dir/src/reader_utils.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/reader.dir/src/reader_utils.i"
	cd /media/ze/Tool/EPM/Debug/Reader && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /media/ze/Tool/EPM/Reader/src/reader_utils.c > CMakeFiles/reader.dir/src/reader_utils.i

Reader/CMakeFiles/reader.dir/src/reader_utils.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/reader.dir/src/reader_utils.s"
	cd /media/ze/Tool/EPM/Debug/Reader && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /media/ze/Tool/EPM/Reader/src/reader_utils.c -o CMakeFiles/reader.dir/src/reader_utils.s

# Object files for target reader
reader_OBJECTS = \
"CMakeFiles/reader.dir/src/reader.o" \
"CMakeFiles/reader.dir/src/reader_utils.o"

# External object files for target reader
reader_EXTERNAL_OBJECTS =

Reader/libreader.a: Reader/CMakeFiles/reader.dir/src/reader.o
Reader/libreader.a: Reader/CMakeFiles/reader.dir/src/reader_utils.o
Reader/libreader.a: Reader/CMakeFiles/reader.dir/build.make
Reader/libreader.a: Reader/CMakeFiles/reader.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/ze/Tool/EPM/Debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libreader.a"
	cd /media/ze/Tool/EPM/Debug/Reader && $(CMAKE_COMMAND) -P CMakeFiles/reader.dir/cmake_clean_target.cmake
	cd /media/ze/Tool/EPM/Debug/Reader && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/reader.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Reader/CMakeFiles/reader.dir/build: Reader/libreader.a
.PHONY : Reader/CMakeFiles/reader.dir/build

Reader/CMakeFiles/reader.dir/clean:
	cd /media/ze/Tool/EPM/Debug/Reader && $(CMAKE_COMMAND) -P CMakeFiles/reader.dir/cmake_clean.cmake
.PHONY : Reader/CMakeFiles/reader.dir/clean

Reader/CMakeFiles/reader.dir/depend:
	cd /media/ze/Tool/EPM/Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/ze/Tool/EPM /media/ze/Tool/EPM/Reader /media/ze/Tool/EPM/Debug /media/ze/Tool/EPM/Debug/Reader /media/ze/Tool/EPM/Debug/Reader/CMakeFiles/reader.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Reader/CMakeFiles/reader.dir/depend

