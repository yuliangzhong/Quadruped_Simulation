# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ubuntu/CMM/a2-yuliangzhong

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/CMM/a2-yuliangzhong/build

# Include any dependencies generated for this target.
include src/libs/utils/CMakeFiles/utils.dir/depend.make

# Include the progress variables for this target.
include src/libs/utils/CMakeFiles/utils.dir/progress.make

# Include the compile flags for this target's objects.
include src/libs/utils/CMakeFiles/utils.dir/flags.make

src/libs/utils/CMakeFiles/utils.dir/src/logger.cpp.o: src/libs/utils/CMakeFiles/utils.dir/flags.make
src/libs/utils/CMakeFiles/utils.dir/src/logger.cpp.o: ../src/libs/utils/src/logger.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/libs/utils/CMakeFiles/utils.dir/src/logger.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/src/logger.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils/src/logger.cpp

src/libs/utils/CMakeFiles/utils.dir/src/logger.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/src/logger.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils/src/logger.cpp > CMakeFiles/utils.dir/src/logger.cpp.i

src/libs/utils/CMakeFiles/utils.dir/src/logger.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/src/logger.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils/src/logger.cpp -o CMakeFiles/utils.dir/src/logger.cpp.s

src/libs/utils/CMakeFiles/utils.dir/src/timer.cpp.o: src/libs/utils/CMakeFiles/utils.dir/flags.make
src/libs/utils/CMakeFiles/utils.dir/src/timer.cpp.o: ../src/libs/utils/src/timer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/libs/utils/CMakeFiles/utils.dir/src/timer.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utils.dir/src/timer.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils/src/timer.cpp

src/libs/utils/CMakeFiles/utils.dir/src/timer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utils.dir/src/timer.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils/src/timer.cpp > CMakeFiles/utils.dir/src/timer.cpp.i

src/libs/utils/CMakeFiles/utils.dir/src/timer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utils.dir/src/timer.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils/src/timer.cpp -o CMakeFiles/utils.dir/src/timer.cpp.s

# Object files for target utils
utils_OBJECTS = \
"CMakeFiles/utils.dir/src/logger.cpp.o" \
"CMakeFiles/utils.dir/src/timer.cpp.o"

# External object files for target utils
utils_EXTERNAL_OBJECTS =

src/libs/utils/libutils.a: src/libs/utils/CMakeFiles/utils.dir/src/logger.cpp.o
src/libs/utils/libutils.a: src/libs/utils/CMakeFiles/utils.dir/src/timer.cpp.o
src/libs/utils/libutils.a: src/libs/utils/CMakeFiles/utils.dir/build.make
src/libs/utils/libutils.a: src/libs/utils/CMakeFiles/utils.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libutils.a"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean_target.cmake
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/utils.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/libs/utils/CMakeFiles/utils.dir/build: src/libs/utils/libutils.a

.PHONY : src/libs/utils/CMakeFiles/utils.dir/build

src/libs/utils/CMakeFiles/utils.dir/clean:
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils && $(CMAKE_COMMAND) -P CMakeFiles/utils.dir/cmake_clean.cmake
.PHONY : src/libs/utils/CMakeFiles/utils.dir/clean

src/libs/utils/CMakeFiles/utils.dir/depend:
	cd /home/ubuntu/CMM/a2-yuliangzhong/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/CMM/a2-yuliangzhong /home/ubuntu/CMM/a2-yuliangzhong/src/libs/utils /home/ubuntu/CMM/a2-yuliangzhong/build /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/utils/CMakeFiles/utils.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/libs/utils/CMakeFiles/utils.dir/depend

