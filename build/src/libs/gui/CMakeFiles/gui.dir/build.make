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
include src/libs/gui/CMakeFiles/gui.dir/depend.make

# Include the progress variables for this target.
include src/libs/gui/CMakeFiles/gui.dir/progress.make

# Include the compile flags for this target's objects.
include src/libs/gui/CMakeFiles/gui.dir/flags.make

src/libs/gui/CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.o: ../src/libs/gui/imgui_widgets/ImGuizmo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/imgui_widgets/ImGuizmo.cpp

src/libs/gui/CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/imgui_widgets/ImGuizmo.cpp > CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.i

src/libs/gui/CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/imgui_widgets/ImGuizmo.cpp -o CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/application.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/application.cpp.o: ../src/libs/gui/src/application.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/application.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/application.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/application.cpp

src/libs/gui/CMakeFiles/gui.dir/src/application.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/application.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/application.cpp > CMakeFiles/gui.dir/src/application.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/application.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/application.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/application.cpp -o CMakeFiles/gui.dir/src/application.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.o: ../src/libs/gui/src/imgui_impl_glfw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/imgui_impl_glfw.cpp

src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/imgui_impl_glfw.cpp > CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/imgui_impl_glfw.cpp -o CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.o: ../src/libs/gui/src/imgui_impl_opengl3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/imgui_impl_opengl3.cpp

src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/imgui_impl_opengl3.cpp > CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/imgui_impl_opengl3.cpp -o CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/inputstate.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/inputstate.cpp.o: ../src/libs/gui/src/inputstate.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/inputstate.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/inputstate.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/inputstate.cpp

src/libs/gui/CMakeFiles/gui.dir/src/inputstate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/inputstate.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/inputstate.cpp > CMakeFiles/gui.dir/src/inputstate.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/inputstate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/inputstate.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/inputstate.cpp -o CMakeFiles/gui.dir/src/inputstate.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/mesh.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/mesh.cpp.o: ../src/libs/gui/src/mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/mesh.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/mesh.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/mesh.cpp

src/libs/gui/CMakeFiles/gui.dir/src/mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/mesh.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/mesh.cpp > CMakeFiles/gui.dir/src/mesh.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/mesh.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/mesh.cpp -o CMakeFiles/gui.dir/src/mesh.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/model.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/model.cpp.o: ../src/libs/gui/src/model.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/model.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/model.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/model.cpp

src/libs/gui/CMakeFiles/gui.dir/src/model.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/model.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/model.cpp > CMakeFiles/gui.dir/src/model.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/model.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/model.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/model.cpp -o CMakeFiles/gui.dir/src/model.cpp.s

src/libs/gui/CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.o: src/libs/gui/CMakeFiles/gui.dir/flags.make
src/libs/gui/CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.o: ../src/libs/gui/src/shadow_map_fbo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/libs/gui/CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.o"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.o -c /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/shadow_map_fbo.cpp

src/libs/gui/CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.i"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/shadow_map_fbo.cpp > CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.i

src/libs/gui/CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.s"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui/src/shadow_map_fbo.cpp -o CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.s

# Object files for target gui
gui_OBJECTS = \
"CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.o" \
"CMakeFiles/gui.dir/src/application.cpp.o" \
"CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.o" \
"CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.o" \
"CMakeFiles/gui.dir/src/inputstate.cpp.o" \
"CMakeFiles/gui.dir/src/mesh.cpp.o" \
"CMakeFiles/gui.dir/src/model.cpp.o" \
"CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.o"

# External object files for target gui
gui_EXTERNAL_OBJECTS =

src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/imgui_widgets/ImGuizmo.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/application.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_glfw.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/imgui_impl_opengl3.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/inputstate.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/mesh.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/model.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/src/shadow_map_fbo.cpp.o
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/build.make
src/libs/gui/libgui.a: src/libs/gui/CMakeFiles/gui.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/CMM/a2-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX static library libgui.a"
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && $(CMAKE_COMMAND) -P CMakeFiles/gui.dir/cmake_clean_target.cmake
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gui.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/libs/gui/CMakeFiles/gui.dir/build: src/libs/gui/libgui.a

.PHONY : src/libs/gui/CMakeFiles/gui.dir/build

src/libs/gui/CMakeFiles/gui.dir/clean:
	cd /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui && $(CMAKE_COMMAND) -P CMakeFiles/gui.dir/cmake_clean.cmake
.PHONY : src/libs/gui/CMakeFiles/gui.dir/clean

src/libs/gui/CMakeFiles/gui.dir/depend:
	cd /home/ubuntu/CMM/a2-yuliangzhong/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/CMM/a2-yuliangzhong /home/ubuntu/CMM/a2-yuliangzhong/src/libs/gui /home/ubuntu/CMM/a2-yuliangzhong/build /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui /home/ubuntu/CMM/a2-yuliangzhong/build/src/libs/gui/CMakeFiles/gui.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/libs/gui/CMakeFiles/gui.dir/depend

