# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /snap/clion/114/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/114/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/tarea2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tarea2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tarea2.dir/flags.make

CMakeFiles/tarea2.dir/strassen.cpp.o: CMakeFiles/tarea2.dir/flags.make
CMakeFiles/tarea2.dir/strassen.cpp.o: ../strassen.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tarea2.dir/strassen.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tarea2.dir/strassen.cpp.o -c "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/strassen.cpp"

CMakeFiles/tarea2.dir/strassen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tarea2.dir/strassen.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/strassen.cpp" > CMakeFiles/tarea2.dir/strassen.cpp.i

CMakeFiles/tarea2.dir/strassen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tarea2.dir/strassen.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/strassen.cpp" -o CMakeFiles/tarea2.dir/strassen.cpp.s

# Object files for target tarea2
tarea2_OBJECTS = \
"CMakeFiles/tarea2.dir/strassen.cpp.o"

# External object files for target tarea2
tarea2_EXTERNAL_OBJECTS =

tarea2: CMakeFiles/tarea2.dir/strassen.cpp.o
tarea2: CMakeFiles/tarea2.dir/build.make
tarea2: CMakeFiles/tarea2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tarea2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tarea2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tarea2.dir/build: tarea2

.PHONY : CMakeFiles/tarea2.dir/build

CMakeFiles/tarea2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tarea2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tarea2.dir/clean

CMakeFiles/tarea2.dir/depend:
	cd "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea2/cmake-build-debug/CMakeFiles/tarea2.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/tarea2.dir/depend

