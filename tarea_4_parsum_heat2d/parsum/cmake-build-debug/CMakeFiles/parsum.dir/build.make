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
CMAKE_SOURCE_DIR = "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/parsum.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/parsum.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/parsum.dir/flags.make

CMakeFiles/parsum.dir/parsum.c.o: CMakeFiles/parsum.dir/flags.make
CMakeFiles/parsum.dir/parsum.c.o: ../parsum.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/parsum.dir/parsum.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/parsum.dir/parsum.c.o   -c "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/parsum.c"

CMakeFiles/parsum.dir/parsum.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/parsum.dir/parsum.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/parsum.c" > CMakeFiles/parsum.dir/parsum.c.i

CMakeFiles/parsum.dir/parsum.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/parsum.dir/parsum.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/parsum.c" -o CMakeFiles/parsum.dir/parsum.c.s

# Object files for target parsum
parsum_OBJECTS = \
"CMakeFiles/parsum.dir/parsum.c.o"

# External object files for target parsum
parsum_EXTERNAL_OBJECTS =

parsum: CMakeFiles/parsum.dir/parsum.c.o
parsum: CMakeFiles/parsum.dir/build.make
parsum: /usr/lib/x86_64-linux-gnu/libmpich.so
parsum: CMakeFiles/parsum.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable parsum"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/parsum.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/parsum.dir/build: parsum

.PHONY : CMakeFiles/parsum.dir/build

CMakeFiles/parsum.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/parsum.dir/cmake_clean.cmake
.PHONY : CMakeFiles/parsum.dir/clean

CMakeFiles/parsum.dir/depend:
	cd "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/tarea_4_parsum_heat2d/parsum/cmake-build-debug/CMakeFiles/parsum.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/parsum.dir/depend

