# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /snap/clion/123/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/123/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/lead_pdc_project.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lead_pdc_project.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lead_pdc_project.dir/flags.make

CMakeFiles/lead_pdc_project.dir/cenatMD.c.o: CMakeFiles/lead_pdc_project.dir/flags.make
CMakeFiles/lead_pdc_project.dir/cenatMD.c.o: ../cenatMD.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/lead_pdc_project.dir/cenatMD.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/lead_pdc_project.dir/cenatMD.c.o   -c "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cenatMD.c"

CMakeFiles/lead_pdc_project.dir/cenatMD.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/lead_pdc_project.dir/cenatMD.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cenatMD.c" > CMakeFiles/lead_pdc_project.dir/cenatMD.c.i

CMakeFiles/lead_pdc_project.dir/cenatMD.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/lead_pdc_project.dir/cenatMD.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cenatMD.c" -o CMakeFiles/lead_pdc_project.dir/cenatMD.c.s

# Object files for target lead_pdc_project
lead_pdc_project_OBJECTS = \
"CMakeFiles/lead_pdc_project.dir/cenatMD.c.o"

# External object files for target lead_pdc_project
lead_pdc_project_EXTERNAL_OBJECTS =

lead_pdc_project: CMakeFiles/lead_pdc_project.dir/cenatMD.c.o
lead_pdc_project: CMakeFiles/lead_pdc_project.dir/build.make
lead_pdc_project: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
lead_pdc_project: CMakeFiles/lead_pdc_project.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable lead_pdc_project"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lead_pdc_project.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lead_pdc_project.dir/build: lead_pdc_project

.PHONY : CMakeFiles/lead_pdc_project.dir/build

CMakeFiles/lead_pdc_project.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lead_pdc_project.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lead_pdc_project.dir/clean

CMakeFiles/lead_pdc_project.dir/depend:
	cd "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug" "/mnt/napster_disk/LEAD University/IIQ - 2020/Computación Paralela y Distribuída/tareas/proyecto/lead_pdc_project/cmake-build-debug/CMakeFiles/lead_pdc_project.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/lead_pdc_project.dir/depend

