# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /home/linuxbrew/.linuxbrew/Cellar/cmake/3.20.1/bin/cmake

# The command to remove a file.
RM = /home/linuxbrew/.linuxbrew/Cellar/cmake/3.20.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cuongng/enzyme2/enzyme

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cuongng/enzyme2/enzyme

# Utility rule file for check-activityanalysis.

# Include any custom commands dependencies for this target.
include test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/compiler_depend.make

# Include the progress variables for this target.
include test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/progress.make

test/ActivityAnalysis/CMakeFiles/check-activityanalysis:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/cuongng/enzyme2/enzyme/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Running enzyme regression tests"
	cd /home/cuongng/enzyme2/enzyme/test/ActivityAnalysis && ../../ -v /home/cuongng/enzyme2/enzyme/test/ActivityAnalysis

check-activityanalysis: test/ActivityAnalysis/CMakeFiles/check-activityanalysis
check-activityanalysis: test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/build.make
.PHONY : check-activityanalysis

# Rule to build all files generated by this target.
test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/build: check-activityanalysis
.PHONY : test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/build

test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/clean:
	cd /home/cuongng/enzyme2/enzyme/test/ActivityAnalysis && $(CMAKE_COMMAND) -P CMakeFiles/check-activityanalysis.dir/cmake_clean.cmake
.PHONY : test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/clean

test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/depend:
	cd /home/cuongng/enzyme2/enzyme && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cuongng/enzyme2/enzyme /home/cuongng/enzyme2/enzyme/test/ActivityAnalysis /home/cuongng/enzyme2/enzyme /home/cuongng/enzyme2/enzyme/test/ActivityAnalysis /home/cuongng/enzyme2/enzyme/test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/ActivityAnalysis/CMakeFiles/check-activityanalysis.dir/depend

