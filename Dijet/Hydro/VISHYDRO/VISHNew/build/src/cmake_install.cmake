# Install script for directory: /Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/VISHNew.e")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code" TYPE EXECUTABLE FILES "/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/build/src/VISHNew.e")
  if(EXISTS "$ENV{DESTDIR}/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/VISHNew.e" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/VISHNew.e")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/VISHNew.e")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/Users/duanpanpan/Desktop/Work/Model/Hydro/VISHYDRO/VISHNew/VISHNew_Code/build/src/CMakeFiles/VISHNew.e.dir/install-cxx-module-bmi-RELEASE.cmake" OPTIONAL)
endif()

