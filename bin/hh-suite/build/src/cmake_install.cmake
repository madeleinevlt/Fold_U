# Install script for directory: /home/mikov/Documents/hh-suite/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/home/mikov/Documents/hh-suite/build/bin/hhblits_omp"
    "/home/mikov/Documents/hh-suite/build/bin/hhsearch_omp"
    "/home/mikov/Documents/hh-suite/build/bin/hhblits_ca3m"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/home/mikov/Documents/hh-suite/build/bin/hhblits"
    "/home/mikov/Documents/hh-suite/build/bin/hhmake"
    "/home/mikov/Documents/hh-suite/build/bin/hhfilter"
    "/home/mikov/Documents/hh-suite/build/bin/hhsearch"
    "/home/mikov/Documents/hh-suite/build/bin/hhalign"
    "/home/mikov/Documents/hh-suite/build/bin/hhconsensus"
    "/home/mikov/Documents/hh-suite/build/bin/a3m_extract"
    "/home/mikov/Documents/hh-suite/build/bin/a3m_database_reduce"
    "/home/mikov/Documents/hh-suite/build/bin/a3m_database_extract"
    "/home/mikov/Documents/hh-suite/build/bin/a3m_database_filter"
    "/home/mikov/Documents/hh-suite/build/bin/cstranslate"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/mikov/Documents/hh-suite/build/src/cs/cmake_install.cmake")

endif()

