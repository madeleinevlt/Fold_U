# Install script for directory: /home/mikov/Documents/hh-suite/scripts

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts" TYPE PROGRAM PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ FILES
    "/home/mikov/Documents/hh-suite/scripts/Align.pm"
    "/home/mikov/Documents/hh-suite/scripts/HHPaths.pm"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts" TYPE PROGRAM FILES
    "/home/mikov/Documents/hh-suite/scripts/addss.pl"
    "/home/mikov/Documents/hh-suite/scripts/create_profile_from_hhm.pl"
    "/home/mikov/Documents/hh-suite/scripts/hhmakemodel.pl"
    "/home/mikov/Documents/hh-suite/scripts/hhsuitedb.py"
    "/home/mikov/Documents/hh-suite/scripts/mergeali.pl"
    "/home/mikov/Documents/hh-suite/scripts/multithread.pl"
    "/home/mikov/Documents/hh-suite/scripts/pdbfilter.pl"
    "/home/mikov/Documents/hh-suite/scripts/renumberpdb.pl"
    "/home/mikov/Documents/hh-suite/scripts/create_profile_from_hmmer.pl"
    "/home/mikov/Documents/hh-suite/scripts/pdb2fasta.pl"
    "/home/mikov/Documents/hh-suite/scripts/reformat.pl"
    "/home/mikov/Documents/hh-suite/scripts/splitfasta.pl"
    "/home/mikov/Documents/hh-suite/scripts/check_a3m.py"
    "/home/mikov/Documents/hh-suite/scripts/ffindex.py"
    "/home/mikov/Documents/hh-suite/scripts/a3m.py"
    "/home/mikov/Documents/hh-suite/scripts/get_a3m_size.py"
    "/home/mikov/Documents/hh-suite/scripts/pdbfilter.py"
    "/home/mikov/Documents/hh-suite/scripts/cif2fasta.py"
    "/home/mikov/Documents/hh-suite/scripts/hhmakemodel.py"
    "/home/mikov/Documents/hh-suite/scripts/hh_reader.py"
    )
endif()

