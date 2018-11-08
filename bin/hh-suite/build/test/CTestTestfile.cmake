# CMake generated Testfile for 
# Source directory: /home/mikov/Documents/hh-suite/test
# Build directory: /home/mikov/Documents/hh-suite/build/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(hhblits "/home/mikov/Documents/hh-suite/build/bin/hhblits" "-h")
set_tests_properties(hhblits PROPERTIES  ENVIRONMENT "HHLIB=/home/mikov/Documents/hh-suite")
add_test(hhsearch "/home/mikov/Documents/hh-suite/build/bin/hhsearch" "-h")
set_tests_properties(hhsearch PROPERTIES  ENVIRONMENT "HHLIB=/home/mikov/Documents/hh-suite")
add_test(hhalign "/home/mikov/Documents/hh-suite/build/bin/hhalign" "-h")
set_tests_properties(hhalign PROPERTIES  ENVIRONMENT "HHLIB=/home/mikov/Documents/hh-suite")
