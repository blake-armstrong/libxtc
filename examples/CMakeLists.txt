# Examples CMakeLists.txt
cmake_minimum_required(VERSION 3.14)

# Basic usage example
add_executable(xtc_example example_usage.cpp)
target_link_libraries(xtc_example PRIVATE xtc::xtc)

# Install examples (optional)
install(TARGETS 
    xtc_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
