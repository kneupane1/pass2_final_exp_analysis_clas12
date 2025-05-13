cmake_minimum_required(VERSION 3.5)
project(pass2_analysis_clas12)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

include_directories(src/include)
# # Use the QADB environment variable for absolute paths
 include_directories($ENV{QADB}/srcC/include)
 include_directories($ENV{QADB}/srcC/rapidjson/include)
 include_directories(/Users/krishnaneupane/Documents/GitHub/clas12-qadb/srcC/include)
 include_directories(/Users/krishnaneupane/Documents/GitHub/clas12-qadb/srcC/rapidjson/include)

# # Print include directories for debugging
 message(STATUS "Include directories: ${CMAKE_INCLUDE_PATH}")
 message(STATUS "QADB include directory: $ENV{QADB}/srcC/include")
 message(STATUS "RapidJSON include directory: $ENV{QADB}/srcC/rapidjson/include")

file(GLOB clas12lib_SRC
    "src/lib/*.cpp"
)
add_library(clas12lib STATIC ${clas12lib_SRC})


file(GLOB clas12_EXE
    "exe/*.cpp"
)

foreach(exe_file IN LISTS clas12_EXE)
    get_filename_component(exe_name ${exe_file} NAME_WE)
    add_executable(${exe_name} ${exe_file})
    target_link_libraries(${exe_name} clas12lib ${ROOT_LIBRARIES})
endforeach()
