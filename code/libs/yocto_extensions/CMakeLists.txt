add_library(yocto_extension 
    yocto_extension_hair.h yocto_extension_hair.cpp
    yocto_extension_montecarlo.h yocto_extension_montecarlo.cpp)

set_target_properties(yocto_extension PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yocto_extension PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yocto_extension yocto yocto_pathtrace)
