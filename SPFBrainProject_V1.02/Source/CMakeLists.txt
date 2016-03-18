cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)
project(LesionSegmentation)

find_package(ITK REQUIRED)
find_package(VTK REQUIRED)

include_directories(HeaderFiles)
include_directories(HeaderFiles/Filters)

set( SPFIO src/Filters/LesionSegmentationCLI.cxx )
set( SPFFilter src/Filters/IsotropicResamplerImageFilter.cxx  src/Filters/SPFLevelSetImageFilter.cxx )

add_library( SPFBrain ${SPFIO} ${SPFFilter} )
#target_link_libraries(SPFBrain ${ITK_LIBRARIES})
#link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)


include(${ITK_USE_FILE})
include(${VTK_USE_FILE})

add_executable(LesionSegmentation src/main.cxx )
target_link_libraries(LesionSegmentation  ${ITK_LIBRARIES} ${VTK_LIBRARIES}  SPFBrain )
