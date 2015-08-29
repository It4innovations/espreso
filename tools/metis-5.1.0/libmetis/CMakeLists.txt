# Add this directory for internal users.
include_directories(.)
# Find sources.
file(GLOB metis_sources *.c)
# Build libmetis.
add_library(metis${IDXTYPEWIDTH} ${METIS_LIBRARY_TYPE} ${GKlib_sources} ${metis_sources})
if(UNIX)
  target_link_libraries(metis${IDXTYPEWIDTH} m)
endif()

if(METIS_INSTALL)
  install(TARGETS metis${IDXTYPEWIDTH}
    LIBRARY DESTINATION libs
    RUNTIME DESTINATION libs
    ARCHIVE DESTINATION libs)
endif()
