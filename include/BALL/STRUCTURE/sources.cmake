### the directory name ###
SET(GROUP STRUCTURE)

FILE(GLOB HEADERS_LIST "include/BALL/${GROUP}/*.h" "include/BALL/${GROUP}/*.iC")	

ADD_BALL_HEADERS("${GROUP}" "${HEADERS_LIST}")

INCLUDE(include/BALL/STRUCTURE/VFLIB/sources.cmake)
INCLUDE(include/BALL/STRUCTURE/FRAGMENTDB/sources.cmake)
