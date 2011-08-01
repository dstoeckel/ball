### list all filenames of the directory here ###
SET(SOURCES_LIST
resourceFileFragmentStorage.C
nameFragmentQuery.C
propertyFragmentQuery.C
nameMapQuery.C
)

ADD_BALL_SOURCES("STRUCTURE/FRAGMENTDB" "${SOURCES_LIST}")
