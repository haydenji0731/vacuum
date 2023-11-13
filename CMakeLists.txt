cmake_minimum_required(VERSION 3.18)
project(vacuum)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fpermissive -DNOCURL=1")
include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")
include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib
        # if this is not specified, need to include the packages manually
        GIT_REPOSITORY https://github.com/samtools/htslib.git
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --without-libdeflate --disable-libcurl --disable-lzma
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(vacuum
        src/vacuum.cpp
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GHashMap.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GVec.hh)

add_executable(junction_extractor 
	src/junction_extractor.cpp 
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.h

        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.h

        ${CMAKE_SOURCE_DIR}/include/gclib/codons.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.h

        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.h

        ${CMAKE_SOURCE_DIR}/include/gclib/GList.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GHashMap.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GBitVec.h
        
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
        )

add_dependencies(vacuum htslib)
target_link_libraries(vacuum ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z bz2 pthread)
install(TARGETS vacuum DESTINATION ${CMAKE_INSTALL_PREFIX})

add_dependencies(junction_extractor htslib)
target_link_libraries(junction_extractor ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z bz2 pthread)
install(TARGETS junction_extractor DESTINATION ${CMAKE_INSTALL_PREFIX})
