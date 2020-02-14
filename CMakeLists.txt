cmake_minimum_required(VERSION 2.8)
project(swamp NONE)


IF(DEFINED ENV{CCP4})
    SET(ccp4_py2   $ENV{CCP4}/lib/py2)
    SET(ccp4_bin   $ENV{CCP4}/bin)
    SET(ccp4_share   $ENV{CCP4}/share)
ELSE()
    MESSAGE(FATAL_ERROR "CCP4 must be installed before using SWAMP, cannot install")
ENDIF()

install (
    DIRECTORY
        swamp
    DESTINATION
        ${ccp4_py2}
    )

install (
    PROGRAMS
        bin/swamp-mr bin/swamp-make-library bin/swamp-results
    DESTINATION
        ${ccp4_bin}
    )

install (
    DIRECTORY
        static
    DESTINATION
        ${ccp4_share}/swamp
    )
