project (StarNEig)
cmake_minimum_required (VERSION 3.3)
enable_language (Fortran)

include (CheckFortranCompilerFlag)

function (enable_fortran_compiler_flag_if_supported flag)
    string(REPLACE "=" "_" flag_name ${flag})
    string (FIND "${CMAKE_Fortran_FLAGS}" "${flag}" flag_already_set)
    if (flag_already_set EQUAL -1)
        check_fortran_compiler_flag (
            "${flag}" Fortran_FLAG_DETECTED${flag_name})
        if (Fortran_FLAG_DETECTED${flag_name})
            set (CMAKE_Fortran_FLAGS
                "${flag} ${CMAKE_Fortran_FLAGS}" PARENT_SCOPE)
        endif ()
    endif ()
endfunction ()

enable_fortran_compiler_flag_if_supported ("-Wall")
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    enable_fortran_compiler_flag_if_supported ("-Wno-maybe-uninitialized")
    enable_fortran_compiler_flag_if_supported ("-Wno-unused-label")
    enable_fortran_compiler_flag_if_supported ("-Wno-unused-variable")
    enable_fortran_compiler_flag_if_supported ("-Wno-unused-dummy-argument")
    enable_fortran_compiler_flag_if_supported ("-Wno-character-truncation")
    enable_fortran_compiler_flag_if_supported ("-Wno-conversion")
    enable_fortran_compiler_flag_if_supported ("-Wno-tabs")

    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        enable_fortran_compiler_flag_if_supported ("-march=native")
        enable_fortran_compiler_flag_if_supported ("-O2")
    endif ()
endif ()

file (GLOB_RECURSE SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]")

add_library (pdtgsen SHARED ${SOURCES})
