SET(PREFIX_ext_exiv2  "${EXTPREFIX}" )

# Note: We explicitly set CMAKE_CXX_STANDARD=14 when compiling exiv2 because
# it uses `auto_ptr`, which has been removed in c++17. Clang 16 switches to
# using -std=c++17 by default, and that causes build failures with libc++.

# FIXME: Remove this setting when exiv2 has been updated to remove usage of
# `auto_ptr`, and also remove the workaround in `FindLibExiv2.cmake`.

if (ANDROID)
    ExternalProject_Add( ext_exiv2
        DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
        URL https://github.com/Exiv2/exiv2/releases/download/v0.27.5/exiv2-0.27.5-Source.tar.gz
        URL_HASH SHA256=35a58618ab236a901ca4928b0ad8b31007ebdc0386d904409d825024e45ea6e2

        INSTALL_DIR ${PREFIX_ext_exiv2}
        CMAKE_ARGS -DCMAKE_CXX_STANDARD=14 -DCMAKE_INSTALL_PREFIX=${PREFIX_ext_exiv2} -DCMAKE_BUILD_TYPE=${GLOBAL_BUILD_TYPE} ${GLOBAL_PROFILE} -DEXIV2_BUILD_SAMPLES=OFF -DICONV_LIBRARY=${PREFIX_ext_exiv2}/lib -DEXIV2_ENABLE_NLS=OFF -DICONV_INCLUDE_DIR=${PREFIX_ext_exiv2}/include -DEXPAT_LIBRARY=$ENV{BUILD_ROOT}/i/lib/libexpat.so -DEXPAT_INCLUDE_DIR=$ENV{BUILD_ROOT}/i/include -DEXIV2_BUILD_EXIV2_COMMAND=OFF
        UPDATE_COMMAND ""
        DEPENDS ext_expat
    )

elseif(WIN32)

    set(WINDOWS_PATCH_COMMAND ${PATCH_COMMAND} -p1 -i ${CMAKE_CURRENT_SOURCE_DIR}/0001-Fix-building-exiv2-with-TDMGCC-10.3.0.patch)
                
    if (ICONV_DEP)
        # our internal ext_iconv has a bit different API
    
        set(WINDOWS_PATCH_COMMAND 
            ${WINDOWS_PATCH_COMMAND}
            COMMAND ${PATCH_COMMAND} -p1 -i ${CMAKE_CURRENT_SOURCE_DIR}/01_fpermissive_error.patch)
    endif()
    

    ExternalProject_Add( ext_exiv2
        DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
        URL https://github.com/Exiv2/exiv2/releases/download/v0.27.5/exiv2-0.27.5-Source.tar.gz
        URL_HASH SHA256=35a58618ab236a901ca4928b0ad8b31007ebdc0386d904409d825024e45ea6e2

        PATCH_COMMAND ${WINDOWS_PATCH_COMMAND}

        INSTALL_DIR ${PREFIX_ext_exiv2}
        CMAKE_ARGS -DCMAKE_CXX_STANDARD=14 -DCMAKE_INSTALL_PREFIX=${PREFIX_ext_exiv2} -DCMAKE_LIBRARY_PATH=${PREFIX_ext_exiv2}/lib -DCMAKE_BUILD_TYPE=${GLOBAL_BUILD_TYPE} ${GLOBAL_PROFILE} -DEXIV2_BUILD_SAMPLES=OFF -DIconv_INCLUDE_DIR=${PREFIX_ext_exiv2}/include -DEXIV2_ENABLE_NLS=OFF -DEXIV2_BUILD_EXIV2_COMMAND=OFF

        UPDATE_COMMAND ""
        DEPENDS ${ICONV_DEP} ${EXPAT_DEP} ${ZLIB_DEP}
    )
else()
    ExternalProject_Add( ext_exiv2

        DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
        URL https://github.com/Exiv2/exiv2/releases/download/v0.27.5/exiv2-0.27.5-Source.tar.gz
        URL_HASH SHA256=35a58618ab236a901ca4928b0ad8b31007ebdc0386d904409d825024e45ea6e2

        INSTALL_DIR ${PREFIX_ext_exiv2}
        CMAKE_ARGS -DCMAKE_CXX_STANDARD=14 -DCMAKE_INSTALL_PREFIX=${PREFIX_ext_exiv2} -DCMAKE_LIBRARY_PATH=${PREFIX_ext_exiv2}/lib -DCMAKE_BUILD_TYPE=${GLOBAL_BUILD_TYPE} ${GLOBAL_PROFILE} -DEXIV2_BUILD_SAMPLES=OFF -DIconv_INCLUDE_DIR=${PREFIX_ext_exiv2}/include -DEXIV2_ENABLE_NLS=OFF -DEXIV2_BUILD_EXIV2_COMMAND=OFF

        UPDATE_COMMAND ""
        DEPENDS ${ICONV_DEP} ext_expat
    )
    
endif()
