SET(PREFIX_ext_freetype "${EXTPREFIX}" )
SET(EXTPREFIX_harfbuzz "${EXTPREFIX}" )

if (POLICY CMP0114)
    cmake_policy(SET CMP0114 NEW)
endif()

set(FREETYPE_URL https://downloads.sourceforge.net/project/freetype/freetype2/2.12.1/freetype-2.12.1.tar.xz)
set(FREETYPE_SHA256 4766f20157cc4cf0cd292f80bf917f92d1c439b243ac3018debf6b9140c41a7f)

ExternalProject_Add(ext_freetype_bootstrap
    DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}

    URL ${FREETYPE_URL}
    URL_HASH SHA256=${FREETYPE_SHA256}

    INSTALL_DIR ${PREFIX_ext_freetype}

    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PREFIX_ext_freetype}
        -DFT_DISABLE_BZIP2=ON
        -DFT_DISABLE_BROTLI=ON
        -DFT_REQUIRE_PNG=ON
        -DFT_DISABLE_HARFBUZZ=ON
        -DCMAKE_BUILD_TYPE=${GLOBAL_BUILD_TYPE}
        -DBUILD_SHARED_LIBS=TRUE
        ${GLOBAL_PROFILE}

    UPDATE_COMMAND ""

    DEPENDS ext_png
)

ExternalProject_Add(
    ext_harfbuzz
    DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
    URL https://github.com/harfbuzz/harfbuzz/releases/download/5.3.1/harfbuzz-5.3.1.tar.xz
    URL_HASH SHA256=4a6ce097b75a8121facc4ba83b5b083bfec657f45b003cd5a3424f2ae6b4434d

    PATCH_COMMAND ${PATCH_COMMAND} -p1 -i ${CMAKE_CURRENT_SOURCE_DIR}/0001-patch-msvc.patch

    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH=${_krita_pythonpath}
        PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1
        ${MESON_BINARY_PATH} setup <BINARY_DIR> <SOURCE_DIR>
            --prefix=${EXTPREFIX_harfbuzz}
            --libdir=lib
            -Dbuildtype=$<IF:$<CONFIG:Debug>,debug,debugoptimized>
            -Dglib=disabled
            -Dgobject=disabled
            -Dtests=disabled
            ${EXTRA_MESON_FLAGS}

    BUILD_COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH=${_krita_pythonpath}
        ${MESON_BINARY_PATH} compile -C <BINARY_DIR> -j${SUBMAKE_JOBS}

    INSTALL_COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH=${_krita_pythonpath}
        ${MESON_BINARY_PATH} install -C <BINARY_DIR>

    UPDATE_COMMAND ""

    DEPENDS ext_freetype_bootstrap ${MESON_DEP} ${PKG_CONFIG_DEP}
)

# Amyspark: just run roughshod!
ExternalProject_Add( ext_freetype
    DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
    URL ${FREETYPE_URL}
    URL_HASH SHA256=${FREETYPE_SHA256}

    INSTALL_DIR ${PREFIX_ext_freetype}

    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DCMAKE_BUILD_TYPE=${GLOBAL_BUILD_TYPE}
        -DFT_DISABLE_BZIP2=ON
        -DFT_DISABLE_BROTLI=ON
        -DFT_REQUIRE_HARFBUZZ=ON
        -DFT_REQUIRE_PNG=ON
        -DBUILD_SHARED_LIBS=TRUE
        ${GLOBAL_PROFILE}

    UPDATE_COMMAND ""
    DEPENDS ext_harfbuzz ext_png
)
