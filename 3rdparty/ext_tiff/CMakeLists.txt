SET(EXTPREFIX_tiff "${EXTPREFIX}" )
ExternalProject_Add( ext_tiff
    DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
    URL https://download.osgeo.org/libtiff/tiff-4.4.0.tar.gz
    URL_HASH SHA256=917223b37538959aca3b790d2d73aa6e626b688e02dcda272aec24c2f498abed

    INSTALL_DIR ${EXTPREFIX_tiff}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTPREFIX_tiff} -DCMAKE_BUILD_TYPE=${GLOBAL_BUILD_TYPE} ${GLOBAL_PROFILE}

    UPDATE_COMMAND ""
)
