set(KRITA_PALETTEDOCKER_SOURCES
    palettedocker.cpp palettedocker_dock.cpp
    )

ki18n_wrap_ui(KRITA_PALETTEDOCKER_SOURCES
    wdgpalettedock.ui
)

kis_add_library(kritapalettedocker MODULE ${KRITA_PALETTEDOCKER_SOURCES})
target_link_libraries(kritapalettedocker kritaui)
install(TARGETS kritapalettedocker  DESTINATION ${KRITA_PLUGIN_INSTALL_DIR})
