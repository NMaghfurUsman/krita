set(kritacolorgenerator_SOURCES
    colorgenerator.cpp
    kis_wdg_color.cpp
    )


ki18n_wrap_ui(kritacolorgenerator_SOURCES
    wdgcoloroptions.ui
    )

kis_add_library(kritacolorgenerator MODULE ${kritacolorgenerator_SOURCES})

target_link_libraries(kritacolorgenerator kritaui)

install(TARGETS kritacolorgenerator  DESTINATION ${KRITA_PLUGIN_INSTALL_DIR})
