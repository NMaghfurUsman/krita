set(KRITA_SPECIFICCOLORSELECTOR_SOURCES 
    specificcolorselector.cc 
    specificcolorselector_dock.cc 
    kis_specific_color_selector_widget.cc 
)

ki18n_wrap_ui(KRITA_SPECIFICCOLORSELECTOR_SOURCES
    wdgSpecificColorSelectorWidget.ui
)


kis_add_library(kritaspecificcolorselector MODULE ${KRITA_SPECIFICCOLORSELECTOR_SOURCES})
target_link_libraries(kritaspecificcolorselector kritaui ${LINK_OPENEXR_LIB})
install(TARGETS kritaspecificcolorselector  DESTINATION ${KRITA_PLUGIN_INSTALL_DIR})
