include(KritaAddBrokenUnitTest)

kis_add_test(KisJP2Test.cpp
    TEST_NAME KisJP2Test
    LINK_LIBRARIES kritaui kritatestsdk
    NAME_PREFIX "plugins-impex-"
)

