include(KritaAddBrokenUnitTest)

kis_add_test(
    KisPdfTest.cpp
    TEST_NAME KisPdfTest
    LINK_LIBRARIES kritaui kritatestsdk
    NAME_PREFIX "plugins-impex-"
    )
