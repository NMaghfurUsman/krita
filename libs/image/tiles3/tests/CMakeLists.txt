
include(KritaAddBrokenUnitTest)

kis_add_tests(
    kis_tiled_data_manager_test.cpp
    kis_low_memory_tests.cpp
    kis_lockless_stack_test.cpp
    kis_chunk_allocator_test.cpp
    kis_memory_window_test.cpp
    kis_store_limits_test.cpp
    kis_swapped_data_store_test.cpp
    kis_tile_data_store_test.cpp
    kis_tile_data_pooler_test.cpp
    LINK_LIBRARIES kritaimage kritatestsdk
    NAME_PREFIX "libs-image-tiles3-"
    )

set_tests_properties(libs-image-tiles3-kis_low_memory_tests PROPERTIES TIMEOUT 180)
