set(STREAMING_ENGINE_FILE_LIST
        global_view_manager.cpp
        local_view_manager.cpp
        order_manager.cpp
        search_engine.cpp
        streaming_engine.cpp
        )

add_executable(HybridFlow.out ${STREAMING_ENGINE_FILE_LIST})
target_link_libraries(HybridFlow.out graph utility_core)