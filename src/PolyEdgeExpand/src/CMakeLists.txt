# Build and install the library
set( OPS_NAME ${PRJ_NAME}Ops)
add_library( ${OPS_NAME} SHARED)
set( source_files
    "PolyEdgeExpand_nodedef.cpp"
    "PolyEdgeExpand_nodedef.h"
)
target_compile_definitions( ${OPS_NAME} PRIVATE TKCM_BUILD_NODEDEF_DLL )
target_sources( ${OPS_NAME} PRIVATE ${source_files})
target_link_libraries( ${OPS_NAME} PUBLIC Amino::Cpp )
target_link_libraries( ${OPS_NAME} PUBLIC Amino::Core )
file(GLOB MY_SOURCE ${BIFROST_LOCATION}/sdk/lib/*.lib)
target_link_libraries( ${OPS_NAME} PUBLIC ${MY_SOURCE} )
bifrost_set_install_rpath(${OPS_NAME})
install( TARGETS ${OPS_NAME} DESTINATION lib )

#Generate and install the json files
set( headers_to_parse 
    PolyEdgeExpand_nodedef.h
)

amino_cpp2json(
    TARGET          ${PRJ_NAME}JsonFiles
    CPP2JSON        ${BIFROST_LOCATION}/bin/cpp2json.exe
    INPUT           ${headers_to_parse}
    INCLUDES        ${CMAKE_CURRENT_SOURCE_DIR}
    INCLUDES        ${BIFROST_LOCATION}/sdk/include
    INCLUDES        ${BIFROST_LOCATION}/sdk/include/Amino
    DEFINES         "_MT"
    PATH            ${BIFROST_LOCATION}/thirdparty/bin/tbb.dll                       
)

# Create config json file to load the bifrost pack
set( CONFIG_PACK_VERSION    "${PROJECT_VERSION}")
set( CONFIG_PACK_NAME       "${PROJECT_NAME}")
set( CONFIG_SHARED_LIBS     "./lib")
set( CONFIG_JSON_FILES      "./json")
set( config_json            ${CMAKE_BINARY_DIR}/${PROJECT_NAME}PackConfig.json )
set( common_config_json     "../../core/PackConfig.json")
configure_file(${common_config_json}.in ${config_json} )
install(FILES ${config_json} DESTINATION . )

# duplicate dll/json to pack dir 
add_custom_command(
    TARGET ${OPS_NAME} POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE_DIR:${OPS_NAME}>/${OPS_NAME}.dll" "${CMAKE_CURRENT_SOURCE_DIR}/../../../pack/bin"
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE_DIR:${OPS_NAME}>/${PRJ_NAME}_nodedef.json" "${CMAKE_CURRENT_SOURCE_DIR}/../../../pack/compounds/${PRJ_NAME}"
)