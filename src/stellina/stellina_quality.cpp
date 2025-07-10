/*
 * stellina_quality.cpp
 * Quality assessment functions for Stellina images
 */

#include "stellina_processor.h"
#include <iostream>
#include <fstream>
#include <cstring>

/**
 * Parse Stellina JSON metadata file
 */
struct stellina_metadata *stellina_parse_json(const char *json_path) {
    if (!json_path || !g_file_test(json_path, G_FILE_TEST_EXISTS)) {
        std::cerr << "JSON file not found: " << json_path << std::endl;
        return nullptr;
    }
    
    struct stellina_metadata *metadata = g_new0(struct stellina_metadata, 1);
    
    // Initialize with defaults
    metadata->quality_accepted = TRUE;
    g_strlcpy(metadata->quality_reason, "No quality info found - assuming accepted", 
              sizeof(metadata->quality_reason));
    metadata->fwhm = 0.0;
    metadata->star_count = 0;
    metadata->used_for_stacking = TRUE;
    g_strlcpy(metadata->original_filename, g_path_get_basename(json_path), 
              sizeof(metadata->original_filename));
    
    JsonParser *parser = json_parser_new();
    GError *error = nullptr;
    
    if (!json_parser_load_from_file(parser, json_path, &error)) {
        std::cerr << "Failed to parse JSON file " << json_path << ": " << error->message << std::endl;
        g_error_free(error);
        g_object_unref(parser);
        stellina_metadata_free(metadata);
        return nullptr;
    }
    
    JsonNode *root = json_parser_get_root(parser);
    if (!JSON_NODE_HOLDS_OBJECT(root)) {
        std::cerr << "Invalid JSON format in " << json_path << std::endl;
        g_object_unref(parser);
        stellina_metadata_free(metadata);
        return nullptr;
    }
    
    JsonObject *root_object = json_node_get_object(root);
    
    // Extract altitude and azimuth
    if (json_object_has_member(root_object, "altitude")) {
        metadata->altitude = json_object_get_double_member(root_object, "altitude");
    } else if (json_object_has_member(root_object, "alt")) {
        metadata->altitude = json_object_get_double_member(root_object, "alt");
    }
    
    if (json_object_has_member(root_object, "azimuth")) {
        metadata->azimuth = json_object_get_double_member(root_object, "azimuth");
    } else if (json_object_has_member(root_object, "az")) {
        metadata->azimuth = json_object_get_double_member(root_object, "az");
    }
    
    // Extract observation time
    if (json_object_has_member(root_object, "date_obs")) {
        const char *date_str = json_object_get_string_member(root_object, "date_obs");
        g_strlcpy(metadata->date_obs, date_str, sizeof(metadata->date_obs));
    } else if (json_object_has_member(root_object, "time")) {
        const char *date_str = json_object_get_string_member(root_object, "time");
        g_strlcpy(metadata->date_obs, date_str, sizeof(metadata->date_obs));
    }
    
    // Check quality indicators (same logic as before)
    gboolean quality_found = FALSE;
    
    if (json_object_has_member(root_object, "quality")) {
        metadata->quality_accepted = json_object_get_boolean_member(root_object, "quality");
        g_strlcpy(metadata->quality_reason, 
                  metadata->quality_accepted ? "Accepted by Stellina" : "Rejected by Stellina",
                  sizeof(metadata->quality_reason));
        quality_found = TRUE;
    }
    
    // ... (rest of the quality checking logic stays the same)
    
    g_object_unref(parser);
    
    std::cout << "Parsed JSON: Alt=" << metadata->altitude << "°, Az=" << metadata->azimuth 
              << "°, Quality=" << (metadata->quality_accepted ? "ACCEPTED" : "REJECTED") << std::endl;
    
    return metadata;
}

/**
 * Check if image meets quality criteria
 */
gboolean stellina_check_quality(const struct stellina_metadata *metadata) {
    if (!metadata) {
        return FALSE;
    }
    
    if (!metadata->quality_accepted) {
        std::cout << "Quality check FAILED: " << metadata->quality_reason << std::endl;
        return FALSE;
    }
    
    std::cout << "Quality check PASSED: " << metadata->quality_reason << std::endl;
    return TRUE;
}

/**
 * Free stellina metadata structure
 */
void stellina_metadata_free(struct stellina_metadata *metadata) {
    if (metadata) {
        g_free(metadata);
    }
}