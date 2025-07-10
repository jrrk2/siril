/*
 * stellina_quality.cpp
 * Quality assessment functions for Stellina images
 */

#include "stellina_processor.h"
#include "core/siril.h"
#include "gui/utils.h"
#include <json-glib/json-glib.h>

/**
 * Parse Stellina JSON metadata file
 */
struct stellina_metadata *stellina_parse_json(const char *json_path) {
    if (!json_path || !g_file_test(json_path, G_FILE_TEST_EXISTS)) {
        siril_log_color_message("JSON file not found: %s\n", "red", json_path);
        return NULL;
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
    GError *error = NULL;
    
    if (!json_parser_load_from_file(parser, json_path, &error)) {
        siril_log_color_message("Failed to parse JSON file %s: %s\n", "red", 
                               json_path, error->message);
        g_error_free(error);
        g_object_unref(parser);
        stellina_metadata_free(metadata);
        return NULL;
    }
    
    JsonNode *root = json_parser_get_root(parser);
    if (!JSON_NODE_HOLDS_OBJECT(root)) {
        siril_log_color_message("Invalid JSON format in %s\n", "red", json_path);
        g_object_unref(parser);
        stellina_metadata_free(metadata);
        return NULL;
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
    
    // Check quality indicators
    gboolean quality_found = FALSE;
    
    // Check for explicit quality field
    if (json_object_has_member(root_object, "quality")) {
        metadata->quality_accepted = json_object_get_boolean_member(root_object, "quality");
        g_strlcpy(metadata->quality_reason, 
                  metadata->quality_accepted ? "Accepted by Stellina" : "Rejected by Stellina",
                  sizeof(metadata->quality_reason));
        quality_found = TRUE;
    }
    
    // Check for used_for_stacking field
    if (json_object_has_member(root_object, "used_for_stacking")) {
        metadata->used_for_stacking = json_object_get_boolean_member(root_object, "used_for_stacking");
        metadata->quality_accepted = metadata->used_for_stacking;
        g_strlcpy(metadata->quality_reason,
                  metadata->used_for_stacking ? "Used in Stellina stacking" : "Not used in Stellina stacking",
                  sizeof(metadata->quality_reason));
        quality_found = TRUE;
    }
    
    // Check for accepted field
    if (json_object_has_member(root_object, "accepted")) {
        metadata->quality_accepted = json_object_get_boolean_member(root_object, "accepted");
        g_strlcpy(metadata->quality_reason,
                  metadata->quality_accepted ? "Accepted by Stellina" : "Rejected by Stellina",
                  sizeof(metadata->quality_reason));
        quality_found = TRUE;
    }
    
    // Check stacking object
    if (json_object_has_member(root_object, "stacking")) {
        JsonObject *stacking = json_object_get_object_member(root_object, "stacking");
        if (json_object_has_member(stacking, "used")) {
            metadata->used_for_stacking = json_object_get_boolean_member(stacking, "used");
            metadata->quality_accepted = metadata->used_for_stacking;
            g_strlcpy(metadata->quality_reason,
                      metadata->used_for_stacking ? "Used in Stellina stacking" : "Not used in Stellina stacking",
                      sizeof(metadata->quality_reason));
            quality_found = TRUE;
        }
    }
    
    // Extract additional quality metrics
    if (json_object_has_member(root_object, "fwhm")) {
        metadata->fwhm = json_object_get_double_member(root_object, "fwhm");
    }
    
    if (json_object_has_member(root_object, "star_count")) {
        metadata->star_count = (int)json_object_get_int_member(root_object, "star_count");
    }
    
    // Check stars object
    if (json_object_has_member(root_object, "stars")) {
        JsonObject *stars = json_object_get_object_member(root_object, "stars");
        if (json_object_has_member(stars, "count")) {
            metadata->star_count = (int)json_object_get_int_member(stars, "count");
        }
        if (json_object_has_member(stars, "fwhm")) {
            metadata->fwhm = json_object_get_double_member(stars, "fwhm");
        }
    }
    
    // Check analysis object
    if (json_object_has_member(root_object, "analysis")) {
        JsonObject *analysis = json_object_get_object_member(root_object, "analysis");
        if (json_object_has_member(analysis, "fwhm")) {
            metadata->fwhm = json_object_get_double_member(analysis, "fwhm");
        }
        if (json_object_has_member(analysis, "star_count")) {
            metadata->star_count = (int)json_object_get_int_member(analysis, "star_count");
        }
    }
    
    g_object_unref(parser);
    
    siril_debug_print("Parsed JSON: Alt=%.2f°, Az=%.2f°, Quality=%s, Stars=%d, FWHM=%.2f\n",
                     metadata->altitude, metadata->azimuth, 
                     metadata->quality_accepted ? "ACCEPTED" : "REJECTED",
                     metadata->star_count, metadata->fwhm);
    
    return metadata;
}

/**
 * Check if image meets quality criteria
 */
gboolean stellina_check_quality(const struct stellina_metadata *metadata) {
    if (!metadata) {
        return FALSE;
    }
    
    // Primary quality check based on Stellina's assessment
    if (!metadata->quality_accepted) {
        siril_debug_print("Quality check FAILED: %s\n", metadata->quality_reason);
        return FALSE;
    }
    
    // Additional quality checks can be added here
    // For example, minimum star count, FWHM thresholds, etc.
    
    if (metadata->star_count > 0 && metadata->star_count < 10) {
        siril_log_message("Warning: Low star count (%d) in %s\n", 
                         metadata->star_count, metadata->original_filename);
    }
    
    if (metadata->fwhm > 0.0 && metadata->fwhm > 5.0) {
        siril_log_message("Warning: High FWHM (%.2f) in %s\n", 
                         metadata->fwhm, metadata->original_filename);
    }
    
    siril_debug_print("Quality check PASSED: %s\n", metadata->quality_reason);
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

/**
 * Print quality summary for an image
 */
void stellina_print_quality_summary(const struct stellina_metadata *metadata) {
    if (!metadata) {
        return;
    }
    
    siril_log_message("Quality Summary for %s:\n", metadata->original_filename);
    siril_log_message("  Status: %s\n", metadata->quality_accepted ? "ACCEPTED" : "REJECTED");
    siril_log_message("  Reason: %s\n", metadata->quality_reason);
    
    if (metadata->star_count > 0) {
        siril_log_message("  Star count: %d\n", metadata->star_count);
    }
    
    if (metadata->fwhm > 0.0) {
        siril_log_message("  FWHM: %.2f pixels\n", metadata->fwhm);
    }
    
    siril_log_message("  Used for stacking: %s\n", metadata->used_for_stacking ? "Yes" : "No");
    siril_log_message("  Coordinates: Alt=%.2f°, Az=%.2f°\n", 
                     metadata->altitude, metadata->azimuth);
}